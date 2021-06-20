#!/usr/bin/env perl

use strict;
use warnings;
use Carp;
use IO::Dir;
use Data::Dumper;
use FindBin qw|$Bin|;
use File::Basename qw|basename dirname|;
use File::Path qw|make_path|;
use Hash::Util qw|lock_keys|;
use List::MoreUtils qw|all uniq|;
use String::ShellQuote;
use Config::Std;
use Perl6::Slurp;
use Getopt::Lucid qw|:all|;
use Utils::File qw|count_line|;
use Utils::File::Iter qw|iter_file|;
use Utils::Hash qw|merge_conf|;
use Utils::Workflow;

use lib "$Bin/../lib";
use Shared qw|check_conf parse_fstr|;
use Variants qw|anno_fields|;


#############################
## Command line interface  ##
#############################

my @spec =  (
	Param("conf|c")->valid(sub { -r }),
	Param("input|in")->valid(sub { -r }),
	Param("outdir|out"),
	Param("engine|eng")->valid(sub { $_ eq 'SGE' || $_ eq 'BASH' }),
	Keypair("param|par"),
	Switch("help|h"),
	Switch("dryrun|dry"),
	Switch("force")
	);

my $opt = Getopt::Lucid->getopt(\@spec);

if ($opt->get_help) {
	print STDERR <<EOF;
Purpose:
	This is a pipeline script to annotate large number of sequence variants.

Usage:
	annovar_pipeline.pl --conf Config [--input InFile] --outdir OutDir

Notes:
	Input to this pipeline can be a variant table file or a directory of variant table files.
	Large input tables will be split into smaller pieces before annotation. Users can specify 
	one or more annotation config files that are used by anno_seqvars.pl, the first config will 
	be used for the primary annotation which will include extra fields from input and collect
	sample level information. Other configs will be used for the secondary variant-only annotations
	used to add extra varinat level info fields that cannot be annotated in the first config. 
	When multiple annotation configs are provided, annotations will be done in parallel and 
	the outputs will be merged by adding new fields from secondary output to the first primary 
	annotation output.

	Before annotation, variants can optionally first be normalized or lifted over to a new genome 
	assembly. If liftover was performed, the annotation configs must be based on the new genome 
	assembly. After annotation, additional variant, gene or sample level information that cannot 
	be gathered by anno_seqvars.pl can be added by one or more post-processing scripts. Post-
	processings scripts take previously annotated variant table as input and output a new variant
	table with additional fields. They will be run sequentially by the order they appear in the 
	config, script appear later can use all the fields from the previous output.

	The final output will be merged from all splits, and sorted by chromosome and position.
	It can also be bgzip'ed and tabix index for position based query.

EOF
	exit 1;
}

$opt->validate({requires => [qw|conf outdir|]});

my %conf = merge_conf($opt->get_conf, $opt->get_param); 
check_conf(\%conf);
$conf{PATH} = { BIN => shell_quote($Bin), MODULE => shell_quote("$Bin/module") };

#############################
## Input files validation  ##
#############################

my (@infiles, @labels);
if ($opt->get_input) {
	$conf{INPUT}{FILE} = $opt->get_input;
}
if(!defined $conf{INPUT}{FILE}) {
	die "Must provide input file!";
}
elsif (-f $conf{INPUT}{FILE}) {
	@infiles = ($conf{INPUT}{FILE});
}
elsif (-d $conf{INPUT}{FILE}) {
	@infiles = map { $conf{INPUT}{FILE}."/".$_ } grep { /\.txt$/ || /\.txt\.gz$/ } 
		IO::Dir->new($conf{INPUT}{FILE})->read();
}
unless(@infiles) {
	die "Cannot find input file: $conf{INPUT}{FILE}!";
}
else {
	if (defined $conf{INPUT}{RENAME}) {
		$conf{INPUT}{RENAME} =~ s/^['"]//; $conf{INPUT}{RENAME} =~ s/['"]$//
	}
	foreach my $infile (@infiles) {
		(my $fbase = basename($infile)) =~ s/\.txt(\.gz)?$//;
		if (defined $conf{INPUT}{RENAME}) {
			my $expr = '$fbase =~ '.$conf{INPUT}{RENAME};
			eval($expr); 
			if ($@) { 
				warn "Error rename $infile!" 
			}
		}
		push @labels, $fbase;
	}
	my @uqlabs = uniq sort @labels;
	unless(scalar(@uqlabs) == scalar(@labels)) {
		die "Input files are not unique!";
	}
}

my $rootdir = $opt->get_outdir;

my $wkf = Utils::Workflow->new($rootdir,
	{ engine => $opt->get_engine, force => $opt->get_force });

# Deterine and check key fields from input file
# Also determine the number of splits
my %groups;
{
	open my $fgrp, ">$rootdir/par/GROUPS" or die "Cannot write to GROUPS file";
	open my $finp, ">$rootdir/par/INPUTS" or die "Cannot write to INPUTS file";

	my $keyfields = parse_fstr($conf{INPUT}{FIELDS}, 1);
	my @keyfields = keys %$keyfields;
	unless(@keyfields == 4) {
		die "Must provide four key fields in input file for Chrom,Position,Ref,Alt!";
	}
	my @xtrafields;
	if (defined $conf{INPUT}{XTRAFIELDS}) {
		my $xtrafields = parse_fstr($conf{INPUT}{XTRAFIELDS}, 1);
		@xtrafields = keys %$xtrafields;
	}
	my @oldvarfields;
	if (defined $conf{LIFTOVER} && $conf{LIFTOVER}{OLDVAR}) {
		@oldvarfields = split(",", $conf{LIFTOVER}{OLDVAR});
	}

	my $ct = 1;
	for(my $ii = 0; $ii < @infiles; $ii ++) {
		my ($it, $fnames) = iter_file($infiles[$ii], { fsep => qr/\t/ });
		foreach my $field (@keyfields, @xtrafields) {
			unless(grep { $_ eq $field } (@$fnames, @oldvarfields) ) {
				die "Cannot find field $field in $infiles[$ii]!";
			}
		}
		if (defined $conf{PARAM}{NSPLIT}) {
			$groups{$labels[$ii]} = [$ct, $ct + $conf{PARAM}{NSPLIT} - 1];
			print $finp $infiles[$ii], "\t", $ct, "\t", $conf{PARAM}{NSPLIT}, "\n";
			$ct += $conf{PARAM}{NSPLIT};
		}
		elsif (defined $conf{PARAM}{NLINE}) {
			my $totline = count_line($infiles[$ii]); $totline --;
			my $nsplit = int($totline/$conf{PARAM}{NLINE}+0.5);
			$groups{$labels[$ii]} = [$ct, $ct + $nsplit - 1];
			print $finp $infiles[$ii], "\t", $ct, "\t", $nsplit, "\n";
			$ct += $nsplit;
		}	
		else {
			die "Must provide number of splits or number of lines per split for input!";
		}
		print $fgrp join("\t", $labels[$ii], @{$groups{$labels[$ii]}}), "\n"; 
	}
	$conf{PARAM}{TOTSPLIT} = $ct - 1;
	$conf{INPUT}{FIELDS} = join(',', @keyfields);
	if (join(',', values %$keyfields) eq 'Chrom,Position,Ref,Alt') {
		$conf{INPUT}{ALIAS} = join(',', map { "$_:$keyfields->{$_}" } @keyfields);
	}
	else {
		my @alias = qw|Chrom Position Ref Alt|;
		$conf{INPUT}{ALIAS} = join(',', map { "$keyfields[$_]:$alias[$_]" } 0..$#keyfields);
	}
	if (defined $conf{PARAM}{BGZIP} && $conf{PARAM}{BGZIP} =~ /^(Y|T)/i) {
		$conf{PARAM}{BGZIP} = "--bgzip";
	}
	else {
		$conf{PARAM}{BGZIP} = "";
	}
}


# Read configs for annotation and determine output field names
# Can used later to check post-processing output files
my @outfields;
{
	# Processing anno config files
	my @annofields;
	if (ref $conf{ANNO}{CONF} eq 'ARRAY') {
		for(my $ii = 0; $ii < @{$conf{ANNO}{CONF}}; $ii ++) {
			my $jj = $ii + 1;
			my %anno;
			read_config $conf{ANNO}{CONF}[$ii] => %anno;
			unless($anno{Input}{HG} eq $conf{ANNO}{HG}) {
				die "Incorrect genome assembly build in config: $anno{Input}{HG} <> $conf{ANNO}{HG}";
			}
			$anno{Input}{Fields} = $conf{INPUT}{ALIAS};
			if ($ii == 0 && defined $conf{INPUT}{XTRAFIELDS}) {
				$anno{Input}{XtraFields} = $conf{INPUT}{XTRAFIELDS};
			}
			else {
				delete $anno{Sample};
				delete $anno{Input}{XtraFields};
			}
			my @conffields = anno_fields(\%anno);
			foreach my $field (@conffields) {
				unless(grep { $field eq $_ } @annofields) {
					push @annofields, $field;
				}
			}
			write_config %anno => "$rootdir/par/Anno_$jj.conf";
		}
		$conf{PARAM}{NUMANNO} = scalar(@{$conf{ANNO}{CONF}});
	}
	else {
		#$conf{ANNO}{CONF} = [$conf{ANNO}{CONF}];
		read_config $conf{ANNO}{CONF} => my %anno;
		unless($anno{Input}{HG} eq $conf{ANNO}{HG}) {
			die "Incorrect genome assembly build in config: $anno{Input}{HG} <> $conf{ANNO}{HG}";
		}
		$anno{Input}{Fields} = $conf{INPUT}{ALIAS};
		if (defined $conf{INPUT}{XTRAFIELDS}) {
			$anno{Input}{XtraFields} = $conf{INPUT}{XTRAFIELDS};
		}
		@annofields = anno_fields(\%anno);
		write_config %anno => "$rootdir/par/Anno.conf"
	}
	$outfields[0] = \@annofields;
	if (defined $conf{ANNO}{POST}) {
		unless(ref $conf{ANNO}{POST}) {
			$conf{ANNO}{POST} = [$conf{ANNO}{POST}];
		}
		unless(ref $conf{ANNO}{FIELDS}) {
			$conf{ANNO}{FIELDS} = [$conf{ANNO}{FIELDS}];
		}
		unless(@{$conf{ANNO}{POST}} == @{$conf{ANNO}{FIELDS}}) {
			die "Different number of post-processing scripts and fields!";
		}
		$conf{PARAM}{NUMPOST} = scalar(@{$conf{ANNO}{POST}});
		for(my $ii = 0; $ii < @{$conf{ANNO}{FIELDS}}; $ii ++) {
			unless(-f $conf{ANNO}{POST}[$ii]) {
				die "Cannot find post-processing script $conf{ANNO}{POST}[$ii]!";
			} 
			my $pfields = parse_fstr($conf{ANNO}{FIELDS}[$ii], 1);
			foreach my $field (values %$pfields) {
				if (grep { $_ eq $field} @annofields) {
					die "Field $field from $conf{ANNO}{POST}[$ii] already exists in annotation!";
				}
			}
			$outfields[$ii+1] = [ @{$outfields[$ii]}, values %$pfields ];
		}
	}
}


#################################
##  Workflow initialization &  ##
##  Working directory setup    ##
#################################

my %opt;
if (defined $conf{LIFTOVER} || defined $conf{INPUT}{NORM}) {
	$wkf->add(split_tab(), { name => "SplitTab", 
							expect => [ map { "tmp/Input.$_.txt" } 1..$conf{PARAM}{TOTSPLIT} ] });
	if (defined $conf{LIFTOVER}) {
		$wkf->add(lift_over(), { name => "LiftOver", depend => "SplitTab", nslots => $conf{PARAM}{TOTSPLIT},
								 expect => [ map { ["wrk/Input.$_.txt", "wrk/Vars.$_.txt"] } 1..$conf{PARAM}{TOTSPLIT} ] });
		$opt{depend} = "LiftOver";
	}
	elsif (defined $conf{INPUT}{NORM}) {
		$wkf->add(norm_vars(), { name => "NormVars", depend => "SplitTab", nslots => $conf{PARAM}{TOTSPLIT},
								 expect => [ map { ["wrk/Input.$_.txt", "wrk/Vars.$_.txt"] } 1..$conf{PARAM}{TOTSPLIT} ] });
		$opt{depend} = "NormVars";
	}
}
else {
	$wkf->add(split_tab(), { name => "SplitTab", 
							 expect => [ map { "wrk/Input.$_.txt" } 1..$conf{PARAM}{TOTSPLIT} ] })
		->add(extract_vars(), { name => "ExtractVars", depend => "SplitTab", nslots => $conf{PARAM}{TOTSPLIT},
								expect => [ map { "wrk/Vars.$_.txt" } 1..$conf{PARAM}{TOTSPLIT} ]  } );
	$opt{depend} = "ExtractVars";

}
if (ref $conf{ANNO}{CONF} eq 'ARRAY') {
	for(my $ii = 0; $ii < @{$conf{ANNO}{CONF}}; $ii ++) {
		my $jj = $ii + 1;
		$wkf->add(anno_vars($jj), { name => "AnnoVars_$jj", %opt, nslots => $conf{PARAM}{TOTSPLIT},
								   	expect => [ map { "wrk/Anno_$jj.$_.txt" } 1..$conf{PARAM}{TOTSPLIT} ],
									callback => \&check_anno })
	}
	$wkf->add(merge_anno(), { name => "MergeAnno", nslots => $conf{PARAM}{TOTSPLIT},
							  depend => [ map { "AnnoVars_$_" } 1..$conf{PARAM}{NUMANNO} ],
							  expect => [ map { "wrk/Anno.$_.txt" } 1..$conf{PARAM}{TOTSPLIT} ],
							  callback => \&check_anno });
	$opt{depend} = "MergeAnno";
}
else {
	$wkf->add(anno_vars(), { name => "AnnoVars", %opt, nslots => $conf{PARAM}{TOTSPLIT},
							 expect => [ map { "wrk/Anno.$_.txt" } 1..$conf{PARAM}{TOTSPLIT} ],
							 callback => \&check_anno });
	$opt{depend} = "AnnoVars";
}
if (defined $conf{ANNO}{POST}) {
	for(my $ii = 0; $ii < @{$conf{ANNO}{POST}}; $ii ++) {
		my $jj = $ii + 1;
		my @exps;
		if ($jj == $conf{PARAM}{NUMPOST}) {
			@exps = map { "wrk/Output.$_.txt" } 1..$conf{PARAM}{TOTSPLIT};
		}
		else {
			@exps = map { "wrk/Output_$jj.$_.txt" } 1..$conf{PARAM}{TOTSPLIT};
		}
		$wkf->add(post_proc($jj), { name => "PostProc_$jj", %opt, nslots => $conf{PARAM}{TOTSPLIT}, 
								 	expect => \@exps, callback => \&check_anno } );
		$opt{depend} = "PostProc_$jj";
	}
}
$wkf->add(merge_final(), { name => "MergeFinal", %opt, 
						   expect => $conf{PARAM}{BGZIP} ? 
						   		[ map { "out/$_.txt.gz" }  @labels] : [ map { "out/$_.txt" }  @labels] });


$wkf->inst(\%conf);
$wkf->run({ conf => $conf{$wkf->{engine}}, dryrun => $opt->get_dryrun });


############################
## Workflow components    ##
############################

# Split input variant table, and prepare site-only table
# => tmp or wrk: Input.XX, Var.XX
sub split_tab {
	my $script;
	if (defined $conf{LIFTOVER} || defined $conf{INPUT}{NORM}) {
		# If liftover or normlization is needed, write to tmpdir first
		$script = <<'EOF';

cat _PARDIR_/INPUTS | while read INFILE START NSPLIT; do
	split_tab.pl --in $INFILE --output _TMPDIR_/Input --nsplit $NSPLIT --start $START --suffix txt
done

EOF
	}
	else {
		$script = <<'EOF';

cat _PARDIR_/INPUTS | while read INFILE START NSPLIT; do
	split_tab.pl --in $INFILE --output _WRKDIR_/Input --nsplit $NSPLIT --start $START --suffix txt
done

EOF
	}
	return $script;
}

sub extract_vars {
	my $script = <<'EOF';

# To be more robust to err, first extract variants to tmpdir, then move.

flatdb.pl --in _WRKDIR_/Input._INDEX_.txt -fsep "\t" --select _INPUT.FIELDS_ | body sort -u  > _TMPDIR_/Vars._INDEX_.txt

mv _TMPDIR_/Vars._INDEX_.txt _WRKDIR_/Vars._INDEX_.txt

EOF
}

sub norm_vars {
	my $script = <<'EOF';

norm_vars.pl --in _TMPDIR_/Input._INDEX_.txt --fields _INPUT.FIELDS_ \
	--seq _INPUT.NORM_ --fsep "\t" --out _WRKDIR_/Input._INDEX_.txt

flatdb.pl --in _WRKDIR_/Input._INDEX_.txt -fsep "\t" --select _INPUT.FIELDS_ | body sort -u  > _TMPDIR_/Vars._INDEX_.txt

mv _TMPDIR_/Vars._INDEX_.txt _WRKDIR_/Vars._INDEX_.txt

EOF
}


# => wrk: Input.XX, Var.XX
# Note: Some variants may not be lifted over to the target assembly.
sub lift_over {
	my $script = <<'EOF';

liftover_vars.pl --in _TMPDIR_/Input._INDEX_.txt --fields _INPUT.FIELDS_ \
	--chain _LIFTOVER.CHAIN_ --seq _LIFTOVER.FASTA_ --fsep "\t" --no-chr \
	--old-var _LIFTOVER.OLDVAR_ --out _WRKDIR_/Input._INDEX_.txt

flatdb.pl --in _WRKDIR_/Input._INDEX_.txt -fsep "\t" --select _INPUT.FIELDS_ | body sort -u > _TMPDIR_/Vars._INDEX_.txt

mv _TMPDIR_/Vars._INDEX_.txt _WRKDIR_/Vars._INDEX_.txt

EOF
}

# Annotate variants 
# The primary config will include samp and external info
# secondary configs will process variants only
sub anno_vars {
	my ($index) = @_;
	my $script = << 'EOF';

perl _PATH.BIN_/anno_seqvars.pl --conf _PARDIR_/NAME.conf --in _WRKDIR_/INPUT._INDEX_.txt \
	--wrkdir _TMPDIR_/NAME._INDEX_ --out _WRKDIR_/NAME._INDEX_.txt

EOF
	if (defined $index) {
		if ($index == 1) {
			$script =~ s/INPUT/Input/;
		}
		else {
			$script =~ s/INPUT/Vars/;
		}
		$script =~ s/NAME/Anno_$index/g;
	}
	else {
		$script =~ s/INPUT/Input/;
		$script =~ s/NAME/Anno/g;
	}
	return $script;
}

# Merge annotation if needed
sub merge_anno {
	my $script = << 'EOF';

perl _PATH.MODULE_/merge_anno.pl _WRKDIR_/Anno_{1.._PARAM.NUMANNO_}._INDEX_.txt _WRKDIR_/Anno._INDEX_.txt

EOF
}

# Run post-processing script
sub post_proc {
	my ($index) = @_;
	my $script = <<'EOF';

bash SCRIPT _WRKDIR_/INPUT._INDEX_.txt _WRKDIR_/OUTPUT._INDEX_.txt

EOF
	if (defined $index) {
		$script =~ s/SCRIPT/$conf{ANNO}{POST}[$index-1]/;
		if ($index == $conf{PARAM}{NUMPOST}) {
			$script =~ s/OUTPUT/Output/;
		}
		else {
			$script =~ s/OUTPUT/Output_$index/;
		}
		if ($index == 1) {
			$script =~ s/INPUT/Anno/;
		}
		else {
			$script =~ s/INPUT/Output_@{[ $index-1 ]}/;
		}
	}
	else {
		$script =~ s/INPUT/Anno/;
		$script =~ s/OUTPUT/Output/;
	}

	return $script;
}

# Merge final results
sub merge_final {
	my $script = <<'EOF';

perl _PATH.BIN_/../utils/merge_var_splits.pl --input _WRKDIR_/INPUT --outdir _OUTDIR_ \
	--group _PARDIR_/GROUPS --tmpdir _TMPDIR_ _PARAM.BGZIP_

EOF
	if (defined $conf{ANNO}{POST}) {
		$script =~ s/INPUT/Output/;
	}
	else {
		$script =~ s/INPUT/Anno/;
	}
	return $script;
}


# Check annotation output have the same number of lines as input
sub check_anno {
	my ($exp) = @_;
	if (-f "$rootdir/$exp") {
		my $fbase = basename($exp);
		my ($input, $index, $suffix);
		if ($fbase =~ /^Anno\.(\d+)\.txt$/) {
			$input = "wrk/Input";
			$index = $1;
		}
		elsif ($fbase =~ /^Anno_(\d+)\.(\d+)\.txt$/) {
			($suffix, $index) = ($1, $2); 
			if ($suffix == 1) {
				$input = "wrk/Input";
			}
			else {
				$input = "wrk/Vars";
			}
		}
		elsif ($fbase =~ /^Output(?:_\d+)?\.(\d+)\.txt$/) {
			$index = $1;
			$input = "wrk/Anno";
		}
		else {
			die "Cannot determine input from output: $exp!";
		}
		my $n_in = count_line("$rootdir/$input.$index.txt");
		my $n_out = count_line("$rootdir/$exp");
		if ($n_in == $n_out) {
			# Also check header for the first file in output, assuming all other files have the same header
			if ($index == 1 && !defined $suffix) {
				my ($it, $fnames) = iter_file("$rootdir/$exp", { fsep => qr/\t/ });
				my $header1;
				my $header2;
				if ($fbase =~ /^Anno/) {
					$header1 = join(",", @$fnames);
					$header2 = join(",", @{$outfields[0]});
				}
				elsif ($fbase =~ /^Output_(\d+)/) {
					my $ii = $1;
					$header1 = join(",", sort @$fnames);
					$header2 = join(",", sort @{$outfields[$ii]});
				}
				elsif ($fbase =~ /^Output/) {
					$header1 = join(",", sort @$fnames);
					$header2 = join(",", sort @{$outfields[-1]});
				}
				else {
					die "Cannot recognize file to check header: $fbase";
				}
				unless($header1 eq $header2) {
					warn "File $rootdir/$exp does not have expected header!";
					print STDERR ">>>Observed: ", $header1, "\n", ">>>Expected: ", $header2, "\n";
				}
			}
			return 1;
		}
		else {
			return 0;
		}
	}
	else {
		return 0;
	}
}


