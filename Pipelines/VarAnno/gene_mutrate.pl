#!/usr/bin/env perl
use strict;
use warnings;
use FindBin qw|$Bin|;
use List::Util qw|sum|;
use List::MoreUtils qw|all|;
use Config::Std;
use IO::Dir;
use File::Temp qw|tempdir|;
use Hash::Util qw|lock_keys|;
use Getopt::Lucid qw|:all|;
use Utils::Parser qw|sql_query|;
use Utils::File::Iter qw|iter_file|;

use lib "$Bin/../lib";
use Shared qw|parse_filters slurp_xref expand_dat check_conf|;
use Variants qw|mutrater|;


############################
## Command line interface ##
############################

my @spec =  (
	Param("conf|c")->valid(sub { -r }),
	Param("input|in|i")->valid(sub { -r }),
	Param("wrkdir|wrk|w"),
	Param("output|out|o"),
	Switch("force"), # Under force mode, previously generated files will not be over-written
	Switch("help|h")
	);

my $opt = Getopt::Lucid->getopt(\@spec);

if ($opt->get_help) {
	print STDERR <<EOF;
Purpose:
	This is a standalone scripts to calculate aggregated haploid mutation rates for different classes 
	of variants of each gene. It can also be used to summarize the variant counts/rates in each gene. 

Usage:
	gene_mutrate.pl --conf Config --in Input [--wrk WrkDir] --out Output

Options:
	--in: (Required) Input file should be annotation of all possible enumeration of SNVs in genes or 
		  exome. It can also be a directory containing one table per gene. The input can be generated 
		  using enum_cdssnvs.pl. 
	--wrk: When input is a directory, the working directory will be used to store intermediate
		   results from each individual input table. When it is not specified, a temp directory will be used.
	--force: When final or intermediate output files already exist, it will be skipped. Turn on this switch 
			 to overwrite existing output.  

Notes:
	The mutation rate for a class of variants is an aggregation of point mutation rates of all possible
	variants of the same class. For SNVs, the point mutation rates depend mainly on local sequence context.
	We can make use of published 3mer (Samocha K et al. 2014) or 7mer (Carlson J et al. 2019) context
	dependent rates. Mutation rates are also known to vary at mega-base scale across the genomic regions,
	and correlated with different genomic annotations including DNA methylation, chromatin modifications,
	and recombination rates. To accomondate region-specific variation, it is possible to provide scaling
	factors for different genomic regions on top of the local mutation rate (Samocha K et al. 2014).
	For more complex mutation rate models (e.g. Carlson J et al. 2019), it is possible to store predicted 
	per-base mutation rate by the model in a tabix database file for position-based querying. All these
	options are available in config file. 

	For indels, their mutation rates are typically derived from the aggregated rates of some class of SNVs.
	For example, we usually assume the rate of frameshift indels to be 1.3x of nonsense (stop gained) SNVs.

	By default, this script groups mutation rates by individual genes identified by GeneID column. Mutation
	rates for gene sets or whole exome can then be calculated by taking summation across genes, because 
	overlapping between genes are generally negligible after excluding read-through genes. 


Input/output:
	The key input to this script is the annotation of all possible enumeration of SNVs.
	An example can be found in Dropbox/Data/AnnoDB/ExomeSNVs/GencodeV19, which includes annotation of
	all coding SNVs in 19944 coding genes. When input is a directory, the script will run in parallel
	to process each table file in the directory as input (using GNU parallel, SGE not supported), then 
	the final results are merged from the output of individual runs. 

	The config file should specify the method to get per-bp mutation rate for each alternative allele at 
	each genomic positions. It also defines variant classes based on the information from input annotation 
	table. The output will have one line for each gene including gene ID, extra gene level information, 
	and mutation rates for different classes of varaints. 

	The input file can also be any list of annotated variants. The mutation rate method section is config 
	file can also be removed, in such case the script will simply count the observed variants and 
	summarize the gene level weighted count for each defined variant class.  

EOF
	exit 1;
}

$opt->validate({ requires => [qw|conf input output|] });


my $infile = $opt->get_input;
my $outfile = $opt->get_output;
my $config = $opt->get_conf;

if (-f $outfile) {
	unless($opt->get_force) {
		print STDERR "Output file $outfile already exists!\n";
		exit 1;
	}
}

read_config $config => my %conf;
check_conf(\%conf);
lock_keys(%conf);

if (-d $infile) {
	my @vartabs = sort grep { !/^\./ && -f "$infile/$_" } IO::Dir->new($infile)->read();
	# Working directory must be provided when doing parallization
	my $wrkdir = $opt->get_wrkdir;
	if ($wrkdir) {
		mkdir $wrkdir unless -d $wrkdir;
	} else {
		$wrkdir = tempdir(CLEANUP => 1);
	}

	my $option = "";
	if (exists $conf{Parallel} && defined $conf{Parallel}{jobs}) {
		$option = join(" ", map { "--$_ $conf{Parallel}{$_}" } keys %{$conf{Parallel}});
	}
	else {
		die "Must provide number of threads in config file for parallel";
	}
	my $args = "";
	if ($opt->get_force) {
		$args = "--force";
	}

	my $command = qq(find $infile/ -not -path '*/\.*' -type f | parallel --eta $option "perl $0 --in {} --conf $config --out $wrkdir/{/}.txt $args 2>$wrkdir/{/}.err");
	print $command, "\n";
	system($command);

	# Then merge the output
	open my $fout, ">$outfile" or die "Cannot write to $outfile";
	for(my $ii = 0; $ii < @vartabs; $ii ++) {
		unless(-f "$wrkdir/$vartabs[$ii].txt") {
			die "Cannot find split $ii output $vartabs[$ii] in the working directory";
		}
		open my $fin, "$wrkdir/$vartabs[$ii].txt" or die "Cannot open $vartabs[$ii]";
		<$fin> if $ii > 0;
		while(<$fin>) {
			print $fout $_;
		}
	}
	exit 0;
}

# Mutation rate calculator
my $rater;
if (exists $conf{MutRate}) {
	$rater = mutrater($conf{MutRate});
}

# Parse filters for varinat classes
my ($varclass, $callbacks, $fields) = parse_filters($conf{Variant}{Class}, $conf{Variant}{Class_Filter});

# Default field for "grouping" will be GeneID
my $grpid = "GeneID";
if (exists $conf{Input} && defined $conf{Input}{GroupID}) {
	$grpid = $conf{Input}{GroupID};
}
# Optional "weight" field in input
my $wtfd;
if (exists $conf{Input} && defined $conf{Input}{Weight}) {
	$wtfd = $conf{Input}{Weight};
}

# Check fields in the input table file
my ($it, $fnames) = iter_file($infile, { fsep => qr/\t/ });
foreach my $field (keys %$fields) {
	unless(grep { $_ eq $field } @$fnames) {
		die "Cannot find field $field in the input file";
	}
}
foreach my $field (qw|Chrom Position Ref Alt|, $grpid) {
	unless(grep { $_ eq $field } @$fnames) {
		die "Cannot find required field $field in the input file";
	}
}
if (defined $wtfd) {
	unless(grep { $_ eq $wtfd } @$fnames) {
		die "Cannot find weight field $wtfd in the input file";
	}
}

my $filter;
if (exists $conf{Input} && defined $conf{Input}{Filter}) {
	my $filtexpr = $conf{Input}{Filter};
	$filtexpr =~ s/^["']//; $filtexpr =~ s/["']$//;
	($filter, my $tokens) = sql_query($filtexpr, 1);
	foreach my $tok (@$tokens) {
		if ($tok->[0] eq 'FIELD') {
			unless(grep { $tok->[1] eq $_} @$fnames) {
				die "Cannot find Filter field $tok->[1] in table $infile";
			}
		}
	}
}

# Tally gene variants to calculate mutation rate
my (%genert, %known);
while(my $dat = $it->()) {
	my $varid = join(":", @{$dat}{qw|Chrom Position Ref Alt|});
	next if defined $known{$varid};
	#my $mrt = $rater->($dat);
	my $mrt;
	if (defined $rater) {
		$mrt = $rater->($dat);
	}
	elsif (defined $wtfd) {
		$mrt = $dat->{$wtfd};
	}
	else {
		$mrt = 1;
	}
	$known{$varid} = 1;
	if ($dat->{$grpid} =~ /;/) {
		foreach my $data (expand_dat($dat, { sep => ';', optional => [$grpid, keys %$fields] })) {
			if (defined $filter) {
				next unless $filter->($data);
			}
			for(my $ii = 0; $ii < @$varclass; $ii ++) {
				if ($callbacks->[$ii]->($data)) {
					$genert{$data->{$grpid}}{$varclass->[$ii]} += $mrt;
				}
			}
		}
	}
	else {
		if (defined $filter){
			next unless $filter->($dat);
		}
		for(my $ii = 0; $ii < @$varclass; $ii ++) {
			if ($callbacks->[$ii]->($dat)) {
				$genert{$dat->{$grpid}}{$varclass->[$ii]} += $mrt;
			}
		}
	}
}

# Parse and calculate combined fields
my (@combs, @combdefs);
if (defined $conf{Variant}{Comb}) {
	if (ref $conf{Variant}{Comb} eq 'ARRAY') {
		@combs = @{$conf{Variant}{Comb}};
		unless( defined $conf{Variant}{Comb_Define} &&
				ref $conf{Variant}{Comb_Define} eq 'ARRAY' &&
				scalar(@{$conf{Variant}{Comb_Define}}) == @combs ) {
			die "Incorrect number of comb field definitions";
		}
		@combdefs = @{$conf{Variant}{Comb_Define}};
	}
	else {
		@combs = ($conf{Variant}{Comb});
		unless (defined $conf{Variant}{Comb_Define} ) {
			die "Cannot find definition for comb field";
		}
		else {
			if (ref $conf{Variant}{Comb_Define} eq 'ARRAY') {
				die "Incorrect number of comb field definition";
			}
		}
		@combdefs = ($conf{Variant}{Comb_Define});
	}
	for(my $ii = 0; $ii < @combs; $ii ++) {
		if (grep { $_ eq $combs[$ii] } @$varclass) {
			die "Combined class $_ was already defined";
		}
		my @terms = split('\+', $combdefs[$ii]);
		my (@coeff, @fields);
		foreach my $term (@terms) {
			my @spterm = split('\*', $term);
			if (@spterm == 1) {
				push @coeff, 1;
				push @fields, $spterm[0];
			}
			elsif (@spterm == 2) {
				unless($spterm[0] =~ /^[0-9\.]+$/) {
					die "Incorrect format of coefficient in combination term: $term";
				}
				push @coeff, $spterm[0];
				push @fields, $spterm[1];
			}
			else {
				die "Cannot split $term to find coefficient and field";
			}
		}
		foreach my $field (@fields) {
			unless (grep { $_ eq $field } @$varclass) {
				die "Cannot find $field from defined variant class";
			}
		}
		foreach my $gid (keys %genert) {
			foreach my $field (@fields) {
				unless(defined $genert{$gid}{$field}) {
					$genert{$gid}{$field} = 0;
				}
			}
			$genert{$gid}{$combs[$ii]} = 
				sum(map { $coeff[$_]*$genert{$gid}{$fields[$_]} } 0..$#fields);
		}
	}
}

if (defined $conf{Output}{FieldPref} && $conf{Output}{FieldPref} =~ /_$/) {
	$conf{Output}{FieldPref} =~ s/_$//;
}

my @header = map { defined $conf{Output}{FieldPref} ? $conf{Output}{FieldPref}."_".$_ : $_ } (@$varclass, @combs);

my ($xtraginfo, @gfields);
if (defined $conf{Output}{GXref}) {
	my ($xgdat, $xgfds) = slurp_xref($conf{Output}{GXref}, $conf{Output}{GXref_Fields});
	$xtraginfo = $xgdat;
	@gfields = @$xgfds;
	foreach my $gfield (@gfields) {
		if (grep { $gfield eq $_ } @header) {
			die "Field $gfield already appear in the genetab output!";
		}
	}
}

open my $fout, ">$outfile" or die "Cannot write to $outfile";
print $fout join("\t", $grpid, @gfields, @header), "\n";
my @allvclass = (@$varclass, @combs);
foreach my $gid (sort keys %genert) {
	my @ginfo;
	if (defined $xtraginfo) {
		@ginfo = map { $xtraginfo->{$gid}{$_} // "." } @gfields;  
	}
	# Genes with all zero-mutation rate will be skipped
	if (all { !defined $genert{$gid}{$_}  ||  $genert{$gid}{$_} == 0 } @allvclass) {
		next;
	}

	if (exists $conf{Output}{Log10} && $conf{Output}{Log10} =~ /^Y/i) {
		print $fout join("\t", $gid, @ginfo, 
			map { !defined $genert{$gid}{$_} || $genert{$gid}{$_} == 0 ? 
				$conf{Output}{NAstring} : log($genert{$gid}{$_})/log(10) } @allvclass), "\n";
	}
	else {
		print $fout join("\t", $gid, @ginfo, map { $genert{$gid}{$_} // 0 } @allvclass), "\n";
	}
}

