#!/usr/bin/env perl

use strict;
use warnings;
use Carp;
use IO::Dir;
use IO::File;
use Data::Dumper;
use Perl6::Slurp;
use FindBin qw|$Bin|;
use File::Copy;
use File::Path qw|make_path|;
use File::Basename;
use Version::Compare;
use Cwd qw|abs_path|;
use Getopt::Lucid qw|:all|;
use List::MoreUtils qw|all uniq|;
use Config::Std;
use String::ShellQuote;
use Genome::Ranges::IntSet;
use Utils::Workflow;
use Utils::Hash qw|merge_conf|;
use Utils::File::Iter qw|iter_file|;


use lib "$Bin/../lib/";
use Shared qw|read_list parse_tabfile|;


############################
## Command line interface ##
############################

my @spec =  (
	Param("conf|c")->valid(sub { -r }),
	Param("list|l")->valid(sub { -r }),
	Param("tab")->valid(sub { -r }),
	Param("ped")->valid(sub { -r }),
	Param("prefix"),
	Param("rename")->valid(sub { -r }),
	Param("remove")->valid(sub { -r }),
	Param("outdir|out"),
	Param("engine|eng")->valid(sub { $_ eq 'SGE' || $_ eq 'BASH' }),
	Keypair("param|par"),
	Switch("help|h"),
	Switch("local"),
	Switch("dryrun|dry"),
	Switch("force")
	);

my $opt = Getopt::Lucid->getopt(\@spec);

if ($opt->get_help) {
	print STDERR <<EOF;
Purpose:
	This is a pipeline script used to extract information from DeepVariant VCFs.

Usage:
	dv_recall2.pl --conf Config --list VCFList [--tab VarTab] --outdir WorkDir

Options:
	--list: VCF file list or directory. The list file can have 1 or 2 columns per line. 
			If the list have 2 columns, the second column will be used as sample ID.
			Otherwise, file basename after stripping off suffix will be used as sample ID.
	--prefix : The output file prefix, default will be the basename of input table file. 
	--tab : The variant table. It overrides the file specified in the config.
	--ped : The pedigree file. It overrides the file specified in the config.
	--rename : Rename original VCF files to match sample IDs in the variants table file.
	--remove : Sample removal list. In case of renaming, new IDs should be used in removal list.

Notes:
	This script is intended to use when large number variants per individual needs DeepVariant validation,
	and we already have variant calls made by DeepVariant in VCFs or gVCFs. In this case, we can quickly 
	extract DeepVariant results for each variant. The input VCF can only have one sample per file. 

EOF
	exit 1;
}

$opt->validate({requires => [qw|conf list outdir|]});

my $rootdir = $opt->get_outdir;

my %conf = merge_conf($opt->get_conf, $opt->get_param); 
$conf{PATH}{MODULE} = shell_quote("$Bin/module");
$conf{PATH}{UTIL} = shell_quote("$Bin/../utils");

if ($opt->get_tab) {
	$conf{INPUT}{FILE} = $opt->get_tab;
}
unless (-f $conf{INPUT}{FILE}) {
	croak "Cannot find input table file: $conf{INPUT}{FILE}";
}

my $outfile;
if ($opt->get_prefix) {
	my $prefix = $opt->get_prefix;
	if ($prefix !~ /\.txt$/) {
		$outfile = $prefix.".txt";
	}
	else {
		$outfile =$prefix;
	}
}
else {
	if ($opt->get_tab) {
		$outfile = basename($opt->get_tab);
	}
	else {
		$outfile = basename($conf{INPUT}{FILE});
	}
}
if ($outfile =~ /\.gz$/) {
	$outfile =~ s/\.gz$//;
}
$conf{PATH}{OUTPUT} = $outfile;


# Input should be DeepVariant VCFs per-individual 
my ($vcfs, $filetype) = read_list($opt->get_list, 
	{ suffix => ['vcf', 'g.vcf.gz', 'vcf.gz', 'gvcf.gz'], rename => $opt->get_rename, remove => $opt->get_remove });
print STDERR "Finished reading $filetype file list\n";

# Get family ID for each sample, and samples in each family
my (%fids, %famsamps);
if ($opt->get_ped) {
	$conf{PED}{FILE} = $opt->get_ped;
}
if (exists $conf{PED}) {
	unless(-f $conf{PED}{FILE}) {
		die "Cannot find PED file $conf{PED}{FILE}";
	}
	$conf{PED}{OPTION} = "--ped $conf{PED}{FILE} ";
	%fids = map { (split)[1,0] } slurp $conf{PED}{FILE};
	while(my ($iid, $fid) = each %fids) {
		push @{$famsamps{$fid}}, $iid;
	}
	if (defined $conf{PED}{IGNORE}) {
		$conf{PED}{OPTION} .= "--ped-ignore $conf{PED}{IGNORE} ";
	}
	if (defined $conf{PED}{TWINS}) {
		$conf{PED}{OPTION} .= "--ped-twins $conf{PED}{TWINS} ";
	}
}
else {
	$conf{PED}{OPTION} = "";
}


#######################################################################
# Parse input file to prepare sample-specific variant calling intervals
#######################################################################

my $wkf = Utils::Workflow->new($rootdir,
	{ engine => $opt->get_engine, force => $opt->get_force, strict_var => 1 });

# Determine the number of families/samps that needs DeepVar validation
my (%sampbygrp, @groups);
{
	my %known;
	if (exists $conf{INPUT}{SAMPLES}) {
		open my $fin, $conf{INPUT}{SAMPLES} or die "Cannot read sample list: $conf{INPUT}{SAMPLES}";
		while(<$fin>) {
			my $iid = (split)[0];
			$known{$iid} = 1;
			if (defined $fids{$iid}) {
				foreach my $sampid (@{$famsamps{$fids{$iid}}}) {
					next if $sampid eq $iid;
					$known{$sampid} = 1;
				}
			}
		}
	}
	else {
		my ($it, $fnames, $input_fields) =
			parse_tabfile($conf{INPUT}{FILE}, $conf{INPUT}{FIELDS}, 5);
		while(my $dat = $it->()) {
			my ($iid, $chr, $pos, $ref, $alt) = @{$dat}{@$input_fields};
			$known{$iid} = 1;
			if (defined $fids{$iid}) {
				foreach my $sampid (@{$famsamps{$fids{$iid}}}) {
					next if $sampid eq $iid;
					$known{$sampid} = 1;
				}
			}
		}
	}
	
	my $nsamp = scalar(keys %known);
	my $nincl = grep { defined $vcfs->{$_} } keys %known; 
	if (exists $conf{PED}) {
		print STDERR "A total of $nsamp samples and their family members are found in the variants table, ",
			"$nincl of them have DeepVariant VCFs .\n";
	}
	else {
		print STDERR "A total of $nsamp samples found in the variants table, ",
			"$nincl of them have DeepVariant VCFs.\n";
	}

	foreach my $iid (keys %known) {
		if (defined $fids{$iid}) {
			next if defined $sampbygrp{"Fam".$fids{$iid}};
			$sampbygrp{"Fam_".$fids{$iid}} = [ grep { defined $known{$_} } @{$famsamps{$fids{$iid}}} ];
		}
		else {
			$sampbygrp{"Indiv_".$iid} = [$iid];
		}
	}

	# Creating working directory
	@groups = sort keys %sampbygrp;
	#open my $flst, ">$rootdir/par/AllGrps.txt" or die "Cannot write to AllGrps.txt"
	for(my $ii = 0; $ii < @groups; $ii ++) {
		my $jj = $ii + 1;
		open my $fout, ">$rootdir/par/GROUP.$jj" or die "Cannot write to GROUP.$jj";
		print $fout $groups[$ii], "\n";
		#print $flst $groups[$ii], "\n";
		make_path "$rootdir/wrk/$groups[$ii]";
		#make_path "$rootdir/wrk/$groups[$ii]/vcf";
		open my $flst, ">$rootdir/wrk/$groups[$ii]/keep.txt" or die "Cannot write to wrk/$groups[$ii]/keep.txt";
		foreach my $iid (@{$sampbygrp{$groups[$ii]}}) {
			if (defined $vcfs->{$iid}) {
				print $flst $iid, "\t", $vcfs->{$iid}, "\n";
			}
			else {
				print $flst $iid, "\t.\n";
			}
		}
	}
}

write_config %conf, "$rootdir/par/run.conf" unless $opt->get_dryrun;


#################################
##  Workflow initialization    ##
#################################
$wkf->add(prep_input(),
		{ name => "PrepInput", expect => [ map { "wrk/$_/input.txt" } @groups ], 
			nslots => scalar(@groups) })
	->add(dv_collect(), 
		{ name => "DVFetch",  expect => [ map { "wrk/$_/output.txt" } @groups ], 
			nslots => scalar(@groups), depend => "PrepInput", deparray => 1 })
	->add(merge_res(),
		{ name => "MergeRes", expect => ["out/$outfile.gz.tbi", "out/$outfile.gz.tbi"], depend => "DVFetch" });

$wkf->inst(\%conf);
$wkf->run({ conf => $conf{$wkf->{engine}}, dryrun => $opt->get_dryrun  });


############################
## Workflow components    ##
############################

# Prepare input files for DeepVariant
sub prep_input {
	my $script = <<'EOF';

read GRPID < _PARDIR_/GROUP._INDEX_

perl _PATH.MODULE_/prep_dvinput.pl _INPUT.FILE_ _INPUT.FIELDS_  _WRKDIR_/$GRPID

EOF
}


# Collect DeepVariant calls from VCFs directly
sub dv_collect {
	my $script = <<'EOF';

read GRPID < _PARDIR_/GROUP._INDEX_

perl _PATH.MODULE_/sum_dvrecall.pl \
	--input _WRKDIR_/$GRPID/input.txt --fields _INPUT.FIELDS_ \
	--fasta _PATH.FASTA_ --dvout _WRKDIR_/$GRPID/vcf \
	--output _WRKDIR_/$GRPID/output.txt _PED.OPTION_ \
	--fields-add _OUTPUT.OUTFIELDS_ --nearby _OUTPUT.NEARBY_

EOF
}

# Merge results
sub merge_res {
	my $script = <<'EOF';

FIRSTTAB=
find _WRKDIR_ -name 'output.txt' | while read TABLE; do
	if [[ -z $FIRSTTAB ]]; then 
		cat $TABLE
		FIRSTTAB=1
	else
		awk 'NR>1' $TABLE
	fi 
done > _OUTDIR_/_PATH.OUTPUT_

perl _PATH.UTIL_/tabix_tabfile.pl _OUTDIR_/_PATH.OUTPUT_

EOF
}



