#!/usr/bin/env perl

use strict;
use warnings;
use Carp;
use Data::Dumper;
use FindBin qw|$Bin|;
use File::Basename qw|dirname|;
use File::Path qw|make_path|;
use Hash::Util qw|lock_keys|;
use List::MoreUtils qw|all|;
use String::ShellQuote;
use Config::Std;
use Perl6::Slurp;
use Getopt::Lucid qw|:all|;
use Utils::File qw|count_line|;
use Utils::Workflow;

use lib "$Bin/../utils";
use Shared qw|check_conf|;


#################################
## Command line interface      ##
## Config parsing & validation ##
#################################

my @spec =  (
	Param("outdir|o"),
	Param("conf|c")->valid(sub { -r }),
	Param("engine|eng")->valid(sub { $_ eq 'SGE' || $_ eq 'BASH' }),
	Switch("help|h"),
	Switch("force"),
	Switch("dryrun|dry"),
	);

my $opt = Getopt::Lucid->getopt(\@spec);

if ($opt->get_help) {
	print STDERR <<EOF;

Purpose:
	This is a pipeline script for enumerating and annotating all possible exome SNVs 
	("in silico saturation mutagenesis").

Usage: 
	enum_cdssnvs.pl --conf Conf --outdir OutDir

Notes:
	The coding territory of exome is defined by a gene table from UCSC genome browser.
	A small part of intronic regions can be included to cover splice sites. All possible
	SNVs are enumerated by mutating each base in the territory to three other nucleotides. 
	The annotations is performed by anno_seqvars.pl with one or more user specified configs.
	When multiple annotation configs are specified, the output will be merged by sequentially
	adding new fields not seen in the previous tables.

	Because the output is typically used in downstream analysis for calculating mutation 
	rates, it is recommended to use the same version of gene table and annotation configs
	as in data processing steps.  

Input/output:
	The output will be a large tabix indexed database file in which each row is one possible
	SNV and its annotations. We will also have one variant table for each individual gene 
	in the out directory useful for downstream parallel processing.

EOF
	exit 1;
}

$opt->validate({ requires => [qw|conf outdir|] });

read_config $opt->get_conf => my %conf;
check_conf(\%conf);

if (!exists $conf{PARAM}{NSPLIT} ||  $conf{PARAM}{NSPLIT} < 2) {
	croak "Incorrect number of split: $conf{PARAM}{NSPLIT}";
}

$conf{PATH}{BINDIR} = shell_quote($Bin);


while (my ($label, $file) = each %{$conf{PATH}}) {
	next if $label eq 'BINDIR';
	croak "Cannot find $label file: $file" unless -f $file;
}

my $rootdir = $opt->get_outdir;

my $wkf = Utils::Workflow->new($rootdir, { engine => $opt->get_engine, force => $opt->get_force });

# Read annotation configs and check input fields
# config should contain standard fields for input table, no extra fields, and no sample level information
read_config  $conf{ANNO}{CONF} => my %anno;
unless($anno{Input}{HG} eq $conf{ANNO}{HG}) {
	die "Incorrect genome assembly build in config: $anno{Input}{HG} <> $conf{ANNO}{HG}";
}

$anno{Input}{Fields} = "Chrom,Position,Ref,Alt";
# Remove sample level fields and other extra-fields
if (defined $anno{Input}{XtraFields}) {
	delete $anno{Input}{XtraFields};
}
if (defined $anno{Sample}) {
	delete $anno{Sample};
}
write_config %anno, "$rootdir/par/Anno.conf";
$conf{ANNO}{CONF} = "$rootdir/par/Anno.conf";

#open my $fgrp, ">$rootdir/par/groups.txt" or die "Cannot write to $rootdir/par/groups.txt";
#print $fgrp "ExomeAnno\t1\t".$conf{PARAM}{NSPLIT}, "\n";
#close $fgrp;


##############################
## Workflow initialization  ##
##############################

# Add jobs
$wkf->add(prep_file(), { name => "PrepFiles", 
		expect => [ map { "wrk/Exome.$_.bed" } 1..$conf{PARAM}{NSPLIT} ] })
	->add(anno_var(),  { name => "AnnoVars", depend => "PrepFiles", nslots => $conf{PARAM}{NSPLIT},
		expect => [ map { ["wrk/Vars.$_.txt", "wrk/Anno.$_.txt"] }  1..$conf{PARAM}{NSPLIT} ],
		callback => \&check_anno  })
	->add(merge_tab(), { name => "MergeTabs",  depend => "AnnoVars", 
		expect => ["ExomeAnno.txt.gz","ExomeAnno.txt.gz.tbi", "ExomeAnno.genes.txt"] })
	->add(gene_tabs(), { name => "GeneTabs",  depend => "AnnoVars", expect => "out/.done",
		callback => \&check_genetabs });

write_config %conf, "$rootdir/par/run.conf" unless $opt->get_dryrun;


################################
## Kickstart workflow engine  ##
################################

$wkf->inst(\%conf);
$wkf->run({ conf => $conf{$wkf->{engine}}, dryrun => $opt->get_dryrun });


############################
## Workflow components    ##
############################


sub prep_file {
	my $script =<< 'EOF';

#cds_bed.pl --dbtab _PATH.GENETAB_ --hgchr --nochr --padding _PARAM.PADDING_ | \
# awk '$1~/[X1-9][0-9]*/"' >  _WRKDIR_/Exome.bed

cds_bed.pl --dbtab _PATH.GENETAB_ --hgchr --nochr --padding _PARAM.PADDING_ --out _WRKDIR_/Exome.bed

#NDIGIT=$(perl -e 'print int(log($ARGV[0])/log(10))+1' _PARAM.NSPLIT_)

# Split then rename
# Note: this require split utility from GNU coreutils v8
#split -n l/_PARAM.NSPLIT_ -a $NDIGIT -d _WRKDIR_/Exome.bed _TMPDIR_/Exome.

#perl -e '$ndigit=$ARGV[0]; for $ii (1..$ARGV[1]){ printf("%d\t%0${ndigit}d\n", $ii, $ii-1) }' $NDIGIT _PARAM.NSPLIT_ | \
#	while read II SUF; do mv _TMPDIR_/Exome.$SUF _WRKDIR_/Exome.$II.bed; done

split_tab.pl --in _WRKDIR_/Exome.bed --out _WRKDIR_/Exome --suffix bed --nsplit _PARAM.NSPLIT_ --noheader

EOF
	return $script;
}

sub anno_var {
	my $script =<< 'EOF';

perl _PATH.BINDIR_/module/enum_snvs_byrng.pl --bed _WRKDIR_/Exome._INDEX_.bed \
	--out _WRKDIR_/Vars._INDEX_.txt --seq _PATH.FASTA_

# We will test if the resulting variant file is empty
NL=$(cat _WRKDIR_/Vars._INDEX_.txt | wc -l | awk '{print $1}')

if [[ $NL == 1 ]]; then
	touch _WRKDIR_/Anno._INDEX_.txt
else
	perl _PATH.BINDIR_/anno_seqvars.pl --conf _ANNO.CONF_ --in _WRKDIR_/Vars._INDEX_.txt \
		--wrkdir _TMPDIR_/Anno._INDEX_ --out _WRKDIR_/Anno._INDEX_.txt
fi

EOF
	return $script;
}

sub check_anno {
	my @exp = @_;
	if (all { -f "$rootdir/$_" } @exp) {
		my $n_in = count_line("$rootdir/$exp[0]");
		my $n_out = count_line("$rootdir/$exp[1]");
		if ($n_in == $n_out || $n_in == 1 && $n_out == 0) {
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

sub merge_tab {
	my $script = <<'EOF';

# Merge all split annotations into a large exome-wide annotation table, indexed by chrom,position
perl _PATH.BINDIR_/../utils/merge_var_splits.pl --input _WRKDIR_/Anno --outdir _OUTDIR_/.. \
	--nsplit _PARAM.NSPLIT_ --group ExomeAnno --bgzip

# Also extract all annotated genes (will be used to check completeness of output)
zcat _OUTDIR_/../ExomeAnno.txt.gz | csvtk cut -t -f GeneID | awk 'NR>1 & $1!="."' | \
	sed -e 's/;/\n/g' | sort -u > _OUTDIR_/../ExomeAnno.genes.txt

EOF
	return $script;
}


sub gene_tabs {
	my $script = << 'EOF';

rm -f _OUTDIR_/.done

# Extract gene-specific annotation tables from the full exome table (in splits)
perl _PATH.BINDIR_/module/extract_anno_pergene.pl --prefix _WRKDIR_/Anno --nsplit _PARAM.NSPLIT_ \
	--outdir _OUTDIR_ --gzip

# Create dummy output file
touch _OUTDIR_/.done

EOF
	return $script;
}

sub check_genetabs {
	my ($dummyexp) = @_;
	my $outdir = "$rootdir/".dirname($dummyexp);
	if (-f "$outdir/../ExomeAnno.genes.txt") {
		my @genes = map { chomp } slurp "$outdir/../ExomeAnno.genes.txt";
		if (all { -s "$outdir/$_.gz" } @genes) {
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


