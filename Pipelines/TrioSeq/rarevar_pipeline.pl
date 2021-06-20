#!/usr/bin/env perl

use strict;
use warnings;
use Carp;
use Tie::IxHash;
use FindBin qw|$Bin|;
use POSIX qw|ceil|;
use Digest::MD5;
use Getopt::Lucid qw|:all|;
use IO::Prompt;
use Data::Dumper;
use Config::Std;
use File::Basename;
use Cwd qw|abs_path|;
use List::MoreUtils qw|all any uniq|;
use Digest::MD5 qw|md5_hex|;
use File::Path qw|make_path|;
use String::ShellQuote;
use Utils::Workflow;
use Utils::Hash qw|merge_conf|;
use Utils::File qw|count_line|;
use Genome::Ranges::IntSet;

use lib "$Bin/../lib";
use Shared qw|parse_fstr|;

#############################
## Command line interface  ##
#############################

my @spec =  (
	Param("conf|c")->valid(sub { -r }),
	Param("list|l")->valid(sub { -r }),
	Param("ped")->valid(sub { -r }),
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
	This is a pipeline script to run variant filtering and annotation in parallel.

Usage: 
	rarevar_pipeline.pl --conf Config --list VCFList [--ped Cohort.ped] --out OutDir 

Options:
	--list: VCF file list (required). It should contain 3~4 columns described below.
	--ped: A standard 6-col PED file describing pedigree relationship of samples in VCF file.
			It overrides the PED file in config.

Notes:
	This is a wrapper of variant filtering from VCF and annotation on the filtered variants. 
	It currently supports one of dnv_finder, var_filter and site_filter for variant filtering.
	It uses anno_seqvars for annotating filtered variants. The config files for filtering should be 
	specified in the config. Because the filter is applied to process all VCF files in the list, 
	users should make sure the parameters in the config file can be applied to all VCF files in the list. 
	All fields in the filtered variant table will automatically be included in annotation output.
	The annotation step can be optional.

	The input list is a list of VCF and region combinations. For large VCF files, it is preferable
	to split the tasks into different genomic regions.   
		List file format: VCF REGION LABEL [NSPLIT]
	REGION can be either a chromosome name, a region specification ("Chr:Start-End"), or a BED file.
	The fourth (optional) column NSPLIT is the number of slots for parallelization, if NSPLIT is 
	provided then REGION must be a BED file. LABEL is used to group differnet VCF-region combinations 
	and used as the file name prefix in output. Different LABEL can be assigned to different sequencing 
	batches or different callers. Typically there should be no sample duplication in the same group, 
	because when merging result, the pipeline does not do dedup. 
	
Output:
	Filtering and annotation will be performed in parallel for each split of VCF-region combination.
	The output from each split will be merged by group. Merged variants of each group will appear
	in the out directory named after LABEL.

EOF
	exit 1;
}


$opt->validate({requires => [qw|list conf outdir|]});


my %conf = merge_conf($opt->get_conf, $opt->get_param); 


unless ($conf{PATH}{PROG} eq 'dnv_finder.pl' || $conf{PATH}{PROG} eq 'var_filter.pl' ||
		$conf{PATH}{PROG} eq 'site_filter.pl' ) {
	die "Currently only support pipelining dnv_finder.pl or var_filter.pl or site_filter.pl";
}
$conf{PATH}{BINDIR} = shell_quote($Bin);

if ($opt->get_ped) {
	$conf{Pedigree}{File} = $opt->get_ped;
}
unless(defined $conf{Pedigree}{File} && -f $conf{Pedigree}{File}) {
	unless($conf{PATH}{PROG} eq 'site_filter.pl') {
		die "Cannot find pedigree file: $conf{Pedigree}{File}";
	}
}

#############################
## Input files preparation ##
#############################

my $rootdir = $opt->get_outdir;

my $wkf = Utils::Workflow->new($rootdir,
	{ engine => $opt->get_engine, force => $opt->get_force });

my (@vcfbeds, %groups);
{
	tie my %inbygrp, 'Tie::IxHash';
	open my $fin, $opt->get_list or die "Cannot open input list";
	while(<$fin>) {
		my @line = split;
		unless(@line == 3 || @line == 4) {
			die "Incorrect number of columns, should give: VCF REGION LABEL [NSPLIT]";
		}
		my ($vcf, $region, $label, $nsplit) = @line;
		unless (defined $nsplit) {
			push @{$inbygrp{$label}}, [$vcf, $region];
		}
		else {
			# Split the region for the group
			# To allow a group contains multiple region files that need split,
			# the split region files will be named after group name and MD5 sum
			# of input region file name.
			unless(-f $region) {
				die "Region $region for split is not a file!";
			}
			make_path "$rootdir/tmp/$label";
			my $digest = md5_hex($region);
			system("split_bed.pl --nsplit $nsplit --bed $region --out $rootdir/tmp/$label/$digest");
			# Validae that split bed was run succesfully,
			my $nlorig = count_line($region);
			my $nlsplit = 0;
			for(my $ii = 1; $ii <= $nsplit; $ii ++) {
				$nlsplit += count_line("$rootdir/tmp/$label/$digest.$ii.bed");
				push @{$inbygrp{$label}}, [$vcf, "$rootdir/tmp/$label/$digest.$ii.bed"];
			}
			unless($nlorig == $nlsplit) {
				die "Split files does not have the same total number of lines as original";
			}
		}
	}
	my $ct = 0;
	my @labels = keys %inbygrp;
	foreach my $label (@labels) {
		foreach my $dat (@{$inbygrp{$label}}) {
			push @vcfbeds, $dat;
			$ct ++;
			push @{$groups{$label}}, $ct;
		}
	}
	# Write DATA files
	for(my $ii = 0; $ii < @vcfbeds; $ii ++) {
		my $jj = $ii + 1;
		open my $fout, ">$rootdir/par/DATA.$jj" or die "Cannot write to DATA.$jj";
		print $fout join("\t", @{$vcfbeds[$ii]}), "\n";
	}
	# Write GROUP files
	open my $fgrp, ">$rootdir/par/GROUPS" or die "Cannot write to GROUPS";
	foreach my $label (@labels) {
		my @index = @{$groups{$label}};
		print $fgrp join("\t", $label, $index[0], $index[$#index]), "\n";
	}
}


# Customize and write the filter and anno config file
{
	read_config $conf{PATH}{FILTER} => my %filter;
	$filter{Pedigree} = $conf{Pedigree} if defined $conf{Pedigree};
	write_config %filter, "$rootdir/par/Filter.conf";
	$conf{PATH}{FILTER} = "$rootdir/par/Filter.conf";

	my @fields;
	my %stfd = ('Chrom' => 1, 'Position' => 1, 'Ref' => 1, 'Alt' => 1);
	my $site_f = parse_fstr($filter{Output}{Site}, 1);
	push @fields, grep { not defined $stfd{$_} } values %$site_f;
	
	my $geno_f = parse_fstr($filter{Output}{Geno}, 1);
	my $siteonly;
	if ($conf{PATH}{PROG} =~ /^dnv_finder/) {
		foreach my $samp (qw|Offspring Father Mother|) {
			push @fields, map { "$_.$samp" } values %$geno_f;
		}
	}
	elsif ($conf{PATH}{PROG} =~ /^var_filter/) {
		if (defined $conf{Pedigree}{Include}) {
			push @fields, values %$geno_f;
			push @fields, qw|FamMembers Relations Phenotypes|;
			push @fields, map { "$_.FamMembers" } values %$geno_f;
		}
		else {
			$siteonly = 1;
		}
	}
	elsif ($conf{PATH}{PROG} =~ /^site_filter/) {
		$siteonly = 1;
	}
	else {
		die "Cannot determine genotype fields: $conf{PATH}{PROG}";
	}
	
	if (exists $conf{PATH}{ANNO}) {
		read_config $conf{PATH}{ANNO} => my %anno;
		$anno{Input}{Fields} = "Chrom,Position,Ref,Alt";
		$anno{Input}{HG} = $filter{Genome}{Build};
		$anno{Input}{XtraFields} = join(',', @fields);
		unless ($siteonly) {
			if ($conf{PATH}{PROG} =~ /^dnv_finder/) {
				$anno{Sample}{IDField} = 'IID,FamID,Father,Mother,Gender,Pheno';
			}
			elsif ($conf{PATH}{PROG} =~ /^var_filter/) {
				$anno{Sample}{IDField} = 'IID,FamID,Gender,Pheno';
			}
		}
		else {
			delete $anno{Sample}{IDField};
		}
		write_config %anno, "$rootdir/par/Anno.conf";
		$conf{PATH}{ANNO} = "$rootdir/par/Anno.conf", 
	}
}
write_config %conf, "$rootdir/par/run.conf" unless $opt->get_dryrun;


#################################
##  Workflow initialization &  ##
##  Working directory setup    ##
#################################

$wkf->add(var_filt(), { name => 'VarFilter', nslots => scalar(@vcfbeds), 
				expect => [ map { "wrk/filter.$_.txt" } 1..scalar(@vcfbeds) ] });
my %deparg;
if (exists $conf{PATH}{ANNO}) {
	$wkf->add(var_anno(), { name => 'VarAnno', nslots => scalar(@vcfbeds), 
				depend => 'VarFilter', deparray => 1,
				expect => [ map { "wrk/anno.$_.txt" } 1..scalar(@vcfbeds) ] });
	$deparg{depend} = 'VarAnno';
}
else {
	$deparg{depend} = 'VarFilter';
}
$wkf->add(var_merge(), { name => 'VarMerge',  %deparg,
			expect => [ map { "out/$_.txt" } sort keys %groups ] });



################################
## Kickstart workflow engine  ##
################################

$wkf->inst(\%conf);
$wkf->run({ conf => $conf{$wkf->{engine}}, dryrun => $opt->get_dryrun });


############################
## Workflow components    ##
############################

sub var_filt {
	my $script = <<'EOF';

read VCF REGION < _PARDIR_/DATA._INDEX_

perl _PATH.BINDIR_/_PATH.PROG_ --conf _PATH.FILTER_ --vcf $VCF --bed $REGION \
	--out _WRKDIR_/filter._INDEX_.txt

EOF
}

sub var_anno {
	my $script = <<'EOF';

NL=$(cat _WRKDIR_/filter._INDEX_.txt | wc -l | awk '{print $1}')

if [[ $NL == 1 ]]; then
	touch  _WRKDIR_/anno._INDEX_.txt
else
	perl _PATH.BINDIR_/../VarAnno/anno_seqvars.pl --in _WRKDIR_/filter._INDEX_.txt \
		--out _WRKDIR_/anno._INDEX_.txt --conf _PATH.ANNO_ --wrkdir _TMPDIR_/anno._INDEX_
fi

EOF
}


sub var_merge {
	my $script;
	if (exists $conf{PATH}{ANNO}) {
		$script = <<'EOF';

perl _PATH.BINDIR_/../utils/merge_var_splits.pl --input _WRKDIR_/anno --outdir _OUTDIR_ --group _PARDIR_/GROUPS

EOF
	}
	else {
		$script = <<'EOF';

perl _PATH.BINDIR_/../utils/merge_var_splits.pl --input _WRKDIR_/filter --outdir _OUTDIR_ --group _PARDIR_/GROUPS

EOF
	}	
}
