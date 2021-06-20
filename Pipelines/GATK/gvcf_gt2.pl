#!/usr/bin/env perl

use strict;
use warnings;
use Carp;
use IO::Dir;
use IO::File;
use Data::Dumper;
use FindBin qw|$Bin|;
use File::Copy;
use File::Path qw|make_path|;
use Cwd qw|abs_path|;
use Getopt::Lucid qw|:all|;
use List::MoreUtils qw|all any uniq|;
use Config::Std;
use Utils::Workflow;
use Utils::Hash qw|merge_conf|;


use lib "$Bin/../lib/";
use Shared qw|read_list split_chrs|;
use Wrapper qw|vcf_sampids|;


############################
## Command line interface ##
############################

my @spec =  (
	Param("conf|c")->valid(sub { -r }),
	Param("list|l")->valid(sub { -r }),
	Param("group|g")->valid(sub { -r }),
	Param("rename")->valid(sub { -r }),
	Param("remove")->valid(sub { -r }),
	Param("outdir|out"),
	Param("engine|eng")->valid(sub { $_ eq 'SGE' || $_ eq 'BASH' }),
	Keypair("param|par"),
	Switch("wgs"),
	Switch("help|h"),
	Switch("dryrun|dry"),
	Switch("force")
	);

my $opt = Getopt::Lucid->getopt(\@spec);


if ($opt->get_help) {
	print STDERR <<EOF;
Purpose:
	This is a pipeline script to perform joint genotyping small number of gVCF files.

Usages:
	gvcf_gt2.pl --conf Config --list GVCF.list --out RootDir [--eng Engine]

Options:
	--list: gVCF file list or directory. The list should have 1 or 2 columns per line. Path to gVCF 
			(required) and sample ID. If only one column is available, sample ID will be taken from 
			file basename after removing suffix. This option can be specified multiple times. Only file 
			name suffix '.g.vcf.gz' is considered valid. 
	--wgs:  Under WGS mode, genotyping task will be split into chromsomes, and results from multiple
			splits will be gathered to create final VCFs.
	--group: Sample groups. It should have 2 columns per line, IID and GroupID.
			 IID must be found in gVCF file list. When sample groups are defined, genotyping will be 
			 performed for each group of samples. Samples are allowed to exist in multiple groups. 
			 Samples without GroupID will be skipped. 
	--rename: Sample rename list.
	--remove: A list of bad samples to be removed.

Notes:
	This is the second version of GATK GVCF genotyping piepline that implements GATK's merge gVCF+
	typing approach suitable for genotyping small number of samples from the same family. When a sample
	group file is provided to define sample groups, genotyping is performed for each group of samples,
	and there will be one VCF output for each group. Otherwise, genotyping is performed for each sample
	and the output is one VCF for each sample. 

Dependencies:
	gatk(v4), picard, bcftools, tabix.

EOF
	exit 1;
}

$opt->validate({requires => [qw|conf list outdir|]});

my $rootdir = $opt->get_outdir;

my %conf = merge_conf($opt->get_conf, $opt->get_param); 


############################################
## Input files validation and collection  ##
############################################

my %rename;

my $gvcfs = read_list($opt->get_list, { suffix => "g.vcf.gz", 
						remove => $opt->get_remove, rename => $opt->get_rename });
while(my ($iid, $gvcf) = each %$gvcfs) {
	my @iids = vcf_sampids($gvcf);
	die "gVCF file $gvcf contains multiple samples" unless @iids == 1;
	if ($iids[0] ne $iid) {
		$rename{$iid} = $iids[0];
		#push @rename, [$iids[0], $iid];
	}
}

my (@samps, %group, @grps, @needmrg, @needrn);
if ($opt->get_group) {
	open my $fin, $opt->get_group or die "Cannot open group file";
	while(<$fin>) {
		next if /^\s*$/;
		my ($iid, $gid) = split;
		unless (defined $gvcfs->{$iid}) {
			warn "Cannot find gVCF file for $iid in group $gid";
		}
		else {
			push @{$group{$gid}} => $iid;
		}
	}
	@grps = sort keys %group;
	@samps = sort map { @$_ } values %group;
	foreach my $gid (@grps) {
		if (@{$group{$gid}} > 1) {
			push @needmrg, $gid;
		}
		if (any { defined $rename{$_} } @{$group{$gid}}) {
			push @needrn, $gid;
		}
	}
}
else {
	@samps = sort keys %$gvcfs;
	@grps = @samps;
	@group{@samps} = map { [$_] } @samps;
	foreach my $iid (@samps) {
		if (defined $rename{$iid}) {
			push @needrn, $iid;
		}
	}
}

#################################
##  Workflow initialization    ##
##  Working directory setup    ##
#################################

my $wkf = Utils::Workflow->new($rootdir,
	{ engine => $opt->get_engine, force => $opt->get_force, strict_var => 1 });


if ($opt->get_wgs) {
	foreach my $gid (@grps) {
		make_path "$rootdir/wrk/$gid";
	}
	my @chrgrp = split_chrs($conf{PATH}{SEQDICT});
	for(my $ii = 0; $ii < @chrgrp; $ii ++) {
		my $jj = $ii + 1;
		open my $fout, ">$rootdir/par/SEQGRP.$jj" or die "Cannot write to SEQGRP";
		print $fout join(",", @{$chrgrp[$ii]}), "\n";
	}
	$conf{PARAM}{NSPLIT} = @chrgrp;
}
else {
	$conf{PARAM}{NSPLIT} = 1;	
}

for(my $ii = 0; $ii < @grps; $ii ++) {
	my $jj = $ii + 1;
	open my $fout, ">$rootdir/par/GRPID4GT.$jj" or die "Cannot write to GRPID4GT";
	if (@{$group{$grps[$ii]}} == 1) {
		my $iid = $group{$grps[$ii]}[0];
		print $fout $grps[$ii], "\t", $gvcfs->{$iid}, "\n";
	}
	else {
		print $fout $grps[$ii], "\t$rootdir/wrk/$grps[$ii].g.vcf.gz\n";
	}
}

if (@needmrg) {
	for(my $ii = 0; $ii < @needmrg; $ii ++) {
		my $jj = $ii + 1;
		open my $fout, ">$rootdir/par/GRPID4COMB.$jj" or die "Cannot write to GRPID4COMB";
		print $fout $needmrg[$ii], "\t", 
			join(',', map { $gvcfs->{$_} } @{$group{$needmrg[$ii]}} ), "\n";
	}
}

if (@needrn) {
	#open my $fout, ">$rootdir/par/RENAME" or die "Cannot write to RENAME list";
	#for(my $ii = 0; $ii < @rename; $ii ++) {
	#	print $fout $rename[$ii][0], "\t", $rename[$ii][1], "\n";
	#}
	for(my $ii = 0; $ii < @needrn; $ii ++) {
		my $jj = $ii + 1;
		open my $fout, ">$rootdir/par/GRPID4RENAME.$jj" or die "Cannot write to GRPID4RENAME.$jj";
		print $fout $needrn[$ii], "\n";
		open my $frn, ">$rootdir/par/RENAME.$jj" or die "Cannot write to RENAME.$jj";
		foreach my $iid (@{$group{$needrn[$ii]}}) {
			if (defined $rename{$iid}) {
				print $frn $rename{$iid}, "\t", $iid, "\n";
			}
		}
	}
}

write_config %conf, "$rootdir/par/run.conf" unless $opt->get_dryrun;


##############################
## Workflow initialization  ##
##############################

my %deparg;
if (@needmrg) {
	$wkf->add(merge_gvcfs(), { name => "MergeGVCFs", 
		expect => [ map { ["wrk/$_.g.vcf.gz", "wrk/$_.g.vcf.gz.tbi"] } @needmrg ],
		nslots => scalar(@needmrg)  });
	$deparg{depend} = "MergeGVCFs";
}

if ($opt->get_wgs) {
	my $unit = $conf{PARAM}{NSPLIT};
	$wkf->add(gt_gvcfs_split(), { name => "GtGVCFs", %deparg,
			expect => get_exp_split($unit), 
			nslots => scalar(@grps)*$unit, step => 1 })
		->add(vcf_merge(), { name => "VcfMerge", depend => "GtGVCFs", deparray => 1,
			expect => [ map { ["out/$_.vcf.gz", "out/$_.vcf.gz.tbi"] } @grps ],
			nslots => scalar(@grps)*$unit, step => $unit });
	$deparg{depend} = "VcfMerge";
}
else {
	$wkf->add(gt_gvcfs(), { name => "GtGVCFs", %deparg,
			expect => [ map { ["out/$_.vcf.gz", "out/$_.vcf.gz.tbi"] } @grps ],
			nslots => scalar(@grps) });
	$deparg{depend} = "GtGVCFs";
}

if (@needrn) {
	$wkf->add(rename_vcf(), { name => "Rename", %deparg,
			expect => [ map { ["out/$_.vcf.gz.bak", "out/$_.vcf.gz.tbi.bak"] } @needrn ],
			nslots => scalar(@needrn) });
}

$wkf->inst(\%conf);
$wkf->run({ conf => $conf{$wkf->{engine}}, dryrun => $opt->get_dryrun  });


############################
## Workflow components    ##
############################

# If a group contain more than one gVCFs, we need to combine gVCFs
# otherwise, they can be genotyped directly
sub merge_gvcfs {
	my $script  = <<'EOF';

read GRPID GVCFs < _PARDIR_/GRPID4COMB._INDEX_

NGVCF=$(echo $GVCFs | awk -F, '{print NF}')

if [[ $NGVCF == 1 ]]; then
	echo "Merge GVCFs need more than one GVCF as input"
	exit 1
fi
 
INPUTS=$(echo $GVCFs | sed -e 's/,/ -V /g')
gatk --java-options "_GATK.JAVAOPT_" \
	CombineGVCFs \
	-R _PATH.FASTA_ \
	-V $INPUTS \
	-O _WRKDIR_/$GRPID.g.vcf.gz 

EOF
}

# Genotype merged gVCFs.
sub gt_gvcfs {
	my $script = <<'EOF';

read GRPID GVCF < _PARDIR_/GRPID4GT._INDEX_

gatk --java-options "_GATK.JAVAOPT_" \
	GenotypeGVCFs \
	-R _PATH.FASTA_ \
	-V $GVCF \
	-O _OUTDIR_/$GRPID.vcf.gz \
	-L _PATH.TARGETS_ \
	_GATK.OPTION[ \
	]_

EOF
}

# For WGS we will genotype by chr splits
sub gt_gvcfs_split {
	my $script = <<'EOF';

II=$(( (_INDEX_-1)/_PARAM.NSPLIT_+1 ))
SPLIT=$(( (_INDEX_-1)%_PARAM.NSPLIT_+1 ))

read GRPID GVCF < _PARDIR_/GRPID4GT.$II

INTERVAL=`cat _PARDIR_/SEQGRP.$SPLIT | sed -e 's/,/ -L /g'`

gatk --java-options "_GATK.JAVAOPT_" \
	GenotypeGVCFs \
	-R _PATH.FASTA_ \
	-V $GVCF \
	-O _OUTDIR_/$GRPID/split.$SPLIT.vcf.gz \
	-L $INTERVAL \
	_GATK.OPTION[ \
	]_

EOF
}

sub get_exp_split {
	my ($unit) = @_;
	my @allexps;
	foreach my $iid (@grps) {
		for(my $jj = 1; $jj <= $unit; $jj ++) {
			my @exp;
			foreach my $suffix (qw|vcf.gz vcf.gz.tbi|) {
				push @exp, "wrk/$iid/split.$jj.$suffix";
			}
			push @allexps, \@exp;
		}
	}
	return @allexps;
}

# Merge split VCFs
sub vcf_merge {
	my $script = <<'EOF';

II=$(( (_INDEX_-1)/_PARAM.NSPLIT_+1 ))

read GRPID GVCF < _PARDIR_/GRPID4GT.$II

INPUT=`perl -e 'print join(" ", map {"INPUT=_WRKDIR_/$ARGV[1]/split.$_.vcf.gz" } 1..$ARGV[0])' _PARAM.NSPLIT_ $GRPID`

picard _PICARD.JAVAOPT_ MergeVcfs \
	$INPUT \
	OUTPUT=_OUTDIR_/$GRPID.vcf.gz

EOF
}


# Rename if necessary
sub rename_vcf {
	my $script = <<'EOF';

read GRPID < _PARDIR_/GRPID4RENAME._INDEX_

mv _OUTDIR_/$GRPID.vcf.gz _OUTDIR_/$GRPID.vcf.gz.bak

mv _OUTDIR_/$GRPID.vcf.gz.tbi _OUTDIR_/$GRPID.vcf.gz.tbi.bak

bcftools reheader -s _PARDIR_/RENAME._INDEX_ -o _OUTDIR_/$GRPID.vcf.gz _OUTDIR_/$GRPID.vcf.gz.bak

tabix -f -p vcf _OUTDIR_/$GRPID.vcf.gz

EOF
}
