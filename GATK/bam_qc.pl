#!/usr/bin/env perl

use strict;
use warnings;
use Carp;
use IO::File;
use Data::Dumper;
use FindBin qw|$Bin|;
use File::Copy;
use File::Path qw|make_path|;
use Cwd qw|abs_path|;
use Getopt::Lucid qw|:all|;
use Config::Std;
use Hash::Util qw|lock_hash_recurse|;
use String::ShellQuote;
use Utils::Workflow;
use Utils::Hash qw|merge_conf|;

use lib "$Bin/../lib/";
use Shared qw|read_list|;
use Wrapper qw|bam_sampids|;

############################
## Command line interface ##
############################

my @spec =  (
	Param("conf|c")->valid(sub { -r }),
	Param("list|l")->valid(sub { -r }),
	Param("bed|b")->valid(sub { -r }),
	Param("rename")->valid(sub { -r }),
	Param("outdir|out"),
	Param("engine|eng")->valid(sub { $_ eq 'SGE' || $_ eq 'BASH' }),
	Keypair("param|par"),
	Switch("wgs"),
	Switch("novb"),
	Switch("help|h"),
	Switch("dryrun|dry"),
	Switch("force")
	);

my $opt = Getopt::Lucid->getopt(\@spec);

if ($opt->get_help) {
	print STDERR <<EOF;
Purpose:
	This is a pipeline script to calculate BAM/CRAM level QC statistics

Usage:
	bam_qc.pl --conf Conf --list BAMList --outdir RootDir

Options:
	--list: (Required) BAM/CRAM file list should have 1 or 2 columns per line: path to the BAM/CRAM file 
			and sample ID. If sample ID is not provided, it will be taken from the basename of file 
			after stripping off suffix. Only one file type can be found in a input list.
    --wgs:  The default workflow is suitable for WES. For WGS data, this option should be switched on.
    		The main difference in WGS and WES is the way depth of coverage is calculated.
    --novb: Skip running verifyBamID. To run verifyBamID (v1) on CRAMs require conversion to BAM at
    		around known SNP sites. It is very time consuming. Turn on this option to skip verifyBamID.  
    --rename: Sample names can be changed by providing a rename list.
	--bed:  Targeted region of interest in BED format. It overrides the target file in config.

Output:
	The following tools are run for each individual BAM/CRAM file in parallel. The intermediate results 
	are stored in individual-specific subdirectory under wrk.

	Picard
		CollectAlignmentSummaryMetrics: alignment_summary_metrics, 
		CollectGcBiasMetrics: BamStats.gc_bias.pdf, BamStats.gc_bias.summary_metrics, BamStats.gc_bias.detail_metrics 
		CollectInsertSizeMetrics: BamStats.insert_size_histogram.pdf, BamStats.insert_size_metric 
		CollectOxoGMetrics(support BAM file only): BamStats.cpcg_metrics
		CollectHsMetrics (for WES):  BamStats.hs_metrics
		CollectWgsMetrics (for WGS): BamStats.wgs_metrics
		 
	MosDepth
		Depth of coverage per target (for WES) or sliding window (for WGS): DoC.regions.bed.gz 
			DoC per target can be used by xHMM CNV calling.
		Callable regions: DoC.quantized.bed.gz 
			Callable region bed files can be used in callable_merge.pl 
		Average and normalized depthes per chromosome: chrom_NormDP.txt, it will be used to infer sex.

	VerifyBamID
		Estimated level cross-sample contamination (FREEMIX) can be found in VBDI.selfSM.
		Note: Currently it runs verifyBamID v1 under sequence only mode. It will be upgraded to v2 in the future.

	The key summary statistics will be gather into a final table "bam_stats.csv" in out directory, columns
	are described in "bam_stats.README.txt". Samples that does not meet depth of coverage criteria or have high 
	level of contmination (defined by config) will be listed in "badsamps.csv".

	We also generated following plots:
		* samples_DepthOfCov.png -- Scatterplot of mean depth of coverage vs percentage of target at 10x/20x
		* samples_chrXYNormDp_byPredSex -- Scatterplot of normalized chrX and chrY depth colored by inferred sample sex.
				Sex is inferred from clusters of chrX and chrY NPDs using k-means. 
				The results can be used for sex check with phenotype records and identify sex chromosome abnormalities.
		* chrom_NormDP_CDF -- Cumulative distribution of normalized depth of coverage for each chromosome. 
				Sample level DoC per chromosome can be found in chrom_depth.csv. 
				The results can be used to identify likely trisomy. 

Dependencies:
	picard, gatk(v4), mosdepth, bcftools, verifyBamID, csvtk and R-ggplot2.
	Note: gatk(v4) is only used to convert CRAM to BAM around known SNPs in the target region 
			(which are extracted by bcftools) for running verifyBamID. 

EOF
	exit 1;
}

$opt->validate({requires => [qw|conf list outdir|]});

#unless($^O =~ /linux/) {
#	print STDERR "This is script can only run on linux\n";
#	exit 1;
#}

# Read in config parameters
my %conf = merge_conf($opt->get_conf, $opt->get_param); 

if (defined $opt->get_bed) {
	$conf{PATH}{TARGETS} = abs_path($opt->get_bed);
}
$conf{PATH}{MODULE} = shell_quote("$Bin/module");

lock_hash_recurse(%conf);

############################
## Input files validation ##
############################

my %opt = exists $conf{AWS} && exists $conf{AWS}{STAGEIN} ? (notest => 1) : ();
my ($bams, $filetype) = read_list($opt->get_list, 
						{ suffix => ['bam', 'cram'], rename => $opt->get_rename, %opt });
my %bams = %$bams;
my $nbams = keys %bams;

exit 1 unless $nbams;

##############################
## Workflow initialization  ##
## Config file parsing      ##
##############################

my $rootdir = $opt->get_outdir;

my $wkf = Utils::Workflow->new($rootdir, 
	{ engine => $opt->get_engine, force => $opt->get_force, strict_var => 1});

# Add jobs to workflow
# NOTE: need module level script to collect sample-level summary statistics
$wkf->add(prep_files(), { name => "PrepFiles", expect => [prep_files_out()] });

if (exists $conf{AWS} && exists $conf{AWS}{STAGEIN}) {
	$wkf->add(bam_qc(), { name => "BamQC", nslots => $nbams, depend => "PrepFiles",
			expect => [ map { [ bam_qc_files($_) ] } sort keys %bams ]  });	
}
else {
	$wkf->add(bam_stats(), { name => "BAMStats", nslots => $nbams,
			expect => [ map { [ bam_stats_files($_) ] } sort keys %bams ] });
	# Note: currently there is no support for CRAM file for OxoG metrics
	if ($filetype eq 'bam') {
		$wkf->add(bam_oxog(), { name => "BAMOxoG", nslots => $nbams,
			expect => [ map { [ bam_oxog_files($_) ] } sort keys %bams ] });
	}
	$wkf->add(depth_cov(), { name => "DepthOfCov", nslots => $nbams, depend => "PrepFiles",
			expect => [ map { [ depth_cov_files($_) ] } sort keys %bams] })
		->add(mos_depth(), { name => "MosDepth", nslots => $nbams, depend => "PrepFiles",
			expect => [ map { [ mos_depth_files($_) ] } sort keys %bams ] });
	unless ($opt->get_novb) {
		$wkf->add(contam_test($filetype), { name => "VerifyBamID", nslots => $nbams, depend => "PrepFiles",
				expect => [ map { [ contam_test_files($_) ] } sort keys %bams ] });		
	}
}
if (exists $conf{AWS} && exists $conf{AWS}{STAGEOUT}) {
	$wkf->add(stage_out_indiv(), { name => "StageOutIndiv", nslots => $nbams,
			depend => exists $conf{AWS}{STAGEIN} ? "BamQC" :
						($opt->get_novb ? [qw|BAMStats DepthOfCov|] : [qw|BAMStats DepthOfCov VerifyBamID|]),
			expect => [ map { "wrk/$_/stagedout" } sort keys %bams ]  })
}

$wkf->add(sum_stats(), { name => "SumStats", 
			depend => exists $conf{AWS} && exists $conf{AWS}{STAGEIN} ? "BamQC" : 
						($opt->get_novb ? [qw|BAMStats DepthOfCov|] : [qw|BAMStats DepthOfCov VerifyBamID|]),
			expect => ["out/bam_stats.csv", "out/chrom_depth.csv"] })
	->add(plot_stats(), { name => "PlotStats", depend => "SumStats", interp => "Rscript",
			expect => ["out/samples_DepthOfCov.png", "out/samples_chrXYNormDp_byPredSex.png", "out/chrom_NormDP_CDF.pdf"] } );
if (exists $conf{AWS} && exists $conf{AWS}{STAGEOUT}) {
	$wkf->add(stage_out(), { name => "StageOut", depend => "PlotStats", expect => "out/stagedout"  });
}


# Write configs to par directory.
write_config %conf, "$rootdir/par/run.conf" unless $opt->get_dryrun;

# For for existance of external files
for my $FILE (qw|SEQDICT FASTA TARGETS|) {
	croak "Cannot find $FILE!" unless -f $conf{PATH}{$FILE};
}
unless($opt->get_novb) {
	croak "Cannot find K1GSNP!" unless -f $conf{PATH}{K1GSNP};
}


##############################
## Working directory setup  ##
##############################

# Write IID2BAM parameter files to par dir.
my $ii = 1;
foreach my $iid (sort keys %bams) {
	make_path "$rootdir/wrk/$iid";
	open my $fout, ">$rootdir/par/IID2BAM.$ii" or die "Cannot write to par/IID2BAM.$ii";
	print $fout $iid, "\t", $bams{$iid}, "\n";
	$ii ++;
}


################################
## Kickstart workflow engine  ##
################################

$wkf->inst(\%conf);
$wkf->run({ conf => $conf{$wkf->{engine}}, dryrun => $opt->get_dryrun });


############################
## Workflow components    ##
############################

# Prepare intervals and K1G SNPs within intervals
sub prep_files {
	my $script =<<'EOF';

picard _RSRC.PICARDOPT_ BedToIntervalList I=_PATH.TARGETS_ \
	O=_WRKDIR_/targets.interval_list SD=_PATH.SEQDICT_

# For HsMetrics, the padded regions will be used as baits
picard _RSRC.PICARDOPT_ IntervalListTools I=_WRKDIR_/targets.interval_list \
	PADDING=_PARAM.BAITPAD_ ACTION=UNION O=_WRKDIR_/targets_padded.interval_list

# keeping only the first three columns for BED file
# Note: bed file is zero-based half-open; interval_list file in 1-based, both end closed
#grep -v '^@' _WRKDIR_/targets_padded.interval_list | \
#	awk 'BEGIN{OFS="\t"}{print $1, $2-1, $3}' > _WRKDIR_/targets_padded.bed

EOF
	unless($opt->get_novb) {
		$script .= <<'EOF';

# Common SNPs (of 1KG) within target regions, used by verifyBamID
bcftools view -m 2 -M 2 -v snps --exclude 'AF<0.01 || AF>0.99' \
	-o _WRKDIR_/targets.snps.vcf -T _PATH.TARGETS_ _PATH.K1GSNP_

EOF
	}
	return $script;
}

sub prep_files_out {
	if ($opt->get_novb) {
		return ("wrk/targets.interval_list", "wrk/targets_padded.interval_list");
	}
	else {
		return ("wrk/targets.interval_list", "wrk/targets_padded.interval_list", "wrk/targets.snps.vcf");
	}
}

# This is used when stage-in is enabled
sub bam_qc {
	my $stagein = << 'EOF';
read IID BAMFILE < _PARDIR_/IID2BAM._INDEX_

TMPDIR=$(stagein.pl --input $BAMFILE --output _TMPDIR_/BAM._INDEX_ --profile _AWS.PROFILEIN_ --tmpdir _AWS.STAGEIN_ --all)

trap "echo Clean up $TMPDIR; rm -fR $TMPDIR; exit 1" EXIT SIGINT SIGTERM

read BAMFILE < _TMPDIR_/BAM._INDEX_

EOF
	my $bam_qc = bam_stats();
	if ($filetype eq 'bam') {
		$bam_qc .= bam_oxog();
	}
	$bam_qc .= depth_cov();
	$bam_qc .= mos_depth();
	unless ($opt->get_novb) {
		$bam_qc .= contam_test($filetype);
	}
	$bam_qc =~ s/read IID BAMFILE < _PARDIR_\/IID2BAM._INDEX_//g;
	return $stagein.$bam_qc;
}

sub bam_qc_files {
	my @outfiles = bam_stats_files(@_);
	if ($filetype eq 'bam') {
		push @outfiles => bam_oxog_files(@_);
	}
	push @outfiles => depth_cov_files(@_);
	push @outfiles => mos_depth_files(@_);
	unless ($opt->get_novb) {
		push @outfiles => contam_test_files(@_);
	}
	return @outfiles;
}

sub stage_out_indiv {
	my $indivfiles = join(',', bam_qc_files('$IID'));
	my $script = <<'EOF';
read IID BAMFILE < _PARDIR_/IID2BAM._INDEX_

stageout.pl --files INDIVFILES --profile _AWS.PROFILEOUT_ --indir _WRKDIR_/.. --outdir _AWS.STAGEOUT_/wrk/$IID 

touch _WRKDIR_/$IID/stagedout

EOF
	$script =~ s/INDIVFILES/$indivfiles/;
	return $script;
}

sub stage_out {
	my $script = <<'EOF';

stageout.pl --files bam_stats.csv,chrom_depth.csv,samples_DepthOfCov.png,samples_chrXYNormDp_byPredSex.png,chrom_NormDP_CDF.pdf \
	 --profile _AWS.PROFILEOUT_ --indir _OUTDIR_ --outdir _AWS.STAGEOUT_/out

touch _OUTDIR_/stagedout

EOF
	return $script;
}


sub bam_stats {
	my $script = <<'EOF';

read IID BAMFILE < _PARDIR_/IID2BAM._INDEX_

#picard _RSRC.PICARDOPT_ CollectMultipleMetrics I=$BAMFILE R=_PATH.FASTA_ \
#	O=_WRKDIR_/$IID/BamStats PROGRAM=CollectGcBiasMetrics

picard _RSRC.PICARDOPT_ CollectAlignmentSummaryMetrics I=$BAMFILE R=_PATH.FASTA_ \
	O=_WRKDIR_/$IID/BamStats.alignment_summary_metrics

picard _RSRC.PICARDOPT_ CollectGcBiasMetrics I=$BAMFILE R=_PATH.FASTA_ \
	O=_WRKDIR_/$IID/BamStats.gc_bias.detail_metrics \
	S=_WRKDIR_/$IID/BamStats.gc_bias.summary_metrics \
	CHART=_WRKDIR_/$IID/BamStats.gc_bias.pdf

picard _RSRC.PICARDOPT_ CollectInsertSizeMetrics I=$BAMFILE R=_PATH.FASTA_ \
	O=_WRKDIR_/$IID/BamStats.insert_size_metrics \
	H=_WRKDIR_/$IID/BamStats.insert_size_histogram.pdf
 
grep -v '^#' _WRKDIR_/$IID/BamStats.alignment_summary_metrics | head -n6 | csvtk transpose -t | \
	csvtk pretty -t > _WRKDIR_/$IID/BamStats.alignment_summary_metrics.txt 

grep -v '^#' _WRKDIR_/$IID/BamStats.gc_bias.summary_metrics | csvtk transpose -t | \
	csvtk pretty -t > _WRKDIR_/$IID/BamStats.gc_bias.summary_metrics.txt 

grep -v '^#' _WRKDIR_/$IID/BamStats.insert_size_metrics | head -n4 | csvtk transpose -t | \
	csvtk pretty -t > _WRKDIR_/$IID/BamStats.insert_size_metrics.txt 

EOF
	
	return $script;
}

sub bam_stats_files {
	my ($iid) = @_;
	my @outfiles = map { "wrk/$iid/BamStats.$_" }
		qw|insert_size_histogram.pdf insert_size_metrics insert_size_metrics.txt
			gc_bias.pdf gc_bias.summary_metrics gc_bias.summary_metrics.txt 
			gc_bias.detail_metrics alignment_summary_metrics alignment_summary_metrics.txt|;
	return @outfiles;
}

sub bam_oxog {
	my $script .= <<'EOF';

read IID BAMFILE < _PARDIR_/IID2BAM._INDEX_

picard _RSRC.PICARDOPT_ CollectOxoGMetrics I=$BAMFILE R=_PATH.FASTA_ \
	O=_WRKDIR_/$IID/BamStats.cpcg_metrics

grep -v '^#' _WRKDIR_/$IID/BamStats.cpcg_metrics | csvtk transpose -t | \
	csvtk pretty -t > _WRKDIR_/$IID/BamStats.cpcg_metrics.txt

EOF

	return $script;
}

sub bam_oxog_files {
	my ($iid) = @_;
	my @outfiles = map { "wrk/$iid/BamStats.$_" } qw|cpcg_metrics cpcg_metrics.txt|;
	return @outfiles;
}

sub depth_cov {
	my $script;
	unless($opt->get_wgs) {
		$script = <<'EOF'

read IID BAMFILE < _PARDIR_/IID2BAM._INDEX_

picard _RSRC.PICARDOPT_ CollectHsMetrics \
	I=$BAMFILE R=_PATH.FASTA_ TI=_WRKDIR_/targets.interval_list \
	BAIT_INTERVALS=_WRKDIR_/targets_padded.interval_list \
	O=_WRKDIR_/$IID/BamStats.hs_metrics \
	COVERAGE_CAP=_PARAM.COVCAP_ \
	MINIMUM_MAPPING_QUALITY=_PARAM.MINMQ_

grep -v '^#' _WRKDIR_/$IID/BamStats.hs_metrics | head -n4 | csvtk transpose -t | \
	csvtk pretty -t > _WRKDIR_/$IID/BamStats.hs_metrics.txt

EOF
	}
	else {
		$script = <<'EOF'

read IID BAMFILE < _PARDIR_/IID2BAM._INDEX_

picard _RSRC.PICARDOPT_ CollectWgsMetrics \
       I=$BAMFILE R=_PATH.FASTA_ O=_WRKDIR_/$IID/BamStats.wgs_metrics \
       INTERVALS=_WRKDIR_/targets.interval_list

grep -v '^#' _WRKDIR_/$IID/BamStats.wgs_metrics | head -n4 | csvtk transpose -t | csvtk pretty -t \
	> _WRKDIR_/$IID/BamStats.wgs_metrics.txt

EOF
	}
}

sub depth_cov_files {
	my ($iid) = @_;
	my @files;
	unless($opt->get_wgs) {
		push @files, "wrk/$iid/BamStats.hs_metrics", "wrk/$iid/BamStats.hs_metrics.txt";
		return @files;
			
	}
	else {
		push @files, "wrk/$iid/BamStats.wgs_metrics", "wrk/$iid/BamStats.wgs_metrics.txt";
		return @files;
	}
}

sub mos_depth {
	my $script;
	unless($opt->get_wgs) {
		$script = <<'EOF'

read IID BAMFILE < _PARDIR_/IID2BAM._INDEX_

export MOSDEPTH_Q0=NO_COVERAGE
export MOSDEPTH_Q1=LOW_COVERAGE
export MOSDEPTH_Q2=CALLABLE

mosdepth -n --quantize _PARAM.DPCUTWES_ --by _PATH.TARGETS_ --mapq _PARAM.MINMQ_ \
	--fasta _PATH.FASTA_ --threads _RSRC.NT_  _WRKDIR_/$IID/DoC $BAMFILE

perl _PATH.MODULE_/calc_chrom_normdp.pl _WRKDIR_/$IID/DoC.regions.bed.gz _WRKDIR_/$IID/chrom_NormDP.txt

EOF
	}
	else {
		$script = <<'EOF'

read IID BAMFILE < _PARDIR_/IID2BAM._INDEX_

export MOSDEPTH_Q0=NO_COVERAGE
export MOSDEPTH_Q1=LOW_COVERAGE
export MOSDEPTH_Q2=CALLABLE
export MOSDEPTH_Q3=HIGH_COVERAGE

mosdepth -n --quantize _PARAM.DPCUTWGS_ --by _PARAM.WINSIZE_ --mapq _PARAM.MINMQ_ --threads _RSRC.NT_  \
	--fasta _PATH.FASTA_ _WRKDIR_/$IID/DoC $BAMFILE

perl _PATH.MODULE_/calc_chrom_normdp.pl _WRKDIR_/$IID/DoC.regions.bed.gz _WRKDIR_/$IID/chrom_NormDP.txt

EOF
	}
	return $script;
}


sub mos_depth_files {
	my ($iid) = @_;
	return ("wrk/$iid/DoC.regions.bed.gz", "wrk/$iid/DoC.quantized.bed.gz",
		"wrk/$iid/DoC.mosdepth.region.dist.txt", "wrk/$iid/DoC.mosdepth.global.dist.txt",
		"wrk/$iid/chrom_NormDP.txt");
}

sub contam_test {
	my ($filetype) = @_;
	my $script;
	if ($filetype eq 'bam') { 
		$script = <<'EOF';

read IID BAMFILE < _PARDIR_/IID2BAM._INDEX_

verifyBamID --vcf _WRKDIR_/targets.snps.vcf --bam $BAMFILE --out _WRKDIR_/$IID/VBID \
	--ignoreRG --chip-none --free-full --maxDepth 1000 --precise 

EOF
	}
	else {
		$script =<<'EOF';

read IID BAMFILE < _PARDIR_/IID2BAM._INDEX_

# Extract targeted region from CRAM file using GATK4
# Note there is a bug in GATK3 for converting CRAMs
gatk --java-options "_RSRC.GATKOPT_" PrintReads \
	-I $BAMFILE  -L _WRKDIR_/targets.snps.vcf \
	-R _PATH.FASTA_  -O _WRKDIR_/$IID.bam

verifyBamID --vcf _WRKDIR_/targets.snps.vcf --bam _WRKDIR_/$IID.bam --out _WRKDIR_/$IID/VBID \
	--ignoreRG --chip-none --free-full --maxDepth 1000 --precise 

sleep 2

rm -f _WRKDIR_/$IID.bam
rm -f _WRKDIR_/$IID.bai

EOF
	}
	return $script;
}

sub contam_test_files {
	my ($iid) = @_;
	my @outfiles = map { "wrk/$iid/VBID.$_" } qw|selfSM depthSM log|;
	return @outfiles;
}

sub sum_stats {
	my $script;
	unless($opt->get_wgs) {
		$script =<<'EOF';

perl _PATH.MODULE_/collect_bamqc_stats.pl _WRKDIR_  _OUTDIR_ 

EOF
	}
	else {
		$script =<<'EOF';

perl _PATH.MODULE_/collect_bamqc_stats.pl _WRKDIR_  _OUTDIR_ WGS

EOF
	}
	return $script;
}

sub plot_stats {
	my $script=<<'EOF';

library(ggplot2)

X<-read.csv("_OUTDIR_/bam_stats.csv", header=T)

Y<-subset(X, DP.Mean<_PARAM.MINDP_ | DP.PCT10x<_PARAM.MINPCT10X_ | ID.FREEMIX>_PARAM.MAXFREEMIX_)

if (nrow(Y) > 0) {
	write.csv(Y, "_OUTDIR_/badsamps.csv", quote=FALSE, row.names=FALSE)
}

Z<-subset(X, ! IID %in% Y$IID)
gender<-kmeans(Z[, c("NDP.chrX","NDP.chrY")], 2, nstart = 20) 

if (gender$centers[1,"NDP.chrX"]>gender$centers[2,"NDP.chrX"]) {
	sex<-ifelse(gender$cluster == 1, "Female", "Male")
} else {
	sex<-ifelse(gender$cluster == 1, "Male", "Female")
} 

g1<-ggplot(Z) + geom_point(aes(x=NDP.chrX, y=NDP.chrY, color=sex)) +
	labs(x="Normalized chrX depth", y="Normalized chrY depth", color="Gender")
ggsave("_OUTDIR_/samples_chrXYNormDp_byPredSex.png", g1, dpi=300)

W<-rbind(data.frame(MeanDP=Z$DP.Mean, PCTx=Z$DP.PCT10x, Cutoff=rep(10,nrow(Z))),
	data.frame(MeanDP=Z$DP.Mean, PCTx=Z$DP.PCT20x, Cutoff=rep(20,nrow(Z))) )
g2<-ggplot(W) + geom_point(aes(x=MeanDP, y=PCTx, color=as.factor(Cutoff))) +
	labs(x="Mean on-target depth", y="Percent of targets at 10x/20x", color="Cutoff")
ggsave("_OUTDIR_/samples_DepthOfCov.png", g2, dpi=300)

V<-read.csv("_OUTDIR_/chrom_depth.csv")
pdf("_OUTDIR_/chrom_NormDP_CDF.pdf")
for(C in seq(1,22)) {
	CHR<-paste0('NDP_chr',C)
	DAT<-V[[CHR]]
	plot(1:length(DAT)/(length(DAT)-1), sort(DAT), main=CHR,
    	xlab="Quantile", ylab=paste0(CHR,"normalized depth"))
}
dev.off()

EOF
	return $script;
}

