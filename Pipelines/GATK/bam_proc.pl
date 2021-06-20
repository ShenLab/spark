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
use List::MoreUtils qw|all none uniq|;
use Config::Std;
use String::ShellQuote;
use Utils::Workflow;
use Utils::Hash qw|merge_conf|;

use lib "$Bin/../lib/";
use Shared qw|read_list split_chrs|;
use Wrapper qw|bam_rgids bam_ismapped|;

############################
## Command line interface ##
############################

my @spec =  (
	Param("conf|c")->valid(sub { -r }),
	Param("list|l")->valid(sub { -r }),
	Param("rename")->valid(sub { -r }),
	Param("outdir|out"),
	Param("engine|eng")->valid(sub { $_ eq 'SGE' || $_ eq 'BASH' }),
	Keypair("param|par"),
	Switch("mapped"),
	Switch("wgs"),
	Switch("novbsr"),
	Switch("help|h"),
	Switch("dryrun|dry"),
	Switch("force")
	);

my $opt = Getopt::Lucid->getopt(\@spec);

if ($opt->get_help) {
	print STDERR <<EOF;
Purpose:
	This is a pipeline script for processing raw sequence data to create seq alignment files ready for
	variant calling.

Usage:
	bam_proc.pl --conf Conf --list BAMList --outdir RootDit [--mapped --wgs --eng Engine]

Options:
	--list: BAM file list typically has one or two columns per line: path to the BAM file and 
			sample ID (optional). Multiple BAM files associated with the same ID is presumed from 
			different read groups. If no sample ID is provided, it will take file base name after 
			removing suffix. 
	--mapped: Indicate that the input BAM or CRAM files contain mapped reads. By default, it will
			check file headers to determine if they contain mapped reads. Note that only one file 
			type is allowed for each run! If this switch is turned on, it will skip the checking and
			pressume all input BAMs contain mapped reads.
	--wgs:	Indicate that workflow should be adjusted to handle WGS data. 
	--novbsr: Turn on this switch to skip the base quality recalibration (BQSR) step. 
			Recommended for processing data generated by the latest sequencing platform or targeted 
			sequencing of small number of genes.

Notes:
	The script implements GATK "best practice" for processing raw sequencing data. Raw sequencing reads
	are mapped by bwa, then sorted and fixed by samtools/picard. After that, likely duplicated reads are 
	marked up by picard. Base quality scores are recalibrated by GATK. Finally, initial variant calls
	are genereated by GATK's HaplotypeCaller.
	
	Although GATK's best practice still recommends using BQSR, DeepVariant's developers recommend skip
	this step because it does improve performance of DeepVariant. BQSR also increase the resulting 
	BAM/CRAM file size to storing both original and calibrated base qualities. The step becomes very 
	time consuming when processing WGS data. 

	The pipeline can process both WES (including any targeted sequencing) and WGS sequence data. Only differences
	are BQSR and HCVarCall steps. For WGS data, these tasks will be split into smaller segments and use
	scatter-and-gather to speed up the process. For WES data, initial variant calls will be made only
	within targeted interavals plus flanking regions specified in config. For wGS data, variant calls will be
	made across the genome.

Input/output:
	The requrired input must be pair-end sequencing data in mapped or unmapped BAM format. Mapped CRAM file 
	with reference FASTA is also supported. We support mapped data because of the need for remapping legacy data.
	The advantage of uBAM over fastq is that it can store reads in a more compact way with read	group information.
	(https://software.broadinstitute.org/gatk/documentation/article.php?id=6484)
	Only one sample is allowed in a BAM file. Input uBAM files must additionally comply with the following 
	requirements:
		1. Files must pass validation by ValidateSamFile
		2. Reads are provided in query-sorted order
		3. Reads should be associated with RG information which must include RG ID, sample name, library name, 
			and platform fileds.
		(Optionally) uBAM can also contain unique molecular identifier (UMI) information for each reads. 
		The custom barcode tag can be used for marking up duplicates. 

	Sample IDs is used as basename for output files. The following files will be generated for each individual:
		IID.cram/bam: Aligned sequence reads in BAM and CRAM formats.
		IID.bam.flagstat/cram.flagstat: Flag stats of aligned sequence reads. 
			Flagstat is useful to verify total number of reads in the final output is the same as input.
		IID.g.vcf.gz: Initial variant calling gVCF. It will be used in joint genotyping across multiple samples.
		IID.dup_metrics: Duplication metrics generated from picard MarkDuplicates.
			Note: this metrics is not available from bam_qc.pl.

Dependencies:
	bwa, samtools, picard, and gatk (v4).

EOF
	exit 1;
}

$opt->validate({requires => [qw|conf list outdir|]});

my $rootdir = $opt->get_outdir;

############################################
## Input files validation and collection  ##
############################################

my %conf = merge_conf($opt->get_conf, $opt->get_param); 

while(my ($label, $file) = each %{$conf{PATH}}) {
	next if $label eq 'PREFIX';
	if (ref $file eq 'ARRAY') {
		unless(all { -f $_ } @$file) {
			croak "Cannot find $label files";
		}
	}
	else {
		unless(-f $file) {
			croak "Cannot find $label file";
		}
	}
}
# Further check BWA related files
if (exists $conf{PATH}{PREFIX}) {
	my @bwasuf = qw|amb ann bwt pac sa alt|;
	unless(all { -f "$conf{PATH}{PREFIX}.$_" } @bwasuf) {
		croak "Cannot find bwa index files for $conf{PATH}{FASTA}";
	}
}
else {
	my @bwasuf = qw|amb ann bwt pac sa|;
	unless(all { -f "$conf{PATH}{FASTA}.$_" } @bwasuf) {
		croak "Cannot find bwa index files for $conf{PATH}{FASTA}";
	}
}


# Ordered sample list
my @samps;
# If mapped, ordered {sample, mBAM, RGs} triplets
my @mbams;
# Ordered {sample, RG} pairs
my @rgs;
# uBAM files: one file for each sample and RG combination
my %ubams;
# mapped flag
my $mapped;
{
	# BAM files and file type
	my %opt = exists $conf{AWS}{STAGEIN} ? (notest => 1) : ();
	my ($bams, $filetype) = read_list($opt->get_list, 
		{suffix => [qw|bam cram|], rename => $opt->get_rename, multi => 1, %opt});

	if ($opt->get_mapped) {
		$mapped = 1;
	}
	else {
		my @ismapped = map { bam_ismapped($_) } map { @$_ } values %$bams;
		if (all { $_ } @ismapped) {
			print STDERR "All data files contain mapped reads\n";
			$mapped = 1;
		}
		elsif (none { $_ } @ismapped) {
			print STDERR "All data files contain unmapped reads\n";
			$mapped = 0;
		}
		else {
			croak "Data files should be either all mapped or all unmapped";
		}
	}

	foreach my $iid (sort keys %$bams) {
		push @samps, $iid;
		if ($mapped) {
			my @samprgs;
			foreach my $bamfile (@{$bams->{$iid}}) {
				my @rgids = bam_rgids($bamfile);
				push @samprgs, @rgids;
				push @mbams, [$iid, $bamfile, [@rgids]];
				foreach my $rgid (@rgids) {
					push @rgs, [$iid, $rgid];
					$ubams{$iid}{$rgid} = "$rootdir/wrk/$iid/$rgid.bam";
				}
			}
			my @uqrgid = uniq sort @samprgs;
			unless(@uqrgid == @samprgs) {
				croak "Read group IDs of $iid are not unique";
			}
		}
		else {
			foreach my $bamfile (@{$bams->{$iid}}) {
				my @rgids = bam_rgids($bamfile);
				croak "uBAM file should contain one and only one RG: $bamfile" unless @rgids == 1;
				push @rgs, [$iid, $rgids[0]];
				if (defined $ubams{$iid}{$rgids[0]}) {
					croak "Read group $rgids[0] for sample $iid has been found before!";
				}
				$ubams{$iid}{$rgids[0]} = $bamfile;
			}
		}
	}
}

print STDERR "Total number of samples: ", scalar(@samps), "\n";
print STDERR "Total number of read groups: ", scalar(@rgs), "\n";

exit 1 unless @samps > 0 && @rgs > 0;

#################################
##  Workflow initialization    ##
##  Working directory setup    ##
#################################

my $wkf = Utils::Workflow->new($rootdir, 
	{ engine => $opt->get_engine, force => $opt->get_force, strict_var => 1});

# Set up parameter files in the working directory files
# for mapped reads, it will store mapped bam per sample and unmapped bam per read group
# for unmapped reads, it will store unmapped bam per read group
if ($mapped) {
	my $prevsamp = "";
	for (my $ii = 0; $ii < @mbams; $ii ++) {
		my $jj = $ii + 1;
		my ($samp, $mbam, $rgs) = @{$mbams[$ii]};
		if ($samp ne $prevsamp) {
			open my $fout, ">$rootdir/par/IID2MBAM.$jj" or die "Cannot write to IID2MBAM";
			print $fout $samp, "\t", $mbam, "\n";
		}
		$prevsamp = $samp;
	}
}
for(my $ii = 0; $ii < @samps; $ii ++) {
	my $jj = $ii + 1;
	open my $fout, ">$rootdir/par/SAMPIID.$jj" or die "Cannot write to SAMPIID";
	print $fout $samps[$ii], "\t", scalar(keys %{$ubams{$samps[$ii]}}), "\n";
	make_path "$rootdir/wrk/$samps[$ii]";
	make_path "$rootdir/tmp/$samps[$ii]";
}
for(my $ii = 0; $ii < @rgs; $ii ++) {
	my $jj = $ii + 1;
	open my $fout, ">$rootdir/par/IID2UBAM.$jj" or die "Cannot write to IID2UBAM";
	print $fout join("\t", @{$rgs[$ii]}, $ubams{$rgs[$ii][0]}{$rgs[$ii][1]}), "\n"; 
	make_path "$rootdir/tmp/$rgs[$ii][0]/$rgs[$ii][1]";
}

# Under WGS mode, create sequencing group
if ($opt->get_wgs) {
	my @chrgrp = split_chrs($conf{PATH}{SEQDICT});
	# Now write to seq group files.
	for(my $ii = 0; $ii < @chrgrp; $ii ++) {
		my $jj = $ii + 1;
		open my $fout, ">$rootdir/par/SEQGRP.$jj" or die "Cannot write to SEQGRP";
		open my $fout2, ">$rootdir/par/SEQGRPUM.$jj" or die "Cannot write to SEQGRPUM";
		print $fout join(" ", map { "-L $_" } @{$chrgrp[$ii]}), "\n";
		if ($ii == $#chrgrp) {
			print $fout2 join(" ", map { "-L $_" } @{$chrgrp[$ii]}), " -L unmapped\n";
		}
		else {
			print $fout2 join(" ", map { "-L $_" } @{$chrgrp[$ii]}), "\n";
		}
	}
	$conf{PARAM}{NSPLIT} = @chrgrp;
}

# Get BWA version
$conf{PARAM}{BWAVER}=`bwa 2>&1 | grep -e '^Version' | sed 's/Version: //'`;
chomp($conf{PARAM}{BWAVER});

# Add module path to config
$conf{PATH}{MODULE} = shell_quote("$Bin/module");

write_config %conf, "$rootdir/par/run.conf" unless $opt->get_dryrun;

# Adding tasks
if ($mapped) {
	$wkf->add(prep_ubam(), { name => "PrepUBams", 
			expect => [get_sprgbam_exp("bam")], nslots => scalar(@mbams),
			callback => \&check_nonempty  })
 		->add(bwa_aln(), { name => "BwaAln", depend => "PrepUBams",
			expect => [get_rgbam_exp("mapped.bam", "mapped.bam.flagstat")], nslots => scalar(@rgs),
			callback => \&check_nonempty });
}
else {
	$wkf->add(bwa_aln(), { name => "BwaAln", 
		expect => [get_rgbam_exp("mapped.bam", "mapped.bam.flagstat")], nslots => scalar(@rgs),
		callback => \&check_nonempty });
}
my %arg;
if (scalar(@samps) == scalar(@rgs)) {
	%arg = (deparray => 1);
}
$wkf->add(merge_bam(), { name => "MergeBam", depend => "BwaAln", %arg,
		expect => [get_rgbam_exp("merged.bam", "merged.bam.flagstat")], nslots => scalar(@rgs),
		callback => \&check_nonempty  })
	->add(sort_fix('RG'), { name => "SortFixRG",  depend => "MergeBam", deparray => 1,
		expect => [get_rgbam_exp("sorted.bam", "sorted.bai", "sorted.bam.flagstat")], 
		nslots => scalar(@rgs), callback => \&check_nonempty  })
	->add(mark_dup(), { name => "MarkDup", depend => "SortFixRG", %arg,
		expect => [get_sampbam_exp("wrk", "dupmarked.bam", "dup_metrics", "dupmarked.bam.flagstat")],
		nslots => scalar(@samps), callback => \&check_nonempty  })
	->add(collect_metrics(), { name => "DupMetrics", depend => "MarkDup", 
		expect => "out/dup_metrics.csv", callback => \&check_nonempty  })
	->add(sort_fix('Sample'), { name => "SortFixSamp", depend => "MarkDup", deparray => 1,
		expect => [get_sampbam_exp("wrk", "sorted.bam", "sorted.bai", "sorted.bam.flagstat")], 
		nslots => scalar(@samps), callback => \&check_nonempty  });

unless($opt->get_novbsr) {
	if ($opt->get_wgs) {
		# For WGS data
		my $unit = $conf{PARAM}{NSPLIT};
		$wkf->add(bq_recal_split(), { name => "BQRecalSplit", depend => "MarkDup",
				expect => [get_sampbam_splitexp($unit, "recal_data.csv")], nslots => @samps*$unit,
				callback => \&check_nonempty })
			->add(bq_recal_gather(), { name => "BQRecalGather", depend => "BQRecalSplit", deparray => 1,
				expect => [get_sampbam_exp("wrk", "recal_data.csv")], step => $unit, nslots => @samps*$unit,
				callback => \&check_nonempty })
			->add(apply_bqsr_split(), { name => "AppBQSRSplit", depend => "BQRecalGather", deparray => 1,
				expect => [get_sampbam_splitexp($unit, "recal.bam")], nslots => @samps*$unit,
				callback => \&check_nonempty  })
			->add(apply_bqsr_gather(), { name => "AppBQSRGather", depend => "AppBQSRSplit", deparray => 1,
				expect => [get_sampbam_exp("out", "bam", "bai", "bam.flagstat")], step => $unit, 
				nslots => @samps*$unit, callback => \&check_nonempty  })
			->add(convert_cram(), { name => "ConvertCram", depend => "AppBQSRGather", 
				expect => [get_sampbam_exp("out", "cram", "cram.crai", "cram.flagstat")], 
				nslots => scalar(@samps), callback => \&check_nonempty  })
			->add(hc_varcall_split(), { name => "HCSplit", depend => "AppBQSRSplit", deparray => 1,
				expect => [get_sampbam_splitexp($unit, "g.vcf.gz", "g.vcf.gz.tbi")], 
				nslots => @samps*$unit, callback => \&check_nonempty  })
			->add(hc_varcall_merge(), { name => "HCGather", depend => "HCSplit", deparray => 1,
				expect => [get_sampbam_exp("out", "g.vcf.gz", "g.vcf.gz.tbi")], step => $unit, 
				nslots => @samps*$unit, callback => \&check_nonempty  });
			if (exists $conf{AWS}{STAGEOUT}) {
				$wkf->add(stage_out_merged(), { name => "StageOut", depend => "HCGather", deparray => 1,
						expect => [get_sampbam_exp("out", "stagedout")], step => $unit,  nslots => @samps*$unit });
			}
	}
	else {
		# For exome data
		$wkf->add(bq_recal(), { name => "BaseQRecal", depend => "SortFixSamp", deparray => 1,
				expect => [get_sampbam_exp("wrk", "recal_data.csv")], nslots => scalar(@samps),
				callback => \&check_nonempty })
			->add(apply_bqsr(), { name => "ApplyBQSR", depend => "BaseQRecal", deparray => 1,
				expect => [get_sampbam_exp("out", "bam", "bai", "bam.flagstat")],
				nslots => scalar(@samps), callback => \&check_nonempty  })
			->add(convert_cram(), { name => "ConvertCram", depend => "ApplyBQSR", deparray => 1,
				expect => [get_sampbam_exp("out", "cram", "cram.crai", "cram.flagstat")], 
				nslots => scalar(@samps), callback => \&check_nonempty  })
			->add(hc_varcall(), { name => "HCVarCall", depend => "ApplyBQSR", deparray => 1,
				expect => [get_sampbam_exp("out", "g.vcf.gz", "g.vcf.gz.tbi")], 
				nslots => scalar(@samps), callback => \&check_nonempty  });
		if (exists $conf{AWS}{STAGEOUT}) {
			$wkf->add(stage_out(), { name => "StageOut", depend => "HCVarCall", deparray => 1,
					expect => [get_sampbam_exp("out", "stagedout")], nslots => scalar(@samps) });
		}
	}
}
else {
	$wkf->add(link_bams(), { name => "LinkBAMs", depend => "SortFixSamp", deparray => 1,
				expect => [get_sampbam_exp("out", "bam", "bai", "bam.flagstat")],
				nslots => scalar(@samps), callback => \&check_nonempty  })
		->add(convert_cram(), { name => "ConvertCram", depend => "LinkBAMs", deparray => 1,
				expect => [get_sampbam_exp("out", "cram", "cram.crai", "cram.flagstat")], 
				nslots => scalar(@samps), callback => \&check_nonempty  })
		->add(hc_varcall(), { name => "HCVarCall", depend => "LinkBAMs", deparray => 1,
				expect => [get_sampbam_exp("out", "g.vcf.gz", "g.vcf.gz.tbi")], 
				nslots => scalar(@samps), callback => \&check_nonempty  });
		if (exists $conf{AWS}{STAGEOUT}) {
			$wkf->add(stage_out(), { name => "StageOut", depend => "HCVarCall", deparray => 1,
					expect => [get_sampbam_exp("out", "stagedout")],
					nslots => scalar(@samps) });
		}
}



################################
## Kickstart workflow engine  ##
################################

$wkf->inst(\%conf);

$wkf->run({ conf => $conf{$wkf->{engine}}, dryrun => $opt->get_dryrun  });


############################
## Workflow components    ##
############################

# Prepare unmapped BAMs from mapped BAMs
# This step may require large memory and large temp directory
sub prep_ubam {
	my $script = "\nread IID MBAM < _PARDIR_/IID2MBAM._INDEX_\n";
	if ($mapped && exists $conf{AWS}{STAGEIN}) {
		$script .= <<'EOF';
TMPDIR=$(stagein.pl --input $MBAM --output _TMPDIR_/MBAM._INDEX_ --profile _AWS.PROFILEIN_ --tmpdir _AWS.STAGEIN_ --all)

trap "echo Clean up $TMPDIR; rm -fR $TMPDIR; exit 1" EXIT SIGINT SIGTERM

read MBAM < _TMPDIR_/MBAM._INDEX_

EOF
	}
	$script .= <<'EOF';
picard _RSRC.PICARDREVERTOPT_ \
	-Djava.io.tmpdir=_TMPDIR_/$IID \
	RevertSam \
	INPUT=$MBAM \
	REFERENCE_SEQUENCE=_PATH.OLDFASTA_ \
	OUTPUT=_WRKDIR_/$IID \
	OUTPUT_BY_READGROUP=true \
	OUTPUT_BY_READGROUP_FILE_FORMAT=bam \
	SANITIZE=true \
	MAX_DISCARD_FRACTION=0.005 \
	ATTRIBUTE_TO_CLEAR=XT \
	ATTRIBUTE_TO_CLEAR=XN \
	ATTRIBUTE_TO_CLEAR=AS \
	ATTRIBUTE_TO_CLEAR=OC \
	ATTRIBUTE_TO_CLEAR=OP \
	SORT_ORDER=queryname \
	RESTORE_ORIGINAL_QUALITIES=true \
	REMOVE_DUPLICATE_INFORMATION=true \
	REMOVE_ALIGNMENT_INFORMATION=true \
	VALIDATION_STRINGENCY=LENIENT \
	SAMPLE_ALIAS=$IID

EOF
	return $script;
}

# Get read-group level expected bam file output
sub get_sprgbam_exp {
	my @exp;
	foreach my $info (@mbams) {
		my ($iid, $bamfile, $rgids) = @$info;
		my @bams;
		foreach my $rgid (@$rgids) {
			foreach my $suffix (@_) {
				push @bams, "wrk/$iid/$rgid.$suffix";
			}
		}
		push @exp, [@bams];
	}
	return @exp;
}

sub bwa_aln {
	my $script = "\nread IID RGID UBAM < _PARDIR_/IID2UBAM._INDEX_\n";
	if (!$mapped && exists $conf{AWS}{STAGEIN}) {
		$script .= <<'EOF';
TMPDIR=$(stagein.pl --input $UBAM --output _TMPDIR_/UBAM._INDEX_ --profile _AWS.PROFILEIN_ --tmpdir _AWS.STAGEIN_ --all)

trap "echo Clean up $TMPDIR; rm -fR $TMPDIR; exit 1" EXIT SIGINT SIGTERM

read UBAM < _TMPDIR_/UBAM._INDEX_

EOF
	}
	$script .= <<'EOF';
picard _RSRC.PICARDOPT_ SamToFastq \
	INPUT=$UBAM \
	FASTQ=/dev/stdout \
	INTERLEAVE=true \
	NON_PF=true | \
  	bwa mem -K 100000000 -p -v 3 \
  	-t _RSRC.BWANT_ \
  	-Y _PATH.FASTA_ /dev/stdin - | \
  	samtools view -1 - > _WRKDIR_/$IID/$RGID.mapped.bam

samtools flagstat _WRKDIR_/$IID/$RGID.mapped.bam \
	> _WRKDIR_/$IID/$RGID.mapped.bam.flagstat

EOF
	# If PREFIX is provided, will switch to BWA index files given by PREFIX
	if (exists $conf{PATH}{PREFIX}) {
		$script =~ s/FASTA/PREFIX/;
	}
	return $script;
}

sub merge_bam {
	return <<'EOF';

read IID RGID UBAM < _PARDIR_/IID2UBAM._INDEX_

picard _RSRC.PICARDOPT_ MergeBamAlignment \
	VALIDATION_STRINGENCY=SILENT \
	EXPECTED_ORIENTATIONS=FR \
	ATTRIBUTES_TO_RETAIN=X0 \
	ALIGNED_BAM=_WRKDIR_/$IID/$RGID.mapped.bam \
	UNMAPPED_BAM=$UBAM \
	OUTPUT=_WRKDIR_/$IID/$RGID.merged.bam \
	REFERENCE_SEQUENCE=_PATH.FASTA_ \
	PAIRED_RUN=true \
	SORT_ORDER="unsorted" \
	IS_BISULFITE_SEQUENCE=false \
	ALIGNED_READS_ONLY=false \
	CLIP_ADAPTERS=false \
	MAX_RECORDS_IN_RAM=2000000 \
	ADD_MATE_CIGAR=true \
	MAX_INSERTIONS_OR_DELETIONS=-1 \
	PRIMARY_ALIGNMENT_STRATEGY=MostDistant \
	PROGRAM_RECORD_ID="bwamem" \
	PROGRAM_GROUP_NAME="bwamem" \
	PROGRAM_GROUP_VERSION=_PARAM.BWAVER_ \
	PROGRAM_GROUP_COMMAND_LINE="bwa mem -K 100000000 -p -v 3" \
	UNMAPPED_READ_STRATEGY=COPY_TO_TAG \
	ALIGNER_PROPER_PAIR_FLAGS=true \
	UNMAP_CONTAMINANT_READS=true

samtools flagstat _WRKDIR_/$IID/$RGID.merged.bam  \
	> _WRKDIR_/$IID/$RGID.merged.bam.flagstat

EOF
}


# Get read-group level expected bam file output
sub get_rgbam_exp {
	my @exp;
	foreach my $info (@rgs) {
		my ($iid, $rgid) = @$info;
		my @bams;
		foreach my $suffix (@_) {
			push @bams, "wrk/$iid/$rgid.$suffix";
		}
		push @exp, [@bams];
	}
	return @exp;
}

# Aggregate aligned+merged flowcell BAM files and mark duplicates
# We take advantage of the tool's ability to take multiple BAM inputs and write out a single output
# to avoid having to spend time just merging BAM files.
sub mark_dup {
	my $script =<< 'EOF';

read IID NRG < _PARDIR_/SAMPIID._INDEX_

# Check if sorted bam files are available for all read groups
NBAMS=$(ls _WRKDIR_/$IID/*.sorted.bam | wc -l)
if [[ $NRG != $NBAMS ]]; then
	echo "Not all sorted bams are available" 2>&1
	exit 1
fi

INPUT=`perl -e 'print join(" ", map {"INPUT=$_" } @ARGV)' _WRKDIR_/$IID/*.sorted.bam`

picard _RSRC.PICARDOPT_ MarkDuplicates \
	$INPUT \
	OUTPUT=_WRKDIR_/$IID.dupmarked.bam \
	TMP_DIR=_TMPDIR_/$IID \
	METRICS_FILE=_WRKDIR_/$IID.dup_metrics \
	VALIDATION_STRINGENCY=SILENT \
	ASSUME_SORT_ORDER="queryname"

samtools flagstat _WRKDIR_/$IID.dupmarked.bam > _WRKDIR_/$IID.dupmarked.bam.flagstat

EOF
	if (exists $conf{PARAM}{BCTAG}) {
		$script =~ s/VALIDATION_STRINGENCY/BARCODE_TAG=_PARAM.BCTAG_ VALIDATION_STRINGENCY/;
	}
	return $script;
}

sub collect_metrics {
	return <<'EOF'

perl _PATH.MODULE_/collect_dup_metrics.pl _WRKDIR_ _OUTDIR_

EOF
}


sub sort_fix {
	my ($level) = @_;
	if ($level eq 'RG') {
		return <<'EOF';

read IID RGID UBAM < _PARDIR_/IID2UBAM._INDEX_

picard _RSRC.PICARDOPT_ SortSam \
	INPUT=_WRKDIR_/$IID/$RGID.merged.bam \
    OUTPUT=/dev/stdout \
    TMP_DIR=_TMPDIR_/$IID/$RGID \
    SORT_ORDER="coordinate" \
    CREATE_INDEX=false \
    CREATE_MD5_FILE=false | \
picard _RSRC.PICARDOPT_ SetNmMdAndUqTags \
    INPUT=/dev/stdin \
    OUTPUT=_WRKDIR_/$IID/$RGID.sorted.bam \
    CREATE_INDEX=true \
    REFERENCE_SEQUENCE=_PATH.FASTA_

samtools flagstat _WRKDIR_/$IID/$RGID.sorted.bam \
	> _WRKDIR_/$IID/$RGID.sorted.bam.flagstat

EOF
	}
	elsif ($level eq 'Sample') {
		return <<'EOF';

read IID NRG < _PARDIR_/SAMPIID._INDEX_

picard _RSRC.PICARDOPT_ SortSam \
	INPUT=_WRKDIR_/$IID.dupmarked.bam \
    OUTPUT=/dev/stdout \
    TMP_DIR=_TMPDIR_/$IID \
    SORT_ORDER="coordinate" \
    CREATE_INDEX=false \
    CREATE_MD5_FILE=false | \
picard _RSRC.PICARDOPT_ SetNmMdAndUqTags \
    INPUT=/dev/stdin \
    OUTPUT=_WRKDIR_/$IID.sorted.bam \
    CREATE_INDEX=true \
    REFERENCE_SEQUENCE=_PATH.FASTA_

samtools flagstat _WRKDIR_/$IID.sorted.bam \
	> _WRKDIR_/$IID.sorted.bam.flagstat

EOF
	}
	else {
		croak "Cannot recognize $level";
	}
}

sub link_bams {
	return <<'EOF';

read IID NRG < _PARDIR_/SAMPIID._INDEX_

cd _OUTDIR_

ln -s ../wrk/$IID.sorted.bam $IID.bam
ln -s ../wrk/$IID.sorted.bai $IID.bai
ln -s ../wrk/$IID.sorted.bam.flagstat $IID.bam.flagstat

EOF
}

sub bq_recal {
	return <<'EOF';

read IID NRG < _PARDIR_/SAMPIID._INDEX_

gatk --java-options "_RSRC.GATKOPT_" BaseRecalibrator \
	-R _PATH.FASTA_ \
	-I _WRKDIR_/$IID.sorted.bam \
	--use-original-qualities \
	-O _WRKDIR_/$IID.recal_data.csv \
    --known-sites _PATH.KNOWNSITES[ --known-sites ]_

EOF
}

sub apply_bqsr {
	return <<'EOF';

read IID NRG < _PARDIR_/SAMPIID._INDEX_

gatk --java-options "_RSRC.GATKOPT_" ApplyBQSR \
	-R _PATH.FASTA_ \
	-I _WRKDIR_/$IID.sorted.bam \
	-O _OUTDIR_/$IID.bam \
	-bqsr _WRKDIR_/$IID.recal_data.csv \
	--static-quantized-quals 10 --static-quantized-quals 20 --static-quantized-quals 30 \
    --add-output-sam-program-record \
    --use-original-qualities

samtools flagstat _OUTDIR_/$IID.bam > _OUTDIR_/$IID.bam.flagstat

EOF
}

sub convert_cram {
	return <<'EOF';

read IID NRG < _PARDIR_/SAMPIID._INDEX_

samtools view -C -T _PATH.FASTA_ -o _OUTDIR_/$IID.cram _OUTDIR_/$IID.bam 

samtools index _OUTDIR_/$IID.cram

samtools flagstat _OUTDIR_/$IID.cram > _OUTDIR_/$IID.cram.flagstat

EOF
}

sub hc_varcall {
	return <<'EOF';

read IID NRG < _PARDIR_/SAMPIID._INDEX_

gatk --java-options "_RSRC.GATKOPT_" HaplotypeCaller \
	-R _PATH.FASTA_ \
	-I _OUTDIR_/$IID.bam \
	-O _OUTDIR_/$IID.g.vcf.gz \
	-L _PATH.TARGETS_ -ip _PARAM.PADDING_ \
	--max-alternate-alleles 3 \
	-ERC GVCF 

EOF
}

sub stage_out {
	return <<'EOF';

read IID NRG < _PARDIR_/SAMPIID._INDEX_

stageout.pl --files $IID.g.vcf.gz,$IID.g.vcf.gz.tbi,$IID.cram,$IID.cram.crai --profile _AWS.PROFILEIN_ \
	--indir _OUTDIR_ --outdir _AWS.STAGEOUT_

touch _OUTDIR_/$IID.stagedout

EOF
}

# Get sample level expected bam file output
sub get_sampbam_exp {
	my $outdir = shift @_;
	my @exp;
	foreach my $iid (@samps) {
		my @bams;
		foreach my $suffix (@_) {
			push @bams, "$outdir/$iid.$suffix";
		}
		push @exp, [@bams];
	}
	return @exp;
}

sub check_nonempty {
	my (@exp) = @_;
	if (all { -s "$rootdir/$_" } @exp) {
		return 1;
	}
	else {
		return 0;
	}
}

sub bq_recal_split {
	return <<'EOF';

II=$(( (_INDEX_-1)/_PARAM.NSPLIT_+1 ))
read IID NRG < _PARDIR_/SAMPIID.$II

SPLIT=$(( (_INDEX_-1)%_PARAM.NSPLIT_+1 ))
INTERVAL=`cat _PARDIR_/SEQGRP.$SPLIT`

gatk --java-options "_RSRC.GATKOPT_" BaseRecalibrator \
	-R _PATH.FASTA_ \
	-I _WRKDIR_/$IID.sorted.bam \
	--use-original-qualities \
	-O _WRKDIR_/$IID/split.$SPLIT.recal_data.csv \
	--known-sites _PATH.KNOWNSITES[ --known-sites ]_ \
	$INTERVAL

EOF
}

sub bq_recal_gather {
	return <<'EOF';

II=$(( (_INDEX_-1)/_PARAM.NSPLIT_+1 ))
read IID NRG < _PARDIR_/SAMPIID.$II

INPUTS=$(perl -e 'print join(" ", map { "-I _WRKDIR_/$ARGV[0]/split.$_.recal_data.csv" } 1.._PARAM.NSPLIT_)' $IID)

gatk --java-options "_RSRC.GATKOPT_" GatherBQSRReports \
	$INPUTS \
	-O _WRKDIR_/$IID.recal_data.csv

EOF
}

sub apply_bqsr_split {
	return <<'EOF';

II=$(( (_INDEX_-1)/_PARAM.NSPLIT_+1 ))
read IID NRG < _PARDIR_/SAMPIID.$II

SPLIT=$(( (_INDEX_-1)%_PARAM.NSPLIT_+1 ))
INTERVAL=`cat _PARDIR_/SEQGRPUM.$SPLIT`


gatk --java-options "_RSRC.GATKOPT_" ApplyBQSR \
	-R _PATH.FASTA_ \
	-I _WRKDIR_/$IID.sorted.bam \
	-O _WRKDIR_/$IID/split.$SPLIT.recal.bam \
	$INTERVAL \
	-bqsr _WRKDIR_/$IID.recal_data.csv \
	--static-quantized-quals 10 --static-quantized-quals 20 --static-quantized-quals 30 \
    --add-output-sam-program-record \
    --create-output-bam-md5 \
    --use-original-qualities

EOF
}

sub apply_bqsr_gather {
	return <<'EOF';

II=$(( (_INDEX_-1)/_PARAM.NSPLIT_+1 ))
read IID NRG < _PARDIR_/SAMPIID.$II

INPUT=`perl -e 'print join(" ", map {"INPUT=_WRKDIR_/$ARGV[1]/split.$_.recal.bam" } 1..$ARGV[0])' _PARAM.NSPLIT_ $IID`

picard _RSRC.PICARDOPT_ GatherBamFiles \
	$INPUT \
	OUTPUT=_OUTDIR_/$IID.bam \
	CREATE_INDEX=true CREATE_MD5_FILE=true

samtools flagstat _OUTDIR_/$IID.bam > _OUTDIR_/$IID.bam.flagstat

EOF
}

sub hc_varcall_split {
	return <<'EOF';

II=$(( (_INDEX_-1)/_PARAM.NSPLIT_+1 ))
read IID NRG < _PARDIR_/SAMPIID.$II

SPLIT=$(( (_INDEX_-1)%_PARAM.NSPLIT_+1 ))
INTERVAL=`cat _PARDIR_/SEQGRP.$SPLIT`

gatk --java-options "_RSRC.GATKOPT_" HaplotypeCaller \
	-R _PATH.FASTA_ \
	-I _WRKDIR_/$IID/split.$SPLIT.recal.bam \
	-O _WRKDIR_/$IID/split.$SPLIT.g.vcf.gz \
	$INTERVAL \
	-ip _PARAM.PADDING_ \
	--max-alternate-alleles 3 \
	-ERC GVCF

EOF
}

sub hc_varcall_merge {
	return <<'EOF';

II=$(( (_INDEX_-1)/_PARAM.NSPLIT_+1 ))
read IID NRG < _PARDIR_/SAMPIID.$II

INPUT=`perl -e 'print join(" ", map {"INPUT=_WRKDIR_/$ARGV[1]/split.$_.g.vcf.gz" } 1..$ARGV[0])' _PARAM.NSPLIT_ $IID`

picard _RSRC.PICARDOPT_ MergeVcfs \
	$INPUT \
	OUTPUT=_OUTDIR_/$IID.g.vcf.gz

EOF
}


sub stage_out_merged {
	return <<'EOF';

II=$(( (_INDEX_-1)/_PARAM.NSPLIT_+1 ))
read IID NRG < _PARDIR_/SAMPIID.$II

stageout.pl --files $IID.g.vcf.gz,$IID.g.vcf.gz.tbi,$IID.cram,$IID.cram.crai --profile _AWS.PROFILEIN_ \
	--indir _OUTDIR_ --outdir _AWS.STAGEOUT_

touch _OUTDIR_/$IID.stagedout

EOF
}

sub get_sampbam_splitexp {
	my $nsplit = shift @_;
	my @exp;
	foreach my $iid (@samps) {
		foreach my $jj (1..$nsplit) {
			my @bams;
			foreach my $suffix (@_) {
				push @bams, "wrk/$iid/split.$jj.$suffix";
			}
			push @exp, [@bams];
		}
	}
	return @exp;
}


