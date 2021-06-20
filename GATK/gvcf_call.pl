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
use List::MoreUtils qw|all uniq|;
use Config::Std;
use Utils::Workflow;
use Utils::Hash qw|merge_conf|;


use lib "$Bin/../lib/";
use Shared qw|read_list split_chrs|;
use Wrapper qw|bam_sampids|;

############################
## Command line interface ##
############################

my @spec =  (
	Param("conf|c")->valid(sub { -r }),
	Param("list|l")->valid(sub { -r }),
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
	This is a pipeline script to call variants from BAM/CRAMs and create per-individual gVCFs.

Usage:
	gvcf_call.pl --conf Conf --list BAMList --outdir OutDir 

Options:
	--list: BAM/CRAM file list or directory. The list should have 1 or 2 columns per line: path to 
			BAM/CRAM file (required) and sample ID. If only one columns is available, sample ID will 
			be taken from the basename after removing suffix. If a directory is provided, we will 
			look for BAM/CRAM files from the directory. Only one file type is allowed.  
	--wgs:  Under WGS mode, variants will be called across the genome instead of within the targeted 
			region. We will split whole genome into chunks containing different chromosomes.
	--rename: Sample rename list.
	--remove: A list of bad samples to be removed.

Notes:
	The script currently support variant calling using GATK (v4, HaplotypeCaller), weCall, and 
	DeepVariant. To accomondate rich options in each variant caller. It's possible to apply different 
	tools and options to generate gVCF from individual BAM/CRAM files.

	In the config file, each section whose name starts with GATK/WECALL/DEEPVARIANT represents a 
	combination of a caller and its options. All combinations of tools and options will be applied
	to each individual BAM/CRAMs. Each combination will be associated with an output directory with 
	the same name as section header (in lower cases) to store per-individual gVCF files. gVCF output 
	for each individual will have sample ID as basename and appear in the VCF header.

	Per-individual gVCF files are used for joint genotyping across multiple samples by GATK (gvcf_gt.pl)
	or GLnexus.

Dependencies:
	gatk, weCall, DeepVariant (docker), tabix.

EOF
	exit 1;
}

$opt->validate({requires => [qw|conf list outdir|]});

my $rootdir = $opt->get_outdir;

my %conf = merge_conf($opt->get_conf, $opt->get_param); 


############################################
## Input files validation and collection  ##
############################################
my ($filetype, %bams, @samps, @rename);

my %opt = exists $conf{AWS} && exists $conf{AWS}{STAGEIN} ? (notest => 1) : ();
if ($opt->get_list) {
	my ($bams, $type) = read_list($opt->get_list,
		{ suffix => ['bam', 'cram'], 
		  rename => $opt->get_rename, remove => $opt->get_remove, %opt });
	$filetype = $type;
	%bams = %$bams;
	while(my ($iid, $bamfile) = each %bams) {
		# Note: htslib support s3 file, but for some reasons I cannot get access on head node
		# as a workaround, I can login to computing note to create script from dry run 
		# then run wkf_run on head node to 
		my @iids = bam_sampids($bamfile);
		croak "Bam file $bamfile contains multiple samples!" unless @iids == 1;
		if ($iids[0] ne $iid) {
			push @rename, [$iids[0], $iid];
		}
	}
}
else {
	croak "Must provide a BAM/CRAM file list or specify a directory containing bam or cram files."
}
@samps = sort keys %bams;


#################################
##  Workflow initialization    ##
##  Working directory setup    ##
#################################
my @callers = grep { /^(GATK|WECALL|DEEPVAR)/ } sort keys %conf;
if (@callers > 0) {
	print STDERR "The following callers and their setting can be found in config:\n";
	print STDERR join("\t", @callers), "\n";
}
else {
	print STDERR "No caller was found in config.";
	exit 1;
}

my $wkf = Utils::Workflow->new($rootdir,
	{ dir => [map { lc($_) } @callers], engine => $opt->get_engine, 
	  force =>  $opt->get_force, strict_var => 1 });

if ($opt->get_wgs) {
	foreach my $iid (@samps) {
		make_path "$rootdir/wrk/$iid";
	}
}

# set up parameter files
if ($filetype eq 'cram') {
	for(my $ii = 0; $ii < @samps; $ii ++) {
		my $jj = $ii + 1;
		open my $fout, ">$rootdir/par/IID2CRAM.$jj" or die "Cannot write to IID2CRAM";
		print $fout $samps[$ii], "\t", $bams{$samps[$ii]}, "\n";
		open $fout, ">$rootdir/par/IID2BAM.$jj" or die "Cannot write to IID2BAM";
		print $fout $samps[$ii], "\t", "$rootdir/out/$samps[$ii].bam", "\n";
	}
}
else {
	for(my $ii = 0; $ii < @samps; $ii ++) {
		my $jj = $ii + 1;
		open my $fout, ">$rootdir/par/IID2BAM.$jj" or die "Cannot write to IID2BAM";
		print $fout $samps[$ii], "\t", $bams{$samps[$ii]}, "\n";
	}
}

if (@rename) {
	for(my $ii = 0; $ii < @rename; $ii ++) {
		my $jj = $ii + 1;
		open my $fout, ">$rootdir/par/RENAME.$jj" or die "Cannot write to RENAME";
		print $fout $rename[$ii][0], "\t", $rename[$ii][1], "\n";
	}
}

if ($opt->get_wgs) {
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

write_config %conf, "$rootdir/par/run.conf" unless $opt->get_dryrun;

##############################
## Workflow initialization  ##
##############################
if (grep { $_ =~ /^weCall/i } @callers) {
	if ($filetype eq 'cram') {
		$wkf->add(cram2bam(), { name => "CRAM2BAM", 
				expect => [ get_exp("out", "bam", "bam.bai") ], 
				nslots => scalar(@samps)*$conf{PARAM}{NSPLIT}, 
				step => $conf{PARAM}{NSPLIT} });
	}
}

foreach my $caller (@callers) {
	$caller = ucfirst(lc($caller));
	my $deparg = {};
	if ($caller =~ /^weCall/i && $filetype eq 'cram') {
		$deparg = { depend => "CRAM2BAM", deparray => 1 };
	}
	my @outsuf = qw|g.vcf.gz g.vcf.gz.tbi|;
	if ($caller =~ /^DeepVar/i) {
		push @outsuf, qw|vcf.gz vcf.gz.tbi|;
	}
	if ($opt->get_wgs) {
		my $unit = $conf{PARAM}{NSPLIT};
		$wkf->add(varcall_split($caller), { name => $caller."Split", %$deparg,
				expect => [ get_exp_split(lc($caller), $unit, @outsuf) ],
				nslots => scalar(@samps)*$unit, step => 1  })
			->add(varcall_merge($caller), { name => $caller, depend => $caller."Split", deparray => 1,
				expect => [ get_exp(lc($caller), @outsuf) ], nslots => scalar(@samps)*$unit, step => $unit });
	}
	else {
		$wkf->add(varcall($caller), { name => $caller, %$deparg,
				expect => [ get_exp(lc($caller), @outsuf) ], nslots => scalar(@samps) });
	}
	if (@rename) {
		$wkf->add(rename_gvcf($caller), { name => "Rename".$caller, depend => $caller,
			expect => [ get_exp_rename(lc($caller), map { ($_, $_.".bak") } @outsuf) ], 
			nslots => scalar(@rename) });
		if (exists $conf{AWS}{STAGEOUT}) {
			$wkf->add(stage_out($caller), { name => "StageOut".$caller, depend => "Rename".$caller,
				expect => [ get_exp(lc($caller), "stagedout") ], nslots => scalar(@samps) });
		}
	}
	else {
		if (exists $conf{AWS}{STAGEOUT}) {
			$wkf->add(stage_out($caller), { name => "StageOut".$caller, depend => $caller, deparray => 1,
				expect => [ get_exp(lc($caller), "stagedout") ], nslots => scalar(@samps) });
		}
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

# If targeted regions are provided
sub varcall {
	my ($caller) = @_;
	my $stagein = <<'EOF';

TMPDIR=$(stagein.pl --input $BAM --output _TMPDIR_/BAM._INDEX_ --profile _AWS.PROFILEIN_ --tmpdir _AWS.STAGEIN_ --all)

trap "echo Clean up $TMPDIR; rm -fR $TMPDIR; exit 1" EXIT SIGINT SIGTERM

read BAM < _TMPDIR_/BAM._INDEX_

EOF
	my $script;
	if ($caller =~ /^gatk/i) {
		if ($filetype eq 'cram') {
			$script = "\nread IID BAM < _PARDIR_/IID2CRAM._INDEX_\n";
		}
		else {
			$script = "\nread IID BAM < _PARDIR_/IID2BAM._INDEX_\n";
		}
		if (exists $conf{AWS}{STAGEIN}) {
			$script .= $stagein;
		}
		$script .=<<'EOF';

gatk --java-options "_SOFTWARE.JAVAOPT_" \
	HaplotypeCaller \
	-R _PATH.FASTA_ \
	-I $BAM \
	-O _SOFTWAREDIR_/$IID.g.vcf.gz \
	-L _PATH.TARGETS_ \
	-ERC GVCF \
	_SOFTWARE.OPTION[ \
	]_

EOF
	}
	elsif ($caller =~ /^wecall/i) {
		$script = "\nread IID BAM < _PARDIR_/IID2BAM._INDEX_\n";
		if (exists $conf{AWS}{STAGEIN} && $filetype eq 'bam') {
			$script .= $stagein;
		}
		$script .= <<'EOF';

# Cleanup the temp working directory
mkdir -p _TMPDIR_/$IID/software
rm -fR _TMPDIR_/$IID/software/*

weCall --inputs $BAM \
	--refFile _PATH.FASTA_ \
	--regions _PATH.TARGETS_ \
	--output _SOFTWAREDIR_/$IID.g.vcf \
	--outputRefCalls 1 \
	--logFilename _TMPDIR_/$IID/software.log \
	--workDir _TMPDIR_/$IID/software  \
	_SOFTWARE.OPTION[ \
	]_

bgzip _SOFTWAREDIR_/$IID.g.vcf
tabix -f -p vcf _SOFTWAREDIR_/$IID.g.vcf.gz

EOF

	}
	elsif ($caller =~ /^DeepVar/i) {
		if ($filetype eq 'cram') {
			$script = "\nread IID BAM < _PARDIR_/IID2CRAM._INDEX_\n";
		}
		else {
			$script = "\nread IID BAM < _PARDIR_/IID2BAM._INDEX_\n";
		}
		if (exists $conf{AWS}{STAGEIN}) {
			$script .= $stagein;
		}
		my $CALLER = uc($caller);
		if (!defined $conf{$CALLER}{Docker} || $conf{$CALLER}{Docker} =~ /^Y|T/) {
			#print STDERR Dumper $caller, \%conf;
			$script .=<<'EOF';

# Run deepvariant using Docker
FADIR=$(dirname _PATH.FASTA_)
TGDIR=$(dirname _PATH.TARGETS_)
MODIR=$(dirname _SOFTWARE.MODEL_)
BAMDIR=$(dirname $BAM)
ROOTDIR=$(dirname _WRKDIR_)

docker run --rm \
    -v $FADIR:$FADIR -v $TGDIR:$TGDIR -v $BAMDIR:$BAMDIR -v $ROOTDIR:$ROOTDIR \
	gcr.io/deepvariant-docker/deepvariant:_SOFTWARE.VERSION_ \
	/opt/deepvariant/bin/make_examples \
	--mode calling \
	--ref _PATH.FASTA_ \
	--reads $BAM --use_ref_for_cram=true \
	--regions _PATH.TARGETS_ \
	--examples _SOFTWAREDIR_/$IID.examples.tfrecords.gz \
	--gvcf _SOFTWAREDIR_/$IID.gvcf.tfrecords.gz

docker run --rm \
    -v $MODIR:$MODIR -v $ROOTDIR:$ROOTDIR \
	gcr.io/deepvariant-docker/deepvariant:_SOFTWARE.VERSION_ \
	/opt/deepvariant/bin/call_variants \
	--outfile _SOFTWAREDIR_/$IID.variants.tfrecords.gz \
	--examples _SOFTWAREDIR_/$IID.examples.tfrecords.gz \
	--checkpoint _SOFTWARE.MODEL_

docker run --rm \
    -v $FADIR:$FADIR -v $ROOTDIR:$ROOTDIR \
	gcr.io/deepvariant-docker/deepvariant:_SOFTWARE.VERSION_ \
	/opt/deepvariant/bin/postprocess_variants \
	--ref _PATH.FASTA_ \
	--infile _SOFTWAREDIR_/$IID.variants.tfrecords.gz \
	--outfile _SOFTWAREDIR_/$IID.vcf \
	--nonvariant_site_tfrecord_path _SOFTWAREDIR_/$IID.gvcf.tfrecords.gz \
	--gvcf_outfile _SOFTWAREDIR_/$IID.g.vcf


bgzip _SOFTWAREDIR_/$IID.vcf
tabix -f -p vcf _SOFTWAREDIR_/$IID.vcf.gz

bgzip _SOFTWAREDIR_/$IID.g.vcf
tabix -f -p vcf _SOFTWAREDIR_/$IID.g.vcf.gz

EOF
		}
		else {
			if (exists $conf{CONDA}) {
			$script .= <<'EOF';
source _CONDA.INIT_
conda activate _CONDA.ENV_
EOF
			}
			$script .= <<'EOF';

python _SOFTWARE.BINPATH_/make_examples.zip \
	--mode calling \
	--ref _PATH.FASTA_ \
	--reads $BAM --use_ref_for_cram=true \
	--regions _PATH.TARGETS_ \
	--examples _SOFTWAREDIR_/$IID.examples.tfrecords.gz \
	--gvcf _SOFTWAREDIR_/$IID.gvcf.tfrecords.gz

python _SOFTWARE.BINPATH_/call_variants.zip \
	--outfile _SOFTWAREDIR_/$IID.variants.tfrecords.gz \
	--examples _SOFTWAREDIR_/$IID.examples.tfrecords.gz \
	--checkpoint _SOFTWARE.MODEL_

python _SOFTWARE.BINPATH_/postprocess_variants.zip \
	--ref _PATH.FASTA_ \
	--infile _SOFTWAREDIR_/$IID.variants.tfrecords.gz \
	--outfile _SOFTWAREDIR_/$IID.vcf \
	--nonvariant_site_tfrecord_path _SOFTWAREDIR_/$IID.gvcf.tfrecords.gz \
	--gvcf_outfile _SOFTWAREDIR_/$IID.g.vcf


bgzip _SOFTWAREDIR_/$IID.vcf
tabix -f -p vcf _SOFTWAREDIR_/$IID.vcf.gz

bgzip _SOFTWAREDIR_/$IID.g.vcf
tabix -f -p vcf _SOFTWAREDIR_/$IID.g.vcf.gz

EOF
		}
	}
	else {
		croak "Caller $caller is not supported in gvcf_call.pl";
	}
	$script =~ s/software/${(\ lc($caller) )}/g;
	$script =~ s/SOFTWARE/${(\ uc($caller) )}/g;
	return $script;
}

sub get_exp {
	my $dir = shift @_;
	my @allexps;
	foreach my $iid (@samps) {
		my @exp;
		foreach my $suffix (@_) {
			push @exp, "$dir/$iid.$suffix";
		}
		push @allexps, \@exp;
	}
	return @allexps;
}

sub stage_out {
	my ($caller) = @_;
	my $script = <<'EOF';
read IID BAM < _PARDIR_/IID2CRAM._INDEX_

stageout.pl --files $IID.g.vcf.gz,$IID.g.vcf.gz.tbi --profile _AWS.PROFILEOUT_ --indir _SOFTWAREDIR_ --outdir _AWS.STAGEOUT_/software

EOF
	if ($caller =~ /^DeepVar/i) {
		$script .= <<'EOF'
stageout.pl --files $IID.vcf.gz,$IID.vcf.gz.tbi --profile _AWS.PROFILEOUT_ --indir _SOFTWAREDIR_ --outdir _AWS.STAGEOUT_/software

touch _SOFTWAREDIR_/$IID.stagedout
EOF
	}
	else {
		$script .= "\ntouch _SOFTWAREDIR_/\$IID.stagedout\n";
	}
	$script =~ s/software/${(\ lc($caller) )}/g;
	$script =~ s/SOFTWARE/${(\ uc($caller) )}/g;
	return $script;
}


sub cram2bam {
	my $script = <<'EOF';

II=$(( (_INDEX_-1)/_PARAM.NSPLIT_+1 ))
read IID CRAM < _PARDIR_/IID2CRAM.$II

EOF
	if (exists $conf{AWS}{STAGEIN}) {
		$script .= <<'EOF';

TMPDIR=$(stagein.pl --input $CRAM --output _TMPDIR_/CRAM._INDEX_ --profile _AWS.PROFILEIN_ --tmpdir _AWS.STAGEIN_ --all)

trap "echo Clean up $TMPDIR; rm -fR $TMPDIR; exit 1" EXIT SIGINT SIGTERM

read CRAM < _TMPDIR_/CRAM._INDEX_

EOF
	}
	if (!$opt->get_wgs && defined $conf{PATH}{TARGETS}) {
		if (exists $conf{CRAM2BAM}{OPTION}) {
			$script .= <<'EOF';

gatk --java-options "_PICARD.GATKOPT_" PrintReads \
	-I $CRAM  -L _PATH.TARGETS_ \
	-R _PATH.FASTA_  -O _OUTDIR_/$IID.bam \
	_CRAM2BAM.OPTION[ \
	]_

sleep 1

mv _OUTDIR_/$IID.bai _OUTDIR_/$IID.bam.bai

EOF
		}
		else {
			$script .= <<'EOF';

gatk --java-options "_PICARD.GATKOPT_" PrintReads \
	-I $CRAM  -L _PATH.TARGETS_ \
	-R _PATH.FASTA_  -O _OUTDIR_/$IID.bam

sleep 1

mv _OUTDIR_/$IID.bai _OUTDIR_/$IID.bam.bai

EOF
		}	
	}
	else {
		$script .= << 'EOF';

samtools view -T _PATH.FASTA_ -b -h \
	-o _OUTDIR_/$IID.bam $CRAM 

samtools index _OUTDIR_/$IID.bam

EOF
	}
	return $script;
}

sub rename_gvcf {
	my ($caller) = @_;
	my $script =<<'EOF';

read OLDID NEWID < _PARDIR_/RENAME._INDEX_

mv _SOFTWAREDIR_/$NEWID.g.vcf.gz _SOFTWAREDIR_/$NEWID.g.vcf.gz.bak

mv _SOFTWAREDIR_/$NEWID.g.vcf.gz.tbi _SOFTWAREDIR_/$NEWID.g.vcf.gz.tbi.bak

bcftools reheader -s _PARDIR_/RENAME._INDEX_ \
	-o _SOFTWAREDIR_/$NEWID.g.vcf.gz _SOFTWAREDIR_/$NEWID.g.vcf.gz.bak

tabix -f -p vcf _SOFTWAREDIR_/$NEWID.g.vcf.gz

EOF
	if ($caller =~ /^DeepVar/i) {
		$script .= <<'EOF';

mv _SOFTWAREDIR_/$NEWID.vcf.gz _SOFTWAREDIR_/$NEWID.vcf.gz.bak

mv _SOFTWAREDIR_/$NEWID.vcf.gz.tbi _SOFTWAREDIR_/$NEWID.vcf.gz.tbi.bak

bcftools reheader -s _PARDIR_/RENAME._INDEX_ \
	-o _SOFTWAREDIR_/$NEWID.vcf.gz _SOFTWAREDIR_/$NEWID.vcf.gz.bak

tabix -f -p vcf _SOFTWAREDIR_/$NEWID.vcf.gz

EOF
	}
	my $CALLER = uc($caller);
	$script =~ s/SOFTWARE/$CALLER/g;
	return $script;
}

sub get_exp_rename {
	my $dir = shift @_;
	my @allexps;
	foreach my $iid (map { $_->[1]  } @rename) {
		my @exp;
		foreach my $suffix (@_) {
			push @exp, "$dir/$iid.$suffix";
		}
		push @allexps, \@exp;
	}
	return @allexps;
}

# For WGS case, we split chromosomes based on rough length
sub varcall_split {
	my ($caller) = @_;
	my $script = <<'EOF';

II=$(( (_INDEX_-1)/_PARAM.NSPLIT_+1 ))
SPLIT=$(( (_INDEX_-1)%_PARAM.NSPLIT_+1 ))

EOF
	if ($caller =~ /^gatk/i) {
		if (exists $conf{AWS}{STAGEIN}) {
			warn "Stage-in is not supported for WGS calling by GATK!";
		}
		if ($filetype eq 'cram') {
			$script .= "\nread IID BAM < _PARDIR_/IID2CRAM.\$II\n";
		}
		else {
			$script .= "\nread IID BAM < _PARDIR_/IID2BAM.\$II\n";
		}
		$script .=<< 'EOF';

INTERVAL=`cat _PARDIR_/SEQGRP.$SPLIT | sed -e 's/,/ -L /g'`

gatk --java-options "_SOFTWARE.JAVAOPT_" \
	HaplotypeCaller \
	-R _PATH.FASTA_ \
	-I $BAM \
	-O _WRKDIR_/$IID/software.$SPLIT.g.vcf.gz \
	-L $INTERVAL \
	-ERC GVCF \
	_SOFTWARE.OPTION[ \
	]_

EOF
	}
	elsif ($caller =~ /^wecall/i) {
		if (exists $conf{AWS}{STAGEIN} && $filetype eq 'bam') {
			warn "Stage-in is not supported for WGS calling by weCall!";
		}
		$script .= << 'EOF';

read IID BAM < _PARDIR_/IID2BAM.$II

REGION=`cat _PARDIR_/SEQGRP.$SPLIT`

mkdir -p _TMPDIR_/$IID.$SPLIT
rm -fR _TMPDIR_/$IID.$SPLIT/*

weCall --inputs $BAM \
	--refFile _PATH.FASTA_ \
	--regions $REGION \
	--output _WRKDIR_/$IID/software.$SPLIT.g.vcf \
	--outputRefCalls 1 \
	--logFilename _WRKDIR_/$IID/software.$SPLIT.log \
	--workDir _TMPDIR_/$IID.$SPLIT  \
	_SOFTWARE.OPTION[ \
	]_

bgzip _WRKDIR_/$IID/software.$SPLIT.g.vcf
tabix -f -p _WRKDIR_/$IID/software.$SPLIT.g.vcf

EOF

	}
	elsif ($caller =~ /^DeepVar/i) {
		if ($filetype eq 'cram') {
			$script .= "\nread IID BAM < _PARDIR_/IID2CRAM.\$II\n";
		}
		else {
			$script .= "\nread IID BAM < _PARDIR_/IID2BAM.\$II\n";
		}
		if (exists $conf{AWS}{STAGEIN}) {
			warn "Stage-in is not supported for WGS calling by DeepVariant!";
		}
		my $CALLER = uc($caller);
		if (!defined $conf{$CALLER}{Docker} || $conf{$CALLER}{Docker} =~ /^Y|T/) {
			$script .= <<'EOF';

REGIONS=`cat _PARDIR_/SEQGRP.$SPLIT | sed -e 's/,/ /g'`

# Run deepvariant using Docker
FADIR=$(dirname _PATH.FASTA_)
TGDIR=$(dirname _PATH.TARGETS_)
MODIR=$(dirname _SOFTWARE.MODEL_)
BAMDIR=$(dirname $BAM)
ROOTDIR=$(dirname _WRKDIR_)

docker run --rm \
    -v $FADIR:$FADIR -v $TGDIR:$TGDIR -v $BAMDIR:$BAMDIR -v $ROOTDIR:$ROOTDIR \
	gcr.io/deepvariant-docker/deepvariant:_SOFTWARE.VERSION_ \
	/opt/deepvariant/bin/make_examples \
	--mode calling \
	--ref _PATH.FASTA_ \
	--reads $BAM --use_ref_for_cram=true \
	--regions "$REGIONS" \
	--examples _WRKDIR_/$IID/software.$SPLIT.examples.tfrecords.gz \
	--gvcf _WRKDIR_/$IID/software.$SPLIT.gvcf.tfrecords.gz

docker run --rm \
    -v $MODIR:$MODIR -v $ROOTDIR:$ROOTDIR \
	gcr.io/deepvariant-docker/deepvariant:_SOFTWARE.VERSION_ \
	/opt/deepvariant/bin/call_variants \
	--outfile _WRKDIR_/$IID/software.$SPLIT.variants.tfrecords.gz \
	--examples _WRKDIR_/$IID/software.$SPLIT.examples.tfrecords.gz \
	--checkpoint _SOFTWARE.MODEL_

docker run --rm \
    -v $FADIR:$FADIR -v $ROOTDIR:$ROOTDIR \
	gcr.io/deepvariant-docker/deepvariant:_SOFTWARE.VERSION_ \
	/opt/deepvariant/bin/postprocess_variants \
	--ref _PATH.FASTA_ \
	--infile _WRKDIR_/$IID/software.$SPLIT.variants.tfrecords.gz \
	--outfile _WRKDIR_/$IID/software.$SPLIT.vcf \
	--nonvariant_site_tfrecord_path _WRKDIR_/$IID/software.$SPLIT.gvcf.tfrecords.gz \
	--gvcf_outfile _WRKDIR_/$IID/software.$SPLIT.g.vcf


bgzip _WRKDIR_/$IID/software.$SPLIT.vcf
tabix -f -p vcf _WRKDIR_/$IID/software.$SPLIT.vcf.gz

bgzip _WRKDIR_/$IID/software.$SPLIT.g.vcf
tabix -f -p vcf _WRKDIR_/$IID/software.$SPLIT.g.vcf.gz

EOF
		}
		else {
			if (exists $conf{CONDA}) {
			$script .= <<'EOF';
source _CONDA.INIT_
conda activate _CONDA.ENV_
EOF
			}
			$script .= <<'EOF';

python _SOFTWARE.BINPATH_/make_examples.zip \
	--mode calling \
	--ref _PATH.FASTA_ \
	--reads $BAM --use_ref_for_cram=true \
	--regions "$REGIONS" \
	--examples _WRKDIR_/$IID/software.$SPLIT.examples.tfrecords.gz \
	--gvcf _WRKDIR_/$IID/software.$SPLIT.gvcf.tfrecords.gz

python _SOFTWARE.BINPATH_/call_variants.zip \
	--outfile _WRKDIR_/$IID/software.$SPLIT.variants.tfrecords.gz \
	--examples _WRKDIR_/$IID/software.$SPLIT.examples.tfrecords.gz \
	--checkpoint _SOFTWARE.MODEL_

python _SOFTWARE.BINPATH_/postprocess_variants.zip \
	--ref _PATH.FASTA_ \
	--infile _WRKDIR_/$IID/software.$SPLIT.variants.tfrecords.gz \
	--outfile _WRKDIR_/$IID/software.$SPLIT.vcf \
	--nonvariant_site_tfrecord_path _WRKDIR_/$IID/software.$SPLIT.gvcf.tfrecords.gz \
	--gvcf_outfile _WRKDIR_/$IID/software.$SPLIT.g.vcf

bgzip _WRKDIR_/$IID/software.$SPLIT.vcf
tabix -f -p vcf _WRKDIR_/$IID/software.$SPLIT.vcf.gz

bgzip _WRKDIR_/$IID/software.$SPLIT.g.vcf
tabix -f -p vcf _WRKDIR_/$IID/software.$SPLIT.g.vcf.gz

EOF
		}
	}
	else {
		croak "Caller $caller is not supported";
	}
	$script =~ s/software/${(\ lc($caller) )}/g;
	$script =~ s/SOFTWARE/${(\ uc($caller) )}/g;
	return $script;
}

sub get_exp_split {
	my $caller = shift @_;
	my $step = shift @_;
	my @allexps;
	foreach my $iid (@samps) {
		for(my $jj = 1; $jj <= $step; $jj ++) {
			my @exp;
			foreach my $suffix (@_) {
				push @exp, "wrk/$iid/$caller.$jj.$suffix";
			}
			push @allexps, \@exp;
		}
	}
	return @allexps;
}


sub varcall_merge {
	my ($caller) = @_;
	my $script =<<'EOF';

II=$(( (_INDEX_-1)/_PARAM.NSPLIT_+1 ))

read IID BAM < _PARDIR_/IID2BAM.$II

INPUT=`perl -e 'print join(" ", map {"INPUT=_WRKDIR_/$ARGV[1]/software.$_.g.vcf.gz" } 1..$ARGV[0])' _PARAM.NSPLIT_ $IID`

picard _PICARD.JAVAOPT_ MergeVcfs \
	$INPUT \
	OUTPUT=_SOFTWAREDIR_/$IID.g.vcf.gz CREATE_INDEX=true

EOF
	if ($caller =~ /^deepvar/i) {
		$script .= <<'EOF'

INPUT=`perl -e 'print join(" ", map {"INPUT=_WRKDIR_/$ARGV[1]/software.$_.vcf.gz" } 1..$ARGV[0])' _PARAM.NSPLIT_ $IID`

picard _PICARD.JAVAOPT_ MergeVcfs \
	$INPUT \
	OUTPUT=_SOFTWAREDIR_/$IID.vcf.gz CREATE_INDEX=true

EOF
	}
	$script =~ s/software/${(\ lc($caller) )}/g;
	$script =~ s/SOFTWARE/${(\ uc($caller) )}/g;
	return $script;
}


