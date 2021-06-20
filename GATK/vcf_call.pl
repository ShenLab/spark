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
use Hash::Util qw|lock_hash_recurse|;
use Utils::Workflow;
use Utils::Hash qw|merge_conf|;


use lib "$Bin/../lib/";
use Shared qw|read_list read_dir split_chrs|;
use Wrapper qw|bam_sampids|;

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
	This is a pipeline script to call variants from BAM/CRAMs for each individual or group of 
	individuals and create VCFs.

Usage:
	vcf_call.pl --conf Config --list BAM.list --outdir RootDir 

Options:
	--list: BAM/CRAM file list or directory. The list should have 1 or 2 columns per line: path to 
			BAM/CRAM file (required) and sample ID. If only one column is available, sample ID will 
			after removbe taken from file basename ing suffix. If a directory is provided, we will 
			look for BAM/CRAM files from the directory. Only one file type is allowed.  
	--wgs:  Under WGS mode, variants will be called across the genome instead of focusing on the targeted 
			region. 
	--group: Sample groups. It should have 2 columns per line, IID and GroupID. IID must be found in 
			 BAM file list. Each sample is allowed to exist in multiple groups. This option is only 
			 available for callers that supports multi-sample joint calling directly from BAM/CRAMs. 
			 A multi-sample VCF will be generated for each group, and samples without GroupID will not
			 be used.
	--rename: Sample rename list.
	--remove: A list of bad samples to be removed.
	
Notes:
	The script currently supports the following variant callers: GATK (v4, HapotypeCaller and MuTect2),
	weCall, bcftools (v1, multi-allelic caller), and DeepVariant. Differnent tools and options can be
	used to generate single or multiple sample VCFs from BAM/CRAMs.

	In the config file, each section whose name starts with GATK/MUTECT/WECALL/BCFTOOLS/DEEPVARIANT 
	represents a combination of a caller and its options. All combinations of tools and options will 
	be applied to each individual BAM/CRAMs. Each combination will be associated with an output directory 
	with the same name as section header (in lower cases) to store VCF files. Per-individual VCF will 
	have sample ID as basename; per-group VCF will have group ID as basename.

Dependencies:
	gatk, weCall, bcftools, DeepVariant (docker), tabix.

EOF
	exit 1;
}

$opt->validate({requires => [qw|conf list outdir|]});

my $rootdir = $opt->get_outdir;

my %conf = merge_conf($opt->get_conf, $opt->get_param); 

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

lock_hash_recurse(%conf);

############################################
## Input files validation and collection  ##
############################################
my ($filetype, %bams, @samps, %rename);

my %opt = exists $conf{AWS} && exists $conf{AWS}{STAGEIN} ? (notest => 1) : ();
if ($opt->get_list) {
	my $list = $opt->get_list;
	unless(-f $list || -d $list) {
		die "The BAM/CRAM list $list does not exist!"
	}
	if (-d $list) {
		$list = abs_path($list);
	}
	my ($bams, $type) = read_list($list, { suffix => ['bam', 'cram'], 
		  rename => $opt->get_rename, remove => $opt->get_remove, %opt });
	$filetype = $type;
	%bams = %$bams;
	while(my ($iid, $bamfile) = each %bams) {
		my @iids = bam_sampids($bamfile);
		croak "Bam file $bamfile contains multiple samples!" unless @iids == 1;
		if ($iids[0] ne $iid) {
			$rename{$iid} = $iids[0];
			#push @rename, [$iids[0], $iid];
		}
	}
}
else {
	croak "Must provide a BAM/CRAM file list or specify a directory of BAM/CRAM files."
}

my (%group, @grps, @needrn);

# when group file is provided, only samples in the group file will be processed
if ($opt->get_group) {
	open my $fin, $opt->get_group or die "Cannot open group file";
	while(<$fin>) {
		next if /^\s*$/;
		my ($iid, $gid) = split;
		unless (defined $bams{$iid}) {
			warn "Cannot find BAM file for $iid in group $gid";
		}
		else {
			push @{$group{$gid}} => $iid;
		}
	}
	@grps = sort keys %group;
	@samps = sort map { @$_ } values %group;
	foreach my $gid (@grps) {
		if (any { defined $rename{$_} } @{$group{$gid}}) {
			push @needrn, $gid;
		}
	}
}
else {
	# If no group files is give, groups will be the same as samples
	@samps = sort keys %bams;
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
my @callers = grep { /^(GATK|WECALL|BCFTOOLS|DEEPVAR|MUTECT)/ } sort keys %conf;
if (@callers > 0) {
	print STDERR "The following callers and their setting can be found in config:\n";
	print STDERR join("\t", @callers), "\n";
	if (grep { /^DEEPVAR/ } @callers) {
		if ($opt->get_group) {
			# Group info will be ignored for callers that do not support mult-sample calling
			warn "DeepVariant does not support multi-sample calling\n";
		}
	}
}
else {
	print STDERR "No caller was found in config!";
	exit 1;
}

# If whatshap is specified, turn on hapflag.
my $hapflag;
if (grep { $_ eq 'WHATSHAP' } keys %conf) {
	$hapflag = 1;
}

my $wkf = Utils::Workflow->new($rootdir,
	{ dir => [map { lc($_) } @callers], engine => $opt->get_engine, force => $opt->get_force, strict_var => 1 });

if ($opt->get_wgs) {
	#foreach my $iid (@samps) {
	#	make_path "$rootdir/wrk/$iid";
	#}
	foreach my $gid (@grps) {
		make_path "$rootdir/wrk/$gid";
	}
}

# Set up parameter files
# IID2BAM: BAM or CRAM for individuals, use original path 
# GRPID2CRAM: CRAM for grouped individuals, use original path
# GRPID2BAM: BAM for grouped individuals (if original data is CRAM format, BAM is converted from them)
if ($filetype eq 'cram') {
	for(my $ii = 0; $ii < @samps; $ii ++) {
		my $jj = $ii + 1;
		open my $fout, ">$rootdir/par/IID2BAM.$jj" or die "Cannot write to IID2BAM";
		print $fout $samps[$ii], "\t", $bams{$samps[$ii]}, "\n";
		#open $fout, ">$rootdir/par/IID2BAM.$jj" or die "Cannot write to IID2BAM";
		#print $fout $samps[$ii], "\t", "$rootdir/out/$samps[$ii].bam", "\n";
	}
	for(my $ii = 0; $ii < @grps; $ii ++) {
		my $jj = $ii + 1;
		open my $fout, ">$rootdir/par/GRPID2CRAM.$jj" or die "Cannot write to GRPID2CRAM";
		print $fout $grps[$ii], "\t", join(",", map { $bams{$_} } @{$group{$grps[$ii]}}), "\n";
		open $fout, ">$rootdir/par/GRPID2BAM.$jj" or die "Cannot write to GID2BAM";
		print $fout $grps[$ii], "\t", join(",", map { "$rootdir/out/$_.bam" } @{$group{$grps[$ii]}}), "\n";
	}
}
else {
	for(my $ii = 0; $ii < @samps; $ii ++) {
		my $jj = $ii + 1;
		open my $fout, ">$rootdir/par/IID2BAM.$jj" or die "Cannot write to IID2BAM";
		print $fout $samps[$ii], "\t", $bams{$samps[$ii]}, "\n";
	}
	for(my $ii = 0; $ii < @grps; $ii ++) {
		my $jj = $ii + 1;
		open my $fout, ">$rootdir/par/GRPID2BAM.$jj" or die "Cannot write to GRPID2BAM";
		print $fout $grps[$ii], "\t", join(",", map { $bams{$_} } @{$group{$grps[$ii]}}), "\n";
	}
}

if (@needrn) {
	for(my $ii = 0; $ii < @needrn; $ii ++) {
		my $jj = $ii + 1;	
		open my $frn, ">$rootdir/par/RENAME.$jj" or die "Cannot write to RENAME.$jj";
		if ($opt->get_group) {
			open my $fout, ">$rootdir/par/NEEDRN.$jj" or die "Cannot write to NEEDRN.$jj";
			print $fout $needrn[$ii], "\n";
			foreach my $iid (@{$group{$needrn[$ii]}}) {
				if (defined $rename{$iid}) {
					print $frn $rename{$iid}, "\t", $iid, "\n";
				}
			}
		}
		else {
			my $origid = $rename{$needrn[$ii]};
			die "Cannot find original ID for $needrn[$ii]" unless defined $origid;
			print $frn $origid, "\t", $needrn[$ii], "\n";
		}
	}
}

write_config %conf, "$rootdir/par/run.conf" unless $opt->get_dryrun;

##############################
## Workflow initialization  ##
##############################

# For software that support BAM input only, CRAM will be converted to BAM and bamflag will be on
# Note: bamflag seems legacy, can be ignored.
my $bamflag;
if (grep { $_ =~ /^weCall/i } @callers) {
	if ($filetype eq 'cram') {
		$wkf->add(cram2bam(), { name => "CRAM2BAM", 
				expect => [ get_ind_exp("out", "bam", "bam.bai") ], 
				nslots => scalar(@samps)*$conf{PARAM}{NSPLIT},
				step => $conf{PARAM}{NSPLIT} });
		unless($opt->get_wgs) {
			$bamflag = 1;
		}
	}
}

foreach my $caller (@callers) {
	$caller = ucfirst(lc($caller));
	my $deparg = {};
	if ($caller =~ /^weCall/i && $filetype eq 'cram') {
		if (scalar(@samps) == scalar(@grps)) {
			$deparg = { depend => "CRAM2BAM", deparray => 1 };
		}
		else {
			$deparg = { depend => "CRAM2BAM" };
		}	
	}
	my $nsamps;
	if ($caller =~ /^deepvar/i) {
		$nsamps = scalar(@samps);
	}
	else {
		$nsamps = scalar(@grps);
	}

	my @suffix = qw(vcf.gz vcf.gz.tbi);
	my @mergesuffix = @suffix;
	if ($caller =~ /^mutect/i) {
		push @suffix, "vcf.gz.stats";
		push @mergesuffix, "vcf.gz.stats";
	}
	elsif ($hapflag) {
		# Note: we do not merge unphased VCF
		push @suffix, qw(unphased.vcf.gz unphased.vcf.gz.tbi);
	}

	if ($opt->get_wgs) {
		my $unit = $conf{PARAM}{NSPLIT};
		$wkf->add(varcall_split($caller, $bamflag, $hapflag), { name => $caller."Split", %$deparg,
				expect => [ get_exp_split("wrk", $unit, @suffix) ],
				nslots => $nsamps * $unit, step => 1   })
			->add(varcall_merge($caller), { name => $caller, depend => $caller."Split", deparray => 1,
				expect => [ get_exp(lc($caller), @mergesuffix) ],
				nslots => $nsamps * $unit, step => $unit });
	}
	else {
		$wkf->add(varcall($caller, $bamflag, $hapflag), { name => $caller, %$deparg,
				expect => [ get_exp(lc($caller), @suffix) ],
				nslots => $nsamps });
	}
	if (@needrn) {
		if ($caller =~ /^deepvar/i || !$opt->get_group) {
			$wkf->add(rename_ind_vcf($caller), { name => "Rename".$caller, depend => $caller,
				expect => [ get_exp_indrename(lc($caller), "orig.vcf.gz", "orig.vcf.gz.tbi", "vcf.gz", "vcf.gz.tbi") ], 
				nslots => scalar(@needrn) });
		}
		else {
			$wkf->add(rename_grp_vcf($caller), { name => "Rename".$caller, depend => $caller,
				expect => [ get_exp_grprename(lc($caller), "orig.vcf.gz", "orig.vcf.gz.tbi", "vcf.gz", "vcf.gz.tbi") ], 
				nslots => scalar(@needrn) });
		}
		if (exists $conf{AWS} && exists $conf{AWS}{STAGEOUT}) {
			$wkf->add(stage_out($caller), { name => "StageOut".$caller, depend => "Rename".$caller,
				expect => [ get_exp(lc($caller), "stagedout") ], nslots => $nsamps });
		}
	}
	else {
		if (exists $conf{AWS} && exists $conf{AWS}{STAGEOUT}) {
			$wkf->add(stage_out($caller), { name => "StageOut".$caller, depend => $caller, deparray => 1,
				expect => [ get_exp(lc($caller), "stagedout") ], nslots => $nsamps });
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
	my ($caller, $bamflag, $hapflag) = @_;
	# Stage-in will change the path of BAM/CRAM files
	my $stagein = <<'EOF';

TMPDIR=$(stagein.pl --input $BAM --output _TMPDIR_/BAM._INDEX_ --profile _AWS.PROFILEIN_ --tmpdir _AWS.STAGEIN_ --all)

trap "echo Clean up $TMPDIR; rm -fR $TMPDIR; exit 1" EXIT SIGINT SIGTERM

read BAM < _TMPDIR_/BAM._INDEX_

EOF
	my $script;
	if ($caller =~ /^(gatk|bcftools|mutect)/i) {
		if ($filetype eq 'cram') {
			#if ($bamflag) {
			#	$script = "\nread GRPID BAM < _PARDIR_/GRPID2BAM._INDEX_\n";
			#}
			#else {
				$script = "\nread GRPID BAM < _PARDIR_/GRPID2CRAM._INDEX_\n";
			#}
		}
		else {
			$script = "\nread GRPID BAM < _PARDIR_/GRPID2BAM._INDEX_\n";
		}
		if (exists $conf{AWS} && exists $conf{AWS}{STAGEIN}) {
			$script .= $stagein;
		}
		if ($caller =~ /^gatk/i) {
			$script .=<<'EOF';

INPUTS=$(echo $BAM | sed -e 's/,/ -I /g')

gatk --java-options "_SOFTWARE.JAVAOPT_" \
	HaplotypeCaller \
	-R _PATH.FASTA_ \
	-I $INPUTS \
	-O _SOFTWAREDIR_/$GRPID.vcf.gz \
	-L _PATH.TARGETS_ \
	_SOFTWARE.OPTION[ \
	]_

EOF
		}
		elsif ($caller =~ /^mutect/i) {
			$script .=<<'EOF';

INPUTS=$(echo $BAM | sed -e 's/,/ -I /g')

gatk --java-options "_SOFTWARE.JAVAOPT_" \
	Mutect2 \
	-R _PATH.FASTA_ \
	-I $INPUTS \
	-O _SOFTWAREDIR_/$GRPID.vcf.gz \
	-L _PATH.TARGETS_ \
	_SOFTWARE.OPTION[ \
	]_

EOF
		}
		else {
			$script .=<<'EOF';

INPUTS=$(echo $BAM | sed -e 's/,/ /g')

bcftools mpileup -f _PATH.FASTA_ \
	-R _PATH.TARGETS_ \
	_SOFTWARE.MPILEUPOPT[ \
	]_ 	$INPUTS | \
	bcftools call -mv -Oz \
	_SOFTWARE.CALLOPT[ \
	]_  -o _SOFTWAREDIR_/$GRPID.vcf.gz

tabix -f -p vcf _SOFTWAREDIR_/$GRPID.vcf.gz

EOF
		}
	}
	elsif ($caller =~ /^wecall/i) {
		$script = "\nread GRPID BAM < _PARDIR_/GRPID2BAM._INDEX_\n";
		if (exists $conf{AWS} && exists $conf{AWS}{STAGEIN} && $filetype eq 'bam') {
			$script .= $stagein;
		}
		$script .= <<'EOF';

# Cleanup the temp working directory
mkdir -p _TMPDIR_/$GRPID/software
rm -fR _TMPDIR_/$GRPID/software/*

weCall --inputs $BAM \
	--refFile _PATH.FASTA_ \
	--regions _PATH.TARGETS_ \
	--output _SOFTWAREDIR_/$GRPID.vcf \
	--outputRefCalls 0 \
	--logFilename _TMPDIR_/$GRPID/software.log \
	--workDir _TMPDIR_/$GRPID/software  \
	_SOFTWARE.OPTION[ \
	]_

bgzip _SOFTWAREDIR_/$GRPID.vcf
tabix -f -p vcf _SOFTWAREDIR_/$GRPID.vcf.gz

EOF
	}
	elsif ($caller =~ /^DeepVar/i) {
		# Although DeepVar does not support multi-sample calling
		# we still read in GRPID
		$script = "\nread GRPID BAM < _PARDIR_/IID2BAM._INDEX_\n";
		if (exists $conf{AWS} && exists $conf{AWS}{STAGEIN}) {
			$script .= $stagein;
		}
		my $CALLER = uc($caller);
		if (!defined $conf{$CALLER}{Docker} || $conf{$CALLER}{Docker} =~ /^Y|T/) {
		#if ($conf{$caller}{Docker} =~ /^Y|T/) {
			$script .= <<'EOF';

IID=$GRPID

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
	--reads $BAM \
	--regions _PATH.TARGETS_ \
	--examples _SOFTWAREDIR_/$IID.examples.tfrecords.gz

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
	--outfile _SOFTWAREDIR_/$IID.vcf

bgzip _SOFTWAREDIR_/$IID.vcf
tabix -f -p vcf _SOFTWAREDIR_/$IID.vcf.gz

EOF
		}
		else {
			if (exists $conf{$CALLER}{CONDAENV}) {
			$script .= <<'EOF';
IID=$GRPID

source _PATH.CONDARC_
conda activate _SOFTWARE.CONDAENV_

EOF
			}
			$script .= <<'EOF';

python _SOFTWARE.BINPATH_/make_examples.zip \
	--mode calling \
	--ref _PATH.FASTA_ \
	--reads $BAM \
	--regions _PATH.TARGETS_ \
	--examples _SOFTWAREDIR_/$IID.examples.tfrecords.gz

python _SOFTWARE.BINPATH_/call_variants.zip \
	--outfile _SOFTWAREDIR_/$IID.variants.tfrecords.gz \
	--examples _SOFTWAREDIR_/$IID.examples.tfrecords.gz \
	--checkpoint _SOFTWARE.MODEL_

python _SOFTWARE.BINPATH_/postprocess_variants.zip \
	--ref _PATH.FASTA_ \
	--infile _SOFTWAREDIR_/$IID.variants.tfrecords.gz \
	--outfile _SOFTWAREDIR_/$IID.vcf

bgzip _SOFTWAREDIR_/$IID.vcf
tabix -f -p vcf _SOFTWAREDIR_/$IID.vcf.gz

EOF
		}		
	}
	else {
		croak "Caller $caller is not supported";
	}
	# Adding whatshap for RB phasing
	if ($hapflag) {
		if ($conf{WHATSHAP}{CONDAENV}) {
			$script .= <<'EOF';
source _PATH.CONDARC_
conda activate _WHATSHAP.CONDAENV_

EOF
		}
		if ($caller =~ /^DeepVar/i) {
			$script .= <<'EOF';

mv _SOFTWAREDIR_/$IID.vcf.gz _SOFTWAREDIR_/$IID.unphased.vcf.gz
mv _SOFTWAREDIR_/$IID.vcf.gz.tbi _SOFTWAREDIR_/$IID.unphased.vcf.gz.tbi

whatshap phase -r _PATH.FASTA_ _WHATSHAP.OPTION[ ]_ -o _SOFTWAREDIR_/$IID.vcf \
	_SOFTWAREDIR_/$IID.unphased.vcf.gz $BAM

bgzip _SOFTWAREDIR_/$IID.vcf
tabix -f -p vcf _SOFTWAREDIR_/$IID.vcf.gz

EOF
		}
		else {
			$script .= <<'EOF';

INPUTS2=$(echo $BAM | sed -e 's/,/ /g')

mv _SOFTWAREDIR_/$GRPID.vcf.gz _SOFTWAREDIR_/$GRPID.unphased.vcf.gz
mv _SOFTWAREDIR_/$GRPID.vcf.gz.tbi _SOFTWAREDIR_/$GRPID.unphased.vcf.gz.tbi

whatshap phase -r _PATH.FASTA_ _WHATSHAP.OPTION[ ]_ -o _SOFTWAREDIR_/$GRPID.vcf.gz \
	_SOFTWAREDIR_/$GRPID.unphased.vcf.gz $INPUTS2

bgzip _SOFTWAREDIR_/$GRPID.vcf
tabix -f -p vcf _SOFTWAREDIR_/$GRPID.vcf.gz

EOF
		}
	}
	$script =~ s/software/${(\ lc($caller) )}/g;
	$script =~ s/SOFTWARE/${(\ uc($caller) )}/g;
	return $script;
}

sub get_ind_exp {
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

sub get_exp {
	my $dir = shift @_;
	my @iids;
	if ($dir =~ /^deepvar/i) {
		@iids = @samps;
	}
	else {
		@iids = @grps;
	}

	my @allexps;
	foreach my $grpid (@iids) {
		my @exp;
		foreach my $suffix (@_) {
			push @exp, "$dir/$grpid.$suffix";
		}
		push @allexps, \@exp;
	}
	return @allexps;
}


sub stage_out {
	my ($caller) = @_;
	my $script;
	if ($caller =~ /^DeepVar/i) {
		$script = <<'EOF';

read GRPID BAM < _PARDIR_/IID2BAM._INDEX_

EOF
	}
	else {
		$script = <<'EOF';

read GRPID BAM < _PARDIR_/GRPID2BAM._INDEX_

EOF
	}
	if ($caller =~ /^mutect/i) {
		$script .=<< 'EOF';

stageout.pl --files $GRPID.vcf.gz,$GRPID.vcf.gz.stats,$GRPID.vcf.gz.tbi --profile _AWS.PROFILEOUT_  \
	--indir _SOFTWAREDIR_ --outdir _AWS.STAGEOUT_/software

EOF
	}
	else {
		$script .=<< 'EOF';

stageout.pl --files $GRPID.vcf.gz,$GRPID.vcf.gz.tbi --profile _AWS.PROFILEOUT_ \
	--indir _SOFTWAREDIR_ --outdir _AWS.STAGEOUT_/software

EOF
	}
	$script .= 'touch _SOFTWAREDIR_/$GRPID.stagedout'."\n";
	$script =~ s/software/${(\ lc($caller) )}/g;
	$script =~ s/SOFTWARE/${(\ uc($caller) )}/g;
	return $script;
}


sub cram2bam {
	my $script = <<'EOF';

II=$(( (_INDEX_-1)/_PARAM.NSPLIT_+1 ))
read IID BAM < _PARDIR_/IID2BAM.$II

EOF
	if (exists $conf{AWS} && exists $conf{AWS}{STAGEIN}) {
		$script .= <<'EOF';

TMPDIR=$(stagein.pl --input $BAM --output _TMPDIR_/BAM._INDEX_ --profile _AWS.PROFILEIN_ --tmpdir _AWS.STAGEIN_ --all)

trap "echo Clean up $TMPDIR; rm -fR $TMPDIR; exit 1" EXIT SIGINT SIGTERM

read BAM < _TMPDIR_/BAM._INDEX_

EOF
	}
	if (!$opt->get_wgs && defined $conf{PATH}{TARGETS}) {
		if (exists $conf{CRAM2BAM}{OPTION}) {
			$script .= <<'EOF';

gatk --java-options "_PICARD.GATKOPT_" PrintReads \
	-I $BAM  -L _PATH.TARGETS_ _CRAM2BAM.OPTION[ ]_ \
	-R _PATH.FASTA_  -O _OUTDIR_/$IID.bam

sleep 1

mv _OUTDIR_/$IID.bai _OUTDIR_/$IID.bam.bai

EOF
		}
		else {
			$script .= <<'EOF';

gatk --java-options "_PICARD.GATKOPT_" PrintReads \
	-I $BAM  -L _PATH.TARGETS_ \
	-R _PATH.FASTA_  -O _OUTDIR_/$IID.bam

sleep 1

mv _OUTDIR_/$IID.bai _OUTDIR_/$IID.bam.bai

EOF
		}	
	}
	else {
		$script .= <<'EOF';

samtools view -T _PATH.FASTA_ -b -h \
	-o _OUTDIR_/$IID.bam $BAM 

samtools index _OUTDIR_/$IID.bam

EOF
	}
	return $script;
}

sub rename_ind_vcf {
	my ($caller) = @_;
	my $script =<<'EOF';

read OLDID NEWID < _PARDIR_/RENAME._INDEX_

mv _SOFTWAREDIR_/$NEWID.vcf.gz _SOFTWAREDIR_/$NEWID.orig.vcf.gz
mv _SOFTWAREDIR_/$NEWID.vcf.gz.tbi _SOFTWAREDIR_/$NEWID.orig.vcf.gz.tbi

bcftools reheader -s _PARDIR_/RENAME._INDEX_ \
	-o _SOFTWAREDIR_/$NEWID.vcf.gz _SOFTWAREDIR_/$NEWID.orig.vcf.gz

tabix -f -p vcf _SOFTWAREDIR_/$NEWID.vcf.gz

EOF
	my $CALLER = uc($caller);
	$script =~ s/SOFTWARE/$CALLER/g;
	return $script;
}

sub get_exp_indrename {
	my $dir = shift @_;
	my @allexps;
	#foreach my $iid (map { $_->[1]  } @rename) {
	foreach my $iid (sort keys %rename) {
		my @exp;
		foreach my $suffix (@_) {
			push @exp, "$dir/$iid.$suffix";
		}
		push @allexps, \@exp;
	}
	return @allexps;	
}

sub rename_grp_vcf {
	my ($caller) = @_;
	my $script =<<'EOF';

read GRPID < _PARDIR_/NEEDRN._INDEX_

mv _SOFTWAREDIR_/$GRPID.vcf.gz _SOFTWAREDIR_/$GRPID.orig.vcf.gz
mv _SOFTWAREDIR_/$GRPID.vcf.gz.tbi _SOFTWAREDIR_/$GRPID.orig.vcf.gz.tbi

bcftools reheader -s _PARDIR_/RENAME._INDEX_ \
	-o _SOFTWAREDIR_/$GRPID.vcf.gz _SOFTWAREDIR_/$GRPID.orig.vcf.gz

tabix -f -p vcf _SOFTWAREDIR_/$GRPID.vcf.gz

EOF
	my $CALLER = uc($caller);
	$script =~ s/SOFTWARE/$CALLER/g;
	return $script;
}

sub get_exp_grprename {
	my $dir = shift @_;
	my @allexps;
	foreach my $grpid (@needrn) {
		my @exp;
		foreach my $suffix (@_) {
			push @exp, "$dir/$grpid.$suffix";
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
	if ($caller =~ /^(gatk|bcftools|mutect)/i) {
		if (exists $conf{PATH}{STAGEIN}) {
			warn "Stage-in is not supported for WGS calling by $caller!";
		}
		if ($filetype eq 'cram') {
			$script .= "\nread GRPID BAM < _PARDIR_/GRPID2CRAM.\$II\n";
		}
		else {
			$script .= "\nread GRPID BAM < _PARDIR_/GRPID2BAM.\$II\n";
		}
		if ($caller =~ /^gatk/i) {
			$script .=<< 'EOF';

INTERVAL=`cat _PARDIR_/SEQGRP.$SPLIT | sed -e 's/,/ -L /g'`

INPUTS=$(echo $BAM | sed -e 's/,/ -I /g')

gatk --java-options "_SOFTWARE.JAVAOPT_" \
	HaplotypeCaller \
	-R _PATH.FASTA_ \
	-I $INPUTS \
	-O _WRKDIR_/$GRPID/software.$SPLIT.vcf.gz \
	-L $INTERVAL \
	_SOFTWARE.OPTION[ \
	]_

EOF
		}
		elsif ($caller =~ /^mutect/i) {
			$script .=<<'EOF';

INTERVAL=`cat _PARDIR_/SEQGRP.$SPLIT | sed -e 's/,/ -L /g'`

INPUTS=$(echo $BAM | sed -e 's/,/ -I /g')

gatk --java-options "_SOFTWARE.JAVAOPT_" \
	Mutect2 \
	-R _PATH.FASTA_ \
	-I $INPUTS \
	-O _WRKDIR_/$GRPID/software.$SPLIT.vcf.gz \
	-L $INTERVAL \
	_SOFTWARE.OPTION[ \
	]_

EOF
		}
		else {
			$script .=<<'EOF';

INTERVAL=`cat _PARDIR_/SEQGRP.$SPLIT`

INPUTS=$(echo $BAM | sed -e 's/,/ /g')

bcftools mpileup -f _PATH.FASTA_ \
	-r $INTERVAL \
	_SOFTWARE.MPILEUPOPT[ \
	]_ 	$INPUTS | \
	bcftools call -mv -Oz \
	_SOFTWARE.CALLOPT[ \
	]_  -o _WRKDIR_/$GRPID/software.$SPLIT.vcf.gz

tabix -f -p vcf _SOFTWAREDIR_/$GRPID.vcf.gz

EOF
		}
	}
	elsif ($caller =~ /^wecall/i) {
		if (exists $conf{PATH}{STAGEIN} && $filetype eq 'bam') {
			warn "Stage-in is not supported for WGS calling by weCall!";
		}
		$script .= << 'EOF';

read GRPID BAM < _PARDIR_/GRPID2BAM.$II

REGION=`cat _PARDIR_/SEQGRP.$SPLIT`

mkdir -p _TMPDIR_/$GRPID.$SPLIT
rm -fR _TMPDIR_/$GRPID.$SPLIT/*

weCall --inputs $BAM \
	--refFile _PATH.FASTA_ \
	--regions $REGION \
	--output _WRKDIR_/$GRPID/software.$SPLIT.vcf \
	--outputRefCalls 0 \
	--logFilename _WRKDIR_/$GRPID/software.$SPLIT.log \
	--workDir _TMPDIR_/$GRPID.$SPLIT  \
	_SOFTWARE.OPTION[ \
	]_

bgzip _WRKDIR_/$GRPID/software.$SPLIT.vcf
tabix -f -p _WRKDIR_/$GRPID/software.$SPLIT.vcf

EOF

	}
	elsif ($caller =~ /^DeepVar/i) {
		if (exists $conf{PATH}{STAGEIN}) {
			warn "Stage-in is not supported for WGS calling by DeepVariant!";
		}
		$script .= << 'EOF';

read GRPID BAM < _PARDIR_/IID2BAM._INDEX_

IID=$GRPID

REGIONS=`cat _PARDIR_/SEQGRP.$SPLIT | sed -e 's/,/ /g'`

EOF
		my $CALLER = uc($caller);
		if (!defined $conf{$CALLER}{Docker} || $conf{$CALLER}{Docker} =~ /^Y|T/) {
		#if ($conf{$caller}{Docker} =~ /^T|Y/) {
			$script .= <<'EOF'
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
	--reads $BAM \
	--regions "$REGIONS" \
	--examples _WRKDIR_/$IID/software.$SPLIT.examples.tfrecords.gz 

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
	--outfile _WRKDIR_/$IID/software.$SPLIT.vcf 


bgzip _WRKDIR_/$IID/software.$SPLIT.vcf
tabix -f -p vcf _WRKDIR_/$IID/software.$SPLIT.vcf.gz

EOF
		}
		else {
			if (exists $conf{$CALLER}{CONDAENV}) {
			$script .= <<'EOF';
source _SOFTWARE.CONDARC_
conda activate _SOFTWARE.CONDAENV_

EOF
			}
			$script .= <<'EOF';

python _SOFTWARE.BINPATH_/make_examples.zip \
	--mode calling \
	--ref _PATH.FASTA_ \
	--reads $BAM \
	--regions "$REGIONS" \
	--examples _WRKDIR_/$IID/software.$SPLIT.examples.tfrecords.gz 

python _SOFTWARE.BINPATH_/call_variants.zip \
	--outfile _WRKDIR_/$IID/software.$SPLIT.variants.tfrecords.gz \
	--examples _WRKDIR_/$IID/software.$SPLIT.examples.tfrecords.gz \
	--checkpoint _SOFTWARE.MODEL_

python _SOFTWARE.BINPATH_/postprocess_variants.zip \
	--ref _PATH.FASTA_ \
	--infile _WRKDIR_/$IID/software.$SPLIT.variants.tfrecords.gz \
	--outfile _WRKDIR_/$IID/software.$SPLIT.vcf 

bgzip _WRKDIR_/$IID/software.$SPLIT.vcf
tabix -f -p vcf _WRKDIR_/$IID/software.$SPLIT.vcf.gz

EOF
		}
	}
	else {
		croak "Caller $caller is not supported";
	}
	if ($hapflag) {
		if ($conf{WHATSHAP}{CONDAENV}) {
			$script .= <<'EOF';
conda activate _WHATSHAP.CONDAENV_

EOF
		}
		if ($caller =~ /^DeepVar/i) {
			$script .= <<'EOF';

mv _WRKDIR_/$IID/software.$SPLIT.vcf.gz _WRKDIR_/$IID/software.$SPLIT.unphased.vcf.gz 
mv _WRKDIR_/$IID/software.$SPLIT.vcf.gz.tbi _WRKDIR_/$IID/software.$SPLIT.unphased.vcf.gz.tbi 

whatshap phase -r _PATH.FASTA_ _WHATSHAP.OPTION[ ]_ -o _WRKDIR_/$IID/software.$SPLIT.vcf  \
	_WRKDIR_/$IID/software.$SPLIT.unphased.vcf.gz $BAM

bgzip _WRKDIR_/$IID/software.$SPLIT.vcf
tabix -f -p vcf _WRKDIR_/$IID/software.$SPLIT.vcf.gz 

EOF
		}
		else {
			$script .= <<'EOF';

INPUTS2=$(echo $BAM | sed -e 's/,/ /g')

mv _WRKDIR_/$GRPID/software.$SPLIT.vcf.gz _WRKDIR_/$GRPID/software.$SPLIT.unphased.vcf.gz
mv _WRKDIR_/$GRPID/software.$SPLIT.vcf.gz.tbi _WRKDIR_/$GRPID/software.$SPLIT.unphased.vcf.gz.tbi

whatshap phase -r _PATH.FASTA_ _WHATSHAP.OPTION[ ]_ -o _WRKDIR_/$GRPID/software.$SPLIT.vcf \
	_WRKDIR_/$GRPID/software.$SPLIT.unphased.vcf.gz $INPUTS2

bgzip _WRKDIR_/$GRPID/software.$SPLIT.vcf
tabix -f -p vcf _WRKDIR_/$GRPID/software.$SPLIT.vcf.gz

EOF
		}
	}
	$script =~ s/software/${(\ lc($caller) )}/g;
	$script =~ s/SOFTWARE/${(\ uc($caller) )}/g;
	return $script;
}

sub get_exp_split {
	my $caller = shift @_;
	my $step = shift @_;
	my @allexps;
	foreach my $iid (@grps) {
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

read GRPID BAM < _PARDIR_/GRPID2BAM.$II

INPUT=`perl -e 'print join(" ", map {"INPUT=_WRKDIR_/$ARGV[1]/software.$_.vcf.gz" } 1..$ARGV[0])' _PARAM.NSPLIT_ $GRPID`

picard _PICARD.JAVAOPT_ MergeVcfs \
	$INPUT \
	OUTPUT=_SOFTWAREDIR_/$GRPID.vcf.gz

EOF
	# Mutect also need to merge stats
	if ($caller =~ /^mutect/i) {
		$script .= <<'EOF';

STATS=`perl -e 'print join(" -stats ", map {"INPUT=_WRKDIR_/$ARGV[1]/software.$_.vcf.gz.stats" } 1..$ARGV[0])' _PARAM.NSPLIT_ $GRPID`

gatk --java-options "_SOFTWARE.JAVAOPT_" -stats $STATS -O _SOFTWAREDIR_/$GRPID.vcf.gz.stats

EOF
	}
	$script =~ s/software/${(\ lc($caller) )}/g;
	$script =~ s/SOFTWARE/${(\ uc($caller) )}/g;
	return $script;
}


