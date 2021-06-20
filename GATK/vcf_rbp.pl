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
use File::Basename;
use File::Path qw|make_path|;
use Cwd qw|abs_path|;
use Getopt::Lucid qw|:all|;
use String::ShellQuote;
use List::MoreUtils qw|all any uniq|;
use Config::Std;
use Hash::Util qw|lock_hash_recurse|;
use Utils::Workflow;
use Utils::Hash qw|merge_conf|;
use Utils::List qw|all_combs|;

use lib "$Bin/../lib/";
use Shared qw|read_list check_conf|;
use Wrapper qw|vcf_sampids|;

############################
## Command line interface ##
############################

my @spec =  (
	Param("conf|c")->valid(sub { -r }),
	Param("vcf|v")->valid(sub { -r }),
	Param("list|l")->valid(sub { -r }),
	Param("group|g")->valid(sub { -r }),
	Param("rename")->valid(sub { -r }),
	Param("remove")->valid(sub { -r }),
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
	This is a pipeline script to do read-backed phasing using whatshap.

Usage:
	vcf_rbp.pl --conf Config --vcf VCFList [--list CRAMList] --outdir outdir

Options:
	--vcf: VCF file, file list or directory. See notes below on how they will be processed. 
	--group: Define individual group membership when a single large joint genotyping VCF is provided.
			It should have 2 columns per line, IID and GroupID. IID should be found in VCF file header.
	--list: BAM/CRAM file list or directory. It will override the cram list specified in the config.
	--rename: Sample rename list, overrides the sample rename list in config
	--remove: A list of bad samples to be removed, overrides the bad sample list in config.
	
Notes:
	If a single VCF file is provided, it is presumed to be from cohort wide joint genotyping.
	Each individual or group of individuals (if groups are defined) will be extracted and processed
	separately. Only individuals with BAM/CRAM files will be selected for processing.

	If a VCF file list is provided, there are two possibilities: 1. VCF files in the list can be 
	segmented large VCFs from cohort wide genotyping split in different regions. In this case
	different VCFs should contain the same set of individuals,  and is sorted by position.
	Each individual or group of individuals (if groups are defined) will be extracted and merged 
	across splits before processing. 2. VCF files in the list can also be one file per group (or sample), 
	an optional second column can be used to specify group ID. In this case, VCF list will be 
	processed in the same way as VCF directory below.

	If a directory of VCF files is provided, all samples in each VCF will be treated as a group and 
	processed together. If group membership is provided by group option, they will be ignored.
 	Only VCF files containing at least one sample in the BAM/CRAM list will be processed. 

 	To look for aligned reads for samples in VCF file, WhatsHap makes use of sample IDs from CRAM/BAM 
 	file header which may be different from VCFs. The pipeline use align_vcf_bam to align samples 
 	in VCF and BAMs accounting for the possible differences in IDs .

Dependencies:
	bcftools, whatshap 	

EOF
	exit 1;
}

$opt->validate({requires => [qw|conf vcf outdir|]});

my $rootdir = $opt->get_outdir;

my %conf = merge_conf($opt->get_conf, $opt->get_param); 
#read_config $opt->get_conf => my %conf;
check_conf(\%conf);

$conf{PATH}{UTIL} = shell_quote("$Bin/utils");


############################################
## Input files validation and collection  ##
############################################

# Read CRAM/BAM files
if ($opt->get_list) {
	$conf{PATH}{CRAMLIST} = $opt->get_list;
}
if ($opt->get_remove) {
	$conf{PATH}{REMOVE} = $opt->get_remove;
}
if ($opt->get_rename) {
	$conf{PATH}{RENAME} = $opt->get_rename;
}
foreach my $FILE (qw|CRAMLIST BAMLIST REMOVE RENAME|) {
	if (defined $conf{PATH}{$FILE}) {
		unless(-f $conf{PATH}{$FILE} || $FILE =~ /LIST$/ && -d $conf{PATH}{$FILE}) {
			die "Cannot find $FILE file: $conf{PATH}{$FILE}";
		}
	}
}
unless(defined $conf{PATH}{CRAMLIST} || defined $conf{PATH}{BAMLIST}) {
	die "Must provide CRAM or BAM file list!";
}

my ($bams, $filetype) = read_list(defined $conf{PATH}{CRAMLIST} ? $conf{PATH}{CRAMLIST} : $conf{PATH}{BAMLIST}, 
								{ suffix => { bam => qr/\.bam$/, cram => qr/\.cram$/ }, 
								  remove => $conf{PATH}{REMOVE}, rename => $conf{PATH}{RENAME} });

# Group membership file
my %group;
if ($opt->get_group) {
	%group = map { (split)[0,1] } slurp $opt->get_group;
}

# VCF files
# @allvcfs: When VCF list of large cohort-wide genotyping is provided: All VCF file paths 
# @allsamps: ..., All sample IDs for cohort-wide genotyping VCFs
# %grp2bam: Map group ID to BAM file path
# %grp2vcf: When VCF directory or a list of family or indiv level VCF files is provided: Map group ID to VCF path
# %sampbygrp: When group file is provided, map sample ID to group ID

my (@allvcfs, @allsamps, %grp2bam, %grp2vcf, %sampbygrp);
my $vcflist = $opt->get_vcf;

# Model: Bulk, BulkSplit, SingleList or SingleDir
my $vcfmode;

if ($vcflist =~ /\.vcf\.gz$/ || $vcflist =~ /\.vcf$/) {
	print STDERR "VCF file: $vcflist\n";
	unless (-f $vcflist) {
		die "Cannot find input VCF file: $vcflist";
	}
	$conf{PATH}{VCF} = $vcflist; 
	@allsamps = grep { defined $bams->{$_} } vcf_sampids($vcflist);
	unless (@allsamps > 0) {
		print STDERR "No sample in VCF is associated with CRAM/BAM files\n";
		exit 1;
	}
	foreach my $iid (@allsamps) {
		# Fall back to IID if GID cannot be found
		my $gid = $group{$iid} // $iid;
		push @{$sampbygrp{$gid}} => $iid;
		push @{$grp2bam{$gid}} => $bams->{$iid};
	}
	$vcfmode = "Bulk";
}
elsif (-f $vcflist) {
	# Do a quick test for first line or first two VCFs
	my $mode;
	{
		open my $fin, $vcflist or die "Cannot open VCF list for reading";
		my $line1= <$fin>;
		my @dat = split(/\s+/, $line1);
		if (@dat >= 2) {
			$mode = "List";
		} 
		else {
			if (@dat == 0) {
				die "No data in $vcflist";
			}
			unless(-f $dat[0]) {
				die "Cannot open VCF file $dat[0]!";
			}
		}
		unless(defined $mode) {
			my $samplist1 = join(",", sort vcf_sampids($dat[0]));
			my $line2 = <$fin>;
			@dat = split(/\s+/, $line2);
			if (@dat >= 1) {
				die "Cannot open VCF file $dat[0]"; 
			}
			my $samplist2 = join(",", sort vcf_sampids($dat[0]));
			if ($samplist1 eq $samplist2) {
				$mode = "Split";
			}
			else {
				$mode = "List";
			}
		}
	}

	if ($mode eq 'Split') {
		print STDERR "Reading VCF splits from list $vcflist\n";
		my @vcfsamps;
		open my $fin, $vcflist or die "Cannot open VCF list for reading";
		while(<$fin>) {
			my $vcffile = (split)[0];
			unless ($vcffile =~ /\.vcf\.gz$/ || $vcffile =~ /\.vcf$/) {
				die "Incorrect VCF file name: $vcffile";
			}
			else {
				unless(-f $vcffile) {
					die "Cannot find VCF file: $vcffile";
				}
			}
			unless (@vcfsamps) {
				@vcfsamps = vcf_sampids($vcffile);
			}
			else {
				#my @vcfsamps2 = vcf_sampids($vcffile);
				my $samplist1 = join(",", sort @vcfsamps);
				my $samplist2 = join(",", sort(vcf_sampids($vcffile)));
				unless($samplist1 eq $samplist2) {
					print STDERR $samplist1, "\n", $samplist2, "\n";
					die "VCF files in list $vcflist do not have identical set of samples!";
				}
			}
			unless (@allsamps) {
				@allsamps = grep { defined $bams->{$_} } @vcfsamps;
				unless (@allsamps > 0) {
					print STDERR "No sample in VCF is associated with CRAM/BAM files\n";
					exit 1;
				}
				foreach my $iid (@allsamps) {
					my $gid = $group{$iid} // $iid;
					push @{$sampbygrp{$gid}} => $iid;
					push @{$grp2bam{$gid}} => $bams->{$iid};
				}
			}
			push @allvcfs, $vcffile;
		}
		$vcfmode = "BulkSplit";
	}
	elsif ($mode eq 'List') {
		print STDERR "Reading VCF files from list $vcflist\n";
		my $count = 0;
		open my $fin, $vcflist or die "Cannot open VCF list for reading";
		while(<$fin>) {
			my ($vcffile, $grpid) = (split)[0,1];
			unless($vcffile =~ /\.vcf\.gz$/ || $vcffile =~ /\.vcf$/) {
				die "Incorrect VCF file name: $vcffile";
			}
			else {
				unless(-f $vcffile) {
					die "Cannot find VCF file: $vcffile";
				}
			}
			my @subsamps = grep { defined $bams->{$_} } vcf_sampids($vcffile);
			next unless @subsamps > 0;
			unless (defined $grpid) {
				$grpid = basename($vcffile); 
				$grpid =~ s/\.vcf(\.gz)?$//;
			}
			if (defined $grp2vcf{$grpid}) {
				die "Sample group $grpid is alreay defined!";
			}
			$grp2vcf{$grpid} = $vcffile;
			foreach my $iid (@subsamps) {
				push @{$grp2bam{$grpid}} => $bams->{$iid};
				push @{$sampbygrp{$grpid}} => $iid;
			}
			$count ++;
			if ($count % 500 == 0) {
				print STDERR "$count VCFs read\n";
			}
		}
		if ($count == 0) {
			print STDERR "No VCF file have overlapping sample with BAM/CRAM files!";
			exit 1;
		}
		$vcfmode = "SingleList";
	}
}
elsif (-d $vcflist) {
	print "Reading VCFs from directory $vcflist\n";
	my @vcffiles = grep { /\.vcf\.gz$/ || /\.vcf$/ } IO::Dir->new($vcflist)->read();
	unless(@vcffiles > 0) {
		print STDERR "No VCF file can be found in $vcflist\n";
		exit 1;
	}
	my $count = 0;
	foreach my $vcffile (@vcffiles) {
		my $gid = basename($vcffile); 
		$gid =~ s/\.vcf(\.gz)?$//;
		my @subsamps = grep { defined $bams->{$_} } vcf_sampids("$vcflist/$vcffile");
		next unless @subsamps > 0;
		if (defined $grp2vcf{$gid}) {
			die "Sample group $gid is alreay defined: $vcffile!";
		}
		$grp2vcf{$gid} = "$vcflist/$vcffile";
		foreach my $iid (@subsamps) {
			push @{$grp2bam{$gid}} => $bams->{$iid};
			push @{$sampbygrp{$gid}} => $iid;
		}
		$count ++;
		if ($count % 500 == 0) {
			print STDERR "$count VCFs read\n";
		}
	}
	if ($count == 0) {
		print STDERR "No VCF file have overlapping sample with BAM/CRAM files!";
		exit 1;
	}
	$vcfmode = "SingleDir";
}
else {
	die "Cannot read from VCF list: $vcflist";
}

#################################
##  Workflow initialization    ##
##  Working directory setup    ##
#################################

my $wkf = Utils::Workflow->new($rootdir, { engine => $opt->get_engine, force => $opt->get_force });

# GRPID2BAM.JJ, GRPID2VCF.JJ, VCFSPLIT.JJ
my @allgrps = sort keys %grp2bam;
for(my $ii = 0; $ii < @allgrps; $ii ++) {
	my $jj = $ii + 1;
	open my $fout, ">$rootdir/par/GRPID2BAM.$jj" or die "Cannot write to GRPID2BAM.$jj";
	print $fout $allgrps[$ii], "\t", join(',', uniq sort @{$grp2bam{$allgrps[$ii]}}), "\n";
}

#if (-f $vcflist && $vcflist !~ /\.vcf(\.gz)?$/) {
if ($vcfmode eq 'BulkSplit') {
	for(my $ii = 0; $ii < @allvcfs; $ii ++) {
		my $jj = $ii + 1;
		open my $fout, ">$rootdir/par/VCFSPLIT.$jj" or die "Cannot write to VCFSPLIT.$jj";
		print $fout $allvcfs[$ii], "\n";
	}
	$conf{PARAM}{NSPLIT} = scalar(@allvcfs);
}
else {
	$conf{PARAM}{NSPLIT} = 1;
}

for(my $ii = 0; $ii < @allgrps; $ii ++) {
	my $jj = $ii + 1;
	open my $fout, ">$rootdir/par/GRPID2VCF.$jj" or die "Cannot write to GRPID2VCF.$jj";
	if ($vcfmode =~ /^Single/) {
		print $fout $allgrps[$ii], "\t", $grp2vcf{$allgrps[$ii]}, "\n";
	}
	else {
		print $fout $allgrps[$ii], "\t", "$rootdir/wrk/$allgrps[$ii].vcf.gz", "\n";
		open my $flst, ">$rootdir/wrk/$allgrps[$ii].keep.txt" or die "Cannot write to $allgrps[$ii].keep.txt";
		print $flst join("\n", @{$sampbygrp{$allgrps[$ii]}}), "\n";
	}
	open my $flst, ">$rootdir/wrk/$allgrps[$ii].bam_list.txt" or die "Cannot write to $allgrps[$ii].bam_list.txt";
	#print $flst join("\n", @{$sampbygrp{$allgrps[$ii]}}), "\n";
	foreach my $iid (@{$sampbygrp{$allgrps[$ii]}}) {
		print $flst $bams->{$iid}, "\t", $iid, "\n";
	}
}

if ($vcfmode =~ /^Single/) {
	for(my $ii = 0; $ii < @allgrps; $ii ++) {
		my $jj = $ii + 1;
		open my $fout, ">$rootdir/par/GRPID2VCF.$jj" or die "Cannot write to GRPID2VCF.$jj";
		print $fout $allgrps[$ii], "\t", $grp2vcf{$allgrps[$ii]}, "\n";
	}
}
else {
	for(my $ii = 0; $ii < @allgrps; $ii ++) {
		my $jj = $ii + 1;
		open my $fout, ">$rootdir/par/GRPID2VCF.$jj" or die "Cannot write to GRPID2VCF.$jj";
		print $fout $allgrps[$ii], "\t", "$rootdir/wrk/$allgrps[$ii].vcf.gz", "\n";
		open my $flst, ">$rootdir/wrk/$allgrps[$ii].keep.txt" or die "Cannot write to $allgrps[$ii].keep.txt";
		print $flst join("\n", @{$sampbygrp{$allgrps[$ii]}}), "\n";
	}
}

##############################
## Workflow initialization  ##
##############################

my %opt;
if ($vcfmode =~ /^Bulk/) {
	if ($vcfmode eq 'Bulk') {
		$wkf->add(vcfsub(), { name => "VCFSub", nslots => scalar(@allgrps), 
							  expect => [map { ["wrk/$_.vcf.gz", "wrk/$_.vcf.gz.tbi"] } @allgrps] });
		%opt = (depend => "VCFSub", deparray => 1);
	}
	else {
		$wkf->add(vcfsub_split(), { name => "VCFSubSplit", 
									nslots => scalar(@allgrps) * $conf{PARAM}{NSPLIT}, 
									expect => [ map { ["wrk/$_->[1].$_->[0].vcf.gz", "wrk/$_->[1].$_->[0].vcf.gz.tbi"] } 
												all_combs([1..$conf{PARAM}{NSPLIT}], \@allgrps) ] })
			->add(vcfmerge(), { name => "VCFMerge", nslots => scalar(@allgrps) * scalar(@allvcfs),
								depend => "VCFSubSplit", deparray => 1, step => $conf{PARAM}{NSPLIT},
								expect => [map { ["wrk/$_.vcf.gz", "wrk/$_.vcf.gz.tbi"] } @allgrps] });
		%opt = (depend => "VCFMerge", deparray => 1, step => $conf{PARAM}{NSPLIT});
	}
}
$wkf->add(whatshap(), { name => "WhatsHap", nslots => scalar(@allgrps) * $conf{PARAM}{NSPLIT}, %opt,
						expect => [map { ["out/$_.vcf.gz", "out/$_.vcf.gz.tbi"] } @allgrps] });


$wkf->inst(\%conf);
$wkf->run({ conf => $conf{$wkf->{engine}}, dryrun => $opt->get_dryrun  });


############################
## Workflow components    ##
############################

# Extract selected samples from large pVCF
sub vcfsub {
	my $script = <<'EOF';

read GRPID BAMLIST < _PARDIR_/GRPID2BAM._INDEX_

bcftools view -S _WRKDIR_/$GRPID.keep.txt -a -c 1 -O z -o _WRKDIR_/$GRPID.vcf.gz _PATH.VCF_

tabix -p vcf _WRKDIR_/$GRPID.vcf.gz 

EOF
}

sub vcfsub_split {
	my $script = <<'EOF';

II=$(( (_INDEX_-1)/_PARAM.NSPLIT_+1 ))
SPLIT=$(( (_INDEX_-1)%_PARAM.NSPLIT_+1 ))

read VCFSPLIT < _PARDIR_/VCFSPLIT.$SPLIT
read GRPID BAMLIST < _PARDIR_/GRPID2BAM.$II

bcftools view -S _WRKDIR_/$GRPID.keep.txt -a -c 1 -O z -o _WRKDIR_/$GRPID.$SPLIT.vcf.gz $VCFSPLIT

tabix -p vcf _WRKDIR_/$GRPID.$SPLIT.vcf.gz

EOF
}

# Merge-sort VCFs from splits
sub vcfmerge {
	my $script = <<'EOF';

II=$(( (_INDEX_-1)/_PARAM.NSPLIT_+1 ))

read GRPID BAMLIST < _PARDIR_/GRPID2BAM.$II

for JJ in $(seq 1 _PARAM.NSPLIT_); do 
	echo _WRKDIR_/$GRPID.$JJ.vcf.gz 
done > _WRKDIR_/$GRPID.vcf_list.txt

bcftools concat -f _WRKDIR_/$GRPID.vcf_list.txt -a -d none -O z -o _WRKDIR_/$GRPID.vcf.gz

tabix -p vcf _WRKDIR_/$GRPID.vcf.gz

EOF
}

# Run whatshap for read backed phasing
# To account for the differences in sample IDs encoded in BAM/CRAM header and VCF header
# we use align_vcf_bam to modify VCF header before and after whatshap run
sub whatshap {
	my $script;
	if (defined $conf{WHATSHAP}{CONDAENV}) {
		$script .= <<'EOF'
source _PATH.CONDARC_
conda activate _WHATSHAP.CONDAENV_

EOF
	}
	$script .= <<'EOF';

II=$(( (_INDEX_-1)/_PARAM.NSPLIT_+1 ))

read GRPID BAMLIST < _PARDIR_/GRPID2BAM.$II
read GRPID VCF < _PARDIR_/GRPID2VCF.$II

INPUTS=$(echo $BAMLIST | sed -e 's/,/ /g')

perl _PATH.UTIL_/align_vcf_bam.pl --vcf $VCF --out _WRKDIR_/$GRPID.vcf.gz --bam _WRKDIR_/$GRPID.bam_list.txt

whatshap phase -r _PATH.FASTA_ _WHATSHAP.OPTION[ ]_ -o _OUTDIR_/$GRPID.vcf \
	_WRKDIR_/$GRPID.vcf.gz $INPUTS

bgzip _OUTDIR_/$GRPID.vcf

tabix -p vcf _OUTDIR_/$GRPID.vcf.gz

perl _PATH.UTIL_/align_vcf_bam.pl --vcf _OUTDIR_/$GRPID.vcf.gz --bam _WRKDIR_/$GRPID.bam_list.txt --reverse

EOF
}

