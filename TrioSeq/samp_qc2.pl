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
use String::ShellQuote;
use List::MoreUtils qw|all any uniq|;
use Config::Std;
use Hash::Util qw|lock_hash_recurse|;
use Utils::Hash qw|merge_conf|;
use Utils::File qw|open_file|;
use Utils::Workflow;

use lib "$Bin/../lib";
use Shared qw|check_conf|;
use Wrapper qw|vcf_sampids|;

############################
## Command line interface ##
############################

my @spec =  (
	Param("conf|c")->valid(sub { -r }),
	Param("vcf|v")->valid(sub { -r }),
	Param("ped|p")->valid(sub { -r }),
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
	This is the non-interactive of sample level QC for family-based sequencing study. 
	It is suitable for large scale study.

Usage:
	samp_qc2.pl --conf Config --vcf VCF [--ped PED] --outdir OutDir

Options:
	--vcf: Cohort level joint genotyping VCF file. It can also be a list of VCF files.
	--ped: Plink style pedigree file. It overrides the PED file defined in the config.
		   Only samples exists in both VCF and PED file will be used in QC.

Notes:
	Sequencing depths encoded in genotype fields from VCF file are extracted for each individual 
	and will be used to verify sample sex (based on normalized depth chrX and chrY) and to identify
	chr21 trisomy.

	SNP genotypes will be converted to plink binary genotype format and used calculate pairwise
	kinship coefficient and validate pedigree relationships. 

Input/output:
	SNPGeno.{bim,bed,fam} - SNP genotypes in binary PED format, used in kinship estimates. 
	King.kin.txt - Annotated within-family relatedness estimates and known pedigree relationships.
	King.kin0.txt - Annotated between-family relatedness estimates.
	King.png - Plot of pairwise kinship estimates vs known pedigree relationships.
	SiteDP_merged.txt - Normalized depth of coverage of chr21, chrX and chrY and all autosomes.
	SiteDP_chrXY.png - Sex check using normalized chrX/chrY depth.
	SiteDP_chr21.png - Cumulative distribution plot of chr21 normalized depth.

EOF
	exit 1;
}


$opt->validate({ requires => [qw|conf vcf outdir|] });

my $rootdir = $opt->get_outdir;

my %conf = merge_conf($opt->get_conf, $opt->get_param); 
#read_config $opt->get_conf => my %conf;

check_conf(\%conf);

# If PED file is provided, we will over-ride the PED file in config
if ($opt->get_ped) {
	$conf{PEDIGREE}{FILE} = $opt->get_ped;
	unless(-f $conf{PEDIGREE}{FILE}) {
		die "Cannot find variant table: $conf{PEDIGREE}{FILE}";
	}
}
#my %known = map { (split)[1] => 1 } $conf{PEDIGREE}{FILE};

if (defined $conf{PEDIGREE}{TWINS}) {
	unless(-f $conf{PEDIGREE}{TWINS}) {
		die "Cannot find twin pairs list: $conf{PEDIGREE}{TWINS}";
	}
	$conf{PEDIGREE}{TWINS} = "--twins $conf{PEDIGREE}{TWINS}";
}
else {
	$conf{PEDIGREE}{TWINS} = "";
}

my $vcflist = $opt->get_vcf;
unless(-f $vcflist) {
	die "Cannot find VCF file or list: $vcflist";
}

$conf{VCFDEPTH}{AUTOSOME} = [map { ($_, "chr$_") } 1..22];
$conf{PATH} = { VCF => $vcflist, MODULE => shell_quote("$Bin/module"), UTIL => shell_quote("$Bin/utils") };



#############################
## Input files validation  ##
#############################

my ($vcfsplit, @allvcfs);
unless ($vcflist =~ /\.vcf\.gz$/ || $vcflist =~ /\.vcf$/) {
	$vcfsplit = 1;
	open my $fin, $vcflist or die "Cannot open VCF list for reading";
	while(<$fin>) {
		my $vcffile = (split)[0];
		unless($vcffile =~ /\.vcf\.gz$/ || $vcffile =~ /\.vcf$/) {
			die "Incorrect VCF file name: $vcffile";
		}
		else {
			unless(-f $vcffile) {
				die "Cannot find VCF file: $vcffile";
			}
		}
		push @allvcfs, $vcffile;
	}
}


# Read VCF or VCF list and determine the splits for calculating depth
# For VCF split, the first chromosome seen will be taken as the chromosome of that VCF.
# Only canonical chromosomes will be kept in depth calculation
my (%chrom2vcf, @filtvcfs);
unless ($vcfsplit) {
	foreach my $chrom (qw|21 X Y auto|) {
		push @{$chrom2vcf{$chrom}} => $vcflist;
	}
}
else {
	foreach my $vcffile (@allvcfs) {
		my $vcfchr = vcf_chrom($vcffile);
		if ($vcfchr ne 'Y' && $vcfchr ne 'chrY') {
			push @filtvcfs => $vcfchr;
		}
		foreach my $chrom (qw|21 X Y|) {
			if ($vcfchr eq $chrom || $vcfchr eq "chr$chrom") {
				push @{$chrom2vcf{$chrom}} => $vcffile;
			}
		}
		if ($vcfchr =~ /^(chr)?\d+$/) {
			push @{$chrom2vcf{auto}} => $vcffile;
		}
	}
	unless(@filtvcfs) {
		print STDERR "No VCF file left?"; exit 1;
	}
}


#################################
##  Workflow initialization    ##
##  Working directory setup    ##
#################################

my $wkf = Utils::Workflow->new($rootdir,
	{ engine => $opt->get_engine, force => $opt->get_force });

my @depthout;
my $depthct = 1;
if ($vcfsplit) {
	foreach my $chrom (qw|auto 21 X Y|) {
		for(my $ii = 0; $ii < @{$chrom2vcf{$chrom}}; $ii ++) {
			my $jj = $ii + 1;
			open my $fout, ">$rootdir/par/VCFDEPTH.$depthct" or die "Cannot write to VCFDEPTH.$depthct";
			print $fout $chrom2vcf{$chrom}[$ii], "\t", $chrom, "\t", $jj, "\n"; 
			push @depthout, $chrom eq 'auto' ? "SiteDP_$chrom.$jj" : "SiteDP_chr$chrom.$jj";
			$depthct ++;
		}
	}
	open my $flst, ">$rootdir/par/SNPGeno.list" or die "Cannot write to SNPGeno.list";
	for(my $ii = 0; $ii < @filtvcfs; $ii ++) {
		my $jj = $ii + 1;
		open my $fout, ">$rootdir/par/VCF.$jj" or die "Cannot write to VCF.$jj";
		print $fout $filtvcfs[$ii], "\n";
		print $flst "$rootdir/wrk/SNPGeno.$jj\n";
	}
}
else {
	foreach my $chrom (qw|auto 21 X Y|) {
		open my $fout, ">$rootdir/par/VCFDEPTH.$depthct" or die "Cannot write to VCFDEPTH.$depthct";		
		print $fout $chrom2vcf{$chrom}[0], "\t", $chrom, "\t1\n";
		push @depthout, $chrom eq 'auto' ? "SiteDP_$chrom.1" : "SiteDP_chr$chrom.1";
		$depthct ++;
	}
}

$wkf->add(vcf_depth(), { name => "VcfDepth", nslots => scalar(@depthout),
 						 expect => [ map { "wrk/$_.idepth" } @depthout ] })
	->add(plot_depth(), { name => "PlotDepth", depend => "VcfDepth",
						   expect => [qw|out/SiteDP_merged.txt out/SiteDP_chr21.png out/SiteDP_chrXY.png|] });
if ($vcfsplit) {
	$wkf->add(vcf2bped_split(), { name => "Vcf2Bped", nslots => scalar(@filtvcfs),
					expect => [ map { ["wrk/SNPGeno.$_.fam", "wrk/SNPGeno.$_.bed", "wrk/SNPGeno.$_.bim" ] } 
								1..scalar(@filtvcfs) ]  });
}
else {
	$wkf->add(vcf2bped(), { name => "Vcf2Bped",
						 	expect => ["wrk/SNPGeno.fam", "wrk/SNPGeno.bed", "wrk/SNPGeno.bim" ] });
}
$wkf->add(merge_geno(), { name => "MergeGeno", depend => "Vcf2Bped",
						  expect => ["out/SNPGeno.fam", "out/SNPGeno.bed", "out/SNPGeno.bim" ] })
	->add(king(), { name => "King", depend => "MergeGeno", 
					expect => ["out/King.kin", "out/King.kin0"] })
	->add(plot_king(), { name => "PlotKing", depend => "King",
						 expect => ["out/King.kin", "out/King.kin0", "out/King.png"] });


$wkf->inst(\%conf);
$wkf->run({ conf => $conf{$wkf->{engine}}, dryrun => $opt->get_dryrun });


############################
## Workflow components    ##
############################

# Convert VCF to plink format, and keep biallelic SNP genotypes
# Note that chrY may be filtered out if the same threshold is applied
# In case of multiple VCF inputs, only autosomes and chrX will be kept
sub vcf2bped {
	my $script = <<'EOF';

_SNPGENO.PLINK_ --vcf _PATH.VCF_ _SNPGENO.OPTION[ ]_  \
	--make-bed --out _WRKDIR_/SNPGeno 
 
EOF
}

sub vcf2bped_split {
	my $script = <<'EOF';

read VCF < _PARDIR_/VCF._INDEX_

_SNPGENO.PLINK_ --vcf $VCF _SNPGENO.OPTION[ ]_ \
	--make-bed --out _WRKDIR_/SNPGeno._INDEX_

EOF
}

# If case of multiple VCFs, need to concat multiple converted SNP geno files
sub merge_geno {
	my $script;
	if ($vcfsplit) {
		$script = <<'EOF';

plink --merge-list _PARDIR_/SNPGeno.list \
	--allow-extra-chr --make-bed --out _WRKDIR_/SNPGeno

EOF
	}
	# Test if we need to remove samples
	$script .= <<'EOF';

table_intersect.pl -a _WRKDIR_/SNPGeno.fam -a-select 1,2 -b _PEDIGREE.FILE_ -b-select 2 \
	--negative --out _WRKDIR_/BadSamps.txt

plink --bfile _WRKDIR_/SNPGeno --allow-extra-chr --remove _WRKDIR_/BadSamps.txt --make-bed --out _OUTDIR_/SNPGeno

table_intersect.pl -a _OUTDIR_/SNPGeno.fam -a-select 2 -b _PEDIGREE.FILE_ -b-noheader | \
	awk '{print $2, $1, $3, $4, $5, $6}' > _OUTDIR_/SNPGeno_Update.fam

NLINE=$(cat _OUTDIR_/SNPGeno.fam | wc -l)
NLINE2=$(cat _OUTDIR_/SNPGeno_Update.fam | wc -l)

if [[ $NLINE == $NLINE2 ]]; then
	mv _OUTDIR_/SNPGeno.fam _OUTDIR_/SNPGeno.fam.bak
	mv _OUTDIR_/SNPGeno_Update.fam _OUTDIR_/SNPGeno.fam
else
	echo "Incorrect number of samples in PED file after update!" 2>&1
fi

EOF
}

# Depth of coverage per sample from VCF DP fields
# In case of multiple VCFs, need to merge VCFs before running this script
# Need to calculate VCF depth separately on each VCFs, then take weighted sums.
sub vcf_depth {
	my $script = <<'EOF';	

read VCF CHROM II < _PARDIR_/VCFDEPTH._INDEX_

if [[ $CHROM == "auto" ]]; then
	vcftools --gzvcf $VCF --chr _VCFDEPTH.AUTOSOME[ --chr ]_ \
		--depth _VCFDEPTH.OPTION[ ]_ --out _WRKDIR_/SiteDP_$CHROM.$II
else
	vcftools --gzvcf $VCF --chr $CHROM --chr chr$CHROM --depth \
		_VCFDEPTH.OPTION[ ]_ --out _WRKDIR_/SiteDP_chr$CHROM.$II
fi

EOF
}

sub vcf_chrom {
	my ($vcf) = @_;
	my $fin = open_file($vcf);
	my $chr;
	while(<$fin>) {
		next if /^#/;
		$chr = (split)[0];
		last;
	}
	return $chr;
}

# Plot VCF depth for sex check and calling 21 trisomy 
sub plot_depth {
	my $script = << 'EOF';

perl _PATH.MODULE_/merge_idepth.pl _WRKDIR_/SiteDP _PEDIGREE.FILE_ _OUTDIR_/SiteDP_merged.txt

Rscript _PATH.MODULE_/plot_idepth.R _OUTDIR_/SiteDP

EOF
}

# Run king to estimate pairwise relatedness
sub king {
	my $script = <<'EOF';

king -b _OUTDIR_/SNPGeno.bed --related _KING.OPTION[ ]_ --prefix _OUTDIR_/King

EOF
}

# Reformat and plot king results
sub plot_king {
	my $script;
	if (defined $conf{SAMPINFO}) {
		$script = <<'EOF';

perl _PATH.UTIL_/anno_king_out.pl --input _OUTDIR_/King --ped _PEDIGREE.FILE_ \
	--ignore _PEDIGREE.IGNORE_ _PEDIGREE.TWINS_ --output _OUTDIR_/King \
	--sxref _SAMPINFO.TABLE[ ]_ \
	--sxref-fields _SAMPINFO.FIELDS[ ]_

EOF
	}
	else {
		$script = <<'EOF';

perl _PATH.UTIL_/anno_king_out.pl --input _OUTDIR_/King --ped _PEDIGREE.FILE_ \
	--ignore _PEDIGREE.IGNORE_ _PEDIGREE.TWINS_ --output _OUTDIR_/King

EOF
	}
	$script .= <<'EOF';

Rscript _PATH.MODULE_/plot_king.R _OUTDIR_/King

EOF
	return $script;
}
