#/usr/bin/env perl

use strict;
use warnings;
use Carp;
use Scalar::Util qw(looks_like_number);
use Graph::Undirected;
use FindBin qw|$Bin|;
use Getopt::Lucid qw|:all|;
use IO::Prompt;
use String::ShellQuote;
use List::MoreUtils qw|all|;
use Perl6::Slurp;
use Data::Dumper;
#use Utils::File::Iter qw|iter_file|;
use Utils::Workflow;
use Utils::Hash qw|prompt_params|;
use Genome::UCSC qw|%PAR|;
use Genet::Ped;

# The default input is vcf, ped, and output directory
my @spec = (
	Param("vcf|v")->valid(sub { -r $_ && -r "$_.tbi" }),
	Param("ped|p")->valid(sub { -r }),
	Param("ignore")->default("_Re"),
	Param("reuse")->valid(sub { -d }),
	Param("outdir|out"),
	Switch("help|h"),
	Switch("dryrun"),
	Switch("force")
	);

my $opt = Getopt::Lucid->getopt( \@spec );

if ($opt->get_help) {
	#print STDERR "Usage: samp_qc.pl --vcf VCF --ped PED --out OUTDIR\n";
	#exit 1;
	print STDERR <<EOF;
Purpose:
	This is an interactive script for sample level QC of joint genotyping VCF from family-based sequencing study.

Usage:
	samp_qc.pl --ped PED --vcf VCF --outdir OutDir 

Options:
	--force: Over-write existing output directory.
	--reuse: Provide previous directory to this option, then re-use the site-depth results from previous run. 	

Notes:
	The aim of sample level QC is to verify that sample sex and familial relationships specified in the PED file is
	consistent with genotypes and sequencing depths in VCF. We may also identify chromosomal abnormalities from 
	QC process. 

	This script is suitable for analyzing small cohort (<5000 samples). There are two major time-consuming steps: 
	calculating average depth of coverage per sample, and calculating/plotting pairwise relatedness across all samples. 
	For large cohort, it is recommended to first convert VCF to BPED by plink2 then	use samp_snp_qc.pl for interactive QC.
	Sequencing depth information can also be collected from the output of bam_qc.pl.

	The script runs following major steps in order. Each step will display an interactive menu of adjustable parameters. 
	After each run, users can examine if desired output files are generated and optionally re-run the step after further
	adjusting parameters. 
		
		1. Convert VCF to plink binary format (BPED).
			It runs vcf2bped.pl utility for conversion. The parameters are generally good for GATK joint VCF from 
			high coverage exome sequencing. Only autosomes and chrX are variants are included in the output.
			The convereted bped file should be checked if it contain enough variants especially sex chromosomes. 
			Thresholds of some parameters can be relaxed to increase the number of variants passing QC.
		
		2. Calculate depth of coverage per sample for autosomes and chrX/Y from VCF DP field.
			This is one of the time-consuming steps, a thinning parameter can be specified to speed up DP calculation.
			If --reuse option is enabled, this step will be skipped.
		
		3. Run Peddy for sex and relatedness check and create interactive plots.
			Sex and relatedness check are also implemented in Peddy. However, Peddy sometimes produced inaccurate results
			because it only uses chrX SNP genotype to infer sex and it only selected common coding SNPs for calculating kinship 
			coefficient. But we still use Peddy as part of the pipeline because it can produce interactive plots useful
			for diagnosis of sample errors.
		
		4. Generating sample QC statistics from SNP genotypes
			The following sample level stats are calculated: Genotype call rate, inbreeding coefficients, heterozygosity rate, 
			Mendelian err, and runs of homozygosity. Not all of them are used in plotting later. But they can be useful for 
			diagnosis of outlier or err samples. 

		5. Run king and plink to calculate pairwise relationship statistics.
			King's kinship coefficient estimation is robust to population strcuture. Peddy implements the king's robust estimates
			but is applied only on known common SNPs which can create biased estimate due to ascertainment of markers. Plink estimates
			the Prob of IBD=0,1,2 assuming that samples are from a homogenous population, which can also be biased in samples 
			have diverse ancestries.
		
		6. Making various QC plots.
			The following QC plots are made:
			sample_DPvsGMissHQ.png, sample_HetDpvsGMissHQ.png
				-- Scatterplot of average depth at all or heterozygotes autosome SNPs vs genotype mmissingness.  

			sample_chrXFHet_byPEDSex.png, sample_chrXFHet_byPredSex.png 
				- Scatterplot of chrX inbreeding coefficient and het/hom ratio colored by pedigree or predicted sex,
				  highlighting samples with inconsistent PED vs predicted sex.
 			sample_chrXYNormDp_byPEDSex.png, sample_chrXYNormDp_byPredSex.png
				- Scatterplot of chrX and chrY normalized depth colored by pedigree or predicted sex,
				  highlighting samples with inconsistent PED vs predicted sex.
				  Prediction of sex in these plots are made based on user defined cutoff for normalized chrY depth (need to exam 
				  chrX/Y depth to determine the cutoff).  If chrY data are not available, we take Peddy's sex prediction which is
				  based on its default cutoff of chrX het/hom ratio.
			sample_PC1-4.png, sample_PC1_PC2.png, sample_PC1_PC3.png, sample_PC1_PC4.png, sample_PC2_PC3.png
				- Projection to first 4 PCs defined by K1G populations colored by predicted ancestry.
				  Results are taken from the Peddy's output, see also pop_pca.pl. 
			relpairs_KinshipIBS0_king.png, relpairs_KinshipIBS0_peddy.png
				- Scatterplot of kinship coefficient and IBS0 estimated by king or peddy for all pairs of samples in the cohort
				   colored by different pairs of known relationships.
			relpairs_IBD1_IBD2_plink.png, relpairs_IBD_Ternary_plink.png
				- Scatterplot and ternary plot of IBD estimates by plink for all pairwise of samples in the cohort colored by 
				  different pairs of known relationships.
				  In the above figures, known pairwise relationships are inferred from PED file. For samples in incomplete families,
				  the pairwise relationships may not be inferred correctly. PED file can also not distinguish sibs from MZ twins. 
				  Sample pairs whose pedigree relationships are inconsistent with inferred relationships are highlighted.
				  Kinship coefficients for all pairs of samples in the cohort are plotted. It is useful to identify cryptic related
				  but become very time consuming when sample size is over 1000.  
				 
		7. Create files for likely sample errors.
			If sample errors are found, they are listed in the following "patch" files s based on best guess.
			The real situation of sample error are complex, and should always be evaluated case-by-case.
			  patch.badsamps.txt - Samples with extremely low depth 
			  patch.sex_update.txt - Samples of incorrect sex
			  patch.nonpar.txt -  parent-offspring pairs that are not related by genotypes
			  patch.swapped.txt - Likely swapped pairs
			  patch.dups.txt  - Duplicated samples or twins
			  patch.rels.txt  - Cryptically related samples 
			The format above patch files are described in "utils/fix_ped_vcf.pl".
			  
Dependencies:
	plink, king, peddy, vcftools, csvkit, R (ggplot2, ggrepel, ggtern, GGally, grid)

EOF
	exit 1;
}

$opt->validate({ requires => [qw|vcf ped outdir|] });
my $rootdir = $opt->get_outdir;

# Create workflow
my $wkf = Utils::Workflow->new($rootdir, {engine => "BASH", force => $opt->get_force});

# Re-use part of previous results: site depth, assuming same sample names
if (my $refdir = $opt->get_reuse) {
	print STDERR "Re-using previous results on depth\n";
	if(all { -f "$refdir/wrk/site.auto.idepth" } qw|auto chrX chrY|) {
		my %samps = map { (split)[1] => 1} slurp $opt->get_ped;
		foreach my $mid (qw|auto chrX chrY|) {
			open my $fin, "$refdir/wrk/site.$mid.idepth" or die "Cannot open ref $mid.idepth";
			open my $fout, ">$rootdir/wrk/site.$mid.idepth" or die "Cannot open $mid.idepth";
			print $fout join("\t", qw|INDV N_SITES MEAN_DEPTH|), "\n";
			while(<$fin>) {
				my @dat = split;
				next unless defined $samps{$dat[0]};
				print $fout join("\t", @dat), "\n";
			}
		}
	}
	else {
		die "Cannot find site depth files from $refdir";
	}
}
else {
	$wkf->add(depth_script(), { name => "SiteDepth", expect => [map { "wrk/site.$_.idepth" } qw|auto chrX chrY|] });
}

# Adding components to workflow
$wkf->add(convert_script(), { name => "VCF2BPed" ,	expect => [map { "wrk/geno.$_" } qw|bim bed fam|],
							  callback => \&check_geno_files })
	->add(peddy_script(), 	{ name => "Peddy", 		expect => [ peddy_files() ]} )
	->add(genoqc_script(), 	{ name => "GenoQC", 	expect => [ genoqc_files() ] })
	->add(kinship_script(), { name => "Kinship", 	expect => [ kinship_files() ] })
	->add(plot_script(),  	{ name => "Plot", 		expect => [ plot_files() ] })
	->add(patch_script(),   { name => "Patch",   	expect => [] });

# Prompt user to confirm parameters
my %conf = (default => { VCF => $opt->get_vcf, PED => $opt->get_ped, MODULE => shell_quote("$Bin/module"), IGNORE => $opt->get_ignore });

my $hgdefault;
if (-f "$rootdir/par/hgbuild") {
	open my $fin, "$rootdir/par/hgbuild" or die "Cannot open hgbuild";
	$hgdefault = <$fin>;
	chomp($hgdefault);
}

my %hgsupport = (hg19 => 1, b37 => 1, hg38 => 1, b38 => 1);
while(my $hgbuild = prompt("Human genome build: ", -default => $hgdefault // "b37")) {
#my $hgbuild = "b37";
	unless (defined $hgsupport{$hgbuild}) {
		print STDERR "$hgbuild is not supported\n";
		next;
	}
	$conf{GENOME}{HGBUILD} = $hgbuild;
	$conf{GENOME}{CHRX} = $hgbuild =~ /^hg/ ? "chrX" : "X";
	$conf{GENOME}{CHRY} = $hgbuild =~ /^hg/ ? "chrY" : "Y";
	$conf{GENOME}{PARX1} = $PAR{$hgbuild}{X1};
	$conf{GENOME}{PARX2} = $PAR{$hgbuild}{X2};
	$conf{GENOME}{PARY2} = $PAR{$hgbuild}{Y2};
	last;
}

# Write hgbuild to param file
{
	open my $fout, ">$rootdir/par/hgbuild" or die "Cannot write to hgbuild";
	print $fout $conf{GENOME}{HGBUILD}, "\n";
}


# Instantiate the workflow
$wkf->inst(\%conf);

# QC Parameters for converting genotypes
print "* Parameters for Converting Genotype\n";
my @params =
 (	SNVTR => "SNV Tranche Threshold" => 99.9, 
	MAF => "Min. Minor Allele Freq" => 0.01,
    VQ 	=> "Variant Quality for Autosome" => 100,
    GQ  => "Min. GQ for Autosome" => 60,
    MISS => "Max. Missing Rate for Autosome" => 0.02,
    HWE => "P-value Threshold for HWE-test for Autosome" => 1.0e-5,
    VQX => "Variant Quality for chrX" => 100,
    GQX => "Min. GQ for chrX" => 40,
    MISSX => "Max. Missing Rate for chrX" => 0.03,
    HWEX => "P-value Threshold for HWE for chrX" => 1.0e-6);

my $params = read_pars(\@params, "convert_params");
while(my $covpar = prompt_params($params)) {
	write_pars($covpar, "convert_params");
	$wkf->run({ tasks => "VCF2BPed", dryrun => $opt->get_dryrun, interact => 1});
	last if $opt->get_dryrun;
	print "Please check the expected output:\n";
	print join("\n", map { "\t".$_ } $wkf->get_expected("VCF2BPed")), "\n";
	my $rerun = prompt("Re-run convertion? ", -yn, -default => 'n');
	if ($rerun) {
		foreach my $outfile ($wkf->get_expected("VCF2BPed")) {
			unlink $outfile if -f $outfile;
			#last;
		}
	}
	else {
		last;
	}
}

@params = (THIN => "Distance (bp) between two sites in VCF for thinning" => 0);
$params = read_pars(\@params, "depth_params");
unless ($opt->get_reuse) {
	while(my $depthpar = prompt_params($params)) {
		write_pars($depthpar, "depth_params");
		$wkf->run({ tasks => "SiteDepth", dryrun => $opt->get_dryrun, interact => 1 });		
		last if $opt->get_dryrun;
		print "Please check the expected output:\n";
		print join("\n", map { "\t".$_ } $wkf->get_expected("SiteDepth")), "\n";
		my $rerun = prompt("Re-run SiteDepth? ", -yn, -default => 'n');
		if ($rerun) {
			foreach my $outfile ($wkf->get_expected("SiteDepth")) {
				unlink $outfile if -f $outfile;
			}
		}
		else {
			last;
		}
	}
}

print "* Options for running Peddy\n";
@params = 
( 	NT => "Number of parallel threads (reduce the number for large data set)" => 4,
	PLOT  => "Plot Peddy's output (disable this option for large data set)" => "--plot",
	SITES => "Customize the sites for PCA and relatedness calculation" => 
			 $conf{GENOME}{HGBUILD} eq 'b38' || $conf{GENOME}{HGBUILD} eq 'hg38' ? "--sites hg38" : "" );
$params = read_pars(\@params, "peddy_opts");
while(my $peddyopt = prompt_params($params)) {
	write_pars($peddyopt, "peddy_opts");
	$wkf->run({ tasks => "Peddy", dryrun => $opt->get_dryrun, interact => 1 });
	last if $opt->get_dryrun;
	print "Please check the expected output\n";
	print join("\n", map { "\t".$_ } $wkf->get_expected("Peddy")), "\n";
	my $rerun = prompt("Re-run Peddy? ", -yn, -default => 'n');
	if ($rerun) {
		foreach my $outfile ($wkf->get_expected("Peddy")) {
			unlink $outfile if -f $outfile;
		}
	}
	else {
		last;
	}
}

# Run GenoQC and King for relatedness test using king
print "### Calculating sample QC and relatedness statistics\n";
foreach my $taskname (qw|GenoQC Kinship|) {
	if (all { -f $_ } $wkf->get_expected($taskname)) {
		my $rerun = prompt("Re-run $taskname? ", -yn, -default => 'n');
		if ($rerun) {
			foreach my $outfile ($wkf->get_expected($taskname)) {
				unlink $outfile if -f $outfile;
			}
		}
	}
	$wkf->run({ tasks => $taskname, dryrun => $opt->get_dryrun, interact => 1 });
}


print "* Parameters for Sample QC\n";
@params = 
( 	MINDP => "Minimal allowable mean depth at heterozygotes" => 10,
	MINMYDP => "Minimal chrY norm depth for males (if 0, using peddy's prediction)" => 0.2,
	KINTOOL  => "Which software do we use to collect output statistics?" => "king",
	IDSEP => "Separator between sample IDs used in creating ID for a pair" => "_",
	MINDUP => "Minimal Kinship Coeff for Duplicated Samples" => 0.45,
	MIN1D => "Minimal Kinship Coeff for 1st-deg Relatives" => 0.177,
	MAX1D => "Maximal Kinship Coeff for 1st-deg Relatives" => 0.354,
	POIBS0 => "Maximal IBS0 Proportion for Parent-Offspring" => 0.005, 
	CRPREL => "Minimal Kinship Coeff for Defining Cryptic Relatedness" => 0.0884,
	NOPLOT => "Do not generate relatedness plots for large data set" => 0);

$params = read_pars(\@params, "relpair_params");
while(my $kinpar = prompt_params($params)) {
	$kinpar->{KINTOOL} = lc($kinpar->{KINTOOL});
	unless($kinpar->{KINTOOL} eq 'king' || $kinpar->{KINTOOL} eq 'peddy') {
		print STDERR "Only support output from king or peddy\n";
		next;
	}
	write_pars($kinpar, "relpair_params");
	$wkf->run({ tasks => "Plot", dryrun => $opt->get_dryrun, interact => 1 });
	last if $opt->get_dryrun;
	print "Please check the expected output\n";
	print join("\n", map { "\t".$_ } $wkf->get_expected("Plot")), "\n";
	my $rerun = prompt("Re-run plotting? ", -yn, -default => 'n');
	if ($rerun) {
		foreach my $outfile ($wkf->get_expected("Plot")) {
			unlink $outfile if -f $outfile;
		}
	}
	else {
		last;
	}
}

$wkf->run({ tasks => "Patch", dryrun => $opt->get_dryrun });

sub write_pars {
	my ($param, $file) = @_;
	my $fout = IO::File->new($wkf->get_subdir("par")."/$file", "w") or croak "Cannot write $file";
	while(my ($key, $val) = each %$param) {
		if (looks_like_number($val)) {
			print $fout $key."=".$val, "\n";
		}
		else {
			print $fout $key."='".$val, "'\n";
		}
	}
}

sub read_pars {
	my ($params, $file) = @_;
	croak "The length of input params list must be a multiple of 3" if scalar(@$params)%3 > 0;
	my $parfile = $wkf->get_subdir("par")."/$file";
	if (-f $parfile) {
		my %pars;
		my $fin = IO::File->new($parfile);
		while(<$fin>) {
			chomp();
			my ($key, $val) = split('=', $_);
			$val =~ s/^'//; $val =~ s/'$//;
			$pars{$key} = $val;
		}
		for(my $ii = 0; $ii < int(@$params/3); $ii ++) {
			croak "Parameter $params->[3*$ii] cannot be found" unless defined $pars{$params->[3*$ii]};
			$params->[3*$ii+2] = $pars{$params->[3*$ii]};
		}
	}
	return $params;
}

## Files and scripts

# Calculate depth of coverage on autosomes and chrX/Y
sub depth_script {
	return <<'EOF';

source _PARDIR_/depth_params

vcftools --gzvcf _VCF_ --not-chr _GENOME.CHRX_ --not-chr _GENOME.CHRY_ --thin $THIN --depth --out _WRKDIR_/site.auto 

vcftools --gzvcf _VCF_ --chr _GENOME.CHRY_ --from-bp _GENOME.PARX1_ --to-bp _GENOME.PARY2_ --thin $THIN --depth --out _WRKDIR_/site.chrY

vcftools --gzvcf _VCF_ --chr _GENOME.CHRX_ --from-bp _GENOME.PARX1_ --to-bp _GENOME.PARX2_ --thin $THIN --depth --out _WRKDIR_/site.chrX

EOF
}


sub convert_script {
	return <<'EOF';

source _PARDIR_/convert_params

vcf2bped.pl --vcf _VCF_ --ped _PED_ --hg _GENOME.HGBUILD_ --out _WRKDIR_/geno --clean --snp-only --vqsr-snp $SNVTR --maf $MAF \
	--vq $VQ --gq $GQ  --geno $MISS --hwe $HWE --vq-chrX $VQX --gq-chrX $GQX --geno-chrX $MISSX --hwe-chrX $HWEX

EOF
}

sub check_geno_files {
	my (@exp) = @_;
	if (all { -f "$rootdir/$_" } @exp) {
		my $bimfile = "$rootdir/$exp[0]";
		open my $fin, $bimfile or die "Cannot open $bimfile";
		my ($n_auto, $n_x, $n_y) = (0) x 3;
		while(<$fin>) {
			my $chr = (split)[0];
			if ($chr <= 22) {
				$n_auto ++;
			}
			elsif ($chr == 23) {
				$n_x ++;
			}
			elsif ($chr == 24) {
				$n_y ++;
			}
		}
		print "Number of markers of autosome and chrX: $n_auto, $n_x\n";
		if ($n_auto > 0 && $n_x > 0) {
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

sub peddy_script {
	return <<'EOF';

source _PARDIR_/peddy_opts

peddy -p $NT $PLOT $SITES --prefix _WRKDIR_/peddy _VCF_ _PED_

if [[ -f _WRKDIR_/peddy.background_pca.json ]]; then
	in2csv _WRKDIR_/peddy.background_pca.json > _WRKDIR_/peddy.background_pca.csv
fi

EOF
}

sub peddy_files {
# ped_check.rel-difference.csv may not always be found in the output
	return map { "wrk/peddy.$_" } qw|
background_pca.json 
background_pca.csv
het_check.csv
html
ped_check.csv
peddy.ped
sex_check.csv|;
}

sub genoqc_script {
	return <<'EOF';

plink --bfile _WRKDIR_/geno --missing --out _WRKDIR_/geno

#plink --bfile _WRKDIR_/geno.auto.pass1 --missing --out _WRKDIR_/geno.auto.pass1

plink --bfile _WRKDIR_/geno --check-sex --out _WRKDIR_/geno

plink --bfile _WRKDIR_/geno --het --out _WRKDIR_/geno

plink --bfile _WRKDIR_/geno --ibc --out _WRKDIR_/geno

plink --bfile _WRKDIR_/geno --mendel-duos --mendel --out _WRKDIR_/geno

plink --bfile _WRKDIR_/geno --homozyg --out _WRKDIR_/geno

perl _MODULE_/merge_samp_qcstats.pl _PED_ _WRKDIR_ _OUTDIR_

EOF
}


sub genoqc_files {
	return ("out/sample.csv",
 		map { "wrk/geno.$_" } qw|
het
ibc
imiss
imendel
sexcheck
hom.indiv
|);
}

sub kinship_script {
	return <<'EOF';

plink --bfile _WRKDIR_/geno --genome --out _WRKDIR_/geno

gzip -f _WRKDIR_/geno.genome

king -b _WRKDIR_/geno.bed --kinship --related --ibs --prefix _WRKDIR_/geno

gzip -f _WRKDIR_/geno.ibs0 _WRKDIR_/geno.kin0

perl _MODULE_/merge_kinship_est.pl _PED_ _WRKDIR_ _OUTDIR_ _IGNORE_

EOF
}

sub kinship_files {
	return ("out/relpairs.csv",
		map { "wrk/geno.$_" } qw|
ibs
ibs0.gz
kin
kin0.gz
genome.gz
|);
}

sub plot_script{
	return 'Rscript _MODULE_/plot_samp_stats.R _OUTDIR_ _WRKDIR_ _PARDIR_/relpair_params';
}

sub plot_files {
	my @files;
	my $fs = IO::File->new("$Bin/module/plot_samp_stats.R") or croak "Cannot open plot script";
	while(<$fs>) {
		if (m{/([\w\-]+\.png)}) {
			push @files, $1;
		}
	}
	return map { "out/$_" } @files;
}

sub patch_script {
	return 'perl _MODULE_/create_patch_files.pl _PED_ _OUTDIR_';
}



