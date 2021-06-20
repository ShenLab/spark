#!/usr/bin/env perl

use strict;
use warnings;
use Carp;
use FindBin qw|$Bin|;
use POSIX qw|ceil|;
use Getopt::Lucid qw|:all|;
use IO::Prompt;
use Data::Dumper;
use Cwd qw|abs_path|;
use List::MoreUtils qw|all any uniq|;
use Perl6::Slurp;
use String::ShellQuote;
#use Utils::File::Iter qw|iter_file|;
use Utils::Workflow;
use Utils::Hash qw|prompt_params|;
#use Utils::File qw|count_line|;
use Genome::UCSC qw|%PAR|;
use Genet::Ped;


#############################
## Command line interface  ##
#############################
# The default input is plink bped format genotypes and a output directory
my @spec = (
	Param("geno|i")->valid(sub { -r "$_.bed" && -r "$_.bim" && -f "$_.fam" }),
	Param("outdir|out"),
	Param("prefix"),
	Switch("help|h"),
	Switch("dryrun|dry"),
	Switch("force")
	);

my $opt = Getopt::Lucid->getopt(\@spec);

if ($opt->get_help) {
	print STDERR <<EOF;
Purpose:
	This is an interactive script for sample and SNP level QC of SNP genotyping data.
	Perform standard sample and SNP level QC of SNP array data.

Usage: 
	samp_snp_qc.pl --geno PREFIX [--prefix PREFIX] --out OUTDIR

Options:
	--force: Over-write existing output directory.
	--prefix: Specify the prefix for cleaned genotype data output (Default: GenoHQ).

Notes: 
	The aim of sample and SNP level QC is to produce a final cleaned SNP genotype data ready for
	genetic association analysis. The script checks sample sex with sex chromosome genotypes. If family 
	members are included in genotyping, it also verifies if the pairwise kinship coefficient is
	consistent with estimated kinship coefficient. Problematic samples with low call rates or unusual 
	heterozygosity rates will be removed as bad samples. Samples with unknown sex, incorrect familial 
	relationships, or cryptic relatedness will be included in the final cleaned data set but not used
	for calculating SNP level QC statistics used for SNP inclusion/exclusion. 

	QC for SNP array genotype data is an iterative process on SNP and sample levels. The script will 
	run following major steps in order. Each step will display an interactive menu of adjustable parameters.
	After each run, users can examine the output files and decide if the step needs re-run after
	changing parameters.

	1. Sex check.
		This step examines sample genotype call rate and check sex based on sex chromosome genotypes.
		It is recommended to remove samples with very low call rate before running this script, because
		they will confound the genotype-based sex check at first step. 

		The following figures are generated:
		pass_1_miss.png - Cumulative distribution plot of genotype call rate per sample, 
						  highlighting samples with low call rate.
		pass_1_het.png - Scatterplot of estimated inbreeding coefficient and missing genotype rate.
						 Sample with unusual low rate of inbreeding may suggest possible contamination.
						 See also array_qc.pl for detecting contamination from BAF signals. 
		pass_1_chrX_FvsMiss.png - Scatterplot of chrX inbreeding coeff vs. missing genotype rate.
						 Thresholds of chrX inbreeding coeff is used to infer genotype sex. Samples with 
						 inconsistent genotype and pedigree sex are highlighted.  
		pass_1_chrXY_FXvsYcnt.png - (Available if chrY SNPs exist) Scatterplot of chrX inbreeding coeff
						 vs. chrY genotype counts. When chrX inbreeding level and chrY geno counts predict
						 different sex, it is a sign of likely have chrX abnormality. For QC steps below,
						 we use genotype sex predicted by chrX inbreeding level.

	2. Relatedness check.
		This step will first perform automated SNP level QC filter after excluding samples with low call 
		rate and unknown sex in the previous step to create an intermediate cleaned of genotype data.  
		Then it runs king on cleaned genotype with "--related" options to estimate kinship coeff and regions
		identical by descent (IBD) for all pairs of samples within families and identify likely cryptic related 
		samples between families. Note that king's IBD segments calling algorithm requires high density 
		SNP genotypes and will become less accurate for genotype data from targted sequencing. 

		To create a cleaned data file for SNP level QC, culprit sample in family causing inconsistent pairwise
		relationships with pedigree file will be identified and removed, and for duplications and cryptic related
		samples, only one with the highest call rate will be kept. Samples removed at this step can be added back 
		andhave their familial relationship manually fixed at later stage before generating the final cleaned data. 

		The following two figures are generated:
		pass_2.relpairs.png -- Scatterplot of kinship coefficient vs IBS0 by king.
		pass_2.relpairs.2.png -- Scatterplot of Prob(IBD=1) vs Prob(IBD=2) by king.
			Both plots use color and shape to indicate inferred and pedigree relationships, and highlight 
			sample pairs with likely relationship error.

	3. SNP level QC filter.
		This step generates a clean data after removing error samples identified from the previous step,
		then generate SNP QC statistics file from cleaned sample genotypes, and determine the appropriate 
		filtering thresholds for marker inclusion. 

		The distribution of QC statistics are displayed in the following plots for determining QC thresholds:
		pass_3_comm_miss.png, pass_3_rare_miss.png
			- Cumulative distribution of missing rate for common or rare SNPs			
		pass_3_hwe_qq.png - QQ plot of HWE p-values against null expectation
		pass_3_hist_hhcount.png - Histogram of haploid heterozygotes counts for SNPs on haploid chromosome
		pass_3_hist_duperr.png - (Available if duplicated samples exist) Histogram of duplicate inconsistency count.
		pass_3_hist_mendel.png - (Available if family trios or duos exist) Histogram of Mendelian inheritance error.

	4. Sex and relationship check again.
		This step generate cleaned genotype data based on SNP level filters derived from the previous step, then re-run
		the sex and relationship checks	to confirm previous results.

		The following figures are generated, they are similar to the figures generated in Step 1 and 2. 
		pass_4_HetvsMiss.png, pass_4_chrX_FXvsMiss.png, pass_4_chrX_F_XvsAUTO.png, pass_4_chrXY_FXvsYcnt.png,
		pass_4.relpairs.png and pass_4.relpairs.2.png.

	5. Patch genotype files manually
		If sample errors are found, users should decide how to fix them by providing following patch files:
		
		 patch.badsamps.txt - Samples to be removed, format: IID
		 patch.dups.txt - Error duplications, format: IID1 [IID2] IID3 ...
		 				  The target sample that was duplicated should appear first, the sample that should 
		 				  be kept as the main data is shown in square bracket.
		 patch.rels.txt - Cryptic related samples, format: IID1 [IID2] IID3 ...
		 				  Same format as duplicates. All samples within the same line are related with each other, 
		 				  but do not reflected in the original PED file.
 		 patch.nonpar.txt:	Parent update, format IID1 DAD/MOM IID2.
 		 					IID2 is optional, if IID2 is not provided, the original parent (DAD or MOM) for IID1 
 		 					will be set to missing; if IID2 is provided, it will be used to replace original parent.   
		 patch.sex_update.txt:	Sex update, format: IID GENDER
		 patch.swapped.txt:	Samples that are swapped, format: IID1 IID2. They will be swapped back.

	6. Final check.
		The final cleaned genotype data will be in out directory (GenoHQ.fam/bim/bed). This step will run sex
		and relationship check again on patched genotype data to confirm patching is correct.
		The following figures are generated, they are similar to the figures generated in Step 1 and 2. 
		GenoHQ_HetvsMiss.png, GenoHQ_chrX_FXvsMiss.png, GenoHQ_chrX_F_XvsAUTO.png, GenoHQ_chrXY_FXvsYcnt.png,
		GenoHQ.relpairs.png and GenoHQ.relpairs.2.png. GenoHQ can be replaced by the prefix provided by --prefix option.

		Optionally, a quick genome-wide association (GWA) scan is performed and qq-plot of genome-wide p-values 
		is created for inspecting if there is global inflation in test statistics. Depending on the data type, 
		either a case-control or family-based test is performed. The affection status in the .fam file is used
		as phenotype for GWA scan. 

Dependencies:
	plink (v2), king, R (ggplot2, plotrix, ggrepel, qqman)
	plink-1.07 is required if family-based association test using --dfam option is used for GWA scan.

EOF
	exit 1;
}

$opt->validate({requires => [qw|geno outdir|]});
my $outprefix;
if ($opt->get_prefix) {
	$outprefix = $opt->get_prefix;
}
else {
	$outprefix = "GenoHQ";
}

###########################
## Input file validation ##
###########################
# Check if input ped file contain families and if chrY markers can be found
my $input = $opt->get_geno;
my ($famflag, $chrYflag, $ccflag);
{
	my %samp = map { my ($fid, $iid) = (split)[0, 1]; $fid.$;.$iid => 1 } slurp "$input.fam";
	open my $fin, "$input.fam" or die "Cannot open $input.fam file";
	my %pheno;
	while(<$fin>) {
		my ($fid, $iid, $dad, $mom, $sex, $phe) = split;
		$pheno{$phe} ++ if $phe == 1 || $phe == 2;
		if ($dad ne '0' && defined $samp{$fid,$dad} ||
		    $mom ne '0' && defined $samp{$fid,$mom}) {
			$famflag = 1;
		}
	}
	if (keys %pheno == 2) {
		$ccflag = 1;
	}
	my %chr;
	open my $fbim, "$input.bim" or die "Cannot open $input.bim file";
	while(<$fbim>) {
		my $chr = (split)[0];
		$chr{$chr} = 1;
	}
	unless(defined $chr{"23"} || defined $chr{X}) {
		croak "Cannot find chrX in the genotype";
	}
	if (defined $chr{"24"} || $chr{Y}) {
		$chrYflag = 1;
	}
}
if ($famflag) {
	print STDERR "The genotype data are from family-based sample.\n";
} 
else {
	print STDERR "The genotype data are from population-based sample.\n";
}
if ($chrYflag) {
	print STDERR "chrY markers are found\n";
} 
else {
	print STDERR "No chrY marker is found\n";
}
if ($ccflag) {
	print STDERR "Sample are composed of both cases and controls\n";
} 
else {
	print STDERR "All samples have the same phenotype\n";
}

####################
# Create workflow ##
####################
my $rootdir = $opt->get_outdir;
my $wkf = Utils::Workflow->new($rootdir, {engine => "BASH", force => $opt->get_force});

$wkf->add(first_pass_run($chrYflag),  { name => "SampQC",  expect => [ first_pass_files($chrYflag) ] })
	->add(first_pass_plot(), { name => "SampQCPlot", expect => [ first_pass_figs($chrYflag) ] })
	->add(second_pass_run(), { name => "Kinship", expect => [ second_pass_files($famflag) ]})
	->add(third_pass_run(), { name => "SNPQC" , expect => "wrk/pass_3.snp_qc_stats.txt" })
	->add(third_pass_filter(), { name => "SNPQCFilter", expect => [ third_pass_files() ] })
	->add(fourth_pass_run($chrYflag), { name => "FinalSampQC", expect => [ fourth_pass_files($chrYflag, $famflag) ] })
	->add(fifth_pass_run($chrYflag, $ccflag, $famflag), 
		{ name => "FinalPlots", expect => [ fifth_pass_files($chrYflag, $famflag) ] });

my %conf = (default => { MODULE => shell_quote("$Bin/module"), 
						 UTILS  => shell_quote("$Bin/utils"), 
						 RAWGENO => abs_path($input), PREFIX => $outprefix });
$wkf->inst(\%conf);

# Optional files in each step
my %optional = ( SampQCPlot => 
	[map { "wrk/pass_1$_"} qw|_miss.badsamps.csv _het.badsamps.csv _chrX.sexerr.csv 
							.remove .sexupdate .rename a.bim a.bed a.fam|],
	Kinship => 
	[map { "wrk/pass_2.$_" } qw|kin0 repairs.remove repairs.dups.txt repairs.cyprel.txt|],
	FinalSampQC =>
	[map { "wrk/pass_4.$_" } qw|badsamp.remove badsamp.csv sexna.remove sexerr.csv
								sexupdate relpairs.remove relpairs.dups.txt relpairs.cyprel.txt|],
	FinalPlots =>
	[map { "out/$outprefix." } qw|kin0 badsamp.remove badsamp.csv sexna.remove sexerr.csv
								sexupdate relpairs.remove relpairs.dups.txt relpairs.cyprel.txt| ]
	);
if ($ccflag) {
	if ($famflag) {
		push @{$optional{FinalPlots}}, map { "out/$outprefix." } qw|dfam dfam.qq.png dfam.mth.png|;
	}
	else {
		push @{$optional{FinalPlots}}, map { "out/$outprefix." } qw|assoc assoc.qq.png assoc.mth.png|;
	}
}


###############################
# Run workflow interactively ##
###############################
# Pass 1: remove bad samples and check sex
print "\n* Generate sample quality control statistics (SampQC)\n";
if (all { -f $_ } $wkf->get_expected("SampQC")) {
	my $rerun = prompt("Re-run SampQC? ", -yn, -default => 'n');
	if ($rerun) {
		foreach my $outfile ($wkf->get_expected("SampQC")) {
			unlink $outfile if -f $outfile;
		}
	}
}
$wkf->run({ tasks => "SampQC", dryrun => $opt->get_dryrun, interact => 1 });


# QC Parameters for filtering bad samples and determine SNP sex
print "\n* Parameters for first pass sample QC (SampQCPlot)\n";
my @params =
 (	MISSMAX => "Max. allowable missing genotype rate on unfiltered SNPs" => 0.1, 
	WGFMIN 	=> "Min. autosome inbreeding coeff using raw SNP genotypes (for detecting contamination)" => -0.3,
    FFMAX	=> "Max. chrX inbreeding coeff for female using raw SNP genotypes" => 0.3,
   	MFMIN  	=> "Min. chrX inbreeding coeff for male using raw SNP genotyes" => 0.7);
interactive_run(\@params, "pass_1", "SampQCPlot");

# Pass 2: relationship check, before that we need to clean SNPs
# Parameters for quick SNP QC and relationshp tests
print "\n* Parameters for quick SNP QC and relationship test (Kinship)\n";
@params =
 (	MINMAF 	=> "Min. minor allele frequency" => 0.02,
 	MAXGMIS => "Max. missing genotype rate" => 0.01,
 	MINHWE 	=> "Min. HWE p-value" => 0.000001,
 	MINDUP	=> "Min. kinship coeff for duplicated sample or MZ twin" => 0.45,
 	MIN1D 	=> "Min. kinship coeff for first-degree relative pair" => 0.177,
 	MAX1D 	=> "Max. kinship coeff for first-degree relative pair" => 0.354,
 	POIBS0 	=> "Max. IBS0 for parent-offspring pair" => 0.01,
 	CRPREL 	=> "Min. kinship coeff for defining cryptic related pairs" => 0.17,
 	IGNORE  => "ID suffix for duplicated samples" => '_Re');
interactive_run(\@params, "pass_2", "Kinship");

# Pass 3: SNP QC filter
if (all { -f $_ } $wkf->get_expected("SNPQC")) {
	my $rerun = prompt("Re-run SNPQC? ", -yn, -default => 'n');
	if ($rerun) {
		foreach my $outfile ($wkf->get_expected("SNPQC")) {
			unlink $outfile if -f $outfile;
		}
	}
}
$wkf->run({ tasks => "SNPQC", dryrun => $opt->get_dryrun, interact => 1 });

print "\n* Parameters of QC filters to generate final high-quality SNPs (SNPQCFilter)\n";
# Check how many sample left in the genotype data
my $ct = count_samp("$rootdir/wrk/pass_2.fam", "$rootdir/wrk/pass_2.relpairs.remove");
@params = 
  (	MINMAF 	=> "Min. minor allele frequency" => 0.01,
  	RAREMAF	=> "Max. MAF to define rare variants" => 0.05,
  	MAXRMIS	=> "Max. missing genotype rate for rare variants" => 0.01,
  	MAXCMIS => "Max. missing genotype rate for common variants" => 0.03,
  	CTRLHWE => "Use HWE in controls only?" => "FALSE",
 	MINHWE 	=> "Min. HWE p-value" => 0.00001); 
if ($ct->{Male} > 0) {
	push @params, "MAXHHCNT" => "Max. haploid heterozygotes on sex chromosome (#Males=$ct->{Male})" 
		=> ceil($ct->{Male}*0.025+0.5);
}
if ($ct->{Trio}+$ct->{Duo} > 0) {
	push @params, "MAXMENDEL" => "Max. Mendel error (#Trios=$ct->{Trio}, #Duos=$ct->{Duo})" 
		=> ceil(($ct->{Trio}+$ct->{Duo})*0.025+0.5);
}
if (-f "$rootdir/wrk/pass_2.relpairs.dups") {
	my $ndups = count_dups("$rootdir/wrk/pass_2.relpairs.dups");
	push @params, "MAXDUPERR"=> "Max. duplication discordance (#Dups=$ndups)"
		=> ceil($ndups*0.05+0.5);
}
interactive_run(\@params, "pass_3", "SNPQCFilter");

# Pass 4: Generate sample summary statistics based on HQ SNP genotypes
print "\n* Parameters for final round of sample QC\n";
@params =
 (  MISSMAX => "Max. allowable missing genotype rate on HQ SNPs" => 0.1, 
	WGFMIN 	=> "Min. autosome inbreeding coeff" => -0.2,
	WGFMAX 	=> "Max. autosome inbredding coeff" => 0.2,
 	FFMAX	=> "Max. chrX inbreeding coeff for female using raw SNP genotypes" => 0.3,
   	MFMIN  	=> "Min. chrX inbreeding coeff for male using raw SNP genotyes" => 0.7, 
 	MINDUP	=> "Min. kinship coeff for duplicated sample or MZ twin" => 0.45,
 	MIN1D 	=> "Min. kinship coeff for first-degree relative pair" => 0.177,
 	MAX1D 	=> "Max. kinship coeff for first-degree relative pair" => 0.354,
 	POIBS0 	=> "Max. IBS0 for parent-offspring pair" => 0.01,
 	CRPREL 	=> "Min. kinship coeff for defining cryptic related pairs" => 0.0884,
 	IGNORE  => "ID suffix for duplicated samples" => '_Re',
 	SEXIMP 	=> "Remove sample with unknown sex after imputation?" => 'FALSE' );
interactive_run(\@params, "pass_4", "FinalSampQC");


# Pass 5: Patch the genotype data based on output from Pass 4 QC
# This need some manual review on the results, and should be iterative.
# (optionally) perform association test to evaluate inflation

my %stdsuffs  = ('badsamps.txt' => 'Samples to be removed', 
		 		'swapped.txt' => 'Samples that are swapped', 
		 		'sex_update.txt' => 'Sex update',
				'nonpar.txt' => 'Parents update', 
				'dups.txt' => 'Error duplications', 
				'rels.txt' => 'Cryptic related samples');
# We will convert those patch files to plink format

my $patchpref = "$rootdir/wrk/patch";
print "* Please review the results from previous step and write patch files to $patchpref\n";
print join("\n", map { ".$_:\t$stdsuffs{$_}" } sort keys %stdsuffs), "\n";
exit if $opt->get_dryrun();

prompt("* Continue after patch files are ready", -tty);

while(1) {
	foreach my $suf (sort keys %stdsuffs) {
		my $desc = $stdsuffs{$suf};
		if (-f "$patchpref.$suf") {
			print "Found $patchpref.$suf\n";
		}
		else {
			print "$desc ($suf) is not found\n";
		}
	}
	if (prompt("Confirm the patch file list?", -yn, -default => 'n')) {
		my $outfile = translate_patch($patchpref, "$rootdir/wrk/pass_5");
		print "* Please review the output $outfile and prepare final clean data file to $rootdir/out/$outprefix\n";
		if (prompt("Skip manual patching?", -default => "n", -yn)) {
			foreach my $suf (qw|bim bed fam|) {
				link("$outfile.$suf", "$rootdir/out/$outprefix.$suf");
			}
		}
		else {
			while(1) {
				prompt("Continue after final genotype data are manually prepared", -tty);
				if (all { -f "$rootdir/out/$outprefix.$_" } qw|bim bed fam|) {
					print "Final genotype data found\n";
					last;
				}
				else {
					print "Cannot find final genotype data: $rootdir/out/$outprefix\n";
				}
			}
		}
		print "\n* Parameters for making plots on final genotype data\n";
		@params =
		(  MISSMAX => "Max. allowable missing genotype rate on HQ SNPs" => 0.1, 
			WGFMIN 	=> "Min. autosome inbreeding coeff" => -0.2,
			WGFMAX 	=> "Max. autosome inbredding coeff" => 0.2,
			FFMAX	=> "Max. chrX inbreeding coeff for female using raw SNP genotypes" => 0.3,
			MFMIN  	=> "Min. chrX inbreeding coeff for male using raw SNP genotyes" => 0.7, 
			MINDUP	=> "Min. kinship coeff for duplicated sample or MZ twin" => 0.45,
			MIN1D 	=> "Min. kinship coeff for first-degree relative pair" => 0.177,
			MAX1D 	=> "Max. kinship coeff for first-degree relative pair" => 0.354,
			POIBS0 	=> "Max. IBS0 for parent-offspring pair" => 0.01,
			CRPREL 	=> "Min. kinship coeff for defining cryptic related pairs" => 0.0884,
			IGNORE  => "ID suffix for duplicated samples" => '_Re',
			GWAS	=> "Perform genome-wide association test?" => 'FALSE' );
		interactive_run(\@params, "pass_5", "FinalPlots");
		if (prompt("Exit QC?", -yn, -default => "y")) {
			last;
		}
	}
	else {
		print "Continue to curate the patch file ...\n";
	}
}

sub translate_patch {
	my ($prefix, $geno) = @_;

	# Please note the order of execution:

	# 1. First remove badsamp and err dups 

	# 2. Update IDs for swapped samples, and update FIDs for crypt related samples

	# 3. Update sex, those already swapped or removed will be skipped, those have changed IDs due to
	# family merging or should be specified with new FIDs.

	# 4. Update parents: those already swapped will be skipped, those changed IDs due to
	# family merging should be specified with new FIDs.

	# Note: currently we only merge family ID for cryptic-related samples
	# If family merging happens, user will be prompted for family name separator
	# If one of cryp-rel sample is to be removed, they can be included in dups file
	# If cryp-rel showed PO relationships, they can be specified in parents update file.

	my $infile = "$rootdir/wrk/pass_4";
	my $outfile = "$geno.step1";
	
	my %famid = map { (split)[1,0] } slurp "$infile.fam";	

	# Step 1:
	my %remove;
	if (-f "$prefix.badsamps.txt" || -f "$prefix.dups.txt") {
		open my $frm, ">$outfile.remove" or die "Cannot write to $outfile.remove";
		if (-f "$prefix.badsamps.txt") {
			open my $fin, "$prefix.badsamps.txt" or die "Cannot open badsamps.txt";
			while(<$fin>) {
				next if /^#/ || /^\s+$/;
				my $iid = (split)[0];
				die "Cannot find family ID for $iid" unless defined $famid{$iid};
				$remove{$iid} = 1;
				print $frm $famid{$iid}, "\t", $iid, "\n";
			}
		}
		if (-f "$prefix.dups.txt") {
			open my $fin, "$prefix.dups.txt" or die "Cannot open dups.txt";
			open my $fup, ">$outfile.update-ids" or die "Cannot write to $outfile.update-ids";
			while(<$fin>) {
				next if /^#/ || /^\s+$/;
				my @iids = split;
				my $keepid;
				unless(any { /^\[[^\[\]]+\]$/ } @iids) {
					print STDERR join("\t", @iids), "\n";
					die "Must specify the sample to retain by [IID]";
				}
				else {
					my @keepid = grep { /^\[[^\[\]]+\]$/ } @iids;
					unless (@keepid == 1) {
						die "Can only specify one sample to retain";
					}
					$keepid = $keepid[0];
					$keepid =~ s/^\[//; $keepid =~ s/\]$//;
					unless (defined $famid{$keepid}) {
						die "Cannot find family ID for kept sample: $keepid"
					}
				}
				my $dupid = shift @iids;
				$dupid =~ s/^\[//; $dupid =~ s/\]$//;
				foreach my $iid (@iids) {
					$iid =~ s/^\[//; $iid =~ s/\]$//;
					$remove{$iid} = 1;
					unless (defined $famid{$iid}) {
						die "Cannot find family ID for dup sample: $iid";
					}
					print $frm $famid{$iid}, "\t", $iid, "\n";
				}
				# Note: remove were perform after rename
				if ($keepid ne $dupid) {
					print $fup join("\t", $famid{$keepid}, $keepid, $famid{$dupid}, $dupid), "\n";
					print $fup join("\t", $famid{$dupid}, $dupid, $famid{$keepid}, $keepid), "\n";
				}
			}
			close $fup;
		}
		close $frm;
		if (-f "$outfile.update-ids") {
			system("plink --bfile $infile --update-ids $outfile.update-ids --remove $outfile.remove --make-bed --out $outfile");
		}
		else {
			system("plink --bfile $infile --remove $outfile.remove --make-bed --out $outfile");
		}
		$infile = $outfile;
	}
	else {
		unlink(glob("$outfile.*"));
		$outfile = $infile;
	}

	# Step 2:
	$outfile = "$geno.step2";
	my (%swap);
	if (-f "$prefix.swapped.txt" || -f "$prefix.rels.txt") {
		open my $fup, ">$outfile.update-ids" or die "Cannot write to $outfile.update-ids";
		if (-f "$prefix.swapped.txt") {
			open my $fin, "$prefix.swapped.txt" or die "Cannot open swapped.txt";
			while(<$fin>) {
				next if /^#/ || /^\s+$/;
				my ($iid1, $iid2) = (split)[0,1];
				$swap{$iid1} = 1; 
				$swap{$iid2} = 1;
				die "Cannot find family ID for $iid1" unless defined $famid{$iid1};
				die "Cannot find family ID for $iid2" unless defined $famid{$iid2};
				print $fup join("\t", $famid{$iid1}, $iid1, $famid{$iid2}, $iid2), "\n";
				print $fup join("\t", $famid{$iid2}, $iid2, $famid{$iid1}, $iid1), "\n";
			}
		}
		if (-f "$prefix.rels.txt") {
			my $concat = prompt("Please specify string to concat family IDs (default '_'): ", -default => "_");
			# Identify all samples in the family that are not removed
			# and update their family IDs
			open my $fin, "$prefix.rels.txt" or die "Cannot open rels.txt";
			while(<$fin>) {
				next if /^#/ || /^\s+$/;
				my @iids = split;
				my @fids = uniq sort map { $famid{$_} // do { croak "Cannot find famid for $_" } } @iids;
				unless (@fids > 1) {
					print STDERR $fids[0], ": ", join("\t", @iids), "\n";
					croak "Cryptic related samples should be from different families";
				}
				my $newfid = join($concat, sort @fids);
				# update FIDs
				foreach my $fid (@fids) {
					foreach my $iid (grep { $famid{$_} eq $fid } keys %famid) {
						next if defined $remove{$iid};
						$famid{$iid} = $newfid;
						print $fup join("\t", $fid, $iid, $newfid, $iid), "\n";
					}
				}
			}
		}
		close $fup;
		system("plink --bfile $infile --update-ids $outfile.update-ids --make-bed --out $outfile");
		$infile = $outfile;
	}
	else {
		unlink(glob("$outfile.*"));
		$outfile = $infile;
	}

	# Step 3:
	$outfile = "$geno.step3";
	if (-f "$prefix.sex_update.txt") {
		open my $fin, "$prefix.sex_update.txt" or die "Cannot open sex_update.txt";
		open my $fup, ">$outfile.update-sex" or die "Cannot write to update-sex";
		while(<$fin>) {
			next if /^#/ || /^\s+$/;
			my ($iid, $sex) = (split)[0,1];
			my $gender;
			if ($sex eq '1' || $sex =~ /^M/i) {
				$gender = 1;
			}
			elsif ($sex eq '2' || $sex =~ /^F/i) {
				$gender = 2;
			}
			else {
				warn "Cannot recognize sex: $sex -- set to unknown";
				$gender = 0;
			}
			die "Cannot find family ID for $iid" unless defined $famid{$iid};
			if (defined $swap{$iid}) {
				warn "Sample $iid has already been swapped";
				next;
			}
			if (defined $remove{$iid}) {
				warn "Sample $iid has been removed";
				next;
			}
			print $fup join("\t", $famid{$iid}, $iid, $gender), "\n";
		}
		close $fup;
		system("plink --bfile $infile --update-sex $outfile.update-sex --make-bed --out $outfile");
		$infile = $outfile;
	}
	else {
		unlink(glob("$outfile.*"));
		$outfile = $infile;
	}
	
	# Step 4:
	# If new parents from other families were found, we need to update FIDs
	# manually later.
	$outfile = "$geno.step4";
	if (-f "$prefix.nonpar.txt") {
		my %dadid = map { (split)[1,2] } slurp "$rootdir/wrk/pass_4.fam";
		my %momid = map { (split)[1,3] } slurp "$rootdir/wrk/pass_4.fam"; 
		my (%pars);
		open my $fin, "$prefix.nonpar.txt" or die "Cannot open nonpar.txt";
		while(<$fin>) {
			next if /^#/ || /^\s+$/;
			my ($iid, $role, $parid) = (split)[0,1,2];
			$parid = '0' unless defined $parid;
			die "Cannot find family ID for $iid" unless defined $famid{$iid};
			if ($parid ne '0') {
				die "Cannot find family ID for $parid" unless defined $famid{$parid};
			}
			if (defined $swap{$iid}) {
				warn "Sample $iid has already been swapped, skip parants update";
			}
			if (uc($role) eq 'DAD' || uc($role) eq 'FATHER') {
				$pars{$iid}{DAD} = $parid;
			} 
			elsif (uc($role) eq 'MOM' || uc($role) eq 'MOTHER') {
				$pars{$iid}{MOM} = $parid;
			}
			else {
				die "Cannot recognize the role: $role";
			}
		}
		open my $fup, ">$outfile.update-parents" or die "Cannot write to update-parents";
		foreach my $iid (keys %pars) {
			print $fup $famid{$iid}, "\t", $iid, "\t", $pars{$iid}{DAD} // $dadid{$iid}, "\t", $pars{$iid}{MOM} // $momid{$iid}, "\n";  
		}
		close $fup;
		system("plink --bfile $infile --update-parents $outfile.update-parents --make-bed --out $outfile");
		$infile = $outfile;
	}
	else {
		unlink(glob("$outfile.*"));
		$outfile = $infile;
	}

	if ($infile eq "$rootdir/wrk/pass_4") {
		print STDERR "No patch file was applied\n";
		return $infile;
	}
	else {
		return $outfile;
	}
}

# Note: in population-based sample, perform case-control association test
# In family-based sample, perform TDT or DFAM test
# The presence of inflation in population-based sample may also be due
# to the population stratification, because we have not performed PCA to
# exclude population outliers. Family-based association are more robust.

##########################
## Workflow Components  ##
##########################

# First pass: calculate missing rate and perform sex check
# Obviously bad samples will be removed based on call rate, 
# Sex is imputed based on chrX het rate
# Note: Final data should have no mssing sex
sub first_pass_run {
	my ($chrYflag) = @_;
	my $script = <<'EOF';

plink --bfile _RAWGENO_ --missing --out _WRKDIR_/pass_1

plink --bfile _RAWGENO_ --het --out _WRKDIR_/pass_1

plink --bfile _RAWGENO_ --check-sex --out _WRKDIR_/pass_1_chrX

EOF
	if ($chrYflag) {
		$script .=<<'EOF';

plink --bfile _RAWGENO_ --check-sex ycount --out _WRKDIR_/pass_1_chrXY

EOF
	}
	return $script;
}

sub first_pass_files {
	my ($chrYflag) = @_;
	my @files = map { "wrk/pass_1$_" } qw|_chrX.sexcheck .imiss .lmiss|;
	if ($chrYflag) {
		push @files, "wrk/pass_1_chrXY.sexcheck";
	}	
	return @files;
}

# Plot summary statistics and generate a bad sample list
# Remove samples with grossly large number of missing genotypes
# Impute genders for gender errors or missing genders
# Note that SNP sex is determined from chrX heterozygosity
sub first_pass_plot {
	my $script = <<'EOF';

Rscript _MODULE_/plot_qc_pass1.R _WRKDIR_/pass_1 _PARDIR_/pass_1.txt

# If removal and sex-update list exist, fix the data before going to second pass
if [[ -f _WRKDIR_/pass_1.remove ]]; then
	REMOVE="--remove _WRKDIR_/pass_1.remove"
fi
if [[ -f _WRKDIR_/pass_1.sexupdate ]]; then
	UPDATE="--update-sex _WRKDIR_/pass_1.sexupdate"
fi

plink --bfile _RAWGENO_ $REMOVE $UPDATE --make-bed --out _WRKDIR_/pass_1

# Fix fam file: zero out parents in the fam file with incompatible sex
# and correct sex for parents whose SNP sex are unknown
perl _MODULE_/fix_fam_file.pl _WRKDIR_/pass_1.fam

if [[ -f _WRKDIR_/pass_1.fam.rename ]]; then
	plink --bfile _WRKDIR_/pass_1 --update-ids _WRKDIR_/pass_1.fam.rename --make-bed --out _WRKDIR_/pass_1a
fi


EOF
	return $script;
}

sub first_pass_figs {
	my ($chrYflag) = @_;
	my @files = map { "wrk/pass_1.$_" } qw|bim bed fam fam.bak|;
	push @files => map { "wrk/pass_1_$_.png" } qw|miss het chrX_FvsMiss|;
	if ($chrYflag) {
		push @files, "wrk//pass_1_chrXY_FXvsYcnt.png";
	}
	return @files;
}

# Second pass: perform a quick SNP cleanup, then estimate relatedness
# Then find incorrect familial relpairs and cryptic related pairs in the population
# For incorrect familial pairs, both pair will be excluded
# For cryptic related pairs, only one with the highest call rate will be kept.
# Incorrect familial samples should be fixed manually.
# Currently we only verify 1d relationship.
sub second_pass_run {
	my $script = <<'EOF';

source _PARDIR_/pass_2.txt

if [[ -f _WRKDIR_/pass_1a.bed ]]; then
	PASS1=_WRKDIR_/pass_1a
else
	PASS1=_WRKDIR_/pass_1
fi

plink --bfile $PASS1 --maf $MINMAF --geno $MAXGMIS --hwe $MINHWE \
	--make-bed --out _WRKDIR_/pass_2

plink --bfile _WRKDIR_/pass_2 --missing --out _WRKDIR_/pass_2

#king -b _WRKDIR_/pass_2.bed --related --degree 2 --prefix _WRKDIR_/pass_2
king -b _WRKDIR_/pass_2.bed --related --prefix _WRKDIR_/pass_2

perl _MODULE_/find_err_relpairs.pl _WRKDIR_/pass_2 _PARDIR_/pass_2.txt

EOF
	return $script;
}

sub second_pass_files {
	my ($famflag) = @_;
	my @files = map { "wrk/pass_2.$_" } qw|bim bed fam imiss lmiss|;
	if ($famflag) {
		push @files, map { "wrk/pass_2.$_" } qw|kin relpairs.R relpairs.png relpairs.2.png relpairs.csv|;
	}
	push @files, "wrk/pass_2allsegs.txt";
	return @files;
}

# After error relationship were fixed, we perform detailed SNP level QC
# on the cleaned data after excluding samples with incorrect relationships. 
sub third_pass_run {
	my $script = <<'EOF';

if [[ -f _WRKDIR_/pass_2.relpairs.remove ]]; then
	REMOVE="--remove _WRKDIR_/pass_2.relpairs.remove"
fi
if [[ -f _WRKDIR_/pass_2.relpairs.dups ]]; then
	perl _UTILS_/plink_dup_discord.pl --input _WRKDIR_/pass_1 --dups _WRKDIR_/pass_2.relpairs.dups
	DUP="--duperr _WRKDIR_/pass_1.duperr"
fi

rm -fR _WRKDIR_/pass_3
perl _UTILS_/create_snp_qc_stats.pl --input _WRKDIR_/pass_1 $REMOVE $DUP \
	--output _WRKDIR_/pass_3.snp_qc_stats.txt --wrkdir _WRKDIR_/pass_3

EOF
	return $script;
} 

# Then select appropriate filtering parameters
sub third_pass_filter {
	my $script = <<'EOF';

# Plot SNP statistics and create bad SNP list
Rscript _MODULE_/plot_qc_pass3.R _WRKDIR_ _PARDIR_/pass_3.txt

if [[ -f _WRKDIR_/pass_1a.bed ]]; then
	PASS1=_WRKDIR_/pass_1a
else
	PASS1=_WRKDIR_/pass_1
fi

plink --bfile $PASS1 --exclude _WRKDIR_/pass_3.exclude --make-bed --out _WRKDIR_/pass_4


EOF
	return $script;
}

sub third_pass_files {
	my @files = map { "wrk/pass_3.$_" } qw|badsnps.txt exclude|;
	push @files, map { "wrk/pass_3_$_" } qw|comm_miss.png rare_miss.png hwe_qq.png hist_hhcount.png hist_duperr.png|;
	push @files, map { "wrk/pass_4.$_" } qw|bed bim fam|;
	if ($famflag) {
		push @files, "wrk/pass_3_hist_mendel.png";
	}
	return @files;
}

# Last pass: generate gender check and relpair test statistics
sub fourth_pass_run {
	my ($chrYflag) = @_;
	my $script =<<'EOF';

# After filter, perform gender check and relatedness inference

plink --bfile _WRKDIR_/pass_4 --missing --out _WRKDIR_/pass_4

plink --bfile _WRKDIR_/pass_4 --het --out _WRKDIR_/pass_4

plink --bfile _WRKDIR_/pass_4 --check-sex --out _WRKDIR_/pass_4_chrX

king -b _WRKDIR_/pass_4.bed --related --degree 2 --prefix _WRKDIR_/pass_4

perl _MODULE_/find_err_relpairs.pl _WRKDIR_/pass_4 _PARDIR_/pass_4.txt _WRKDIR_/pass_4

EOF
	if ($chrYflag) {
		$script .=<<'EOF';
plink --bfile _WRKDIR_/pass_4 --check-sex ycount --out _WRKDIR_/pass_4_chrXY

EOF
	}
	$script .=<<'EOF';

if [[ -f _WRKDIR_/pass_1a.bed ]]; then
	PASS1=_WRKDIR_/pass_1a
else
	PASS1=_WRKDIR_/pass_1
fi

grep -v '^IGNORE' _PARDIR_/pass_4.txt > _TMPDIR_/pass_4.txt
Rscript _MODULE_/plot_qc_pass4.R _WRKDIR_/pass_4 $PASS1.fam _TMPDIR_/pass_4.txt _WRKDIR_/pass_4

EOF

	return $script;
}

sub fourth_pass_files {
	my ($chrYflag, $famflag) = @_;
	my @files = map { "wrk/pass_4$_" } qw|.lmiss .imiss .het _chrX.sexcheck _HetvsMiss.png
										_sample.csv _chrX_FXvsMiss.png _chrX_F_XvsAUTO.png|;
	if ($chrYflag) {
		push @files, "wrk/pass_4_chrXY.sexcheck", "wrk/pass_4_chrXY_FXvsYcnt.png";
	}
	if ($famflag) {
		push @files, "wrk/pass_4.kin";
		push @files, map { "wrk/pass_4.relpairs.$_" } qw|R csv png 2.png|; 
	}
	push @files, "wrk/pass_4allsegs.txt";
	return @files;
}

# Final step: QC check on the cleaned data
sub fifth_pass_run {
	my ($chrYflag, $ccflag, $famflag) = @_;
	my $script =<<'EOF';

# After filter, perform gender check and relatedness inference

plink --bfile _OUTDIR_/_PREFIX_ --missing --out _OUTDIR_/_PREFIX_

plink --bfile _OUTDIR_/_PREFIX_ --het --out _OUTDIR_/_PREFIX_

plink --bfile _OUTDIR_/_PREFIX_ --check-sex --out _OUTDIR_/_PREFIX_'_chrX'

king -b _OUTDIR_/_PREFIX_.bed --related --degree 2 --prefix _OUTDIR_/_PREFIX_

perl _MODULE_/find_err_relpairs.pl _OUTDIR_/_PREFIX_ _PARDIR_/pass_5.txt _OUTDIR_/_PREFIX_

EOF
	if ($chrYflag) {
		$script .=<<'EOF';
plink --bfile _OUTDIR_/_PREFIX_ --check-sex ycount --out _OUTDIR_/_PREFIX_'_chrXY'

EOF
	}
	$script .=<<'EOF';

grep -v '^IGNORE' _PARDIR_/pass_5.txt > _TMPDIR_/pass_5.txt
Rscript _MODULE_/plot_qc_pass4.R _OUTDIR_/_PREFIX_ _OUTDIR_/_PREFIX_.fam _TMPDIR_/pass_5.txt _OUTDIR_/_PREFIX_

EOF
	if ($ccflag) {
		$script .= 'source _PARDIR_/pass_5.txt'."\n";
		$script .= 'if [[ $GWAS =~ ^T ]]; then'."\n";
		if ($famflag) {
			$script .=  "  plink-1.07 --noweb --bfile _OUTDIR_/_PREFIX_ --dfam --out _OUTDIR_/_PREFIX_\n".
						"  Rscript _UTILS_/plink_qqman.R _OUTDIR_/_PREFIX_.dfam _OUTDIR_/_PREFIX_.bim\n";
		}
		else {
			$script .=  "  plink --bfile _OUTDIR_/_PREFIX_ --assoc --out _OUTDIR_/_PREFIX_\n".
						"  Rscript _UTILS_/plink_qqman.R _OUTDIR_/_PREFIX_.assoc _OUTDIR_/_PREFIX_.bim\n";
		}
		$script .= "fi\n\n";
	}
	return $script;	
}

sub fifth_pass_files {
	my ($chrYflag, $famflag) = @_;
	my @files = map { "out/$outprefix$_" } qw|.lmiss .imiss .het _chrX.sexcheck _HetvsMiss.png
										_sample.csv _chrX_FXvsMiss.png _chrX_F_XvsAUTO.png|;
	if ($chrYflag) {
		push @files, "out/${outprefix}_chrXY.sexcheck", "out/${outprefix}_chrXY_FXvsYcnt.png";
	}
	if ($famflag) {
		push @files, "out/$outprefix.kin";
		push @files, map { "out/$outprefix.relpairs.$_" } qw|R csv png 2.png|; 
	}
	else {
		push @files, "out/$outprefix.kin0";
	}
	return @files;
}

######################
## Helper functions ##
######################

sub write_pars {
	my ($param, $file) = @_;
	my $fout = IO::File->new($wkf->get_subdir("par")."/$file", "w") 
		or croak "Cannot write $file";
	while(my ($key, $val) = each %$param) {
		print $fout $key."=".$val,"\n";
	}
}

sub read_pars {
	my ($params, $prefix) = @_;
	croak "The length of input params list must be a multiple of 3" if scalar(@$params)%3 > 0;

	my $parfile = $wkf->get_subdir("par")."/$prefix.txt";
	if (-f $parfile) {
		my %pars;
		my $fin = IO::File->new($parfile);
		while(<$fin>) {
			chomp();
			my ($key, $val) = split('=', $_);
			$pars{$key} = $val;
		}
		for(my $ii = 0; $ii < int(@$params/3); $ii ++) {
			croak "Parameter $params->[3*$ii] cannot be found" unless defined $pars{$params->[3*$ii]};
			$params->[3*$ii+2] = $pars{$params->[3*$ii]};
		}
	}

	return $params;
}

sub interactive_run {
	my ($params, $prefix, $taskname) = @_;
	$params = read_pars($params, $prefix);

	while(my $par = prompt_params($params)) {
		write_pars($par, "$prefix.txt");
		$wkf->run({ tasks => $taskname, dryrun => $opt->get_dryrun, interact => 1 });
		last if $opt->get_dryrun;
		print "Please check the expected output: \n";
		print join("\n", map { "\t".$_ } $wkf->get_expected($taskname)), "\n";
		my $rerun = prompt("Re-run $taskname? ", -yn, -default => 'n');
		if ($rerun) {
			foreach my $outfile ($wkf->get_expected($taskname)) {
				unlink $outfile if -f $outfile;
				#last;
			}
			if (defined $optional{$taskname}) {
				foreach my $file (@{$optional{$taskname}}) {
					my $path = $rootdir."/".$file;
					unlink $path if -f $path;
				}
			}
		}
		else {
			last;
		}
	}
}

sub count_samp {
	my ($famfile, $rmfile) = @_;
	my (%fiids);
	open my $fin, $famfile or die "Cannot open $famfile";
	while(<$fin>) {
		my ($fid, $iid, $sex) = (split)[0,1,4];
		$fiids{$fid,$iid} = $sex;
	}
	if (-f $rmfile) {
		open my $frm, $rmfile or die "Cannot open $rmfile";
		while(<$fin>) {
			my ($fid, $iid) = (split)[0,1];
			delete $fiids{$fid,$iid};
		}
	}
	my $n_m = grep { $_ == 1 } values %fiids;
	my $n_f = grep { $_ == 2 } values %fiids;
	my $n_tot = keys %fiids;
	my ($n_trio, $n_duo) = (0, 0);
	open $fin, $famfile or die "Cannot open $famfile";
	while(<$fin>) {
		my ($fid, $iid, $dad, $mom) = (split)[0,1,2,3];
		if ($dad ne '0' && $mom ne '0') {
			if (all { defined $fiids{$fid,$_} } ($iid, $dad, $mom)) {
				$n_trio ++;
			}
			elsif (defined $fiids{$fid,$iid} && defined $fiids{$fid,$dad}) {
				$n_duo ++;
			}
			elsif (defined $fiids{$fid,$iid} && defined $fiids{$fid,$mom}) {
				$n_duo ++;
			}
		}
		elsif ($dad ne '0') {
			if (defined $fiids{$fid,$iid} && defined $fiids{$fid,$dad}) {
				$n_duo ++;
			}
		}
		elsif ($mom ne '0') {
			if (defined $fiids{$fid,$iid} && defined $fiids{$fid,$mom}) {
				$n_duo ++;
			}
		}
	}
	return { Male => $n_m, Female => $n_f, Total => $n_tot, Trio => $n_trio, Duo => $n_duo };
}

sub count_dups {
	my ($dupfile) = @_;
	my $npairs;
	open my $fin, $dupfile or die "Cannot open dup file $dupfile";
	while(<$fin>) {
		my @samps = split;
		my $nsamp = @samps / 2;
		unless ($nsamp >= 2) {
			croak "Incorrect number of samples in the line: $_";
		}
		$npairs += $nsamp*($nsamp-1)/2;
	}
	return $npairs;
}


