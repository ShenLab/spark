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
use IO::Prompt;
use Perl6::Slurp;
use List::MoreUtils qw|all|;
use String::ShellQuote;
use Utils::Workflow;
use Utils::Hash qw|merge_conf|;
use Genet::File::VCF qw|parse_vcf_header|;

use lib "$Bin/../lib/";
use Shared qw|read_list|;
use Wrapper qw|bam_sampids|;

############################
## Command line interface ##
############################

my @spec =  (
	Param("conf|c")->valid(sub { -r }),
	Param("vcf|v")->valid(sub { -r $_ || -r "$_.bed" && -r "$_.fam" && -r "$_.bim" }),
	Param("list|l")->valid(sub { -r $_  }),
	Param("rename")->valid(sub { -r $_ }),
	Param("outdir|out"),
	Param("engine|eng")->valid(sub { $_ eq 'SGE' || $_ eq 'BASH' }),
	Keypair("param|par"),
	Switch("search"),
	Switch("help|h"),
	Switch("dryrun|dry"),
	Switch("force")
	);


my $opt = Getopt::Lucid->getopt(\@spec);

if ($opt->get_help) {
	print STDERR <<EOF;
Purpose:
	This is pipeline script to verify sample ID of BAM files based on known SNP genotypes.

Usage:
	bam_id.pl --conf Configs --vcf Genotypes --list BAMList --outdir RootDir

Options:
	--list: BAM file list. It has two columns per line: absolute path to the BAM file and sample ID. 
			In case only one column is provided, it will be taken from the basename of file after 
			stripping of suffix. ID should be matched to SNP genotypes in the VCF file. 
	--rename: A rename list for BAM files to match the IDs in genotypes.
	--vcf:  Multi-sample VCF file of SNP genotypes on selected sites. AN and AC will be used to 
			determine population allele frequency.
	--search: The default test is to verify the sample CRAM/BAM file matched to the same sample in VCF.
			  Under search mode, each CRAM/BAM file will also be tested against all samples in VCF
			  to find the best matching sample.

Notes:
	The matching of BAM/CRAM to known SNP genotypes is performed by verifyBamID (v1). Only samples in
	VCF will be used for testing. "Sequence only" mode test of contamination is part of bam_qc pipeline.
	CRAM files are supported but it involves convertion of CRAM to BAM at known SNP sites.

	When array-based genotypes are used, VCF file can be created from plink bped format. But reference
	alleles must be corrected to match the reference sequence. And AN/AC fields should be added to the VCF. 
	
Output:
	VerifyBamID will be run on each sample BAM file against VCF in parallel. Intermediate results per 
	individual are stored in wrk directory. Final results for matching the same individual genotype 
	for each BAM/CRAM will be written to summary.txt under out directory. 

	Under search mode, the best matching individual will be written to wrk/IID.bestSM for each individual.  
	It can be used for diagnosis when sample mismatches are found.

EOF
	exit 1;
}

$opt->validate({requires => [qw|conf vcf list outdir|]});

# Read in config parameters
my %conf = merge_conf($opt->get_conf, $opt->get_param); 
# Path to find module scripts
$conf{PATH}{MODULE} = shell_quote("$Bin/module");


############################
## Input files validation ##
############################
my %opt = exists $conf{AWS}{STAGEIN} ? (notest => 1) : ();
my ($bams, $filetype) = read_list($opt->get_list, { suffix => ['bam', 'cram'],
						rename => $opt->get_rename});
my %bams = %$bams;

my $rootdir = $opt->get_outdir;

# We will only consider IIDs within genotype/VCF file
my (%vcfids, $genoflag);
{
	my $geno = $opt->get_vcf;
	if ($geno =~ /\.vcf\.gz$/) {
		my ($header, $info, $sampids) = parse_vcf_header($geno);
		@vcfids{@$sampids} = @$sampids;
		if (defined $conf{PATH}{OVERCHAIN}) {
			$conf{PATH}{GENO} = abs_path($geno);
			$conf{PATH}{VCF} = "$rootdir/tmp/cohort.vcf.gz";
		}
		else {
			$conf{PATH}{VCF} = abs_path($geno);
		}
	}
	elsif (all { -f "$geno.$_" } qw|bed bim fam|) {
		$conf{PATH}{GENO} = abs_path($geno);
		$conf{PATH}{VCF} = "$rootdir/tmp/cohort.vcf.gz";
 		%vcfids = map { (split)[1] => 1 } slurp "$geno.fam";
		$genoflag = 1;
	}
	else {
		die "Cannot determine the type of genotype";
	}
}

my $nvcf = keys %vcfids;
my $nbams = keys %bams;
my $ncomm = grep { defined $vcfids{$_} } keys %bams;
print STDERR "A total of $nvcf samples with known genotypes, $nbams in $filetype file list;  $ncomm are in common and will be used in ID test\n";

exit 1 unless $ncomm;

#unless (prompt("Continue to run the workflow: ", -yes_no, -default => 'y')) {
#	exit 1;
#}


##############################
## Workflow initialization  ##
## Config file parsing      ##
##############################

my $wkf = Utils::Workflow->new($rootdir, 
	{ engine => $opt->get_engine, force => $opt->get_force});

my %deparg;
if ($genoflag) {
	$wkf->add(prep_vcf(), { name => "PrepVCF",  
		expect => ["tmp/cohort.vcf.gz", "tmp/cohort.vcf.gz.tbi"]});
	$deparg{depend} = "PrepVCF";
}
else {
	if (defined $conf{PATH}{OVERCHAIN}) {
		$wkf->add(lift_over(), name => "LiftOver",
			expect => ["tmp/cohort.vcf.gz", "tmp/cohort.vcf.gz.tbi"]);
		$deparg{depend} = "LiftOver";
	}
}

# Add jobs to workflow
$wkf->add(test_bamid($filetype), { name => "TestBamID", nslots => $ncomm, %deparg,
	expect => [ map { ["wrk/$_.selfSM", "wrk/$_.depthSM", "wrk/$_.log", "wrk/$_.vcf"] } 
					   grep { defined $vcfids{$_} } sort keys %bams ] })
	->add(sum_results(), { name => "Summary", depend => "TestBamID", 
		expect => "out/summary.csv" });

# Check external file
croak "Cannot find FASTA reference!" unless -f "$conf{PATH}{FASTA}";

# Write run-time config to par directory
write_config %conf, "$rootdir/par/run.conf" unless $opt->get_dryrun;


##############################
## Working directory setup  ##
##############################

# Write IID2BAM parameter files to par dir.
my $ii = 0;
foreach my $iid (sort keys %bams) {
	if (defined $vcfids{$iid}) {
		$ii ++;
		open my $fout, ">$rootdir/par/IID2BAM.$ii" or die "Cannot write to par/IID2BAM.$ii";
		print $fout $iid, "\t", $bams{$iid}, "\n";
	}
}


################################
## Kickstart workflow engine  ##
################################

$wkf->inst(\%conf);
$wkf->run({ conf => $conf{$wkf->{engine}}, dryrun => $opt->get_dryrun });


############################
## Workflow components    ##
############################

# Note: VCF converted from plink may not have correct REF and ALT alleles!
# After liftover, variants whose REF allele does match the reference seq will be dropped.
sub prep_vcf {
	my $script = <<'EOF';

TMPVCF=$(echo _PATH.VCF_ | sed -e 's/\.vcf\.gz$//')

plink --chr 1-22 --bfile _PATH.GENO_ --recode vcf-iid --out $TMPVCF

EOF
	if (chr_flag($conf{PATH}{FASTA})) {
		$script =~ s/vcf-iid/vcf-iid --output-chr chrM /;
	}
	if (defined $conf{PATH}{OVERCHAIN}) {
		$script .= <<'EOF';

liftover_vcf.pl --invcf $TMPVCF.vcf --chain _PATH.OVERCHAIN_ --seq _PATH.FASTA_ \
	--hgchr --output $TMPVCF.lifted --wrkdir _TMPDIR_/liftover

zcat $TMPVCF.lifted.vcf.gz | fill-an-ac | bgzip -c > _PATH.VCF_
tabix -p vcf _PATH.VCF_

EOF

	}
	else {
		$script .= <<'EOF';

fill-an-ac < $TMPVCF.vcf | bgzip -c > _PATH.VCF_
tabix -p vcf _PATH.VCF_

EOF
	}
	return $script;
}

sub chr_flag {
	my ($fasta) = @_;
	open my $fin, $fasta or die "Cannot open FASTA: $fasta";
	my $header = <$fin>;
	if ($header =~ /^>chr/) {
		return 1;
	}
	else {
		return 0;
	}
}

sub lift_over {
	my $script = <<'EOF';

liftover_vcf.pl --invcf _PATH.GENO_ --chain _PATH.OVERCHAIN_ --seq _PATH.FASTA_ \
	--hgchr --output $TMPVCF.lifted --wrkdir _TMPDIR_/liftover

zcat $TMPVCF.lifted.vcf.gz | fill-an-ac | bgzip -c > _PATH.VCF_
tabix -p vcf _PATH.VCF_

EOF
	return $script;
}

sub test_bamid {
	my ($filetype) = @_;
	my $script = "\nread IID BAMFILE < _PARDIR_/IID2BAM._INDEX_\n";
	if (exists $conf{AWS}{STAGEIN}) {
		$script .= <<'EOF';
TMPDIR=$(stagein.pl --input $BAMFILE --output _TMPDIR_/BAM._INDEX_ --profile _AWS.PROFILEIN_ --tmpdir _AWS.STAGEIN_ --all)

trap "echo Clean up $TMPDIR; rm -fR $TMPDIR; exit 1" EXIT SIGINT SIGTERM

read BAMFILE < _TMPDIR_/BAM._INDEX_

EOF
	}
	$script = <<'EOF';
if (( _PARAM.MAXDEPTH_ > 50  )); then
	PRECISE="--precise"
else 
	PRECISE=""
fi

EOF
	unless($opt->get_search) {
		$script .= <<'EOF';
GENOTYPE=_WRKDIR_/$IID.vcf

bcftools view -s $IID -I -v snps -o _WRKDIR_/$IID.vcf _PATH.VCF_ 

EOF
	}
	else {
		$script .= <<'EOF';
GENOTYPE=_PATH.VCF_

EOF
	}
	if ($filetype eq 'bam') {
		$script .= <<'EOF';
verifyBamID --vcf $GENOTYPE --bam $BAMFILE --out _WRKDIR_/$IID \
	--best --ignoreRG --maxDepth _PARAM.MAXDEPTH_ $PRECISE

EOF
	}
	else {
		$script .= <<'EOF';
# Extract targeted region from CRAM file using GATK
# we use GATK printreads module to extract aligned reads around common snps
# GATK v3 contains a bug, so this step is implemented in GATK4

gatk --java-options "_RSRC.GATKOPT_" PrintReads \
	-I $BAMFILE  -L $GENOTYPE \
	-R _PATH.FASTA_  -O _WRKDIR_/$IID.bam

verifyBamID --vcf $GENOTYPE --bam _WRKDIR_/$IID.bam --out _WRKDIR_/$IID \
	--best --ignoreRG --maxDepth _PARAM.MAXDEPTH_ $PRECISE 

sleep 2

rm -f _WRKDIR_/$IID.bam
rm -f _WRKDIR_/$IID.bai

EOF
	}
	return $script;
}

sub sum_results {
	my $script = <<'EOF';

perl _PATH.MODULE_/collect_bamid_res.pl _WRKDIR_  _OUTDIR_ 

EOF
	unless($opt->get_search) {
		$script .= <<'EOF'

csvtk filter2 -f '$FREEMIX>_PARAM.MINFREEMIX_ || $CHIPMIX>_PARAM.MINCHIPMIX_' _OUTDIR_/summary.csv \
	> _OUTDIR_/badsamps.csv

LC=$(wc -l _OUTDIR_/badsamps.csv | awk '{print $1}')
if [[ $LC == 0 || $LC == 1 ]]; then
	rm -f _OUTDIR_/badsamps.csv
fi

EOF
	}
	return $script;
}

