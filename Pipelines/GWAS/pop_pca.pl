#!/usr/bin/env perl

use strict;
use warnings;
use FindBin qw|$Bin|;
use Config::Std;
use File::Basename;
use String::ShellQuote;
use Hash::Util qw|lock_hash_recurse|;
use Getopt::Lucid qw|:all|;
use List::MoreUtils qw|all none|;
use Utils::Workflow;
use Utils::Hash qw|merge_conf|;


#############################
## Command line interface  ##
#############################
# The default input is plink bped format genotypes and a output directory
my @spec = (
	Param("conf|c")->valid(sub { -r }),
	Param("geno|i")->valid(sub { -r "$_.bed" && -r "$_.bim" && -f "$_.fam" }),
	Param("outdir|out"),
	Param("engine|eng")->valid(sub { $_ eq 'SGE' || $_ eq 'BASH' }),
	Param("prefix"),
	Keypair("param|par"),
	Switch("help|h"),
	Switch("dryrun|dry"),
	Switch("force")
	);

my $opt = Getopt::Lucid->getopt(\@spec);

if ($opt->get_help) {
	print STDERR <<'EOF';
Purpose:
	This is a pipeline script to run PCA on SNP genotypes and infer population ancestry.

Usage:
	pop_pca.pl --conf Config [--geno SNPGeno --prefix Prefix] --out OutDir

Options:
	--geno: SNP genotype file prefix in plink binary PED format. 
			It overrides the genotype file prefix in config.
	--prefix: Output file prefix. Default will be the prefix of input gneotype file.

Notes:
	The input to PCA analysis include SNP genotype data from one or more reference populations and 
	a test population. All genotype should be stored in binary PED format, and assuming sample and 
	SNPs passed QC.

	The script will run following steps:
	1. Extract genotypes of (presumed unrelated) founders from reference populations.
	2. Merge genotypes from reference panel with the testing cohort.
	   All genotypes will be lifted over to the same genome assembly if necessary before merging.
	3. Infer PC axes using reference panel and project test cohort onto those PCs.
	   If samples in the test cohort are unrelated, direct PCA can also be performed on merged 
	   reference+test populations. 
	4. Calculate ancestral probabilities for each individual in the test cohort.
	   We implemented the geometric algorithm of SNPweight (Chen 2013) for ancestry inference.
	   The populations in the reference should represent different ancestral populations
	   and should not include admixed populations.

	   Chen et al. 2013: Improved ancestry inference using weights from external reference panels
 	   Bioinformatics 29:1399–1406.

Output:
	Prefix.PCs.txt - First N PCs. The default N is 10.
	Prefix_PopRefs.PCs1-$Num.png - Pairwise scatterplots of first $Num (<=N) of PCs of reference populations
	Prefix.PCs1-${Num}_Prob${Prob}.png - Pairwise scatterplots of first $Num (<=N) of PCs of the test population.
				Samples are color coded by predicted ancestral population using probability threshold $Prob.  
	Prefix.Ancestry.txt - (Available if ancestry prediction is switched on) Predicted population ancestry.
				For each test sample, it gives probabilities of different ancestral populations represented 
				in the population references.

Dependencies:
	smartpca, plink
	There is a technical issue in smartpca that sample ID can only be no more than 39 characters.
	In merged genotpes, cohort label + ":" + original ID in genotype data is used as sample ID.
	Make sure IDs used by smartpca do not exceed smartpca's max length. 

EOF
	exit 1;
}

$opt->validate({requires => [qw|conf outdir|]});

my $rootdir = $opt->get_outdir;

my %conf = merge_conf($opt->get_conf, $opt->get_param); 
$conf{default}{MODULE} = shell_quote("$Bin/module");


if ($opt->get_geno) {
	# overriding the test genotype in the config file
	$conf{TEST}{GENO} = $opt->get_geno;
}
unless(all { -f $conf{TEST}{GENO}.".$_" } qw|bim bed fam|) {
	die "Not all test genotype files can be found: $conf{TEST}{GENO}";
}
# Out prefix will be used as label for output, test file, and all test samples
if ($opt->get_prefix) {
	$conf{TEST}{LABEL} = $opt->get_prefix;
}
else {
	unless(defined $conf{TEST}{LABEL}) {
		$conf{TEST}{LABEL} = (split(q|\.|,basename($conf{TEST}{GENO})))[0];
	}
}

lock_hash_recurse(%conf);

my ($nofamilyname, $nopcaproj, $noancpred);

sub check_sampids {
	my ($prefix, $label) = @_;
	my $lablen = length($label);
	open my $fin, "$prefix.fam" or die "Cannot read $prefix.fam";
	while(<$fin>) {
		my ($fid, $iid) = (split)[0,1];
		my $fidlen = length($fid);
		my $iidlen = length($iid);
		if ($lablen + $iidlen + 1 > 39) {
			if ($iidlen + 1 + 1 <= 39) {
				print STDERR "$label:$fid:$label:$iid\n";
				die "SmartPCA only support sample ID up to max. 39 character, try to change another label";
			}
			else {
				die "SmartPCA only support sample ID up to max. 39 character."
			}
		}
		elsif ($lablen*2 + $fidlen + $iidlen + 3 > 39) {
			$nofamilyname = 1;
		}
	}
}

check_sampids($conf{TEST}{GENO}, $conf{TEST}{LABEL});


#################################
##  Workflow initialization    ##
##  Working directory setup    ##
#################################

my $wkf = Utils::Workflow->new($rootdir, 
	{ engine => $opt->get_engine, force => $opt->get_force});

my @refs = grep {  /^REF_/ } sort keys %conf;
my @prefs;
if (@refs > 0) {
	# Writing configs for reference panel data
	my %known = ($conf{TEST}{LABEL} => 1);
	open my $fout, ">$rootdir/par/POPREF" or die "Cannot write to POPREF";
	open my $flst, ">$rootdir/par/PROCESS" or die "Cannot write to PROCESS";
	foreach my $ref (@refs) {
		unless(all { -f $conf{$ref}{GENO}.".$_" } qw|bim bed fam|) {
			die "Not all genotype related files can be found for $ref";
		}
		unless (-f $conf{$ref}{LABEL}) {
			die "Cannot find population label file for $ref"
		}
		my $prefix = (split('_', $ref))[1];
		check_sampids($conf{$ref}{GENO}, $prefix);
		if (defined $known{$prefix}) {
			die "Population label $prefix has been defined!";
		}
		print $fout join("\t", @{$conf{$ref}}{qw|GENO LABEL POP|}, $prefix), "\n";
		push @prefs, $prefix;
		if (exists $conf{$ref}{CHAIN}) {
			unless (-f $conf{$ref}{CHAIN}) {
				die "Liftover chain file $conf{$ref}{CHAIN} does not exist";
			}
			print $flst $prefix, "\t", $conf{$ref}{CHAIN}, "\n"; 
		}
		else {
			print $flst $prefix, "\tNA\n";
		}
	}
	print STDERR "The following population reference data will be used:\n";
	print STDERR join(" ", @prefs), "\n";
	
	# Also adding test into process list
	if (exists $conf{TEST}{CHAIN}) {
		unless (-f $conf{TEST}{CHAIN}) {
			die "Liftover chain file $conf{TEST}{CHAIN} does not exist";
		}
		print $flst $conf{TEST}{LABEL}, "\t", $conf{TEST}{CHAIN}, "\n";
	}
	else {
		print $flst $conf{TEST}{LABEL}, "\tNA\n";
	}
}
else {
	print STDERR "No population reference data was found in the config";
	exit 1;
}

# Prepare smartpca parameter file
{
	open my $fout, ">$rootdir/par/SMARTPCA" or die "Cannot write to SMARTPCA";
	print $fout <<EOF;
genotypename:       $rootdir/out/$conf{TEST}{LABEL}_Merged.bed
snpname:            $rootdir/out/$conf{TEST}{LABEL}_Merged.bim
indivname:          $rootdir/out/$conf{TEST}{LABEL}_Merged.fam
evecoutname:        $rootdir/out/$conf{TEST}{LABEL}.evec
evaloutname:        $rootdir/out/$conf{TEST}{LABEL}.eval
EOF
	# Default is to perform PCA projection, 
	# If PCAPROJ is specified and set to FALSE or NO, we will do PCA on merged ref + test data
	# make sure all samples in the test data are unrelated
	unless (exists $conf{TEST}{PCAPROJ} && $conf{TEST}{PCAPROJ} =~ /^F|N/i) {
		print $fout "poplistname:	$rootdir/out/$conf{TEST}{LABEL}_PopRefs.list\n";
	}
	else {
		$nopcaproj = 1;
	}
	if (exists $conf{TEST}{INFERANC} && $conf{TEST}{INFERANC} =~ /^F|N/i) {
		$noancpred = 1;
	}

	if ($nofamilyname) {
		print $fout "familynames:	NO\n";
	}
	# Additional parameters
	foreach my $parname (keys %{$conf{SMARTPCA}}) {
		print $fout $parname, ":\t", $conf{SMARTPCA}{$parname}, "\n";
	}
}

###########################################
## Workflow initialization and kickstart ##
###########################################

$wkf->add(extract_popref(), { name => "ExtractPopRef", expect => [popref_output()] })
 	->add(extract_test(),   { name => "ExtractTest", expect => [testdat_output()] })
 	->add(proc_geno(), { name => "ProcGeno", expect => [geno_output()],
 						 depend => [qw|ExtractPopRef ExtractTest|] } )
 	->add(merge_reftest(), { name => "MergeGeno", expect => [map { "out/$conf{TEST}{LABEL}_Merged.$_" } qw|bim bed fam fam.bak|],
 							 depend => "ProcGeno"})
 	->add(run_smartpca(), { name => "SmartPCA", expect => [map { "out/$conf{TEST}{LABEL}.$_"} qw|evec eval| ],
 							depend => "MergeGeno" })
 	->add(proc_pca(), { name => "ProcPCA", depend => "SmartPCA",
 						expect => $noancpred ? 
 						[map { "out/$conf{TEST}{LABEL}$_" } (".PCs.txt", ".PCs1-$conf{PLOT}{NUMPC}_Prob$conf{PLOT}{PROB}.png", "_PopRefs.PCs1-$conf{PLOT}{NUMPC}.png")] :
 						[map { "out/$conf{TEST}{LABEL}$_" } (".PCs.txt",  ".Ancestry.txt", ".Ancestry.coeffs.txt", ".PCs1-$conf{PLOT}{NUMPC}_Prob$conf{PLOT}{PROB}.png", 
 																"_PopRefs.PCs1-$conf{PLOT}{NUMPC}.png")]
 						 });

write_config %conf, "$rootdir/par/run.conf" unless $opt->get_dryrun;

$wkf->inst(\%conf);
$wkf->run({ conf => $wkf->{engine} eq 'BASH' ? undef : $conf{$wkf->{engine}}, dryrun => $opt->get_dryrun  });


############################
## Workflow components    ##
############################

# Extract subset of founders from population reference panel
sub extract_popref {
	my $script = <<'EOF';

cat _PARDIR_/POPREF | while read GENO LABEL POPSTR PREFIX; do
	perl _MODULE_/extract_founders.pl $GENO $LABEL $POPSTR _WRKDIR_/$PREFIX
done

find _WRKDIR_ -name '*.label' | xargs cut -f2 | sort | uniq > _OUTDIR_/_TEST.LABEL_"_PopRefs.list"

EOF
}

sub popref_output {
	my @exps;
	foreach my $pref (@prefs) {
		push @exps, map { "wrk/$pref.$_" } qw|bim bed fam log keep label|;
	}
	push @exps, "out/$conf{TEST}{LABEL}_PopRefs.list";
	return @exps;
}


sub extract_test {
	my $option = "";
	foreach my $action (qw|exclude remove keep|) {
		my $ACTION = uc($action);
		if (exists $conf{TEST}{$ACTION}) {
			unless (-f $conf{TEST}{$ACTION}) {
				die "Cannot find $action file: $conf{TEST}{$ACTION}";
			}
			if ($action eq 'exclude') {
				$option .= "--$action range $conf{TEST}{$ACTION} ";
			}
			else {
				$option .= "--$action $conf{TEST}{$ACTION} ";
			}
		}
	}
	my $script = <<'EOF';

plink --bfile _TEST.GENO_ OPTION --make-bed --out _WRKDIR_/_TEST.LABEL_

awk -v LABEL=_TEST.LABEL_ '{print $2, LABEL}' _WRKDIR_/_TEST.LABEL_.fam > _WRKDIR_/_TEST.LABEL_.label

EOF
	$script =~ s/OPTION/$option/;
	return $script;
}

sub testdat_output {
	my @exps = map { "wrk/$conf{TEST}{LABEL}.$_" } qw|bim bed fam log label|;
	return @exps;
}

# Process reference panel and test genotypes to the same genome assembly
# and using the same marker name nomenclature. 
sub proc_geno {
	my $script = <<'EOF';

cat _PARDIR_/PROCESS | while read PREFIX CHAIN; do
	rm -fR _TMPDIR_/$PREFIX
	if [[ $CHAIN == "NA" ]]; then
		merge_bped.pl --input _WRKDIR_/$PREFIX --output _WRKDIR_/$PREFIX.Proc --tmpdir _TMPDIR_/$PREFIX
	else 
		merge_bped.pl --input _WRKDIR_/$PREFIX --output _WRKDIR_/$PREFIX.Proc --tmpdir _TMPDIR_/$PREFIX --chain $CHAIN
	fi
done

EOF
}

sub geno_output {
	my @exps;
	foreach my $pref (@prefs, $conf{TEST}{LABEL}) {
	 	push @exps => map { "wrk/$pref.Proc.$_" } qw|bim bed fam|;
	}
	return @exps;
}

# Merge extracted reference panel with test genotypes
sub merge_reftest {
	my $script = <<'EOF';

INPUT=$(cut -f1 _PARDIR_/PROCESS | perl -e '@a=map { chomp; "_WRKDIR_/$_.Proc" } <>; print join(" ", @a)')
TAG=$(cut -f1 _PARDIR_/PROCESS | perl -MFile::Basename -e '@a=map { chomp; $_ } <>; print join(",", @a)')

rm -fR _TMPDIR_/Merged
merge_bped.pl --input $INPUT --processed --output _OUTDIR_/_TEST.LABEL_"_Merged" --tags $TAG --tmpdir _TMPDIR_/_TEST.LABEL_"_Merged"

# Fix the merged fam file to adding population labels
perl _MODULE_/fix_fam_4pca.pl _OUTDIR_/_TEST.LABEL_"_Merged.fam" _WRKDIR_/*.label 

EOF
}


sub run_smartpca {
	my $script = <<'EOF';

smartpca -p _PARDIR_/SMARTPCA > _OUTDIR_/_TEST.LABEL_.log

EOF
}

# Post-process smart PCA output
sub proc_pca {
	my $script = <<'EOF';

# PCs.txt will be used for plot and ancestral inference
perl _MODULE_/refmt_smartpca_out.pl _OUTDIR_/_TEST.LABEL_.evec _OUTDIR_/_TEST.LABEL_.PCs.txt

EOF
	unless ($noancpred) {
		$script .= <<'EOF'
# Ancestry.txt is the final output of predicted ancestry
perl _MODULE_/infer_ancestry.pl _OUTDIR_/_TEST.LABEL_.evec _OUTDIR_/_TEST.LABEL_.Ancestry.txt _TEST.LABEL_

EOF
	}
	$script .= <<'EOF';
# Plot PCs
Rscript _MODULE_/plot_pop_pca.R _OUTDIR_/_TEST.LABEL_ _PLOT.NUMPC_ _PLOT.PROB_

EOF

}



