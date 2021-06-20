#!/usr/bin/env perl

use strict;
use warnings;
use File::Copy;
use Data::Dumper;
use Data::Table;
use Utils::File::Iter qw|data_tab|;

my ($pedfile, $wrkdir, $outdir) = @ARGV;

unless(@ARGV == 3) {
	print STDERR <<EOF

Purpose: Used by samp_qc to merge sample level information into one file

Usage: merge_samp_qcstats PEDFILE WRKDIR OUTDIR

Pedigree structure, sample sex and phenotypes are defined in pedigree file.
Raw output from vcftools, plink, and peddy should be found in working dir.
Output will be csv format, contain all information about sample quality,
gender, etc.

EOF
}
else {
	die "Cannot find ped file" unless -f $pedfile;
	die "Cannot find working directory" unless -d $wrkdir;
	die "Cannot find output directory" unless -d $outdir;
}

# First apply simple patch to imendel file
if (-f "$wrkdir/geno.imendel") {
	copy "$wrkdir/geno.imendel", "$wrkdir/geno.imendel.bak";
	my %totct;
	open my $fin, "$wrkdir/geno.imendel.bak" or die "Cannot open geno.imendel.bak"; <$fin>;
	while(<$fin>) {
		my ($fid, $iid, $count) = split;
		$totct{"$fid\t$iid"} += $count;
	}
	# write to new imendel file
	open my $fout, ">$wrkdir/geno.imendel" or die "Cannot write to geno.imendel";
	print $fout join("\t", qw|FID IID N|), "\n";
	foreach my $fiid (keys %totct) {
		print $fout $fiid, "\t", $totct{$fiid}, "\n";
	}
}

my %opts = (
	'site.auto.idepth' => { alias => { INDV => "IID", MEAN_DEPTH => "DP_auto" },
							subset => [qw|IID DP_auto|] },
	'site.chrX.idepth' => { alias => { INDV => "IID", MEAN_DEPTH => "DP_chrX" }, 
							subset => [qw|IID DP_chrX|] },
	'site.chrY.idepth' => { alias => { INDV => "IID", MEAN_DEPTH => "DP_chrY" }, 
							subset => [qw|IID DP_chrY|] },
	'geno.het' =>         { alias => { F => 'F_auto' }, subset => [qw|IID F_auto|] },
	'geno.ibc' => 		  { alias => { Fhat1 => 'F1_auto', Fhat2 => 'F2_auto', Fhat3 => 'F3_auto' },
							subset => [qw|IID F1_auto F2_auto F3_auto|] },
	#'geno.auto.pass1.imiss' =>  { alias => { F_MISS => "GMiss_all"  }, subset => [qw|IID GMiss_all|] },
	'geno.imiss' =>       { alias => { F_MISS => "GMiss_HQ"  },  subset => [qw|IID GMiss_HQ|] },
	'geno.sexcheck' =>    { alias => { F => 'F_chrX' }, subset => [qw|IID F_chrX|] },
	'geno.imendel' 	=>    { alias => { N => 'MendelErr'}, subset => [qw|IID MendelErr|] },
	'geno.hom.indiv' =>   { alias => { KB => 'HomKB' } , subset => [qw|IID HomKB|] },
	'peddy.het_check.csv' => { alias => { sample_id => "IID", ,het_ratio => 'HetRt_auto', 
										   mean_depth => 'DP_het', 'ancestry-prediction' => 'Popu_pred' },
 							   subset => [qw|IID PC1 PC2 PC3 PC4 DP_het Popu_pred HetRt_auto|] },
	'peddy.sex_check.csv' => { alias => { sample_id => "IID", het_ratio => 'HetRt_chrX', 
 										  ped_sex => 'Sex_ped', predicted_sex => 'Sex_pred' },
 							    subset => [qw|IID Sex_ped Sex_pred HetRt_chrX|] },
 	);

my $tab = data_tab($pedfile, { header => 0,  alias => {1 => "FID", 2 => "IID", 6 => "Pheno"}, 
	subset => [qw|FID IID Pheno|]  });

foreach my $file (sort keys %opts) {
	my $opt = $opts{$file};
	unless (-f "$wrkdir/$file") {
		if ($file eq 'geno.imendel') {
			warn "Cannot find $wrkdir/$file!";
			next;
		}
		else {
			die "Cannot find $wrkdir/$file!";
		}
	}
	print "Slurp $file\n";
	die "Cannot find $wrkdir/$file" unless -f "$wrkdir/$file";
	$tab = $tab->join(data_tab("$wrkdir/$file", { %$opt, chomp => 3 }), 1, ["IID"], ["IID"]);
}

# Add NDP_chrX and NDP_chrY
my (@NdpchrX, @NdpChrY);
for(my $ii = 0; $ii < $tab->nofRow; $ii ++) {
	my $DpChrX = $tab->elm($ii, "DP_chrX");
	my $DpChrY = $tab->elm($ii, "DP_chrY");
	my $DpAuto = $tab->elm($ii, "DP_auto");
	unless(defined $DpAuto) {
		warn "Cannot find autosome depth for ".$tab->elm($ii, "IID");
		print Dumper $tab; exit 1;
		push @NdpchrX, "NA";
		push @NdpChrY, "MA";
	}
	else {
		if ($DpChrX =~ /nan/i) {
			push @NdpchrX, "NA";
		}
		else {
			push @NdpchrX, sprintf("%.3f", $DpChrX/$DpAuto);
		}
		if ($DpChrY =~ /nan/i) {
			push @NdpChrY, "NA";
		}
		else {
			push @NdpChrY, sprintf("%.3f", $DpChrY/$DpAuto);
		}
	}
}

$tab->addCol(\@NdpchrX, "NDP_chrX");
$tab->addCol(\@NdpChrY, "NDP_chrY");

open my $fout, ">", "$outdir/sample.csv" or die "Cannot write to sample file";
print $fout $tab->csv;

open my $fdoc, ">$outdir/sample.README.txt" or die "Cannot write to sample README file";
print $fdoc <<EOF;
DP_auto		Mean depth of coverage of all autosome variants in VCF file (by vcftools)
DP_chrX/chrY 	Mean depth of coverage of all chrX/chrY variants that are not on PAR (by vcftools)
NDP_chrX/chrY 	Normalized depth of coverage of all chrX/chrY variants
DP_het		Mean depth on a subset of well genotyped heterzygotes (by peddy)
F_auto/chrX 	The moment estimate of inbreeding coefficient using autosomal/chrX SNPs (by plink)
F1/F2/F3_auto	Three different version of F_auto estimates based on GCTA (by plink2)
GMiss_HQ 	Genotyping missing rate for all high-quality SNPs that passed genotype QC (by plink)
MendelErr 	Number of Mendel errors for individuals that are part of parent-offspring pairs or trios (by plink2)
HomKB		Total length of runs of homozygosity (in KB) identified in the individual (by plink with default parameters)
HetRt_chrX 	Observed heterzygotes to homozygotes ratio on chrX (by peddy, using a subset of well genotyped exonic SNPs)
Sex_ped 	Sample Gender given in PED file
Sex_pred 	Predicted gender based on HetRt_chrX (by peddy)
HetRt_auto 	Observed heterzygotes to homozygotes ratio on autosomes (by peddy, using a subset of well genotyped exonic SNPs)
PC1-4 	 	First four principle component axes calculated by projecting genotype vector onto the space defined by 1000 genomes populations (by peddy) 
Popu_pred 	Predicted ancestry of origin (one of the major populations in 1000 genomes or unknown (by peddy, using SVM classification)
EOF

