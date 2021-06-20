#!/usr/bin/env perl

use strict;
use warnings;
use Perl6::Slurp;
use Data::Dumper;
use Genet::Ped;
use List::Util qw|min max|;
use Utils::List qw|all_pairs|;
use Utils::File qw|count_line|;
use Utils::File::Iter qw|iter_file slurp_file|;


my ($pedfile, $wrkdir, $outdir, $suffix, $thres) = @ARGV;

unless(@ARGV == 3 || @ARGV == 4 || @ARGV == 5) {
	print STDERR <<EOF

Purpose: Used by samp_qc to merge pairwise kinship coefficient

Usage: merge_kinship_est PEDFILE WRKDIR OUTDIR [SUFFIX] [THRES]

Pedigree structure, sample sex and phenotypes are defined in pedigree file.
Raw output from king, plink, and peddy should be found in working dir.
Output will be csv format, contain all information about sample quality,
gender, etc.

Users can optionally provide a SUFFIX string to define duplicated sample.

Users can also provide a Kinship threshold to filter kin0 file, useful for large data set.

EOF
}
else {
	die "Cannot find ped file" unless -f $pedfile;
	die "Cannot find working directory" unless -d $wrkdir;
	die "Cannot find output directory" unless -d $outdir;
}

# First searching for duplicated sample
# They will also be ignored when reading PED file
my $ignore;
my %dups;
if (defined $suffix) {
	$ignore = qr/$suffix\d*$/;
	my %samps = map { (split)[1] => 1 } slurp $pedfile;
	my @matchdup = grep { /$ignore/ } keys %samps;
	foreach my $iid (@matchdup) {
		(my $prefix = $iid) =~ s/$ignore//;
		if (defined $samps{$prefix}) {
			push @{$dups{$prefix}}, $iid;
		}
		else {
			print STDERR "The original sampple of $iid does not exist in PED file\n";
		}
	}
}

# Merge all useful information about relative pairs into one file.

my $ped = Genet::Ped->new($pedfile, {ignore => $ignore});
my %relpairs;
while(my ($iid, $dupcat) = each %dups) {
	foreach my $pair (all_pairs($iid, @$dupcat)) {
		$relpairs{$pair->[0],$pair->[1]} = "Duplicate";
		$relpairs{$pair->[1],$pair->[0]} = "Duplicate";
	}
}
foreach my $pairs (values %{ $ped->get_po_pairs() }) {
	foreach my $pair (@$pairs) {
		my @samp = ([$pair->[0]], [$pair->[1]]);
		for my $ii (0,1) {
			if (defined $dups{$pair->[$ii]}) {
				push @{$samp[$ii]}, @{$dups{$pair->[$ii]}};
			}
		}
		for(my $jj = 0; $jj < @{$samp[0]}; $jj ++) {
			for(my $kk = 0; $kk < @{$samp[1]}; $kk ++) {
				$relpairs{$samp[0][$jj],$samp[1][$kk]} = "Parent-Offspring";
				$relpairs{$samp[1][$kk],$samp[0][$jj]} = "Parent-Offspring";
			}
		}
	}
}
foreach my $pairs (values %{ $ped->get_sib_pairs() }) {
	foreach my $pair (@$pairs) {
		my @samp = ([$pair->[0]], [$pair->[1]]);
		for my $ii (0,1) {
			if (defined $dups{$pair->[$ii]}) {
				push @{$samp[$ii]}, @{$dups{$pair->[$ii]}};
			}
		}
		for(my $jj = 0; $jj < @{$samp[0]}; $jj ++) {
			for(my $kk = 0; $kk < @{$samp[1]}; $kk ++) {
				$relpairs{$samp[0][$jj],$samp[1][$kk]} = "Sib-Pair";
				$relpairs{$samp[1][$kk],$samp[0][$jj]} = "Sib-Pair";
			}
		}
	}
}

#my $totline = count_line("$wrkdir/geno.kin0.gz");
#my $maxline = 1_000_000; # Capped at 1M lines
#my $step = max(int($totline/$maxline), 1);
#if ($step > 1) {
#	print STDERR "Trim between-family data by a factor of $step\n";
#}

my %kin;
{
	my $iter = iter_file("$wrkdir/geno.kin");
	while(my $dat = $iter->()) {
		$dat->{FID1} = $dat->{FID};
		$dat->{FID2} = $dat->{FID};
		$kin{$dat->{ID1},$dat->{ID2}} = $dat;
	}
 
	$iter = iter_file("$wrkdir/geno.kin0.gz");
	while(my $dat = $iter->()) {
		if (defined $thres) {
			next unless $dat->{Kinship} > $thres;
		}
		$kin{$dat->{ID1},$dat->{ID2}} = $dat; 
	}
}


foreach my $file (qw|geno.ibs geno.ibs0.gz|) {
	my $it = iter_file("$wrkdir/".$file);
	while(my $dat = $it->()) {
		next unless defined $kin{$dat->{ID1}, $dat->{ID2}} || defined $kin{$dat->{ID2}, $dat->{ID1}};
		my $info  = $kin{$dat->{ID1}, $dat->{ID2}} // $kin{$dat->{ID2}, $dat->{ID1}};
		unless (defined $dat->{FID1}) {
			if (defined $relpairs{$dat->{ID1}, $dat->{ID2}}) {
				$info->{TYPE} = $relpairs{$dat->{ID1}, $dat->{ID2}};
			}
			else {
				$info->{TYPE} = "Within-Family";
			}	
		}
		else {
			$info->{TYPE} = "Between-Family";
		}
		foreach my $k (qw|N_SNP N_IBS0 N_IBS1 N_IBS2|) {
			$info->{$k} = $dat->{$k};
		}
	}
}

my $it = iter_file("$wrkdir/geno.genome.gz");
while(my $dat = $it->()) {
	next unless defined $kin{$dat->{IID1}, $dat->{IID2}} || defined $kin{$dat->{IID2}, $dat->{IID1}};
	my $info  = $kin{$dat->{IID1}, $dat->{IID2}} // $kin{$dat->{IID2}, $dat->{IID1}};
	next unless defined $info;
	foreach my $k (qw|PI_HAT Z0 Z1 Z2|) {
		$info->{$k} = $dat->{$k};
	} 
}

my %peddy = (n => 'n_snp', ibs0 => 'n_ibs0', ibs2 => 'n_ibs2');
$it = iter_file("$wrkdir/peddy.ped_check.csv");
while(my $dat = $it->()) {
	next unless defined $kin{$dat->{sample_a},$dat->{sample_b}} || defined $kin{$dat->{sample_b},$dat->{sample_a}};
	my $info = $kin{$dat->{sample_a},$dat->{sample_b}} // $kin{$dat->{sample_b},$dat->{sample_a}};
	while(my ($k, $a) = each %peddy) {
		$info->{$a} = $dat->{$k};
	}
	$info->{n_ibs1} = $info->{n_snp}-$info->{n_ibs0}-$info->{n_ibs2};
	$info->{ibs0} = sprintf("%.4f", $info->{n_ibs0}/$info->{n_snp});
	$info->{kinship} = $dat->{rel};
}

my @fields = qw|FID1 ID1 FID2 ID2 TYPE Kinship IBS0 N_SNP N_IBS0 N_IBS1 N_IBS2
				PI_HAT Z0 Z1 Z2 kinship ibs0 n_snp n_ibs0 n_ibs1 n_ibs2|;
open my $fout, ">$outdir/relpairs.csv" or die "Cannot write to relpairs.csv.gz";
print $fout join(",", @fields), "\n";
foreach my $dat (values %kin) {
	print $fout join(',', @{$dat}{@fields}), "\n";
}

open my $fdoc, ">$outdir/relpairs.README.txt" or die "Cannot write to relpairs.README.txt";
print $fdoc <<EOF;
FID1/ID1 	Family ID and sample ID for individual 1
FID2/ID2 	Family ID and sample ID for individual 2
TYPE		Type of the individual pair (Parent-Offspring, Sib-Pair, Within-Family, and Between-Family)
Kinship 	Estimated kinship coefficient (by king v2 using HQ genotypes, NOTE: this is roughly half of plink's PI_HAT)
IBS0 		Proportion of markers did not share any allele identical-by-state betwen the pair (by king)
N_SNP 		Total number of SNPs in HQ genotype data set
N_IBS0/1/2 	The number of SNPs in HQ genotype data set that share 0, 1, and 2 alleles identical-by-state betwen the pair (by king)
PI_HAT 	 	Estimated coefficient of relatedness (by plink using HQ genotypes)
Z0/Z1/Z2	Estimated proportion of markers that share 0, 1, and 2 alleles identical-by-descent
kinship 	Estimated kinship coefficient * 2 (by a modified version of king implemented by peddy, using a subset of HQ genotypes)
ibs0 		Proportion of markers did not share any allele identical-by-state betwen the pair (by peddy, using a subset of HQ genotypes)
n_ibs0/1/2	The number of SNPs in HQ genotype data set that share 0, 1, and 2 alleles identical-by-state betwen the pair (by peddy, using a subset of HQ genotypes)
EOF

