#!/usr/bin/env perl

use strict;
use warnings;
use Data::Dumper;
use List::MoreUtils qw|any|;
use Utils::File::Iter qw|iter_file|;

# Summarize per-target depth into normalized sex chromosome depth.
my ($infile, $outfile) = @ARGV;

unless(@ARGV == 2) {
	print STDERR "$0 PER_TARGET_DEPTH OUTPUT\n";
	exit 1;
}
unless (-f $infile) {
	die "Cannot find input file $infile";
}

my $fin;
if ($infile =~ /\.gz$/) {
	open $fin, "zcat $infile |" or die "Cannot open gzip file";
}
else {
	open $fin, $infile or die "Cannot open text file";
}


my %stats;
while(<$fin>) {
	my @dat = split;
	my ($chr, $start, $end, $meandp);
	if (@dat == 4) {
		($chr, $start, $end, $meandp) = @dat;
	}
	elsif (@dat == 5) {
		($chr, $start, $end, undef, $meandp) = @dat;
	}
	else {
		die "Incorrect number of fields in the line: $_";
	}
	if ($chr =~ /^\d+$/ || $chr =~ /^chr\d+$/) {
		$stats{auto}{size} += $end - $start;
		$stats{auto}{cumdp} += $meandp * ($end - $start);
		my $chrom = $chr =~ /^\d+$/ ? "chr$chr" : $chr;
		$stats{$chrom}{size} += $end - $start;
		$stats{$chrom}{cumdp} += $meandp * ($end - $start);
	}
	elsif ($chr eq 'X' || $chr eq 'chrX') {
		$stats{chrX}{size} += $end - $start;
		$stats{chrX}{cumdp} += $meandp * ($end - $start);
	}
	elsif ($chr eq 'Y' || $chr eq 'chrY' ) {
		$stats{chrY}{size} += $end - $start;
		$stats{chrY}{cumdp} += $meandp * ($end - $start);
	}
}

my @chroms = map { "chr$_" } 1..22, 'X', 'Y';

if (any { !defined $stats{$_} } qw|auto chrX chrY|) {
	# print Dumper \%stats;
	print STDERR "Target coverage does not include sex chromosomes\n";
#	exit 1;	
}


my $dp_auto = $stats{auto}{cumdp}/$stats{auto}{size};

open my $fout, ">$outfile" or die "Cannot write to $outfile";
print $fout join("\t", qw|DP_auto|, map { "NDP_".$_ } @chroms), "\n";
print $fout $dp_auto;
foreach my $chrom (@chroms) {
	if ($dp_auto > 0 && defined $stats{$chrom}{cumdp}) {
		print $fout "\t", $stats{$chrom}{cumdp}/$stats{$chrom}{size}/$dp_auto;
	}
	else {
		print $fout "\t", "NA";
	}
}
print $fout "\n";


