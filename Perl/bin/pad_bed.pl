#!/usr/bin/env perl

use strict;
use warnings;
use List::Util qw|max|;
use Getopt::Lucid qw|:all|;
use Genome::Ranges;
use Genome::Ranges::IntSet;

my ($infile, $padlen, $outfile) = @ARGV;

unless (@ARGV == 2 || @ARGV == 3) {
	print STDERR "Usage:\n";
	print STDERR "$0 INFILE PADLEN [OUTFILE]\n";
	exit 1;
}

unless ($infile eq '-' || -f $infile) {
	die "Cannot find input bed file\n";
}
unless ($padlen =~ /^\d+$/) {
	die "The padding length is not an integer\n";
}

#my @chroms;
my %ranges;
my $fin;
if ($infile eq '-') {
	$fin = \*STDIN;
}
else {
	open $fin, $infile or die "Cannot open bed file";
}
while(<$fin>) {
	my ($chrom, $start, $end) = (split)[0,1,2];
	$start = max(0, $start-$padlen);
	$end = $end + $padlen;
	push @{$ranges{$chrom}}, [$start+1, $end];
	#push @chroms, $chrom;
}


my $rngpad = Genome::Ranges::IntSet->new(\%ranges)->to_ranges();

my $fout;
if ($outfile) {
	open $fout, ">$outfile" or die "Cannot write to output file";
}
else {
	$fout = \*STDOUT;
}

foreach my $chr (sort %$rngpad) {
	foreach my $intv (@{$rngpad->{$chr}}) {
		print $fout join("\t", $chr, $intv->[0]-1, $intv->[1]), "\n";
	}
}
