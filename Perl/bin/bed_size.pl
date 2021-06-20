#!/usr/bin/env perl

use strict;
use warnings;
use Genome::Ranges qw|validate_elem|;
use Utils::Number qw|commafy|;
use Utils::File qw|open_file|;

unless(@ARGV) {
	print STDERR "Usage: bed_size BEDFILES ...\n";
}

foreach my $bedfile (@ARGV) {
	if (-f $bedfile || $bedfile eq '-') {
		my $fin;
		if (-f $bedfile) {
			#open $fin, $bedfile or die "Cannot open input file";
			$fin = open_file($bedfile) or die "Cannot open input $bedfile";
		}
		else {
			$fin = \*STDIN;
			$bedfile = "STDIN";
		}
		my $size;
		while(<$fin>) {
			my ($chrom, $st, $ed) = (split)[0,1,2];
			validate_elem($chrom, $st+1, $ed);
			#unless ($chrom =~ /^\w+$/ && $st =~ /^\d+$/ && $ed =~ /^\d+$/) {
			#	next;
			#}
			$size += $ed - $st;
		}
		print commafy($size), "\t", $bedfile, "\n";
	}
	else {
		warn "Cannot find $bedfile";
	}
}
