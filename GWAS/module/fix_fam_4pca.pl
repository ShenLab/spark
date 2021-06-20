#!/usr/bin/env perl

use strict;
use warnings;
use File::Basename;
use Perl6::Slurp;
use File::Copy qw|copy|;


my ($famfile) = shift @ARGV;

my %fids = map { (split)[1,0] } slurp $famfile;

# Sore population label for each individual in the merged famfile
my %poplabel;
foreach my $label (@ARGV) {
	my $prefix = basename($label);
	$prefix =~ s/\.label$//;
	open my $fin, $label or die "Cannot open $label";
	while(<$fin>) {
		my ($iid, $popu) = split;
		my $newid = "$prefix:$iid";
		unless (defined $fids{$newid}) {
			warn "Cannot find family ID for $newid";
		}
		$poplabel{$newid} = $popu;
	}
}

foreach my $iid (keys %fids) {
	if (!defined $poplabel{$iid}) {
		die "Cannot find population label for $iid";
	}
}

copy($famfile, $famfile.".bak");
open my $fin, $famfile.".bak" or die "Cannot open $famfile.bak for reading";
open my $fout, ">$famfile" or die "Cannot write to $famfile";
while(<$fin>) {
	my @dat = split;
	$dat[5] = $poplabel{$dat[1]};
	print $fout join("\t", @dat), "\n";
} 

