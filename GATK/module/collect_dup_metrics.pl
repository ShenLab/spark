#!/usr/bin/env perl

use strict;
use warnings;
use Carp;
use IO::Dir;
use File::Basename;

# Collect BamQC statistics
my ($wrkdir, $outdir) = @ARGV;
unless(@ARGV == 2) {
	print STDERR "$0 WRKDIR OUTDIR\n";
	exit 1;
}
unless (-d $wrkdir && -d $outdir) {
	die "Cannot read working dir or output dir";
}

my (%metrics);
my @iids = grep { !/^\./ && -d "$wrkdir/$_" } IO::Dir->new($wrkdir)->read();
foreach my $iid (@iids) {
	croak "Cannot find metrics file for $iid"
		unless -f "$wrkdir/$iid.dup_metrics";
	open my $fin, "$wrkdir/$iid.dup_metrics" or die "Cannot open dup_metrics";
	my @fields;
	while(<$fin>) {
		if (/^\w+/) {
			chomp();
			@fields = split /\t/;
			last;
		}
	}
	my %dat;
	while(<$fin>) {
		chomp();
		my @values = split /\t/;
		@dat{@fields} = @values;
		last;
	}
	croak "Cannot find dup metrics" unless defined $dat{PERCENT_DUPLICATION};
	$metrics{$iid} = $dat{PERCENT_DUPLICATION};
}

open my $fout, ">$outdir/dup_metrics.csv";
print $fout join(',', qw|IID PCT_DUP|), "\n";
foreach my $iid (@iids) {
	print $fout join(',', $iid, $metrics{$iid}), "\n";
}