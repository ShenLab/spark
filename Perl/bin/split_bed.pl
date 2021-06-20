#!/usr/bin/env perl

use strict;
use warnings;
use Carp;
use List::Util qw|min max|;
use Utils::File qw|open_file|;
use Utils::Number qw|commafy|;
use Utils::List qw|which_min|;
use Genome::Ranges::IntSet;
use Getopt::Euclid;

# First calculate the expected length for each split;
my ($totlen, %ranges);
my $fin = open_file($ARGV{'--bed'}) or die "Cannot open $ARGV{'--bed'}";
while(<$fin>) {
	my ($chrom, $start, $end) = (split)[0,1,2];
	$start = max(0, $start-$ARGV{'--pad'});
	$end = $end + $ARGV{'--pad'};
	push @{$ranges{$chrom}}, [$start+1, $end];
}

my $rng = Genome::Ranges::IntSet->new(\%ranges)->to_ranges();

my $nreg = $rng->count;
my $nsplit = $ARGV{'--nsplit'};
if ($nreg < $nsplit) {
	print STDERR "Number of intervals ($nreg) is less than number of splits ($nsplit), consider change nsplit)\n";
	exit 1;
}

my $avelen = $rng->size/$nsplit;
print STDERR "Average length of each split is: ", commafy(sprintf("%.2f", $avelen)), "\n";

my @splits;
my @sizes = (0) x $nsplit;

my $ii = 0;
my $it = $rng->iter();
while(my ($chrom, $start, $end) = $it->()) {
	while($sizes[$ii] > $avelen) {
		$ii ++;
		if ($ii >= $nsplit) {
			$ii = $ii % $nsplit;
		}
	}
	push @{$splits[$ii]}, [$chrom, $start, $end];
	$sizes[$ii] = $end-$start+1;
	$ii ++;
	if ($ii >= $nsplit) {
		$ii = $ii % $nsplit;
	}
}

# Now write to results
my $prefix = $ARGV{'--output'};
for(my $ii = $ARGV{'--start'}; $ii < $nsplit + $ARGV{'--start'}; $ii ++) {
	open my $fout, ">$prefix.$ii.bed" or die "Cannot write to $prefix.$ii.bed";
	foreach my $intv (@{$splits[$ii-1]}) {
		print $fout join("\t", $intv->[0], $intv->[1]-1, $intv->[2]), "\n";
	}
}


__END__

=head1 NOTE

Split genomic interval list into roughly evenly distributed splits.

The script will determine the length of intervals to be appeared in the final interval list based 
on the number of splits. The order of intervals in the original bed file will NOT be preserved.

Output file will be "Prefix.N.bed" with N starting from 1 or otherwise a specified starting index.

=head1 REQUIRED ARGUMENTS

=over

=item -[-]bed [=] <file>

The input BED file, gzipped BED file is supported.

=for Euclid:
	file.type: readable

=item -[-]out[put] [=] <prefix>

The output file prefix.

=back

=head1 OPTIONS

=over

=item -[-]start [=] <number>

The starting index for numbering split files.

=for Euclid:
	number.default: 1
	number.type: integer > 0

=item -[-]pad [=] <length>

Length of interval padding. Interval will be merged if overlap after padding.

=for Euclid:
	length.default: 0

=item -[-]nsplit [=] <number>

Split the table into N small files with similar number of lines.

=back

=cut



