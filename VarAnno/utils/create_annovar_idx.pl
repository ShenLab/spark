#!/bin/usr/env perl

use strict;
use warnings;
use Carp;
use List::MoreUtils qw|all|;
use Getopt::Euclid;

# Create file index
open my $fin, $ARGV{'<dbtab>'} or croak "Cannot open database table";

my $fout;
my $idxfile = "$ARGV{'<dbtab>'}.idx";

if (-f $idxfile && ! $ARGV{'--force'}) {
	croak "Index file $idxfile already exists!";
}


my @fields = split(',', $ARGV{'--fields'});
unless(@fields == 2) {
	die "Must provide chrom and position!";
}

my $fsep = qr/$ARGV{'--fsep'}/;


my @index;
if (all { /^\d+/ } @fields) {
	@index = map { $_-1 } @fields;
}

my $header;
while($header = <$fin>) {
	if ($header !~ /^#/) {
		# If index already defined: no header line
		if (@index) {
			seek($fin, -length($header), 1);
		}
		else {
			my @fnames = split($fsep, $header);
			for(my $ii = 0; $ii < @fields; $ii ++) {
				$index[$ii] =  first { $fnames[$_-1] eq $fields[$ii] } 1..@fnames;
				unless(defined $index[$ii]) {
					die "Cannot find field $fields[$ii] from input file";
				}
			}
		}
		last;
	}
}

my $BINSZ = $ARGV{'--binsize'};
my $DBSZ = (stat "$ARGV{'<dbtab>'}")[7];

open $fout, ">", $idxfile or croak "Cannot open index file for writing";
print $fout join("\t", "#BIN", $BINSZ, $DBSZ), "\n";

my ($offset0, $offset1);
$offset0 = tell $fin;

my $prev_bin = -1;
my $prev_chrom = "";
my ($chrom, $pos);
while(<$fin>) {
	my $linelen = length($_);
 	($chrom, $pos) = (split($fsep))[@index];
	my $bin = $pos - ($pos % $BINSZ);
	if ($chrom ne $prev_chrom || $bin != $prev_bin) {
		if ($prev_bin < 0) {
			$offset0 = tell $fin; $offset0 -= $linelen;
			$prev_bin= $bin;
		}
		else {
			$offset1 = tell $fin;
			$offset1 -= $linelen;
			print $fout join("\t", $prev_chrom, $prev_bin, $offset0, $offset1), "\n";
			$prev_bin = $bin;
			$offset0 = $offset1;
			$offset1 = undef;
		}
		$prev_chrom = $chrom if $chrom ne $prev_chrom;
	}
}

$offset1 = tell $fin;
print $fout join("\t", $chrom, $prev_bin, $offset0, $offset1), "\n";

__END__

=head1 NAME

create_annovar_idx -- Create ANNOVAR database file index for fast searching.

=head1 DESCRIOTION

The standard ANNOVAR db table will have Chr/Start/End/Ref/Alt as the first five columns,
and sorted based on Chr/Start. We build index based these on Chr/Start.
The results are validated by comparing with idx files provided by ANNOVAR.

We extend this script to any index any line-based database file with or without header.
User should specify the field names for chromosome and position.

=head1 REQUIRED ARGUMENTS

=over

=item <dbtab>

The database table file.

=for Euclid:
  dbtab.type: readable

=back

=head1 OPTIONS

=over

=item -[-]bin[size] [=] <size>

The bin size.

=for Euclid:
  size.default: 100

=item -[-]fields [=] <string>

Field names for chrom and position. If no header in dbfile, use line number (1-based) as field name.

=for Euclid:
	string.default: "1,2"

=item -[-]fsep [=] <regex>

Field separator. Default is "\t".

=for Euclid:
	regex.default: "\t"

=item -[-]skip [=] <comment>

Skip lines beginning with comment string.

=for Euclid:
	comment.default: '#'

=item -[-]f[orce]

Over-write existing index.

=back

=cut
