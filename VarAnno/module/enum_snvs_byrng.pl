#!/usr/bin/env perl

use strict;
use warnings;
use Carp;
use Getopt::Euclid;
use FaSlice;
use Genome::UCSC::TwoBit;

# Reference genomc sequence for fetching sequence context
my $sq;
if ($ARGV{'--seq'} =~ /\.2bit$/) {
	$sq = Genome::UCSC::TwoBit->new($ARGV{'--seq'});
}
else {
    # use index fasta sequence otherwise
	$sq = FaSlice->new(file => $ARGV{'--seq'});
}

open my $fout, ">$ARGV{'--output'}" or die "Cannot write to output";
print $fout join("\t", qw|Chrom Position Ref Alt|), "\n";

my $fin = IO::File->new($ARGV{'--bed'});
while(<$fin>) {
	my ($chrom, $st, $ed) = split;
	(my $chr = $chrom) =~ s/^chr//; 
	foreach my $pos ($st+1..$ed) {
		my $ref = $sq->get_slice($chrom, $pos, $pos);
		next unless $ref =~ /^[ACGT]$/;
		foreach my $alt (qw|A C G T|) {
			next if $alt eq $ref;
			print $fout join("\t", $chr, $pos, $ref, $alt), "\n";
		}
	}
}


__END__

=head1 NAME

enum_snvs_by_rng.pl -- Enumerate SNVs on given genomic regions.

=head1 DESCRIPTION

This script will generate all possible SNVs in the given genomic
regions along with their sequence context for further annotation.

=head1 REQUIRED ARGUMENTS

=over 

=item -[-]bed [=] <file>

Bed file containing genomic regions.

=item -[-]out[put] [=] <file>

Output file name prefix, will create a vcf file and an ANNOVAR input file.

=item -[-]seq [=] <refseq>

Reference sequence file. Chromosome name must match with those in the bed file.

=back

=head1 OPTIONS

=over

=item -[-]len [=] <length>

Length of sequence context (default: 3).

=for Euclid:
	length.default: 3

=cut