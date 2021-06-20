#!/usr/bin/env perl

use strict;
use warnings;
use Carp;
use Getopt::Euclid;
use Iterator::Simple qw|igrep|;
use Set::IntSpan;
use Genome::UCSC::GeneTable qw|iter_geneTab gdat_list_coding_exons gdat_list_introns|;


my $padding;
if (defined $ARGV{'--padding'}) {
	die "Padding length must be a number" unless $ARGV{'--padding'} =~ /^\d+$/;
	$padding = $ARGV{'--padding'};
}
else {
	$padding = 0;
}

my ($db, $tab) = split('\.', $ARGV{'--dbtab'});

my $iter;
if (-f $ARGV{'--dbtab'}) {
	$iter = igrep { $_->{cdsStart} < $_->{cdsEnd} } iter_genePred($ARGV{'--dbtab'});
}
else {
	$iter = igrep { $_->{cdsStart} < $_->{cdsEnd} } iter_geneTab($db, $tab, 'ucsc');
}
my $fout;
if (defined $ARGV{'--output'}) {
	open $fout, ">$ARGV{'--output'}" or die "Cannot write to output";
}
else {
	$fout = \*STDOUT;
}

while(my $dat = $iter->()) {
	if ($ARGV{'--hgchr'}) {
		next unless $dat->{chrom} =~ /^chr[0-9XY]+$/;
	}
	my @intervals;
	foreach my $intv (gdat_list_coding_exons($dat)) {
		push @intervals, [$intv->[0], $intv->[1]];
	}
  	foreach my $intv (gdat_list_introns($dat)) {
  		push @intervals, [$intv->[0], $intv->[0]+$padding-1];
  		push @intervals, [$intv->[1]-$padding+1, $intv->[1]];
  	}
  	my $rng = Set::IntSpan->new(\@intervals);
  	print $fout $dat->{name}, "\t", $rng->size(), "\n";
}


__END__

=head1 NAME

cds_len : Calculate coding sequence length for each transcript.

=head1 USAGE

cds_len.pl --dbtab hg19.refGene --pad 8 --out refGene.cdsLen.txt

=head1 REQUIRED ARGUMENTS

=over 

=item -[-]dbtab [=] <database>

UCSC genome database and gene table, separated by dot.
Can also be database file in genePred format.

=back

=head1 OPTIONS

=item -[-]out[put] [=] <bedfile>

The output file name, writing to stdout if not provided.
The output intervals will be 0-based half-open, following BED file format.

=item -[-]pad[ding] [=] <number>

Length of extension of coding exon into intronic regions (contain splicing signals).

=for Euclid:
  number.default: 2

=item -[-]hgchr

Focus on canonical chromosomes only.

=cut
