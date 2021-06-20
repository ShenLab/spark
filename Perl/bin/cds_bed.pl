#!/usr/bin/env perl

use strict;
use warnings;
use Carp;
use Perl6::Slurp;
use Getopt::Euclid;
use Iterator::Simple qw|igrep|;
use Utils::Number qw|commafy|;
use Genome::UCSC::GeneTable qw|iter_geneTab iter_genePred gdat_list_coding_exons gdat_list_introns|;
use Genome::Ranges::IntSet;

my $padding;
if (defined $ARGV{'--padding'}) {
	die "Padding length must be a number" unless $ARGV{'--padding'} =~ /^\d+$/;
	$padding = $ARGV{'--padding'};
}
else {
	$padding = 0;
}

my (%gincl, %gexcl, %tincl, %texcl);
if ($ARGV{'--include-genes'}) {
	%gincl = map { (split)[0] => 1 } slurp $ARGV{'--include-genes'};
}
if ($ARGV{'--include-trans'}) {
	%tincl = map { (split)[0] => 1 } slurp $ARGV{'--include-trans'};
}
if ($ARGV{'--exclude-genes'}) {
	%gexcl = map { (split)[0] => 1 } slurp $ARGV{'--exclude-genes'};
}
if ($ARGV{'--exclude-trans'}) {
	%texcl = map { (split)[0] => 1 } slurp $ARGV{'--exclude-trans'};
}


my ($db, $tab) = split('\.', $ARGV{'--dbtab'});

my %intervals;
my $iter;
if (-f $ARGV{'--dbtab'}) {
	$iter = igrep { $_->{cdsStart} < $_->{cdsEnd} } iter_genePred($ARGV{'--dbtab'});
}
else {
	$iter = igrep { $_->{cdsStart} < $_->{cdsEnd} } iter_geneTab($db, $tab, 'ucsc');
}
while(my $dat = $iter->()) {
	if ($ARGV{'--hgchr'}) {
		next unless $dat->{chrom} =~ /^chr[0-9XY]+$/;
	}
	if (%tincl) {
		next unless defined $tincl{$dat->{name}};
	}
	if (%texcl) {
		next if defined $texcl{$dat->{name}};
	}
	if (%gincl) {
		next unless defined $gincl{$dat->{name2}};
	}
	if (%gexcl) {
		next if defined $gexcl{$dat->{name2}};
	}

	if ($ARGV{'--nochr'}) {
		$dat->{chrom} =~ s/^chr//;	
	}
	my @cds = gdat_list_coding_exons($dat);
	foreach my $intv (@cds) {
		push @{$intervals{$dat->{chrom}}}, [$intv->[0], $intv->[1]];
	}
  	my @introns = gdat_list_introns($dat);
  	foreach my $intv (@introns) {
  		push @{$intervals{$dat->{chrom}}}, [$intv->[0], $intv->[0]+$padding-1];
  		push @{$intervals{$dat->{chrom}}}, [$intv->[1]-$padding+1, $intv->[1]];
  	}
}

my $extrng = Genome::Ranges::IntSet->new(\%intervals);
my $totlen = $extrng->size - $extrng->count;
print STDERR "Total size of merged intervals: ", commafy($totlen), "\n";

my $fout;
if (defined $ARGV{'--output'}) {
	open $fout, ">$ARGV{'--output'}" or die "Cannot write to output";
}
else {
	$fout = \*STDOUT;
}

my $it = $extrng->to_ranges()->iter();
while(my $dat = $it->()) {
	my ($chrom, $start, $end) = @{$dat}[0,1,2];
	print $fout join("\t", $chrom, $start-1, $end), "\n";
}


__END__

=head1 USAGE

cds_bed.pl --dbtab hg19.refGene --pad 8 --nochr --out refGene.b37.bed

=head1 NOTES

Coding regions are be defined as CDS part of coding exons plus first and last N
bases of introns that presumed to contain splice signals (N=2 by default, but can)

ome extension
from coding exons into flanking intronic regions (defined by a padding parameter
default is 2). 

**: UTR
||: coding sequence
--: introns

    ******-------***|||||||||-------||||||||-------||||********
                    ...........   ............   ......
          ..   ..   ...........   ............   ......


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

=item -[-]nochr

Strip off the chr prefix of chromosome names.

=item -[-]hgchr

Focus on canonical chromosomes only.

=item -[-]include-trans [=] <translist>

A list of transcript IDs to be included.

=item -[-]exclude-trans [=] <translist>

A list of transcript IDs to be included.

=item -[-]include-genes [=] <genelist>

A list of gene IDs to be included.

=item -[-]exclude-genes [=] <genelist>

A list of gene IDs to be included.

=cut
