#!/usr/bin/env perl

use strict;
use warnings;
use Carp;
use File::Basename;
use File::Random qw|random_line|;
use Getopt::Euclid;
use FaSlice;
use Genome::UCSC qw|hg_chrom|;
use Genome::UCSC::TwoBit;
use Utils::Seq qw|is_sym is_snv|;

my %types = (bim => [0, 3, 4, 5]);

my @infds;
if ($ARGV{'--input'} =~ /\.(\w+)$/) {
	my $intype = $1;
	if (!defined $types{$intype}) {
		croak "File type $intype is not supported!";
	}
	@infds = @{$types{$intype}};
}
else {
	croak "Cannot determine the input marker file type: $ARGV{'--input'}";
}


my (@seqfiles, @sqs);
foreach my $ref (@{$ARGV{'--seq'}}) {
	my $sq;
	if ($ref =~ /\.2bit$/) {
		$sq = Genome::UCSC::TwoBit->new($ref);
	}
	elsif ($ref =~ /\.fasta$/ || $ref =~ /\.fa$/) {
		$sq = FaSlice->new(file => $ref);
	}
	else {
		croak "Cannot recognize the type of sequence file: $ref";
	}
	push @sqs, $sq;
	push @seqfiles, basename($ref);
}

my ($tot, %match);
open my $fin, $ARGV{'--input'} or die "Cannot open input marker file";
while(<$fin>) {
	my ($chr, $pos, $a1, $a2) = (split(/\s+/))[@infds];
	my $chrom = hg_chrom($chr);
	(my $chr_ = $chrom) =~ s/^chr//;
	if (is_snv($a1, $a2) && is_sym($a1, $a2)) {
		$tot ++;
		for(my $ii = 0; $ii < @seqfiles; $ii ++) {
			my $refal;
			if ($sqs[$ii]->exists($chrom, $pos)) {
				$refal = $sqs[$ii]->get_slice($chrom, $pos, $pos);
			}
			elsif ($sqs[$ii]->exists($chr_, $pos)) {
				$refal = $sqs[$ii]->get_slice($chr_, $pos, $pos);
			}
			else {
				next;
			}
			if ($refal eq $a1 || $refal eq $a2) {
				$match{$seqfiles[$ii]} ++;
			}
		}
	}
}

if (defined $tot && $tot > 0) {
	print "Total number of symmetric SNPs tested: $tot\n";
	for(my $ii = 0; $ii < @seqfiles; $ii ++) {
		print $match{$seqfiles[$ii]}, " matches alleles of $seqfiles[$ii]\n";
	}
}
else {
	print "No asymmetric SNP was found\n";
}



__END__

=head1 NAME

test_hgbuild -- Quick check the genome reference build for a genetic marker file.

=head1 USAGE

test_hgbuild [options] -in BIMFILE -seq /path/to/hg18.2bit /path/to/hg19.fasta ...

=head1 REQUIRED ARGUMENTS

=over 

=item -[-]in[put] [=] <file>

Input marker map to be checked. Format will be determined from file suffix.

=for Euclid:
	file.type: readable

=item -[-]seq [=] <files>...

One or more reference genome sequence files. 2bit or fasta format.

=for Euclid:
	files.type: readable

=back

=cut

