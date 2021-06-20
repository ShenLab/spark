#!/usr/bin/env perl
use strict;
use warnings;
use FindBin qw|$Bin|;
use Utils::File::Iter qw|iter_file|;

use lib "$Bin/../../lib";
use Variants qw|var_type|;


# Prepare de novo near input
# Input should be a standard annotation file, nearby variants in the same individual should have been removed.
my ($input, $output) = @ARGV;

# List of effects that need fixing so denovoear can recognize
my @effpatch = qw|missense coding_sequence protein_altering splice_acceptor splice_donor frameshift initiator_codon conserved_exon_terminus|;

my ($it, $fnames) = iter_file($input, { fsep => qr/\t/ });

foreach my $fd (qw|Chrom Position Ref Alt GeneID GeneEff|) {
	unless(grep { $fd eq $_ } @$fnames) {
		die "Cannot find field $fd from input file!";
	}
}

open my $fout, ">$output" or die "Cannot write to $output";
print $fout join("\t", qw|gene_name chr pos consequence snp_or_indel|), "\n";
while(my $dat = $it->()) {
	if ($dat->{GeneID} =~ /;/ || $dat->{GeneEff} =~ /;/) {
		die "Gene variants has not been split!";
	}
	my $vartype = var_type($dat->{Ref}, $dat->{Alt});
	if (grep { $_ eq $dat->{GeneEff}  } @effpatch) {
		$dat->{GeneEff} .= "_variant";
	}
	print $fout join("\t", @{$dat}{qw|GeneID Chrom Position GeneEff|}, $vartype), "\n";
}

