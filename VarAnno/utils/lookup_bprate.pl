#!/usr/bin/env perl
use strict;
use warnings;
use FindBin qw|$Bin|;
use Utils::File::Iter qw|iter_file|;
use Getopt::Euclid;

use lib "$Bin/../../lib";
use Shared qw|parse_tabfile|;
use Variants qw|mutrater|;

my $rater;
{
	my @fields = qw|method lookup fasta rate weight scale chunk|;
	my %deftab = ('3merDenovoNear' => "$ENV{HOME}/Dropbox/Data/Genetics/MutRate_3merDenovoNear.txt", 
				  '7merMrEel' => "$ENV{HOME}/Dropbox/Data/Genetics/MutRate_7merMrEel.txt" );
	my %opt;
	foreach my $fd (@fields) {
		my $Fd = ucfirst($fd);
		if ($fd eq 'fasta') {
			$opt{$Fd} = $ARGV{'--fasta'} // $ARGV{'--seq'};
		}
		else {
			$opt{$Fd} = $ARGV{"--$fd"};
		}
	}
	unless(defined $deftab{$opt{Method}}) {
		die "Local mutation rate method $opt{Method} is not supported!";
	}
	unless(defined $opt{Lookup}) {
		$opt{Lookup} = $deftab{$opt{Method}};
	}
	$rater = mutrater(\%opt);
}

my ($it, $fnames, $keyfields) = parse_tabfile($ARGV{'--input'}, $ARGV{'--in-fields'}, 4, 5);
if (@$keyfields == 4) {
	unless($ARGV{'--seq'}) {
		die "Must provide genome refseq if context if not available in the input!";
	}
}

my @outfields = @$fnames;
my @fieldsadd;
{
	@fieldsadd = split(',', $ARGV{'--field-add'});
	unless(@fieldsadd == 1 || @fieldsadd == 2) {
		die "Can only added one or two additional fields!";
	}
	foreach my $fd (@fieldsadd) {
		if (grep { $fd eq $_ } @outfields) {
			die "The added field $fd already exists in the input!";
		}
	}
	push @outfields, @fieldsadd;
}

my $fout;
if ($ARGV{'--output'}) {
	open $fout, ">$ARGV{'--output'}" or die "Cannot write to $ARGV{'--output'}";
}
else {
	$fout = \*STDOUT;
}
print $fout join("\t", @outfields), "\n";

my ($prev_chrom, $prev_pos);
while(my $dat = $it->()) {
	my ($chrom, $pos, $ref, $alt, $context) = @{$dat}{@$keyfields};
	if ( defined $ARGV{'--chunk'} && defined $prev_chrom && $prev_chrom eq $chrom) {
		if ($pos < $prev_pos) {
			die "Input file is not sorted by Chrom,Position when chunk is enabled";
		}
	}
	if ($ref =~ /^[ACGT]$/ && $alt =~ /^[ACGT]$/) {
		my %site = (Chrom => $chrom, Position => $pos, Ref => $ref, Alt => $alt);
		if (defined $context) {
			$site{Context} = $context;
		}
		my ($prob, $fromto) = $rater->(\%site);
		if (@fieldsadd == 2) {
			print $fout join("\t", @{$dat}{@$fnames}, $fromto, $prob), "\n";
		}
		else {
			print $fout join("\t", @{$dat}{@$fnames}, $prob), "\n";
		}
	}
	else {
		if (@fieldsadd == 2) {
			print $fout join("\t", @{$dat}{@$fnames}, ".", "."), "\n";
		}
		else {
			print $fout join("\t", @{$dat}{@$fnames}, "."), "\n";
		}
	}
	$prev_chrom = $chrom;
	$prev_pos = $pos;
}


__END__


=head1 NAME

lookup_bprate.pl -- Looking up per-bp mutation rate from sequence context.

=head1 NOTES

This utility can be used to fetch per-bp allele-specific mutation rate for each given variant.
Currently, we only have models for calculating SNV rate based on sequence context and scaling factor.
For indels/MNVs, a missing value will be assigned.

=head1 REQUIRED ARGUMENTS

=over

=item -[-]in[put] [=] <table>

Input table. It should be tab separated containing genomic cooridnates and ref/alt alleles.

=back

=head1 OPTIONS

=over

=item -[-]in-fields [=] <fields>

Standard input file fields, should correspond to Chrom,Position,Ref,Alt,(Context).
No rename is needed!

=for Euclid:
	fields.default: "Chrom,Position,Ref,Alt,Context"

=item -[-]out[put] [=] <table>

Output table file name. If not provided, will dump to STDOUT. Fields will be tab separated.

=item -[-]field[s]-add [=] <string>

Additional field(s) added to the input. If two fields are specified, should corresponds to Context and MutRate.

=for Euclid:
	string.default: "MutProb"

=item -[-]method [=] <string>

Method for local mutation rates from sequence context. Two methods are supported:
3merDenovoNear and 7merMrEel.

=for Euclid:
	string.default: "3merDenovoNear"

=item -[-]lookup [=] <table>

Lookup table for local mutation rates, the default will be automatically set based on method. 

=item -[-]seq [=] <seq>

Genome reference sequence file used to fetch sequence context if context cannot be found from input.
Can be either fasta or 2bit format.

=item -[-]fasta [=] <seq>

Alias to --seq.

=item -[-]rate [=] <db>

Tabix indexed database table of per-bp allele specific mutation rate.

=item -[-]weight [=] <db>

Tabix indexed database table of per-bp weight.

=item -[-]scale [=] <factor>

Either a constant or a bedgraph file to scale the mutation rate.

=for Euclid:
	factor.default: 1

=item -[-]chunk [=] <size>

Chunk size for storing per-bp weight and rates. When this is enabled, instead of querying weight/rate
at each base pair, we can directly fetch the weight/rate from RAM. It also requires that the input
file should be sorted by Chrom and Position.

=back




