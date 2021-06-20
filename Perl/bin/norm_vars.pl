#!/usr/bin/env perl

use strict;
use warnings;
use Carp;
use IO::Prompt;
use List::MoreUtils qw|all none|;
use FaSlice;
use Utils::Seq qw|rev_comp|;
use Utils::Number qw|commafy|;
use Utils::List qw|parse_fields|;
use Utils::File::Iter qw|iter_file|;
use Genome::UCSC qw|hg_chrom is_hgchr|;
use Genome::UCSC::TwoBit;
use Genome::UCSC::Liftover;
use Genet::Var qw|normalize|;
use Getopt::Euclid;


my @fields = parse_fields($ARGV{'--fields'});
unless(@fields == 4) {
	croak "Must provide field names for chr, pos, ref, alt";
}

print STDERR "The following fields will be parsed for chr, pos, ref, alt:\n", 
	join("\t", @fields), "\n";

# Iterator to variant file
my $noheader;
if (all { /^\d+$/ } @fields) {
	$noheader = 1;
}
my $fsep = defined $ARGV{'--fsep'} ? qr/$ARGV{'--fsep'}/ : qr/\s+/;
my ($it, $fnames) = iter_file($ARGV{'--input'} eq '-' ? \*STDIN : $ARGV{'--input'},
		{ skip => $ARGV{'--skip'}, sheet => $ARGV{'--sheet'},
 	  	  header => ! $noheader, fsep => $fsep	});

# Check field names
foreach my $field (@fields) {
	if (none { $field eq $_ } @$fnames) {
		print STDERR "Field $field cannot be found in the data file\n";
		exit 1;
	}
}
my ($f_chr, $f_pos, $f_ref, $f_alt) = @fields;


my ($varchr, $refchr);

# Sequence query object
my $sq;
if ($ARGV{'--seq'} =~ /\.2bit$/) {
	$sq = Genome::UCSC::TwoBit->new($ARGV{'--seq'});
	my $firstchr = (sort keys %{$sq->{SEQLEN}})[0];
	if ($firstchr =~ /^chr/) {
		$refchr = 1;
	}
	else {
		$refchr = 0;
	}
}
else {
	#croak "Currently only support .2bit file";
	$sq = FaSlice->new(file => $ARGV{'--seq'});
	open my $fin, $ARGV{'--seq'} or die "Cannot open $ARGV{'--seq'}";
	my $firstline = <$fin>;
	if ($firstline =~ /^>chr/) {
		$refchr = 1;
	}
	else {
		$refchr = 0;
	}
}

if ($ARGV{'--out'} =~ /\.txt$/) {
	$ARGV{'--out'} =~ s/\.txt$//;
}
open my $fout, ">$ARGV{'--out'}.txt" or croak "Cannot write to $ARGV{'--out'}.txt";
print $fout join("\t", @$fnames), "\n";
open my $flog, ">$ARGV{'--out'}.log" or croak "Cannot write to $ARGV{'--out'}.log";

my ($total, $normed) = (0, 0);
while(my $dat = $it->()) {
	my ($chr, $pos, $ref, $alt) = @{$dat}{@fields};
	unless(defined $varchr) {
		if ($chr =~ /^chr/) {
			$varchr = 1;
		}
		else {
			$varchr = 0;
		}
		unless($varchr == $refchr) {
			warn "Chromosome nomenclatures in variant table and sequence file are different";
		}
	}

	$total ++;

	unless ($pos =~ /^[0-9,]+$/ && $ref =~ /^[ACGT]+$/i && $alt =~ /^[ACGT]+$/i) {
		print $flog "$total\tINVALID\t$chr:$pos:$ref:$alt\n";
		next;
	}
	$pos =~ s/,//g if $pos =~ /,/;
	if ($ref =~ /[acgt]/) {
		$ref = uc($ref);
		$dat->{$f_ref} = $ref;
	}
	if ($alt =~ /[acgt]/) {
		$alt = uc($alt);
		$dat->{$f_alt} = $alt;
	}
	
	if ($varchr == 1 && $refchr == 0) {
		$chr =~ s/^chr//; $chr = 'MT' if $chr eq 'M';
	}
	elsif ($varchr == 0 && $refchr == 1) {
		$chr = hg_chrom($chr);
	}

	my $var = { CHROM => $chr, POS => $pos, REF => $ref, ALT => $alt };
	normalize($var, $sq);
	unless ($var->{POS} == $pos && $var->{REF} eq $ref && $var->{ALT} eq $alt) {
		$normed ++;
		print $flog "$total\tNORMALIZED\t$chr:$pos:$ref:$alt => ", join(":", @{$var}{qw|CHROM POS REF ALT|}), "\n";
		# Noramlization does not change chr
		#$dat->{$f_chr} = $var->{CHROM};
		$dat->{$f_pos} = $var->{POS};
		$dat->{$f_ref} = $var->{REF};
		$dat->{$f_alt} = $var->{ALT};
	}
	print $fout join("\t", @{$dat}{@$fnames}), "\n";
}
close $flog;

if ($normed > 0) {
	print STDERR "Total ", commafy($total), " variants in the output, ", 
		commafy($normed), " changed allelic representation\n";
}
else {
	print STDERR "Total ", commafy($total), " variants in the output, no change in allelic representation\n";
	unless (-s "$ARGV{'--out'}.log") {
		unlink("$ARGV{'--out'}.log");
		print STDERR "Removed log file\n";
	}
}



__END__

=head1 NAME

norm_vars -- Normalize genetic variants

=head1 USAGE

norm_vars [options] -in INPUT -seq HGSEQ -out OUTPUT

=head1 DESCRIPTION

NOTE: Variant is represented by the four mandatory fields for CHR:POS:REF:ALT.
Both REF and ALT alleles should match regex /^[ACGT]+$/. Deletions/Insertions
should follow the representation used in VCF file. 

=head1 REQUIRED ARGUMENTS

=over 

=item -[-]in[put] [=] <file>

Input data file for liftover. Use '-' to read from STDIN.

=item -[-]field[s] [=] <string>

Field names for chromosome, position, reference allele and alternative allele. 
The fields should be comma separated. 

=item -[-]seq [=] <file>

Reference sequence for the target genome, 2bit or fasta format. 

=for Euclid:
  file.type: readable

=item -[-]out[put] [=] <file>

Output file name prefix. 

The main output is the original input file with ref/alt fields modified if a variant is normalized.
It will also output a list of variants whose ref/alt alleles has been changed (.log).

=for Euclid:
  file.type: writable

=back

=head1 OPTIONS

=over

=item -[-]sheet [=] <num>

If input is Excel workbook, specify the sheet number.

=item -[-]skip [=] <nrow>

Number of rows to skip at the begniing.

=item -[-]fsep [=] <regex>

Field separator.

=back

=cut


