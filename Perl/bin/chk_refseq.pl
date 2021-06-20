#!/usr/bin/env perl

use strict;
use warnings;
use Carp;
use Getopt::Euclid;
use IO::File;
use IO::Prompt;
use File::Basename;
use List::MoreUtils qw|any none all|;
use FaSlice;
use File::Copy qw|move|;
use Data::Dumper;
use Genome::UCSC qw|hg_chrom|;
use Genome::UCSC::TwoBit;
use Utils::Number qw|commafy|;
use Utils::Seq qw|rev_comp|;
use Utils::List qw|parse_fields|;
use Utils::Hash qw|peek_hash|;
use Utils::File::Iter qw|iter_file|;


my @fields = parse_fields($ARGV{'--fields'});
unless(@fields == 3) {
	croak "Must provide field names for chr, pos, ref";
}

print STDERR "The following fields will be parsed for chr, pos, ref:\n", 
	join("\t", @fields), "\n";

# Iterator to the file
my %opt;
if (defined $ARGV{'--fsep'}) {
	$opt{fsep} = qr/$ARGV{'--fsep'}/;
}
if (all { /^\d+$/ } @fields) {
	$opt{header} = 0;
}
else {
	$opt{header} = 1;
}
#my $fsep = defined $ARGV{'--fsep'} ? qr/$ARGV{'--fsep'}/ : qr/\s+/;
my ($it, $fnames) = iter_file($ARGV{'--input'} eq '-' ? \*STDIN : $ARGV{'--input'},
	{ skip => $ARGV{'--skip'}, sheet => $ARGV{'--sheet'},maxcol => $ARGV{'--maxcol'}, maxrow => $ARGV{'--maxrow'}, %opt });

# Check field names
foreach my $field (@fields) {
	if (none { $field eq $_ } @$fnames) {
		print STDERR "Field $field cannot be found in the data file\n";
		exit 1;
	}
}
my ($f_chr, $f_pos, $f_ref) = @fields;


my ($varchr, $refchr);

# Reference seq file handler
my $sq;
if ($ARGV{'--seq'} =~ /\.2bit$/) {
	$sq = Genome::UCSC::TwoBit->new($ARGV{'--seq'});
	my $firstchr = peek_hash($sq->{SEQLEN}); # (sort keys %{$sq->{SEQLEN}})[0];
	if ($firstchr =~ /^chr/) {
		$refchr = 1;
	}
	else {
		$refchr = 0;
	}
}
elsif ($ARGV{'--seq'} =~ /\.fasta$/ || $ARGV{'--seq'} =~ /\.fa$/) {
	# use index fasta sequence otherwise
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
else {
	croak "Cannot recognize the type of sequence file: $ARGV{'--seq'}";
}
print STDERR "Start checking reference seq\n";


# Now parse the file, the first four columns will assumed to be
# chromsome, position, ref allele, and alt allele
my ($total, $right, $lineno) = (0, 0, 0);
my %badsites;
while(my $dat = $it->()) {
	$lineno ++;
	my @data = @{$dat}{@fields};
	my ($chr, $pos, $ref) = @data[0,1,2];
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

	# check if desired data fields are captured
	unless ($chr =~ /^\w+$/ && $pos =~ /^[0-9,]+$/ && $ref =~ /^[ACGT]+$/) {
		$badsites{$chr,$pos,$ref} = [$lineno, "INVALID_LINE"];
		if ($ARGV{'--verbose'}) {
			print STDERR "Ignore invalid line: ", join("\t", $chr, $pos, $ref), "\n";
			exit 1 if $ARGV{'--strict'};
		}
		next;
	}
	(my $loc = $pos) =~ s/,//g;

	$total ++;

	# We need to check if chr exist in the reference sequence file
	if ($varchr == 1 && $refchr == 0) {
		$chr =~ s/^chr//; $chr = 'MT' if $chr eq 'M';
	}
	elsif ($varchr == 0 && $refchr == 1) {
		$chr = hg_chrom($chr);
	}

	unless ($sq->exists($chr, $loc, $loc+length($ref)-1)) {
		$badsites{$chr,$pos,$ref} = [$lineno, "OUT_OF_RANGE"];
		next;
	}

	my $refsq = $sq->get_slice($chr, $loc, $loc+length($ref)-1);
	unless ($ref eq $refsq) {
		if ($ref eq rev_comp($refsq)) {
			$badsites{$chr,$pos,$ref} = [$lineno, "REVCOMP_MATCH"];
		}
		else {
			$badsites{$chr,$pos,$ref} = [$lineno, "REFSEQ_MISMATCH"];
		}
		next;
	}

	$right ++;
}

print STDERR "Total ", commafy($total), " valid positions checked: ", commafy($right), " are ref matched.\n";
if (%badsites) {
	my $fout;
	if ($ARGV{'--output'}) {
		open $fout, ">$ARGV{'--output'}" or die "Cannot write to output file";
	}
	else {
		$fout = \*STDOUT;
	}
	foreach my $loc (sort { $badsites{$a}[0] <=> $badsites{$b}[0] } keys %badsites) {
		my ($chr, $pos, $ref) = split($;, $loc);
		my ($lineno, $reason) = @{$badsites{$loc}};
		print $fout join("\t", "LINE$lineno: ".$reason, $chr, $pos, $ref), "\n";
	}
}



__END__

=head1 NAME

chk_refseq -- Check reference sequence given genomic positions.

=head1 USAGE

chk_refseq [options] -in INPUT -out FILE -seq HGSEQ -fields CHR,POS,REF ...

=head1 DESCRIPTION

The required input file should contain information about chromosome name and
starting position (1-based) for a reference allele. The utility will check if 
the reported reference allele matches that in the  reference genome. Usually, 
this will be the first step for variant annotation.

Indels are represented differently by different tools. Here we took the same
format as used by VCF standard.

The following types of errors will be reported:
INVALID - the input file does not have valid chr/pos/ref fields.
MISMATCH - the reported ref does match ref seq.

When error is found, the input file will be backed up, then only correct lines
will be written to the new file, errors will be written in a separate file.

=head1 REQUIRED ARGUMENTS

=over 

=item -[-]in[put] [=] <file>

Input data file to be checked.

=for Euclid:
  file.type: readable

=item -[-]field[s] [=] <string>

Field names for chromosome, position, and reference allele.
The fields should be comma separated. When commma is part of field name
it can be escaped by '\'; then '\' itself can also be escaped by another '\'.
Under VCF to bim file mode, those are extra fields to be extracted and appended.

=item -[-]seq [=] <file>

Reference sequence file, 2bit or fasta format. 

=for Euclid:
  file.type: readable

=back

=head1 OPTIONS

=over

=item -[-]out[put] [=] <prefix>

Output file name. If not provided, will dump to STDOUT.

=for Euclid:
  prefix.type: writable

=item -[-]sheet [=] <num>

If the input is an Excel workbook, specify the sheet number.

=item -[-]maxrow [=] <num>

=item -[-]maxcol [=] <num>

For spreadsheet: the max range of columns and rows.

=item -[-]skip [=] <nrow>

Number of rows to skip at the begniing.

=item -[-]fsep [=] <regex>

Field separator, the default is white space.

=item -[-]verbose

Verbose mode

=item -[-]strict

Strict mode, die on invalid line.

=back



