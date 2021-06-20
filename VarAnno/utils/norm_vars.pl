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
use Utils::Number qw|commafy|;
use Utils::List qw|parse_fields|;
use Utils::File::Iter qw|iter_file|;


my @fields = parse_fields($ARGV{'--fields'});
unless(@fields == 4) {
	croak "Must provide field names for chr, pos, ref, alt";
}

my $header = 1;
if (all { /^\d+$/ } @fields) {
	$header = 0;
}

my $fsep = defined $ARGV{'--fsep'} ? qr/$ARGV{'--fsep'}/ : qr/\s+/;
my ($it, $fnames) = iter_file($ARGV{'--input'},
	{ skip => $ARGV{'--skip'}, sheet => $ARGV{'--sheet'},
 	  header => $header, fsep => $fsep });

my $outfile = $ARGV{'--output'};
if ($outfile =~ /\.txt$/) {
	$outfile =~ s/\.txt$//;
}

open my $fvcf, "| vcf-sort -c > $outfile.vcf" or die "Cannot write VCF";
print $fvcf "##fileformat=VCFv4.0\n";
print $fvcf '#'. join("\t", qw|CHROM POS ID REF ALT QUAL FILTER INFO|), "\n";

while(my $dat = $it->()) {
	my ($chr, $pos, $ref, $alt) = @{$dat}{@fields};
	$pos =~ s/,// if $pos =~ /,/;
	my $varid = sprintf("%s:%d:%s:%s", $chr, $pos, $ref, $alt);
	print $fvcf join("\t", $chr, $pos, $varid, $ref, $alt, 100, 'PASS', 100), "\n";
}

close $fvcf;

my $genomeref = $ARGV{'--seq'};

# normalize variants
system(qq|vt normalize $outfile.vcf -r $genomeref -o $outfile.norm.vcf|);

unless(-f "$outfile.norm.vcf") {
	print STDERR "Cannot normalize variants correctly, check input files!";
	exit 1;
}

# collect results
my %norm;
{
	open my $fin, "$outfile.norm.vcf" or die "Cannot open normalized vcf";
	while(<$fin>) {
		next if /^#/;
		my ($chr, $pos, $varid, $ref, $alt) = (split)[0,1,2,3,4];
		my $newvid = sprintf("%s:%d:%s:%s", $chr, $pos, $ref, $alt);
		if ($newvid ne $varid) {
			$norm{$varid} = $newvid;
		}
	}
}

# Writing to new result file if not empty
if (%norm) {
	print STDERR "Total number of normalized variants: ", scalar(keys %norm), "\n";
	my ($infile, $fout, $flog);
	$infile = $ARGV{'--input'};
	open $fout, ">$outfile.txt" or die "Cannot write to output file";
	open $flog, ">$outfile.log" or die "Cannot write to error file";
	my ($it, $fnames) = iter_file($infile, { skip => $ARGV{'--skip'}, 
		sheet => $ARGV{'--sheet'}, header => $header, fsep => $fsep });
	print $fout join("\t", @$fnames), "\n";
	print $flog join("\t", "Variant", "Normalized"), "\n";
	while(my $dat = $it->()) {
		my ($chr, $pos, $ref, $alt) = @{$dat}{@fields};
		(my $loc = $pos) =~ s/,//;
		my $varid = sprintf("%s:%d:%s:%s", $chr, $loc, $ref, $alt);
		if (defined $norm{$varid}) {
			my @data = split(':', $norm{$varid});
			unless(@data == 4) {
				croak "Incorrect format for normalized variant: $norm{$varid}";
			}
			# Commafy if original position has comma
			if ($pos =~ /,/) {
				$data[1] = commafy($data[1]);
			}
			@{$dat}{@fields} = @data;
			print $fout join("\t", @{$dat}{@$fnames}), "\n"; 
			print $flog join("\t", $varid, $norm{$varid}), "\n";
		}
		else {
			print $fout join("\t", @{$dat}{@$fnames}), "\n"; 
		}
	}
}
else {
	print STDERR "No variant need normalization\n";
}



__END__

=head1 NAME

norm_vars.pl -- Normalize variants (Wrapper around vt norm).

=head1 USAGE

norm_vars.pl [options] -in INPUT -out FILE -fasta HGSEQ -fields CHR,POS,REF,ALT...

=head1 NOTES

This is a wrapper of vt norm. We used this wrapper to validate our own perl-based
implementation.

=head1 REQUIRED ARGUMENTS

=over 

=item -[-]in[put] [=] <file>

Input data file to be checked.

=for Euclid:
  file.type: readable

=item -[-]field[s] [=] <string>

Field names for chromosome, position, reference allele and an optional alternative
allele. The fields should be comma separated. When commma is part of field name
it can be escaped by '\'; then '\' itself can also be escaped by another '\'.
Under VCF to bim file mode, those are extra fields to be extracted and appended.

=item -[-]seq [=] <file>

Reference sequence file, must be in fasta format. 

=for Euclid:
  file.type: readable

=item -[-]out[put] [=] <prefix>

Do not clob the original input file. Write output to new files.

=for Euclid:
  prefix.type: writable

=back

=head1 OPTIONS

=over

=item -[-]sheet [=] <num>

If the input is an Excel workbook, specify the sheet number.

=item -[-]skip [=] <nrow>

Number of rows to skip at the begniing.

=item -[-]fsep [=] <regex>

Field separator.


=back

=cut


