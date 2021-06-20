use strict;
use warnings;
use List::Util qw|first|;
use List::MoreUtils qw|uniq|;
use Getopt::Euclid;
use Utils::List qw|parse_fields|;
use Utils::Hash qw|peek_hash|;
use Utils::File::Iter qw|iter_file|;
use Genome::UCSC qw|hg_chrom|;
use Genome::UCSC::TwoBit;
use FaSlice;

# Sequence query object
my $sq;
if ($ARGV{'--seq'} =~ /\.2bit$/) {
	$sq = Genome::UCSC::TwoBit->new($ARGV{'--seq'});
}
else {
	#croak "Currently only support .2bit file";
	$sq = FaSlice->new(file => $ARGV{'--seq'});
}
my $refchr;
{
	my $firstchr = peek_hash($sq->{SEQLEN});
	if ($firstchr =~ /^chr/) {
		$refchr = 1;
	}
	else {
		$refchr = 0;
	}
}

my ($it, $fnames) = iter_file($ARGV{'--input'}, { fsep => qr/\t/ });

my (@infields, $infmt);
{
	@infields = parse_fields($ARGV{'--in-fields'}); 
	if (@infields == 5) {
		$infmt = 'ANNOVAR';
	}
	else {
		die "Incorrect number of input fields for a variant";
	}
	foreach my $fd (@infields) {
		unless (grep { $_ eq $fd } @$fnames) {
			die "Cannot find input field $fd from $ARGV{'--input'}!";
		}
	}
}

my @outfields = @$fnames;
my @stdfields;
{
	# Remove in fields, replace standard fields
	# Keep only the first of input fields to splice in standard fields
	for(my $ii = 1; $ii < @infields; $ii ++) {
		my $index = first { $infields[$ii] eq $outfields[$_] } 0..$#outfields;
		splice @outfields, $index, 1;
	}
	#print STDERR join("\t", @outfields), "\n";

	@stdfields = parse_fields($ARGV{'--out-fields'});
	my $index = first { $infields[0] eq $outfields[$_] } 0..$#outfields;
	splice @outfields, $index, 1, @stdfields;
	#print STDERR join("\t", @outfields), "\n";

	my @uqoutfields = uniq sort @outfields;
	unless(scalar(@outfields) == scalar(@uqoutfields)) {
		die "Output file fields are not unique!";
	}
}

my $fout;
if ($ARGV{'--output'}) {
	$fout = IO::File->new($ARGV{'--output'}, "w");
}
else {
	$fout = \*STDOUT;
}
print $fout join("\t", @outfields), "\n";

my $inchr;
while(my $dat = $it->()) {
	if ($infmt eq 'ANNOVAR') {
		my ($chr, $start, $end, $ref, $alt) = @{$dat}{@infields};
		unless(defined $inchr) {
			$inchr = $chr =~ /^chr/ ? 1 : 0;
		}

		if ($ref =~ /^[ACGTacgt]+$/ && $alt =~ /^[ACGTacgt]+$/) {
			$dat->{$stdfields[0]} = $chr;
			$dat->{$stdfields[1]} = $start;
			$dat->{$stdfields[2]} = uc($ref);
			$dat->{$stdfields[3]} = uc($alt);
		}
		else {
			my $chrom = $chr;
			if ($inchr == 1 && $refchr == 0) {
				$chrom =~ s/^chr//; $chrom = 'MT' if $chrom eq 'M';
			}
			elsif ($inchr == 0 && $refchr == 1) {
				$chrom = hg_chrom($chr);
			}
			if ($ref =~ /^[ACGTacgt]+$/ && $alt eq '-') { # deletion
				$start --;
				my $bp = $sq->get_base($chrom, $start);
				$dat->{$stdfields[0]} = $chr;
				$dat->{$stdfields[1]} = $start;
				$dat->{$stdfields[2]} = $bp.uc($ref);
				$dat->{$stdfields[3]} = $bp;
			}
			elsif ($ref eq '-' && $alt =~ /^[ACGTacgt]+$/) { # insertion
				my $bp = $sq->get_base($chrom, $start);
				$dat->{$stdfields[0]} = $chr;
				$dat->{$stdfields[1]} = $start;
				$dat->{$stdfields[2]} = $bp;
				$dat->{$stdfields[3]} = $bp.uc($alt);
			}
			else {
				warn "Cannot find input variant: $chr, $start, $end, $ref, $alt!";
			}
		}
	}
	else {
		die "Cannot recognize input format";
	}
	print $fout join("\t", @{$dat}{@outfields}), "\n";
}



__END__

=head1 NAME

fix_var_chrpos -- Fix variant representation to standard format.

=head1 NOTE

The standard format to represent a variant is to use four columns which represent Chrom, Positon, Ref, Alt.
Ref and Alt must be /^[ACGT]+$/. This is the same format used by VCF, in which variants are represented by
CHROM, POS, REF, ALT.

Variants are represented differently in the output of other software:
	ANNOVAR uses five columns: CHR, START, END, REF, ALT; and for indels REF or ALT is "-".  
	VEP: uses location to store chr/pos and allele for alt allele.

The script will determine the input format based on the number of fields specified for the input, and 
convert the variant to standard representation.

=head1 REQUIRED ARGUMENTS

=over 

=item -[-]in[put] [=] <file>

Input variant table. It must be tab separated.

=item -[-]in-fields [=] <string>

Fields to represent a variant in the input file.

=item -[-]seq [=] <file>

Reference sequence file, in fasta or 2bit format. 

=back

=head1 OPTIONS

=over

=item -[-]out[put] [=] <file>

Output variant table.

=item -[-]out-fields [=] <string>

Standard fields to represent a variant in the output file. 
Field name should not be in conflict with other fields from input.

=for Euclid:
	string.default: "Chrom,Position,Ref,Alt"

=cut
