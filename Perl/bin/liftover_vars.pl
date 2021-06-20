#!/usr/bin/env perl

use strict;
use warnings;
use Carp;
use List::MoreUtils qw|all none|;
use FaSlice;
use Utils::Seq qw|rev_comp|;
use Utils::Number qw|commafy|;
use Utils::List qw|parse_fields|;
use Utils::Hash qw|str2hash|;
use Utils::File::Iter qw|iter_file|;
use Genet::Var qw|normalize|;
use Genome::UCSC qw|hg_chrom is_hgchr|;
use Genome::UCSC::TwoBit;
use Genome::UCSC::Liftover qw|check_chain_chrpref|;
use Getopt::Euclid;


# Parse input file fields
my @fields = parse_fields($ARGV{'--fields'});
unless(@fields == 4) {
	croak "Must provide field names for chr, pos, ref, alt";
}
print STDERR "The following fields will be parsed for chr, pos, ref, alt:\n", 
	join("\t", @fields), "\n";
#unless (prompt("Confirm and proceed? ", -yn, -default => 'y')) {
#	print "Exiting...\n";
#	exit 1;
#}

# Iterator to variant file
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
if (defined $ARGV{'--alias'}) {
	$opt{alias} = str2hash($ARGV{'--alias'}, { psep => ',', kvsep => ':' });
}
#my $fsep = defined $ARGV{'--fsep'} ? qr/$ARGV{'--fsep'}/ : qr/\s+/;
my ($it, $fnames) = iter_file($ARGV{'--input'} eq '-' ? \*STDIN : $ARGV{'--input'},
	{ skip => $ARGV{'--skip'}, sheet => $ARGV{'--sheet'}, maxcol => $ARGV{'--maxcol'}, maxrow => $ARGV{'--maxrow'}, %opt });

# Check field names
foreach my $field (@fields) {
	if (none { $field eq $_ } @$fnames) {
		print STDERR "Field $field cannot be found in the data file\n";
		exit 1;
	}
}
my ($f_chr, $f_pos, $f_ref, $f_alt) = @fields;

# Check if old variant field exists
my @outfields = @$fnames;
my $oldvarflag;
if ($ARGV{'--oldvar'}) {
	if ($ARGV{'--oldvar'} =~ /,/) {
		my @oldvarfds = split(',', $ARGV{'--oldvar'});
		unless(@oldvarfds == 4) {
			die "Must provide four fields for oldvar!"
		}
		foreach my $fd (@oldvarfds) {
			if (grep { $_ eq $fd } @$fnames) {
				print STDERR "$fd field already exists!";
				exit 1;
			}
		}
		$oldvarflag = 4;
		unshift @outfields, @oldvarfds;
	}
	else {
		if (grep { $_ eq $ARGV{'--oldvar'} } @$fnames) {
			print STDERR "$ARGV{'--oldvar'} field already exists!\n";
			exit 1;
		}
		$oldvarflag = 1;
		unshift @outfields, $ARGV{'--oldvar'};
	}
}

# Sequence query object
my $sq;
if ($ARGV{'--seq'} =~ /\.2bit$/) {
	$sq = Genome::UCSC::TwoBit->new($ARGV{'--seq'});
}
else {
	#croak "Currently only support .2bit file";
	$sq = FaSlice->new(file => $ARGV{'--seq'});
}

# Determine chr prefix in chain file and create lift chain obj
my ($inchr, $outchr) = check_chain_chrpref($ARGV{'--chain'});

my $lift = Genome::UCSC::Liftover->new($ARGV{'--chain'});
#my ($fromdb, $todb) = $lift->get_fromto;
#if (defined $todb && $ARGV{'--seq'} !~ /$todb/i) {
#	carp "Ref seq file name ($ARGV{'--seq'}) does not match $todb?";
#}

if ($ARGV{'--output'} =~ /\.txt$/) {
	$ARGV{'--output'} =~ s/\.txt$//;
}
open my $fout, ">$ARGV{'--output'}.txt" or croak "Cannot write to $ARGV{'--output'}.txt";
if ($opt{header} == 1) {
	print $fout join("\t", @outfields), "\n";
}
open my $flog, ">$ARGV{'--output'}.log" or croak "Cannot write to $ARGV{'--output'}.log";

my $varchr; # <- flag if chr prefix exist for variants in the table
my ($total, $mapped) = (0, 0);
while(my $dat = $it->()) {
	my ($chr, $pos, $ref, $alt) = @{$dat}{@fields};
	unless(defined $varchr) {
		$varchr = $chr =~ /^chr/ ? 1 : 0;
		unless($varchr == $inchr) {
			warn "Chromosome nomenclatures in variant table and chain files are different";
		}
	}

	my $varid = join(':', $chr, $pos, $ref, $alt);
	#if ($ARGV{'--oldvar'}) {
	#	$dat->{$ARGV{'--oldvar'}} = $varid;
	#}
	if ($oldvarflag) {
		if ($oldvarflag == 1) {
			$dat->{$outfields[0]} = $varid;
		}
		else {
			$dat->{$outfields[0]} = $chr;
			$dat->{$outfields[1]} = $pos;
			$dat->{$outfields[2]} = $ref;
			$dat->{$outfields[3]} = $alt;
		}
	}
	$total ++;

	unless ($chr =~ /^\w\S*$/ && $pos =~ /^[0-9][0-9,]+$/ && $ref =~ /^[ACGTN]+$/ && $alt =~ /^[ACGTN]+$/) {
		print $flog "$total\tINVALID\t$varid\n";
		next;
	}
	$pos =~ s/,//g if $pos =~ /,/;

	#$chr = "chr".$chr if $chr !~ /^chr/;

	if ($varchr == 1 && $inchr == 0) {
		$chr =~ s/^chr//; $chr = 'MT' if $chr eq 'M';
	}
	elsif ($varchr == 0 && $inchr == 1) {
		$chr = hg_chrom($chr);
	}
	
	my ($tg_chrom, $tg_start, $tg_strand, $tg_score) = $lift->query($chr, $pos);
	my ($tg_chrom2, $tg_end, $tg_strand2, $tg_score2);
	if (length($ref) > 1) {
		($tg_chrom2, $tg_end, $tg_strand2, $tg_score2) = $lift->query($chr, $pos+length($ref)-1)
	}
	else {
		($tg_chrom2, $tg_end, $tg_strand2, $tg_score2) = ($tg_chrom, $tg_start, $tg_strand, $tg_score);
	}

	my ($outflag, @tags, $mapinfo);
	if (defined $tg_start && defined $tg_end) {
		my $tg_len = abs($tg_end - $tg_start) + 1;
		# Require: start and end should be on the same chromsome and same strand
		if ( $tg_chrom eq $tg_chrom2 && $tg_strand eq $tg_strand2 ) {
			if ($tg_len == length($ref) || $tg_len == length($alt)) {
				if ($tg_strand eq '+') {
					unless($tg_end >= $tg_start) {
						push @tags, "INVERSION";
						$mapinfo = "$tg_chrom:$tg_start:$tg_end";
					}
					else {
						my $refsq = $sq->get_slice($tg_chrom, $tg_start, $tg_end);
						$mapinfo = "$tg_chrom:$tg_start:$refsq";
						if ($refsq eq $ref) {
							$dat->{$f_chr} = $tg_chrom;
							$dat->{$f_pos} = $tg_start;
							if (is_hgchr($tg_chrom)) {
								$outflag = 1;
								if ($tg_chrom ne $chr) {
									push @tags, "DIFF_CHROM";
								}
							}
							else {
								if ($ARGV{'--no-hgchr'}) {
									$outflag = 1;
								}
								push @tags, "NOT_CANONCHR";
							}
						}
						elsif ($refsq eq $alt) {
							# allele switch after liftover
							$dat->{$f_chr} = $tg_chrom;
							$dat->{$f_pos} = $tg_start;
							$dat->{$f_ref} = $alt;
							$dat->{$f_alt} = $ref;
							if (is_hgchr($tg_chrom)) {
								$outflag = 1;
								push @tags, 'REFALT_SWITCH'
							}
							else {
								if ($ARGV{'--no-hgchr'}) {
									$outflag = 1;
								}
								push @tags, 'NOT_CANONCHR', 'REFALT_SWITCH';
							}
						}
						else {
							push @tags, 'REFALT_MISMATCH';
						}
					}
				}
				else {
					# on negative strand
					unless($tg_end <= $tg_start) {
						push @tags, "INVERSION" ,"NEGATIVE_STRAND";
						$mapinfo = "$tg_chrom:$tg_end:$tg_start";
					}
					else {
						my $refsq = $sq->get_slice($tg_chrom, $tg_end, $tg_start);
						$mapinfo = "$tg_chrom:$tg_end:$refsq";
						if (rev_comp($refsq) eq $ref) {
							$dat->{$f_chr} = $tg_chrom;
							$dat->{$f_pos} = $tg_end;
							$dat->{$f_ref} = rev_comp($ref);
							$dat->{$f_alt} = rev_comp($alt);
							if (is_hgchr($tg_chrom)) {
								$outflag = 1;
								if ($tg_chrom ne $chr) {
									push @tags, "DIFF_CHROM" ,"NEGATIVE_STRAND";
								}
							}
							else {
								if ($ARGV{'--no-hgchr'}) {
									$outflag = 1;
								}
								push @tags, 'NOT_CANONCHAR', 'NEGATIVE_STRAND';
							}
						}
						elsif (rev_comp($refsq) eq $alt) {
							$dat->{$f_chr} = $tg_chrom;
							$dat->{$f_pos} = $tg_end;
							$dat->{$f_ref} = rev_comp($alt);
							$dat->{$f_alt} = rev_comp($ref);
							if (is_hgchr($tg_chrom)) {
								$outflag = 1;
								push @tags, 'SWITCH_REFALT', 'NEGATIVE_STRAND'
							}
							else {
								if ($ARGV{'--no-hgchr'}) {
									$outflag = 1;
								}
								push @tags, 'NOT_CANONCHR', 'REFALT_SWITCH', 'NEGATIVE_STRAND';
							}
						}
						else {
							push @tags, 'REFALT_MISMATCH', 'NEGATIVE_STRAND';
						}
					}
				}
			}
			else {
				if ($tg_strand eq '+') {
					$mapinfo = "$tg_chrom:$tg_start:$tg_end";
				}
				else {
					$mapinfo = "$tg_chrom:$tg_end:$tg_start";
				}
				push @tags, 'LENGTH_MISMATCH';
			}
		}
		else {
			$mapinfo = "$tg_chrom:$tg_strand:$tg_chrom2:$tg_strand2";
			push @tags, 'TWOEND_DFBLOCK';
		}
	}
	else {
		push @tags, 'UNMAPPED';
	}
	if ($outflag) {
		unless ($ARGV{'--no-norm'}) {
			if (length($dat->{$f_ref}) > 1 || length($dat->{$f_alt}) > 1) {
				my $var = { CHROM => $dat->{$f_chr}, POS => $dat->{$f_pos}, 
							REF => $dat->{$f_ref}, ALT => $dat->{$f_alt} };
				normalize($var, $sq);
				$dat->{$f_pos} = $var->{POS}; 
				$dat->{$f_ref} = $var->{REF};
				$dat->{$f_alt} = $var->{ALT};
			}
		}
		if ($ARGV{'--nochr'}) {
			$dat->{$f_chr} =~ s/^chr//; 
			$dat->{$f_chr} = 'MT' if $dat->{$f_chr} eq 'M'; 
		}
		print $fout join("\t", map { $_ // $ARGV{'--nastr'} } @{$dat}{@outfields}), "\n";
		$mapped ++;
	}
	if (@tags) {		
		print $flog $total, "\t", join(',', @tags), "\t", $varid, "\t", $mapinfo // "", "\n";
	}
}

print STDERR "Total ", commafy($total), " variants examined: ", commafy($mapped), " are lifted over\n";


__END__

=head1 NAME

liftover_vars -- Liftover genetic variants

=head1 USAGE

liftover_vars [options] -in INPUT -chain CHAIN -seq HGSEQ -out OUTPUT

=head1 DESCRIPTION

Don't use this script for VCF file. The intended use of the script is for biallelic 
annotated sequence variants.

NOTE: Variant is represented by the four mandatory fields: Chrom,Position,Ref,Alt.
Both Ref and Alt alleles should match regex /^[ACGTN]+$/. Deletions/Insertions
should follow the representation used in VCF file. 

=head1 REQUIRED ARGUMENTS

=over 

=item -[-]in[put] [=] <file>

Input data file for liftover.

=for Euclid:
  file.type: readable

=item -[-]field[s] [=] <string>

Field names for chromosome, position, reference allele and an optional alternative
allele. The fields should be comma separated. 

=item -[-]chain [=] <file>

Liftover chain file. When contig/chromsome name nomenclature does no match the one used in VCF, 
it still works for canonical chromosomes but results on non-chr contigs may not be accurate.

=for Euclid:
  file.type: readable

=item -[-]seq [=] <file>

Reference sequence for the target genome, 2bit or fasta format, should match the target (new)
assembly in the liftover chain file. This file is used to check ref seq in the new assembly. 
And We assume the original input has already been checked for ref seq in the orignal assembly.

=for Euclid:
  file.type: readable

=item -[-]out[put] [=] <file>

Output file name prefix. Two files prefix.txt and prefix.log will be in the output. 

The main output is the original input file with related fields subsituted.
A list of variants that cannot be mapped and those strand has been flipped can be found in the log file.

=for Euclid:
  file.type: writable

=back

=head1 OPTIONS

=over

=item -[-]sheet [=] <num>

If input is Excel workbook, specify the sheet number.

=item -[-]maxrow [=] <num>

=item -[-]maxcol [=] <num>

For spreadsheet: the max range of columns and rows.

=item -[-]skip [=] <nrow>

Number of rows to skip at the begniing.

=item -[-]fsep [=] <regex>

Field separator of the input file.

=item -[-]alias [=] <str>

Rename input file header.

=item -[-]no[-]chr

Strip off the 'chr' prefix in the output.

=item -[-]old[-]var [=] <field>

Provide the name for new columns to store original variant. 
It can be one column: "Chrom:Pos:Ref:Alt", or comma separated four columns (for chrom, pos, ref, alt).

=item -[-]no-norm

By default, the variant representation will be normalized. This switch will turn off normalization (not recommended).

=item -[-]no-hgchr

By default, only canonical chromosome will be in the output. Turning on this switch will output everything,
including those remapped to non-chr contigs.

=item -[-]nastr [=] <string>

Missing value string for empty fields.

=for Euclid:
	string.default: "."

=back

=cut


