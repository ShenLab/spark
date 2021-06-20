#!/usr/bin/env perl

use strict;
use warnings;
use Carp;
use IO::Prompt;
use File::Temp qw|tempdir|;
use List::MoreUtils qw|all none|;
use FaSlice;
use Utils::Seq qw|rev_comp|;
use Utils::Number qw|commafy|;
use Utils::List qw|parse_fields|;
use Utils::File::Iter qw|iter_file|;
use Genet::Var qw|normalize|;
use Genome::UCSC qw|hg_chrom is_hgchr|;
use Genome::UCSC::TwoBit;
use Genome::UCSC::Liftover qw|check_chain_chrpref|;
use Getopt::Euclid;


# Parse input file fields
my @fields = parse_fields($ARGV{'--fields'});
unless(@fields == 3 || @fields == 4) {
	croak "Must provide field names for chrom, start, end and optionally strand";
}
if (@fields == 3) {
	print STDERR "The following fields will be parsed for chrom, start, end :\n";
}
else {
	print STDERR "The following fields will be parsed for chrom, start, end, and strand :\n"; 
}
print STDERR join("\t", @fields), "\n";

# Iterator to table containing intervals
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
	{ skip => $ARGV{'--skip'}, sheet => $ARGV{'--sheet'}, maxcol => $ARGV{'--maxcol'}, maxrow => $ARGV{'--maxrow'}, %opt });

# Check field names
foreach my $field (@fields) {
	if (none { $field eq $_ } @$fnames) {
		print STDERR "Field $field cannot be found in the data file\n";
		exit 1;
	}
}
my ($f_chr, $f_st, $f_ed, $f_strand) = @fields;

# Check if old interval field exists
my @outfields = @$fnames;
my $oldintflag;
if (defined $ARGV{'--old-int'}) {
	if ($ARGV{'--old-int'} =~ /,/) {
		my @oldintfds = split(',', $ARGV{'--old-int'});
		unless(@oldintfds == 3 || @oldintfds == 4) {
			die "Must provide three or four fields for old-int!"
		}
		foreach my $fd (@oldintfds) {
			if (grep { $_ eq $fd } @$fnames) {
				print STDERR "$fd field already exists!";
				exit 1;
			}
		}
		$oldintflag = scalar(@oldintfds);
		unshift @outfields, @oldintfds;
	}
	else {
		unless (none { $_ eq $ARGV{'--old-int'} } @$fnames) {
			print STDERR "$ARGV{'--old-int'} field already exists!\n";
			exit 1;
		}
		$oldintflag = 1;
		unshift @outfields, $ARGV{'--old-int'};
	}
}

# Determine chr prefix in chain file and create lift chain obj
my ($inchr, $outchr) = check_chain_chrpref($ARGV{'--chain'});
my $lift = Genome::UCSC::Liftover->new($ARGV{'--chain'});

# Working directory
my $wrkdir;
if ($ARGV{'--wrkdir'}) {
	$wrkdir = $ARGV{'--wrkdir'};
	mkdir $wrkdir;
	croak "Cannot create working directory" unless -d $wrkdir;
}
else {
    $wrkdir = tempdir(CLEANUP => 1);
}

# Validate entries in input file and prepare input BED6 file 
my (%invalid, $tabchr, $ndigits);
{
	my $it2 = iter_file($ARGV{'--input'} eq '-' ? \*STDIN : $ARGV{'--input'},
		{ skip => $ARGV{'--skip'}, sheet => $ARGV{'--sheet'}, %opt });
	open my $fout, ">$wrkdir/input.bed" or die "Cannot write to $wrkdir/input.bed";

	my %known;
	while(my $dat = $it2->()) {
		my ($chrom, $start, $end, $strand) = @{$dat}{@fields};
		my $rangeid = @fields == 4 ? "$chrom:$start-$end($strand)" : "$chrom:$start-$end";

		# Validate data
		my $errflag;
		if (defined $strand) {
			unless ($strand eq '+' || $strand eq '-') {
				$invalid{$rangeid} = "Invalid strand: $strand";
				$errflag = 1;
			}
		}
		unless ($chrom =~ /^\w\S*$/ && ($start =~ /^[0-9][0-9,]*$/ || $ARGV{'--log10-scale'} && $start =~ /^[0-9][0-9,\.]*$/) &&
						 			   ($end =~ /^[0-9][0-9,]*$/   || $ARGV{'--log10-scale'} && $end =~ /^[0-9][0-9,\.]*$/)) {
			$invalid{$rangeid} = "Incorret position: $chrom:$start-$end";
			$errflag = 1;			 			 
		}
		$start =~ s/,//g if $start =~ /,/; $end =~ s/,//g if $end =~ /,/;
		
		unless(defined $tabchr) {
			$tabchr = $chrom =~ /^chr/ ? 1 : 0;
			unless($tabchr == $inchr) {
				warn "Chromosome nomenclatures in input and chain files are different";
			}
		}

		if ($ARGV{'--log10-scale'}) {
			if ($start =~ /\.([0-9]+)$/) {
				my $traildigit = $1;
				if(!defined $ndigits || length($traildigit) > $ndigits) {
					$ndigits = length($traildigit);
				}	
			}
			if ($f_ed ne $f_st && $end =~ /\.([0-9]+)$/) {
				my $traildigit = $1;
				if(!defined $ndigits || length($traildigit) > $ndigits) {
					$ndigits = length($traildigit);
				}
			}
		}

		# Prepare input
		if ($tabchr == 1 && $inchr == 0) {
			$chrom =~ /^chr/;; $chrom = 'MT' if $chrom eq 'M';
		}
		elsif ($tabchr == 0 && $inchr == 1) {
			$chrom = hg_chrom($chrom);
		}

		if ($ARGV{'--log10-scale'}) {
			$start = int($start * 10 ** $ARGV{'--log10-scale'});
			$start -= 1;
			if ($f_ed ne $f_st) {
				$end = int($end * 10 ** $ARGV{'--log10-scale'});
			}
		}
		else {
			unless ($ARGV{'--bed'}) {
				$start -= 1 if $start > 0;
			}
		}

		if (@fields == 3) {
			print $fout join("\t", $chrom, $start, $end, $rangeid), "\n";
		}
		else {
			print $fout join("\t", $chrom, $start, $end, $rangeid, 100, "+"), "\n";
		}
	}
}

# Run liftOver on input BED file
system(qq|liftOver -minMatch=$ARGV{'--fraction'} $wrkdir/input.bed $ARGV{'--chain'} $wrkdir/output.bed $wrkdir/output.unmapped|);

# Collect new genomic positions (+strand) from output BED, collect failure reason for unmapped intervals
my (%newrng, %unmapped);
{
	open my $fbed, "$wrkdir/output.bed" or die "Cannot open $wrkdir/output.bed";
	while(<$fbed>) {
		my @dat = split;
		if ($ARGV{'--log10-scale'} || !defined $ARGV{'--bed'}) {
			$dat[1] += 1;
		}
		if (@dat == 4) {
			$newrng{$dat[3]} = join("\t", @dat[0,1,2]);
		}
		elsif (@dat == 6) {
			$newrng{$dat[3]} = join("\t", @dat[0,1,2,5]);
		}
		else {
			die "Incorrect number of columns in output.bed : $_";
		}	 
	}
	open my $ferr, "$wrkdir/output.unmapped" or die "Cannot open $wrkdir/output.unmapped";
	while(<$ferr>) {
		my $reason;
		if (/^#/) {
			chomp;
			$reason = $_; $reason =~ s/^#//;
		}
		$_ = <$ferr>;
		unless($_) {
			die "Cannot find original entry that failed due to $reason";
		}
		my @dat = split;
		unless(@dat == 4 || @dat == 6) {
			die "Incorrect number of columns of original entry: $_";
		}
		$unmapped{$dat[3]} = $reason;
	}
}


# Reformat output 
if ($ARGV{'--output'} =~ /\.txt$/) {
	$ARGV{'--output'} =~ s/\.txt$//;
}
open my $fout, ">$ARGV{'--output'}.txt" or croak "Cannot write to $ARGV{'--output'}.txt";
if ($opt{header}) {
	print $fout join("\t", @outfields), "\n";
}
open my $flog, ">$ARGV{'--output'}.log" or croak "Cannot write to $ARGV{'--output'}.log";

my ($total, $remapped) = (0, 0);
while(my $dat = $it->()) {
	my ($chrom, $start, $end, $strand) = @{$dat}{@fields};
	my $rangeid = @fields == 4 ? "$chrom:$start-$end($strand)" : "$chrom:$start-$end";
	$total ++;

	if (defined $invalid{$rangeid}) {
		print $flog $rangeid, "\t", $invalid{$rangeid}, "\n";
		next;
	}
	elsif (defined $unmapped{$rangeid}) {
		print $flog $rangeid, "\t", $unmapped{$rangeid}, "\n";
		next;
	}

	#if (defined $ARGV{'--old-int'}) {
	#	$dat->{$ARGV{'--old-int'}} = $rangeid;
	#}
	if ($oldintflag) {
		if ($oldintflag == 1) {
			$dat->{$outfields[0]} = $rangeid;
		}
		elsif ($oldintflag >= 3) {
			$dat->{$outfields[0]} = $chrom;
			$dat->{$outfields[1]} = $start;
			$dat->{$outfields[2]} = $end;
			$dat->{$outfields[3]} = $strand // "." if $oldintflag == 4;
		}
	}

	unless(defined $newrng{$rangeid}) {
		die "Cannot find new location for a remapped range: $rangeid";
	}
	my @newloc = split(/\t/, $newrng{$rangeid});
	$remapped ++;

	# Updating existing genomic positions
	$dat->{$f_chr} = $newloc[0];
	if ($ARGV{'--hgchr'}) {
		next unless is_hgchr($dat->{$f_chr});
	}
	if ($ARGV{'--nochr'}) {
		$dat->{$f_chr} =~ s/^chr//; 
		$dat->{$f_chr} = 'MT' if $dat->{$f_chr} eq 'M'; 
	}

	if ($ARGV{'--log10-scale'}) {
		my $fmt = $ndigits ? "%.${ndigits}f" : "%d";
		$dat->{$f_st} = sprintf($fmt, $newloc[1]/(10**$ARGV{'--log10-scale'}));
		if ($f_ed ne $f_st) {
			$dat->{$f_ed} = sprintf($fmt, $newloc[2]/(10**$ARGV{'--log10-scale'}));
		}
	}
	else {
		$dat->{$f_st} = $newloc[1];
		if ($f_ed ne $f_st) {
			$dat->{$f_ed} = $newloc[2];
		}
	}	
	if (@fields == 4) {
		$dat->{$f_strand} = $newloc[3];
	}

	print $fout join("\t", map { $_ // $ARGV{'--nastr'} } @{$dat}{@outfields}), "\n";
}

print STDERR "Total ", commafy($total), " intervals examined: ", commafy($remapped), " are lifted over\n";


close $fout;
close $flog;


__END__

=head1 NAME

liftover_intervals -- Liftover genomic intervals

=head1 USAGE

liftover_intervals [options] -in INPUT -chain CHAIN -seq HGSEQ -out OUTPUT

=head1 DESCRIPTION

This is a wrapper of UCSC liftOver utility for lifting over genomic intervals.
We can use this utility to liftOver any genomic intervals with flexible format.

=head1 REQUIRED ARGUMENTS

=over 

=item -[-]in[put] [=] <file>

Input data file for liftover.

=for Euclid:
  file.type: readable

=item -[-]field[s] [=] <string>

Field names for chromosome, start, end and an optional strand. The fields should be comma separated. 
If end is not available, it can be the same as start.

=item -[-]chain [=] <file>

Liftover chain file. When contig/chromsome name nomenclature does no match the one used in VCF, 
it still works for canonical chromosomes but results on non-chr contigs may not be accurate.

=for Euclid:
  file.type: readable

=item -[-]out[put] [=] <file>

Output file name prefix. 

The main output is the original input file with related fields subsituted.
It will also output a list of intervals that cannot be mapped (.log).

=for Euclid:
  file.type: writable

=back

=head1 OPTIONS

=over

=item -[-]frac[tion] [=] <factor>

Minimal fraction of base pairs in the original interval that are preseved after liftover (default 0.95).

=for Euclid:
	factor.default: 0.95
	factor.type: number, factor >= 0 && factor <= 1

=item -[-]wrkdir [=] <dir>

Working directory.

=item -[-]sheet [=] <num>

If input is an Excel workbook, specify the sheet number.

=item -[-]maxrow [=] <num>

=item -[-]maxcol [=] <num>

For spreadsheet: the max range of columns and rows.

=item -[-]skip [=] <nrow>

Number of rows to skip at the begniing.

=item -[-]fsep [=] <regex>

Field separator of the input file.

=item -[-]log10[-scale] [=] <scale>

The log10 scale of original positions.
If original positions are in Mb, set this factor to 6; and 3 for Kb.
The Output will be in the same scale and have the same number of effective digits.

=for Euclid:
	scale.type: int > 0

=item -[-]bed

Indicate that input intervals are BED style (0-based half-open).

=item -[-]old[-int] [=] <field>

Provide the name for new column(s) to store original interval. 
It can be one column: "Chrom:Start-End_Strand", three columns or four columns.

=item -[-]no[-]chr

Strip off the 'chr' prefix in the output.

=item -[-]hgchr

By default, all remapped chromosomes and contigs will be in the output. 
Turning on this switch will only caonincal chromosomes.

=item -[-]nastr [=] <string>

Missing value string for empty fields.

=for Euclid:
	string.default: "."

=back

=cut


