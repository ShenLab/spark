use strict;
use warnings;
use Getopt::Euclid;
use List::MoreUtils qw|all none|;
use List::BinarySearch qw|binsearch_pos|;
use Utils::List qw|insert_after parse_fields|;
use Utils::File::Iter qw|iter_file|;
use Genome::UCSC qw|hg_chrom is_hgchr|;

# Parse input file fields
my @fields = parse_fields($ARGV{'--fields'});
unless(@fields == 2) {
	die "Must provide field names for chrom and position";
}
print STDERR "The following fields will be parsed for Chrom,Position: ", join(",", @fields), "\n";

# Prepare input table iterator
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
	{ skip => $ARGV{'--skip'}, sheet => $ARGV{'--sheet'}, 
	  maxcol => $ARGV{'--maxcol'}, maxrow => $ARGV{'--maxrow'}, %opt });

foreach my $field (@fields) {
	if (none { $field eq $_ } @$fnames) {
		print STDERR "Field $field cannot be found in the data file\n";
		exit 1;
	}
}

# Iterator over variant file and find out the genetic map positions
my @outfields = @$fnames;
if (grep { $_ eq $ARGV{'--field-add'} } @$fnames) {
	unless($ARGV{'--over-write'}) {
		die "The new field $ARGV{'--field-add'} already exist in input file!";
	}
	else {
		warn "The original field $ARGV{'--field-add'} will be over-written";
	}
}
else {
	insert_after(\@outfields, $fields[1], $ARGV{'--field-add'});
}

# Lookup genetic positions for physical locations
my (%pmap, %gmap, $m_chr);
{
	open my $fmap, $ARGV{'--map'} or die "Cannot open $ARGV{'--map'}";
	my ($prev_chr, $prev_pos, $prev_gpos);
	while(<$fmap>) {
		my ($chr, $pos, $gpos) = split;
		unless(defined $m_chr) {
			if ($chr =~ /^chr/) {
				$m_chr = 1;
			}
			else {
				$m_chr = 0;
			}
		}
		if (defined $prev_chr && $prev_chr eq $chr) {
			unless($pos >= $prev_pos && $gpos >= $prev_gpos) {
				die "Genet map is not position sorted in chromsome $chr";
			}
		}
		
		push @{$pmap{$chr}}, $pos;
		push @{$gmap{$chr}}, $gpos;

		($prev_chr, $prev_pos, $prev_gpos) = ($chr, $pos, $gpos);
	}
}


my $fout;
if ($ARGV{'--output'}) {
	open $fout, ">$ARGV{'--output'}" or die "Cannot write to $ARGV{'--output'}";
}
else {
	$fout = \*STDOUT;
}
print $fout join("\t", @outfields), "\n" if $opt{header};

my $t_chr;
while(my $dat = $it->()) {
	my ($chr, $pos) = @{$dat}{@fields};
	unless(defined $t_chr) {
		if ($chr =~ /^chr/) {
			$t_chr = 1;
		}
		else {
			$t_chr = 0;
		}
		if ($t_chr != $m_chr) {
			warn "Chromosome nomenclautures in variant table and recomb map file are different!";
		}
	}
	if ($t_chr == 0 && $m_chr == 1) {
		$chr = hg_chrom($chr);
	}
	elsif ($t_chr == 1 && $m_chr == 0) {
		$chr =~ s/^chr//;
	}
	if (defined $pmap{$chr}) {
		my $index = binsearch_pos {$a <=> $b} $pos, @{$pmap{$chr}};
		if ($index < 0 || $index > @{$pmap{$chr}}) {
			die "Index out of bound: $index"
		}
		elsif ($index == 0) {
			$dat->{$ARGV{'--field-add'}} = $gmap{$chr}[$index];
		}
		elsif ($index == @{$pmap{$chr}}) {
			$dat->{$ARGV{'--field-add'}} = $gmap{$chr}[$index-1];
		}
		else {
			unless($pos >= $pmap{$chr}[$index-1] && $pos <= $pmap{$chr}[$index]) {
				die "Position $pos is out of range of bounding interval: ".
					$pmap{$chr}[$index-1]."~".$pmap{$chr}[$index];
			}
			# linear interp 
			my $factor = ($pos-$pmap{$chr}[$index-1])/($pmap{$chr}[$index]-$pmap{$chr}[$index-1]);
			$dat->{$ARGV{'--field-add'}} = $gmap{$chr}[$index-1] +
					$factor * ($gmap{$chr}[$index]-$gmap{$chr}[$index-1]);
		}
	}
	else {
		$dat->{$ARGV{'--field-add'}} = '.';
	}
	print $fout join("\t", @{$dat}{@outfields}), "\n";
}


__END__

=head1 NAME

lookup_genetpos -- looking up genetic positions for physical locations

=head1 REQUIRED ARGUMENTS

=over 

=item -[-]in[put] [=] <file>

Input data table.

=for Euclid:
  file.type: readable

=item -[-]field[s] [=] <string>

Field names for chromosome, position. The fields should be comma separated. 
For file without header, it can be a list of column numbers.

=item -[-]map [=] <file>

The recombination map file. Format: Chr, Position, Genet position.

=for Euclid:
  file.type: readable

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

=item -[-]field-add [=] <field>

The new field name to be added to the output. Default: "GenetPos".

=for Euclid:
	field.default: "GenetPos"

=item -[-]over-write

Over-write an existing column if additional field name already exist in input.

=item -[-]out[put] [=] <file>

The output file name.

=item -[-]scale [=] <factor>

Scaling factor to transform from cM to M or from M to cM.

=item -[-]nastr [=] <string>

Missing value string for empty fields.

=for Euclid:
	string.default: "."

=back

=cut
