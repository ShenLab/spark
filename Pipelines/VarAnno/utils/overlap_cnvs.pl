#!/usr/bin/env perl

use strict;
use warnings;
use Data::Dumper;
use File::Basename;
use FindBin qw|$Bin|;
use Sort::Maker;
use List::MoreUtils qw|uniq all|;
use Getopt::Euclid;
use Utils::List qw|parse_fields|;
use Utils::Hash qw|peek_hash|;
use Utils::File::Iter qw|iter_file|;
use Genome::UCSC qw|hg_chrom|;
use Genome::Ranges qw|validate_elem|;
use Genome::UCSC::BinKeeper;


use lib "$Bin/../../lib";
use Shared qw|parse_fstr parse_tabfile read_probes find_proberng region_bk|;
use Variants qw|cnv_type parse_operations|;

# Validate field specifications for input and dbfile
my $dbfields = parse_fstr($ARGV{'--db-fields'}, 1);
my $matchtype;
{
	# Note: Orders of the fields are important for input
	my @infields = parse_fields($ARGV{'--in-fields'});
	unless(@infields == 4) {
		die "Incorrect number of input file fields: $ARGV{'--in-fields'}";
	}
	foreach my $stdfield (qw|Chrom Start End|) {
		unless(grep { $_ eq $stdfield } values %$dbfields) {
			die "Cannot find standard field $stdfield from dbfile";
		}
	}
	if (grep { $_ eq 'Type' } values %$dbfields) {
		$matchtype = 1;
	}
}

# Parse and create operators
my $ops = parse_operations($ARGV{'--operation'}, 
			{ orsep => $ARGV{'--or-sep'}, dfsep => $ARGV{'--df-sep'}, dbfields => $dbfields });

# Check additional fields are absent from input file
my ($it, $fnames, $keyfields) = parse_tabfile($ARGV{'--input'}, $ARGV{'--in-fields'}, 4, 4);
my @sumfields = @$fnames;
foreach my $opname (keys %$ops) {
	if (grep { $_ eq $opname } @$fnames) {
		unless($ARGV{'--over-write'}) {
			die "The operation output column $opname already exists in the input";
		}
	}
	else {
		push @sumfields, $opname;
	}
}

# Indicate chromsome name prefix
my ($dbchr, $cnvchr, $probechr);

# If provided, probes files will be used as backbone coordinates of all intervals
my $probes;
if ($ARGV{'--probes'}) {
	$probes = read_probes($ARGV{'--probes'}, $ARGV{'--probe-fields'});
	my $firstpbchr = peek_hash($probes);
	if ($firstpbchr =~ /^chr/) {
		$probechr = 1;
	}
	else {
		$probechr = 0;
	}
}

# Parse database file name
my ($dbfile, $dboutfd) = split(":", $ARGV{'--dbfile'});
unless(-f $dbfile) {
	die "Cannot find dbfile: $dbfile";
}
unless(defined $dboutfd) {
	$dboutfd = (split(q|\.|, basename($dbfile)))[0];
}

# Check detail fields are absent from input file
my @detfields = @$fnames;
my @dbxfds = values %$dbfields;
my @calcfds = qw|iLen T Q TxQ|;
foreach my $field (@dbxfds, @calcfds) {
	if (grep { $dboutfd."_".$field eq $_ } @$fnames) {
		die "Field ${dboutfd}_$field already exist in the input file";
	}
	push @detfields, $dboutfd."_".$field;
}


# Prepare database binkeeper for querying
my $bk = region_bk($dbfile, $dbfields, 
					{ probes => $probes, poverlap => $ARGV{'--probe-overlap'},
					  callback => sub {
					  	 my $dat = shift @_;
					  	 if (defined $dat->{Type}) {
					  	 	$dat->{Type} = cnv_type($dat->{Type}, { strict => 1 });
					  	 }
					  	 return $dat;
					  }	 });
my $firstdbchr = peek_hash($bk);
if ($firstdbchr =~ /^chr/) {
	$dbchr = 1;
}
else {
	$dbchr = 0;
}

# Prepare output files
my ($fsum, $fdet);
if ($ARGV{'--output'}) {
	my $outprefix = $ARGV{'--output'};
	if ($outprefix =~ /\.txt$/) {
		$outprefix =~ s/\.txt$//;	
	}
	open $fsum, ">$outprefix.txt" or die "Cannot write to $outprefix.txt";
	print $fsum join("\t", @sumfields), "\n";
	open $fdet, ">$outprefix.details.txt" or die "Cannot write to $outprefix.details.txt";
	print $fdet join("\t", @detfields), "\n";
}
else {
	$fsum = \*STDOUT;
	print $fsum join("\t", @sumfields), "\n";
}

# Sorter will determine the order of overlapping intervals
my $sorter = sort_maker($ARGV{"--sort"});

# Testing overlap and write to output
while(my $dat = $it->()) {
	my ($chrom, $start, $end, $type) = @{$dat}{@$keyfields};
	validate_elem($chrom, $start, $end);

	unless(defined $cnvchr) {
		if ($chrom =~ /^chr/) {
			$cnvchr = 1;
		}
		else {
			$cnvchr = 0;
		}
		unless($dbchr == $cnvchr) {
			warn "Chromosome nomenclatures in variant table and database are different";
		}
		if ($ARGV{'--probes'}) {
			unless($probechr == $dbchr) {
				warn "Chromosome nomenclatures in probe file and database are different";
			}
		}
	}

	# Make chr name consistent between tab and db
	if ($cnvchr == 1 && $dbchr == 0) {
		$chrom =~ s/^chr//; 
	}
	elsif ($cnvchr == 0 && $dbchr == 1) {
		$chrom = hg_chrom($chrom);
	}

	my @overlaps;
	if ($ARGV{'--probes'}) {
		# First, for each CNV region find overlapping probe ranges from the probe list
		my $chrom_probe = $chrom;
		if ($dbchr == 1 && $probechr == 0) {
			$chrom_probe =~ s/^chr//; 
		}
		elsif ($dbchr == 0 && $probechr == 1) {
			$chrom_probe = hg_chrom($chrom_probe);
		}
		my ($lo_pos, $hi_pos) = find_proberng($probes, $chrom_probe, $start, $end, 
											{ overlap => $ARGV{'--probe-overlap'} });
		# Second, find overlapping probes ranges of database entries
		if (defined $lo_pos && $hi_pos) {
			@overlaps = $bk->find_range($chrom, $lo_pos, $hi_pos, 
				{ tCover => $ARGV{'--targ-overlap'}, qCover => $ARGV{'--query-overlap'}, calc => 1 });
		}
	}
	else {
		@overlaps = $bk->find_range($chrom, $start, $end, 
			{ tCover => $ARGV{'--targ-overlap'}, qCover => $ARGV{'--query-overlap'}, calc => 1 });
	}

	# If type is specified, we need to match the CNV type 
	if ($matchtype) {
		my $cnv_type = cnv_type($type, { strict => 1 });
		@overlaps = $sorter->(grep { $_->[2]{Type} eq $cnv_type } @overlaps);
	}
	else {
		@overlaps = $sorter->(@overlaps);
	}

	# Write each overlapping regions to details file
	if ($fdet) {
		foreach my $overlap (@overlaps) {
			print $fdet join("\t", @{$dat}{@$fnames}, @{$overlap->[2]}{@dbxfds}, @{$overlap->[3]}{@calcfds}), "\n";
		}
	}
	
	# Collapse the results to summary file
	if (@overlaps) {
		while(my ($outfield, $callback) = each %$ops) {
			$dat->{$outfield} = $callback->(map { $_->[2] } @overlaps);
			if (!defined $dat->{$outfield} || $dat->{$outfield} eq "") {
				$dat->{$outfield} = $ARGV{'--nastr'};
			}
		}
	}
	print $fsum join("\t", map { $_ // $ARGV{'--nastr'} } @{$dat}{@sumfields}), "\n";
}


# Make a custom sort subroutine from specification
# Example: TxQ(number,descending);Start(number,ascending)
#        => (number => { code => '$_->[3]{TxQ}', descending => 1 })
sub sort_maker {
	my ($fstr) = @_;
	my @params;
	foreach my $conf (split(";", $fstr)) {
		my ($field, $type, $order) = ($conf =~ /^(\w+)\((\w+),(\w+)\)$/);
		unless(defined $field && defined $type && defined $order) {
			die "Cannot parse field, type and order from $conf";
		}
		unless($type eq "number" || $type eq "string") {
			die "Does not support field type: $type";
		}
		unless($order eq "ascending" || $order eq "descending") {
			die "Does not support sort order: $order";
		}
		if (grep { $field eq $_ } values %$dbfields) {
			push @params, $type, { code => '$_->[2]{'.$field.'}', $order => 1 };
		}
		else {
			if (grep { $field eq $_ } qw|iLen T Q TxQ|) {
				push @params, $type, { code => '$_->[3]{'.$field.'}', $order => 1 };
			}
			else {
				die "Cannot recognize sorting field: $field";
			}
		}
	}
	my $sorter = make_sorter(plain => 1, @params);
	return $sorter;
}

__END__

=head1 NAME

overlap_cnvs.pl -- Find ovelapping features for CNVs.

=head1 NOTES

To find and summarize features overlapping with CNVs, we need following two extension to bedtools
operation. 
1. We need to take into account the CNV type or alleles when looking for overlapping CNVs.
2. To account for imprecision of CNV boundaries, we may use target regions (array probes)
to define overlapping proportions.

This script basically is a combined action of "intersect" and "groupby" modules of bedtools.
We will first generate a defailed overlapping table. It contains columns from the original input
plus the specified fields in the dbfile (Chrom,Start,End,Type,Label,...), followed by calculated
fields for overlaps (T,Q,iLen,TxQ). The additional fields will have column name prefix specified
as part of dbfile (See below).

Then one or more group operations can be specified to create summary output. We selected operations
from bedtools' groupby module that are relevant for CNVs. Valid operations are:

	* count - number of overlapping intervals 
	* count_distinct - distinct number of overlapping intervals based on Label column
		<-- count and count_distinct are useful for determining CNV carrier frequencies
	* collapse - print a comma separated list of Labels
	* first - print the first Label of all overlapping intervals
	* distinct - print a comma separated list of distinct Labels
	* freq - print a comma separated list of Labels and the number of times they appear

If some additional operations not defined by this script. We can apply bedtools groupby directly
on defailed overlappig output.

We also provide ways to sort overlapping intervals in the detailed output for each input CNV.
The default is to sort overlapping intervals by TxQ from largest to smallest. Some summary operation
may depends on the sorting.


=head1 REQUIRED ARGUMENTS

=over

=item -[-]in[put] [=] <file>

Input CNV table. It must be tab separated with required fields.

=for Euclid:
	file.type: readable

=item -[-]db[file] [=] <dbfile>

The database file of intervals-based genomic features. It must be TAB separated with required fields 
specified below. Intervals from database will be stored in BinKeeper for range-based query. 

Optionally, we can specify (after the file name) the column name prefix that appear in the detailed 
overlapping file. Default will be dbfile name prefix are removing suffix. Note: Column names from dbfile
cannot contain special characters.

=item -[-]db[-]fields [=] <fstr>

Specify fields from DB file. Fields that specifies a region or CNV should be renamed to standard names: 
"Chrom,Start,End,Type,Label". The first three are required, Type is optional. If Type is provided, 
we will only look for genomic regions that match the Type with CNVs. Other fields are optional, and
will appear as they are in the detailed overlapping output and used by operation.

=item -[-]op[eration] [=] <fstr>

One or more operations and their alias as output column name. 
Example: count_distinct(FamID,IID):CarrierFreq;distinct(Syndrome):KnownPathoCNVs

Valid operations are: count,count_distinct,collapse,first,distinct,freq
These operations were taken from bedtools groupby utility. 
See bedtools document for details: https://bedtools.readthedocs.io/en/latest/content/tools/groupby.html

=back

=head1 OPTIONS

=over

=item -[-]out[put] [=] <file>

Output file name prefix. It not provided, will write to STDOUT; 
if provided, also write out a details table (.details.txt) will include all overlapping intervals.

=item -[-]in[-]fields [=] <fstr>

The input file field names for Chrom,Start,End,Type. Must be in the specified order without alias.
Default are taken from the output of anno_cnvs pipeline.

=for Euclid:
	fstr.default: "Chrom,Start,End,Type"

=item -[-]over[-]write

Ovwr-write existing columns of the same name in input.

=item -[-]or-sep [=] <string>

Separator for joining different record labels in the output field. 
It will be used when the specified operation is to list multiple items.

=for Euclid:
	string.default: ','

=item -[-]df[-]sep [=] <char>

The db field separator, used to join multiple db fields.

=for Euclid:
	char.default: "|"

=item -[-]query-over[lap] [=] <perc>

Minimal overlap required as fraction of CNV, default: 0 (any overlap).

=for Euclid:
	perc.default: 1e-10

=item -[-]targ-over[lap] [=] <perc>

Minimal overlap required as fraction of database entry, default 0 (any overlap).

=for Euclid:
	perc.default: 1e-10

=item -[-]sort [=] <string>

Sort the overlapping intervals by the fields spcified in string. It should be a field name,
followed by field type (number or string) and sort order (ascending or descending).

=for Euclid:
	string.default: "TxQ(number,descending);Start(number,ascending)"

=item -[-]probes [=] <probfile>

Genomic positions for SNP array probes or exome targets.
When probe file is specified, overlap length and proportion will be calculated based
on the number of overlapping probes. Database regions with no probe coverage will
not be used for overlap testing. Input CNVs with no probe coverage will find no overlaps.

=item -[-]probe-fields [=] <fstr>

Specify columns of chrom,start and end in the probes file. End can be optional.
Default will be inferred based on file name suffix. For example, "1:Chrom,2:Start,3:End"
for BED file, "Chr:Chrom,Position:Start" for PFB file. It must be provided for other
file types.

=item -[-]probe-over[lap] [=] <perc>

If probe file is provided, this is the percentage of a probe that are covered by the
query interval for it to be included as a backbone for that query range, default; 0.5.

=for Euclid:
	perc.default: 0.5

=item -[-]nastr [=] <string>

Output value for missing data, default ".".

=for Euclid:
	string.default: "."

=back



