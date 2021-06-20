#!/usr/bin/env perl

use strict;
use warnings;
use Data::Dumper;
use FindBin qw|$Bin|;
use Getopt::Euclid;
use Utils::Parser qw|sql_query|;
use List::MoreUtils qw|all any uniq|;
use Utils::File::Iter qw|iter_file|;

use lib "$Bin/../../lib";
use Shared qw|parse_fstr|;

# First determine which fields will be used as key, which will be uses as additional output
# We used the following rule:
# Field names that matche to the input and not specified by over-write will be used
# as keys for gene ID or gene name to fetch additional gene level information. The primary key
# will be the ID/name field appear at the first.
my @genefd; # this is the gene ID/name field(s) in the input
my @outfds; # this is the additional field(s) added to input (appear at the beginning of output)
my @keyfd; # this is the linked gene ID/name field(s) in Gxref
my @valfds; # this is the value field(s) to fetch gene level information from Gxref
my ($it, $fnames); # Iterator and field names of Gxref
my ($it2, $fnames2); # Iterator and field names of input 
my $fields; # parsed gxref fields
{
	# Parse gxref field specification and validate with gxref file
	$fields = parse_fstr($ARGV{'--gxref-fields'}, 1);
	my %opt;
	if (defined $ARGV{'--gxref-fsep'}) {
		$opt{fsep} = qr/$ARGV{'--gxref-fsep'}/;
	}
	else {
		unless($ARGV{'--gxref'} =~ /(csv|xlsx|ods|xls)$/) {
			$opt{fsep} = qr/\t/;
		}
	}
	if (all { /^\d+$/ } keys %$fields) {
		$opt{header} = 0;
	}

	($it, $fnames) = iter_file($ARGV{'--gxref'}, { skip => $ARGV{'--gxref-skip'}, %opt, 
		sheet => $ARGV{'--gxref-sheet'}, maxcol => $ARGV{'--gxref-maxcol'}, maxrow => $ARGV{'--gxref-maxrow'} });
	foreach my $infd (keys %$fields) {
		unless (grep { $infd eq $_ } @$fnames) {
			#print join("\n", @$fnames), "\n";
			die "Cannot find field $infd from gxref: $ARGV{'--gxref'}";
		}
	}

	($it2, $fnames2) = iter_file($ARGV{'--input'} eq '-' ? \*STDIN : $ARGV{'--input'}, { fsep => qr/\t/ });

	my %overwrite;
	if ($ARGV{'--over-write'}) {
		my @owfds = split(',', $ARGV{'--over-write'});
		foreach my $fname (@owfds) {
			$overwrite{$fname} = 1;
		}
	}
	# Determine various fields
	while(my ($field, $alias) = each %$fields) {
		if (grep { $alias eq $_ } @$fnames2) {
			if (defined $overwrite{$alias}) {
				warn "Field $alias in the input will be over-written!";
				delete $overwrite{$alias};
				push @valfds, $field;
			}
			else {
				push @genefd, $alias;
				push @keyfd, $field;	
			}
		}
		else {
			push @outfds, $alias;
			push @valfds, $field;
		}
	}
	if (%overwrite) {
		die "The following over-write fields can not found in input file: ", join(',', keys %overwrite);
	}
}

unless(@keyfd > 0 && @valfds > 0) {
	print STDERR "Empty key or value field(s)!\n";
	exit 1;
}

print STDERR "The following field(s) from Gxref will be used as ID or name for genes: ", join(',', @keyfd), "\n";
print STDERR "... matched to the field(s) from input: ", join(',', @genefd), "\n";
print STDERR "The following field(s) from Gxref will be used to fetch info: ", join(',', @valfds), "\n";
print STDERR "... and will appear in the output as: ", join(',', map { $fields->{$_} } @valfds), "\n";


# Gxref filter
my ($callback, $tokens);
if ($ARGV{'--gxref-filter'}) {
	($callback, $tokens) = sql_query($ARGV{'--gxref-filter'}, 1);
	foreach my $tok (@$tokens) {
		if ($tok->[0] eq 'FIELD') {
			unless (grep { $tok->[1] eq $_ } @$fnames) {
				die "Cannot find filtering field $tok->[1] from gxref file: $ARGV{'--gxref'}";
			}
		}
	}
}

# Fetch gene level information
# For each ID/name field, there will be a %ginfo hash stored in the array 
my @ginfo;
while(my $dat = $it->()) {
	if ($ARGV{'--gxref-filter'}) {
		next unless $callback->($dat);
	}
	for(my $ii = 0; $ii < @keyfd; $ii ++) {
		if (defined $ginfo[$ii]{$dat->{$keyfd[$ii]}}) {
			if ($ARGV{'--all'}) {
				foreach my $valfd (@valfds) {
					$ginfo[$ii]{$dat->{$keyfd[$ii]}}{$valfd} .= $ARGV{'--multi-sep'}.$dat->{$valfd};
				}
			}
			else {
				warn "Information for $dat->{$keyfd[$ii]} has already been found";
			}
		}
		else {
			foreach my $valfd (@valfds) {
				$ginfo[$ii]{$dat->{$keyfd[$ii]}}{$valfd} = $dat->{$valfd};
			}
		}
	}
}


my $fout;
if ($ARGV{'--output'}) {
	open $fout, ">$ARGV{'--output'}";
}
else {
	$fout = \*STDOUT;
}
print $fout join("\t", @outfds, @$fnames2), "\n";

while(my $dat = $it2->()) {
	# Expand gene ID/name fields
	my @geneids; # [G1,N1,C1],[G2,N2,C2]...
	foreach my $gid (split(';', $dat->{$genefd[0]})) {
		push @geneids, [$gid];
	} 
	# First field will be the primary ID, all others are secondary
	# they should have the same length as primary ID or 1 if collapsed
	for(my $ii = 1; $ii < @genefd; $ii ++) {
		my @gnames = split(';', $dat->{$genefd[$ii]});
		unless (@gnames == @geneids || @gnames == 1) {
			die "Incorrect secondary gene ID/name field $genefd[$ii]: $dat->{$genefd[$ii]}";
		}
		if (@gnames == 1 && @geneids > 1) {
			@gnames = ($gnames[0]) x scalar(@geneids);
		}
		for(my $jj = 0; $jj < @geneids; $jj ++) {
			push @{$geneids[$jj]}, $gnames[$jj];
		}
	}

	# Collect gene information for each gene starting from matching the primary ID
	my %outvals;
	foreach my $geneid (@geneids) {
		my $flag;
		for (my $ii = 0; $ii < @$geneid; $ii ++) {
			if (defined $ginfo[$ii]{$geneid->[$ii]}) {
				my $gdat = $ginfo[$ii]{$geneid->[$ii]};
				foreach my $valfd (@valfds) {
					my $alias = $fields->{$valfd};
					push @{$outvals{$alias}}, $gdat->{$valfd};
				}
				$flag = 1;
				last;
			}
		}
		unless($flag) {
			foreach my $valfd (@valfds) {
				my $alias = $fields->{$valfd};
				push @{$outvals{$alias}}, $ARGV{'--nastr'};
			}
		}
	}

	# Update the current data structure in-situ
	foreach my $alias (map { $fields->{$_} } @valfds) {
		my @uqvals = uniq sort @{$outvals{$alias}};
		if (@uqvals == 1) {
			$dat->{$alias} = $uqvals[0];
		}
		else {
			$dat->{$alias} = join(";", @{$outvals{$alias}});
		}
	}

	print $fout join("\t", @{$dat}{@outfds}, @{$dat}{@$fnames2}), "\n";
}



__END__

=head1 NAME

anno_genes.pl -- Adding additional gene level information.

=head1 NOTES

The default way to add gene information used by standard in anno_seqvar pipeline is to lookup information
by GeneID from an external table file. However, GeneID is not always available in other sources, and we may
want to add additional gene level information after we have done the variant annotation. It is also often
the case that matching by a single gene ID or name will miss a few genes when IDs have changed. 

This script is written to append gene level information later to an annotated variant table file. It supports 
looking up genes by any ID or name as long as they can be matched to the query file. Miltiple ID/name matching
is also supported: the first ID will be the primary ID, followed by secondary IDs. Priority will be given
to matching with ID with higher precedence. 

By default, the overlapping fields between gxref and input are considered as gene ID or name fields, they
will be used for query gene level information with the priority given to the matched one that appear first.
And fields not appear in the input will be used as new value fields. There is an option to update existing 
field in the input, they will be removed from overlapping fields used for gene matching.

When multiple genes exist in a row (separated by ";"), information for each gene will be collected.
When multiple entries for a gene exist, there is an option to output then all. The additional columns will 
appear at the beginning of the input table. 

The format of external gene info can be text, csv or excel files with felixble format. The limitation is you can 
only specify one lookup table file at a time. 

=head1 REQUIRED ARGUMENTS

=over

=item -[-]in[put] [=] <table>

Input table. It should be tab separated with gene ID and/or name field.
Use '-' if input is piped from STDIN.

=item -[-]gxref [=] <file>

Lookup table of gene level extra information.

=for Euclid:
	file.type: readable

=item -[-]gxref-fields [=] <string>

Fields in the lookup table used for matching and output. The first one or more fields should always 
be gene ID or name and should be renamed to match the same column name as input. All remaining
fields will appear as additional columns (appearing at the first) in the output table.
The new column fields should not exist in the input unless specified by --over-write.

An example: EnsemblID:GeneID,GeneName:Symbol,AutismDscore:Dscore,EnsemblScore:Dscore_Ens

When multiple entries corresponding to the same gene exist in the lookup table, the default is
to used the same convention as in the annotation pipeline, only the *first one* will be used unless
--all is turned on.

=back

=head1 OPTIONS

=over

=item -[-]over[-]write [=] <fields>

Specify the fields from the input to be over-written. 

=item -[-]out[put] [=] <table>

Output table file name. If not provided, will dump to STDOUT. Fields will be tab separated.

=item -[-]all

Slurp all information for each gene.

=item -[-]multi-sep [=] <sepchar>

The separator for multuple records of the same gene. Only used when --all option is turned on.

=for Euclid:
	sepchar.default: "|"

=item -[-]nastr [=] <string>

Output value for missing data, default ".".

=for Euclid:
	string.default: "."

=item -[-]gxref-filter [=] <string>

Filter the gxref table before slurping the data.

=item -[-]gxref-fsep [=] <string>

Gxref file: field separator. It defailts to tab for text file; "," for csv file.

=item -[-]gxref-skip [=] <number>

Gxref file: number of rows to skip at the begnning.

=item -[-]gxref-sheet [=] <number>

Gxref file: if we lookup gene information from an Excel workbook, specify the sheet name or number.

=item -[-]gxref-maxcol [=] <number>

=item -[-]gxref-maxrow [=] <number>

Gxref file: if we lookup gene information from an Excel workbook, the max range of columns and rows.

=back



