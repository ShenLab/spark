use strict;
use warnings;
use Perl6::Slurp;
use FindBin qw|$Bin|;
use File::Basename;
use List::MoreUtils qw|uniq all|;
use Getopt::Euclid;
use Utils::List qw|insert_after|;
use Utils::Parser qw|sql_query|;
use Utils::File::Iter qw|iter_file|;
use Genome::UCSC qw|hg_chrom|;
use Genome::UCSC::GeneTable qw|iter_genePred iter_geneTab gdat_is_coding gdat_list_coding_exons|;
use Genome::UCSC::BinKeeper;

use lib "$Bin/../../lib";
use Shared qw|parse_tabfile parse_fstr|;
use Variants qw|cnv_type|;

my $fields = parse_fstr($ARGV{'--gxref-fields'}, 1);

# Reading gene level information from gene lookup table
my (%genes, $field_flag, $new_column);
# %genes will store gene information
# field_Flag will be 1 if additional fields from gxref are available
# new_column will be the name of new column added to the existing columns in input
{	
	my @lfields = keys %$fields;
	my $idcol = shift @lfields; # this is the linke gene ID/name field in Gxref
	$field_flag = 1 if @lfields > 0;
	
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

	my ($gxref, $suffix) = split(':', $ARGV{'--gxref'});
	if (!defined $suffix) {
		$new_column = (split(q|\.|, basename($gxref)))[0];
	}
	else {
		$new_column = $suffix;
	}

	my ($it, $fnames) = iter_file($gxref, { skip => $ARGV{'--gxref-skip'}, %opt, 
		sheet => $ARGV{'--gxref-sheet'}, maxcol => $ARGV{'--gxref-maxcol'}, maxrow => $ARGV{'--gxref-maxrow'} });
	foreach my $infd (keys %$fields) {
		unless (grep { $infd eq $_ } @$fnames) {
			die "Cannot find field $infd from gxref: $ARGV{'--gxref'}";
		}
	}

	my ($callback, $tokens);
	if ($ARGV{'--gxref-filter'}) {
		($callback, $tokens) = sql_query($ARGV{'--gxref-filter'}, 1);
		foreach my $tok (@$tokens) {
			if ($tok->[0] eq 'FIELD') {
				unless (grep { $tok->[1] eq $_ } @$fnames) {
					die "Cannot find filtering field $tok->[1] from gxref: $ARGV{'--gxref'}";
				}
			}
		}
	}

	while(my $dat = $it->()) {
		if ($ARGV{'--gxref-filter'}) {
			next unless $callback->($dat);
		}
		if (@lfields > 0) {
			if (defined $genes{$dat->{$idcol}}) {
				warn "Information for $dat->{$idcol} has already been found";
				next;
			}
			$genes{$dat->{$idcol}} = join($ARGV{'--gf-sep'}, @{$dat}{@lfields});	
		}
		else {
			$genes{$dat->{$idcol}} = 1;
		}
	}
}

# Gene and CDS binkeeper if we want to find overlapping genes defined by UCSC genome browser track
my ($genebk, $cdsbk);
if (defined $ARGV{'--genetab'}) {
	$genebk = Genome::UCSC::BinKeeper->new();
	$cdsbk = Genome::UCSC::BinKeeper->new();

	my $iter;
	if (-f $ARGV{'--genetab'}) {
		$iter = iter_genePred($ARGV{'--genetab'});
	}
	else {
		unless($ARGV{'--genetab'} =~ /^\w+\.\w+$/) {
			die "Incorrect specification of genetab: $ARGV{'--genetab'}";
		}
		$iter = iter_geneTab(split('\.', $ARGV{'--genetab'}));
	}

	# Put all CDS exons to binkeeper, also put a CDS start/end to gene binkeeper
	while(my $dat = $iter->()) {
		next unless gdat_is_coding($dat);
		$genebk->add($dat->{chrom}, $dat->{cdsStart}, $dat->{cdsEnd}, $dat->{name2});
		foreach my $cdexon (gdat_list_coding_exons($dat)) {
			$cdsbk->add($dat->{chrom}, $cdexon->[0], $cdexon->[1], $dat->{name2});
		}
	}
}

# Prepare output
my ($it, $fnames) = iter_file($ARGV{'--input'}, { fsep => qr/\t/ });
my $f_gene = (values %$fields)[0];
unless(grep { $_ eq $f_gene } @$fnames) {
	unless ($ARGV{'--genetab'}) {
		die "Cannot find gene ID/name field $f_gene from input file";
	}
}
my @outfields = @$fnames;
if (grep { $_ eq $new_column } @$fnames) {
	unless($ARGV{'--over-write'}) {
		die "New field $new_column was already defined in the input";
	}
	else {
		warn "The field $new_column from input file will be over-written";
	}
}
else {
	# The new column of highlighted gene will appear after full gene list 
	if (grep { $_ eq $f_gene } @$fnames) {
		insert_after(\@outfields, $f_gene, $new_column);
	}
	else {
		push @outfields, $new_column;
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

# Now iterate through input file to highlight genes
if (defined $ARGV{'--genetab'}) {
	my @keyfields = split(',', $ARGV{'--cnv-fields'});
	unless(@keyfields == 4) {
		die "Incorrect number of CNV fields in input file: $ARGV{'--cnv-fields'}";
	}
	foreach my $field (@keyfields) {
		unless(grep { $_ eq $field } @$fnames) {
			die "Cannot find CNV field $field from input file: $ARGV{'--cnv-fields'}";
		}
	}
	while(my $dat = $it->()) {
		my ($chr, $st, $ed, $type) = @{$dat}{@keyfields};
		my $cntype = cnv_type($type);
		$chr = hg_chrom($chr) if $chr !~ /^chr/;
		my @cngenes = grep { defined $genes{$_} } uniq sort map { $_->[2] }
			$genebk->find_range($chr, $st, $ed, { tCover => $ARGV{'--gene-overlap'} });
		my @disruptexons = grep { defined $genes{$_} } uniq sort map { $_->[2] }
			$cdsbk->find_range($chr, $st, $ed, { tCover => $ARGV{'--exon-overlap'} });
		if (@cngenes > 0 || @disruptexons > 0) {
			if ($cntype eq 'Del') {
				if ($field_flag) {
					$dat->{$new_column} = 
					join($ARGV{'--or-sep'}, map { "$_($genes{$_})" } @disruptexons);
				}
				else {
					$dat->{$new_column} = 
					join($ARGV{'--or-sep'}, @disruptexons);
				}
			}
			else {
				my %known = map { $_ => 1 } @cngenes;
				if ($field_flag) {
					$dat->{$new_column} = 
					join($ARGV{'--or-sep'}, 
						map { $known{$_} ? "$_($genes{$_})" : "${_}*($genes{$_})" } @disruptexons);
				}
				else {
					$dat->{$new_column} = 
					join($ARGV{'--or-sep'}, map { $known{$_} ? $_ : "${_}*" } @disruptexons);
				}
			}
		}
		print $fout join("\t", map { $dat->{$_} // "." } @outfields), "\n";
	}
}
else {
	while(my $dat = $it->()) {
		if ($dat->{$f_gene} ne '.') {
			my @highlights = grep { defined $genes{$_} } split(";", $dat->{$f_gene});
			if (@highlights) {
				if ($field_flag) {
					$dat->{$new_column} = join($ARGV{'--or-sep'}, map { "$_($genes{$_})" } @highlights);
				}
				else {
					$dat->{$new_column} = join($ARGV{'--or-sep'}, map { "$_" } @highlights);
				}		
			}
		}
		print $fout join("\t", map { $dat->{$_} // "." } @outfields), "\n";
	}
}




__END__

=head1 NAME

highlight_cnvgenes.pl -- highlight genes affected by CNVs.

=head1 NOTES

The default way to higlight genes in standard anno_cnv pipeline is to define gene sets
and higlight genes within these sets.

This script supplement the existing pipeline and also provide a few extensions:
1. We can lookup gene sets by gene names or other non-ensembl IDs.
2. We can select additional fields from gene table to show together with highlighted genes. 
3. We can also query the CNV overlapping genes from UCSC database table instead of
using genes from existing annotations. 
4. When we query genes from database table by ourselves, if highlighted genes is not 
entirely duplicated, it will have an asterisk after gene name. E.g. NRXN1*(SFARI1or2)

In this script, we only support looking up genes and relevant fields from a gene table.
This is similar to define gene set from two columns of gene table and form a superset 
of all gene sets fromt this table. But we now support multiple columns and not limited to
categorical variables. Gene table can be text, csv, or excel format, and can be filtered 
by SQL expression. 

=back

=head1 REQUIRED ARGUMENTS

=over

=item -[-]in[put] [=] <table>

Input CNV table file to be annotated. It should be tab separated with required fields.

=item -[-]gxref [=] <file>

Lookup table to define gene sets to be highlighted. Optionally, we can specify the column name
of highlighted genes in the output after ":". Default will be the basename of gxref file after 
removing suffix. E.g. /Path/To/SFARI/Gene/Table:SFARIGenes. The highlighted genes will appear
after gene IDs or names that are used to query the gxref file.

=item -[-]gxref-fields [=] <string>

Provide gene name and other columns for the gxref file. 

The first field must be gene ID or name used to look up the gene set. If we rely on previously annotated 
genes for CNVs in the input, it should be renamed to match the field name in input. For example: 
GeneName:Symbols,pLI,lof_z matches GeneName in gxref file to Symbols in input.

If we query genes by ourselves from UCSC genome browser, make sure that gene names or IDs used by the 
databasae (name2 field) match the ones used in gxref table.

The values of remaining fields will be concatenated and appended to the genes higlighted in the bracket,
e.g. GATA4(0.7|2.5). If no other field is provided, only gene names will be shown.

=back

=head1 OPTIONS

=over

=item -[-]gxref-filter [=] <expression>

SQL like filtering expression that applies to gxref table.

=item -[-]out[put] [=] <file>

The output file name. If not provided, results will be written to STDOUT.

=item -[-]over-write

If output field name already exist in input. Turn on this switch to over-write existing field in the input file.

=item -[-]gf[-]sep [=] <char>

The gene field separator. It is used to separate multiple additional values for a gene, default "|".

=for Euclid:
	char.default: "|"

=item -[-]or[-]sep [=] <char>

The output record separator, used to separate multiple genes

=for Euclid:
	char.default: ","

=item -[-]cnv[-]fields [=] <string>

Provide field names for chromosome, start, end and type in the input CNV table.
They will only be used when we query gene table based on CNV coordinates.

=for Euclid:
	string.default: "Chrom,Start,End,Type"

=item -[-]genetab [=] <dbtab>

Database and gene table used to query genes within a range, e.g. "hg19.refGene".
Or it can also be a UCSC database dump file. Gene name column (name2, 13th column) will be used
to identify each gene. The gene name nomenclature should match the one used in defining gene set. 

=item -[-]gene-over[lap] [=] <perc>

Minimial overlap percentage of a gene to be considered as fully deleted or duplicated (default: 0.99).

=for Euclid:
	perc.default: 0.99

=item -[-]exon-over[lap] [=] <perc>

Minimal overlap percentage of an exon to be considered as affected (default: 1e-10).

=for Euclid:
	perc.default: 1e-10

=item -[-]gxref-fsep [=] <string>

Gxref file: field separator. It defailts to tab for text file; "," for csv file.

=item -[-]gxref-skip [=] <number>

Gxref file: number of rows to skip at the begnning

=item -[-]gxref-sheet [=] <number>

Gxref file: if we lookup gene information from an Excel workbook, specify the sheet name or number.

=item -[-]gxref-maxcol [=] <number>

=item -[-]gxref-maxrow [=] <number>

Gxref file: if we lookup gene information from an Excel workbook, the max range of columns and rows.

=back


