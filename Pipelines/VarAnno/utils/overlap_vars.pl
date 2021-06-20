use strict;
use warnings;
use Data::Dumper;
use FindBin qw|$Bin|;
use List::Util qw|max first|;
use List::MoreUtils qw|any all uniq|;
use File::Which;
use File::Basename;
use Getopt::Euclid;
use Utils::List qw|parse_fields|;
use Utils::File::Iter qw|iter_file|;
use Genome::UCSC qw|hg_chrom|;
use Genome::Ranges qw|range_overlap|;
use Genome::UCSC::BinKeeper;
use Bio::DB::BigFile;

use lib "$Bin/../../lib";
use Shared qw|region_bk parse_fstr parse_tabfile|;
use Variants qw|var_type parse_operations|;


# Tabix and bigwig parsing library or utiliy flag
my ($tabixflag, $tabixlib, $tabixbin, $bigwigflag, $bigwiglib, $bigwigbin);
if ($ARGV{'--dbfile'} =~ /\.gz$/ && -f "$ARGV{'--dbfile'}.tbi") {
	$tabixflag = 1;
	eval {
		require Bio::DB::HTS::Tabix;
	};
	unless($@){
		$tabixlib = 1;
	}
	else {
		$tabixbin = which("tabix");
		die "Cannot find tabix utility" unless $tabixbin;
	}
}
if ($ARGV{'--dbfile'} =~ /\.bw$/ || $ARGV{'--dbfile'} =~ /\.bigwig$/) {
	$bigwigflag = 1;
	eval {
		require Bio::DB::BigFile;
	};
	unless($@) {
		$bigwiglib = 1;
	}
	else {
		$bigwigbin = which("bigWigToBedGraph");
		die "Cannot find bigWigToBedGraph utility" unless $bigwigbin;
	}
}

# Validate field specifications for input and database file
my $dbfields; # database fields and alias
my (@keycols, @valcols); # Columns of key and value fields in database file (for text/bed files only)
my @valnames; # Field name alias for values to be collected (key columns are standardized)
my ($dbchr, $tabchr); # Chromosome nomenclature in database and input files
{
	my @infields = parse_fields($ARGV{'--in-fields'});
	unless(@infields == 4) {
		die "Incorrect number of input file fields: $ARGV{'--in-fields'}";
	}
	
	my %stdfields = (Chrom => 1, Start => 1, End => 1);
	$dbfields = parse_fstr($ARGV{'--db-fields'}, 1);
	
	@valnames = grep { !defined $stdfields{$_} } values %$dbfields;
	if (@valnames == 0) {
		die "Cannot find value field name for $ARGV{'--dbfile'}";
	}

	if ($bigwigflag) {
		if (@valnames > 1) {
			die "Can only specify one value field for bigwig file!";
		}
		if ($bigwiglib) {
			my $wig = Bio::DB::BigFile->bigWigFileOpen($ARGV{'--dbfile'});
			my $list = $wig->chromList();
			my $first = $list->head();
			$dbchr = $first->name =~ /^chr/ ? 1 : 0;
		}
		else {
			open my $fin, "$bigwigbin $ARGV{'--dbfile'} stdout | head |" or die "Cannot open $bigwigbin pipe for input";
			my $firstline = <$fin>;
			$dbchr = $firstline =~ /^chr/ ? 1 : 0;
		}
	}
	else {
		unless(join(",", (values %{$dbfields})[0,1,2]) eq "Chrom,Start,End") {
			die "The first three fields from dbfile must be Chrom,Start,End!";
		}

		my %opt;
		if (all { /^\d+$/ } keys %$dbfields) {
			$opt{header} = 0;
			# Truning on fastgz may truncate the bgzipped file when parsing,
			# but it's safe if we only look at first line
			if ($ARGV{'--dbfile'} =~ /\.gz$/) {
				$opt{fastgz} = 1;
			}
		}
	
		my ($dbit, $dbfnames) = iter_file($ARGV{'--dbfile'}, { fsep => qr/\t/, fastgz => 1, %opt });
		my $f_chr;
		foreach my $dbfield (keys %$dbfields) {
			my $ii = first {  $dbfnames->[$_] eq $dbfield } 0..$#$dbfnames;
			unless(defined $ii) {
				die "Cannot find dbfield $dbfield from dbfile: $ARGV{'--dbfile'}";
			}
			my $dbalias = $dbfields->{$dbfield};
			if ($dbalias eq 'Chrom') {
				$f_chr = $dbfield;
			}
			if (defined $stdfields{$dbalias}) {
				push @keycols, $ii; 
			}
			else {
				push @valcols, $ii;
			}
		}
		my $first = $dbit->();
		if ($first->{$f_chr} =~ /^chr/) {
			$dbchr = 1;
		}
		else {
			$dbchr = 0;
		}
	}
}

my ($it, $fnames, $keyfields) = parse_tabfile($ARGV{'--input'}, $ARGV{'--in-fields'}, 4, 4);

my $dbh;
{
	if ($tabixflag) {
		$dbh = create_tabixfile_dbh($ARGV{'--dbfile'});
	}
	elsif ($bigwigflag) {
		$dbh = create_bigwig_dbh($ARGV{'--dbfile'});
	}
	else {
		$dbh = create_binkeeper_dbh($ARGV{'--dbfile'});
	}
}


my $ops;
if (defined $ARGV{'--operation'}){
	$ops = parse_operations($ARGV{'--operation'}, {dbfields => $dbfields });
}
else {
	$ops = parse_operations(join(";", map { "first($_):$_" } @valnames), { dbfields => $dbfields });
}

my @outfields = @$fnames;
foreach my $opname (keys %$ops) {
	if (grep { $_ eq $opname } @$fnames) {
		unless($ARGV{'--over-write'}) {
			die "The operation output column $opname already exist in the output";
		}
	}
	else {
		push @outfields, $opname;
	}
}


my $fout;
if($ARGV{'--output'}) {
	open $fout, ">$ARGV{'--output'}" or die "Cannot write to $ARGV{'--output'}";
}
else {
	$fout = \*STDOUT;
}
print $fout join("\t", @outfields), "\n";


while(my $dat = $it->()) {
	my ($chrom, $pos, $ref, $alt) = @{$dat}{@$keyfields};
	unless(defined $tabchr) {
		if ($chrom =~ /^chr/) {
			$tabchr = 1;
		}
		else {
			$tabchr = 0;
		}
		unless($dbchr == $tabchr) {
			warn "Chromosome nomenclature in DB and input are different!";
		}
	}
	if ($dbchr == 1 && $tabchr == 0) {
		$chrom = hg_chrom($chrom);
	}
	elsif ($dbchr == 0 && $tabchr == 1) {
		$chrom =~ s/^chr//; $chrom = 'MT' if $chrom eq 'M';
	}

	my $vartype = var_type($ref, $alt);

	my ($start, $end);
	if ($vartype eq 'SNV') {
		$start = $pos;
		$end = $pos;
	}
	else {
		if ($ARGV{'--snv-only'}) {
			print $fout join("\t", map { $dat->{$_} // $ARGV{'--nastr'} } @outfields), "\n";
			next;
		}
		if ($vartype eq 'MNV' || length($ref) > length($alt)) {
			$start = $pos;
			$end = $pos + length($ref) - 1;
		}
		elsif (length($ref) < length($alt)) {
			$start = $pos;
			$end = $pos + length($ref);
		}
		else {
			die "Cannot determine the range for $chrom:$start:$ref>$alt";
		}
	}

	my @vals = $dbh->($chrom, $start, $end);

	if (@vals) {
		while(my ($outfield, $callback) = each %$ops) {
			$dat->{$outfield} = $callback->(@vals);
			if (!defined $dat->{$outfield} || $dat->{$outfield} eq "") {
				$dat->{$outfield} = $ARGV{'--nastr'};
			}
		}
	}
	print $fout join("\t", map { $_ // $ARGV{'--nastr'} } @{$dat}{@outfields}), "\n";
}

#
# Prepare database query handler
#
# -- binkeeper ==
sub create_binkeeper_dbh {
	my ($dbfile) = @_;
	# region_bk  will automatically adjust for the 0-based coordinate system based on file name suffix
	my $bk = region_bk($dbfile, $dbfields, { bed => $ARGV{'--db-bed'} });
	my $dbh = sub {
		my ($chrom, $start, $end) = @_;
		my @vals = map { $_->[2] } $bk->find_range($chrom, $start, $end, { qCover => $ARGV{'--overlap'} });
		return @vals;
	}
}

# -- line-based databasae file with tabix (including bed file) --
sub create_tabixfile_dbh {
	my ($dbfile) = @_;
	unless($dbfile =~ /\.gz$/ && -f "$dbfile.tbi") {
		die "Database file must be bgzipped and have tabix index";
	}
	unless(@keycols == 3) {
		die "Database file must have three key columns: Chrom,Start,End";
	}
	unless(@valnames == @valcols) {
		die "Incorrect number of value columns";
	}
	my $bedflag = $ARGV{'--db-bed'};
	if ($dbfile =~ /\.bed\.gz$/ || $dbfile =~ /\.bedgraph\.gz$/) {
		$bedflag = 1;
	}
	my $dbh;
	if ($tabixflag) {
		my $tabix = Bio::DB::HTS::Tabix->new(filename => $dbfile, use_tmp_dir => 1);
		$dbh = sub {
			my ($chrom, $start, $end) = @_;
			my @vals;
			my $region = sprintf("%s:%d-%d", $chrom, max(0, $start-$ARGV{'--padding'}), $end+$ARGV{'--padding'});
			my $it = $tabix->query($region);
			while(my $line = $it->next) {
				chomp($line);
				my @dat = split(/\t/, $line);
				my $dbstart = $bedflag ? $dat[$keycols[1]]+1 : $dat[$keycols[1]];
				my $dbend = $dat[$keycols[2]];
				my $intersect = range_overlap($dbstart, $dbend, $start, $end);
				if ($intersect >= $ARGV{'--overlap'}*($end-$start+1)) {
					my %val;
					@val{@valnames} = @dat[@valcols];
					push @vals, \%val;
				}
			}
			return @vals;
		};
	}
	else {
		$dbh = sub {
			my ($chrom, $start, $end) = @_;
			my @vals;
			my @cpos = map { $_+1 } @keycols;
			my $region = sprintf("%s:%d-%d", $chrom, max(0, $start-$ARGV{'--padding'}), $end+$ARGV{'--padding'});
			open my $fin, "$tabixbin -s $cpos[0] -b $cpos[1] -e $cpos[2] $dbfile $region |"
				or die "Cannot open tabix pipe for input";
			while(<$fin>) {
				chomp;
				my @dat = split(/\t/);
				my $dbstart = $bedflag ? $dat[$keycols[1]]+1 : $dat[$keycols[1]];
				my $dbend = $dat[$keycols[2]];
				my $intersect = range_overlap($dbstart, $dbend, $start, $end);
				if ($intersect >= $ARGV{'--overlap'}*($end-$start+1)) {
					my %val;
					@val{@valnames} = @dat[@valcols];
					push @vals, \%val;
				}
			}
			return @vals;
		};
	}
	return $dbh;
}
	

# -- Query overlapping intervals from bigwig file --
# The intervals collected from bigwig should by default be 0-based (bedflag=1 by default)
sub create_bigwig_dbh {
	my ($dbfile) = @_;
	if (@valnames > 1) {
		die "Can only collect one value field from a bigwig file!";
	}
	my $dbh;
	if ($bigwiglib) {
		my $wig = Bio::DB::BigFile->bigWigFileOpen($dbfile);
		$dbh = sub {
			my ($chrom, $start, $end) = @_;
			my $intervals = $wig->bigWigIntervalQuery($chrom, max(0, $start-$ARGV{'--padding'}), $end+$ARGV{'--padding'});
			my @vals;
			my $intv = $intervals->head;
			if ($intv) {
				my $intersect = range_overlap($intv->start+1, $intv->end, $start, $end);
				if ($intersect >= $ARGV{'--overlap'}*($end-$start+1)) {
					push @vals, { $valnames[0] => $intv->value };
				}
			}
			return @vals;
		};
	}
	else {
		$dbh = sub {
			my ($chrom, $start, $end) = @_;
			my $qstart = max(0, $start - $ARGV{'--padding'});
			my $qend = $end + $ARGV{'--padding'};
			open my $fin, "$bigwigbin -chrom=$chrom -start=$qstart -end=$qend $dbfile stdout |"
				or die "Cannot open $bigwigbin pipe for input";
			my @vals;
			while(<$fin>) {
				my @dat = split;
				my $intersect = range_overlap($dat[1]+1, $dat[2], $start, $end);
				if ($intersect >= $ARGV{'--overlap'}*($end-$start+1)) {
					push @vals, { $valnames[0] => $dat[3] };
				}
			}
			return @vals;
		};
	}
	return $dbh;
}


__END__

=head1 NAME

overlap_vars.pl -- Find and summarize overlapping features for sequence variants.

=head1 NOTES

Typical usage of this script is to look for overlapping genomic features of a sequence variants (e.g.
evolutionary conserved region or functional protein domains), or to query a score for the genomic
position of a sequence variant. Although precise allele is not directly used, the allele length
will be used in determining the genomic interval defined by the variant.

We have the following rule for defining the genomic position of indels, insertion is defined by the
two flanking bases of insertion point, deletion or MNV is defined by the span of deleted sequence. 
There is an option to ignore indels/MNVs.

UCSC BED file, bgzipped BED file or bigwig can be used as database files, we also support general line-based
text file with tabix index. Similar to overlap_cnv utlity, one or more group operations can be specified to 
create summary output appended to input.

Genomic coordinate in BED file is zero-based. And in all other format, it should be 1-based.
Note bedGraph converted from bigWig format should also have 0-based coordinates.

=head1 REQUIRED ARGUMENTS

=over

=item -[-]in[put] [=] <file>

Input sequence variant table. It must be tab-separated with required fields.

=item -[-]db[file] [=] <dbfile>

The database file of interval-based genomic features. File name suffix .bed,.bedgreph,.bw will be automatically recognzed.

=item -[-]db[-]fields [=] <fstr>

Database fields. The first three are required and must be corresponds to Chrom,Start,End (not needed to bigWig format).
The remaining fields are optional and will be used by operation.

=back

=head1 OPTIONS

=over

=item -[-]op[eration] [=] <fstr>

One or more operations to db fields and their alias as output column name.
Default is "first(FieldName):FieldName" for all non-standard field names from database.
The operation syntax can be found in documentation of overlap_cnvs.pl 
It supports all groupby operations of BEDtools: https://bedtools.readthedocs.io/en/latest/content/tools/groupby.html

=item -[-]out[put] [=] <file>

Output file name. If not provided, will be written to STDOUT.

=item -[-]in[-]fields [=] <fstr>

The input file field name for Chrom,Position,Ref,Alt. Must be in the specified order without using alias.

=for Euclid:
	fstr.default: "Chrom,Position,Ref,Alt"

=item -[-]over[-]write

Ovwr-write existing columns of the same name in input.

=item -[-]overlap [=] <frac>

Min required fraction of variant range that are overlapped by the genomic feature. Must be greater than 0.

=for Euclid:
	frac.default: 0.5

=item -[-]snv-only

Only annotate SNVs.

=item -[-]db-bed 

Indicate that genomic intervals in database are BED-style (i.e., half-open 0-based).
By default, we infer coordinate system based on file name suffix. 

=item -[-]pad[ding] [=] <length>

Add padding length to the query interval when using tabix or bigwig tools to account for edge cases.

=for Euclid:
	length.default: 1

=item -[-]nastr [=] <string>

Output value for missing data, default ".".

=for Euclid:
	string.default: "."

=back




