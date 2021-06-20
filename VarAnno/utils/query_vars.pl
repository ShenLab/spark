use strict;
use warnings;
use File::Which;
use FindBin qw|$Bin|;
use Perl6::Slurp;
use Data::Dumper;
use List::Util qw|min max first|;
use List::MoreUtils qw|all uniq|;
use Getopt::Euclid;
use Iterator::Simple qw|iterator iter imap|;
use Utils::File::Iter qw|iter_file|;
use Genet::File::VCF qw|parse_vcf_header|;
use Genet::Var qw|normalize|;
use FaSlice;
use Genome::UCSC qw|hg_chrom|;
use Genome::UCSC::TwoBit;


use lib "$Bin/../../lib";
use Shared qw|parse_fstr|;
use Variants qw|get_fieldsinfo expand_site|;

my ($tabixflag, $tabixlib, $tabixbin);
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


# Find out fields for matching (Should always include Chrom,Position when using index)

my $dbfields; # Hash ref of all selected database fields and aliases (matching input)
my (@keyfds, @valfds); # The ordered field names in *database* for matching and collecting information
#my (@keynames, @valnames); # The ordered field alias for the above fields
my (@keycols, @valcols); # the ordered column numbers (0-based) in *database* for matching and collecting information 
my @addnames; # Additional field names to appear at the end of column names from input file
my (@keynames, @valnames); # Names of fields in the input file that will be used as key to query database and also values
my ($dbchr, $tabchr); # Chromosome nomenclature in database and input files
my ($it, $fnames);
{
	$dbfields = parse_fstr($ARGV{'--db-fields'}, 1);
	# Look for matching fields between database and input
	my ($init, $infnames) = iter_file($ARGV{'--input'} eq '-' ? \*STDIN : $ARGV{'--input'}, { fsep => qr/\t/ });
	($it, $fnames) = ($init, $infnames);
	my ($dbit, $dbfnames);
	if ($ARGV{'--dbfile'} =~ /\.vcf$/ || $ARGV{'--dbfile'} =~ /\.vcf\.gz$/) {
		$dbfnames = [qw|CHROM POS ID REF ALT QUAL FILTER|];
		my ($header, $info, undef) = parse_vcf_header($ARGV{'--dbfile'});
		push @$dbfnames, keys %{$info->{INFO}};
	}
	else {
		my %opt;
		if (all { /^\d+$/ } keys %$dbfields) {
			$opt{header} = 0;
			# Truning on fastgz may truncate the bgzipped file when parsing,
			# but it's safe if we only look at first line
			if ($ARGV{'--dbfile'} =~ /\.gz$/) {
				$opt{fastgz} = 1;
			}
		}
		($dbit, $dbfnames) = iter_file($ARGV{'--dbfile'}, { fsep => qr/\t/, %opt });
	}

	my %overwrite;
	if ($ARGV{'--over-write'}) {
		%overwrite = map { (split)[0] => 1 } split(',', $ARGV{'--over-write'});
	}
	while(my ($field, $alias) = each %$dbfields) {
		my $fcol = first { $field eq $dbfnames->[$_] } 0..$#$dbfnames;
		unless(defined $fcol) {
			die "Cannot find field $field from database file";
		}
		my $flagkey;
		if (grep { $alias eq $_ } @$infnames) {
			if (defined $overwrite{$alias}) {
				delete $overwrite{$alias};
				push @valnames, $alias;
			}
			else {
				$flagkey = 1;
				push @keynames, $alias;
			}
		}
		else {
			push @addnames, $alias;
			push @valnames, $alias;
		}
		if ($flagkey) {
			push @keyfds, $field;
			push @keycols, $fcol;
		}
		else {
			push @valfds, $field;
			push @valcols, $fcol;
		}
	}
	unless(@keyfds > 0 && @valfds > 0) {
		die "Cannot find key or value columns from database file: ", join(',', keys %overwrite);
	}
	if (%overwrite) {
		die "The following over-write fields cannot be found in the input: ", join(", ", keys %overwrite);
	}
	unless (defined $dbit) {
		open my $fdb, $ARGV{'--dbfile'} or die "Cannot open $ARGV{'--dbfile'}";
		while(<$fdb>) {
			unless(/^#/) {
				$dbchr = /^chr/ ? 1 : 0;
				last;
			}
		}
	}
	else {
		my $dat = $dbit->();
		my $chrom = $dat->{$keyfds[0]};
		$dbchr = $chrom =~ /^chr/ ? 1 : 0;
	}
}
print STDERR "The following field(s) from DB will be used as keys: ", join(',', @keyfds), "\n";
print STDERR "The following field(s) from DB will be used to fetch values: ", join(',', @valfds), "\n";


my $seq;
if ($ARGV{'--seq'}) {
	if ($ARGV{'--seq'} =~ /\.2bit$/) {
		$seq = Genome::UCSC::TwoBit->new($ARGV{'--seq'});
	}
	else {
		$seq = FaSlice->new(file => $ARGV{'--seq'});
	}
}

# Create database query hanler
# Database type will be determined by the file name suffix
my $dbh;
{
	my $db = $ARGV{'--dbfile'};
	my ($vcf, $sitetype);
	if ($db =~ /\.vcf$/ || $db =~ /\.vcf\.gz/) {
		$vcf = Genet::File::VCF->new($db);
		($sitetype, undef) = get_fieldsinfo($db);
	}
	my %inkeys;
	# When annovar db file is used, input will be slurp first, because database is very large
	if (-f "$db.idx") {
		if ($ARGV{'--input'} eq '-') {
			die "Stream input from STDIN is not supported for database with ANNOVAR index";
		}
		my $it = iter_file($ARGV{'--input'}, { fsep => qr/\t/, subset => [@keynames] });
		unless(@keynames >= 2) {
			die "Should have at least two key fields when querying database file with ANNOVAR index!";
		}
		while(my $dat = $it->()) {
			unless(defined $tabchr) {
				my $chrom = $dat->{$keynames[0]};
				$tabchr = $chrom =~ /^chr/ ? 1 : 0;
				unless($dbchr == $tabchr) {
					warn "Chromosome nomenclature in input file and DB are different!";
				}
			}

			my @keys = @{$dat}{@keynames};
			if ($dbchr == 1 && $tabchr == 0) {
				$keys[0] = hg_chrom($keys[0]);
			}
			elsif ($dbchr == 0 && $tabchr == 1) {
				$keys[0] =~ s/^chr//; $keys[0] = 'MT' if $keys[0] eq 'M';
			}
			my $key = join("\t", @keys);
			# First two key fields will be assumed to be chr and pos
			$inkeys{$key} = $keys[0].$;.$keys[1];
		}
	}

	if ($tabixflag) {
		unless(@keycols >= 2) {
			die "Incorrect number of key columns for index based database query!";
		}
		if (defined $vcf) {
			$dbh = create_tabixvcf_dbh($vcf, $sitetype);
		}
		else {
			$dbh = create_tabixfile_dbh($db);
		}
	}
	elsif (-f "$db.idx") {
		if (defined $vcf) {
			$dbh = create_idxvcf_dbh($vcf, $sitetype, \%inkeys);
		}
		else {
			$dbh = create_idxfile_dbh($db, \%inkeys);
		}
	}
	else {
		if (defined $vcf) {
			$dbh = create_vcf_dbh($vcf, $sitetype);
		}
		else {
			$dbh = create_file_dbh($db);
		}
	}
}


# Iterator the input files to add annotations
#my ($it, $fnames) = iter_file($ARGV{'--input'} eq '-' ? \*STDIN : $ARGV{'--input'}, {fsep => qr/\t/});
my $fout;
if (defined $ARGV{'--output'}) {
	$fout = IO::File->new($ARGV{'--output'}, "w");
}
else {
	$fout = \*STDOUT;	
}
print $fout join("\t", @$fnames, @addnames), "\n";


while(my $dat = $it->()) {
	my @keys = @{$dat}{@keynames};
	unless(defined $tabchr) {
		$tabchr = $keys[0] =~ /^chr/ ? 1 : 0;
		unless($dbchr == $tabchr) {
			warn "Chromosome nomenclature in input file and DB are different!";
		}
	}
	if ($dbchr == 1 && $tabchr == 0) {
		$keys[0] = hg_chrom($keys[0]);
	}
	elsif ($dbchr == 0 && $tabchr == 1) {
		$keys[0] =~ s/^chr//; $keys[0] = 'MT' if $keys[0] eq 'M';
	}

	# Iterator of DB data with matched keys
	my $iter = $dbh->(@keys);

	my %outvals;
	my $hitflag;
	while(my $value = $iter->()) {
		unless(@$value == @valnames) {
			die "Incorrect number of values for: ", join("\t", @keys);
		}
		for(my $ii = 0; $ii < @valnames; $ii ++) {
			push @{$outvals{$valnames[$ii]}}, $value->[$ii];
		}
		$hitflag = 1;
		last unless $ARGV{'--all'};
	}

	unless($hitflag) {
		foreach my $alias (@valnames) {
			$dat->{$alias} = $ARGV{'--nastr'};
		}
	}
	else {
		foreach my $alias (@valnames) {
			if ($ARGV{'--all'}) {
				$dat->{$alias} = join($ARGV{'--multi-sep'}, @{$outvals{$alias}});
			}
			else {
				$dat->{$alias} = $outvals{$alias}[0];
			}
		}
	}

	print $fout join("\t", @{$dat}{@$fnames}, @{$dat}{@addnames}), "\n";
}


#
# Prepare database query handler
#

# If tabix index is available, we query the database file for chrom and position, and slurp
# all data matched for the desginated key.
# Alternatively, relevant data from database will be slurpped first to store in a hash for
# lookup by keys.

# --- line-based database file with tabix ---
sub create_tabixfile_dbh {
	my ($dbfile) = @_;
	my $dbh;
	if ($tabixlib) {
		my $tabix = Bio::DB::HTS::Tabix->new(filename => $dbfile, use_tmp_dir => 1);
		$dbh = sub {
			my $querykey = join("\t", @_);
			my ($chr, $pos) = @_[0,1];
			my $tabixiter = $tabix->query("$chr:$pos-$pos");
			return slurp_iter(iter($tabixiter), $querykey);
		};
	}
	else {
		$dbh = sub {
			my $querykey = join("\t", @_);
			my ($chr, $pos) = @_[0,1];
			my $c_seq = $keycols[0] + 1;
			my $c_pos = $keycols[1] + 1;
			my $region = sprintf("%s:%d-%d", $chr, $pos, $pos);
			open my $fin, "$tabixbin -s $c_seq -b $c_pos -e $c_pos $dbfile $region |"
				or die "Cannot open tabix pipe for input";
			return slurp_iter($fin, $querykey);
		};
	}
	return $dbh;
}

# Helper functions
sub slurp_iter {
	my ($iter, $querykey) = @_;
	my @vals;
	while(my $line = <$iter>) {
		#print STDERR $querykey, "\t", $line, "\n";
		my ($key, $val) = line2keyval($line);
		if ($key eq $querykey) {
			push @vals, $val;
		}
	}
	return iter([@vals]);
}

sub line2keyval {
	my ($line) = @_;
	chomp($line);
	my @dat = split(/\t/, $line);
	my $key = join("\t", @dat[@keycols]);
	my @val = @dat[@valcols];
	return ($key, \@val);
}

# -- line based database file with ANNOVAR index --
sub create_idxfile_dbh {
	my ($dbfile, $inkeys) = @_;
	my ($bb, $BIN, $chrpos) = load_idx(@_);
	my %db;
	open my $fin, $dbfile or die "Cannot open $dbfile";
	foreach my $cp (sort keys %$chrpos) {
		my ($chr, $pos) = split($;, $cp);
		my $bin = $pos - ($pos%$BIN);
		if (defined $bb->{$chr,$bin}) {
			my $fh = bounded_fh($fin, $bb->{$chr,$bin});
			while(my $line = <$fh>) {
				my ($key, $val) = line2keyval($line);
				next unless defined $inkeys->{$key};
				push @{$db{$key}}, $val;
			}
		}
	}
	return hash_dbh(\%db);
}

# Helper functions
# Return a handler for querying from bash
sub hash_dbh {
	my ($dbhash) = @_;
	my $dbh = sub {
		my $querykey = join("\t", @_);
		return hash_iter($dbhash, $querykey);
	};
	return $dbh;
}

# Iterate over packed values given hash and key
sub hash_iter {
	my ($hash, $key) = @_;
	if (defined $hash->{$key}) {
		my $vals = $hash->{$key};
		return iter($vals);
	}
	else {
		return iter([]);
	}
}

# Load ANNOVAR index on candidate sites
sub load_idx {
	my ($dbfile, $inkeys) = @_;
	my %chrpos;
	foreach my $cp (values %$inkeys) {
		$chrpos{$cp} = 1;
	}
	my %db;
	my ($bb, $BIN) = read_annovaridx($dbfile, \%chrpos);
	return ($bb, $BIN, \%chrpos);
}

sub read_annovaridx {
	my ($dbfile, $chrpos) = @_;
	unless(-f $dbfile) {
		die "Cannot find DB file: $dbfile";
	}
	unless(-f "$dbfile.idx") {
		die "Cannot find DB index file $dbfile.idx"
	}
	open my $fidx, "$dbfile.idx" or die "Cannot open $dbfile.idx!";
	my $line = <$fidx>;
	my ($BIN, $DBSIZE);
	if ($line =~ m/BIN\t(\d+)\t(\d+)/) {
		($BIN, $DBSIZE) = ($1, $2);
	}
	else {
		die "Malformed database index header line: $line";
	}

	my %keep;
	foreach my $cp (keys %$chrpos) {
		my ($chr, $pos) = split($;, $cp);
		my $bin = $pos - ($pos % $BIN);
		$keep{$chr,$bin} = 1;
	}

	my %bb;
	while(<$fidx>) {
		my ($chr, $bin, $offset0, $offset1) = split;
		next unless defined $keep{$chr,$bin};
		$bb{$chr,$bin} = "$offset0,$offset1";
	}
	return (\%bb, $BIN);
}

# Return a bounded input file handle
sub bounded_fh {
	my ($fin, $range) = @_;
	my ($cmin, $cmax) = split(',', $range);
	seek($fin, $cmin, 0);
	my $chere = $cmin;
	return iterator {
		return if $chere > $cmax;
		my $line = <$fin>;
		$chere += length($line);
		return $line;
	};
}


# --- line-based database file without any index --
sub create_file_dbh {
	my ($dbfile) = @_;
	my %db;
	my $it = iter_file($dbfile, { fsep => qr/\t/ });
	while(my $dat = $it->()) {
		my $key = join("\t", @{$dat}{@keyfds});
		push @{$db{$key}}, [@{$dat}{@valfds}];
	}
	return hash_dbh(\%db);
}


# --- bgzipped VCF file with tabix ---
sub create_tabixvcf_dbh {
	my ($vcf, $sitetype) = @_;
	my $dbh;
	if ($tabixlib) {
		my $tabix = Bio::DB::HTS::Tabix->new(filename => $vcf->{file}, use_tmp_dir => 1);
		$dbh = sub {
			my ($chr, $pos) = @_[0,1];
			my $querykey = join("\t", @_);
			my $region = sprintf("%s:%d-%d", $chr, max(1, $pos-$ARGV{'--padding'}), $pos+$ARGV{'--padding'});
			my $tabixiter = $tabix->query($region);
			$vcf->{ftype} = 1;
			$vcf->{file} = iter($tabixiter);
			my $vcfit = $vcf->iter({ samp => [] });
			my $dbhash = slurp_vcf($vcfit, $sitetype);
			return hash_iter($dbhash, $querykey);
		};
	}
	else {
		$dbh = sub {
			my ($chr, $pos) = @_[0,1];
			my $querykey = join("\t", @_);
			my $region = sprintf("%s:%d-%d", $chr, max(1, $pos-$ARGV{'--padding'}), $pos+$ARGV{'--padding'});
			my $vcfit = $vcf->iter({ samp => [], region => $region });
			my $dbhash = slurp_vcf($vcfit, $sitetype);
			return hash_iter($dbhash, $querykey);
		};
	}
	return $dbh;
}

# Helper function
sub slurp_vcf {
	my ($vcfit, $sitetype, $dbref, $inkeys) = @_;
	if (ref $vcfit eq 'Genet::File::VCF') {
		$vcfit = $vcfit->iter({ samp => [] });
	}
	my %db;
	unless(defined $dbref) {
		$dbref = \%db;
	}
	while(my ($site, $geno) = $vcfit->()) {
		# Only for some public available VCF that may not conform to standard format.
		next unless defined $site->{ALT}; 
		my @sites = expand_site($site, $sitetype);
		for(my $ii = 0; $ii < @sites; $ii ++) {
			next if $sites[$ii]{ALT} eq '*';

			if (defined $seq && $sites[$ii]{REF} =~ /^[ACGT]+$/ && $sites[$ii]{ALT} =~ /^[ACGT]+$/) {
				normalize($sites[$ii], $seq);
			}
			foreach my $site (@sites) {
				my $key = join("\t", @{$site}{@keyfds});
				if (defined $inkeys) {
					next unless defined $inkeys->{$key};
				}
				my @val = map { ref $_ eq 'ARRAY' ? join(',', @$_) : $_ } @{$site}{@valfds};
				push @{$dbref->{$key}} => \@val;
			}
		}
	}	
	return $dbref;
}

# --- VCF file with ANNOVAR index --
sub create_idxvcf_dbh {
	my ($vcf, $sitetype, $inkeys) = @_;
	my ($bb, $BIN, $chrpos) = load_idx($vcf->{file}, $inkeys);
	my %db;
	open my $fin, $vcf->{file} or die "Cannot open $vcf->{file}";
	foreach my $cp (sort keys %$chrpos) {
		# first create an bounded file handle
		my ($chr, $pos) = split($;, $cp);
		my $bin = $pos - ($pos%$BIN);
		if (defined $bb->{$chr,$bin}) {
			my $fh = bounded_fh($fin, $bb->{$chr,$bin}); 
			$vcf->{file} = $fh;
			$vcf->{ftype} = 1;
			my $vcfit = $vcf->iter({ samp => [] });
			slurp_vcf($vcf, $sitetype, \%db, $inkeys);
		}
	}
	return hash_dbh(\%db);
}


# --- VCF file without any index ---
sub create_vcf_dbh {
	my $dbhash = slurp_vcf(@_);
	return hash_dbh($dbhash);
}

__END__

=head1 NAME

query_var.pl -- Annotating variant specific scores.

=head1 NOTES

The script queries variant level information from database files.

It extends the ANNOVAR's database table index scheme to support header line fields.
It also extends the query_genevar utility to support VCF file that can match any extra database fields. 
Note that query_genevar only matches GeneID, but also support matching AA changes.

The support for VCF is similar in principle as used by VEP, but we have an option perform internal 
variant normalization that works better for matching indels. And we also "expand" INFO fields
that are of type R,G,or number >2. In such cases, the corresponding database field name specified 
should also be the "expanded" version.

For range-based annotation, use range_intersect utility for BED file. We also have overlap_vars utility
for more functions including the support for bigWig formats.

When index to the database file is not available, it will slurp all the data to store in a hash 
for hash-based query (suitable for small data set). And this approach can be used to query information
not just for variants. In this sense, it also extend the anno_genes utility for matching multiple
fields between input and database. But this script does NOT "expand" packed data within a row (e.g. 
for overlapping genes), nor does it have priority based matching scheme in anno_genes.


=head1 REQUIRED ARGUMENTS

=over 

=item -[-]in[put] [=] <table>

Input table file. Must be tab-separated. Use "-"" to stream from STDIN.

=item -[-]db[file] [=] <file>

The database file used for query variant-specific information. 
The following types of database files are supported (Must be tab-separated):
Line-based text file with or without ANNOVAR index, line-based bgziped database file with tabix index,
or VCF file with or without index.

=item -[-]db[-]fields [=] <string>

Fields in the database used for matching and output. Alias to the database field is used for matching
fields from input file. If index file are available, first two fields will be used as chrom and position.
All fields that shared between database file and input (after aliasing) minus those to be over-written will 
be used for matching.

=back

=head1 OPTIONS

=over

=item -[-]over[-]write [=] <fields>

Specify the fields from the input to be over-written. 

=item -[-]out[put] [=] <table>

Output table file name. If not provided, will dump to STDOUT. Fields will be tab separated.

=item -[-]all

Extract and output all matched terms. Default is report the first one that matches.

=item -[-]multi-sep [=] <string>

Separator for multiple output records, useful when --all option is on.

=for Euclid:
	string.default: ","

=item -[-]nastr [=] <string>

Output value for missing data, default ".".

=for Euclid:
	string.default: "."

=item -[-]pad[ding] [=] <length>

Length padded to both sides of the query position. Default is 0, suitable for lookup SNVs.
Set the padding to greater than 0 if VCF dbfile also includes indels. 

=for Euclid:
	length.default: 0


=back



