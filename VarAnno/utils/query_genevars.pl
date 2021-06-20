use strict;
use warnings;
use File::Which;
use FindBin qw|$Bin|;
use List::Util qw|first|;
use List::MoreUtils qw|all any uniq|;
use Getopt::Euclid;
use Data::Dumper;
use Iterator::Simple qw|iter|;
use Utils::Hash qw|chk_default|;
use Utils::File qw|open_file|;
use Utils::File::Iter qw|iter_file|;
use Genome::UCSC qw|hg_chrom|;

use lib "$Bin/../../lib";
use Shared qw|parse_fstr expand_dat|;

# We will test if Bio::DB::HTS::Tabix  is installed on the machine
# If not, we will fall back to calling tabix utility
my ($tabix_lib, $tabix_bin) = ('N', '');
eval {
	require Bio::DB::HTS::Tabix;
};
unless($@){
	$tabix_lib = "Y";
}
else {
	$tabix_bin = which("tabix");
	die "Cannot find tabix utility" unless $tabix_bin;
}

my ($dbchr, $tabchr);

# Parse database fields and decide which columns to use for tabix, and which fields to appear in output
# Also determine if chr prefix exist in db file.
my $fields = parse_fstr($ARGV{'--db-fields'}, 1);
# When parsing tab separated dbfiles, fnum will store the column number for given field alias
my %fnum;
# All non-standard fields will be used as output fields
my @outfields;
{	
	my $fin = open_file($ARGV{'--dbfile'});
	my $header = <$fin>; chomp($header);
	my @infields = split(/\t/, $header);
	if (all { /^\d+$/ } keys %$fields) {
		while(my ($ii, $alias) = each %$fields) {
			chk_default(\%fnum, $alias, $ii);
			if ($ii > @infields) {
				die "Cannot find column $ii from dbfile";
			}
		}
	}
	else {
		while(my ($field, $alias) = each %$fields) {
			my $ii = first { $infields[$_-1] eq $field } 1..@infields;
			unless(defined $ii) {
				die "Cannot find field $field from dbfile";
			}
			else {
				chk_default(\%fnum, $alias, $ii);
			}
		}
	}
	# Now check if all standard fields can be found
	foreach my $alias (qw|Chrom Position Ref Alt|) {
		unless(defined $fnum{$alias}) {
			die "Cannot find standard field $alias after aliasing";
		}
	}
	my %stdfields = (Chrom => 1, Position => 1, Ref => 1, Alt => 1, GeneID => 1, AAChg => 1);
	@outfields = grep { !defined $stdfields{$_} } values %$fields;

	unless(defined $fnum{GeneID} || defined $fnum{AAChg}) {
		die "Cannot find at least one of GeneID or AAChg column";
	}
	unless(@outfields) {
		die "Cannot find output fields";
	}
	my $line = <$fin>;
	my $chrom = (split(/\t/, $line))[$fnum{Chrom}-1];
	if ($chrom =~ /^chr/) {
		$dbchr = 1;
	}
	else {
		$dbchr = 0;
	}
}	

my ($it, $fnames) = iter_file($ARGV{'--input'} eq '-' ? \*STDIN : $ARGV{'--input'}, {fsep => qr/\t/});
my @expandfds = qw|GeneID|;
{
	my @stdinfields = qw|Chrom Position Ref Alt GeneID|;
	if (defined $fnum{AAChg} || $ARGV{'--misonly'}) {
		push @stdinfields, 	qw|GeneEff TransIDs TransEffs AAChg|;
		push @expandfds, qw|GeneEff TransIDs TransEffs AAChg|;
	}
	foreach my $field (@stdinfields) {
		unless(any { $_ eq $field } @$fnames) {
			die "Cannot find standard field $field from input file";
		}
	}
}

# Appending additional columns
my @colnames = @$fnames;
foreach my $field (@outfields) {
	if (grep { $field eq $_ } @$fnames) {
		unless($ARGV{'--over-write'}) {
			die "Field $field already exist in the input file, please rename to a different alias!";
		}
		else {
			warn "Field $field in the input file will be over-written.";
		}
	}
	else {
		push @colnames, $field;
	}
}

my $fout;
if (defined $ARGV{'--output'}) {
	$fout = IO::File->new($ARGV{'--output'}, "w");
}
else {
	$fout = \*STDOUT;	
}
print $fout join("\t", @colnames), "\n";

my $dbh;
if ($tabix_lib eq 'Y') {
	my $tabix = Bio::DB::HTS::Tabix->new(filename => $ARGV{'--dbfile'}, use_tmp_dir => 1);
	$dbh = sub {
		my ($chr, $pos) = @_;
		#my $region = sprintf("%s:%d-%d", $chr, $pos, $pos);
		return $tabix->query("$chr:$pos-$pos");
	};
}
else {
	$dbh = sub {
		my ($chr, $pos) = @_;
		my $region = sprintf("%s:%d-%d", $chr, $pos, $pos);
		open my $fin, "$tabix_bin -s $fnum{Chrom} -b $fnum{Position} -e $fnum{Position} $ARGV{'--dbfile'} $region |"
			or die "Cannot open tabix pipe for input";
		return iter($fin);
	};
}

my @f_keys = keys %fnum;
my @f_idx = map { $_ - 1 } values %fnum;

# If gene ID is provided, we will use GeneID, otherwise matching by AAChg
# To query annotations other than missense predictions, we only use GeneID for matching

while(my $dat = $it->()) {
	my %outvals;
	unless(defined $tabchr) {
		if($dat->{Chrom} =~ /^chr/) {
			$tabchr = 1;
		}
		else {
			$tabchr = 0;
		}
		unless($tabchr == $dbchr) {
			warn "Chromosome nomenclatures in variant table and database are different";
		}
	}

	my $chrom = $dat->{Chrom};
	if ($tabchr == 1 && $dbchr == 0) {
		$chrom =~ s/^chr//; $chrom = 'MT' if $chrom eq 'M';
	}
	elsif ($tabchr == 0 && $dbchr == 1) {
		$chrom = hg_chrom($chrom);
	}

	foreach my $data (expand_dat($dat, { fields => \@expandfds, sep => ';' })) {
		if ($ARGV{'--misonly'} && $data->{GeneEff} !~ /^missense/ && $data->{TransEffs} !~ /missense/) {
			foreach my $field (@outfields) {
				push @{$outvals{$field}}, $ARGV{'--nastr'};
			}
			next;
		}
		# Parse AAChg if the field was specified and exists in dbfile
		my $aachg;
		if (defined $fnum{AAChg}) {
			my @transeffs = split(',', $data->{TransEffs});
			my @aachgs = split(',', $data->{AAChg});	
			my $jj = first { $transeffs[$_] eq 'missense' } 0..$#transeffs;
			#unless(defined $jj) {
				#my $varid = join(":", @{$dat}{qw|Chrom Position Ref Alt|});
				#die "Cannot find missense variants for $varid";
			#}
			if (defined $jj) {
				if ($aachgs[$jj] =~ /^p\.([A-Z]+)\d+([A-Z]+)$/) {
					$aachg = $1."/".$2;
				}
				elsif ($aachgs[$jj] =~ /\d+\:([A-Z]+\/[A-Z]+)$/) {
					$aachg = $1;
				}
				else {
					warn "Cannot parse $aachgs[$jj] to find AA change";
				}
			}
		}
		# Search the database
		my $iter = $dbh->($chrom, $dat->{Position});
		my $hitflag;
		my %hits;
		while(my $line = $iter->next) {
			chomp($line);
			my @values = split(/\t/, $line);
			my %info;
			@info{@f_keys} = @values[@f_idx];		
			if ($info{Chrom} eq $chrom && $info{Position} == $dat->{Position} && 
				$info{Ref} eq $dat->{Ref} && $info{Alt} eq $dat->{Alt} && 
				(defined $info{GeneID} && $info{GeneID} eq $data->{GeneID} ||
				 defined $info{AAChg} && defined $aachg && $info{AAChg} eq $aachg)) {
				foreach my $field (@outfields) {
					push @{$hits{$field}}, $info{$field};
				}
				$hitflag = 1;
				last unless $ARGV{'--all'}; 
			}
		}
		unless($hitflag) {
			foreach my $field (@outfields) {
				push @{$outvals{$field}}, $ARGV{'--nastr'};
			}
		}
		else {
			foreach my $field (@outfields) {
				if ($ARGV{'--all'}) {
					push @{$outvals{$field}}, $hits{$field};
				}
				else {
					if (@{$hits{$field}} > 1) {
						my $varid = join(":", @{$dat}{qw|Chrom Position Ref Alt|});
						die "More than one values of $field for $varid were collected for under first hit mode";
					} 
					push @{$outvals{$field}}, $hits{$field}[0];
				}
			}
		}
	}
	# Appending to the output
	foreach my $field (@outfields) {
		if (defined $outvals{$field}) {
			my @geneids = split(';', $dat->{GeneID});
			unless(scalar(@{$outvals{$field}}) == scalar(@geneids)) {
				print STDERR scalar(@{$outvals{$field}}), "<>", scalar(@geneids), "\n";
				my $varid = join(":", @{$dat}{qw|Chrom Position Ref Alt|});
				die "Output values for $field does not have the same length as geneids: $varid";
			}
			my @outval;
			if ($ARGV{'--all'}) {
				@outval = map { ref $_ eq 'ARRAY' ? join($ARGV{'--multi-sep'}, @$_) : $_ } @{$outvals{$field}};
			}
			else {
				if (any { ref $_ eq 'ARRAY' } @{$outvals{$field}}) {
					die "More than one predictions were collected for $field";
				}
				@outval = @{$outvals{$field}};
			}
			my @uqvals = uniq sort @outval;
			if (@uqvals == 1) {
				$dat->{$field} = $uqvals[0];
			}
			else {
				$dat->{$field} = join(";", @outval);
			}
 		}
		else {
			$dat->{$field} = $ARGV{'--nastr'}; #join(";", (".") x scalar(@geneids));
		}
	}
	print $fout join("\t", @{$dat}{@colnames}), "\n";
}


__END__

=head1 NAME

query_genevar.pl -- Annotating gene-specific variant predictions/scores.

=head1 NOTES

To annotate the predicted pathogenicity or deleteriousness of missense variants from pre-computed results,
we need to match variants with not only genomic coordinates, but also additional information including
gene/transcript ID or amino acid changes. 

The default behavior of ANNOVAR is take the first record in its database matching genomic coordinates.
This is not correct in case of overlapping genes. So we initially developed this script to assign gene-specific 
annotations for each missense variant, and later generalized to look for any gene-specific variant prediction
scores.

In principle, missense annotation should be transcript-specific which is adopted in VEP's built in annotation
of PolyPhen and SIFT. In reality, it is not very informative to have different predictions for the same AA change 
in different transcripts. So we will only give gene-specific predictions for each variant. We have two approaches
for matching genes: 
1. GeneID: this is quick and accurate when the gene set used to generating the lookup table and that used
in the annotation are the same. But it's not always the case that tool developer and users of annotation pipeline
can use the same set of GeneIDs.
2. AAChg: to account for possible gene ID changes, we can also match by amino acid changes. This will be reasonably
accurate in most cases.

For the second approach, there may be some very specicial cases when the same missense variant in one gene have 
different AA changes because the reading frames of different transcripts are not in frame (not sure if it really exist?).
If this happens, we will only fetch the prediction for the *first* matched amino acid change.

The database files will be tabix indexed for position-based query and will include gene ID and/or amino acid changes 
for matching against querying variants.

Although the missense prediction should be gene-specific, some overlapping genes have highly similar sequences
and also may have the same gene name. In overlapping gene regions, if prediction for the same missense variant
are the same across all genes, only a single prediction result will be shown.

=head1 REQUIRED ARGUMENTS

=over

=item -[-]in[put] [=] <table>

Input annotated variant table. Should be the tabular format collected from VEP outputs, 
including standard fields of Chrom,Pos,Ref,Alt,GeneID,(GeneEff,TransIDs,TransEffs,AAChg).
We will use GeneID for ID-based match, and AAChg for amino acids-based match.
GeneID must be present in the input even if it was not used for matching.
Note: use '-' as file name to stream from STDIN.

=item -[-]db[file] [=] <file>

The database file used for querying annotations. 
Note: we only support line-based bgzipped database table files indexed by tabix. VCF is not supported!

=item -[-]db[-]fields [=] <string>

Fields in the database used for matching and output.
Example1: chrom:Chrom,pos:Position,ref:Ref,alt:Alt,ENSG:GeneID,Amino_acids:AAChg,MPC
Example2: 1:Chrom,2:Position,3:Ref,4:Alt,5:AAChg,6:REVEL

Fields with alias Chrom,Position,Ref,Alt,GeneID,AAChg will be used to matching variants.
Different approach of variant matching will be chosen depending on the presence of GeneID or AAChg
fields. If both fields are present, we will give priority to geneID-based match.
Always check the version of gene set used in annotation to decided if GeneID-based matching is appropriate.  

AAChg field in the database file should always have "R/T" format.
 
All other fields will appear as additional columns in the output variant table. Additional output fields 
should not exist in input file.

=back

=head1 OPTIONS

=over

=item -[-]over[-]write

Under over-writing mode, the existing columns from the input file can be over-written.

=item -[-]out[put] [=] <outfile>

The output file name. Additional columns will be appended.
If not provided, results will be written to STDOUT.

=item -[-]mis[-]only

Restrict to variants whose gene-based annotations are missense. Note: it is possible that some variants
were annotated as missense in some transcript, but LGD in other transcript. We consider a variant as missense
if its functional consequence on *any one* transcript was annotated as missense. 
Truning on this may reduce the runtime by skipping some variants. The default is to lookup database for all 
SNVs, and can be used to generally look up any gene-specific variant predictions.

=item -[-]all

The default is to extract prediction related fields from the first matched record.
Turning on this option will extract and output all predictions.
When this option is turned on, we will not merge identical predicionts.

=item -[-]multi-sep [=] <string>

Separator for multiple output records, useful when --all option is on.

=for Euclid:
	string.default: ","

=item -[-]nastr [=] <string>

Output value for missing data, default ".".

=for Euclid:
	string.default: "."

=back

