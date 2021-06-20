use strict;
use warnings;
use List::Util qw|sum|;
use File::Basename;
use FindBin qw|$Bin|;
use List::MoreUtils qw|notall|;
use Getopt::Euclid;
use Data::Dumper;
use Perl6::Slurp;
use Utils::Stat qw|mean|;
use Utils::Parser qw|sql_query|;
use Utils::File::Iter qw|iter_file|;

use lib "$Bin/../../lib";
use Shared qw|expand_dat|;

# Slurp gene expression values across selected tissues
my (%tx, %gx, @tissues);
{
	my ($it, $fnames) = iter_file($ARGV{'--exprmat'}, { fsep => qr/\t/ });	
	unless (@$fnames >= 3) {
		die "Expression matrix must have 3 or more columns!";
	}
	if (defined $ARGV{'--select'} || defined $ARGV{'--select-pattern'}) {
		if (defined $ARGV{'--select'}) {
			foreach my $tissue (split(',', $ARGV{'--select'})) {
				unless(grep { $tissue eq $_ } @$fnames) {
					die "Cannot find tissue $tissue in exprssion matrix";
				}
				push @tissues, $tissue;
			}
		}
		if (defined $ARGV{'--select-pattern'}) {
			my $regex = qr/$ARGV{'--select-pattern'}/;
			foreach my $tissue (grep { /$regex/ } @$fnames[2..$#$fnames]) {
				unless(grep { $tissue eq $_ } @tissues) {
					push @tissues, $tissue;
				}
			}
		}
		unless(@tissues) {
			die "No tissue is selected by --select/--select-pattern option";
		}
	}
	else {
		@tissues = @$fnames[2..$#$fnames];
	}
	# Now read tissue-specific transcript expression
	my %g2t;
	while(my $dat = $it->()) {
		my ($tid, $gid) = @{$dat}{@$fnames[0,1]};
		push @{$g2t{$gid}} => $tid;
		$tx{$tid} = [@{$dat}{@tissues}];
	}
	# And calculate gene expression by tanking the summation over all transcripts of a gene
	while(my ($gid, $trans) = each %g2t) {
		foreach my $tid (@$trans) {
			for(my $ii = 0; $ii < @tissues; $ii ++) {
				$gx{$gid}[$ii] += $tx{$tid}[$ii];
			} 
		}
	}
}


# Now processing input
my ($it, $fnames) = iter_file($ARGV{'--input'} eq '-' ? \*STDIN : $ARGV{'--input'}, {fsep => qr/\t/});
my @expandfds = qw|GeneID GeneEff TransIDs TransEffs|;
foreach my $field (@expandfds) {
	unless(grep { $field eq $_ } @$fnames) {
		die "Cannot find standard field $field from input!";
	}
}

my ($transfilt, @filtfds);
if (defined $ARGV{'--trans-filter'}) {
	($transfilt, my $tokens) = sql_query($ARGV{'--trans-filter'}, 1);
	foreach my $tok (@$tokens) {
		if ($tok->[0] eq 'FIELD') {
			unless(grep { $tok->[1] eq $_ } @$fnames) {
				die "The field $tok->[1] in filter expression cannot be found in the input";	
			}
			unless(grep { $tok->[1] eq $_ } @expandfds) {
				push @expandfds, $tok->[1];
			}
			unless(grep { $tok->[1] eq $_ } @filtfds) {
				push @filtfds, $tok->[1];
			}
		}
	}
}

my %genetrans;
if (defined $ARGV{'--gene-trans'}) {
	my %exclude;
	if (defined $ARGV{'--exclude-trans'}) {
		%exclude = map { (split)[0] => 1 } slurp $ARGV{'--exclude-trans'};
	}
	open my $fin, $ARGV{'--gene-trans'} or die "Cannot open gene-trans list";
	while(<$fin>) {
		my ($gid, $tid) = split;
		next if defined $exclude{$tid};
		push @{$genetrans{$gid}} => $tid;
	}
}
else {
	if (defined $ARGV{'--exclude-trans'}) {
		warn "--exclude-trans option will be ignored";
	}
}

my $outfd;
if ($ARGV{'--add-field'}) {
	$outfd = $ARGV{'--add-field'};	
}
else {
	my $fbase = basename($ARGV{'--exprmat'});
	$outfd = 'pExt_'.(split('\.', $fbase))[0];
}

my @outfields = @$fnames;
if (grep { $outfd eq $_ } @$fnames) {
	unless($ARGV{'--over-write'}) {
		die "Output field $outfd has been defined in the input!";
	}
	else {
		warn "Output field $outfd has been defined in the input!";
	}
}
else {
	push @outfields, $outfd;
}


my $fout;
if (defined $ARGV{'--output'}) {
	$fout = IO::File->new($ARGV{'--output'}, "w");
}
else {
	$fout = \*STDOUT;	
}
print $fout join("\t", @outfields), "\n";

my $ndigit = $ARGV{'--ndigit'};

while(my $dat = $it->()) {
	if ($dat->{TransIDs} eq '.') {
		$dat->{$outfd}  = ".";
		print $fout join("\t", @{$dat}{@outfields}), "\n";
		next;
	}
	my @pext;
	foreach my $data (expand_dat($dat, { fields => \@expandfds, sep => ';' })) {
		# Find out transcripts that are annotated the same as GeneEff
		my @transIDs = split(',', $data->{TransIDs});
		if (notall { defined $tx{$_} } @transIDs) {
			die "Not all transcript can be found in the expression matrix $ARGV{'--expr'}: $data->{TransIDs}";
		}
		my @transEffs = split(',', $data->{TransEffs});

		# select all transcripts whose TransEff equal GeneEff
		# If trans-filter is provided, further filter can be performed	
		my @transSub;
		if (defined $transfilt) {
			my @dat4filt;
			foreach my $field (@filtfds) {
				my @transInfo = split(',', $data->{$field});
				unless (@transInfo == @transIDs) {
					print STDERR "NTrans=", scalar(@transIDs), "\n", "$field=", $dat->{$field}, "\n";
					die "Filter field $field does not have the same length as TransIDs";
				}
				for (my $ii = 0; $ii < @transIDs; $ii ++) {
					$dat4filt[$ii]{$field} = $transInfo[$ii];
				}
			}
			@transSub = map { $transIDs[$_] } grep { $transfilt->($dat4filt[$_]) && $transEffs[$_] eq $data->{GeneEff} } 0..$#transIDs;
		}
		else {
			@transSub = map { $transIDs[$_] } grep {  $transEffs[$_] eq $data->{GeneEff} } 0..$#transIDs;
			unless(@transSub) {
				die "No transcript has the same annotated effect as $data->{GeneEff}";
			}
		}

		# Calculate normalized expression for affected tissues
		# Note the default is to use the sum of transcripts listed in TransID
		my @tissueNormExprs;
		if (@transSub) {
			my @transAll;
			if ($ARGV{'--gene-trans'}) {
				unless(defined $genetrans{$data->{GeneID}}) {
					warn "Cannot find transcripts for gene $dat->{GeneID}!";
					@transAll = @transIDs;
				}
				else {
					@transAll = @{$genetrans{$data->{GeneID}}};
				}
			}
			else {
				@transAll = @transIDs;
			}
			for(my $ii = 0; $ii < @tissues; $ii ++) {
				my $Esub = sum(map { $tx{$_}[$ii] } @transSub);
				my $Etot = sum(map { $tx{$_}[$ii] } @transAll);
				if ($Etot > 0) {
					push @tissueNormExprs, $Esub/$Etot;
				}
			}
		}
		if (@tissueNormExprs) {
			push @pext, sprintf("%.${ndigit}f", mean(@tissueNormExprs));
		}
		else {
			push @pext, $ARGV{'--nastr'};
		}
	}
	$dat->{$outfd} = join(";", @pext);
	print $fout join("\t", @{$dat}{@outfields}), "\n";
}



__END__

=head1 NAME

calc_pext.pl -- Calculate pext measure using tissue-specific expressions.

=head1 NOTE

The Pext metrics (Cummings B et al 2020) is defined for a coding variant as the
"[p]roportion of [ex]pressed [t]ranscripts" that share the same functional effects. 

Our implementation of pext takes annotated variant table and transcript-specific expressions 
as inputs, and calculate the proportion of total transcriptional output from a gene that 
are affected by the variant with the most severe functional effect (see lib/Variants.pm
for variant effect ranking)

We will make use of existing GeneID,GeneEff,TransIDs,TransEffs fields from the annotation
pipeline to find all transcript for each gene and transcript with the most severe functional
consequences caused by the variant. 

We need to make sure that gene and transcripts version used for annotation are the same
as for expression quantification.

=head1 REQUIRED ARGUMENTS
 
=over
 
=item -[-]in[put] [=] <table>

Input variant table. It should be tab separated and containing "GeneID,GeneEff,TransIDs" fields 
that are standard output from variant annotation pipeline. Use '-' if input is piped from STDIN.

=item -[-]expr[mat] [=] <table>

Gene expression matrix file. The first two columns should be transcript ID and gene ID, followed 
by expressions on different tissues. Some example gene expression matrix files can be found in 
`/home/x2680/Dropbox/Data/AnnoDB/GeneExpr/`

=back

=head1 OPTIONS

=over

=item -[-]select [=] <list>

Selecting list of tissue (comma separated) for calculating pext.

=item -[-]select-pattern [=] <regex>

Selecting tissues by regex pattern to include for calculation.

=item -[-]gene-trans [=] <file>

Full list of transcript for each gene, should have two columns: GeneID, TransID. We will use the sum of all 
transcripts that mapped to one gene as the denominator for normalized expression. If this file is not provided,
default is to use TransIDs field, but transcripts listed in TransID field may not be complete.

=item -[-]trans-filter [=] <filter>

Filter transcript. The field used for filter must have the same length as TransIDs.
E.g. for LoF variants, we can require LGD variants must have HC flag.

=item -[-]exclude-trans [=] <list>

A list of transcripts to be removed from gene-transcript table, in order to keep consistent with the
full transcript set used in annotation.

=item -[-]out[put] [=] <table>

Output table file name. If not provided, will it will be dumped to STDOUT. Fields will be tab separated.

=item -[-]add-field [=] <field>

The additional field name. Default will be "pExt_fbase", where fbase is the basename of expression
matrix file after stripping off suffices. It will appear by the end of input fields.

=item -[-]over-write

Over-write existing field.

=item -[-]ndigit [=] <number>

Number of effective digit for pExt values.

=for Euclid:
	number.default: 3

=item -[-]nastr [=] <string>

String for missing value in the output.

=for Euclid:
	string.default: "."


=back

