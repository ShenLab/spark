use strict;
use warnings;
use FindBin qw|$Bin|;
use Tie::IxHash;
use List::Util qw|first|;
use List::BinarySearch qw|binsearch_pos|;
use Utils::List qw|insert_after|;
use Utils::File::Iter qw|iter_file|;
use Getopt::Euclid;
use Config::Std;
use Data::Dumper;

use lib "$Bin/../../lib";
use Shared qw|parse_fstr expand_dat parse_filters parse_bins|;


# If conf is provided, will store class labels and class filters
# NewField => { Label => [qw|ClsA ClsB ...|], Calback => [qw|Expr1 Expr2 ...|] }
# We need to verify that NewField does not exist in the input and all filtering fields exist in the input
# When applied to the input, when multiple filters will pass, only the first label will be applied
# Variant that cannot be classified will have missing value. 
# New columns will also appear in the output
tie my %filter, 'Tie::IxHash';
my %filtfds;
if (defined $ARGV{'--conf'}) {
	my ($it2, $fnames2) = iter_file($ARGV{'--input'}, { fsep => qr/\t/ });
	read_config $ARGV{'--conf'} => my %conf;
	foreach my $newfd (keys %conf) {
		if (grep { $newfd eq $_ } @$fnames2) {
			die "New field in config $newfd already exists in the input!";
		}
		unless(defined $conf{$newfd}{Class} && defined $conf{$newfd}{Class_Filter}) {
			die "Must provide class label and filter for $newfd";
		}
		my ($class, $filters, $fields) = parse_filters($conf{$newfd}{Class}, $conf{$newfd}{Class_Filter});
		foreach my $fd (keys %$fields) {
			unless(grep { $fd eq $_ } @$fnames2) {
				die "Filter field $fd cannot be found in the input!";
			}
			$filtfds{$fd} = 1;
		}
		$filter{$newfd} = { Label => $class, Callback => $filters };
	}
}

# Determining the matching fields between input and weight table.
# @keyfd : categorical variables in the weight table
# @idfd: categorical variables in the input table
# $binfd: score bin field in the weight table
# scorefd: score field in the input table
# @wtfds: weight fields in the weight table
# @outfds: output weight fields in the output table.
my (@keyfds, @idfds, $binfd, $scorefd, @wtfds, @outfds);
{
	my $fields = parse_fstr($ARGV{'--weight-fields'}, 1);
	unless(scalar(keys %$fields) >= 3) {
		die "Should specify at least 3 columns from the lookup table: $ARGV{'--weight-fields'}";
	}
	my $realnum = qr/-?\d[\.0-9]*/;
	# Check that all specified fields can be found in from input/weight tables
	# newly defined columns from config will be added to input fields
	my ($it, $fnames) = iter_file($ARGV{'--weight'}, { fsep => qr/\t/ });
	my ($it2, $fnames2) = iter_file($ARGV{'--input'}, { fsep => qr/\t/ });
	my $dat = $it->();
	my $dat2 = $it->();
	foreach my $tabfd (keys %$fields) {
		unless(grep { $tabfd eq $_ } @$fnames) {
			die "Cannot find field $tabfd from weight: $ARGV{'--weight'}";
		}
		my $infd = $fields->{$tabfd};	
		# If score/bin field has been defined
		if (defined $scorefd && defined $binfd) {
			# Then check output field
			if (grep { $infd eq $_ } @$fnames2) {
				unless($ARGV{'--over-write'}) {
					die "The output field $infd already exists in the input: $ARGV{'--input'}!";
				}
				else {
					warn "The input file field $infd will be over-written";
				}
			}
			elsif (defined $filter{$infd}) {
				die "The output field $infd cannot be the same as one of derived fields!";
			}
			else {
				push @wtfds, $tabfd;
				push @outfds, $infd;
			}
		}
		else {
			unless( (grep { $infd eq $_ } @$fnames2) || defined $filter{$infd} ) {
				die "Cannot find $infd from the input file or as a derived field!";
			}
			# Determine if the field is the score field
			if (defined $ARGV{'--score'}) {
				if ($tabfd eq $ARGV{'--score'}) {
					$binfd = $tabfd;
					$scorefd = $infd;
					#if (defined $filter{$infd}) {
					#	die "The score field $scorefd cannot be one of the derived field";
					#}
				}
				else {
					push @keyfds, $tabfd;
					push @idfds,  $infd;
				}		
			}
			else {
				# If score field is not explcity specified, try to make a guess based o nformat
				if ($dat->{$tabfd} eq "*" || $dat->{$tabfd} =~ /^(\.,|${realnum}[,:])+$realnum$/) {
					unless (@keyfds) {
						die "Score field $tabfd:$infd must come after fields of categorical variables!";
					}
					$binfd = $tabfd;
					$scorefd = $infd;
				}
				else {
					push @keyfds, $tabfd;
					push @idfds,  $infd;
				}
			}
		}
	}
}

unless(@keyfds > 0 && @idfds > 0 && defined $binfd && defined $scorefd) {
	print Dumper \@keyfds, \@idfds, $binfd, $scorefd;
	die "Cannot find matching categorical variable and score fields!";
}
unless(@wtfds > 0 && @outfds > 0) {
	print Dumper \@keyfds, \@idfds, $binfd, $scorefd;
	die "Cannot find weight and output fields!";
}

printf STDERR "The following field(s) of catagorical variable from the input \"%s\" will be used to match field(s) %s in the weight table\n", 
	join(",", @idfds), join(",", @keyfds); 
printf STDERR "The score field %s from the input will be used to lookup range field \"%s\" in the weight table to determine weight\n", $scorefd, $binfd;


my ($it, $fnames) = iter_file($ARGV{'--weight'}, { fsep => qr/\t/ });
	
# Slurp lookup table, also store score bins
# Gene1 1:70:0.05,100 ...
# Then varsbin{Gene1} = "1:70:0.05,100", scorebin{"1:70:0.05,100"} = [1,2,...70,100]
# Missing data is allowed: Gene3 .,1:70:1 ...
# Missing data can only be part of the bin range, and appear at the beginning of bin
# Or use a catch-all bin *, so that any score including missing data will match: Gene2 * 100
# Weight for missing score or catch-all will be stored separately
my (%varsbin, %scorebin, %varwt, %varwt2);
while(my $dat = $it->()) {
	my $tabkey = join("\t", map { $dat->{$_} } @keyfds);
	if (defined $varwt{$tabkey}) {
		warn "Variant weights for [$tabkey] has already been found";
		next;
	}

	my $sbin = $dat->{$binfd};
	$varsbin{$tabkey} = $sbin;

	my @sbins;
	if ($sbin ne '*') {
		@sbins = parse_bins($sbin);
		unless (defined $scorebin{$sbin}) {
			if ($sbins[0] eq '.') {
				$scorebin{$sbin} = [@sbins[1..$#sbins]];	
			}
			else {
				$scorebin{$sbin} = \@sbins;			
			}
		}
	}

	foreach my $wtfd (@wtfds) {
		my @wts = split(',', $dat->{$wtfd});
		if ($sbin eq "*") {
			unless(@wts == 1) {
				die "Incorrect number of weights for $wtfd for * bin";
			}
			$varwt2{$tabkey}{$wtfd} = $wts[0];
		}
		else {
			unless(scalar(@wts) == scalar(@sbins)) {
				die "Incorrect number of weights for $wtfd: ".scalar(@wts)."!=".scalar(@sbins);
			}
			if ($sbins[0] eq ".") {
				$varwt2{$tabkey}{$wtfd} = $wts[0];
				$varwt{$tabkey}{$wtfd} = [@wts[1..$#wts]];
			}
			else {	
				$varwt{$tabkey}{$wtfd} = \@wts;
			}
		}
	}
}

# Annotate input file
my ($it2, $fnames2) = iter_file($ARGV{'--input'}, { fsep => qr/\t/ });

# The added column will be next to the original score field
my @fields = @$fnames2;
#my $f_score = first { $fields[$_] eq $scorefd } 0..$#fields;
#splice @fields, $f_score+1, 0, $outfd;
my @outfds2;
for(my $ii = 0; $ii < @wtfds; $ii ++) {
	my $wtfd = $wtfds[$ii];
	my $outfd = $outfds[$ii];
	unless (grep { $_ eq $outfd } @$fnames2) {
		push @outfds2, $outfd;
	}
}
# If config file is provided, we should further add additional fields defined in config
if (%filter) {
	push @fields, keys %filter;
	if (@outfds2 > 0) {
		push @fields, @outfds2;
	}
}
else {
	if (@outfds2 > 0) {
		insert_after(\@fields, $scorefd, @outfds2);
	}
}


my $fout;
if ($ARGV{'--output'}) {
	open $fout, ">$ARGV{'--output'}";
}
else {
	$fout = \*STDOUT;
}
print $fout join("\t", @fields), "\n";

my @expandfds = grep { !defined $filter{$_} } (@idfds, $scorefd);
if (%filtfds) {
	push @expandfds, sort keys %filtfds;
}

while(my $dat = $it2->()) {
	my %wts;
	foreach my $data (expand_dat($dat, { optional => \@expandfds, sep => ';'})) {
		# First derive custom fields if conf is provided
		if (%filter) {
			foreach my $newfd (keys %filter) {
				$data->{$newfd} = ".";
				for(my $ii = 0; $ii < @{$filter{$newfd}{Callback}}; $ii ++) {
					if ( $filter{$newfd}{Callback}[$ii]->($data) ) {
						$data->{$newfd} = $filter{$newfd}{Label}[$ii];
						last;
					}
				}
				push @{$wts{$newfd}}, $data->{$newfd};
			}
		}

		my $inkey = join("\t", map { $data->{$_} } @idfds); 
		my $score = $data->{$scorefd};
		foreach my $wtfd (@wtfds) {
			# Missing data if no categorical variables are found
			unless (defined $varsbin{$inkey}) {
				push @{$wts{$wtfd}}, $ARGV{'--nastr-out'};
				next;
			}

			my $sbin = $varsbin{$inkey};
			# Catch-all or missing value: put in the backup weight
			if ($sbin eq "*"  || $score eq $ARGV{'--nastr-in'}) {
				push @{$wts{$wtfd}}, $varwt2{$inkey}{$wtfd} // $ARGV{'--nastr-out'};
				next;
			}
			# Determine the numerical ranges in the score bin			
			unless(defined $scorebin{$sbin}) {
				die "Cannot find score bins for $sbin!";
			}
			my @sbins = @{$scorebin{$sbin}};

			# Now lookup the bin slot
			if ($score >= $sbins[-1] || $score < $sbins[0]) {
				if ($ARGV{'--strict'}) {
					die "Score for $inkey is outside range";
				}
				else {
					# If score is larger than the largest bin, use the last weight
					if ($score >= $sbins[-1]) {
						push @{$wts{$wtfd}}, $varwt{$inkey}{$wtfd}[-1];
					}
					else {
						push @{$wts{$wtfd}}, $varwt{$inkey}{$wtfd}[0];
					}
					next;
				}
			}
			my $index = binsearch_pos { $a <=> $b } $dat->{$scorefd}, @sbins;
			if ( $sbins[$index] == $score) {
				push @{$wts{$wtfd}}, $varwt{$inkey}{$wtfd}[$index];
			}
			else {
				unless($index > 0) {
					die "Incorrect zero index found for score $score";
				}
				if ($index > 0) {
					push @{$wts{$wtfd}}, $varwt{$inkey}{$wtfd}[$index-1];
				}
				else {
					push @{$wts{$wtfd}}, $varwt{$inkey}{$wtfd}[$index];
				}
			}
		}
	}
	for(my $ii = 0; $ii < @wtfds; $ii ++) {
		my $outfd = $outfds[$ii];
		my $wtfd = $wtfds[$ii];	
		$dat->{$outfd} = join(';', @{$wts{$wtfd}});
	}
	foreach my $newfd (keys %filter) {
		$dat->{$newfd} = join(";", @{$wts{$newfd}});
	}
	print $fout join("\t", @{$dat}{@fields}), "\n";
}



__END__

=head1 NAME

lookup_varwt.pl -- Lookup gene or category specific variant weights.

=head1 NOTES

We have annotated variant level information like conservation, predicted pathogenicity in the standard
annotation pipeline. But some variant level annotations are gene-specific. For example, Pop-score of 
PSAP method assign a significance level for every CADD score in different genes. It can be derived from 
CADD score and gene ID from existing annotation output using a gene-based lookup table. Another example
is empirically determined variant weight. One such weight is postive predictive value (PPV) that can be
defined for a group variant on some gene set based on the burden analysis. It has been used in DDD study
as weight for DNV enrichment test (DenovoWEST). In that study, they stratified genes based on the sHet, 
stratified variants based on CADD scores, and for missense variants based on presence or absence of missense
constrained regions. PPV was transformed from DNV fold enrichment in different categories, and can be
smoothed and used as weights.

In the examples above, we need to classify variant based on one or more categorical variable (like gene or
gene set, functional class, or presence or absence of MCR). Then we will need an original continuous score
for each variant. The lookup table should contain several columns of categorical variable(s) that matches
the corresponding field in the variant table, followed by a column that defines the end points to group original
scores into ranges (like Start:End:Step or V1,V2,V3... or Min,Start:End:Step,Max), and then one or more types 
of weights associated with different score ranges. The ranges are half-open intervals demarcated by the end points
[V1, V2), [V2, V3),... 
Lookup table should be tab separated, weights of one type should be packed in one column and comma separated.
An example is given in "Dropbox/Data/AnnoDB/PopScores/CADD13_gnomADALL.txt".

To deal with situations when some variants may not have a score. E.g., missense pathogenicity scores are only 
defined for missense variants, indels does not have CADD score. We also support two special types of ranges:
missing value "." and" a catch-all "*".  "." can only used as a part of end points list and must appear at the 
beginning. It matches the score that have missing value in the input. "*" can only used on its own, it matches
any score that cannot be found in the range including missing value. If "*" or "." is not used, missing values 
in the original score will have missing value in weights.

=head1 REQUIRED ARGUMENTS

=over

=item -[-]in[put] [=] <table>

Input variant table. It should be tab separated with necessary fields used to match variant weight lookup table.

=item -[-]w[eigh]t [=] <file>

Lookup table of variant weights.

=for Euclid:
	file.type: readable

=item -[-]w[eigh]t-fields [=] <string>

Define relevant fields in the weight lookup table and their matches to the input.

The first one or more columns should be catagorical fields used to match between lookup table and input,
then comes the column of score ranges specified in the lookup table and its continuous score counterpart in the input,
after that should be one or more columns of weights associated with score ranges and field names for output. 

By default, the script will determine the score ranges field from weight table by examining the content of first row.

Example: Gene:HGNC,CADD13,PopScore_Het:CADD13_PopScore
Gene,CADD13 are fields names from lookup table; HGNC,CADD13 are the corresponding fields in input.
CADD13 in the lookup table defines the score ranges, and in the input it is the original score.
PopScore_Het gives the weights for the score ranges defined in CADD13, CADD13_PopScore will be the additional 
column in the output giving weight for each variant.

=back

=head1 OPTIONS

=over

=item -[-]out[put] [=] <table>

Output table file name. If not provided, will dump to STDOUT. Columns will be tab separated.

=item -[-]score [=] <field>

Manually specify the score-range field in weight file, to avoid ambiguity.

=item -[-]conf[ig] [=] <config>

(Advanced use only) Provide config file to define new fields of input from existing columns.
For example, we can define variant classes based on GeneEff and CADD score. Section header will be
new column name. Each section should include multiple pairs of filter defintion and associated name: 
"Class=XXX" and "Class_Filter='Expr'".
The label for the first filter expression that is passed will be used as the value for the derived field.

=item -[-]strict

By default, if original score is outside of the specified range, minimal value will be assigned for score below
lower bound, and max value for score above upper bound. Under strict mode, program will exit with err.

=item -[-]over-write

By default, we require additional column of weight output should not exist in input.
Turn on this option to over-write an existing weight field in the input file.

=item -[-]nastr-in [=] <string>

NA string for the score from input. Weight table cannot have missing value.

=for Euclid:
	string.default: "."

=item -[-]nastr-out [=] <string>

NA string to appear in the output.

=for Euclid:
	string.default: "."

=back



