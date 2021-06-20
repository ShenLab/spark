#!/usr/bin/env perl
use strict;
use warnings;
use FindBin qw|$Bin|;
use List::Util qw|sum|;
use List::MoreUtils qw|all|;
use Config::Std;
use IO::Dir;
use File::Temp qw|tempdir|;
use Hash::Util qw|lock_keys|;
use Getopt::Lucid qw|:all|;
use Bio::DB::HTS::Tabix;
use FaSlice;
use Genome::UCSC::TwoBit;
use Utils::List qw|all_combs|;
use Utils::File qw|count_line|;
use Utils::File::Iter qw|iter_file|;

use lib "$Bin/../lib";
use Shared qw|parse_filters slurp_xref expand_dat|;
use Variants qw|mutrater|;


############################
## Command line interface ##
############################

my @spec =  (
	Param("conf|c")->valid(sub { -r }),
	Param("input|in|i")->valid(sub { -r }),
	Param("wrkdir|wrk|w"),
	Param("output|out|o"),
	Switch("force"), # Under force mode, previously generated files will not be over-written
	Switch("help|h")
	);

my $opt = Getopt::Lucid->getopt(\@spec);

if ($opt->get_help) {
	#print STDERR "gene_mutrate.pl --conf CONFIG [--wrk WRKDIR] --in INPUT -out OUTPUT\n";
	print STDERR <<EOF;
Purpose:
	This is a standalone script to create gene size matrix used by dnenrich.

Notes:
	The script will have the same command line option as gene_mutrate and can re-use gene_mutrate's
	config file define mutation classes. The output is "transposed" gene size matrix such that 
	each column is a gene and each row represents a mutation class and sequence context combination. 
	The gene size matrix gives total number of possible counts for different classes of SNVs in
	differetn sequence context. These counts are used by dnenrich in permutation.

	One difference from the original implementation is we do not model overlapping genes, because such
	situation can be rare after removing read through genes, and also because functional effect are likely
	different for the same variants in overlapping genes.

EOF
	exit 1;
}

$opt->validate({ requires => [qw|conf input output|] });


my $infile = $opt->get_input;
my $outfile = $opt->get_output;
my $config = $opt->get_conf;

if (-f $outfile) {
	unless($opt->get_force) {
		print STDERR "Output file $outfile already exists!\n";
		exit 1;
	}
}

read_config $config => my %conf;
lock_keys(%conf);


if (-d $infile) {
	my @vartabs = grep { !/^\./ && -f "$infile/$_" } IO::Dir->new($infile)->read();
	# working directory must be provided when parallization
	my $wrkdir = $opt->get_wrkdir;
	if ($wrkdir) {
		mkdir $wrkdir unless -d $wrkdir;
	} else {
		$wrkdir = tempdir(CLEANUP => 1);
	}

	my $option = "";
	if (exists $conf{Parallel} && defined $conf{Parallel}{jobs}) {
		$option = join(" ", map { "--$_ $conf{Parallel}{$_}" } keys %{$conf{Parallel}});
	}
	else {
		die "Must provide number of threads in config file for parallel";
	}
	my $args = "";
	if ($opt->get_force) {
		$args = "--force";
	}

	my $command = qq(find $infile/ -not -path '*/\.*' -type f | parallel --eta $option "perl $0 --in {} --conf $config --out $wrkdir/{/}.txt $args 2>$wrkdir/{/}.err");
	#print $command, "\n"; exit 1;
	system($command);
	
	# Write the file list
	open my $flst, ">$wrkdir/.vartabs" or die "Cannot write to $wrkdir/.vartabs";
	foreach my $vartab (@vartabs) {
		print $flst $vartab, "\n";
	}
	close $flst;

	# For each result file, build an index for faster access
	my @index;
	my $nlines = count_line("$wrkdir/$vartabs[0].txt");
	for(my $ii = 0; $ii < @vartabs; $ii ++) {
		open my $fin, "$wrkdir/$vartabs[$ii].txt" or die "Cannot open $vartabs[$ii]";
		my $ll = 0;
		while(<$fin>) {
			my $offset = tell($fin);
			$index[$ii][$ll++] = $offset;
		}
		unless($ll = $nlines) {
			die "Result file $wrkdir/$vartabs[$ii].txt does not have correct number of lines!";
		}
	}

	# Merge and transpose the output to create final gene size matrix
	open my $fout, ">$outfile" or die "Cannot write to $outfile";
	for(my $ll = 0; $ll < $nlines; $ll ++) {
		for(my $ii = 0; $ii < @vartabs; $ii ++) {
			open my $fin, "$wrkdir/$vartabs[$ii].txt" or die "Cannot open $vartabs[$ii]";
			seek($fin, $index[$ii][$ll-1], 0) if $ll > 0;
			my $line = <$fin>; chomp($line);
			if ($ii == 0) {
				print $fout $line;
			}
			else {	
				my @dat = split(/\t/, $line); 
				die "Incorrect number of columns: $wrkdir/$vartabs[$ii].txt" unless @dat > 2;
				print $fout "\t", join("\t", @dat[2..$#dat]);
				if ($ii == @vartabs - 1) {
					print $fout "\n";	
				}
			}
		}
	}
	exit 0;	
}


# Mutation rate counter
$conf{MutRate}{Scale} = 1;
$conf{MutRate}{Count} = 1;
my $counter = mutrater($conf{MutRate});

# Parse filters 
my ($varclass, $callbacks, $fields) = parse_filters($conf{Variant}{Class}, $conf{Variant}{Class_Filter});

# Check fields in the input table file
my ($it, $fnames) = iter_file($infile);
foreach my $field (keys %$fields) {
	unless(grep { $_ eq $field } @$fnames) {
		die "Cannot find field $field in the input file";
	}
}
foreach my $field (qw|Chrom Position Ref Alt GeneID|) {
	unless(grep { $_ eq $field } @$fnames) {
		die "Cannot find required field $field in the input file";
	}
}

# Tally mutation counts
my (%mutct, %known);
while(my $dat = $it->()) {
	my $varid = join(":", @{$dat}{qw|Chrom Position Ref Alt|});
	next if defined $known{$varid};
	$known{$varid} = 1;
	my ($count, $context) = $counter->($dat);
	#next unless $context =~ /^[ACGT]+\>[ACGT]+$/;

	if ($dat->{GeneID} =~ /;/) {
		foreach my $data (expand_dat($dat, { sep => ';', optional => ['GeneID', keys %$fields] })) {
			for(my $ii = 0; $ii < @$varclass; $ii ++) {
				if ($callbacks->[$ii]->($data)) {
					$mutct{$data->{GeneID}}{$varclass->[$ii]}{$context} += $count if $context =~ /^[ACGT]+\>[ACGT]+$/;
					$mutct{$data->{GeneID}}{$varclass->[$ii]}{'*'} += $count;
				}
			}
		}
	}
	else {
		for(my $ii = 0; $ii < @$varclass; $ii ++) {
			if ($callbacks->[$ii]->($dat)) {
				$mutct{$dat->{GeneID}}{$varclass->[$ii]}{$context} += $count if $context =~ /^[ACGT]+\>[ACGT]+$/;
				$mutct{$dat->{GeneID}}{$varclass->[$ii]}{'*'} += $count;
			}
		}
	}
}


my @genes = sort keys %mutct;
my @context = ('*');
{
	my ($motiflen) = ($conf{MutRate}{Method} =~ /(\d+)mer/);
	my $adjlen = ($motiflen-1)/2;
	my @acgt = ([qw|A C G T|]) x $motiflen; 
	foreach my $from (all_combs(@acgt)) {
		my @to = @$from;
		foreach my $alt (qw|A C G T|) {
			next if $alt eq $from->[$adjlen];
			$to[$adjlen] = $alt; 
			push @context, join("", @$from).">".join("", @to);
		}
	}
}

open my $fout, ">$outfile" or die "Cannot write to $outfile";
print $fout join("\t", qw|Group	Type|, @genes), "\n";
foreach my $context (@context) {
	foreach my $class (@$varclass) {
		print $fout join("\t", '*', $context.":".$class,
			 	map { my $count = $mutct{$_}{$class}{$context};
			 		 defined $count ? int($count) : 0 } @genes), "\n";
	}
}


__END__