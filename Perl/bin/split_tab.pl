#!/usr/bin/env perl

use strict;
use warnings;
use POSIX qw|ceil|;
use Getopt::Euclid;
use Utils::File qw|count_line open_file|;


my $totlines = count_line($ARGV{'--input'});
unless($ARGV{'--no-header'}) {
	$totlines -= 1;
}

my ($nline, $nsplit);
if (defined $ARGV{'--nline'}) {
	$nsplit = int($totlines/$ARGV{'--nline'}+0.5);
	$nline = ceil($totlines/$nsplit);
}
elsif (defined $ARGV{'--nsplit'}) {
	$nsplit = $ARGV{'--nsplit'};
	$nline = ceil($totlines/$nsplit);
}
else {
	die "Should specify the number of lines per split or number of splits";
}

# Decide number of files with nline-1, and number with nline.
# Example: put 20 in 3 part; nline=7, nsplit=3; three parts are 7, 7, 6, nlarger == 2 (i.e. first two will have 7 lines) 
my $nlarger = $totlines % $nsplit;
$nlarger = $nsplit if $nlarger == 0;

if ($nlarger == $nsplit) {
	print STDERR "The input will be split into $nsplit files with $nline lines in each split\n";
}
else {
	print STDERR "The input will be split into $nsplit files with @{[ $nline-1 ]}~$nline lines in each split\n";
}
if ($nline < 2) {
	warn "Number of lines per file is too small!";
	exit 1;	
}

#my $ndigit = int(log($nsplit)/log(10)) + 1;
my $fin = open_file($ARGV{'--input'}) or die "Cannot open $ARGV{'--input'}";
my $header;
unless($ARGV{'--no-header'}) {
	$header = <$fin>;
}

my $ct = 0;
my $ii = $ARGV{'--start'};
my $fout;
while(<$fin>) {
	if ($ct % $nline == 0) {
		if ($ii == $nlarger + $ARGV{'--start'}) {
			$ct = 0;
			$nline --;
		}
		my $outfile = sprintf("%s.%d", $ARGV{'--output'}, $ii);
		if ($ARGV{'--suffix'}) {
			$outfile .= ".".$ARGV{'--suffix'};
		}
		open $fout, ">$outfile" or die "Cannot write to $outfile"; 
		print $fout $header if defined $header;
		$ii ++;
	}
	print $fout $_;
	$ct ++;
}

unless ($ii == $nsplit + $ARGV{'--start'}) {
	die "Incorrect number of splits: $ii!";
	for(my $jj = $ARGV{'--start'}; $jj <= $ii; $jj ++) {
		my $outfile = sprintf("%s.%d", $ARGV{'--output'}, $jj);
		$outfile .= ".".$ARGV{'--suffix'} if $ARGV{'--suffix'};
		unlink $outfile;
	}
}


__END__

=head1 NOTE

Split table file into multiple small er tables, each of which will have a header.

Output file will be "Prefix.N.Suffix" with N starting from 1 or otherwise a specified starting index.

=head1 REQUIRED ARGUMENTS

=over

=item -[-]in[put] [=] <file>

The input table file, gzipped files are supported.

=for Euclid:
	file.type: readable

=item -[-]out[put] [=] <prefix>

The output file prefix.

=back

=head1 OPTIONS

=over

=item -[-]start [=] <number>

The starting index for numbering split files.

=for Euclid:
	number.default: 1
	number.type: integer > 0

=item -[-]nline [=] <number>

Put fixed number of lines in each split, the real number of lines will differ
because adjustment will be made to ensure similar number of lines in each split.

=item -[-]nsplit [=] <number>

Split the table into N small files with similar number of lines.

=item -[-]suffix [=] <string>

Appending suffix to the file name.

=item -[-]no[-]header

Indicate that the input file has no header line.

=back

=cut
