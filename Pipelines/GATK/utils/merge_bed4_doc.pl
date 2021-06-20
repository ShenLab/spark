#!/usr/bin/env perl

use strict;
use warnings;
use IO::Dir;
use FindBin qw|$Bin|;
use Perl6::Slurp;
use File::Temp qw|tempdir|;
use Getopt::Euclid;
use List::MoreUtils qw|uniq|;
use Utils::File qw|open_file count_line|;
use Utils::Stat qw|mean std|;

use lib "$Bin/../../lib";
use Shared qw|read_dir|;

my $files = read_dir($ARGV{'--indir'}, { suffix => $ARGV{'--suffix'}, 
										remove => $ARGV{'--remove'}, rename => $ARGV{'--rename'} });

my @samps;
if ($ARGV{'--keep'}) {
	my %keep = map { (split)[0] => 1 } slurp $ARGV{'--keep'};
 	@samps = grep { defined $keep{$_} } sort keys %$files;
}
else {
	@samps = sort keys %$files;
}
if (@samps > 0) {
	print STDERR "A total of ", scalar(@samps), " samples found\n";
}
else {
	print STDERR "No sample left, exit\n";
	exit 1;
	
}
if ($ARGV{'--batch'}) {
	if (-f $ARGV{'--output'}) {
		die "The output should be a directory under batch mode";
	}
	elsif (-d $ARGV{'--output'}) {
		warn "The output directory already exist, and content will be over-written";
	}
	else {
		mkdir $ARGV{'--output'};
	}
	open my $fin, "$ARGV{'--batch'}" or die "Cannot read batch file: $ARGV{'--batch'}";
	# first orgianize samples in batches, prepare one file per batch
	my $wrkdir = tempdir( CLEANUP => 1 );
	my %bybatch;
	while(<$fin>) {
		my ($iid, $batch) = split;
		if (grep { $iid eq $_ } @samps) {
			push @{$bybatch{$batch}} => $iid;
		}
	}
	open my $flst, ">$wrkdir/Batches.txt" or die "Cannot write to $wrkdir/Batches.txt";
	foreach my $batch (sort keys %bybatch) {
		print $flst $batch, "\n";
		open my $fout, ">$wrkdir/$batch.keep.txt" or die "Cannot write to $wrkdir/$batch.txt";
		print $fout join("\n", @{$bybatch{$batch}}), "\n";
	}
	my $option = "";
	foreach my $opt (qw|stat remove suffix|) {
		if ($ARGV{"--$opt"}) {
			if ($opt eq 'stat') {
				$option .= "--stat ";
			}
			else {
				$option .= "--$opt ".$ARGV{"--$opt"}." ";
			}
		}
	}
	my $cmd = "cat $wrkdir/Batches.txt | parallel --eta --jobs $ARGV{'--jobs'} ".
			  "'perl $0 --indir $ARGV{'--indir'} --output $ARGV{'--output'}/{}.bed.gz ".
			  "--keep $wrkdir/{}.keep.txt $option 2>$ARGV{'--output'}/{}.err'";
	print $cmd, "\n";
	system($cmd);
	exit 1;
}

# Create file handles, also check the number of lines are the same
my (%fin, @nlines);
foreach my $iid (keys %$files) {
	push @nlines => count_line($files->{$iid});
	$fin{$iid} = open_file($files->{$iid});
}
unless(scalar(uniq @nlines) == 1) {
	die "Not all input files have the same number of lines!";
}


my $fout;
my $outfile = $ARGV{'--output'};
unless ($outfile =~ /\.bed\.gz$/) {
	$outfile .= ".bed.gz";
}
open $fout, "| bgzip -c > $outfile" or die "Cannot write bgzip pipe";

if ($ARGV{'--stat'}) {
	print $fout join("\t", qw|Chrom Start End Mean SD|), "\n";
}
else {
	print $fout join("\t", qw|Chrom Start End|, @samps), "\n";
}

while(my $line = $fin{$samps[0]}->getline()) {
	chomp($line);
	my @data = (split(/\t/, $line))[0,1,2,3];
	my $origreg = join(":", @data[0,1,2]);
	for(my $ii = 1; $ii < @samps; $ii ++) {
		my $line2 = $fin{$samps[$ii]}->getline();
		chomp($line2);
		my ($chr, $start, $end, $depth) = (split(/\t/, $line2))[0,1,2,3];
		my $region = join(":", $chr, $start, $end);
		unless($region eq $origreg) {
			# When this happen, the previous output file should be remove
			close $fout;
			unlink $outfile;
			die "Region unmatched: $region ne $origreg";
		}
		push @data, $depth;
	}
	if ($ARGV{'--stat'}) {
		print $fout join("\t", @data[0,1,2], mean(@data[3..$#data]), std(@data[3..$#data])), "\n";
	}
	else {
		print $fout join("\t", @data), "\n";
	}
}


__END__

=head1 NAME

merge_bed4_doc.pl -- Merge DoC in BED4 format to a coverage matrix BED file.

=head1 REQUIRED ARGUMENTS
 
=over

=item -[-]indir [=] <indir>

Input directory of RD data. Prepare input directory to remove bad samples.

=item -[-]out[put] [=] <outfile>

Name for output file in bgzip'ed BED format. If file name does not end with .bed.gz, the suffix will be added. 
When batch file is provided, the output should be a directory.

=back

=head1 OPTIONS

=over

=item -[-]stat 

Output summary statistics instead of full DoC (mean and std).

=item -[-]remove [=] <list>

Samples to be removed.

=item -[-]rename [=] <list>

Rename samples. The provided list should be two column file of old and new ID.

=item -[-]keep [=] <list>

Samples to be included. Number of samples to be processed at once should not exceed the max.
number of possible open files.

=item -[-]suffix [=] <string>

File name suffix.

=for Euclid:
	string.default: 'cov.bed.gz'

=item -[-]batch [=] <file>

Batch file for parallel processing. Should be a list of two columns: IID and Batch.
When batch file is provided, batches will be processed in multiple jobs using GNU parallel.
Output should be a directory and file names will be Batch.bed.gz.

=item -[-]jobs [=] <number>

Number of parallel jobs.

=for Euclid:
	number.default: 4

=back




