use strict;
use warnings;
use Sort::Versions;
use Utils::File qw|open_file|;
use Genome::Ranges::IntSet;

my ($infile, $cutoffs, $outprefix) = @ARGV;


foreach my $cutoff (split(',', $cutoffs)) {
	my $fin = open_file($infile);

	my $fout;
	if (defined $outprefix) {
		open $fout, ">${outprefix}_ge$cutoff.bed" or die "Cannot write to ${outprefix}_ge$cutoff.bed";
	}
	else {
		$fout = \*STDOUT;
	}

	my ($prev_flag, $chr, $pos, $score);
	my @runlist;
	while(<$fin>) {
		($chr, $pos, $score) = split;
		if ($score >= $cutoff) {
			if (@runlist > 0 && $pos != $runlist[-1] + 1) {
				print $fout join("\t", $chr, $runlist[0]-1, $runlist[-1]), "\n";
				@runlist = ();
			}
			push @runlist, $pos;
			$prev_flag = 1;
		}
		else {
			if ($prev_flag) {
				print $fout join("\t", $chr, $runlist[0]-1, $runlist[-1]), "\n";
				@runlist = ();
			}
			$prev_flag = 0;
		}
	}
	if (@runlist > 0) {
		print $fout join("\t", $chr, $runlist[0]-1, $runlist[-1]), "\n";
	}
}


