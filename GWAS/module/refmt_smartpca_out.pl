#/usr/bin/env perl
use strict;
use warnings;
use Data::Dumper;
use List::Util qw|sum|;
use List::MoreUtils qw|all|;
use Math::Matrix;
use Utils::Stat::Running;

my ($input, $output) = @ARGV;

# First pass, determine number of PCs, number of reference populations
# Also calculated mean eigen vectors for each population among first numpc PCs
my ($numpc, @refpops);
{
	my %known;
	open my $fin, $input or die "Cannot open $input";
	# Skip header
	<$fin>;
	while(<$fin>) {
		my @dat = split;
		unless(defined $numpc) {
			$numpc = scalar(@dat)-2;
		}
		my $sampid = shift @dat;
		my $label = pop @dat;
		$known{$label} = 1;
	}
	@refpops = sort keys %known;
	if ($numpc < scalar(@refpops)) {
		die "Number PCs should be at least the same as number of populations";
	}
}


open my $fout, ">$output" or die "Cannot write to $output";
print $fout join("\t", "Label", "FamID", "IID", (map { "PC_$_" } 1..$numpc)), "\n";
{
	open my $fin, $input or die "Cannot re-open $input"; <$fin>;
	while(<$fin>) {
		my @dat = split;
		my $sampid = shift @dat;
		my @idparts = split(':', $sampid);
		my ($famid, $iid);
		if (@idparts == 4) {
			$famid = $idparts[1];
			$iid = $idparts[3];
		}
		elsif (@idparts == 2) {
			$iid = $idparts[1];
		}
		else {
			die "Cannot recognize composite sample ID: $sampid";
		}
		#my ($famid, $iid) = (split(':', $sampid))[1,3];
		my $label = pop @dat;
		print $fout join("\t", $label, $famid // $iid, $iid, @dat), "\n";
	}
}




