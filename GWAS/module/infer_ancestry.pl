#/usr/bin/env perl
use strict;
use warnings;
use Data::Dumper;
use List::Util qw|sum|;
use List::MoreUtils qw|all|;
use Math::Matrix;
use Utils::Stat::Running;
#use Utils::File::Iter qw|iter_file|;

my ($input, $output, $testpop) = @ARGV;

# First pass, determine number of PCs, number of reference populations
# Also calculated mean eigen vectors for each population among first numpc PCs
my ($numpc, @refpops, $numpop, %meanevec);
{
	my %stats;
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
		next if $label eq $testpop;
		for(my $ii = 1; $ii <= $numpc; $ii ++) {
			unless (defined $meanevec{$label}{$ii}) {
				$stats{$label}{$ii} = Utils::Stat::Running->new();
			}
			$stats{$label}{$ii}->push($dat[$ii-1]);
		}
	}
	@refpops = sort keys %meanevec;
	$numpop = scalar(@refpops);
	if ($numpc < scalar(@refpops)-1) {
		die "Number PCs should be at least the same as number of populations - 1";
	}
	foreach my $pop (@refpops) {
		$meanevec{$pop} = [map { $stats{$pop}{$_}->mean() } 1..($numpop-1)];
		push @{$meanevec{$pop}}, 1.0;
	}
}

# Calculate the transformation coefficient
my %coeff;
{	
	my $A = Math::Matrix->new(map { $meanevec{$_} } @refpops);
	foreach my $pop (@refpops) {
		my $x = Math::Matrix->new([map { $_ eq $pop ? 1 : 0 } @refpops]);
		my $E = $A->concat($x->transpose);
		my $s = $E->solve;
		$coeff{$pop} = [map { $_->[0] } @$s];
	}
}

if ($output =~ /\.txt$/) {
	$output =~ s/\.txt$//;
}

# Also output coeffs
open my $fpar, ">$output.coeffs.txt" or die "Cannot write to $output.coeffs.txt";
print $fpar join("\t", "Population", map { "C$_" } 1..($numpop-1)), "\n";
foreach my $pop (@refpops) {
	print $fpar join("\t", $pop, @{$coeff{$pop}}), "\n";
}

# Output predicted ancestry
open my $fout, ">$output.txt" or die "Cannot write to $output.txt";
print $fout join("\t", "Label", "FamID", "IID", (map { "Prob$_" } @refpops), (map { "PC_$_" } 1..$numpc)), "\n";
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
		my $label = pop @dat;
		unless(scalar(@dat) == $numpc) {
			die "Incorrect number of columns in $input";
		}
		#next unless $label eq $testpop;
		my @probs;
		foreach my $pop (@refpops) {
			my $prob = sum(map { $dat[$_-1]*$coeff{$pop}[$_-1] } 1..($numpop-1));
			$prob += $coeff{$pop}[$numpop-1]; 
			if ($prob < 0) {
				$prob = 0;
			}
			elsif ($prob > 100) {
				$prob = 100;
			}
			push @probs, $prob;
 		}
 		if (all { $_ == 0 } @probs) {
 			warn "All probes are 0: $sampid!";
 		}
 		else {
 			my $norm = sum(@probs);
 			for(my $ii = 0; $ii < @probs; $ii ++) {
 				$probs[$ii] = sprintf("%.3f", $probs[$ii]/$norm);
 			}
 		}
		print $fout join("\t", $label, $famid // $iid, $iid, @probs, @dat), "\n";
	}
}




