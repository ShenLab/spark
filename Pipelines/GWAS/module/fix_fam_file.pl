use strict;
use warnings;
use Perl6::Slurp;
use List::MoreUtils qw|uniq|;
use File::Copy;

my ($input, $output) = @ARGV;

# Fix pedigree file after sex update:
# 1. We will zero out parents with unreasonable sex
# 2. We will add back sex code for parents with unknown sex
unless(@ARGV) {
	print STDERR "Fix pedigree file due to incompatible sex issues.\n";
	print STDERR "$0 INPUT [OUTPUT]\n";
}
else {
	unless (-f $input) {
		die "Cannot find input file $input";
	}
	unless (defined $output) {
		print STDERR "Replacing output to $input\n";
		copy($input, $input.'.bak');
		$output = $input;
		$input = $input.'.bak';
	}
}

my %sex = map { my ($fid, $iid, $sex) = (split)[0,1,4]; ($fid.$;.$iid => $sex) } slurp $input;

# identify samples require sex update update
# 1. update from unknown to male/female
# 2. update from male to female or female to male due to sample swap
my %update;
{
	open my $fin, $input or die "Cannot open $input for reading";	
	while(<$fin>) {
		my ($fid, $iid, $dad, $mom, $gender, $phe) = split;
		if ($dad ne '0' && defined $sex{$fid,$dad} && 
			$sex{$fid,$dad} ne '1' && $sex{$fid,$dad} ne '2') {
			print STDERR "Update sex for dad: $fid $dad\n";
			$update{$fid,$dad} = 1;
			$sex{$fid,$dad} = 1;
		}
		if ($mom ne '0' && defined $sex{$fid,$mom} && 
			$sex{$fid,$mom} ne '1' && $sex{$fid,$mom} ne '2') {
			print STDERR "Update sex for mom: $fid $mom\n";
			$update{$fid,$mom} = 2;
			$sex{$fid,$mom} = 2;
		}
	}
}

my %rename;
open my $fin, $input or die "Cannot open $input for reading";
open my $fout, ">$output" or die "Cannot write to $output";
while(<$fin>) {
	my ($fid, $iid, $dad, $mom, $gender, $phe) = split;
	if (defined $update{$fid,$iid}) {
		$gender = $update{$fid,$iid};
	}
	if ($dad ne '0' && $mom ne '0' && 
		defined $sex{$fid,$dad} && defined $sex{$fid,$mom} && 
		$sex{$fid,$dad} eq '2' && $sex{$fid,$mom} eq '1') {
		print STDERR "Parents swapped: $fid $dad <--> $mom\n";
		# parents swapped
		$rename{$fid."\t".$dad} = $fid."\t".$mom;
		$rename{$fid."\t".$mom} = $fid."\t".$dad;
	}
	else {
		if ($dad ne '0' && defined $sex{$fid,$dad} && $sex{$fid,$dad} eq '2') {
			print STDERR "Zero out dad: $fid $dad\n";
			$dad = 0;
		}
		if ($mom ne '0' && defined $sex{$fid,$mom} && $sex{$fid,$mom} eq '1') {
			print STDERR "Zero out mom: $fid $mom\n";
			$mom = 0;
		}
	}
	print $fout join("\t", $fid, $iid, $dad, $mom, $gender, $phe), "\n";
}

if (keys %rename) {
	open my $frn, ">$output.rename" or die "Cannot write rename list";
	while(my ($fiid1, $fiid2) = each %rename) {
		print $frn $fiid1, "\t", $fiid2, "\n";
	}
}
