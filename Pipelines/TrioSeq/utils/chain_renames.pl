use strict;
use warnings;
use Carp;
use List::MoreUtils qw|all|;
use FindBin qw|$Bin|;

use lib "$Bin/../../lib";
use Shared qw|chain_renames|;


unless (@ARGV) {
	print STDERR "Usage:\n";
	print STDERR "$0 RENAME1 RENAME2 .... > RENAME_ALL.txt\n";
	print STDERR "Rename files should be listed from older to newer ones\n";
	exit 1;
}
else {
	unless (all { -f $_ } @ARGV) {
		die "Cannot find all input files";
	}
}

my @idpairs = chain_renames(@ARGV);

for my $pair (@idpairs) {
	print join("\t", @$pair), "\n";
}
