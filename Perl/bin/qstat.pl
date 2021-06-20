#!/usr/bin/env perl
use strict;
use warnings;
use Carp;
use List::Util qw|sum max|;
use File::Which;
use Data::Traverse qw|traverse|;

my $server = $ARGV[0] // "localhost";

# Determine job scheduling engine
my $qstat;
if ($server eq 'localhost') {
	$qstat = which('qstat');
}
else {
	$qstat = `ssh $server "which qstat"`;
}

my $engine;
unless ($qstat) {
	croak "Cannot find qstat on $server";
}
else {
	if ($qstat =~ /gridengine|sge/) {
		$engine = 'SGE';
	}
	else {
		croak "Unsupported engine: $qstat";
	}
}

my $qs;
if ($server eq 'localhost') {
	open $qs, qq{qstat -u '*' |} or croak "Cannot execute qstat on server";
}
else {
	open $qs, qq{ssh $server "qstat -u '*'" |} or croak "Cannot execute qstat on $server";
}
<$qs>; <$qs>;

my (%stats, %jobnum);
my (%cont);
while(<$qs>) {
	my ($jobid, $user, $status, $queue, $jaid) = (split)[0,3,4,7,9];
	my $qname;
	if ($queue) {
		($qname) = ($queue =~ /^(\S+)\@/);
		$qname = "" unless defined $qname;
	}
	else {
		$qname = "";
	}
	$stats{$user}{$qname}{$status} ++;
	push @{$cont{User}}, length($user);
	push @{$cont{Queue}}, length($qname);
	push @{$cont{Status}}, length($status);
	$jobnum{$user} ++;
}

#use Data::Dumper;
#print Dumper \%stats;

traverse { push @{$cont{Count}} => length($b) if /HASH/ } \%stats;
# Determine the appropriate format

my %fmt = (User => 8, Queue => 8, Status => 7, Count => 7);
foreach my $col (keys %fmt) {
	my $maxlen = max(@{$cont{$col}}) // 0 + 1;
	if ( $maxlen > $fmt{$col}) {
		$fmt{$col} = $maxlen;
	}
}

my @fields = qw|User Queue Status Count|;

print "|".join("|", map { '-' x $fmt{$_} } @fields), "|\n";
print "|".join("|", map { sprintf("%$fmt{$_}s", $_) } @fields), "|\n";
#print "|-------|--------|-------|-------|\n";
#print "|    User|   Queue| Status|  Count|\n";
foreach my $user (sort { $jobnum{$a} <=> $jobnum{$b} } keys %stats) {
	my $qinfo = $stats{$user};
	my $ntot = sum( map { scalar(keys %{$qinfo->{$_}}) } keys %$qinfo);
	my $noto = 0;
	print "|".join("|", map { '-' x $fmt{$_} } @fields), "|\n";
	while(my ($queue, $hashref) = each %$qinfo) {
		while(my ($status, $count) = each %$hashref) {
			$noto ++;
			if ($noto == 1) {
				print sprintf("|%$fmt{User}s", $user), sprintf("|%$fmt{Queue}s", $queue),
				sprintf("|%$fmt{Status}s|", $status), sprintf("%$fmt{Count}d|", $count), "\n";
			}
			else {
				print sprintf("|%$fmt{User}s", ""), sprintf("|%$fmt{Queue}s", $queue),
				sprintf("|%$fmt{Status}s|", $status), sprintf("%$fmt{Count}d|", $count), "\n";
			}
		}
	}
}
print "|".join("|", map { '-' x $fmt{$_} } @fields), "|\n";

