#!/usr/bin/env perl 

#use warnings;
#use strict;
use Getopt::Std;
use File::Find;
use File::stat;
use User::pwent;
use User::grent;
use Utils::Number qw|pretty_bytes|;

=head1 NAME
 
lst - list sorted directory contents (depth first)

=head1 DESCRIPTION

Have you ever wondered what the newest or biggest files within a directory are? 
The standard ls program has options for listing out directories sorted in time order 
(the -t flag) and for recursing into subdirectories (the -R flag). However, it pauses
at each directory to display the sorted contents of just that directory.  It doesn't
descend through all subdirectories first and then sort everything it found.
The lst program does just that.

=cut
  
getopts("lusrcmi")                            or die <<DEATH;
Usage: $0 [-mucsril] [dirs ...]
or    $0 -i [-mucsrl] < filelist

Input format:
-i  read pathnames from stdin
Output format:
-l  long listing
Sort on:
-m  use mtime (modify time) [DEFAULT]
-u  use atime (access time)
-c  use ctime (inode change time)
-s  use size for sorting
Ordering:
-r  reverse sort
NB: You may only use select one sorting option at a time.
DEATH

unless ($opt_i || @ARGV) { @ARGV = (".") }

if ($opt_c + $opt_u + $opt_s + $opt_m > 1) {
	die "can only sort on one time or size";
}

$IDX = "mtime";
$IDX = "atime" if $opt_u;
$IDX = "ctime" if $opt_c;
$IDX = "size"  if $opt_s;

$TIME_IDX = $opt_s ? "mtime" : $IDX;

*name = *File::Find::name;  # forcibly import that variable

# the $opt_i flag tricks wanted into taking
# its filenames from ARGV instead of being
# called from find.

if ($opt_i) {
	*name = *_;  # $name now alias for $_
    while (<>) { 
    	chomp; &wanted; 
    }   # ok, not stdin really
}
else {
	find(\&wanted, @ARGV);
}

# sort the files by their cached times, youngest first
@skeys = sort { $time{$b} <=> $time{$a} } keys %time;
  
# but flip the order if -r was supplied on command line
@skeys = reverse @skeys if $opt_r;
  
for (@skeys) {
	unless ($opt_l) {  # emulate ls -l, except for permissions
		print "$_\n";
		next;
	}
	$now = localtime $stat{$_}->$TIME_IDX( );
	printf "%04o\t%3d\t%6s\t%6s\t%s\t%s\t%s\t%s\n",
		$stat{$_}->mode( ) & 07777,
		$stat{$_}->nlink( ),
		user($stat{$_}->uid( )), group($stat{$_}->gid( )),
		pretty_bytes($stat{$_}->size( )),
		$now, $_;
}

# get stat info on the file, saving the desired
# sort criterion (mtime, atime, ctime, or size)
# in the %time hash indexed by filename.
# if they want a long list, we have to save the
# entire stat object in %stat.  yes, this is a
# hash of objects
sub wanted {
	my $sb = stat($_);  # XXX: should be stat or lstat?
	return unless $sb;
	$time{$name} = $sb->$IDX( );  # indirect method call
    $stat{$name} = $sb if $opt_l;
}
  
# cache user number to name conversions; don't worry
# about the apparently extra call, as the system caches the
# last one called all by itself
sub user {
	my $uid = shift;
  	$user{$uid} = getpwuid($uid) ? getpwuid($uid)->name : "#$uid"
  	unless defined $user{$uid};
  	return $user{$uid};
}
  
# cache group number to name conversions; ditto on unworryness
sub group {
  	my $gid = shift;
  	$group{$gid} = getgrgid($gid) ? getgrgid($gid)->name : "#$gid"
  	unless defined $group{$gid};
  	return $group{$gid};
}


