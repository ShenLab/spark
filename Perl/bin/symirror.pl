#!/usr/bin/env perl 


=head1 NAME

symirror - build spectral forest of symlinks

=head1 DESCRIPTION

Recursively duplicates a directory tree, making a shadow forest full of 
symlinks pointing back at the real files.

=head1 USAGE

symirror SRC_DIR DST_DIR

=cut

use warnings;
use strict;
use Pod::Usage;
use Cwd         qw(realpath);
use File::Find  qw(find);

unless(@ARGV == 2) {
    pod2usage(1);
}

# die "usage: $0 realdir mirrordir" unless @ARGV == 2;

my $SRC = realpath $ARGV[0];
my $DST = realpath $ARGV[1];

my $oldmask = umask 077;        # in case was insanely uncreatable
chdir $SRC                      or die "can't chdir $SRC: $!";
unless (-d $DST) {
    mkdir($DST, 0700)           or die "can't mkdir $DST: $!";
}
find { 
    wanted      => \&shadow,
    postprocess => \&fixmode,
} => ".";
umask $oldmask;

sub shadow {
    (my $name = $File::Find::name) =~ s!^\./!!;     # correct name
    return if $name eq ".";
    if (-d) { # make a real dir; we'll copy mode later
        mkdir("$DST/$name", 0700)       
                        or die "can't mkdir $DST/$name: $!";
    } else {  # all else gets symlinked                      
        symlink("$SRC/$name", "$DST/$name")
                        or die "can't symlink $SRC/$name to $DST/$name: $!";
    }
}

sub fixmode {
    my $dir = $File::Find::dir;
    my $mode = (stat("$SRC/$dir"))[2] & 07777;
    chmod($mode, "$DST/$dir") 
        or die "can't set mode on $DST/$dir: $!";
}

