#!/usr/bin/env perl 


=head1 NAME

dangle - detect dangle symbolic links.

=head1 DESCRIPTION

Recursively walk through a directory tree, reporting all symlinks that are dangle.
 * The script can handle the relative path created for the symlink.
 * But It does not trace recursively to links to links.

=head1 USAGE

dangle SRC_DIR ...

=cut

use warnings;
use strict;
use File::Basename qw|dirname|;
use Utils::Dir qw|dir_walk|;

unless(@ARGV >= 1) {
	print STDERR "Usage: $ARGV[0] SRC_DIR ...\n";
  exit 1;
}

sub dangles {
  my $file = shift;
  if ( -l $file ) {
    my $target = readlink $file;
    if ($target =~ /^\//) {
      if ( ! -e $target && ! -l $target ) {
        print "$file -> $target broken\n";
      }
    }
    else {
      my $srcdir = dirname($file);
      my $tgfile = "$srcdir/$target";
      if ( ! -e $tgfile && ! -l $tgfile ) {
        print "$file -> $tgfile broken\n";
      }
    }
  }
}

foreach my $basedir (@ARGV) {
	die "$basedir is not a directory" unless -d $basedir;
	dir_walk($basedir, \&dangles);
}

