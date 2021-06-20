#!/usr/bin/env perl

use strict;
use warnings;
use File::Basename;
use File::Temp qw|tempdir|;
use File::Path qw|remove_tree|;
use File::Copy qw|copy|;
use Utils::File qw|find_linkpath|;
use List::MoreUtils qw|uniq|;
use Getopt::Euclid;

my @inputs = split($ARGV{'--sep'}, $ARGV{'--input'});

# Assuming input file have different base names, so they will not crash when copied to the same dir
my @inbase = map { basename($_) } @inputs;
my @inbaseuq = uniq sort @inbase;
unless(@inbase == @inbaseuq) {
	die "Input files do not have unique basenames!";
}

my $tmpdir = tempdir(CLEANUP => 0, DIR => $ARGV{'--tmpdir'});

my $flag;
foreach my $infile (@inputs) {
	# AWS s3
	if ($infile =~ /^s3:/) {
		my $aws = "aws s3";
		if ($ARGV{'--profile'}) {
			$aws .= " --profile $ARGV{'--profile'}";
		}
		print STDERR "Copy $infile to $tmpdir\n";
		my $code = system(qq|$aws cp $infile $tmpdir --no-progress > /dev/null|);
		die "Copy $infile failed" if $code;
		if ($ARGV{'--all'}) {
			my $indir = dirname($infile);
			my $fbase = basename($infile);
			print STDERR "Copy $infile.* to $tmpdir, if any\n";
			my $code = system(qq|$aws cp $indir $tmpdir --recursive --exclude '*' --include '$fbase.*' --no-progress >/dev/null|);
			die "Copy $infile.* failed" if $code;
		}
	}
	else {
		# ordinary files
		if (-l $infile) {
			my $path = find_linkpath($infile, 1);
			print STDERR "Trace the link: $infile => $path\n";
			my $fbase = basename($infile);
			print STDERR "Copy $path to $tmpdir/$fbase\n";
			copy($path, "$tmpdir/$fbase") or die "Copy $path failed!";
		}
		else {
			print STDERR "Copy $infile to $tmpdir\n";
			copy($infile, $tmpdir) or die "Copy $infile failed!";
		}
		if ($ARGV{'--all'}) {
			my @extra = glob("$infile.*");
			print STDERR "Extra files found: ", join(" ", @extra), "\n";
			foreach my $xtrafile (@extra) {
				if (-l $infile) {
					my $path = find_linkpath($xtrafile, 1);
					print STDERR "Trace the link: $xtrafile => $path";
					my $fbase = basename($xtrafile);
					print STDERR "Copy $path to $tmpdir/$fbase\n";
				}
				else {
					print STDERR "Copy $xtrafile to $tmpdir\n";
					copy($xtrafile, $tmpdir) or die "Copy $xtrafile to $tmpdir";
				}
			}
		}
	}
}
# If all files are copied correctly
$flag = 1;

# Print out new file paths
# Can only reach this point when all input files are copied correctly
my $outlist = join($ARGV{'--sep'}, map { $tmpdir."/".$_ } @inbase);
if ($ARGV{'--output'}) {
	open my $fout, ">$ARGV{'--output'}" or die "Cannot write to $ARGV{'--output'}";
	print $fout $outlist, "\n";
	print $tmpdir, "\n";
}
else {
	print $outlist, "\n";
}


# End block will run even some prior step failed
END {
	unless($flag) {
		print STDERR "Remove $tmpdir\n";
		remove_tree($tmpdir);
	}
}

__END__

=head1 NAME

stagein.pl -- Stage in files from remote directory

=head1 NOTE

A temp dir will be created and for each input files, the utility will try to copy the file to the
tempdir. If all copy operations are successful, it will return the list of copied location. 
We will not check if input file exists on the remote dir, because otherwise the copy operation
will fail. If input is symlink, we should follow the link to find actual files. For file system
mounted by sshfs, this is automatically done by sshfs.

Update 202011: now adding support for AWS s3

=head1 REQUIRED ARGUMENTS

=over

=item -[-]in[put] [=] <files>

One or more input files to stage in.
The input must be files not directories! But we do not check input.

=back
 
=head1 OPTIONS

=item -[-]out[put] [=] <outfile>

Write the output file names to a file.

=item -[-]tmp[dir] [=] <dir>

Base directory for creating temp directory.

=for Euclid:
	dir.default: "/tmp"

=item -[-]profile [=] <name>

For AWS s3, the profile name of the credential to access the bucket.

=item -[-]all

Also stage in additional files that have input file matching the name as prefix.

=item -[-]sep [=] <char>

Input file separator.

=for Euclid:
	char.default: ","

=back






