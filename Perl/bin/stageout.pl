#!/usr/bin/env perl

use strict;
use warnings;
use File::Basename;
use File::Copy qw|copy|;
use File::Path qw|make_path|;
use List::MoreUtils qw|all uniq|;
use Cwd qw|abs_path|;
use Getopt::Euclid;

my @inputs = split($ARGV{'--sep'}, $ARGV{'--files'});
if ($ARGV{'--indir'}) {
	@inputs = map { $ARGV{'--indir'}."/".$_ } @inputs;
}

unless(all { -f $_ } @inputs) {
	print join("\n", @inputs), "\n";
	die "Not all inputs are files!"
}

# Assuming input file have different base names, so they will not crash when copied to the same dir
my @inbase = map { basename($_) } @inputs;
my @inbaseuq = uniq sort @inbase;
unless(@inbase == @inbaseuq) {
	die "Input files do not have unique basenames!";
}

my $s3flag;
my $outdir = $ARGV{'--outdir'};

if ($outdir =~ /^s3:/) {
	$s3flag = 1;
}
else {
	unless(-d $outdir) {
		make_path $outdir;
	}
}

foreach my $infile (@inputs) {
	print STDERR "Copying $infile to $outdir\n";
	if ($s3flag) {
		my $aws = "aws s3";
		if ($ARGV{'--profile'}) {
			$aws .= " --profile $ARGV{'--profile'}";
		}
		my $fbase = basename($infile);
		print STDERR "Copy $infile to $outdir\n";
		my $code = system(qq|$aws cp $infile $outdir/$fbase|);
		die "Copy $infile failed" if $code;		
	}
	else {
		print STDERR "Copy $infile to $outdir\n";
		copy($infile, $outdir) or die "Copy $infile failed!";
	}
}

# Replace input files with symlinks to remote files
# Note for s3 we do not create symlinks
unless($s3flag) {
	foreach my $infile (@inputs) {
		print STDERR "Replace $infile with symlink\n";
		my $inbase = basename($infile);
		my $target = abs_path("$outdir/$inbase");
		unlink($infile);
		symlink($target, $infile);
	}
}


__END__

=head1 NAME

stageout.pl -- Stage out files to remote directory

=head1 NOTE

After files are staged out, original files will be removed and replaced by a symlink.
Directory in remote directory will be created if they do not exist.

Update 202011: now adding support for AWS s3.
We will not create symlinks after stage out to S3.


=head1 REQUIRED ARGUMENTS

=over

=item -[-]files [=] <files>

One or more input files to stage out.
The input must be files not directories!

=item -[-]out[dir] [=] <outdir>

The output root directory. 

=back
 
=head1 OPTIONS

=item -[-]profile [=] <name>

For AWS s3, the profile name of the credential to access the bucket.

=item -[-]indir [=] <indir>

Input directory, such that all input file paths are relative to this directory.

=item -[-]sep [=] <char>

Input file separator.

=for Euclid:
	char.default: ","

=back







