use strict;
use warnings;
use File::Basename qw|dirname|;
use File::Path qw|make_path|;
use File::Spec;
use Cwd qw|abs_path|;
use IO::Dir;
use String::ShellQuote;
use Getopt::Euclid;
use Utils::Dir qw|make_tree list_files|;

my $indir = $ARGV{'--src'};
my $outdir = $ARGV{'--dest'};

unless(-d $indir) {
	die "Cannot find input directory $indir";
}

if (-d $outdir) {
	if ($ARGV{'--force'}) {
		warn "Over-writing existing directory $ARGV{'--dest'}";
	}
	else {
		warn "The destination directory $ARGV{'--dest'} already exists, exit";
		exit 1;
	}
}

if ($ARGV{'--recursive'}) {
	make_tree($indir, $outdir) unless $ARGV{'--dry'};
}
else {
	make_path($outdir) unless $ARGV{'--dry'};
}

my %mirror;

my ($include, $exclude);
if ($ARGV{'--include'}) {
	$include = qr/$ARGV{'--include'}/;
}
if ($ARGV{'--exclude'}) {
	$exclude = qr/$ARGV{'--exclude'}/;
}
foreach my $file (list_files($indir, { type => "f", recur => $ARGV{'--recursive'},
		include => $include, exclude => $exclude })) {
	rsync_file($file);
}

sub rsync_file {
	my ($file) = @_;
	my $fbase = File::Spec->abs2rel($file, $indir);
	my $path;
	if (-l $file) {
		$path = find_linkpath($file, $ARGV{'--tracelink'});
		#if ($ARGV{'--tracelink'}) {
		#	while(1) {
		#		$file = find_linkpath($file);
		#		unless(-l $file) {
		#			$path = $file;
		#			last;
		#		}
		#	}
		#}
		#else {
		#	$path = find_linkpath($file);
		#}
	}
	else {
		$path = abs_path($file);
	}
	$path = shell_quote($path);
	my $outdirq = shell_quote($outdir);

	unless(-f $path || -l $path) {
		if ($ARGV{'--strict'}) {
			die "Path $path does not link to a real file or another symlink";
		}
		else {
			warn "Path $path does not link to a real file or another symlink";
		}
	}
	else {
		if ($ARGV{'--dedup'}) {
			if (defined $mirror{$path}) {
				my $relpath = File::Spec->abs2rel($mirror{$path}, dirname("$outdir/$fbase"));
				if ($ARGV{'--dry'}) {
					print qq|ln -s $relpath $outdir/$fbase\n|;
				}
				else {
					symlink $relpath, "$outdir/$fbase";
				}
			}
			else {
				if ($ARGV{'--dry'}) {
					print qq|rsync -azvP $path $outdirq/$fbase\n|;
				}
				else {
					system(qq|rsync -azvP $path $outdirq/$fbase|);
				}
				$mirror{$path} = "$outdir/$fbase";
			}
		}
		else {
			if ($ARGV{'--dry'}) {
				print qq|rsync -azvP $path $outdirq/$fbase\n|;
			}
			else {
				system(qq|rsync -azvP $path $outdirq/$fbase|);
			}
		}
	}
}



__END__

=head1 NAME

deep_mirror.pl -- Create deep copy of a directory.

=head1 USAGE

deep_mirror.pl --src SRCDIR --dest DESTDIR

=head1 NOTE

The script will do "deep mirror" of a source directory to the target directory.

Any symlinks in the source directory will be resolved, and if '--tracelink' option
is turned on, will trace link to the final original file.

By default, it only mirrors all files under the source directory. If '--recursive'
option is turned on, it will recursively do the same for all sub-directories.

If a symlink points to a directory, then the directory will be mirrored when '--recursive'
is on, otherwise it will be skipped.

We assume no symlinks are broken or circular. The former can be tested by dangle.pl utility.
But it is possible that multiple symlinks will point to the same file. In such cases,
by default the same file will be mirrored multiple times. To avoid duplicated files
occupying large disk spaces, there is an option '--dedup'. When turned on, we will
keep track of original files that have been mirrored, and create symlinks to the first 
one that have been mirrored.

To transfer data between servers, use sshfs to mount remote server as local FS.

=head1 REQUIRED ARGUMENTS
 
=over

=item -[-]src [=] <directory>

The source directory.

=item -[-]dest [=] <directory>

The destination directory. 

=back

=head1 OPTIONS
 
=over

=item -[-]trace[link]

Trace the symlink until it reaches to the final file or directory.

=item -[-]recur[sive]

Recursively mirror all sub-directories under the source directory.

=item -[-]include [=] <pattern>

=item -[-]exclude [=] <pattern>

File name patterns for inclusion or exclusion.

=item -[-]dedup

Avoid mirroring duplicated files.

=item -[-]dry

Print out the rsync commandlines without executing them.

=item -[-]strict

Under strict mode, script will die if any symlink is dangle.

=item -[-]force

Overwrite existing destination directory.

=back
