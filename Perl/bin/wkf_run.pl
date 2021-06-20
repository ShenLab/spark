#!/usr/bin/env perl

use strict;
use warnings;
use Carp;
use Cwd qw|abs_path|;
use Storable qw|retrieve|;
use Config::Std;
use IO::Prompt;
use Getopt::Euclid;
use Utils::Workflow;


# retrieve perl data and verify directory struvture.
croak "Cannot find .perldat under rootdir" unless -f "$ARGV{'--rootdir'}/.perldat";

#local $Storable::Eval = 1;
my $wkf = retrieve "$ARGV{'--rootdir'}/.perldat";

my %conf;
if (defined $ARGV{'--conf'}) {
	my $engine = $wkf->{engine};
	read_config $ARGV{'--conf'} => my %dat;
	croak "Cannot find $engine configuration" unless defined $dat{$engine};
	%conf = %{$dat{$engine}};
}
else {
	my $conf = retrieve "$ARGV{'--rootdir'}/.confdat";
	%conf = %$conf;
}

$ARGV{'--rootdir'} =~ s/\/$//;
croak "rootdir does not match!" unless abs_path($ARGV{'--rootdir'}) eq $wkf->{rootdir};
foreach my $subdir ($wkf->get_subdir) {
	croak "Cannot find subdir $subdir" unless -d $subdir;
}
my $srcdir = $wkf->get_subdir("src");
foreach my $taskname ($wkf->get_all_tasks) {
	croak "Cannot find script for $taskname!" unless -f "$srcdir/$taskname.sh";
}

# Determine the jobs to be executed.
my @all_tasks;
if (defined $ARGV{'--tasks'}) {
	push @all_tasks, split(/[,;]/, $ARGV{'--tasks'});
	foreach my $task (@all_tasks) {
		if (!defined $wkf->{tasks}{$task}) {
			croak "Cannot find task $task in the workflow!";
		}
	}
}
else {
	@all_tasks = $wkf->get_all_tasks;
}

my $logdir = $wkf->get_subdir("log");
if ($ARGV{'--force'}) {
	unless (prompt("By using --force, all previously generated data file will be erased, continue?",
		-yes_no, -default => 'n')) {
		exit 1;
	} 
	foreach my $taskname (@all_tasks) {
		foreach my $outfile ($wkf->get_expected($taskname)) {
			unlink $outfile if -f $outfile;
		}
	}
}

# Now run the workflow
$wkf->run({ conf => \%conf, tasks => \@all_tasks,
			nochain => $ARGV{'--nochain'}, dryrun => $ARGV{'--dryrun'} });


__END__

=head1 NAME

wkf_run -- Submit/execute tasks from existing workflows.

=head1 NOTE

Using default options with dryrun, this script can also be used
to check if all expected outputs of the workflow can be found.

Note: if customized callback is defined in workflow to check the output files
those custom checks will not be implemented by this workflow. 

=head1 REQUIRED ARGUMENTS

=over

=item -[-]root[dir] <dir>

The root directory with default layout created by the workflow manager.

=for Euclid:
    dir.type: readable

=back

=head1 OPTIONS

=over

=item -[-]conf[ig] <file>

Config file for the pipeline. Should be compatible with the format that can be read
by C<Config::Std> module. Engine specific configurations should be in a section
with the engine name.

=for Euclid:
  file.type: readable

=item -[-]task[s] <tasks>

A list of comma separated tasks to be executed.

=item -[-]force

By default, the execution engine is searching for unfinished tasks or slots
within job array by examinig expected outputs. This option would force execution
regardless of the expected outputs. The existing output files for the task to be 
executed will be removed. When the task list is given, only execute those specified tasks.

=item -[-]no[-]chain

By default, if one task is executed, all other tasks that depend on this one
will also be added to the queue list. This option turn off this chaining behavior.

=item -[-]dry[run]

Only print out the command without execution.

=back

=cut
