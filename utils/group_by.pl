use strict;
use warnings;
use Getopt::Euclid;
use FindBin qw|$Bin|;
use Data::Dumper;
use File::Temp qw|tempdir|;
use File::Path qw|make_path|;
use List::Util qw|first|;
use List::MoreUtils qw|all|;
use Utils::Parser qw|sql_query|;
use Utils::File::Iter qw|iter_file|;

use lib "$Bin/../lib";
use Shared qw|parse_fstr parse_tabfile|;
use Variants qw|parse_operations|;


# First determine the column names for the grouping fields
my @gfields; # field names for grouping
my $wrkdir;
if ($ARGV{'--wrkdir'}) {
	unless(-d $ARGV{"--wrkdir"}) {
		mkdir $ARGV{"--wrkdir"};
	}	
	$wrkdir = $ARGV{'--wrkdir'};
}
else {
	$wrkdir = tempdir(CLEANUP => 1);
}
make_path "$wrkdir/tmp";
{
	my ($noheader, %opt);
	if ($ARGV{'--alias'}) {
		$opt{alias} = parse_fstr($ARGV{'--alias'});
		if (all { /^\d+$/ } keys %{$opt{alias}}) {
			$opt{header} = 0;
			$noheader = 1;
		}
	}
	my ($it, $fnames) = iter_file($ARGV{'--input'} eq '-' ? \*STDIN : $ARGV{'--input'}, { fsep => qr/\t/, %opt });

	#@gfields = split(',', $ARGV{'--group'});
	my $gfields = parse_gfields($ARGV{'--group'});
	@gfields = keys %$gfields;
	if (all { /^\d+$/ } @gfields) {
		die "Group fields can not be all numbers, must provide alias to input file without header!";
	}
	foreach my $gfield (@gfields) {
		unless(grep { $_ eq $gfield } @$fnames) {
			die "Group field $gfield cannot be found in input!";
		}
	}

	my $callback;
	if ($ARGV{'--filter'}) {
		my ($cb, $tokens) = sql_query($ARGV{'--filter'}, 1);
		foreach my $tok (@$tokens) {
			if ($tok->[0] eq 'FIELD') {
				unless(grep { $tok->[1] eq $_ } @$fnames) {
					die "The field $tok->[1] in filter expression cannot be found in the input";
				}
			}
		}
		$callback = $cb;
	}

	# Create sort arguments from grouping field or custom sort subrountine
	my $args = sort_args($gfields, $fnames);

	print STDERR "Writing sorted input to a temp file: sort $args\n";
	open my $fsort, "| body sort -T $wrkdir/tmp $args > $wrkdir/sorted_input.txt" or die "Cannopt open sort pipe";
	print $fsort join("\t", @$fnames), "\n";
	while(my $dat = $it->()) {
		if (defined $callback) {
			next unless $callback->($dat);
		}
		print $fsort join("\t", @{$dat}{@$fnames}), "\n";
	}
}

# Now open sorted file
my ($it, $fnames) = iter_file("$wrkdir/sorted_input.txt", { fsep => qr/\t/ });
# Parse operations
my $ops = parse_operations($ARGV{'--operation'}, { dbfields => $fnames });
# Need to check that operation names do not conflict with grouping fields
foreach my $opname (keys %$ops) {
	if (grep { $opname eq $_ } @gfields) {
		die "Operation column $opname conflict with grouping field: $ARGV{'--group'}";
	}
}

my $fout;
if ($ARGV{'--output'}) {
	open $fout, ">$ARGV{'--output'}" or die "Cannot write to $ARGV{'--output'}";
}
else {
	$fout = \*STDOUT;
}
print $fout join("\t", @gfields, keys %$ops), "\n";

# Write grouped summary
print STDERR "Iterating over sorted input and applying group operations\n";
my ($prevgrpid, @groupdat) ;
while(my $dat = $it->()) {
	my $grpid = join("\t", @{$dat}{@gfields});
	if (defined $prevgrpid && $prevgrpid ne $grpid) {
		# Apply operation
		my @outvals;
		foreach my $callback (values %$ops) {
			my $outval = $callback->(@groupdat);
			if (!defined $outval || $outval eq "") {
				$outval = ".";
			}
			push @outvals, $outval;
		}
		print $fout join("\t", $prevgrpid, @outvals), "\n";
		@groupdat = ();
	}
	push @groupdat, $dat;
	$prevgrpid = $grpid;
}
# last group
my @outvals;
foreach my $callback (values %$ops) {
	my $outval = $callback->(@groupdat);
	if (!defined $outval || $outval eq "") {
		$outval = ".";
	}
	push @outvals, $outval;
}
print $fout join("\t", $prevgrpid, @outvals), "\n";


# Make sort command line from specification
# Example: Chrom(V);Start(nr) => -k2,2V -k3,3nr
sub sort_args {
	my ($gfields, $fnames) = @_;
	my @args;
	while(my ($field, $type) = each %$gfields) {
		my $colnum = first { $fnames->[$_-1] eq $field } 1..@{$fnames};
		unless(defined $colnum) {
			die "Cannot find column number for $field";
		}
		push @args, "-k$colnum,${colnum}$type";
	}
	return join(" ", @args);
}

sub parse_gfields {
	my ($fstr) = @_;
	tie my %gfields, "Tie::IxHash";
	foreach my $conf (split(',', $fstr)) {
		my ($field, $type);
		if ($conf =~ /^(\w+)\((\w+)\)$/) {
			$field = $1;
			$type = $2;
		}
		elsif ($conf =~ /^(\w+)$/) {
			$field = $conf;
			$type = "";
		}
		else {
			die "Cannot recognize sort field: $conf";
		}
		if(defined $gfields{$field}) {
			die "Field $field appear multiple times!";
		}
		$gfields{$field} = $type;
	}
	return \%gfields;
}


__END__

=head1 NAME

group_by.pl -- “group by” clause for line-based database file.

=head1 NOTE

Given a input file and a set of "grouping fields", grou_by will compute the summary statistics
on other columns from the same file. This operation has been used as part of overlap_cnvs and
overlap_vars utilities. The script is a standalone implementation of group_by function to create 
summary stats files.

We will sort the input file by the grouping fields then apply the grouping functions on the fly.

=head1 REQUIRED ARGUMENTS

=over

=item -[-]in[put] [=] <file>

Input data table. It must be tab-separated with required fields. Use "-"" to stream from STDIN.

=item -[-]group [=] <fields>

The grouping fields in the input file, and optional their type in sorting.
Example: Chrom(V),Position(n),Ref,Alt

=item -[-]op[eration] [=] <fstr>

One or more operations and their alias as output columns.

=back

=head1 OPTIONS

=over

=item -[-]out[put] [=] <file>

Output file name.

=item -[-]wrk[dir] [=] <dir>

Temp directory for sorting large files.

=item -[-]alias [=] <fstr>

Rename input fields. If input has noheader, then it must be renamed. Otherwise alias is optional. 
Grouping and operation should also use renamed fields. 

=item -[-]filter [=] <expr>

A SQL filter expression that will be applied to input table.

=back

