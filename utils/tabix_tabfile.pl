use strict;
use warnings;
use Getopt::Euclid;
use Data::Dumper;
use List::Util qw|first|;
use List::MoreUtils qw|all|;
use Utils::File::Iter qw|iter_file|;
use FindBin qw|$Bin|;


my $args = index_tabfile($ARGV{'<table>'}, $ARGV{'--fields'});

my @fields = split(',', $ARGV{'--fields'});
unless(@fields == 2 || @fields == 3) {
	die "Incorrect number of positional fields";	
}

unless($ARGV{'<table>'} =~ /\.gz$/) {
	if ($ARGV{'--no-sort'}) {
		system("bgzip $ARGV{'<table>'}");
	}
	else {
		system("csvtk sort -t -k $fields[0] -k $fields[1]:n $ARGV{'<table>'} | bgzip -c > $ARGV{'<table>'}.gz");	
	}
	$ARGV{'<table>'} .= ".gz";
}

my $fopt = "";
if ($ARGV{'--force'}) {
	$fopt = "-f"
}

my $cmd = qq|tabix $args $fopt $ARGV{'<table>'}|;

print $cmd, "\n";
system($cmd);


sub index_tabfile {
	my ($table, $fieldstr) = @_;

	my @fields = split(',', $fieldstr);
	unless(@fields == 2 || @fields == 3) {
		die "Incorrect number of positional fields";
	}

	my $header = 1;
	if (all { /^\d+$/ } @fields) {
		$header = 0;
		#die "Must contain a file header";
	}

	my ($it, $fnames) = iter_file($table, { fsep => qr/\t/, header => $header, fastgz => 1 });

	my @index;
	for(my $ii = 0; $ii < @fields; $ii ++) {
		$index[$ii] = first { $fnames->[$_-1] eq $fields[$ii] } 1..@{$fnames};
		unless(defined $index[$ii]) {
			print Dumper $fnames;
			die "Cannot find column $fields[$ii] from input file";
		}
	}
	if (@index == 2) {
		$index[2] = $index[1];
	}

	my $args = "-s $index[0] -b $index[1] -e $index[2]";

	if ($header) {
		$args = "-S 1 $args";
		if ($ARGV{'--comment'}) {
			$args = "-c $fnames->[0] $args";
		}
	}
	return $args;
}


__END__

=head1 REQUIRED ARGUMENTS
 
=over

=item <table>

The input table file.

=head1 OPTIONS
 
=over

=item -[-]fields [=] <string>

The field names for Chrom,Pos or Chrom,Start,End

=for Euclid:
	string.default: "Chrom,Position"

=item -[-]c[omment]

Turn on tabix -c option to use first column name as comment.
This is useful if header is needed when using tabix for slicing.
Caution is to make sure that first column name will never appear in subsequent data rows.

=item -[-]no[-]sort

If file has been pre-sorted, turn on nosort swithc to skip sorting.

=item -[-]f[orce]

Force overwriting existing index.

=back