use strict;
use warnings;
use Perl6::Slurp;
use List::MoreUtils qw|all|;
use List::Util qw|first|;
use Utils::File::Iter qw|iter_file|;
use Getopt::Euclid;


# Merge variant tables that split into chunks and sort the output by Chrom,Position
# We use GNU sort utility to sort large files.
my $inprefix = $ARGV{'--input'};
my $outdir = $ARGV{'--outdir'};
my $tmpdir = $ARGV{'--tmpdir'};
my $suffix = $ARGV{'--suffix'};

mkdir $outdir unless -f $outdir;

# For multiple splits within each group, field names should be the same
my @groups;
if (-f $ARGV{'--group'}) {
 	@groups = map { [split] } slurp $ARGV{'--group'};
 	foreach my $group (@groups) {
 		unless(all { -f "$inprefix.$_.$suffix"  } $group->[1]..$group->[2]) {
			die "Not all input files $inprefix.{$group->[1]..$group->[2]}.$suffix can be found!";
		}
 	}
}
else {
	if (defined $ARGV{'--nsplit'}) {
		@groups = ([$ARGV{'--group'}, 1, $ARGV{'--nsplit'}]);
		unless(all { -f "$inprefix.$_.$suffix" } 1..$ARGV{'--nsplit'}) {
			die "Not all input files $inprefix.{1..$ARGV{'--nsplit'}}.$suffix can be found!";
		}
	}
	else {
		my $ct = 1;
		while(1) {
			if (-f "$inprefix.$ct.$suffix") {
				$ct ++;
			}
			else {
				last;
			}
		}
		if ($ct <= 1) {
			die "Cannot find input file $inprefix.1.$suffix!";
		}
		else {
			@groups = ([$ARGV{'--group'}, 1, $ct-1]);
		}
	}
}

foreach my $dat (@groups) {
	my ($label, $start, $end) = @$dat;
	my (@fields, $fout, $f_chr, $f_pos);
	foreach my $ii ($start..$end) {
		next if -z "$inprefix.$ii.$suffix";
		my ($it, $fnames) = iter_file("$inprefix.$ii.$suffix", { fsep => qr/\t/ });
		unless(@fields) {
			@fields = @$fnames;
			$f_chr = first { $fields[$_-1] eq 'Chrom' } 1..scalar(@fields);
			$f_pos = first { $fields[$_-1] eq 'Position' } 1..scalar(@fields);
			unless(defined $f_chr && defined $f_pos) {
				die "Cannot find column Chrom and Positon";
			}
			if ($ARGV{'--bgzip'}) {
				open $fout, "| body sort -T $tmpdir -V -k $f_chr,$f_chr -k $f_pos,${f_pos}n | uniq | bgzip -c > $outdir/$label.$suffix.gz" 
					or die "Cannot write open body sort pipe to bgzip compressed file $outdir/$label.$suffix.gz";
			}
			else {
				open $fout, "| body sort -T $tmpdir -V -k $f_chr,$f_chr -k $f_pos,${f_pos}n | uniq > $outdir/$label.$suffix" 
					or die "Cannot write open body sort pipe to write to $outdir/$label.$suffix";
			}
			print $fout join("\t", @fields), "\n";
		}
		else {
			unless(join("\t", @fields) eq join("\t", @$fnames)) {
				die "Input file $inprefix.$ii.$suffix does not have the same fields as $inprefix.$start.$suffix";
			}
		}
		while(my $dat = $it->()) {
			print $fout join("\t", @{$dat}{@fields}), "\n";
		}
	}
	close $fout;
	if ($ARGV{'--bgzip'}) {
		system("tabix -f -s $f_chr -b $f_pos -e $f_pos -c $fields[0] $outdir/$label.$suffix.gz");
	}
}


__END__

=head1 NAME 

merge_var_spits.pl -- Merge variants that split in multiple tables.

=head1 NOTES

Input file will be named as Prefix.N.suffix, where N=1,2....

They will be grouped together by a group file that has three columns: Group Start End

such that Prefix.Start.suffix ~ Prefix.End.suffix will be merged into Group.suffix

=head1 REQUIRED ARGUMENTS
 
=over
 
=item -[-]in[put] [=] <prefix>

Input variant table prefix.

=item -[-]out[dir] [=] <dir>

Output directory.
 
=item -[-]group [=] <file>

Group file or name for output file prefix.

=back
 
=head1 OPTIONS
 
=over

=item -[-]suffix [=] <str>

Split file name suffix, default txt.

=for Euclid:
	str.default: "txt"

=item -[-]nsplit [=] <number>

Number of split in the input, used when only one group in the output.

=item -[-]bgzip 

Sort and compress by bgzip, and index by tabix.

=item -[-]tmp[dir] [=] <dir>

Temp directory used for sorting.

=for Euclid:
	dir.default: "/tmp"

=back

