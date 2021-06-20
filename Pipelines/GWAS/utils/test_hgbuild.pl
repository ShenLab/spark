use strict;
use warnings;
use FindBin qw|$Bin|;
use Utils::File qw|open_file|;
use Utils::File::Iter qw|iter_file|;
use List::MoreUtils qw|all|;
use Getopt::Euclid;

use lib "$Bin/../../lib/";
use Shared qw|parse_fstr merge_tabfiles slurp_xref|;


my @fields;
if ($ARGV{'--fields'}) {
 	@fields = split(',', $ARGV{'--fields'});
}
else {
	if ($ARGV{'--input'} =~ /\.bim$/) {
		@fields = (1, 4, 2);
	}
	elsif ($ARGV{'--input'} =~ /\.pfb$/) {
		@fields = qw|Chr Position Name|;
	}
	else {
		die "Cannot recognize file format, must provide field names for Chrom,Position,MarkerID manually!";
	}
}
unless(@fields == 3) {
	die "Must provide field names for Chrom,Position,MarkerID";
}
print STDERR "The following fields will be parsed for Chrom, Position, MarkerID:\n",
	join("\t", @fields), "\n";


my %opt;
if (defined $ARGV{'--fsep'}) {
	$opt{fsep} = qr/$ARGV{'--fsep'}/;
}
if (all { /^\d+$/ } @fields) {
	$opt{header} = 0;
}
else {
	$opt{header} = 1;
}
#my $fsep = defined $ARGV{'--fsep'} ? qr/$ARGV{'--fsep'}/ : qr/\s+/;
my ($it, $fnames) = iter_file($ARGV{'--input'} eq '-' ? \*STDIN : $ARGV{'--input'},
	{ skip => $ARGV{'--skip'}, %opt });

foreach my $field (@$fnames) {
	unless(grep { $field eq $_ } @$fnames) {
		print "Field $field cannot be found in the input";
		exit 1;
	}
}

# Slurp chr/pos from input
my (%chrpos, $inchr);
while(my $dat = $it->()) {
	my ($chrom, $pos, $mkid) = @{$dat}{@fields};
	next unless $mkid =~ /^rs/;
	unless(defined $inchr) {
		$inchr = $chrom =~ /^chr/ ? 1 : 0;
	}
	$chrpos{$mkid} = "$chrom\t$pos";
}

unless(keys %chrpos) {
	print STDERR "No rs IDs can be found from input marker file!";
	exit 1;
}

my $flog;
if ($ARGV{'--output'}) {
	open $flog, ">$ARGV{'--output'}" or die "Cannot write to $ARGV{'--output'}";
	print $flog join("\t", qw|Marker  Chrom Position dbChr dbPos|), "\n";
}

# Now go-through dbSNP VCF files.
my (%match, $dbchr);
foreach my $dbsnp (@{$ARGV{'--dbsnp'}}) {
	print STDERR "Checking SNP positions from $dbsnp\n";
	my $fin = open_file($dbsnp);
	while(<$fin>) {
		next if /^#/;
		my ($chrom, $pos, $mkid) = (split)[0,1,2];
		unless(defined $dbchr) {
			$dbchr = $chrom =~ /^chr/ ? 1 : 0;
		}
		if (defined $chrpos{$mkid}) {
			if ($inchr == 1 && $dbchr == 0) {
				$chrom = hg_chrom($chrom);
			}
			elsif ($inchr == 0 && $dbchr == 1) {
				$chrom =~ s/^chr//; $chrom = 'MT' if $chrom eq 'M';
			}
			$match{$dbsnp}{Tot} ++;
			if ($chrpos{$mkid} eq "$chrom\t$pos") {
				$match{$dbsnp}{Match} ++;
			}
			else {
				#if ($ARGV{'--verbose'}) {
				#	print STDERR "Inconsistent position for $mkid: $chrpos{$mkid} ne $chrom\t$pos\n";
				#}
				if (defined $flog) {
					print $flog join("\t", $mkid, $chrom, $pos, $chrpos{$mkid}), "\n";
				}
			}
			last if $match{$dbsnp}{Tot} >= $ARGV{'--max'};
		}
	}
}

print "Total number of SNPs read from input file: ", scalar(keys %chrpos), "\n";
foreach my $dbsnp (sort keys %match) {
	printf "%d can be found in %s, %d have matched position\n", 
			$match{$dbsnp}{Tot} // 0, $dbsnp, $match{$dbsnp}{Match} // 0; 
}


__END__

=head1 NAME

test_hgbuild -- Quick check the genome reference build for SNP markers .bim file.

=head1 NOTE

It requires some markers in the bim file have RS IDs and checks against dbSNP database in VCF format.

=head1 REQUIRED ARGUMENTS

=over

=item -[-]in[put] [=] <file>

Input file to be checked.

=for Euclid:
	file.type: readable

=item -[-]dbsnp[s] [=] <files>...

One or more VCF files of dbSNP database.

=for Euclid:
	files.type: readable

=back

=head1 OPTIONS

=over

=item -[-]fields [=] <fstr>

Comma separated list of field names for Chrom,Position and Marker ID. Must be in this order.

=item -[-]skip [=] <nrow>

Number of rows to skip at the begniing.

=item -[-]fsep [=] <regex>

Field separator of the input file.

=item -[-]output [=] <file>

Output log file recording markers with inconsistent positions.

=item -[-]max [=] <number>

Max. number of SNPs to check.

=for Euclid:
	number.default: 1000

=item -[-]verbose

Print out err positions.

=back

=cut
