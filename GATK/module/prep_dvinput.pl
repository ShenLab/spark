use strict;
use warnings;
use Cwd qw|abs_path|;
use File::Path qw|make_path|;
use List::MoreUtils qw|all none|;
use Perl6::Slurp;
use Utils::File::Iter qw|iter_file|;

my ($intab, $keycols, $outdir) = @ARGV;

# List of samples to keep
# Including path to their original VCF.
my %keep = map { (split)[0,1] } slurp "$outdir/keep.txt";

# IID column should appear first in the key columns
my $idcol = (split(',', $keycols))[0];

my ($it, $fnames) = iter_file($intab, { fsep => qr/\t/ });
unless(grep { $idcol eq $_ } @$fnames) {
	die "Cannot find IID column: $idcol";
}

# Prepare links to VCFs
make_path "$outdir/vcf";
while(my ($iid, $vcf) = each %keep) {
	next if $vcf eq '.';
	unless(-f $vcf) {
		die "Cannot find VCF for $iid: $vcf";
	}
	symlink abs_path($vcf), "$outdir/vcf/$iid.vcf.gz";
}

# First extract variants
open my $fout, ">$outdir/input.txt" or die "Cannot write to input.txt";
print $fout join("\t", @$fnames), "\n";
while(my $dat = $it->()) {
	if (defined $keep{$dat->{$idcol}}) {
		print $fout join("\t", @{$dat}{@$fnames}), "\n";
	}
}


