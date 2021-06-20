#!/usr/bin/env perl

use strict;
use warnings;
use Carp;
use IO::Dir;
use File::Basename;

# Collect BAM ID test results
my ($wrkdir, $outdir, $suff) = @ARGV;
unless(@ARGV == 2 || @ARGV == 3) {
	print STDERR "$0 WRKDIR OUTDIR\n";
	exit 1;
}
unless (-d $wrkdir && -d $outdir) {
	die "Cannot read working dir or output dir";
}

$suff = "selfSM" unless defined $suff;

my @resfiles = grep { /\.$suff$/ } IO::Dir->new($wrkdir)->read();

my %res;
foreach my $file (@resfiles) {
	(my $iid = $file) =~ s/\.$suff$//;
	open my $fin, "$wrkdir/$file" or die "Cannot open $wrkdir/$file";
	my $header = <$fin>; $header =~ s/^#//;
	my @fields = split(/\s+/, $header);
	my $line = <$fin>;
	my @values = split(/\s+/, $line);
	die "Number of fields does no match values" unless @fields == @values;
	my %dat = (SEQ_NAME => $iid);
	@dat{@fields} = @values;
	$res{$iid} = \%dat;
}

# Generating output
my @fields = qw|SEQ_NAME SEQ_ID CHIP_ID #SNPS #READS AVG_DP FREEMIX CHIPMIX|;
open my $fout, ">$outdir/summary.csv" or die "Cannot write to summary";
print $fout join(",", @fields), "\n";
while(my ($iid, $dat) = each %res) {
	print $fout join(",", @{$dat}{@fields}), "\n";
}
open $fout, ">$outdir/summary.README.txt" or die "Cannot write to README file";
print $fout <<EOF;
SEQ_NAME	Base file name of the sequence file after removing suffix
SEQ_ID		Sample ID obtaine from sequence file header (may be different from SEQ_NAME)
CHIP_ID		Sample ID in the genotype file to which the sequence data is compared, usually the same as SEQ_NAME but can be different if renamed 
#SNPS		Number of SNPs passing the criteria from VCF file
#READS		Total number of reads loaded from BAM file
AVG_DP		Average sequence depth at sites in VCF file
FREEMIX		Sequence-only estimate of contamination
CHIPMIX		Sequence+genotype estimate of contamination
EOF
