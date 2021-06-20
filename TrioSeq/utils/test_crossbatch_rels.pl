#!/usr/bin/env perl

use strict;
use warnings;
use File::Tempdir;
use List::MoreUtils qw|all|;
use File::Path qw|make_path|;
use Getopt::Euclid;
use Perl6::Slurp;
use File::Copy qw|copy|;
use Utils::File::Iter qw|iter_file|;

# check overlap samples between two batches
my $wrkdir = $ARGV{'--outdir'};
if (!-d $wrkdir) {
	make_path $wrkdir;
	die "Working directory $wrkdir cannot be created" unless -d $wrkdir;
}
else {
	if ($ARGV{'--force'}) {
		warn "Over-writing existing results in $wrkdir";
	}
	else {
		warn "Working directory $wrkdir already exist, exit";
		exit 1;
	}
}

if ($ARGV{'--b1'} eq $ARGV{'--b2'}) {
	die "Input files for two batches cannot be the same!";
}

# First note for the same FID, sequenced in different batches.
# Then for different FID but are related across batches.
if (all { -f "$ARGV{'--b1'}.$_" } qw|bim bed fam|) {
	system(qq|plink --bfile $ARGV{'--b1'} --not-chr x y --snps-only --make-bed --out $wrkdir/Batch_1|);
}
elsif (all { -f "$ARGV{'--b1'}.$_" } qw|vcf.gz ped|) {
	my $opt = "";
	if ($ARGV{'--pass-only'} && $ARGV{'--pass-only'} eq '1') {
		$opt = '--pass-only';
	}
	system(qq|vcf2bped.pl --vcf $ARGV{'--b1'}.vcf.gz --ped $ARGV{'--b1'}.ped --no-chrX --snp-only $opt --out $wrkdir/Batch_1 --clean|);
}
else {
	die "Cannot find genotype files for Batch 1";
}
# Remove duplicates
my @dups = find_dup_samps("$wrkdir/Batch_1.fam", $ARGV{'--ignore'});
if (@dups) {
	open my $frm, ">$wrkdir/Batch_1.removed" or die "Cannot write to removed list for Batch_1";
	print $frm join("\n", map { "$_->[0] $_->[1]" } @dups), "\n";
	close $frm;
	system(qq|plink --bfile $wrkdir/Batch_1 --remove $wrkdir/Batch_1.removed --make-bed --out $wrkdir/Batch_1|);
}


if (all { -f "$ARGV{'--b2'}.$_" } qw|bim bed fam|) {
	system(qq|plink --bfile $ARGV{'--b2'} --not-chr x y --snps-only --make-bed --out $wrkdir/Batch_2|);
}
elsif (all { -f "$ARGV{'--b2'}.$_" } qw|vcf.gz ped|) {
	my $opt = "";
	if ($ARGV{'--pass-only'} && $ARGV{'--pass-only'} eq '2') {
		$opt = '--pass-only';
	}
	system(qq|vcf2bped.pl --vcf $ARGV{'--b2'}.vcf.gz --ped $ARGV{'--b2'}.ped --no-chrX --snp-only $opt --out $wrkdir/Batch_2 --clean|);
}
else {
	die "Cannot find genotype files for Batch 2";
}
@dups = find_dup_samps("$wrkdir/Batch_2.fam", $ARGV{'--ignore'});
if (@dups) {
	open my $frm, ">$wrkdir/Batch_2.removed" or die "Cannot write to removed list for Batch_2";
	print $frm join("\n", map { "$_->[0] $_->[1]" } @dups), "\n";
	close $frm;
	system(qq|plink --bfile $wrkdir/Batch_2 --remove $wrkdir/Batch_2.removed --make-bed --out $wrkdir/Batch_2|);
}

if ($ARGV{'--align-fid'}) {
	my %fid1 = map { (split)[1,0] } slurp "$wrkdir/Batch_1.fam";
	my %fid2 = map { (split)[1,0] } slurp "$wrkdir/Batch_2.fam";
	my $flag;
	foreach my $iid (keys %fid1) {
		if (defined $fid2{$iid} && $fid1{$iid} ne $fid2{$iid}) {
			$flag = 1;
			print STDERR "Align FIDs for samples with identical IIDs\n";
			last;
		}
	}
	if ($flag) {
		my %updatefid2;
		foreach my $iid (keys %fid2) {
			if (defined $fid1{$iid} && $fid1{$iid} ne $fid2{$iid}) {
				$updatefid2{$fid2{$iid}} = $fid1{$iid}
			}
		}
		copy("$wrkdir/Batch_2.fam", "$wrkdir/Batch_2.fam.bak");
		open my $fin, "$wrkdir/Batch_2.fam.bak" or die "Cannot open Batch_2.fam.bak";
		open my $fout, ">$wrkdir/Batch_2.fam" or die "Cannot write to Batch_2.fam";
		while(<$fin>) {
			my @dat = split;
			if (defined $updatefid2{$dat[0]}) {
				$dat[0] = $updatefid2{$dat[0]};
			}
			print $fout join("\t", @dat), "\n";
		}
	}
}


unless($ARGV{'--t1'} ne $ARGV{'--t2'}) {
	die "Must assign different prefix labels to different batches.";
}
system(qq|merge_bped.pl --input $wrkdir/Batch_1 $wrkdir/Batch_2 --output $wrkdir/Merged --tags $ARGV{'--t1'},$ARGV{'--t2'} --tag-iid --tag-sep $ARGV{'--sep'}|);

#if ($ARGV{'--fast'}) {
#	system(qq|king --cpus $ARGV{'--nt'} -b $wrkdir/Merged.bed --duplicate --prefix $wrkdir/Merged|);
#}
#else {
	system(qq|king --cpus $ARGV{'--nt'} -b $wrkdir/Merged.bed --related --degree 2 --prefix $wrkdir/Merged|);
#}

# Parse the results files
# 1. between batch within family relatedness
# 2. between batch across family cryptic relatedness
my $sep = $ARGV{'--sep'};
my (@relpairs, @errpairs);
my $it = iter_file("$wrkdir/Merged.kin");
while(my $dat = $it->()) {
	my $ii = index($dat->{ID1}, $sep);
	my $jj = index($dat->{ID2}, $sep);
	$dat->{BATCH1} = substr($dat->{ID1}, 0, $ii);
	$dat->{BATCH2} = substr($dat->{ID2}, 0, $jj);
	$dat->{IID1} = substr($dat->{ID1}, $ii+1);
	$dat->{IID2} = substr($dat->{ID2}, $jj+1);
	if ($dat->{BATCH1}  ne $dat->{BATCH2}) {
		if ($dat->{IID1} eq $dat->{IID2}) {
			if ($dat->{Kinship} < 0.45) {
				push @errpairs, join("\t", map { $_ // "." } @{$dat}{qw|FID IID1 BATCH1 BATCH2 IBS0 Kinship IBD2Seg IBD1Seg InfType|});
			}
		}
		else {
			if ($dat->{Kinship} > 0.0884) {
				push @relpairs, join("\t", "Within-Family", map { $_ // "." }
					@{$dat}{qw|BATCH1 FID IID1 BATCH2 FID IID2 IBS0 Kinship IBD2Seg IBD1Seg InfType|});
			}
		}
	}
}
if (-f "$wrkdir/Merged.kin0") {
	my $it = iter_file("$wrkdir/Merged.kin0");
	while(my $dat = $it->()) {
		 my $ii = index($dat->{ID1}, $ARGV{'--sep'});
		 my $jj = index($dat->{ID2}, $ARGV{'--sep'});
		$dat->{BATCH1} = substr($dat->{ID1}, 0, $ii);
		$dat->{BATCH2} = substr($dat->{ID2}, 0, $jj);
		$dat->{IID1} = substr($dat->{ID1}, $ii+1);
		$dat->{IID2} = substr($dat->{ID2}, $jj+1);
		if ($dat->{BATCH1} ne $dat->{BATCH2}) {
			#if ($dat->{FID1}.":".$iid1 ne $dat->{FID2}.":".$iid2) {
				push @relpairs, join("\t", "Between-Family", map { $_ // "." }
					@{$dat}{qw|BATCH1 FID1 IID1 BATCH2 FID2 IID2 IBS0 Kinship IBD2Seg IBD1Seg InfType|});
			#	}
			#}
		}
		else {
			warn "Unexpected between-family relatedness within the same cohort\n";
			print STDERR join("\t", map { $_ // "." } @{$dat}{qw|FID1 ID1 FID1 ID2 IBS0 Kinship IBD2Seg IBD1Seg InfType|}), "\n";
		}
	}
}

if (@relpairs) {
	print STDERR "Writing related pairs to $wrkdir/relpairs.txt\n";
	open my $fout, ">$wrkdir/relpairs.txt" or die "Cannot write relpairs";
	print $fout join("\t", qw|TYPE BATCH1 FID1 IID1 BATCH2 FID2 IID2 IBS0 Kinship IBD2Seg IBD1Seg InfType|), "\n";
	print $fout join("\n", @relpairs), "\n";
}
else {
	print STDERR "No cross-batch relatedness was found.\n";
}

if (@errpairs) {
	print STDERR "Writing error pairs to $wrkdir/errpairs.txt\n";
	open my $ferr, ">$wrkdir/errpairs.txt" or die "Cannot write errpairs";
	print $ferr join("\t", qw|FID1 IID1 BATCH1 BATCH2 IBS0 Kinship IBD2Seg IBD1Seg InfType|), "\n";
	print $ferr join("\n", @errpairs), "\n";
}
else {
	print STDERR "No error pair was found.\n";
}

if ($ARGV{'--cleanup'}) {
	foreach my $suffix (qw|bim bed fam log|) {
		foreach my $prefix (qw|Batch_1 Batch_2 Merged|) {
			unlink("$wrkdir/$prefix.$suffix");
			if (-f "$wrkdir/$prefix.$suffix~") {
				unlink("$wrkdir/$prefix.$suffix~");
			}
		}
	}
}

sub find_dup_samps {
	my ($pedfile, $pattern) = @_;
	my $patregex = qr/$pattern/;
	my @dupids;
	my $fin = IO::File->new($pedfile) or die "Cannot reading $pedfile";
	while(<$fin>) {
		my ($fid, $iid) = (split)[0,1];
		if ($iid =~ /$patregex/) {
			#print STDERR $iid, "\n";
			push @dupids, [$fid, $iid];
		}
	}
	return @dupids;
}

__END__

=head1 NAME

crossbatch_rels -- Check sample relatedness between sequencing batches.

=head1 DESCRIPTION

This is a utlity script to check sample relatedness across batches.

For genotypes of each batch, we assume either BED/BIM/FAM files exist, or
VCF/PED files exist. For the latter, convertion with default parameters
will be done.

Samples with *identical FID and IID* across batches will be assumed to 
be identical, and if their genotypes is conflict of this identity,
error will be recorded. If in some data FID may be different for samples
with identical IIDs, there is an option to align FIDs so such sample
can also be checked.

Currently only support comparison between two batches, and for relatedness
closer than 2 degrees. Because distant relatedness can not be reliably detected
using common SNPs even within exome region.

=head1 REQUIRED ARGUMENTS

=over

=item -[-]b1 [=] <prefix>

Prefix for Batch 1.

=item -[-]t1 [=] <string>

Tag for Batch 1, will be used as prefix to IID.

=item -[-]b2 [=] <prefix>

Prefix for Batch 2.

=item -[-]t2 [=] <string>

Tag for Batch 2.

=item -[-]out[dir] [=] <dir>

Output and working directory.

=back

=head1 OPTIONS

=over

=item -[-]ignore [=] <string>

Ignore sample IDs matching the regex string.

=for Euclid:
	string.default: '_Re\d*$'

=item -[-]sep [=] <string>

ID prefix seperator string.

=for Euclid:
  string.default: ':'

=item -[-]force

Force writing to existing directory.

=item -[-]pass-only [=] <batch>

Use only variants tagged "PASS" in the specified batch.

=item -[-]align-fid

Assuming samples with identical IDs are the same sample, align FIDs
between two data set.

=item -[-]clean[up]

Clean up intermediate genotype files.

=item -[-]nt [=] <number>

Number of parallel threads for running king.

=for Euclid:
	number.default: 4

=back


