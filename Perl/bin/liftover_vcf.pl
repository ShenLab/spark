#!/usr/bin/env perl

use strict;
use warnings;
use Carp;
use Data::Dumper;
use List::MoreUtils qw|any all|;
use File::Tempdir;
use IO::Uncompress::Gunzip;
use Getopt::Euclid;
use FaSlice;
use Genet::Var qw|normalize|;
use Utils::Seq qw|rev_comp|;
use Utils::File::Iter qw|iter_file|;
use Genome::UCSC qw|is_hgchr hg_chrom|;
use Genome::UCSC::Liftover qw|check_chain_chrpref|;
use Genet::File::VCF qw|check_vcf_chrpref|;

# Check pre-requisite
my $helpline = `CrossMap.py 2>&1 | head -n1`;
if ($helpline =~ /command not found/i) {
	die "Cannot find CrossMap command line tool";
}
else {
	my ($ver) = ($helpline =~ /(v[\.0-9]+)/);
	unless($ver eq 'v0.2.7') {
		print $helpline, "\n";
		die "Only compatible with CrossMap v0.2.7";
	}
}

my ($wrkdir, $tmpdir);
if ($ARGV{'--wrkdir'}) {
	$wrkdir = $ARGV{'--wrkdir'};
	mkdir $wrkdir;
	croak "Cannot create working directory" unless -d $wrkdir;
}
else {
	$tmpdir = File::Tempdir->new();
    $wrkdir = $tmpdir->name;
}

my ($inchr, $outchr) = check_chain_chrpref($ARGV{'--chain'});
my $vcfchr = check_vcf_chrpref($ARGV{'--invcf'});
# When chromsome names are different, we will issue a warning
# the results on canonicla chromosomes should be almost the same, but non-chr contigs may not be
# processed correctly. If those non-chr contigs are not of primary interest, the warning can be
# safely ignored.
unless($vcfchr == $inchr) {
	warn "Chromosome nomenclatures in VCF and chain files are different";
}

# First extract all variants and create interval BED6 file for each of them
# for testing strandness after liftover
{
	my $it = iter_file($ARGV{'--invcf'}, { ignore => qr/^##/, 
		subset => [qw|#CHROM POS ID REF ALT QUAL FILTER|] });
	open my $fout, ">$wrkdir/input.bed" or die "Cannot write to input.bed";
	open my $fvcf, ">$wrkdir/input.vcf" or die "Cannot write to input.vcf";
	my %known;
	if ($vcfchr == $inchr) {
		while(my $dat = $it->()) {
			my $locid = join(':', @{$dat}{qw|#CHROM POS REF|});
			next if $known{$locid};
			print $fout join("\t", $dat->{'#CHROM'}, $dat->{POS}-1, 
				$dat->{POS}+length($dat->{REF}), $locid, 100, '+'), "\n";
			print $fvcf join("\t", @{$dat}{qw|#CHROM POS ID REF ALT QUAL FILTER|}, $locid), "\n";
			$known{$locid} = 1;
		}
	}
	elsif ($vcfchr == 0 && $inchr == 1) {
		while(my $dat = $it->()) {
			my $locid = join(':', @{$dat}{qw|#CHROM POS REF|});
			next if $known{$locid};
			$dat->{'#CHROM'} = hg_chrom($dat->{'#CHROM'});
			print $fout join("\t", $dat->{'#CHROM'}, $dat->{POS}-1, 
				$dat->{POS}+length($dat->{REF}), $locid, 100, '+'), "\n";
			print $fvcf join("\t", @{$dat}{qw|#CHROM POS ID REF ALT QUAL FILTER|}, $locid), "\n";
			$known{$locid} = 1;
		}
	}
	elsif ($vcfchr == 1 && $inchr == 0) {
		while(my $dat = $it->()) {
			my $locid = join(':', @{$dat}{qw|#CHROM POS REF|});
			next if $known{$locid};
			$dat->{'#CHROM'} =~ s/^chr//; $dat->{'#CHROM'} = 'MT' if $dat->{'#CHROM'} eq 'M';
			print $fout join("\t", $dat->{'#CHROM'}, $dat->{POS}-1, 
				$dat->{POS}+length($dat->{REF}), $locid, 100, '+'), "\n";
			print $fvcf join("\t", @{$dat}{qw|#CHROM POS ID REF ALT QUAL FILTER|}, $locid), "\n";
			$known{$locid} = 1;
		}
	}
	else {
		die "Should not reach here!";
	}
}


# Perform liftover on input bed and VCF
system(qq|CrossMap.py bed $ARGV{'--chain'} $wrkdir/input.bed $wrkdir/output.bed|);
system(qq|CrossMap.py vcf $ARGV{'--chain'} $wrkdir/input.vcf $ARGV{'--seq'} $wrkdir/output.vcf|);

# Collecting CrossMap Results:
# 1. Mark variant positions mapped to negative strand
# 2. Mark variant positions (including 1bp to the right) mapped to >1 new positions 
#    also mark with inconsistent new strands
#    Note: CrossMap split the intervals, if there are gaps. We will make use of this.
# 3. Collect new mapping position and new reference alleles

# Variant positions with different reference alleles will be skipped.
# So variants with Ref/Alt flipped in the new genome assembly will not be in the output.
# Most of such cases will be common variants.

print STDERR "Collecting liftover results.\n";
my (%neg, %nremap, %newpos);
{
	open my $fbed, "$wrkdir/output.bed" or die "Cannot open $wrkdir/output.bed";
	while(<$fbed>) {
		my ($chrom, $start, $end, $name, $score, $strand) = split;
		$nremap{$name} ++;
		if ($strand eq '-') {
			$neg{$name} = 1;
		}
	}
	open my $fvcf, "$wrkdir/output.vcf" or die "Cannot open $wrkdir/output.vcf";
	open my $ferr, ">$wrkdir/output.err.txt" or die "Cannot write to output.err.txt";
	while(<$fvcf>) {
		next if /^#/;
		my @dat = split;
		unless(defined $nremap{$dat[7]}) {
			warn "The locus $dat[7] was remapped in BED but not VCF" if $ARGV{'--verbose'};
			print $ferr join("\t", @dat[0,1,3], "NRemap=0"), "\n";
			next;
		}
		if ($nremap{$dat[7]} != 1) {
			print $ferr join("\t", @dat[0,1,3], "NRemap=".$nremap{$dat[7]}), "\n";
			next;
		}

		# Check refseq for variants mapped to neg strand
		my $refseq_old = (split(':', $dat[7]))[2];
		my $refseq_new = $dat[3];
		if (defined $neg{$dat[7]}) {
			if ($refseq_new eq rev_comp($refseq_old)) {
				$newpos{$dat[7]} = join(":", @dat[0,1,3]);
			}
			else {
				# skip variants with unmatched refseq
				warn "The locus $dat[7] is mapped to negative strand, but but the new refseq($refseq_new) is different from rev comp of the old($refseq_old)"
					if $ARGV{'--verbose'};
				print $ferr join("\t", @dat[0,1,3], "OirgRef=".rev_comp($refseq_old)), "\n";
			}
		}
		else {
			if ($refseq_new eq $refseq_old) {
				$newpos{$dat[7]} = join(":", @dat[0,1,3]);
			}
			else {
				# skip variants with unmatched refseq
				warn "The locus $dat[7] is mapped to plus strand, but the new refseq($refseq_new) is different from old($refseq_old)"
					if $ARGV{'--verbose'};
				print $ferr join("\t", @dat[0,1,3], "OirgRef=".$refseq_old), "\n";
			}
		}
	}
}


my $sq = FaSlice->new(file => $ARGV{'--seq'});

print STDERR "Remapping the variants in input file.\n";

my $outvcf = $ARGV{'--output'};
unless ($outvcf =~ /\.vcf\.gz$/) {
	$outvcf .= ".vcf.gz";
}
open my $fgz, "| vcf-sort -c | bgzip -c > $outvcf" or die "Cannot open pipe";

my $fin;
if ($ARGV{'--invcf'} =~ /\.vcf\.gz$/) {
	$fin = IO::Uncompress::Gunzip->new($ARGV{'--invcf'}, MultiStream => 1);
}
else {
	open $fin, $ARGV{'--invcf'} or die "Cannot open $ARGV{'--invcf'}";
}

while(<$fin>) {
	if (/^#/) {
		print $fgz $_;
	}
	else {
		my @dat = split;
		# Collect new mapping position and new reference sequence
		my $name = join(":", @dat[0,1,3]);
		next unless defined $newpos{$name};

		# Unpack new positions and new reference allele
		($dat[0], $dat[1], $dat[3]) = split(':', $newpos{$name});
		if ($ARGV{'--hgchr'}) {
			next unless is_hgchr($dat[0]);
		}

		# Processing alt alleles
		my @alts;	
		if (defined $neg{$name}) {
			@alts = map { rev_comp($_) } map { uc($_) } split(',', $dat[4]);	
		}
		else {
			@alts = map { uc($_) } split(',', $dat[4]);
		}
		#next unless all { /^[ACGTN\*]+$/ } @alts;

		if (length($dat[3]) == 1 && (all { length($_) == 1 } @alts)) {
			$dat[4] = join(',', @alts);
		}
		else {
			my $var = { CHROM => $dat[0], POS => $dat[1], REF => $dat[3], ALT => [ @alts ] };
			normalize($var, $sq);
			$dat[1] = $var->{POS};
			$dat[3] = $var->{REF};
			$dat[4] = join(',', @{$var->{ALT}});
		}

		if ($ARGV{'--nochr'} && $outchr) {
			$dat[0] =~ s/^chr//;
		}
		print $fgz join("\t", @dat), "\n";
	}
}
close $fgz;

# Then index the resulting VCF file
if (-f $outvcf) {
	system("tabix -p vcf $outvcf");
}


__END__

=head1 NAME

liftover_vcf.pl -- Lift VCF file with genotypes.

=head1 USAGE

liftover_vcf.pl -in INPUT -chain CHAIN -seq HGSEQ -out OUTPUT

=head1 DESCRIPTION

To lift a VCF file, there are several new complications as compared with simply liftover
over bialleli variants. The existance of multi-allelic variants, changes in ref allele
also causes collateral changes in genotype fields, etc. 

The following approach was used to liftover a VCF file. We will store the genomic position of
each reference allele to a BED6 file and used CrossMap to find the new genomic position and
strand. We trim the original VCF to site only, and used CrossMap to find the new position and
new reference allele. We will then compare with new reference allele with the original one.
When a variant was mapped to a different strand, we will rev comp all original alleles.
If a variant's new reference allele remains the same after accounting for strand flip, we will only
update the position and keep other VCF entries as they are.

CrossMap after v0.3.4 also account for the reference allele changes for strand flip, but it does
not post-process the indels, nor does it check if reference allele remains the same.

=head1 NOTE

The script is only compatible with CrossMap v0.2.7. CrossMap v.2.8 failed to map variants using
UCSC's chain files for unclear reason. 

The script only support liftover sequence variants in VCF. Structual variants will be ignored.

We do not update sequence contig information in the new VCF header.

=head1 REQUIRED ARGUMENTS

=over

=item -[-]in[vcf] [=] <file>

Input data file for liftover.

=for Euclid:
  file.type: readable

=item -[-]chain [=] <file>

Liftover chain file. When contig/chromsome name nomenclature does no match the one used in VCF, 
it still works for canonical chromosomes but results on non-chr contigs may not be accurate.

=for Euclid:
  file.type: readable

=item -[-]seq [=] <file>

Reference sequence for the target genome, in FASTA format; should be indexed.
Must match the target assembly in chain file. 

This file is used to check refseq in the new assembly. We assume the original
input has been checked for refseq in the orignal assembly.

=for Euclid:
  file.type: readable

=item -[-]out[put] [=] <file>

Output file name or prefix. The output format will be big-zipped VCF.

=back

=head1 OPTIONS

=over

=item -[-]hgchr

Keep only canonical human chromosome.

=item -[-]no[-]chr

Remove 'chr' prefix in chromosome names in the output VCF file.

=item -[-]wrkdir [=] <dir>

Working directory, if not provided, will use tempdir.

=item -[-]verbose

Output all warning messages.

=back

=cut

