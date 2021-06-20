#!/usr/bin/env perl

use strict;
use warnings;
use Carp;
use IO::File;
use File::Copy;
use File::Which;
use Perl6::Slurp;
use Getopt::Euclid;
use List::MoreUtils qw|all|;
use Genome::UCSC qw|hg_par|;
use Genet::File::VCF;

my $vcf = Genet::File::VCF->new($ARGV{'--vcf'});

# Read pedfile lines
my %pedline;
if ($ARGV{'--ped'}) {
	%pedline = map { my @row = split; $row[1] => [@row] } slurp $ARGV{'--ped'};
}
else {
	%pedline = map { $_ =>  [$_, $_, 0, 0, 0, 0] } $vcf->get_sampids;
}

# Read header information to determine filter flags.
my @incflags = $vcf->list_vqsr_filter({ SNP => $ARGV{'--vqsr-snp'}, INDEL => $ARGV{'--vqsr-indel'} });
print STDERR "The following filters will used to include variants\n";
print STDERR join("\t", @incflags), "\n";

# Create intermediate file
my $prefix = $ARGV{'--out'};

my $plink = $ARGV{'--plink'} // which("plink");
croak "Cannot find plink" unless defined $plink;


foreach my $set (qw|auto chrX|) {
	# Construct PLINK2 command line
	if ($set eq 'chrX') {
		next if $ARGV{'--no-chrX'};
	}
	# We will process autosome and chrX separately
	my $plinkcmd = "$plink --vcf $ARGV{'--vcf'} --const-fid --allow-extra-chr 0 ";
	not defined $ARGV{'--multi'} and $plinkcmd .= "--biallelic-only ";
	$plinkcmd .= "--vcf-filter ".join(" ", grep { $_ ne 'PASS' } @incflags);
	if ($set eq 'chrX') {
		$plinkcmd .= " --chr X --vcf-min-qual $ARGV{'--vq-chrX'} --vcf-min-gq $ARGV{'--gq-chrX'} ";
	}
	else {
		$plinkcmd .= " --chr 1-22 --vcf-min-qual $ARGV{'--vq'} --vcf-min-gq $ARGV{'--gq'} ";
	}
	$plinkcmd .= "--make-bed --out $ARGV{'--out'}.$set.pass1";
	
	print ucfirst($set), " Pass1: convert VCF to binary PED format\n";
	print "-------------\n";
	print $plinkcmd,  "\n";
	print "-------------\n";
	system($plinkcmd);

	# Patch the resulting output files
	copy("$ARGV{'--out'}.$set.pass1.fam", "$ARGV{'--out'}.$set.pass1.fam.bak");
	my $fin = IO::File->new("$ARGV{'--out'}.$set.pass1.fam.bak");
	my $fout = IO::File->new("$ARGV{'--out'}.$set.pass1.fam", "w");
	while(<$fin>) {
		my $iid = (split)[1];
		croak "pedline for $iid cannot be found" unless defined $pedline{$iid};
		print $fout join("\t", @{$pedline{$iid}}), "\n";
	}
	$fout->close;

	# Patch variant IDs in bim file, and variants with '*' alleles.
	# Also update chrX name for markers in in PAR regions
	my %known;
	my @badvars;
	copy("$ARGV{'--out'}.$set.pass1.bim", "$ARGV{'--out'}.$set.pass1.bim.bak");
	$fin = IO::File->new("$ARGV{'--out'}.$set.pass1.bim.bak");
	$fout = IO::File->new( "$ARGV{'--out'}.$set.pass1.bim", "w");
	while(<$fin>) {
		my ($chr, $mkid, $gpos, $pos, $a1, $a2) = split;
		# Patch chrX
		my $exclude;
		if ($chr eq '23') {
			if ( hg_par('X', $pos, $ARGV{'--hg'}) > 0 ) {
				$chr = '25';
			}
		}
		if ($mkid eq '.') {
			$mkid = varid($chr, $pos, $a1, $a2);	
		}
		# Patch duplicated marker ID
		# A few duplication should be ok, when large number of duplicated markers
		# are found, it is possible due to the error in VCF file processing
		if(defined $known{$mkid}) {
			$mkid = $mkid.".$known{$mkid}";
		}
		# Identify other bad variants for exclusion
		if ($chr eq '0' || $a1 !~ /^[ACGT]+$/ || $a2 !~ /^[ACGT]+$/) {
			$exclude = 1;
		}
		if ($ARGV{'--snp-only'}) {
			unless ($a1 =~ /^[ACGT]$/ && $a2 =~ /^[ACGT]$/) {
				$exclude = 1;
			}
		}
		print $fout join("\t", $chr, $mkid, $gpos, $pos, $a1, $a2), "\n";
		push @badvars, $mkid if $exclude;
		$known{$mkid} ++;
	}
	$fout->close;

	print ucfirst($set), " Pass2: Fix split chromosome, and remove bad markers\n";
	$plinkcmd = "$plink --bfile $ARGV{'--out'}.$set.pass1 --make-bed --out $ARGV{'--out'}.$set.pass2 ";
	if (@badvars) {
		my $fexl = IO::File->new("$ARGV{'--out'}.$set.pass2.excl", "w");
		print $fexl join("\n", @badvars), "\n";
		$plinkcmd .= "--exclude $ARGV{'--out'}.$set.pass2.excl ";
	}
	print "-------------\n";
	print $plinkcmd, "\n";
	print "-------------\n";
	system($plinkcmd);

	# Second round plink QC
	print ucfirst($set), " Genotype QC\n";
	$plinkcmd = "$plink --bfile $ARGV{'--out'}.$set.pass2 --make-bed --out $ARGV{'--out'}.$set " ;
	if ($set eq 'chrX') {
		$plinkcmd .= "--maf $ARGV{'--maf'} --geno $ARGV{'--geno-chrX'} --hwe $ARGV{'--hwe-chrX'}";
	}
	else {
		$plinkcmd .= "--maf $ARGV{'--maf'} --geno $ARGV{'--geno'} --hwe $ARGV{'--hwe'}";	
	}
	print "-------------\n";
	print $plinkcmd, "\n";
	print "-------------\n";
	system($plinkcmd);
}

# NOTE: A small number of duplicated variants should not influence downstream analysis
# to fully patch the duplication, use merge_bped.pl 

# Finally merge auto and chrX to create the final genotype file
unless($ARGV{'--no-chrX'}) {
	if (all { -f "$ARGV{'--out'}.chrX.$_" } qw|bim bed fam|) {
		my $plinkcmd = "$plink --bfile $ARGV{'--out'}.auto --bmerge $ARGV{'--out'}.chrX --allow-no-sex --make-bed --out $ARGV{'--out'}";
		print "-------------\n";
		print $plinkcmd, "\n";
		print "-------------\n";
		system($plinkcmd);
	}
}
else {
	foreach my $suf (qw|bim bed fam|) {
		copy("$ARGV{'--out'}.auto.$suf", "$ARGV{'--out'}.$suf");
	}
}


if ($ARGV{'--clean'}) {
	foreach my $set (qw|auto chrX|) {
		foreach my $suf (qw|bim bed fam hh nosex bim.bak fam.bak excl log nosex|) {
			unlink $ARGV{'--out'}.".$set.$suf";
			foreach my $ii (1,2) {
				unlink $ARGV{'--out'}.".$set.pass$ii.$suf";
			}
		}
	}
}

sub varid {
	my ($chr, $pos, $a1, $a2) = @_;
	my $mkid;
	if ($a1 =~ /^[ACGT]+$/i && $a2 =~ /^[ACGT]+$/i && length($a1) < 5 && length($a2) < 5 ) {
		$mkid = "$chr-$pos-".join("-", sort($a1, $a2));
	}
	else {
		$mkid = "$chr-$pos";
	}
	return $mkid;

}

__END__

=head1 NAME

vcf2bped -- Convert VCF to PLINK bped format.

=head1 USAGE

vcf2bped [options] -vcf GENOFILE -ped FAMFILE -out PREFIX

=head1 DESCRIPTION

The utility wrap up PLINK2's function and properly deal with PED file, variant
IDs, and allele issues.

It accept parameters for variant and genotype level QC, but do not remove
individuals. The default parameters are generally good for genotypes
called by GATK.

It runs PLINK2 twice, first convert VCF to bped, then perform initial QC
on bped.

Only autosome and chrX will be processed!

=head1 REQUIRED ARGUMENTS

=over

=item -[-]vcf [=] <file>

Input VCF file.

=for Euclid
  file.type: readable

=item -[-]out[put] [=] <prefix>

Output file prefix.

=back

=head1 OPTIONS

=over

=item -[-]ped [=] <file>

Input pedigree file. Create a ped file with minimal information when not provided.

=for Euclid:
  file.type: readable

=item -[-]snp-only

Only keep SNVs.

=item -[-]no[-]chrX

Do not process chrX variants.

=item -[-]multi

Default is to keep only biallelic sites.

The number of biallelic sites will reduce with increasing number of samples.
Use this to switch to read multi-allelic sites, but keep two most common alleles.
The genotypes involving remaining alleles will be set to missing.

=item -[-]maf [=] <thres>

PLINK's MAF filter threshold (default: 0.01).

=for Euclid:
  thres.default: 0.01

=item -[-]vqsr-snp [=] <thres>

GATK's VQSR filter threshold for SNP. Will convert to FILTER tags (default: 99.8).

=for Euclid:
  thres.default: 99.8

=item -[-]vqsr-indel [=] <thres>

GATK's VQSR filter threshold for INDEL. Will convert to FILTER tags (default: 99.0).

=for Euclid:
  thres.default: 99.0

=item -[-]gq [=] <thres>

Minimal GQ threshold (default: 40).

=for Euclid:
  thres.default: 40

=item -[-]gq-chrX [=] <thres>

Minimal GQ threshold for chrX (default: 30).

=for Euclid:
  thres.default: 30

=item -[-]vq [=] <thres>

Minmial variant level QUAL threshold (default: 100).

=for Euclid:
  thres.default: 100

=item -[-]vq-chrX [=] <thres>

Minmial variant level QUAL threshold for chrX (default: 100).

=for Euclid:
  thres.default: 100

=item -[-]geno [=] <thres>

PLINK's genotype missingness filter threshold (default: 0.02).

=for Euclid:
  thres.default: 0.02

=item -[-]geno-chrX [=] <thres>

PLINK's genotype missingness filter threshold for chrX (default: 0.03).

=for Euclid:
  thres.default: 0.03

=item -[-]hwe [=] <thres>

PLINK's HWE filter threshold (default: 0).

=for Euclid:
  thres.default: 0

=item -[-]hwe-chrX [=] <thres>

PLINK's HWE filter threshold for chrX (default: 0).

=for Euclid:
  thres.default: 0

=item -[-]hg [=] <build>

Human genome assembly version (default: b37).

=for Euclid:
  build.default: "b37"

=item -[-]plink [=] <path>

Provide an alternative path to plink2 executable.

=item -[-]clean

Clean up the intermdiate file.

=cut



