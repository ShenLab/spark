#!/usr/bin/env perl

use strict;
use warnings;
use Carp;
use IO::File;
use Cwd qw|abs_path|;
use List::MoreUtils qw|any|;
use File::Copy;
use File::Which;
use File::Tempdir;
use File::Copy qw|copy|;
use File::Basename qw|basename dirname|;
use Perl6::Slurp;
use Utils::List qw|which_min|;
use Utils::Seq qw|is_snv is_asym rev_comp|;
use Utils::File::Iter qw|slurp_file|;
use Genome::UCSC qw|hg_chrom hg_chr|;
use Genome::UCSC::TwoBit;
use Genome::UCSC::Liftover;
use Getopt::Euclid;

my $plink;
if ($ARGV{'--plink'}) {
	$plink = $ARGV{'--plink'};
}
else {
	$plink = which("plink") or croak "Cannot find plink";
}
$plink .= " --allow-no-sex";

my @infiles = @{$ARGV{'--input'}};
my @labels;
if ($ARGV{'--tags'}) {
	@labels = split(/[,;]/, $ARGV{'--tags'});
	croak "Different number of labels and infiles" unless @infiles == @labels;
}

my ($dir, $tmpdir);
if ($ARGV{'--tmpdir'}) {
 	$tmpdir = $ARGV{'--tmpdir'};
 	mkdir $tmpdir unless -d $tmpdir;
}
else {
	$dir = File::Tempdir->new();
	$tmpdir = $dir->name;
	if (!-d $tmpdir) {
		croak "$tmpdir does not exist!";
	}
}

my @procfiles;
for (my $ii = 0; $ii < @infiles; $ii ++) {
	my $ind = $ii + 1;
	push @procfiles, proc_one_bped($infiles[$ii], "$tmpdir/in-$ind", $labels[$ii], $ARGV{'--processed'});
}


my $merged;
if (@infiles > 1) {
	$merged = merge_bpeds(\@procfiles, $ARGV{'--no-flip'});
}
else {
	$merged = $procfiles[0];
}

my $lifted;
if ($ARGV{'--chain'}) {
	$lifted = liftover_bped($merged, $ARGV{'--chain'});
}
else {
	$lifted = $merged;
}

my $aligned;
if ($ARGV{'--seq'}) {
	$aligned = refalign_bped($lifted, $ARGV{'--seq'});
}
else {
	$aligned = $lifted;
}

my $plinkcmd = "$plink --bfile $aligned --make-bed --out $ARGV{'--output'}";
system($plinkcmd); # or croak "Plink runtime error when recoding final output";


# Standard procedure for one bped:
# 1. Keep only biallelic asymmetric SNPs (unless specified otherwise);
# 2. Exclude multi-allelic SNPs, i.e those mapped to the same location but different alleles;
# 3. Exclude one of the duplicated SNPs that mapped to the same location and having same alleles;
# 4. Standardize SNP names to CHR-POS-sort(A1-A2).
# 5. [Optionally], add prefix to sample FID/IID
# All the operations will be performed in a tempdir
# Return output file prefix, which should be "outfile-num"
sub proc_one_bped {
	my ($infile, $outfile, $tag, $processed) = @_;

	symlink(abs_path($infile.".bed"), $outfile.".bed");
	
	if (defined $tag) {
		open my $fin, $infile.".fam" or croak "Cannot open $infile.fam";
		open my $fout, ">", $outfile.".fam" or croak "Cannot open $outfile.fam";
		while(<$fin>) {
			my @row = split;
			if ($ARGV{'--tag-iid'}) {
				croak "Can only specify one of --tag-iid/--tag-fid" if $ARGV{'--tag-fid'};
				$row[1] = "$tag".$ARGV{'--tag-sep'}.$row[1];
			}
			elsif ($ARGV{'--tag-fid'}) {
				croak "Can only specify one of --tag-iid/--tag-fid" if $ARGV{'--tag-iid'};
				$row[0] = "$tag".$ARGV{'--tag-sep'}.$row[0];
			}
			else {
				$row[0] = "$tag".$ARGV{'--tag-sep'}.$row[0];
				$row[1] = "$tag".$ARGV{'--tag-sep'}.$row[1];
			}
			print $fout join("\t", @row), "\n";
		}
	}
	else {
		symlink(abs_path($infile.".fam"), $outfile.".fam");
	}

	if ($processed) {
		symlink(abs_path($infile.".bim"), $outfile.".bim");
		return $outfile;
	}
	else {
		open my $fmap, $infile.".bim" or croak "Cannot open $infile.bim";
		open my $fbim, ">", $outfile.".bim" or croak "Cannot open $outfile.bim";
		while(<$fmap>) {
			my ($chr, $rsid, $gpos, $pos, $a1, $a2) = split;
			# ($a1, $a2) = sort($a1, $a2);  # The sort here causes error!
			my $varid = "$chr-$pos-".join('-', sort($a1, $a2));
			print $fbim join("\t", $chr, $varid, $gpos, $pos, $a1, $a2), "\n";
		}
		$fbim->close;
		return dedup_bped($outfile);
	}	
}

# Remove duplicated markers
# Assume infile has marker ID standardized.
# For markers mapped at the same location but different alleles, all of them will be removed
# For markers mapped at the same location with the same alleles, keep on with higher genotyping rate.
# Also remove markers not on standard chromosomes and non-SNP variants
sub dedup_bped {
	my ($infile) = @_;

	my (%known, %extract, %dups, %snppos);

	# Search for duplicated markers and patch marker IDs
	copy($infile.".bim", $infile.".bim.bak");
	open my $fin, $infile.".bim.bak" or croak "Cannot open $infile.bim.bak";
	open my $fout, ">", $infile.".bim" or croak "Cannot open $infile.bim";
	while(<$fin>) {
		my ($chr, $mkid, $gpos, $pos, $a1, $a2) = split;
		my $varid;
		if ($known{$mkid}) {
			$varid = $mkid.".$known{$mkid}";
		}
		else {
			$varid = $mkid;
		}
		print $fout join("\t", $chr, $varid, $gpos, $pos, $a1, $a2), "\n";

		if (is_snv($a1, $a2) && $chr > 0 && $chr < 26) {
			$snppos{"$chr-$pos"} ++ unless $known{$mkid};
			# Also require a valid chromosome code
			if ($ARGV{'--keep-sym'} || is_asym($a1, $a2)) {
				$extract{$varid} = 1;
				if ($varid =~ /\.\d+$/) {
					push @{$dups{$mkid}}, $varid;
				}
			}
		}
		$known{$mkid} ++;
	}
	$$fout->close;

	# After patching, extract all bi-allelic asym SNPs that are uniquely mapped
	open my $fext, ">", "$infile-snps.extract" or croak "Cannot open $infile-snps.extract";
	foreach my $varid (keys %extract) {
		my ($chr, $pos) = (split('-', $varid))[0, 1];
		if ($snppos{"$chr-$pos"} == 1) {
			print $fext $varid, "\n";
		}
	}
	$fext->close;
	system("$plink --bfile $infile --extract $infile-snps.extract --make-bed --out $infile-snps");

	# Check for SNPs in the same position with same alleles
	if (keys %dups) {
		foreach my $mkid (keys %dups) {
			push @{$dups{$mkid}}, $mkid;
		}
		my $fdup = IO::File->new("$infile-snps.dups", "w");
		print $fdup join("\n", map { @$_ } values %dups), "\n";
		$fdup->close;
		system("$plink --bfile $infile-snps --extract $infile-snps.dups --missing --out $infile-snps");
		# pick the low missing rate for duplicated markers
		my %lmiss = map { $_->{SNP} => $_->{F_MISS} } slurp_file("$infile-snps.lmiss");
		# Then for each duplicated markers, select one with highest genotyping rate
		my @excl;
		foreach my $dupmk (values %dups) {
			my @dupmks = grep { defined $lmiss{$_} } @$dupmk;
			next unless @dupmks > 1;
			my @gmiss = @lmiss{@dupmks};
			my $ii_min = which_min(@gmiss);
			for(my $ii = 0; $ii < @dupmks; $ii ++) {
				next if $ii == $ii_min;
				push @excl, $dupmks[$ii];
			}
		}
		my $fexl = IO::File->new("$infile-dedup.exclude", "w");
		print $fexl join("\n", @excl), "\n";
		$fexl->close;
		my $cmd = "$plink --bfile $infile-snps --exclude $infile-dedup.exclude --make-bed --out $infile-dedup";
		system($cmd); # or croak "dedup_bped: plink runtime error when processing $infile";

		# Patch the marker name
		copy("$infile-dedup.bim", "$infile-dedup.bim.bak");
		my $fin = IO::File->new("$infile-dedup.bim.bak");
		my $fout = IO::File->new("$infile-dedup.bim", "w");
		while(<$fin>) {
			my ($chr, $mkid, $gpos, $pos, $a1, $a2) = split;
			if ($mkid =~ /\.\d+$/) {
				$mkid =~ s/\.\d+$//;
			}
			print $fout join("\t", $chr, $mkid, $gpos, $pos, $a1, $a2), "\n";
		}
		$fout->close;
		return "$infile-dedup";
	}
	else {
		return "$infile-snps";
	}
}

# Align multiple bped files
# flip strand to align with the first file and keep variants that are common to all
# Return the outfile prefix, defaults to first file prefix + "-merged".
sub merge_bpeds {
	my ($infiles, $noflip) = @_;
	my $refbped = shift @$infiles;
	my @infiles = @$infiles;
	my %refsnp;
	{
		open my $fin, "$refbped.bim" or croak "Cannot open ref bimfile for reading";
		while(<$fin>) {
			my ($chr, $mkid, $gpos, $pos, $a1, $a2) = split;
			#croak "Incorrect marker ID: $_" unless $mkid eq "$chr-$pos-$a1-$a2" || $mkid eq "$chr-$pos-$a2-$a1" ;
			if (defined $refsnp{"$chr-$pos"}) {
				croak "Multiple variants found at $chr-$pos in $refbped";
			}
			$refsnp{"$chr-$pos"} = join("-", sort($a1, $a2));
		}
	}
	my %varcts;
	foreach my $infile (@infiles) {
		if ($noflip) {
			my $fin = IO::File->new("$infile.bim");
			while(<$fin>) {
				my ($chr, $mkid, $gpos, $pos, $a1, $a2) = split;
				my $a1a2 = join("-", sort($a1, $a2));
				if (defined $refsnp{"$chr-$pos"} && $refsnp{"$chr-$pos"} eq $a1a2 ) {
					$varcts{"$chr-$pos-$a1a2"} ++;
				}
			}
		}
		else {
			copy("$infile.bim", "$infile.bim.bak");
			my $fin = IO::File->new("$infile.bim.bak");
			my $fout = IO::File->new("$infile.bim", "w");
			while(<$fin>) {
				my ($chr, $mkid, $gpos, $pos, $a1, $a2) = split;
				#croak "Incorrect marker ID: $_" unless $mkid eq "$chr-$pos-$a1-$a2" || $mkid eq "$chr-$pos-$a2-$a1";
				my $a1a2 = join("-", sort($a1, $a2));
				my $a1a2rc = join("-", sort map { rev_comp($_) } ($a1, $a2));
				#my $a1a2rc2 = join("-", map { rev_comp($_) } reverse sort($a1, $a2));
				if (defined $refsnp{"$chr-$pos"}) {
					if ($refsnp{"$chr-$pos"} eq $a1a2) {
						$varcts{"$chr-$pos-$a1a2"} ++;
					}
					elsif ($refsnp{"$chr-$pos"} eq $a1a2rc) {
						$a1 = rev_comp($a1);
						$a2 = rev_comp($a2);
						$mkid = join('-', $chr, $pos, sort($a1, $a2));
						$varcts{$mkid} ++;
					}
				}
				print $fout join("\t", $chr, $mkid, $gpos, $pos, $a1, $a2), "\n";
			}
		}
	}
	# After going through each infile, generate a list of common SNPs to keep
	my $indir = dirname($refbped);
	open my $fcom, ">", "$indir/common_snps.extract" or croak "Cannot open common_snps";
	foreach my $mkid (keys %varcts) {
		if ($varcts{$mkid} == @infiles) {
			print $fcom $mkid, "\n";
		}
	}
	$fcom->close;

	foreach my $infile ($refbped, @infiles) {
		my $cmd = "$plink --bfile $infile --extract $indir/common_snps.extract --make-bed --out $infile-align";
		system($cmd); #or croak "merge_bpeds: plink runtime error when processing $infile";
	}
	open my $flst, ">", "$indir/merge_list.txt" or croak "Cannot write to merge_list";
	print $flst join("\n", map { "$_-align" } @infiles), "\n";
	$flst->close;

	my $cmd = "$plink --bfile $refbped-align --merge-list $indir/merge_list.txt --make-bed --out $refbped-merged";
	system($cmd); # or croak "merge_bpeds: plink runtime error in merging";

	return "$refbped-merged";
}


# Liftover processed bped
# Return outfile prefix, which should be "infile-targetdb"
sub liftover_bped {
	my ($infile, $chain) = @_;
	my $indir = dirname($infile);

	my $lift = Genome::UCSC::Liftover->new($chain);
	my ($fromdb, $Todb) = ($chain =~ /(hg\d+)To(Hg\d+)/);
	my $todb = lcfirst($Todb);

	# Liftover chr/pos and xclude unmapped SNPs
	open my $fin, "$infile.bim" or croak "Cannot open $infile.bim";
	open my $fexl, ">", "$infile-$todb.exclude" or croak "Cannot open $infile-$todb.exclude";

	my %bimline;
	while(<$fin>) {
		my ($chr, $mkid, $gpos, $pos, $a1, $a2) = split;
		my $chrom = hg_chrom($chr);
		my ($tg_chrom, $tg_pos, $tg_strand, $score) = $lift->query($chrom, $pos);
		my $tg_chr;
		if (defined $tg_chrom) {
			$tg_chr = hg_chr($tg_chrom, $tg_pos, $todb);
			my ($tg_a1, $tg_a2);
			if ($tg_strand eq '-') {
				$tg_a1 = rev_comp($a1);
				$tg_a2 = rev_comp($a2);
			}
			else {
				$tg_a1 = $a1;
				$tg_a2 = $a2;
			}
			# determine marker ID
			my $varid = join('-', $tg_chr, $tg_pos, sort($tg_a1, $tg_a2));
			$bimline{$mkid} = join("\t", $tg_chr, $varid, 0, $tg_pos, $tg_a1, $tg_a2); 
		}
		else {
			print $fexl $mkid, "\n";
		}
	}
	my $cmd = "$plink --bfile $infile --exclude $infile-$todb.exclude --make-bed --out $infile-$todb";
	system($cmd); #or croak "liftover_bped: plink runtime error";

	# Then patch bim file position and alleles
	copy("$infile-$todb.bim", "$infile-$todb.bim.bak");
	open my $fmap, "$infile-$todb.bim.bak" or croak "Cannot open $infile-$todb.bim.bak";
	open my $fout, ">", "$infile-$todb.bim" or croak "Cannot open $infile-$todb.bim";
	while(<$fmap>) {
		my $mkid = (split)[1];
		croak "Cannot find bimline for $mkid" unless defined $bimline{$mkid};
		print $fout $bimline{$mkid}, "\n";
	}
	$fout->close;
	# Then de-duplicate the resulting file.
	return dedup_bped("$infile-$todb");
} 

# Align alleles to reference sequence
sub refalign_bped {
	my ($infile, $seq) = @_;
	my ($todb) = ($infile =~ /\-(\w+)$/);
	croak "Cannot determine todb from $infile" unless defined $todb;
	#if ($seq !~ /$todb/i) {
	#	carp "Ref seq file name ($seq) does not match $todb?";
	#}

	# Verify todb
	my $sq = Genome::UCSC::TwoBit->new($seq);

	copy("$infile.bim", "$infile.bim.bak");
	open my $fin, "$infile.bim.bak" or croak "Cannot open $infile.bim.bak";
	open my $fout, ">", "$infile.bim" or croak "Cannot open $infile.bim";
	my @exclude;
	while(<$fin>) {
		my ($chr, $mkid, $gpos, $pos, $a1, $a2) = split;
		my $chrom = hg_chrom($chr);
		my $refbase = $sq->get_base($chrom, $pos);
		my $varid;
		if ($refbase eq $a1) {
			$varid = "$chr-$pos-$refbase-$a2";
		}
		elsif ($refbase eq $a2) {
			$varid = "$chr-$pos-$refbase-$a1";
		}
		else {
			if ($refbase eq rev_comp($a1)) {
				$varid = "$chr-$pos-$refbase-".rev_comp($a2);
				($a1, $a2) = (rev_comp($a1), rev_comp($a2));
			}
			elsif ($refbase eq rev_comp($a2)) {
				$varid = "$chr-$pos-$refbase-".rev_comp($a1);
				($a1, $a2) = (rev_comp($a1), rev_comp($a2));
			}
			else {
				$varid = $mkid;
				push @exclude, $varid;
			}
		}
		print $fout join("\t", $chr, $varid, $gpos, $pos, $a1, $a2), "\n";
	}
	if (@exclude) {
		open my $fexl, ">", "$infile-refalign.exclude" or croak "Cannot open $infile-refalign.exclude";
		print $fexl join("\n", @exclude), "\n";
		my $cmd = "$plink --bfile $infile --exclude $infile-refalign.exclude --make-bed --out $infile-refalign";
		system($cmd); # or croak "refalign_bped: plink runtime error";

		return "$infile-refalign";
	}
	else {
		return $infile;
	}
}


__END__

=head1 NAME

merge_bped -- Merge multiple bped files in a standardized way.

=head1 DESCRIOTION

We assume that samples from different files are distinct. (i.e., we are not merging
samples genotyped multiple times)

The following standard operations will be applied to each bped:

* keep only asymmetric biallelic SNPs, including symmetric SNPs if specified.
(Asym: A/C <=> T/G, A/G <=> T/C, Sym: A/T, C/G)

* name of SNPs to CHR-POS-sort(A1-A2) (where A1/A2 is based on the
first file or CHR-POS-REF-ALT, where REF/ALT is based on seq file)

* remove duplicated SNPs (keep one with higher genotyping rate)

* [optional] add prefix to sample names and family IDs

Then, merge multiple genotype files, this step will keep only SNPs appearing
in all files and align all SNP alleles to the first file. There is an option for 
asymmetric SNPs to be flipped or not. If only one file is provided, simply continue 
with following optional steps:

* liftOver genomic coordinates

* align SNP alleles to reference sequence

Currently, it only support SNP genotypes, other types of variants are ignored.

=head1 USAGE

merge_bped [options] PREFIX1 PREFIX2 ...

=head1 REQUIRED ARGUMENTS

=over

=item -[-]out[put] [=] <prefix>

Output file prefix

=item -[-]in[put] [=] <infiles>...

Input files (prefix only).

=back

=head1 OPTIONS

=over 

=item -[-]keep-sym

Keep symmetric SNPs, assuming they are on the same strand.
Use this option when only when you sure the strand for those SNPs are correct.

=item -[-]chain [=] <file>

Chain file for liftover. Liftover will be done for merged file.
Must be a UCSC chain file.

=for Euclid:
  file.type: readable

=item -[-]seq [=] <file>

Reference genome sequence file (.2bit). When chain file is provided
This will be the target genome reference. This will be applied to
the merged file to flip allele and adjust marker names (useful for imputation).
This step can be time comsuming.

=for Euclid:
  file.type: readable

=item -[-]processed

Assuming input genotype files have all been processed.

=item -[-]no-flip

Do not flip alleles for asymmetric SNPs when merging.

=item -[-]tmpdir [=] <dirname>

Provide a specified tempdir name. Otherwise, use automatically created tempdir
which will be destoryed after program exits.

=item -[-]tag[s] [=] <labels>

Commad separated list labels to tag different input files.
They will be used to add prefix to sample names and family IDs.

=item -[-]tag-sep [=] <string>

=for Euclid:
  string.default: ':'

=item -[-]tag-fid

Only tag FIDs.

=item -[-]tag-iid

Only tag IIDs.

=item -[-]plink [=] <path>

Use a custom version of plink

=back

=cut