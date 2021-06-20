use strict;
use warnings;
use Data::Dumper;
use Getopt::Euclid;
use File::Temp qw|tempdir|;
use List::MoreUtils qw|all|;
use Utils::Seq qw|rev_comp|;
use File::Copy qw|copy|;
use File::Path qw|make_path|;
use Genome::UCSC qw|hg_chr hg_chrom|;
use Genome::UCSC::Liftover;


unless(all { -f "$ARGV{'--input'}.$_" } qw|bim bed fam|) {
	die "Cannot find bim/bed/fam files"
}

my $lift = Genome::UCSC::Liftover->new($ARGV{'--chain'});
my ($fromdb, $Todb) = ($ARGV{'--chain'} =~ /(hg\d+)To(Hg\d+)/);
unless($fromdb && defined $Todb) {
	die "Cannot parse fromdb and Todb from chain file: $ARGV{'--chain'}!";
}
my $todb = lcfirst($Todb);

my $wrkdir;
if ($ARGV{'--wrkdir'}) {
	$wrkdir = $ARGV{'--wrkdir'};
	make_path $wrkdir unless -d $wrkdir;
}
else {
	$wrkdir = tempdir(CLEANUP => 1);
}

open my $fin, "$ARGV{'--input'}.bim" or die "Cannot open $ARGV{'--input'}.bim";
open my $fexl, ">", "$wrkdir/plink.exclude" or die "Cannot write to plink.exclude";

my %bimline;
while(<$fin>) {
	my ($chr, $mkid, $gpos, $pos, $a1, $a2) = split;
	my $chrom = hg_chrom($chr);
	my ($tg_chrom, $tg_pos, $tg_strand, $score) = $lift->query($chrom, $pos);
	my $tg_chr;
	if (defined $tg_chrom) {
		$tg_chr = hg_chr($tg_chrom, $tg_pos, $todb);
		if ($tg_chr eq '0') {
			print $fexl $mkid, "\n";
		}
		else {
			my ($tg_a1, $tg_a2);
			if ($tg_strand eq '-') {
				$tg_a1 = rev_comp($a1);
				$tg_a2 = rev_comp($a2);
			}
			else {
				$tg_a1 = $a1;
				$tg_a2 = $a2;
			}
			if(defined $bimline{$mkid}) {
				die "Marker ID $mkid is duplicated!";
			}
			$bimline{$mkid} = join("\t", $tg_chr, $mkid, 0, $tg_pos, $tg_a1, $tg_a2);
		}
	}
	else {
		print $fexl $mkid, "\n";
	}
}

my $cmd = "plink --bfile $ARGV{'--input'} --exclude $wrkdir/plink.exclude --make-bed --out $wrkdir/plink";
system($cmd);

copy("$wrkdir/plink.bim", "$wrkdir/plink.bim.bak");
open my $fmap, "$wrkdir/plink.bim.bak" or die "Cannot open $wrkdir/plink.bim.bak";
open my $fout, ">$wrkdir/plink.bim" or die "Cannot write to $wrkdir/plink.bim";
while(<$fmap>) {
	my $mkid = (split)[1];
	die "Cannot find bimline for $mkid" unless defined $bimline{$mkid};
	print $fout $bimline{$mkid}, "\n";
}
close $fout;

system(qq|plink --bfile $wrkdir/plink --make-bed --out $ARGV{'--output'}|);


__END__

=head1 NAME

liftover_bped.pl -- Liftover binary PED file.

=head1 NOTE

SNP map file (.bim) does not regard alleles, So only position will be lifted and reference base
check should be done later. If target position is on the reverse strand, alleles will be flipped.
Liftover is not one-to-one mapping, it is also likely that different markers will be mapped to the
same location after liftover.

=head1 REQUIRED ARGUMENTS

=over 

=item -[-]in[put] [=] <prefix>

Input genotype data file prefix. Marker ID must be unique!

=item -[-]chain [=] <file>

Liftover chain file. Must be UCSC's chain file with chr name prefix. 

=for Euclid:
  file.type: readable

=item -[-]out[put] [=] <prefix>

Output file name prefix. 

=back

=head1 OPTIONS

=over

=item -[-]wrkdir [=] <dir>

Working directory.

=back

=cut
