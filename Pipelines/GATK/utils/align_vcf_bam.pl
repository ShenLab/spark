use strict;
use warnings;
use Cwd qw|abs_path|;
use FindBin qw|$Bin|;
use File::Copy qw|move|;
use File::Basename;
use Data::Dumper;
use Utils::Hash qw|chk_default|;
use List::MoreUtils qw|uniq all|;
use Getopt::Euclid;

use lib "$Bin/../../lib";
use Shared qw|read_list|;
use Wrapper qw|vcf_sampids bam_sampids|;

my ($inprefix, $insuffix);
$inprefix = $ARGV{'--vcf'};
if ($ARGV{'--vcf'} =~ /\.vcf$/) {
	$inprefix =~ s/\.vcf$//;
	$insuffix = "vcf";
}
elsif ($ARGV{'--vcf'} =~ /\.vcf\.gz$/) {
	$inprefix =~ s/\.vcf\.gz$//;
	$insuffix = "vcf.gz";
	unless(-f $ARGV{'--vcf'}.".tbi") {
		die "Cannot find index file for input VCF: $ARGV{'--vcf'}.tbi!";
	}
}
else {
	die "Incorrect name for input VCF: $ARGV{'--vcf'}";
}

my @vcfids = vcf_sampids($ARGV{'--vcf'});

my ($outprefix, $outsuffix);
if (defined $ARGV{'--output'}) {
	$outprefix = $ARGV{'--output'};
	if ($ARGV{'--output'} =~ /\.vcf$/) {
		$outprefix =~ s/\.vcf$//;
		unless ($insuffix eq 'vcf') {
			warn "Output file suffix is not the same as input file suffix!";
		}
		$outsuffix = $insuffix;
	}
	elsif ($ARGV{'--output'} =~ /\.vcf\.gz$/) {
		$outprefix =~ s/\.vcf\.gz$//;
		unless($insuffix eq 'vcf.gz') {
			warn "Output file suffix is not the same as input file suffix!";
		}
		$outsuffix = $insuffix;
	}
	else {
		die "Incorrect name for output VCF: $ARGV{'--output'}";
	}
}
else {
	$outprefix = $inprefix;
	$outsuffix = $insuffix;
}

# ReaID => SeqID, reverse: SeqID => RealID
my %idmap;
{
	my %bams;
	if (@{$ARGV{'--bams'}} == 1 && $ARGV{'--bams'}[0] !~ /\.(cram|bam)$/) {
		my $bams = read_list($ARGV{'--bams'}[0], { suffix => ['bam', 'cram'] });
		%bams = %$bams;
	}
	else {
		unless(all { /\.(cram|bam)$/ } @{$ARGV{'--bams'}}) {
			die "Not all BAM/CRAM files have correct file name suffix!";
		}
		foreach my $bamfile (@{$ARGV{'--bams'}}) {
			my $bambase = basename($bamfile);
			my $iid = $bambase; $iid =~ s/\.(cram|bam)$//;
			$bams{$iid} = $bamfile;
		}
	}

	# Each bam file can have only one ID!
	while(my ($realid, $bamfile) = each %bams) {
		my @iids = bam_sampids($bamfile);
		die "File $bamfile contains multiple samples!" unless @iids == 1;
		unless($ARGV{'--reverse'}) {
			chk_default(\%idmap, $realid, $iids[0]);
		}
		else {
			chk_default(\%idmap, $iids[0], $realid);	
		}
	}
}

unless(all { defined $idmap{$_} } @vcfids) {
	if ($ARGV{'--strict'}) {
		die "Not all VCF samples are assocciated with BAM/CRAM files!";
	}
}


# Prepare rename list
my @rename;
foreach my $vcfid (@vcfids) {
	if (defined $idmap{$vcfid} && $idmap{$vcfid} ne $vcfid) {
		push @rename, [$vcfid, $idmap{$vcfid}];
	}
}

if (@rename) {	
	open my $fout, ">$outprefix.renamed.txt" or die "Cannot write to $outprefix.renamed.txt";
	foreach my $pair (@rename) {
		print $fout join("\t", @$pair), "\n";
	}
	close $fout;

	my $oldinvcf;
	if (abs_path($outprefix) eq abs_path($inprefix)) {
		$oldinvcf = "$inprefix.orig.$insuffix";
		move($ARGV{'--vcf'}, $oldinvcf);
		if ($insuffix eq 'vcf.gz') {
			move($ARGV{'--vcf'}.".tbi", $oldinvcf.".tbi");
		}
	}
	else {
		$oldinvcf = $ARGV{'--vcf'};
	}

	my $code = system(qq|bcftools reheader -s $outprefix.renamed.txt -o $outprefix.$outsuffix $oldinvcf|);
	if ($code) {
		print STDERR "bcftool reheader for $oldinvcf failed!";
		unlink("$outprefix.$outsuffix");
		exit 1;
	}

	if ($outsuffix eq 'vcf.gz') {
		my $code = system(qq|tabix -p vcf $outprefix.$outsuffix|);
		if ($code) {
			print STDERR "tabix index of $outprefix.vcf.gz failed!\n";
			exit 1;
		}
	}
}
else {
	if ($outprefix eq $inprefix) {
		print STDERR "No sample needs rename, keep input file intact\n";
	}
	else {
		print STDERR "No sample needs rename, create symlink to input file\n";
		symlink abs_path($ARGV{'--vcf'}), "$outprefix.$outsuffix";
		if ($outsuffix eq 'vcf.gz') {
			symlink abs_path($ARGV{'--vcf'}.".tbi"), "$outprefix.$outsuffix.tbi";
		}
	}
}

__END__

=head1 NAME

align_vcf_bam.pl -- align sample IDs of VCF with BAM/CRAM files.

=head1 NOTES

Most pipeline tools look for sample ID of BAM/CRAM file from second column a cram list or from the 
basename of the file after removing suffix, which are usually different from sample IDs encoded in
the file header. Sample IDs in VCF headers that generated from the data processing pipeline are 
always renamed to the real ID. It causes discrepancies of sample IDs between VCF and BAMs. 

For applications that make use of IDs store in the VCF or BAM file headers to match samples, this script
will aligned IDs between VCF and BAM files. The default is to rename VCF header (typically smaller) to
keep them consistent with BAM file headers. After processing, rename VCF header again to keep consistent
with real IDs.

=head1 REQUIRED ARGUMENTS

=over

=item -[-]vcf [=] <file>

Input VCF file.

=for Euclid:
	file.type: readable

=item -[-]bam[s] [=] <list>...

Can be one or more bam files or a list/directory of BAM files.

=back

=head1 OPTIONS

=over

=item -[-]out[put] [=] <file>

Output VCF file, it must be different from input. 
If no rename is required, a symbolic link will be created to the input.
If no output file is provided, we will use input file prefix.

=item -[-]reverse

The default is to rename VCF sample to align with sample IDs in BAM header.
Reverse is to rename them to keep VCF sample to real sample IDs.

=item -[-]strict

Require all samples in VCF to be found in at BAM list.

=cut
