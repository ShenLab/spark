#!/usr/bin/env perl

use strict;
use warnings;
use Carp;
use IO::Dir;
use File::Basename;

# Collect BamQC statistics
my ($wrkdir, $outdir, $wgs) = @ARGV;
unless(@ARGV == 2 || @ARGV == 3) {
	print STDERR "$0 WRKDIR OUTDIR [WGS]\n";
	exit 1;
}
unless (-d $wrkdir && -d $outdir) {
	die "Cannot read working dir or output dir";
}


# Specification for statistics to be collected
my %spec = ( 
	"BamStats.alignment_summary_metrics" =>
		{
			PCT_PF_READS_ALIGNED => ["AL.PCTalnPF","The percentage of PF reads pairs that aligned to the reference sequence."],
			PCT_READS_ALIGNED_IN_PAIRS => ["AL.PCTalnPair", "The fraction of reads whose mate pair was also aligned to the genome."],
			PCT_CHIMERAS => ["AL.Chimeras", "The fraction of reads that map outside of a maximum insert size (100k) or that have the two ends mapping to different chromosomes."],
			MEAN_READ_LENGTH => ["AL.AveReadLen", "The mean read length of the set of reads examined."],
			PF_MISMATCH_RATE => ["AL.Mismatch", "The rate of bases mismatching the reference for all bases aligned to the reference sequence."]
		},
	"BamStats.insert_size_metrics" => 
		{ 
			MEDIAN_INSERT_SIZE => ["IS.Median", "The MEDIAN insert size of all paired end reads where both ends mapped to the same chromosome."],
			MEAN_INSERT_SIZE   => ["IS.Mean", "The mean insert size of the \"core\" (after removing outliers) of the distribution."],
			STANDARD_DEVIATION => ["IS.SD", "Standard deviation of insert sizes over the \"core\" of the distribution."]
	 	},
	"BamStats.hs_metrics" =>
	 	{ 
	 		PCT_PF_UQ_READS => ["AL.PCTUniqPF", "The fraction of reads that passed vendor's filter (PF) and not marked as duplicates from total reads."],
	 		MEAN_TARGET_COVERAGE => ["DP.Mean", "The mean coverage of the target region."],
	 		MEDIAN_TARGET_COVERAGE => ["DP.Median", "The median coverage of a target region."],
	 		PCT_SELECTED_BASES => ["HS.PCTSelected", "The fraction of aligned PF-bases located on or near a baited/targeted region"],
	 		PCT_OFF_BAIT => ["HS.PCTOffBait", "The fraction of aligned PF-bases that are mapped away from any baited/targeted region"],
	 		ON_BAIT_VS_SELECTED => ["HS.OnVsSelected", "The fraction of bases on or near baits/targets that are covered by baits/targets"],
	 		FOLD_ENRICHMENT => ["HS.FoldEnrichment", "The fold by which the baited/targeted region has been amplified above genomic background."],
	 		ZERO_CVG_TARGETS_PCT => ["DP.PCT0x", "The fraction of targets that did not reach coverage=1 over any base."],
	 		PCT_TARGET_BASES_1X => ["DP.PCT1x", "The fraction of all target bases achieving 1X or greater coverage."],
	 		PCT_TARGET_BASES_2X => ["DP.PCT2x", "The fraction of all target bases achieving 2X or greater coverage."],
	 		PCT_TARGET_BASES_10X => ["DP.PCT10x", "The fraction of all target bases achieving 10X or greater coverage."],
	 		PCT_TARGET_BASES_20X => ["DP.PCT20x", "The fraction of all target bases achieving 20X or greater coverage."],
	 		PCT_TARGET_BASES_30X => ["DP.PCT30x", "The fraction of all target bases achieving 30X or greater coverage."],
	 		PCT_TARGET_BASES_50X => ["DP.PCT50x", "The fraction of all target bases achieving 50X or greater coverage."],
	 		AT_DROPOUT => ["GC.ATDropout", "A measure of how undercovered <= 50% GC regions are relative to the mean."],
			GC_DROPOUT => ["GC.GCDropout", "A measure of how regions of high GC content (>= 50% GC) are undercovered relative to the mean coverage value."],
	 	},
	"BamStats.gc_bias.summary_metrics" =>
		{
			AT_DROPOUT => ["GC.ATDropout", "A measure of how undercovered <= 50% GC regions are relative to the mean."],
			GC_DROPOUT => ["GC.GCDropout", "A measure of how regions of high GC content (>= 50% GC) are undercovered relative to the mean coverage value."],
		},
	"BamStats.wgs_metrics" =>
		{
			MEAN_COVERAGE => ["DP.Mean", "The mean coverage in bases in target regions, after all filters are applied."],
			MEDIAN_COVERAGE => ["DP.Median", "The median coverage in target regions, after all filters are applied."],
			PCT_EXC_MAPQ => ['AL.PCTExclMQ', 'The fraction of aligned bases that were filtered out because they were in reads with low mapping quality'],
			PCT_EXC_DUPE => ['AL.PCTExclDup', 'The fraction of aligned bases that were filtered out because they were in reads marked as duplicates.'],
			PCT_EXC_TOTAL => ['AL.PCTExclTot', 'The total fraction of aligned bases excluded due to low MAPQ, duplication, without mate pair, low base quality, capped coverage, and secondary observation from an insert with overlapping reads.'],
			PCT_1X => ["DP.PCT1x",	"The fraction of bases that attained at least 1X sequence coverage in post-filtering bases."],
			PCT_5X => ["DP.PCT5x", "The fraction of bases that attained at least 5X sequence coverage in post-filtering bases."],
			PCT_10X => ["DP.PCT10x", "The fraction of bases that attained at least 10X sequence coverage in post-filtering bases."],
			PCT_15X => ["DP.PCT15x", "The fraction of bases that attained at least 15X sequence coverage in post-filtering bases."],
			PCT_20X => ["DP.PCT20x", "The fraction of bases that attained at least 20X sequence coverage in post-filtering bases."],
			PCT_25X	=> ["DP.PCT25x", "The fraction of bases that attained at least 25X sequence coverage in post-filtering bases."],
			PCT_30X => ["DP.PCT30x", "The fraction of bases that attained at least 30X sequence coverage in post-filtering bases."],
			PCT_50X => ["DP.PCT50x", "The fraction of bases that attained at least 50X sequence coverage in post-filtering bases."],
		},	
	"chrom_NormDP.txt" =>
		{
			NDP_chrX => ["NDP.chrX", "Normalized chrX depth = Mean chrX target depth / mean autosome target depth."],
			NDP_chrY => ["NDP.chrY", "Normalized chrY depth = Mean chrY target depth / mean autosome target depth."],
		},
	"VBID.selfSM" =>
		{
			FREEMIX => ["ID.FREEMIX", "Sequence-only estimate of contamination (0-1 scale)"],
		},
	"BamStats.cpcg_metrics" =>
		{
			OXIDATION_Q => ["OXOG.Q", "Phred-scaled oxoG error rate"],
		}
	);


my @files = qw|BamStats.alignment_summary_metrics BamStats.insert_size_metrics
			   chrom_NormDP.txt VBID.selfSM BamStats.cpcg_metrics|;
if ($wgs) {
	push @files, "BamStats.wgs_metrics", "BamStats.gc_bias.summary_metrics";
}
else {
	push @files, "BamStats.hs_metrics"; 
}

my (%stats, %depth);
my @iids = grep { !/^\./ && -d "$wrkdir/$_" } IO::Dir->new($wrkdir)->read();
foreach my $iid (@iids) {
	foreach my $file (@files) {
		my $dat = collect_stat("$wrkdir/$iid/$file");
		next unless $dat;
		while(my ($key, $val) = each %{$spec{$file}}) {
			$stats{$iid}{$val->[0]} = $dat->{$key};
		}
		if ($file eq "chrom_NormDP.txt") {
			$depth{$iid} = $dat;
		}
	}
}

# Output results and README
my @fields = sort { $a->[0] cmp $b->[0] } map { values %{$spec{$_}} } @files;
open my $frdm, ">$outdir/bam_stats.README.txt" or die "Cannot write README";
foreach my $info (@fields) {
	print $frdm join("\t", @$info), "\n";
}

open my $fout, ">$outdir/bam_stats.csv" or die "Cannot write bam stats csv";
print $fout join(",", "IID", map { $_->[0] } @fields), "\n";
foreach my $iid (@iids) {
	next unless defined $stats{$iid};
	print $fout join(",", $iid, map { $stats{$iid}{$_->[0]} // "NA" } @fields), "\n"; 
}

# Also output per-chromsome level normalized depth
open my $fdp, ">$outdir/chrom_depth.csv" or die "Cannot write depth csv file";
my @dpfields = map { "NDP_chr$_" } 1..22, "X", "Y";
unshift @dpfields, "DP_auto"; 
print $fdp join(",", "IID", @dpfields), "\n";
foreach my $iid (@iids) {
	next unless defined $depth{$iid};
	print $fdp join(",", $iid, @{$depth{$iid}}{@dpfields}), "\n";
}


sub collect_stat {
	my ($infile) = @_;
	if (basename($infile) =~ /^BamStats/) {
		if ($infile =~ /alignment_summary_metrics$/) {
			return collect_picard_stat($infile, 0, "PAIR");
		}
		elsif ($infile =~ /cpcg_metrics$/) {
			return collect_picard_stat($infile, 2, "CCG");
		}
		else {
			return collect_picard_stat($infile);
		}
	}
	else {
		return collect_other_stat($infile);
	}
}

sub collect_picard_stat {
	my ($infile, $ncol, $colval) = @_;
	return unless -f $infile;
	my (%data, @fields);
	open my $fin, $infile or die "Cannot read $infile";
	while(<$fin>) {
		if (/^\w+/) {
			chomp();
			@fields = split /\t/;
			last;
		}
	}
	while(<$fin>) {
		chomp();
		my @values = split /\t/;
		if (defined $ncol) {
			next unless $values[$ncol] eq $colval;
		}
		if (@values < @fields) {
			push @values, ("") x (@fields - @values);
		}
		@data{@fields} = @values;
		last;
	}
	return \%data;
}


sub collect_other_stat {
	my ($infile) = @_;
	return unless -f $infile;
	my %data;
	open my $fin, $infile or die "Cannot read $infile";
	$_= <$fin>;
	my @fields = split;
	$_ = <$fin>;
	my @values = split;
	@data{@fields} = @values;
	return \%data;
}
