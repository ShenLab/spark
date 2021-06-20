# Create summary of mutation types
use strict;
use warnings;
use FindBin qw|$Bin|;
use Getopt::Euclid;
use File::Basename;
use File::Temp qw|tempdir|;
#use Genome::UCSC::BinKeeper;
use Utils::Parser qw|sql_query|;

use lib "$Bin/../../lib";
use Shared qw|parse_tabfile|;
use Variants qw|var_type|;


my $wrkdir;
if ($ARGV{'--wrkdir'}) {
	$wrkdir = $ARGV{'--wrkdir'};
	mkdir $wrkdir unless -d $wrkdir;
}
else {
	$wrkdir = tempdir(CLEANUP => 1);
}


my (%vtcount, %ctbysamp);
my ($it, $fnames, $keyfields) = parse_tabfile($ARGV{'--input'}, $ARGV{'--fields'}, 4, 5);
my ($filter, $tokens);
if (defined $ARGV{'--filter'}) {
	($filter, $tokens) = sql_query($ARGV{'--filter'}, 1);
	foreach my $tok (@$tokens) {
		if ($tok->[0] eq 'FIELD') {
			unless(grep { $tok->[1] eq $_ } @$fnames) {
				die "Cannot find Filter field $tok->[1] in variant table $$ARGV{'--input'}";
			}
		}
	}
}

my $iidflag = @$keyfields == 5 ? 1 : 0;
{
	my %known;
	# Write SNVs into a temp VCF for plot
	my $label = (split(q|\.|, basename($ARGV{'--input'})))[0];
	open my $fvcf, ">$wrkdir/SNVs.vcf" or die "Cannot write SNVs.vcf";
	print $fvcf "##fileformat=VCFv4.1\n";
	print $fvcf '#'.join("\t", qw|CHROM POS ID REF ALT QUAL FILTER INFO|), "\n";
	# Write indel length to a data table
	open my $find, ">$wrkdir/IndelLens.txt" or die "Cannot write to Indels.txt";
	print $find join("\t", qw|CHROM POS REF ALT TYPE LEN|), "\n";
	while(my $dat = $it->()) {
		if (defined $filter) {
			next unless $filter->($dat);
		}
		my ($iid, $chrom, $pos, $ref, $alt);
		if ($iidflag) {
			($iid, $chrom, $pos, $ref, $alt) = @{$dat}{@$keyfields};
		}
		else {
			($chrom, $pos, $ref, $alt) = @{$dat}{@$keyfields};
		}
		my $varid = join(":", $chrom, $pos, $ref, $alt);
		if ($iidflag) {
			next if defined $known{$iid,$varid};
		}
		else {
			next if defined $known{$varid};
		}
		
		my $vt = var_type($ref, $alt);
		$vtcount{$vt} ++;
		if ($ref =~ /^[ACGT]$/ && $alt =~ /^[ACGT]$/) {
			print $fvcf join("\t", $chrom, $pos, '.', $ref, $alt, 100, '.', '.'), "\n";
		}
		if ($vt eq 'Indel') {
			my ($type, $length);
			if (length($ref) > length($alt)) {
				$type = "Deletion";
				$length = length($ref)-length($alt);
			}
			else {
				$type = "Insertion";
				$length = length($alt) - length($ref);
			}
			print $find join("\t", $chrom, $pos, $ref, $alt, $type, $length), "\n";
		}

		if ($iidflag) {
			$known{$iid,$varid} = 1;
			$ctbysamp{$iid} ++;
		}
		else {
			$known{$varid} = 1;
		}
	}
	close $fvcf; close $find;

	if ($vtcount{SNV} > 0) {
		open my $fs, ">$wrkdir/SNVPattern.R" or die "Cannot write to SNVPattern.R";
		print $fs <<EOF;
library(MutationalPatterns)

ref_genome<-"BSgenome.Hsapiens.UCSC.$ARGV{'--hgbuild'}"
library(ref_genome, character.only=TRUE)

vcfs<-read_vcfs_as_granges("$wrkdir/SNVs.vcf", "$label", ref_genome)

png("$ARGV{'--output'}.type_occur.png", width=5, height=4, units="in", res=300)
type_occurences <- mut_type_occurrences(vcfs, ref_genome)
plot_spectrum(type_occurences, CT=TRUE)
dev.off()

png("$ARGV{'--output'}.mut_matrix.png", width=10, height=4, units="in", res=300)
mut_mat<-mut_matrix(vcfs, ref_genome)
plot_96_profile(mut_mat, condensed=TRUE)
dev.off()

EOF
	close $fs;
	system(qq|Rscript $wrkdir/SNVPattern.R|);
	}

	# Indel length histogram
	if ($vtcount{Indel} > 0) {
		open my $fs, ">$wrkdir/IndelLenHist.R" or die "Cannot write to IndelLenHist.R";
		print $fs <<EOF;
library(ggplot2)
library(dplyr)

Indels<-read.table("$wrkdir/IndelLens.txt", header=TRUE)

Indels\$LEN<-ifelse(Indels\$LEN>=$ARGV{'--maxlen'}, $ARGV{'--maxlen'}, Indels\$LEN)

LensCount<-count(Indels, TYPE, LEN)

for(Type in c("Deletion", "Insertion")) {
	for(Len in seq(1, $ARGV{'--maxlen'})) {
		if(nrow(filter(LensCount, LEN==Len & TYPE==Type))==0) {
			add_row(LensCount, TYPE=Type, LEN=Len, n=0)
		}
	}
}

for(Type in c("Deletion", "Insertion")) {
	for(Len in seq(1, $ARGV{'--maxlen'})) {
		if(nrow(filter(LensCount, LEN==Len & TYPE==Type))==0) {
			LensCount<-add_row(LensCount, TYPE=Type, LEN=Len, n=0)
		}
	}
}
LensCount\$LEN<-as.factor(LensCount\$LEN)

hist <- ggplot(LensCount, aes(x=LEN, y=n, fill=TYPE)) + 
		geom_bar(stat="identity", position="dodge") +
		labs(x="Length", y="Count", fill="Indel types") +
		scale_x_discrete(labels=c(seq(1,$ARGV{'--maxlen'}-1), ">=$ARGV{'--maxlen'}")) +
		theme_classic() 

ggsave("$ARGV{'--output'}.indel_len.png", hist, width=5, height=4, units="in")

EOF
		close $fs;
		system(qq|Rscript $wrkdir/IndelLenHist.R|);
	}

	# Variant type pie chart
	if (keys %vtcount > 1) {
		open my $fout, ">$wrkdir/VarTypes.txt" or die "Cannot write to VarTypes.txt";
		print $fout join("\t", qw|Type Count|), "\n";
		foreach my $type (keys %vtcount) {
			print $fout $type, "\t", $vtcount{$type}, "\n";
		}
		close $fout;
		open my $fs, ">$wrkdir/VarTypeChart.R" or die "Cannot write to VarTypeChart.R";
		print $fs <<EOF;
library(ggplot2)

VT <- read.table("$wrkdir/VarTypes.txt", header=TRUE)
VT\$Perc<-100*VT\$Count/sum(VT\$Count)
VT\$Desc<-with(VT, paste0(Type, " (", Count, ", ", round(Perc, 1), "%)"))
pie <- ggplot(VT, aes(x="", y=Count, fill=Desc)) + 
		geom_bar(stat="identity", width=1) + 
		coord_polar("y", start=0) + 
		labs(x = NULL, y = NULL, fill = NULL) +
		scale_fill_discrete(name = "Variant types") +
		theme_classic() + 
		theme(axis.line = element_blank(),
        	axis.text = element_blank(), axis.ticks = element_blank())
ggsave("$ARGV{'--output'}.var_type.png", pie, width=6, height=4, units="in")

EOF
		close $fs;
		system(qq|Rscript $wrkdir/VarTypeChart.R|);
	}

	# Check if nsamp and iidflag are both positive
	if ($ARGV{'--nsamp'} && $iidflag) {
		my $nsampwvar = scalar(keys %ctbysamp);
		if ($ARGV{'--nsamp'} < $nsampwvar) {
			die "Incorrect of number of samples!";
		}
		open my $fout, ">$wrkdir/CountPerSamp.txt" or die "Cannot write to CountPerSamp.txt";
		print $fout join("\t", qw|IID Count|), "\n";
		while(my ($iid, $count) = each %ctbysamp) {
			print $fout $iid, "\t", $count, "\n";
		}
		for(my $ii = $nsampwvar; $ii < $ARGV{'--nsamp'}; $ii ++) {
			print $fout "Samp$ii\t0\n";
		}
		close $fout;
		open my $fs, ">$wrkdir/CountFitDistr.R" or die "Cannot write to CountFitDistr.R";
		print $fs <<EOF;
library(vcd)

Counts<-read.table("$wrkdir/CountPerSamp.txt", header=T)
gf<-goodfit(Counts\$Count, type="poisson")
png("$ARGV{'--output'}.count_persamp.png", width=5, height=5, units="in", res=300)
plot(gf)
dev.off()

sink("$ARGV{'--output'}.count_persamp.txt")
print(gf)
summary(gf)
sink()

EOF
		close $fs;
		system(qq|Rscript $wrkdir/CountFitDistr.R|);
	}
}



__END__

=head1 NAME

sum_mutype.pl -- Summarize de novo or rare mutations.

=head1 NOTES

This is a utility to summarize the mutation types observed in a variant table.

The following files will be generated:
1. Use mutational pattern to plot mutation patterns of SNVs
2. Pie chart for SNV,indel,MNV proportions
3. Barplot for indel length
4. Histogram for number of mutation per sample if IID is provided.

We will not assess the allelic balance of variants, which is the purpose of another utility.

=head1 REQUIRED ARGUMENTS

=over

=item -[-]in[put] [=] <table>

Input variant table.

=item -[-]out[put] [=] <prefix>

Output file name prefix.

=back

=head1 OPTIONS

=over 

=item -[-]fields [=] <fstr>

Key fields in the input table, default: IID,Chrom,Position,Ref,Alt.
Must be in the specified order, IID can be optional.

=for Euclid:
	fstr.default: "IID,Chrom,Position,Ref,Alt"

=item -[-]hgbuild [=] <hgbuild>

Human genome build version, default: hg19.

=for Euclid:
	hgbuild.default: "hg19"

=item -[-]filter [=] <expr>

Filter expression.

=item -[-]wrkdir [=] <dir>

Working directory, use temp dir if not provided.

=item -[-]nsamp [=] <number>

Number of samples.

=item -[-]maxlen [=] <length>

Max. indel length to be capped.

=for Euclid:
	length.default: 7

=back


