#!/usr/bin/env perl
use strict;
use warnings;
use Carp;
use IO::Dir;
use IO::File;
use Data::Dumper;
use File::Path qw|make_path|;
use File::Basename;
use List::MoreUtils qw|all none any|;
use Data::Dumper;
use Getopt::Euclid;
use Perl6::Slurp;
use Genet::File::VCF;
use Genome::UCSC qw|hg_chrom|;
use Genome::UCSC::BinKeeper;


my $vcf = Genet::File::VCF->new($ARGV{'--vcf'});

my @samps;
if ($ARGV{'--samp'}) {
	if (-f $ARGV{'--samp'}) {
		my %keep = map { (split)[0] => 1 } slurp $ARGV{'--samp'};
		print STDERR "Reading sample list from $ARGV{'--samp'}\n";
		@samps = grep { defined $keep{$_} } $vcf->get_sampids();
	}
	else {
		@samps = split(',', $ARGV{'--samp'});
		my @allsamps = $vcf->get_sampids();
		foreach my $samp (@samps) {
			if (none { $_ eq $samp } @allsamps) {
				carp "Cannot find sample $samp in VCF!";
			}
		}
	}
}
else {
	@samps = $vcf->get_sampids();
}
unless (@samps) {
	print STDERR "No sample left require processing\n";
	exit 1;
}


my $outdir = $ARGV{'--outdir'};
make_path $outdir unless -d $outdir;
$outdir .= '/' if $outdir !~ /\/$/;


my @filters = qw|PASS|;
if ($ARGV{'--filters'}) {
	push @filters, split(',', $ARGV{'--filters'});
}
else {
	push @filters, grep { /SNP/ } $vcf->list_vqsr_filter({SNP => $ARGV{'--tranche'}});
}
print STDERR "The SNPs with following filter flags will be included:\n";
print STDERR join("\t", @filters), "\n";

my $targets = Genome::UCSC::BinKeeper->new();
if ($ARGV{'--range'}) {
	my $fin = IO::File->new($ARGV{'--range'}) or croak "Cannot read range file";
	while(<$fin>) {
		my ($chr, $st, $ed, $name) = (split)[0,1,2,3];
		$targets->add($chr, [$st, $ed, $name]);
	}
}

my $wrkdir = $outdir;
{
	my $fs = IO::File->new("$wrkdir/tally_abhist.R", "w");
	print $fs <<EOF;
library(mclust)

args = commandArgs(TRUE)
infile = args[1]
tag = args[2]

X<-read.table(infile, header=T)

# Fit normal mixture models
ClstBest<-Mclust(X\$AB, modelNames = c("V"), G = 1:$ARGV{'--cmax'})

# Save data structure
save(ClstBest, file=paste("$outdir",tag,".model.dat",sep=""))

# Output model parameters
sink(paste("$outdir",tag,".best_model.summary",sep=""))
cat("#NCOMP=", ClstBest\$G, sep="", fill=TRUE)
cat("PROP", "MEAN", "VARIANCE", fill=TRUE)
if(ClstBest\$G>1) {
  for(jj in 1:ClstBest\$G) {
     cat(ClstBest\$parameter\$pro[jj],
         ClstBest\$parameter\$mean[jj],
         ClstBest\$parameter\$variance\$sigmasq[jj], fill=TRUE)
  }
} else {
  cat(ClstBest\$parameter\$pro, ClstBest\$parameter\$mean,
      ClstBest\$parameter\$variance\$sigmasq, fill=TRUE)
}
sink()

# Calculate sample statistics
sink(paste("$outdir",tag,".stats.summary",sep=""))
if(require(e1071)) {
	cat("MEAN", "MEDIAN", "STDEV", "VAR", "SKEWNESS", "KURTOSIS", fill=TRUE)
	cat(mean(X\$AB), median(X\$AB), sd(X\$AB), var(X\$AB),
    	skewness(X\$AB), kurtosis(X\$AB), fill=TRUE)
} else {
	cat("MEAN", "MEDIAN", "STDEV", "VAR", fill=TRUE)
	cat(mean(X\$AB), median(X\$AB), sd(X\$AB), var(X\$AB), fill=TRUE)
}
sink()

# Make BIC plot
png(paste("$outdir",tag,".bic.png",sep=""), width=5, height=5, units="in", res=300)
plot(ClstBest, what="BIC")
dev.off()

# Estimate mixture density
DensityBest<-densityMclust(X\$AB, modelNames=c("V"), G=c(ClstBest\$G))

# Make density+histogram plot
png(paste("$outdir",tag,".hist.png",sep=""), width=5, height=5, units="in", res=300)
plot(DensityBest, what="density", X\$AB, breaks=50, hist.col="lightblue", hist.border="white",
     xlab="Allele Balance")
dev.off()

# Make CDF and QQ plot
png(paste("$outdir",tag,".qq.png",sep=""), width=5, height=5, units="in", res=300)
#densityMclust.diagnostic(DensityBest, X\$AB, what="qq")
plot(DensityBest, what="diagnostic", type="qq")
dev.off()

png(paste("$outdir",tag,".cdf.png",sep=""), width=5, height=5, units="in", res=300)
#densityMclust.diagnostic(DensityBest, X\$AB, what="cdf", col=c("red", "blue"), legend=FALSE)
#legend(x="bottomright", legend=c("Est CDF", "Emp CDF"),
# col=c("red", "blue"), lty=c("solid", "dashed"), cex=0.75)
plot(DensityBest, what="diagnostic", type="cdf", col=c("red", "blue"))
dev.off()

# Perform KS test?
pfit<-function(x){
    y=0
    for(ii in 1:ClstBest\$G){
        y=y+ClstBest\$parameters\$pro[ii]*pnorm(x,ClstBest\$parameters\$mean[ii],sqrt(ClstBest\$parameters\$variance\$sigmasq[ii]))
    }
    return(y)
}

if (ClstBest\$G>1) {
  sink(paste("$outdir",tag,".kstest.summary",sep=""))
  print(ks.test(X\$AB, "pfit"))
  sink()
}
EOF
	$fs->close();
}


my (%varcounts, %fout);

foreach my $iid (@samps) {
	$fout{$iid} = IO::File->new("$outdir$iid.ab.txt", "w");
}
foreach my $iid (keys %fout) {
	$fout{$iid}->print(join("\t", qw|CHR POS REF ALT AB DEPTH|), "\n");
}

my $it = $vcf->iter({ strict => 0, samp => \@samps });
while(my $dat = $it->()) {
	my ($site, $geno) = @$dat;
	next unless defined $site->{ALT};
	if ($ARGV{'--range'}) {
	  next unless $targets->any_overlap($site->{CHROM}, $site->{POS}, $site->{POS});
	}

	unless($ARGV{'--no-rsid'}) {
		next unless defined $site->{ID};
	}
	
	next if hg_chrom($site->{CHROM}) eq 'chrX';

	# We will only consider biallelic SNPs
	for(my $ii = 0; $ii < @{$site->{ALT}}; $ii ++) {
		my $jj = $ii + 1;

		next unless $site->{REF} =~ /^[ACGT]$/ && $site->{ALT}[$ii] =~ /^[ACGT]$/;
		next unless any { $site->{FILTER}[0] // "." eq $_ } @filters;

		my $vid = join("\t", @{$site}{qw|CHROM POS REF|}, $site->{ALT}[$ii]);

		foreach my $iid (@samps) {
			my $gdat = $geno->{$iid};
			next unless defined $gdat->{GT};
			$varcounts{$iid}{HOM} ++ if $gdat->{GT} eq "$jj/$jj";
			next unless $gdat->{GT} eq "0/$jj" || $gdat->{GT} eq "$jj/0";
			$varcounts{$iid}{HET} ++;
			my $depth = $gdat->{AD}[0]+$gdat->{AD}[$jj];
			next unless $depth >= $ARGV{'--depth'};
			my $ab = $gdat->{AD}[0]/$depth;
			$fout{$iid}->print(join("\t", $vid, $ab, $depth), "\n");
		}
	}
}

foreach my $iid (keys %fout) {
	$fout{$iid}->close;
	print STDERR "\nPlotting $iid, NumHets=$varcounts{$iid}{HET}\n";
	system(qq|Rscript $wrkdir/tally_abhist.R $outdir$iid.ab.txt $iid|)
		unless $ARGV{'--dry'};
}


# Give the number of het/hom summary
my $fsum = IO::File->new("$outdir/varcounts.txt", "w");
print $fsum join("\t", qw|IID HET HOM|), "\n";
foreach my $iid (sort keys %varcounts) {
	print $fsum join("\t", $iid, $varcounts{$iid}{HET}, $varcounts{$iid}{HOM}), "\n";
}

__END__

=head1 NAME

ab_het.pl -- Allele balance of heterozygotes SNPs.

=head1 DESCRIPTION

Distribution of AB can be used to test cross-sample contamination.
This script will analyze all QC-passed SNPs with rs-IDs.

=head1 REQUIRED ARGUMENTS

=over

=item -[-]vcf [=] <file>

Input VCF file. Should be processed by standard GATK pipeline.
Known variants should be annotated. We will examine the AB of heterozygotes
at known biallelic SNPs to test allelic imbalance.

=for Euclid:
  file.type: readable

=item -[-]out[dir] [=] <outdir>

Output directory.

=back

=head1 OPTIONS

=over

=item -[-]samp[le] [=] <list>

A list of samples to be processed. Otherwise, will process all samples in the VCF (slow).
Also be careful: if all samples are processed, very large number of files needs to be opened.

=item -[-]range[s] [=] <intervals>

Restrict to a list of target intervals.

=item -[-]tranche [=] <threshold>

SNP tranche threshold. It will be used to derive filter tags to keep variants.
When parsing VCF files, only the first tag in FILTER fields will be used.

=for Euclid:
  threshold.default: 99.8

=item -[-]filter[s] [=] <tag>

Manually specify filters to keep sites (in addition to PASS), comma separated list.
It will override VQSR filters. Only the first tag in FILTER fields will be used.

=item -[-]depth [=] <depth>

Minimal depth of the het SNPs to be included in AB histogram.

=for Euclid:
  depth.default: 12

=item -[-]cmax [=] <ncomp>

The max. number of mixture Gaussian components to fit.

=for Euclid:
  ncomp.default: 10

=item -[-]no-rsid

Do not restrict the analysis to variants with rsIDs.

=item -[-]dry

Dry run mode.

=back
