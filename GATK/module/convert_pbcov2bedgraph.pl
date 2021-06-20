use strict;
use warnings;
use File::Which;
use Utils::File qw|open_file|;

# Convert per-base coverage to BED graph 

my ($prefix) = @ARGV;
$prefix =~ s/\.txt\.gz$// if $prefix =~ /\.txt\.gz$/;

my $input = "$prefix.txt.gz";
my $output = "$prefix.bedgraph";

my $fin = open_file($input);
#my $fbed = open_file($output, {write => 1});
print STDERR "Converting per-base coverage to bedGraph\n";
open my $fbed, "| sort -k1,1 -k2,2n > $prefix.bedgraph" or die "Cannot write to output pipe for sort";

my %maxpos;
my ($prevchr, $preval, @prevpos) = ("");
while(<$fin>) {
	chomp;
	my ($chr, $pos, $val) = split /\t/;
	if ($chr ne $prevchr || $val != $preval || $pos != $prevpos[-1]+1) {
		print $fbed join("\t", $prevchr, $prevpos[0]-1, $prevpos[-1], $preval), "\n" if @prevpos;
		$maxpos{$prevchr} = $prevpos[-1]+5000 if @prevpos;
		@prevpos = ();
	}
	$prevchr = $chr;
	$preval = $val;
	push @prevpos, $pos;
}
print $fbed join("\t", $prevchr, $prevpos[0]-1, $prevpos[-1], $preval), "\n";
close $fbed;
$maxpos{$prevchr} = $prevpos[-1]+5000;

#system(qq|bedGraphToBigWig $output /tmp/b37.chrom.sizes BAMQC/Callable/perbase_cov.bw|)
open my $fsz, ">$prefix.chrom.sizes.txt" or die "Cannot write to $prefix.chrom.sizes.txt";
foreach my $chr (sort keys %maxpos) {
	print $fsz $chr, "\t", $maxpos{$chr}+5000, "\n";
}
close $fsz;

my $bgtobw = which("bedGraphToBigWig");
if ($bgtobw) {
	print STDERR "Converting bedGraph to bigWig\n";
	system("$bgtobw $prefix.bedgraph $prefix.chrom.sizes.txt $prefix.bw");
}

my $bigzip = which("bgzip");
my $tabix = which("tabix");
if ($bigzip && $tabix) {
	print STDERR "Compress and index bedGraph\n";
	system("bgzip $prefix.bedgraph");
	system("tabix -p bed $prefix.bedgraph.gz");
}

