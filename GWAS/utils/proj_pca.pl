use strict;
use warnings;
use Perl6::Slurp;
use String::ShellQuote;
use FindBin qw|$Bin|;
use File::Path qw|make_path|;
use Getopt::Euclid;
use File::Basename;
use List::MoreUtils qw|all|;

my $refpanel = $ARGV{'--ref'};
my $rootdir = $ARGV{'--outdir'};
my $input = $ARGV{'--input'};

# Check related files
if (-d $rootdir) {
	unless($ARGV{'--force'}) {
		die "The output directory $rootdir already exist!";
	}
}
foreach my $subdir (qw|wrk out par tmp|) {
	make_path "$rootdir/$subdir";
}

unless(-f "$ARGV{'--path'}/snpwt.$refpanel") {
	die "Cannot find SNP weights for $refpanel";
}
unless(all { "-f $input.$_" } qw|bim bed fam|) {
	die "Not all input genotype related files can be found";
}

# Parse reference panel file to get information about ancestral populations
my (@pops, $npop, %refpc);
{
	open my $fin, "$ARGV{'--path'}/snpwt.$refpanel" or die "Cannot open snpwt.$refpanel";
	<$fin>;
	my $popline = <$fin>;
	@pops = split(/\s+/, $popline);
	$npop = scalar(@pops);
	<$fin>;
	my $pcline = <$fin>;
	my @dat = split(/\s+/, $pcline);
	unless(scalar(@dat) == $npop ** 2) {
		die "Incorrect number of fields on PC line:\n$pcline\n";
	}
	for(my $ii = 0; $ii < @pops; $ii ++) {
		my @popdat = @dat[$ii*$npop..($ii+1)*$npop-1];
		unless($popdat[0] eq $pops[$ii]) {
			die "Inconsistent population label at PC line: $popdat[0]<>$pops[$ii]";
		}
		$refpc{$popdat[0]} = \@popdat;
	}
}


my $label;
if ($ARGV{'--label'}) {
	$label = $ARGV{'--label'};
}
else {
	$label = basename($input);
}

# Convert input file to eigenstrat format
if ($ARGV{'--keep'} || $ARGV{'--remove'}) {
	my %fids = map { (split)[1,0] } slurp "$input.fam";
	my $option;
	foreach my $action (qw|keep remove|) {
		if ($ARGV{"--$action"}) {
			$option .= " --$action $rootdir/tmp/$action.txt";
			open my $fout, ">$rootdir/tmp/$action.txt" or die "Cannot write to $action.txt";
			foreach my $iid (map { (split)[0] } slurp $ARGV{"--$action"}) {
				if (defined $fids{$iid}) {
					print $fout $fids{$iid}, "\t", $iid, "\n";
				}
			}
		}
	}
	system(qq|plink --bfile $input $option --make-bed --out $rootdir/wrk/$label|);
	unless(all { -f "$rootdir/wrk/$label.$_" } qw|bim bed fam|) {
		die "Cannot find all genotype related file after keep/remove";
	}

	open my $fpar, ">$rootdir/par/convertf.conf" or die "Cannot write to convertf.conf";
	print $fpar <<EOF;
genotypename:    $rootdir/wrk/$label.bed
snpname:         $rootdir/wrk/$label.bim
indivname:       $rootdir/wrk/$label.fam
outputformat:    EIGENSTRAT
genotypeoutname: $rootdir/wrk/$label.eigenstratgeno
snpoutname:      $rootdir/wrk/$label.snp
indivoutname:    $rootdir/wrk/$label.ind
familynames:     NO
EOF
	close $fpar;

}
else {
	open my $fpar, ">$rootdir/par/convertf.conf" or die "Cannot write to convertf.conf";
	print $fpar <<EOF;
genotypename:    $input.bed
snpname:         $input.bim
indivname:       $input.fam
outputformat:    EIGENSTRAT
genotypeoutname: $rootdir/wrk/$label.eigenstratgeno
snpoutname:      $rootdir/wrk/$label.snp
indivoutname:    $rootdir/wrk/$label.ind
familynames:     NO
EOF
	close $fpar;
}
system(qq|convertf -p $rootdir/par/convertf.conf|);
unless(all { "$rootdir/wrk/$label.$_" } qw|eigenstratgeno snp ind|) {
	die "Cannot find all output from convertf";
}

# Run SNP weight projection
{
	open my $fpar, ">$rootdir/par/inferanc.conf" or die "Cannot write to inferanc.conf";
	print $fpar <<EOF;
geno:   $rootdir/wrk/$label.eigenstratgeno
snp:    $rootdir/wrk/$label.snp
ind:    $rootdir/wrk/$label.ind
snpwt:  $ARGV{'--path'}/snpwt.$refpanel
predpcoutput:   $rootdir/wrk/${refpanel}_Proj.txt

EOF
	close $fpar;
	system(qq|inferanc -p $rootdir/par/inferanc.conf|);
	unless(-f "$rootdir/wrk/${refpanel}_Proj.txt") {
		die "Cannot find output from inferanc";
	}
}

# Reformat the output
{
	open my $fin, "$rootdir/wrk/${refpanel}_Proj.txt" or die "Cannot open ${refpanel}_Proj.txt";

	open my $fpc, ">$rootdir/out/$label.PCs.txt" or die "Cannot write to $label.PCs.txt";
	print $fpc join("\t", "Label", "IID", map { "PC_$_" } 1..($npop-1)), "\n";
	foreach my $pop (@pops) {
		print $fpc join("\t", $pop, @{$refpc{$pop}}), "\n";
	}

	open my $fout, ">$rootdir/out/$label.Ancestry.txt" or die "Cannot write to $label.Ancestry.txt";
	print $fout join("\t", "IID", "NumSNPs", (map { "PC_$_" } 1..($npop-1)), (map { "Prob$_" } @pops)), "\n";
	while(<$fin>) {
		my @dat = split;
		unless(scalar(@dat) == $npop*2+2) {
			die "Incorrect number of columns in the projection output file: $_";
		}
		my $iid = shift @dat;
		shift @dat;
		print $fpc join("\t", $label, $iid, @dat[1..($npop-1)]), "\n";
		print $fout join("\t", $iid, @dat), "\n";
	}
}

# Plot PCs
my $moduledir = shell_quote("$Bin/module");
system(qq|Rscript $moduledir/plot_pop_proj.R $rootdir/out/$label $ARGV{'--prob'}|);
unless(-f "$rootdir/out/$label.PCs1-".($npop-1).".png") {
	die "Cannot make PCA plot!";
}


__END__

=head1 NAME

proj_pca.pl -- Wrapper of SNPweights to project population samples to PCA space.

=head1 NOTE

It requires convertf from EigenSoft and inferanc from SNPweights.

=head1 REQUIRED ARGUMENTS

=over

=item -[-]in[put] [=] <prefix>

Input genotype data file prefix. It is required that input map file contain SNP rs IDs.

=item -[-]ref [=] [=] <label>

One of the reference panels to use: CO,AA,NA,EA

=item -[-]out[dir] [=] <dir>

The output directory.

=back

=head1 OPTIONS

=over

=item -[-]label [=] <string>

The input data label. Default will be the basename after stripping off suffix.

=item -[-]keep [=] <goodsamp>

A list of samples to be included in the analysis.

=item -[-]remove [=] <badsamp>

A list of samples to be remove in analysis.

=item -[-]path [=] <path>

Path to SNP weights files.

=for Euclid:
	path.default: "$ENV{HOME}/Software/SNPweights/wt"

=item -[-]prob [=] <cutoff>

Probability threshold for making ancestry prediction (only used in plot)

=for Euclid:
	cutoff.default: 0.9

=item -[-]force

=back
