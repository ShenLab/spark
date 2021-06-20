#!/usr/bin/env perl

use strict;
use warnings;
use Carp;
use FindBin qw|$Bin|;
use POSIX qw|ceil|;
use Getopt::Lucid qw|:all|;
use IO::Prompt;
use Data::Dumper;
use Config::Std;
use Cwd qw|abs_path|;
use List::MoreUtils qw|all any uniq|;
use Perl6::Slurp;
use String::ShellQuote;
use File::Path qw|make_path|;
use Utils::Workflow;
use Genome::UCSC qw|%PAR|;
use Genet::Ped;
use Utils::Hash qw|merge_conf|;
use Genome::Ranges::IntSet;


use lib "$Bin/../lib/";
use Shared qw|read_list parse_tabfile|;


#############################
## Command line interface  ##
#############################
# The default input is plink bped format genotypes and a output directory
my @spec = (
	Param("geno|i")->valid(sub { -r "$_.bed" && -r "$_.bim" && -f "$_.fam" }),
	Param("outdir|out"),
	Param("conf")->valid(sub { -r }),
	Param("engine|eng")->valid(sub { $_ eq 'SGE' || $_ eq 'BASH' }),
	Keypair("param|par"),
	Switch("help|h"),
	Switch("dryrun|dry"),
	Switch("force")
	);


my $opt = Getopt::Lucid->getopt(\@spec);

if ($opt->get_help) {
	print STDERR <<EOF;
Purpose:
	This is a pipeline script to call uniparental disomy (UPD) from SNP genotypes.

Usage: 
	upd_test.pl --conf CONF --geno PREFIX --out OUTDIR 

Input/output: 
	The script uses UPDio on trio-based SNP genotypes to infer UPD events. The input to this script should 
	be QC-passed genotypes in plink binary PED format. The script will extract all available trios from the 
	genotype data, then apply UPDio to test if genotypes informative of UPD is significantly enriched in 
	each offspring chromosome. Significant UPD events (by pre-defined p-value cutoff) and inferred UPD type 
	will be written to UPD_Events.txt. Then for each UPD event, the longest segment from the chromsome that 
	contain most of informative SNPs (by pre-defined fraction) will be called as UPD segments in UPD_Segs.txt.
	Plots supporting each UPD events can be found in the wrk directory.

	Although segments that are close to chromosome length are most likely UPD events and smaller segments
	are likely copy number deletions, using SNP genotyoes alone cannot distinguish copy number deletions from 
	segmental UPD event. When CNV calls are available, providing regions of CNV calls can remove UPD calls 
	overlapping these regions. Alternatively, the smaller UPD calls can be annotated for LRR/BAF (cnv_lrrbaf.pl) 
	to test if they are consistent with copy number deletions.

Notes:
	UPDio currently only support calling UPD events on autosomes.

Dependencies:
	UDPio, R(ggplot2, quantsmooth)

EOF
	exit 1;
}

$opt->validate({requires => [qw|geno outdir|]});

my $rootdir = $opt->get_outdir;

my %conf = merge_conf($opt->get_conf, $opt->get_param); 
$conf{PATH}{GENO} = $opt->get_geno;
$conf{PATH}{MODULE} = shell_quote("$Bin/module");

############################
## Input files validation ##
############################

# Identify total number trios.
my @alltrios;
{
	my $ped = Genet::Ped->new($opt->get_geno.".fam");
	my $trios = $ped->get_trios();
	foreach my $fid (sort keys %$trios) {
		foreach my $triad (sort { $a->[0] cmp $b->[0] } @{$trios->{$fid}}) {
			push @alltrios, [$fid, @$triad];
		}
	}
}

unless (@alltrios) {
	print STDERR "No trio was found from the pedigree file";	
	exit 1;
}
else {
	my @allchd = uniq sort map { $_->[1] } @alltrios;
	unless (@allchd == @alltrios) {
		croak "The number of unique children in trios does not equal the number of trios";
	}
	else {
		print "Found a total number of ", scalar(@alltrios), " trios in ped file\n";
	}
}

# Slurp CNVs
my %allcnvs;
if (exists $conf{CNV}) {
	unless(exists $conf{CNV}{Table} && -f $conf{CNV}{Table}) {
		die "Cannot find CNV table: $conf{CNV}{Table}";
	}
	print "CNV data are available\n";
	my ($iter, $fnames, $infields) = 
		parse_tabfile($conf{CNV}{Table}, $conf{CNV}{Fields}, 4, 4);
	while(my $dat = $iter->()) {
		my ($iid, $chr, $start, $end) = @{$dat}{@$infields};
		$chr =~ s/^chr//;
		push @{$allcnvs{$iid}}, [$chr, $start, $end];
	}
}

#################################
##  Workflow initialization &  ##
##  Working directory setup    ##
#################################

my $wkf = Utils::Workflow->new($rootdir,
	{ strict_var => 1, engine => $opt->get_engine, force => $opt->get_force });

# We will ignore the union of CNV regions in parent offspring trios
for(my $ii = 0; $ii < @alltrios; $ii ++) {
	my $jj = $ii + 1;
	my ($fid, $child, $dad, $mom) = @{$alltrios[$ii]};
	#my $trioid = join($conf{PARAM}{IDSep}, $fid, $child, $dad, $mom);
	#push @{$alltrios[$ii]}, $trioid;
	make_path "$rootdir/wrk/$child";
	open my $fout, ">$rootdir/par/TRIO.$jj" or die "Cannot write to TRIO.$jj";
	print $fout join("\t", $fid, $child, $dad, $mom), "\n";
	if (defined $allcnvs{$child}) {
		open my $fbed, ">$rootdir/par/$child.cnvs.bed" or die "Cannot write to child CNVs";
		foreach my $dat (@{$allcnvs{$child}}) {
			print $fbed join("\t", @$dat), "\n";
		}
	}
}

$wkf->add(upd_test(),  { name => "UPDTest", nslots => scalar(@alltrios),
						 expect => [ upd_test_files() ] })
	->add(collect_res(), { name => "CollectRes", depend => "UPDTest", 
						 expect => ["out/UPD_Events.txt", "out/UPD_Segs.txt"] });

$wkf->inst(\%conf);
write_config %conf, "$rootdir/par/run.conf" unless $opt->get_dryrun;

$wkf->run({ conf => $conf{$wkf->{engine}}, dryrun => $opt->get_dryrun });


############################
## Workflow components    ##
############################

sub upd_test {
	my $script = <<'EOF';

read FID CHILD DAD MOM < _PARDIR_/TRIO._INDEX_

echo -e "$FID $CHILD\n$FID $DAD\n$FID $MOM" > _WRKDIR_/$CHILD/samps.txt

plink --bfile _PATH.GENO_ --keep _WRKDIR_/$CHILD/samps.txt \
	--recode vcf-iid --out _WRKDIR_/$CHILD/genos

if [[ -f _PARDIR_/$CHILD.cnvs.bed ]]; then
	CHILDCNV="--child_cnv_data _PARDIR_/$CHILD.cnvs.txt"
fi

perl _PATH.UPDIO_/UPDio.pl --multisample_vcf _WRKDIR_/$CHILD/genos.vcf \
	--childID $CHILD --dadID $DAD --momID $MOM --output_path _WRKDIR_/$CHILD \
	--path_to_R _PATH.RSTAT_ --significance_level _PARAM.SIGNIF_ \
	--include_MI $CHILDCNV

rm -f _WRKDIR_/_INDEX_/genos.vcf

EOF
	return $script;
}

sub upd_test_files {
	my @exps;
	foreach my $jj (1..@alltrios) {
		my $trio = $alltrios[$jj-1];
		my ($fid, $child, $dad, $mom) = @$trio;
		push @exps, [map { "wrk/$child/$child.$_" } qw|events_list table upd log|];
	}
	return @exps;
}

sub collect_res {
	my $script = <<'EOF';

echo -e "IID\tP-value\tType\tChr" > _OUTDIR_/UPD_Events.txt

find _WRKDIR_ -name '*.upd' | \
	while read outfile; do 
		dir=$(dirname $outfile)
		if [[ -s $outfile ]]; then 
			fbase=$(basename $outfile)
			awk -v ID=${fbase%.upd} 'BEGIN{OFS="\t"}{print ID, $1, $2, $3}' $outfile
		fi
done >> _OUTDIR_/UPD_Events.txt

# Then extract UPD segments from events list
perl _PATH.MODULE_/extract_upd_segs.pl --in _OUTDIR_/UPD_Events.txt --out _OUTDIR_/UPD_Segs.txt \
	--wrk _WRKDIR_ --chromsize _PATH.CHROMSIZE_

EOF
	return $script;
}


