use strict;
use warnings;
use IO::Dir;
use IO::File;
use Perl6::Slurp;

# Collect results of DNenrich

# Following file will be created: 
# Enrich.txt -- enrichment of variants hits
# Contra.txt -- contrast between groups

my ($pardir, $wrkdir, $outdir) = @ARGV;

# Find all analysis files
my @analyses = grep { /^analysis\.\d+$/ } IO::Dir->new($pardir)->read();

unless(@analyses) {
	print STDERR "No analysis file was found!\n";
	exit 1;
}

# Also read gene name aliases
my %alias;
if (-f "$pardir/alias") {
	print STDERR "Reading gene name alias\n";
	%alias = map { (split)[0,1] } slurp "$pardir/alias";
}

# Output files
my %fields = ( Enrich => [qw|VarClass SampGroup GeneSet NGenesInSet BackGround Test 
						  	 NMuts NGenesInMuts ObsStat ExpStat FoldEnrich Pval Genes|],
			   Contra => [qw|VarClass SampGrp1 SampGrp2 GeneSet NGenesInSet BackGround Test 
				   		     NMuts NGenesInMuts ObsStat ExpStat FoldEnrich Pval GenesGrp1 GenesGrp2|] );

my %desc = (VarClass 	=> "Variant class or class combination", 
			SampGroup 	=> "Sample group",  SampGrp1 => "Sample group 1", SampGrp2 => "Sample group 2",
			GeneSet 	=> "Gene set used in enrichment test", 
			NGenesInSet => "Total number of genes in the set (restricting to genes in the background)",
			BackGround 	=> "Label of background gene set",
			Test => "The number of DNVs in gene for it to be included in the test (* denotes 1 or more for enrichment, 2 would denote a test of recurrence within that gene set)",
			NMuts 	=> { Enrich => "Total number of mutations of VarClass (restricting to genes in the background)",
						 Contra => "Total number of mutations of VarClass for each sample group (restricting to genes in the background)" } ,
			NGenesInMuts => { Enrich => "Total number of genes associated mutations of VarClass (restricting to genes in the background)",
							  Contra => "Total number of genes associated mutations of VarClass for each sample group (restricting to genes in the background)" },
			ObsStat => { Enrich => "Observed hit statistics (weighted number) of VarClass mutations in GeneSet",
						 Contra => "Observed relative hit proportion for SampGrp 1 of VarClass mutations in GeneSet" },
			ExpStat =>  { Enrich => "Expected hit statistics (weighted number) of mutations of VarClass in GeneSet averaged over simulations",
						  Contra => "Expected relative hit proportion for SampGrp 1 of VarClass mutations in GeneSet" },
			FoldEnrich => "Fold enrichment of ObsStat over ExpStat",
			Pval  => "One-sided p-value of observed statistics",
			Genes => "Name of genes hit by mutations of VarClass (and times of hits if recurrent)",
			GenesGrp1 => "Name of genes hit by mutations of VarClass in sample group 1 (and times of hits if recurrent)",
			GenesGrp2 => "Name of genes hit by mutations of VarClass in sample group 2 (and times of hits if recurrent)");

my %fouts;
foreach my $label (qw|Enrich Contra|) {
	$fouts{$label} = IO::File->new("$outdir/$label.txt", "w");
	$fouts{$label}->print(join("\t", @{$fields{$label}}), "\n");
	open my $fout, ">$outdir/$label.README.txt" or die "Cannot write to $outdir/$label.README.txt"; 
	foreach my $field (@{$fields{$label}}) {
		print $fout $field, "\t", ref $desc{$field} ? $desc{$field}{$label} : $desc{$field}, "\n";
	}
}


for(my $ii = 1; $ii <= @analyses; $ii ++) {
	my $anal = $analyses[$ii-1];
	open my $fin, "$pardir/$anal" or die "Cannot open $pardir/$anal";
	my $line = <$fin>; chomp($line);
	my @dat = split(/\t/, $line);
	my $varcls = shift @dat;
	my $bglabel = shift @dat;
	my @sgrp = @dat;
	my $outfile;
	my @fields;
	if (@sgrp == 1) {
		$outfile = "${varcls}_$sgrp[0]_$bglabel";
	}
	elsif (@sgrp == 2) {
		$outfile = "${varcls}_$sgrp[0]vs$sgrp[1]_$bglabel";
	}
	else {
		die "Incorrect number of sample groups!";
	}
	my $results = collect_res("$wrkdir/dnenrich_$outfile", scalar(@sgrp));

	foreach my $res (@$results) {
		$res->{VarClass} = $varcls;
		$res->{BackGround} = $bglabel;
		if (@sgrp == 1) {
			$res->{SampGroup} = $sgrp[0];
			$fouts{Enrich}->print(join("\t", @{$res}{@{$fields{Enrich}}}), "\n");
		}
		else {
			$res->{SampGrp1} = $sgrp[0];
			$res->{SampGrp2} = $sgrp[1];
			$fouts{Contra}->print(join("\t", @{$res}{@{$fields{Contra}}}), "\n");
		}
	}	
}


sub collect_res {
	my ($resfile, $mode) = @_;
	open my $fin, $resfile or die "Cannot open $resfile";
	# Test for enrichment
	# KEGG_WNT_SIGNALING_PATHWAY	*	0.0514379	3.332	1.39004	149	257	248	.
	# geneset, test type, p-value, observed hit stat, exp hit stat, tot num genes in gs, num mut, num gene in muts, placehold
	my (@res, %genes1, %genes0);
	while(<$fin>) {
		chomp;
		if (/^__GRP(\d+)/) {
			my $grpid = $1;
			unless($grpid == 0 || $grpid == 1) {
				die "Incorrect group ID: $grpid";
			}
			# __GRP1  CHD_BAD_HITS    04-0042 1       1       1       2       745     1
			# Group,  gene set, sampid, mut weight, index_1, index_2, tot num mut in samp, gene ID, gene weight  
			my ($sgrp, $set, $sampid, $mutwt, $index_1, $index_2, $nmut_insamp, $geneid, $genewt) = split(/\t/, $_);
			my $gname = $alias{$geneid} // $geneid;
			if ($grpid == 1) {
				$genes1{$set}{$gname} ++;
			}
			else {
				$genes0{$set}{$gname} ++;
			}
		}
		else {
			my ($set, $test, $pval, $obs_stat, $exp_stat, $ngene_inset, $nmut, $ngene_inmut, $placehold) = split(/\t/, $_);
			next unless $obs_stat > 0;
			push @res => { GeneSet => $set, Test => $test, Pval => $pval, 
						   ObsStat => $obs_stat, ExpStat => $exp_stat, FoldEnrich => sprintf("%.3f", $obs_stat/$exp_stat),
						   NGenesInSet => $ngene_inset, NMuts => $nmut, NGenesInMuts => $ngene_inmut  }; 
		}
	}
	foreach my $res (@res) {
		my (@genes_grp1, @genes_grp0);
		my $minct;
		if ($res->{Test} eq '*') {
			$minct = 1;
		}
		else {
			$minct = $res->{Test};
		}
		my $genes1 = $genes1{$res->{GeneSet}};
		push @genes_grp1, map {  $genes1->{$_} > 1 ? "$_(x$genes1->{$_})" : $_ } 
			sort { $genes1->{$b} <=> $genes1->{$a} || $a cmp $b } 
			grep { $genes1->{$_} >= $minct } keys %{$genes1};
		my $genes0 = $genes0{$res->{GeneSet}};
		if ($mode == 2) {
			push @genes_grp0, map { $genes0->{$_} > 1 ? "$_(x$genes0->{$_})" : $_ } 
				sort { $genes0->{$b} <=> $genes0->{$a} || $a cmp $b } 
				grep { $genes0->{$_} >= $minct } keys %{$genes0};
			$res->{GenesGrp1} = @genes_grp1 > 0 ? join(",", @genes_grp1) : ".";
			$res->{GenesGrp2} = @genes_grp0 > 0 ? join(",", @genes_grp0) : ".";
		}
		else {
			$res->{Genes} = @genes_grp1 > 0 ? join(",", @genes_grp1) : ".";
		}
	}
	return \@res;
}
