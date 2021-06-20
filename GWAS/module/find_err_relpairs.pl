#!/usr/bin/env perl

use strict;
use warnings;
use Carp;
use Graph;
use Data::Dumper;
use Config::Std;
use List::MoreUtils qw|all|;
use Genet::Ped;
use Utils::File::Iter qw|iter_file|;

unless (@ARGV >= 2) {
	print STDERR "This script is aimed for a quick automated soultion to cleanup sample errs\n";
	print STDERR "in family-based GWAS sample for further marker-based QC.\n";
	print STDERR "$0 INPUT_PREFIX PARAMS [OUTPUT_PREFIX]\n";
	exit 1;
}


my ($input, $param, $output) = @ARGV;
$output = $input unless defined $output;

unless(all { -f "$input.$_" } qw|fam bed bim imiss|) {
	print STDERR "Not all required input files (.fam, .bed, .bim and .imiss) can be found!\n";
	exit 1;
}
unless(-f "$input.kin" || -f "$input.kin0") {
	print STDERR "No kinship estimation (.kin and .kin0) can be found!\n";
	exit 0;
}

# Missing genotype rates
my %miss;
{
	my $it = iter_file("$input.imiss", {chomp => 3});
	while(my $dat = $it->()) {
		$miss{$dat->{FID},$dat->{IID}} = $dat->{F_MISS};
	}
}

read_config $param => my %conf or die "Cannot read config file $param";
my %param = %{$conf{""}};

my (@remove, @dups, @rels);

# For family-based sample, 
# 1. Identify all incorrect 1-d rel pairs, and remove appropriate one
# 2. Make a plot of kinship vs IBS0 colored by ped relationships
# (cryptic relatedness between parents can be identified in 2, but no sample will be removed)
my $gped = Graph::Undirected->new;
open my $ftab, ">$output.relpairs.csv" or die "Cannot write to relpairs table";
print $ftab join(",", qw|FID1 ID1 FID2 ID2 Kinship IBS0 IBD1 IBD2 PedType InfType|), "\n";
my $plotflag;

my %pairs;
# First look for duplicated pairs
{
	my $pattern = qr/$param{IGNORE}$/;
	open my $fin, "$input.fam" or die "Cannot open $input.fam";
	while(<$fin>) {
		my ($fid, $iid) = (split)[0,1];
		if ($iid =~ /$pattern$/) {
			(my $oldid = $iid) =~ s/$pattern$//;
			my @iids = sort($oldid, $iid);
			$pairs{$fid,$oldid,$iid} = 'Dup/MZ';
		}
	}
}

if (-f "$input.kin") {
	# Read pedigree file to find 1d pairs
	print STDERR "Verify relationships within families\n";	
	my $ped = Genet::Ped->new("$input.fam", { verbose => 0 });
	my $popairs = $ped->get_po_pairs();
	while (my ($fid, $pairs) = each %{$popairs}) {
		foreach my $pair (@$pairs) {
			my @iids = sort($pair->[0],$pair->[1]);
			$pairs{$fid,$iids[0],$iids[1]} = 'PO';
		}
	}
	my $sbpairs = $ped->get_sib_pairs();
	while (my ($fid, $pairs) = each %{$sbpairs}) {
		foreach my $pair (@$pairs) {
			my @iids = sort($pair->[0],$pair->[1]);
			next if defined $pairs{$fid,$iids[0],$iids[1]};
			$pairs{$fid,$iids[0],$iids[1]} = 'Sib';
		}
	}

	my $gerr = Graph::Undirected->new;

	my %known;
	# Read king output to compare inferred relationships.
	# A plot will be made to show the pedigree and inferred relationships
	my $it = iter_file("$input.kin");
	while(my $dat = $it->()) {
		my @iids = sort ($dat->{ID1},$dat->{ID2});
		$dat->{FIID1} = $dat->{FID}."\t".$iids[0];
		$dat->{FIID2} = $dat->{FID}."\t".$iids[1];

		my $inftype;		
		if ($dat->{Kinship} >= $param{MIN1D} && $dat->{Kinship} <= $param{MAX1D}) {
			if ($dat->{IBS0} <= $param{POIBS0}) {
				$inftype = 'PO';
			}
			else {
				$inftype = 'Sib';
			}
		}
		else {
			if ($dat->{Kinship} >= $param{MINDUP}) {
				$inftype = 'Dup/MZ';
			}
			elsif ($dat->{Kinship} >= $param{CRPREL}) {
				$inftype = "CypRel";
			}
			else {
				$inftype = "Other";
			}
		}
		if (defined $pairs{$dat->{FID},$iids[0],$iids[1]}) {
			if ($pairs{$dat->{FID},$iids[0],$iids[1]} eq $inftype) {
				$gped->add_edge($dat->{FIID1}, $dat->{FIID2});
			}
			else {
				# We can keep special this case as correct?
				if ($inftype eq 'Dup/MZ') {
					push @dups, [@{$dat}{qw|FID ID1 FID ID2|}];
				}
				$gerr->add_edge($dat->{FIID1}, $dat->{FIID2});
			}
		}
		print $ftab join(",", @{$dat}{qw|FID ID1 FID ID2 Kinship IBS0 IBD1Seg IBD2Seg|}, 
			$pairs{$dat->{FID},$iids[0],$iids[1]} // "Other", $inftype), "\n";
		$plotflag = 1;
		$known{$dat->{FID},$iids[0],$iids[1]} = 1;
	}

	# Check that all pedigree 1d-pairs have kinship estimates
	unless(all { $known{$_} } keys %pairs) {
		croak "Not all pedigree 1d pairs have estimated kinship";
	}

	# Now need to decide which samples to exclude to eliminate incorrect relations.
	push @remove, find_removal($gerr, $gped, 1);
}

# For between family samples
# 1. Look for unexpected duplicates or cryptic rel pairs, and remove appropriate samples
# Note: we should set threshold for cryptic rel higher to avoid removing too many sample
if (-f "$input.kin0") {
	print STDERR "Identify unexpected duplicates or 1-d related pairs between families\n";
	my $gerr = Graph::Undirected->new;
	my $it = iter_file("$input.kin0");
	while(my $dat = $it->()) {
		my $inftype;
		if ($dat->{Kinship} >= $param{CRPREL}) {
			
			$dat->{FIID1} = $dat->{FID1}."\t".$dat->{ID1};
			$dat->{FIID2} = $dat->{FID2}."\t".$dat->{ID2};
			$gerr->add_edge($dat->{FIID1}, $dat->{FIID2});

			if ($dat->{Kinship} >= $param{MINDUP}) {
				$inftype = 'Dup/MZ';
				push @dups, [@{$dat}{qw|FID1 ID1 FID2 ID2|}];
			}
			elsif ($dat->{Kinship} >= $param{MIN1D} && $dat->{Kinship} <= $param{MAX1D}) {
				if ($dat->{IBS0} <= $param{POIBS0}) {
					$inftype = 'PO';
				}
				else {
					$inftype = 'Sib';
				}
				push @rels, [@{$dat}{qw|FID1 ID1 FID2 ID2|}];
			}
			else {
				$inftype = 'CypRel';
				push @rels, [@{$dat}{qw|FID1 ID1 FID2 ID2|}];
			}
		}

		# CypRel will not be used for plotting to avoid cluttering
		if (defined $inftype && $inftype ne 'CypRel') {
			print $ftab join(",", @{$dat}{qw|FID1 ID1 FID2 ID2 Kinship IBS0 IBD1Seg IBD2Seg|}, 
				"Unrelated", $inftype), "\n";
			$plotflag = 1;
		}
	}
	push @remove, find_removal($gerr, $gped, 1);
}

# Make plot for inspecting error pairs within family
# and for incorrect 1-d pair between family
if ($plotflag) {
	close $ftab;
	my $script = <<'EOF';
library(ggplot2)
library(ggrepel)
REL<-read.csv("_INFILE_.csv", header=TRUE, as.is=TRUE)
rownames(REL)<-with(REL, ifelse(FID1==FID2, paste(FID1,paste(ID1,ID2,sep=","),sep=":"), 
					paste(paste(FID1,ID1,sep=":"), paste(FID2,ID2,sep=":"), sep=",")))
ERR<-subset(REL, PedType!=InfType)
g<-ggplot(REL) + geom_point(aes(x=IBS0, y=Kinship, color=InfType, shape=PedType)) +
	labs(y="Kinship coefficient", color="Inferred\nRelationship", shape="Pedigree\nRelatinship")
if (nrow(ERR) > 0 & nrow(ERR)< 20) {
	g<-g+geom_text_repel(data=ERR, aes(x=IBS0, y=Kinship), label=rownames(ERR), size=2)
}
ggsave("_INFILE_.png", g, unit="in", width=6, height=5)

f<-ggplot(REL) + geom_point(aes(x=IBD1, y=IBD2, color=InfType, shape=PedType)) +
	labs(x="Prob(IBD=1)", y="Prob(IBD=2)", color="Inferred\nRelationship", shape="Pedigree\nRelatinship")
if (nrow(ERR) > 0 & nrow(ERR)< 20) {
	f<-f+geom_text_repel(data=ERR, aes(x=IBD1, y=IBD2), label=rownames(ERR), size=2)
}
ggsave("_INFILE_.2.png", f, unit="in", width=6, height=5)


EOF
	$script =~ s/_INFILE_/$output.relpairs/g;
	open my $fs, ">$output.relpairs.R" or die "Cannot write R script";
	print $fs $script;
	close $fs;
	system("Rscript $output.relpairs.R");
}

if (@remove) {
	open my $frm, ">$output.relpairs.remove" or die "Cannot write removal list";
	foreach my $fiid (@remove) {
		print $frm $fiid, "\n";
	}
}

if (@dups) {
	open my $fdup, ">$output.relpairs.dups.txt" or die "Cannot write to duplicated pairs";
	foreach my $dat (@dups) {
		print $fdup join("\t", @$dat), "\n";
	}
}

if (@rels) {
	open my $frel, ">$output.relpairs.cyprel.txt" or die "Cannot write to cryptic related pairs";
	foreach my $dat (@rels) {
		print $frel join("\t", @$dat), "\n";
	}
}


# Given graphs of correct and error pedigree relationships
# Determine the minimal list of samples to be removed
# Basic algorithm is to start with nodes that only in error graph
# Nodes in the error graph will be ranked by degree in the error graph and
# the degree ped graph (reverse). When there is a tie, will be ranked by missing
# genotype rates. Then nodes will then be sequentially removed until no edge
# left in the error graph.
sub find_removal {
	my ($gerr, $gped, $verbose) = @_;
	my @remove;
	# Find all edges
	my $count = 1;
	while($gerr->edges()) {
		my @errnd;
		foreach my $node ($gerr->vertices()) {
			my ($fid, $iid) = split(/\t/, $node);
			my $dgerr = $gerr->degree($node);
			my $dgped;
			if (defined $gped) {
				$dgped = $gped->degree($node) // 0;
			}
			else {
				$dgped = 0;
			}
			push @errnd, [$node, $dgerr, $dgped, $miss{$fid,$iid}];
		}
		@errnd = sort {  $b->[1] <=> $a->[1] or $a->[2] <=> $b->[2] or 
			$b->[3] <=> $a->[3] } @errnd;
		if ($verbose) {
			print STDERR "Round $count : --------------------\n";
			print STDERR join("\n", map { join("\t", @$_) } @errnd), "\n";
		}
		$gerr->delete_vertex($errnd[0][0]);
		push @remove, $errnd[0][0];
		$count ++;
	}
	if ($verbose) {
		my @remain = $gerr->vertices();
		print STDERR "Final : --------------------\n";
		print STDERR join("\n", @remain), "\n";
	}
	return @remove;
}


