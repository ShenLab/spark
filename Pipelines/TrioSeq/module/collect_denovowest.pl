use strict;
use warnings;
use IO::Dir;
use Config::Std;
use Data::Dumper;
use FindBin qw|$Bin|;
use List::Util qw|first min|;
use List::MoreUtils qw|all none|;
use Utils::File qw|count_line|;
use Utils::File::Iter qw|iter_file|;
use Statistics::R;

use lib "$Bin/../../lib";
use Shared qw|slurp_xref|;

unless(@ARGV == 3) {
	print STDERR "$0 PARDIR WRKDIR OUTDIR\n";
	exit 1;
}

my ($pardir, $wrkdir, $outdir) = @ARGV;

my %results;
my @analyses = grep { /^genes_bygrp\.\d+$/ } IO::Dir->new($pardir)->read();

# For each gene, values for the following fields of each sample group will be collected:
# AllObserved, AllExpected, pAllEnrich, MisObserved, MisExpected, pMisEnrich, MisEvents, MisDist, pMisCluster
# Then pMisComb will be calculated by Fisher's method, and pDenovoWEST will be calculated
# Missing values will be coded as "."
for(my $ii = 1; $ii <= @analyses; $ii ++) {
	my $anal = $analyses[$ii-1];
	open my $fin, "$pardir/$anal" or die "Cannot open $pardir/$anal";
	my $line = <$fin>; chomp($line);
	my ($gene, $group) = (split(/\t/, $line))[0, 1];
	slurp_enrich($gene, $group);
	slurp_cluster($gene, $group);
}

# Additional gene level information
my ($xtraginfo, @gfields);
if (-f "$pardir/GeneInfo.conf") {
	read_config "$pardir/GeneInfo.conf" => my %conf;
	my ($xgdat, $xgfds) = slurp_xref($conf{GeneTab}{GXref}, $conf{GeneTab}{GXref_Fields});
	$xtraginfo = $xgdat;
	@gfields = @$xgfds;
}

#
open my $fout, ">$outdir/README.txt" or die "Cannot write to README.txt";
print $fout <<EOF;
AllObserved	Observed weighted sum statistics for all DNVs
AllExpected	Expected weighted sum statistics for all DNVs from simulation using baseline mutation rates
pAllEnrich 	One-sided p-value for enrichment of all DNVs 
MisObserved	Observed weighted sum statistics for missense DNVs
MisExpected	Expected weighted sum statistics for missense DNVs from simulation using baseline mutation rates
pMisEnrich	One-sided p-value for enrichment of missense DNVs
MisEvents 	Total number of missense mutational events used by DenovoNear
MisDist 	The observed geometric mean of coding distances for all pairs of missense DNVs 
pMisCluster	One-side p-value for clustering of missense DNVs
pMisComb 	Combined p-value for enrichment and clustering of missense DNVs 
pDenovoWEST	The final DenovoWEST p-value as the minimum of pAllEnrich and pMisComb
EOF

foreach my $group (keys %results) {
	my %res = %{$results{$group}};

	my (@stdfields, $sortkey);
	if (none { defined $_->{pMisCluster} } values %res) {
		if (none { defined $_->{pMisEnrich} } values %res) {
			@stdfields = qw|AllObserved AllExpected pAllEnrich|;
			$sortkey = "pAllEnrich";
		}
		else {
			@stdfields = qw|AllObserved AllExpected pAllEnrich MisObserved MisExpected pMisEnrich|;	
			foreach my $gid (keys %res) {
				$res{$gid}{pDenovoWEST} = min($res{$gid}{pAllEnrich}, $res{$gid}{pMisEnrich});
			}
			$sortkey = "pDenovoWEST";
		}
	}
	else {
		@stdfields = qw|AllObserved AllExpected pAllEnrich MisObserved MisExpected pMisEnrich MisEvents MisDist pMisCluster pMisComb pDenovoWEST|;
		# Calculate pMisComb and pDenovoWEST...
		my $R = Statistics::R->new();
		foreach my $gid (keys %res) {
			if (defined $res{$gid}{pMisEnrich} && defined $res{$gid}{pMisCluster}) {
				if ($res{$gid}{pMisCluster} eq '.') {
					$res{$gid}{pMisComb} = $res{$gid}{pMisEnrich};
				}
				else {
					$R->set('x', -2*(log($res{$gid}{pMisEnrich})+log($res{$gid}{pMisCluster})));
					$R->run('pval<-pchisq(x, 4, lower.tail=F)');
					$res{$gid}{pMisComb} = $R->get('pval');
				}
				if (defined $res{$gid}{pAllEnrich}) {
					$res{$gid}{pDenovoWEST} = min($res{$gid}{pMisComb}, $res{$gid}{pAllEnrich});
				}
				else {
					warn "Cannot find pAllEnrich for $gid!";
					$res{$gid}{pDenovoWEST} = $res{$gid}{pMisComb};
				}
			}
			elsif (defined $res{$gid}{pMisEnrich}) {
				$res{$gid}{pMisComb} = $res{$gid}{pMisEnrich};
				$res{$gid}{pDenovoWEST} = min($res{$gid}{pMisComb}, $res{$gid}{pAllEnrich});
			}
			else {
				$res{$gid}{pMisComb} = ".";
				$res{$gid}{pDenovoWEST} = $res{$gid}{pAllEnrich};
			}
		}
		$sortkey = "pDenovoWEST";
	}

	open my $fout, ">$outdir/DenovoWEST_${group}.txt" or die "Cannot write to $outdir/DenovoWEST_${group}.txt";
	print $fout join("\t", "GeneID", @gfields, @stdfields), "\n";
	# Sort the results based on pDenovoWEST?
	# Map . to 1
	my @genes = sort { $res{$a}{$sortkey} <=> $res{$b}{$sortkey}  } 
				grep { defined $res{$_}{$sortkey} && $res{$_}{$sortkey} ne '.' } keys %res;
	push @genes, grep { !defined $res{$_}{$sortkey} || $res{$_}{$sortkey} eq '.' } sort keys %res;

	foreach my $gid (@genes) {
		my @ginfo;
		if (defined $xtraginfo) {
			@ginfo = map { $xtraginfo->{$gid}{$_} // "." } @gfields;  
		}
		my %gres;
		foreach my $fd (@stdfields) {
			if ($fd =~ /^p/ && defined $res{$gid}{$fd} && $res{$gid}{$fd} ne ".") {
				$gres{$fd} = $res{$gid}{$fd} < 0.01 ? sprintf("%.2E", $res{$gid}{$fd}) : sprintf("%.2f", $res{$gid}{$fd});
			}
			else {
				$gres{$fd} = $res{$gid}{$fd};
			}
		}
		print $fout join("\t", $gid, @ginfo, map { $gres{$_} // "." } @stdfields), "\n";
	}
}



#
# Slurp results files
#
# Format of {all/mis}enrich.txt
# GeneID  Expected        Observed        P-value Info
# ENSG00000000457 0.3692423245755087      0.982   0.06771492872249747     629722811|probability of observing >= 17 mutations is too small
sub slurp_enrich {
	my ($gene, $group) = @_;
	foreach my $typ (qw|all mis|) {
		#print "$wrkdir/$group/${gene}_${typ}enrich.txt\n";
		next unless -f "$wrkdir/$group/${gene}_${typ}enrich.txt";
		my $it = iter_file("$wrkdir/$group/${gene}_${typ}enrich.txt", 
			{ fsep => qr/\t/, alias => {'P-value' => 'pEnrich'} });
		while(my $dat = $it->()) {
			unless($gene eq $dat->{GeneID}) {
				die "Incorret gene ID: $gene ne $dat->{GeneID}";
			}
			foreach my $fd (qw|Expected Observed pEnrich|) {
				my $Typ = ucfirst($typ);
				my $key = $fd; $key =~ s/(?=[A-Z])/$Typ/;
				$results{$group}{$dat->{GeneID}}{$key} = $dat->{$fd};
			}
		}
	}
}

# Format of dncluster.txt: 
# gene_id	mutation_category	events_n	dist	probability
# ENSG00000104888	missense	2	501.0	0.5080324919675079
# ENSG00000104888	nonsense	0	nan	nan
# Also collect missense variants from input
# missense = ["missense_variant", "stop_lost"]

sub slurp_cluster {
	my ($gene, $group) = @_;
	if (-f "$wrkdir/$group/${gene}_dncluster.txt") {
		my $nline = count_line("$wrkdir/$group/${gene}_dncluster.txt");
		if ($nline > 1) {
			my $it = iter_file("$wrkdir/$group/${gene}_dncluster.txt", 
				{ fsep => qr/\t/, alias => {'gene_id' => 'GeneID', 'events_n' => 'Events', 'dist' => 'Dist', 'probability' => 'pCluster'} });
			while(my $dat = $it->()) {
				next unless $dat->{mutation_category} eq 'missense';
				unless($gene eq $dat->{GeneID}) {
					die "Incorrect gene ID in ${gene}_dncluster: $gene ne $dat->{GeneID}";
				}
				if ($dat->{Dist} =~ /,/) {
					my @dist = split(',', $dat->{Dist});
					if (all { $_ eq 'nan' } @dist) {
						$dat->{Dist} = 'nan';
					}
					else {
						$dat->{Dist} = first { $_ ne 'nan' } @dist;
					}
				}
				foreach my $fd (qw|Events Dist pCluster|) {
					my $key = $fd; $key =~ s/(?=[A-Z])/Mis/;
					$results{$group}{$dat->{GeneID}}{$key} = $dat->{$fd} eq 'nan' ? '.' : $dat->{$fd};
				}
			}
		}
		else {
			unless (-f "$wrkdir/$group/${gene}_dnclstin.txt") {
				die "Cannot find $wrkdir/$group/${gene}_dnclstin.txt";
			}
			my $it = iter_file("$wrkdir/$group/${gene}_dnclstin.txt", 
				{ fsep => qr/\t/, alias => { 'gene_name' => 'GeneID' } });
			my $nmis = 0;
			while(my $dat = $it->()) {
				unless($gene eq $dat->{GeneID}) {
					die "Incorrect gene ID in ${gene}_dnclstin: $gene ne $dat->{GeneID}";
				}
				if ($dat->{consequence} =~ /missense|stop_lost|coding_sequence_variant/ && 
					$dat->{snp_or_indel} =~ /SNV/) {
					$nmis ++;
				}
			}
			$results{$group}{$gene}{MisEvents} = $nmis;
			$results{$group}{$gene}{MisDist} = ".";
			$results{$group}{$gene}{pCluster} = ".";
		}
	}
}
