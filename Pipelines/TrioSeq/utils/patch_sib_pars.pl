#!/usr/bin/env perl

use strict;
use warnings;
use Perl6::Slurp;
use List::MoreUtils qw|all none uniq|;
use Getopt::Euclid;
use Graph::Undirected;
use Data::Dumper;
use String::Random;
use Utils::List qw|all_pairs|;

unless(defined $ARGV{'--sibs'} || defined $ARGV{'--half-sibs'}) {
	die "Must provide either sib pair list or half-sib pair list";
}

# Slurp FamID/Sex/Pheno for each individual in the PED file
my %famid = map { (split)[1,0] } slurp $ARGV{'--ped'};
my %sex = map { (split)[1,4] } slurp $ARGV{'--ped'};
my %phe = map { (split)[1,5] } slurp $ARGV{'--ped'};

my $ignore;
if (defined $ARGV{'--ignore'}) {
	$ignore = qr/$ARGV{'--ignore'}/;
}

my (%famsamps, %dad, %mom, %rep, %knownids, %newids);
# %famsamps: available individuals in each family
# %dad/%mom: father/mother of each individual, must be non-0, but can be those not exist in PED
# %rep: original sample ID for tech dups
# %knownids: any sample IDs that already appear in IID column or as parents 
# %newids: parents IDs that are added later to fix the PED file, IID=>[FamID, Sex], those parents are assumed to be founders
{
	open my $fin, $ARGV{'--ped'} or die "Cannot open PED file";
	while(<$fin>) {
		my ($fid, $iid, $dad, $mom, $sex) = (split)[0,1,2,3,4];
		$knownids{$iid} = [$fid, $sex];
		push @{$famsamps{$fid}}, $iid;
		if (defined $ignore && $iid =~ /$ignore/) {
			(my $origid = $iid) =~ s/$ignore//;
			if (defined $famid{$origid}) {
				die "Sample $origid and duplicate $iid do not have the same FamID and sex"
					unless $famid{$origid} eq $fid && $sex{$origid} eq $sex;
				$rep{$iid} = $origid;
				next;
			}
		}
		if ($dad ne '0') {
			$dad{$iid} = $dad;
			$knownids{$dad} = [$fid, 1];
		}
		if ($mom ne '0') {
			$mom{$iid} = $mom;
			$knownids{$mom} = [$fid, 2];
		}
	}
}

# Create a list of parents that are not present in the ped file
my %mispars;
# %mispar: known parents in each family that does appear in IID column of PED file ("missing parents")
{
	foreach my $iid (keys %knownids) {
		unless(defined $famid{$iid}) {
			my $fid = $knownids{$iid}[0];
			push @{$mispars{$fid}} => $iid;
		}
	}
}

my $idgen = String::Random->new;


# Individuals who are part of sib pairs
# IID => { SibID => Sentinel ID, Father => ID, Mother => ID}
my %sibs; 

# Create graphs for full sib pairs
# Full sibs who share the same parents should form a complete subgraph
my $gsib = Graph::Undirected->new();
{
	# First include all sib pairs who already shared both parents in PED file
	foreach my $fid (keys %famsamps) {
		foreach my $pair (all_pairs(@{$famsamps{$fid}})) {
			next if $gsib->has_edge($pair->[0], $pair->[1]);
			if (defined $mom{$pair->[0]} && defined $mom{$pair->[1]} &&
				defined $dad{$pair->[0]} && defined $dad{$pair->[1]} &&
				$mom{$pair->[0]} eq $mom{$pair->[1]} && $dad{$pair->[0]} eq $dad{$pair->[1]}) {
				$gsib->add_edge($pair->[0], $pair->[1]);
			}
		}
	}
	# Then adding additional sib pairs from the specified list
	# After that we should assume all sib pairs should have been included in the graph.
	if ($ARGV{'--sibs'}) {
		open my $fin, $ARGV{'--sibs'} or die "Cannot open $ARGV{'--sibs'}";
		while(<$fin>) {
			my @pair = split;
			unless(@pair == 2) {
				die "Incorrect number of columns in sibs pair: $_";
			}
			unless(all { defined $famid{$_} } @pair) {
				print STDERR join("\t", @pair), "\n";
				die "Cannot find both sib pair in the PED file";
			}
			else {
				unless($famid{$pair[0]} eq $famid{$pair[1]}) {
					print STDERR "$pair[0]:$famid{$pair[0]} and $pair[1]:$famid{$pair[1]}\n";
					die "Not all sib pairs share the same family ID";
				}
			}
			unless($gsib->has_edge(@pair)) {
				$gsib->add_edge(@pair);
			}
		}
	}
	# Find all connected components. If all sib pairs have been added to the graph
	# each component should be a complete subgraph (but we will not check for this). 
	# The first sample from this (sorted) component can be the sentinel sample and everyone in 
	# the same component should share the same pair of parents.
	foreach my $cc ($gsib->connected_components()) {
		my @sibsamp = sort @$cc;
		my @father = uniq sort grep { $_ } map { $dad{$_} } @sibsamp; 
		my @mother = uniq sort grep { $_ } map { $mom{$_} } @sibsamp; 
		if (@father > 1 || @mother > 1) {
			print STDERR "Sibs: ", join(", ", @sibsamp), "\n";
			print STDERR "Father: ", join(", ", @father), "\n";
			print STDERR "Mother: ", join(", ", @mother), "\n";
			die "Multiple fathers or mothers found for a full sib pair";
		}
		# For some sample, his/her father or mother may not be defined at this point.
		foreach my $iid (@sibsamp) {
			$sibs{$iid} = { SibID => $sibsamp[0], Father => $father[0], Mother => $mother[0] };
		}
	}
}

# HalfSibs graph, half sib pairs will be connected
my $ghalfsib = Graph::Undirected->new();
{
	# First include all sib pairs who share one and only one parent in PED file
	# Those share one parent but the other is missing are not included
	foreach my $fid (keys %famsamps) {
		foreach my $pair (all_pairs(@{$famsamps{$fid}})) {
			next if $gsib->has_edge($pair->[0], $pair->[1]);
			next if $ghalfsib->has_edge($pair->[0], $pair->[1]);
			if (defined $mom{$pair->[0]} && defined $mom{$pair->[1]} &&
				defined $dad{$pair->[0]} && defined $dad{$pair->[1]} &&
				($mom{$pair->[0]} eq $mom{$pair->[1]} && $dad{$pair->[0]} ne $dad{$pair->[1]} || 
				 $mom{$pair->[0]} ne $mom{$pair->[1]} && $dad{$pair->[0]} eq $dad{$pair->[1]}) ){
				my $iid1 = defined $sibs{$pair->[0]} ? $sibs{$pair->[0]}{SibID} : $pair->[0];
				my $iid2 = defined $sibs{$pair->[1]} ? $sibs{$pair->[1]}{SibID} : $pair->[1];
				# If one of the half sib also have full sibs, then the sample should be represented
				# by the sentinel sample of the full sibs.
				# We assume that all full sib pairs has been found above.
				if ($iid1 eq $iid2) {
					die "Half-sib pair $pair->[0],$pair->[1] belongs to a known full-sib par";
				}
				$ghalfsib->add_edge($pair->[0], $pair->[1]);
			}
		}
	}
	# We will also assume that all half sibs should be included after further slurping in pairs
	# provided in the list below.
	if (defined $ARGV{'--half-sibs'}) {
		open my $fin, $ARGV{'--half-sibs'} or die "Cannot open $ARGV{'--half-sibs'}";
		while(<$fin>) {
			my @pair = split;
			unless(@pair == 2) {
				die "Incorrect number of columns in half-sibs pair list: $_";
			}
			unless(all { defined $famid{$_} } @pair) {
				die "Cannot find both half-sib pair in the PED file";
			}
			else {
				unless($famid{$pair[0]} eq $famid{$pair[1]}) {
					die "Not all half-sib pairs share the same family ID";
				}
			}	
			my @uv = map { defined $sibs{$_} ? $sibs{$_}{SibID} : $_ } @pair;
			if ($uv[0] eq $uv[1]) {
				die "Half-sib pair $pair[0],$pair[1] belongs to a known full-sib par";
			}
			unless($ghalfsib->has_edge(@uv)) {
				$ghalfsib->add_edge(@uv);
			}
		}
	}
}

# Check that samples share at least one parents should be included in either sib graph
# or half-sib graph at this point
foreach my $fid (keys %famsamps) {
	foreach my $pair (all_pairs(@{$famsamps{$fid}})) {
		if (defined $dad{$pair->[0]} && defined $dad{$pair->[1]} && $dad{$pair->[0]} eq $dad{$pair->[1]} ||
			defined $mom{$pair->[0]} && defined $mom{$pair->[1]} && $mom{$pair->[0]} eq $mom{$pair->[1]}) {
			my @uv = map { defined $sibs{$_} ? $sibs{$_}{SibID} : $_ } @$pair;
			unless($gsib->has_edge($pair->[0],$pair->[1]) || $ghalfsib->has_edge($uv[0],$uv[1])) {
				die "Sample $pair->[0],$pair->[1] share at least one parent but cannot be found in sib or half-sib graph";
			}
		}
	}
}

# Individuals who are part of half-sib pairs
# IID => { Father => ID, Mother => ID}
my %halfsibs;

# We first process half-sib graphs to fill in missing parents
if ($ARGV{'--half-sibs'}) {
	foreach my $cc ($ghalfsib->connected_components()) {
		my @sibsamp = sort @$cc;
		foreach my $iid (@sibsamp) {
			if (defined $sibs{$iid}) {
				$halfsibs{$iid} = { Father => $sibs{$iid}{Father}, Mother => $sibs{$iid}{Mother} };
			}
			else {
				$halfsibs{$iid} = { Father => $dad{$iid}, Mother => $mom{$iid} };
			}
		}
		# We will fill in parental data by the following Algorithm:
		# As long as there are samples with missing parent(s), Loop:
		# we sort those samples based on the completness of parents in graph neighbors
		# (from largest to smallest amount of available parents) then fill in parents in 
		# the first sample with smallest missingness .... until all parents are filled in.
		# Validate parents info in this connected component

		# Find samples with incomplete parents info
		my @sampmisspar = grep { !defined $halfsibs{$_}{Father} || !defined $halfsibs{$_}{Mother} } @sibsamp;
		while(@sampmisspar > 0) {
			# Calculate missingness and neighbor information
			my %metrics;
			foreach my $iid (@sampmisspar) {
				my @par_miss = grep { !defined $halfsibs{$iid}{$_} } qw|Father Mother|;
				die "Both parents of $iid are available?" unless @par_miss > 0;
				my $ct = 0;
				foreach my $nn ($ghalfsib->neighbors($iid)) {
					die "Cannot find halfsib info for $nn" unless defined $halfsibs{$nn};
					foreach my $par_miss (@par_miss) {
						if (defined $halfsibs{$nn}{$par_miss}) {
							$ct ++;
						}
					}
				}
				$metrics{$iid}{Miss} = scalar(@par_miss);
				$metrics{$iid}{Avail} = $ct/scalar(@par_miss);
			}
			@sampmisspar = 
				sort { $metrics{$b}{Avail} <=> $metrics{$a}{Avail} ||
					   $metrics{$a}{Miss}  <=> $metrics{$b}{Miss} } @sampmisspar;
			# Start with the sample ranked at the first
			my $samp = shift @sampmisspar;
			# Find neighbors in connected components also keep list of all non-neighbors
			my @neighbors = $ghalfsib->neighbors($samp);
			# Fill in missing parents based on the constriaints of neigbors
			if (@neighbors == 1) {
				_fill_hspar_nn(\%halfsibs, $samp, $neighbors[0]);
			}
			else {	
				warn "The sample $samp has multiple neighbors in half-sib graph\n";
				_fill_hspar_nn(\%halfsibs, $samp, $neighbors[0]);
				#die "Currently does not support parental inference in such cases";
			}
		}
	}
}

# After fill in half sib parents, revist full sibs to fill in remaining missing parents
foreach my $cc ($gsib->connected_components()) {
	my @sibsamp = sort @$cc;
	unless(all { $sibs{$_}{SibID} eq $sibsamp[0] } @sibsamp) {
		die "Inconsistent SibID for sibs: ".join(", ", @sibsamp);
	}
	my @Par = qw|Father Mother|;
	for my $ii (0..1) {
		if (none { defined $sibs{$_}{$Par[$ii]} } @sibsamp) {
			if (defined $halfsibs{$sibsamp[0]}) {
				die "Sample $sibsamp[0] is part of half-sib but $Par[$ii] is unknown"
					unless defined $halfsibs{$sibsamp[0]}{$Par[$ii]};
				foreach (@sibsamp) {
					$sibs{$_}{$Par[$ii]} = $halfsibs{$sibsamp[0]}{$Par[$ii]};
				}
			}
			else {
				my $randid = _gen_randid($famid{$sibsamp[0]}, $Par[$ii]);
				foreach (@sibsamp) {
					$sibs{$_}{$Par[$ii]} = $randid;
				}
			}
		}
	}
}

foreach my $iid (keys %sibs) {
	if (!defined $mom{$iid}) {
		$mom{$iid} = $sibs{$iid}{Mother} // do { die "Cannot find mother for sib sample $iid" };
	}
	if (!defined $dad{$iid}) {
		$dad{$iid} = $sibs{$iid}{Father} // do { die "Cannot find father for sib sample $iid" };
	}
}
foreach my $iid (keys %halfsibs) {
	if (!defined $mom{$iid}) {
		$mom{$iid} = $halfsibs{$iid}{Mother} // do { die "Cannot find mother for sib sample $iid" };
	}
	if (!defined $dad{$iid}) {
		$dad{$iid} = $halfsibs{$iid}{Father} // do { die "Cannot find father for sib sample $iid" };
	}
}

# Final validation
# Enumerate all edges in sib and half-sib graphs and test their parents
foreach my $pair ($gsib->edges) {
	unless (defined $dad{$pair->[0]} && defined $mom{$pair->[0]} &&
			defined $dad{$pair->[1]} && defined $mom{$pair->[1]} &&
			$dad{$pair->[0]} eq $dad{$pair->[1]} && $mom{$pair->[0]} eq $mom{$pair->[1]}) {
		die "The sib pair $pair->[0],$pair->[1] does not share the same parents after patch";
	}
}
foreach my $pair ($ghalfsib->edges) {
	unless (defined $dad{$pair->[0]} && defined $mom{$pair->[0]} &&
			defined $dad{$pair->[1]} && defined $mom{$pair->[1]} &&
			($dad{$pair->[0]} eq $dad{$pair->[1]} && $mom{$pair->[0]} ne $mom{$pair->[1]} ||
			 $dad{$pair->[0]} ne $dad{$pair->[1]} && $mom{$pair->[0]} eq $mom{$pair->[1]}) ) {
		die "The half-sib pair $pair->[0],$pair->[1] does not share one and only one parent after patch";
	}
}


# Now create final output
open my $fout, ">$ARGV{'--output'}" or die "Cannot write to $ARGV{'--output'}";
foreach my $fid (sort keys %famsamps) {
	# Adding newly created parents and "missing parents" for this family
	unless($ARGV{'--no-strict'}) {
		foreach my $iid (grep { $newids{$_}[0] eq $fid } sort keys %newids) {
			print $fout join("\t", $fid, $iid, 0, 0, $newids{$iid}[1], 0), "\n";
		}
		if (defined $mispars{$fid}) {
			foreach my $iid (@{$mispars{$fid}}) {
				print $fout join("\t", $fid, $iid, 0, 0, $knownids{$iid}[1], 0), "\n";
			}
		}
	}
	foreach my $iid (@{$famsamps{$fid}}) {
		# For tech dup, write out samp info of original sample
		if (defined $rep{$iid}) {
			my $iid2 = $rep{$iid};
			print $fout join("\t", map { $_ // "0" } ($fid, $iid, $dad{$iid2}, $mom{$iid2}, $sex{$iid2}, $phe{$iid2})), "\n";
		}
		else {
			print $fout join("\t", map { $_ // "0" } ($fid, $iid, $dad{$iid}, $mom{$iid}, $sex{$iid}, $phe{$iid})), "\n";
		}
	}
}


# Fill in half-sib's parents data given one neighbor
sub _fill_hspar_nn {
	my ($hsdat, $sampid, $nnid) = @_;
	my @Par = qw|Father Mother|;
	for my $ii (0..1) {
		if (!defined $hsdat->{$sampid}{$Par[$ii]}) {
			# Check if the other field posed a constraint
			if ( defined $hsdat->{$sampid}{$Par[1-$ii]} && defined $hsdat->{$nnid}{$Par[1-$ii]} &&
				 $hsdat->{$sampid}{$Par[1-$ii]} ne $hsdat->{$nnid}{$Par[1-$ii]} &&
				 defined $hsdat->{$nnid}{$Par[$ii]} ) {
				$hsdat->{$sampid}{$Par[$ii]} = $hsdat->{$nnid}{$Par[$ii]};
			}
			else {
				$hsdat->{$sampid}{$Par[$ii]} = _gen_randid($famid{$sampid}, $Par[$ii]);
			}
		}	
	}
	return $hsdat;
}

# When generating new random IDs, make sure to avoid name crash.
sub _gen_randid {
	my ($FamID, $Par) = @_;
	my $randid = $idgen->randregex($ARGV{'--pattern'});
	while(defined $knownids{$randid} || defined $newids{$randid}) {
		$randid = $idgen->randregex($ARGV{'--pattern'});
	}
	$newids{$randid} = [$FamID, $Par eq 'Father' ? 1 : 2];
	return $randid;
}



__END__

=head1 NAME

patch_sib_pars -- Patch parents for siblings of nuclear families in a PED file

=head1 USAGE

patch_sib_pars [options] -in PEDIN -sib LIST1 -out PEDOUT

=head1 DESCRIPTION

For nuclear families, there are several courses of incomplete pedigree. One is missing parents, so the 
relationship between offspring become ambiguous or unknown based on PED file. For example, if three lines
of the PED file represent a single mother two children family but father's column for the two children
are missing, we cannot determine if two children are full or half sib pairs. Because sibships can be 
confidently verified by genotypes, it provides a solution to fix those incomplete pedigrees.

The script takes input a PED file of nuclear families, and a list of full sib pairs and optionally a second
list of half-sib pairs. It verify and fix parents in the PED file, such that full sibs share both parents
and half-sibs share one and only one parents. And no parent will be mssing for sibs. 

=head1 REQUIRED ARGUMENTS

=over 

=item -[-]ped [=] <pedfile>

Pedigree file. 

=for Euclid:
	pedfile.type: readable

=item -[-]out[put] <outfile>

The output patched PED file name.

=back

=head1 OPTIONS

=over

=item -[-]sibs [=] <list>

The full list of full sib pairs.
If sib pairs are from different families, FamID needs to be fixed by fix_ped_vcf providing
those sib pairs and cryptic related pairs.

=item -[-]half-sibs [=] <list>

The full list of half-sib pairs.

=item -[-]pattern [=] <string>

Pattern for IDs of newly added samples, used by randregex() of String::Randome module.

=for Euclid:
	string.default: '[A-Z]{8}'

=item -[-]ignore [=] <string>

Regex pattern for duplicated sample IDs. They will not be considered as sibling.

In case when the original sample does not exist, the dup sample will be renamed.

=for Euclid:
	string.default: '_Re\d*$'

=item -[-]no-strict

Do not require all parents must also present in the IID column of PED file.
Otherwise, each newly added parent will add a new line to the PED output.

=back

