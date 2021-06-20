#!/usr/bin/env perl

use strict;
use warnings;
use Graph::Undirected;
use List::MoreUtils qw|any|;
use Utils::List qw|which_max|;
use Utils::File qw|count_line|;
use Utils::File::Iter qw|iter_file slurp_file|;
use Genet::Ped;

my ($pedfile, $outdir) = @ARGV;

unless(@ARGV == 2) {
	print STDERR <<EOF

Purpose: Used by samp_qc to create patch files that for fix_vcf_ped

Usage: create_patchs PEDFILE OUTDIR

It will create the following files in OUTDIR:

patch.badsamps.txt - bad samples
patch.sex_update.txt - sex updates
patch.nonpar.txt - incorrect parents
patch.swapped.txt - swapped samples
patch.dups.txt   - duplications
patch.rels.txt   - cryptic relatedness

Processing order for PED: remove badsampes, update sex/nonpar, rename/remove duplicates, patch relatedness

Processing order for VCF: remove badsampes and duplicates, rename swapped samples, rename duplicates if not removed

The resulting files should be reviewed manually to create the final patch files.

EOF
}
else {
	die "Cannot find ped file" unless -f $pedfile;
	die "Cannot find output directory" unless -d $outdir;
}

my $ped = Genet::Ped->new($pedfile, {verbose => 0});

my %out;
if (-f "$outdir/sample_err.csv") {
	my $it = iter_file("$outdir/sample_err.csv");
	while(my $dat = $it->()) {
		if ($dat->{LOWDP} eq 'TRUE') {
			$out{badsamps}{$dat->{IID}} = 1;
		}
		elsif ($dat->{SEXERR}) {
			$out{sex_update}{$dat->{IID}} = $dat->{Sex_pred};
		}
	}
}
else {
	print STDERR "Not error samples found\n";
}

my %dp;
{
	my $it = iter_file("$outdir/sample.csv");
	while(my $dat = $it->()) {
		# Note: previously we used DP_het, but now realized that GATK's downsamping causes inaccurate estimate of DP
		# DP_auto also has this issue but less severe
		$dp{$dat->{IID}} = $dat->{DP_auto};
	}
}


foreach my $fid ($ped->get_famids()) {
	my %pars = $ped->get_parents($fid);
	while(my ($iid, $pair) = each %pars) {
		if (@$pair == 2) {
			if (exists $out{sex_update}) {
				if (defined $out{sex_update}{$pair->[0]} && defined $out{sex_update}{$pair->[1]} &&
					$out{sex_update}{$pair->[0]} ne $out{sex_update}{$pair->[1]}) {
					# first type of swapping error: parent swap, both gender error in father & mother
					push @{$out{swapped}} => [ "Family $fid: sample $iid 's Parents swap", $pair->[0], $pair->[1] ];
				}
			}
		}
	}
}
if (defined $out{badsamps}) {
	my $fout = IO::File->new("$outdir/patch.badsamps.txt", "w");
	foreach my $iid (keys %{$out{badsamps}}) {
		print $fout $iid, "\n";
	}
}
else {
	unlink "$outdir/patch.badsamps.txt";
}
if (defined $out{sex_update}) {
	my $fout = IO::File->new("$outdir/patch.sex_update.txt", "w");
	while(my ($iid, $sex) = each %{$out{sex_update}}) {
		print $fout $iid, "\t", $sex, "\n";
	}
}
else {
	unlink "$outdir/patch.sex_update.txt";
}
if (defined $out{swapped}) {
	my $fout = IO::File->new("$outdir/patch.swapped.txt", "w");
	foreach my $pair (@{$out{swapped}}) {
		print $fout '#'.$pair->[0],"\n";
		print $fout $pair->[1], " ", $pair->[2], "\n";
	}
}
else {
	unlink "$outdir/patch.swapped.txt";
}



if (-f "$outdir/relpairs_err.csv") {
	my $dup = Graph::Undirected->new();
	my $rel = Graph::Undirected->new();
	my @datline = sort { $a->{FID1} cmp $b->{FID1} }
		slurp_file("$outdir/relpairs_err.csv");
	
	my @fields = qw|FID1 ID1 TYP1 FID2 ID2 TYP2 Kinship IBS0 PI_HAT Z0 Z1 Z2|;
	my $fdup = IO::File->new("$outdir/relpairs_err_dup.txt", "w");
	print $fdup join("\t", @fields), "\n";
	my $fpo  = IO::File->new("$outdir/relpairs_err_po.txt", "w");
	print $fpo join("\t", 'ERRTYPE', @fields), "\n";
	my $fsib = IO::File->new("$outdir/relpairs_err_sib.txt", "w");
	print $fsib join("\t", 'ERRTYPE', @fields), "\n";
	my $frel = IO::File->new("$outdir/relpairs_err_rel.txt", "w");
	print $frel join("\t", @fields), "\n";
	my $ffat = IO::File->new("$outdir/relpairs_err_fatal.txt", "w");
	print $ffat join("\t", @fields), "\n";

	my (%errpairs, %errindv, %nonpar, %dupcat);
	foreach my $dat (@datline) {
		foreach my $II (qw|1 2|) {
			my $isfd = $ped->is_founder($dat->{"FID$II"}, $dat->{"ID$II"});
			my $sex = $ped->get_sex($dat->{"FID$II"}, $dat->{"ID$II"});
			if ($isfd) {
				if (!defined $sex) {
					$dat->{"TYP$II"} = "Parent";
				}
				elsif ($sex eq '1') {
					$dat->{"TYP$II"} = "Father";
				}
				else {
					$dat->{"TYP$II"} = "Mother";
				}
			}
			else {
				if (!defined $sex) {
					$dat->{"TYP$II"} = "Offspring";
				}
				elsif ( $sex eq '1') {
					$dat->{"TYP$II"} = "Son";
				}
				else {
					$dat->{"TYP$II"} = "Daughter";
				}
			}
		}
		if (any { $dat->{$_} eq 'TRUE' } qw|POERR WFPO BFPO SIBERR WFSIB BFSIB FATAL|){
			$errindv{$dat->{"ID1"}} = 1;
			$errindv{$dat->{"ID2"}} = 1;
		}

		if ($dat->{DUPERR} eq 'TRUE') {
			print $fdup join("\t", @{$dat}{@fields}), "\n";
			push @{$errpairs{DUP}} => [$dat->{ID1}, $dat->{ID2}];
		}
		if (any { $dat->{$_} eq 'TRUE' } qw|POERR WFPO BFPO|) {
			my $errtype = join(",", grep { $dat->{$_} eq 'TRUE' } qw|POERR WFPO BFPO|);
			print $fpo join("\t", $errtype, @{$dat}{@fields}), "\n";
			if ($errtype eq 'POERR') {
				if ($dat->{TYP1} eq 'Father' || $dat->{TYP1} eq 'Mother') {
					my $slot = $dat->{TYP1} eq 'Father' ? 'DAD' : 'MOM';
					$nonpar{$dat->{ID2}}{$slot} = '0';
				}
				elsif ($dat->{TYP2} eq 'Father' || $dat->{TYP2} eq 'Mother') {
					my $slot = $dat->{TYP2} eq 'Father' ? 'DAD' : 'MOM';
					$nonpar{$dat->{ID1}}{$slot} = '0';
				}
			}
		}
		if (any { $dat->{$_} eq 'TRUE' } qw|SIBERR WFSIB BFSIB|) {
			my $errtype = join(",", grep { $dat->{$_} eq 'TRUE' } qw|SIBERR WFSIB BFSIB|);
			print $fsib join("\t", $errtype, @{$dat}{@fields}), "\n";
		}
		if (any { $dat->{$_} eq 'TRUE' } qw|WFREL BFREL|) {
			print $frel join("\t", @{$dat}{@fields}), "\n";
			push @{$errpairs{REL}} => [$dat->{ID1}, $dat->{ID2}];
		}
		if ($dat->{FATAL} eq 'TRUE') {
			print $ffat join("\t", @{$dat}{@fields}), "\n";
		}
	}

	# All duplicated pair and all cryptic related pairs
	if (defined $errpairs{DUP}) {
		foreach my $pair (@{$errpairs{DUP}}) {
			$dup->add_edge($pair->[0], $pair->[1]);
		}
		my $fout = IO::File->new("$outdir/patch.dups.txt", "w");
		foreach my $cc ($dup->connected_components) {	
			my @noerr = grep { !defined $errindv{$_} } @$cc;
			if (@noerr == 1) {
				my @outlst = @noerr;
				push @outlst, grep { $_ ne $noerr[0] } @$cc;
				my @dps = map { $dp{$_} } @outlst;
				# find sample with highest depth
				my $ii = which_max(@dps);
				$outlst[$ii] = "[".$outlst[$ii]."]";
				print $fout join("\t", @outlst), "\n";
			}
			else {
				print $fout join("\t", sort { $dp{$b} <=> $dp{$a} } @$cc), "\n";
			}
			foreach (@$cc) {
				$dupcat{$_} = 1;
			}
		}
	}
	else {
		unlink "$outdir/patch.dups.txt";
	}

	if (defined $errpairs{REL}) {
		foreach my $pair (@{$errpairs{REL}}) {
			$rel->add_edge($pair->[0], $pair->[1]);
		}
		my $fout = IO::File->new("$outdir/patch.rels.txt", "w");
		foreach my $cc ($rel->connected_components) {
			print $fout join("\t", @$cc), "\n";
		}
	}
	else {
		unlink "$outdir/patch.rels.txt";
	}

	foreach my $suf (qw|dup po sib rel fatal|) {
		my $file = "$outdir/relpairs_err_$suf.txt";
		unlink $file if count_line($file) == 1;
	}

	# Analyze all PO errors
	if (%nonpar) {
		my $fout = IO::File->new("$outdir/patch.nonpar.txt", "w");
		while(my ($iid, $info) = each %nonpar) {
			# We will not output the cases with dup err?
			next if defined $dupcat{$iid};
			while(my ($pa, $iid2) = each %$info) {
				print $fout join("\t", $iid, $pa), "\n";
			}
		}
	}
	else {
		unlink "$outdir/patch.nonpar.txt";
	}
}




	


