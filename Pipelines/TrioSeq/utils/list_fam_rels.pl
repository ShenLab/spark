#!/usr/bin/env perl

use strict;
use warnings;
use FindBin qw|$Bin|;
use Perl6::Slurp;
use List::MoreUtils qw|uniq all|;
use Getopt::Euclid;

use lib "$Bin/../../lib";
use Shared qw|fam_rels reltype fam_kins|;


my $twins;
if ($ARGV{'--twins'}) {
	my @twins;
	open my $fin, "$ARGV{'--twins'}" or die "Cannot open twins";
	while(<$fin>) {
		next if /^#/ || /^\s*$/;
		my @dat = split;
		unless (@dat == 2) {
			die "Incorrect number of columns in twins file: $_";
		}
		push @twins, \@dat;
	}
	$twins = \@twins;
}

my %famid = map { (split)[1,0] } slurp $ARGV{'--ped'};

my (@samps, $pairflag);
if ($ARGV{'--samps'}) {
	if (-f $ARGV{'--samps'} || $ARGV{'--samps'} eq '-') {
		# @samps = map { (split)[0] } slurp $ARGV{'--samps'}; 
		my $fin;
		if ($ARGV{'--samps'} eq '-') {
			$fin = \*STDIN;
		}
		else {
			open $fin, $ARGV{'--samps'} or die "Cannot open samp list";
		}
		my $ncol;
		while(<$fin>) {
			my @a = split;
			unless (@a == 1 || @a == 2) {
				die "Incorrect number of columns in sample list: $_";
			}
			unless (defined $ncol) {
				$ncol = scalar(@a);
			}
			else {
				unless(@a == $ncol) {
					die "Inconsistent number of columns in sample list: $_";
				}
			}
			if ($ncol == 1) {
				if (defined $famid{$a[0]}) {
					push @samps, $a[0];
				}
				else {
					warn "Sample $a[0] cannot be found in the pedigree, skip\n";
				}
			}
			else {
				if (defined $famid{$a[0]} && defined $famid{$a[1]}) {
					if ($famid{$a[0]} eq $famid{$a[1]}) {
						push @samps, \@a;
					}
					else {
						warn "Sample pair $a[0],$a[1] are not in the same family, skip\n";
					}
				}
				else {
					warn "Sample pair $a[0],$a[1] cannot be found in the pedigree, skip\n";
				}
				
			}
		}
		$pairflag = 1 if $ncol == 2;
	}
	else {
		foreach my $iid (split(',', $ARGV{'--samps'})) {
			if (defined $famid{$iid}) {
				push @samps, $iid;
			}
			else {
				warn "Sample $iid cannot be found in the pedigree, skip\n";
			}
		}
	}
}
else {
	@samps = map { (split)[1] } slurp $ARGV{'--ped'};
}

#if ($pairflag) {
#	unless(all { defined $famid{$_->[0]} && defined $famid{$_->[1]}  } @samps) {
#		my @sampswofid = uniq sort grep { !defined $famid{$_} } map { @$_ } @samps;
#		print STDERR "Samples w/o family IDs: \n", join("\n", @sampswofid), "\n";
#		die "Not all sample pairs in the list can be found in PED file!";
#	}
#}
#else {
#	unless (all { defined $famid{$_} } @samps) {
#		die "Not all samples in the list can be found in the PED file";
#	}
#}
my %opt;
if ($ARGV{'--samps'}) {
	if ($pairflag) {
		$opt{select} = [uniq sort map { ($famid{$_->[0]}, $famid{$_->[1]}) } @samps];
	}
	else {
		$opt{select} = [uniq sort map { $famid{$_} } @samps];
	}
}

my ($famsamp, $famrels) = fam_rels($ARGV{'--ped'}, 
	{ twins => $twins, strict => $ARGV{'--strict'}, ignore => $ARGV{'--ignore'}, 
	  shorten => 1, fullsib => 1, verbose => $ARGV{'--verbose'}, %opt });


my $idcoefs;
if ($ARGV{'--calc'}) {
	$idcoefs = fam_kins($ARGV{'--ped'}, 
		{ twins => $twins, strict => $ARGV{'--strict'}, ignore => qr/$ARGV{'--ignore'}/,
		  full => undef, verbose => $ARGV{'--verbose'}, wrkdir => $ARGV{'--wrkdir'}, %opt });
}

my $fout;
if (defined $ARGV{'--output'}) {
	open $fout, ">$ARGV{'--output'}" or die "Cannot write to $ARGV{'--output'}";
}
else {
	$fout = \*STDOUT;
}

my @outfields = qw|FamID IID1 IID2 Relation|;
if ($ARGV{'--reltype'}) {
	push @outfields, "RelType";
}
if ($ARGV{'--calc'}) {
	push @outfields, qw|Kinship ProbIBD2 ProbIBD1|;
}

print $fout join("\t", @outfields), "\n";
if ($pairflag) {
	foreach my $pair (@samps) {
		next unless $famid{$pair->[0]} eq $famid{$pair->[1]};
		my $fid = $famid{$pair->[0]};
		my @outline = ($fid, $pair->[0], $pair->[1], $famrels->{$pair->[0]}{$pair->[1]});
		if ($ARGV{'--reltype'}) {
			push @outline, reltype($famrels->{$pair->[0]}{$pair->[1]});
		}
		if ($ARGV{'--calc'}) {
			if (defined $idcoefs->{$pair->[0]}{$pair->[1]}) {
				push @outline, @{$idcoefs->{$pair->[0]}{$pair->[1]}};
			}
			else {
				push @outline, (".") x 3;
			}
		}
		print $fout join("\t", @outline), "\n";
	}
}
else {
	foreach my $iid (@samps) {
		my $rels = $famrels->{$iid};
		foreach my $relid (sort keys %$rels) {
			my @outline = ($famid{$iid}, $iid, $relid, $rels->{$relid});
			if ($ARGV{'--reltype'}) {
				push @outline, reltype($rels->{$relid});
			}
			if ($ARGV{'--calc'}) {
				my $fid = $famid{$iid};
				if (defined $idcoefs->{$iid}{$relid}) {
					push @outline, @{$idcoefs->{$iid}{$relid}};
				}
				else {
					push @outline, (".") x 3;
				}
			}
			print $fout join("\t", @outline), "\n";
		}
	}
}



__END__

=head1 NAME

list_fam_rels -- Listing all family members for each individual.

=head1 USAGE

list_fam_rels [options] -in PED -out OUTPUT -list SAMPS ...

=head1 DESCRIPTION

The script is mainly used for testing fam_rels function.

Given a pedigree file, and a list of samples that exist in the pedigree.
For each sample, list all his/her family members and show their relationship.
The output will contain four columns: FamID, IID1, IID2, Relation (Kinship, ProbIB2, ProbIBD1)
Relation is the familial relationship of IID2 w.r.t to IID1.

=head1 REQUIRED ARGUMENTS

=over 

=item -[-]ped [=] <pedfile>

Pedigree file.

=for Euclid:
	pedfile.type: readable

=back

=head1 OPTIONS

=over

=item -[-][rel]type

Also report reltype up to 2nd degree.

=item -[-]calc

Also calculate kinship coeffs.

=item -[-]out[put] <outfile>

The output file name, write to STDOUT if not providd.

=item -[-]samp[s] [=] <list>

The sample list, otherwise use all samples in the PED file.
Sample list can be a comma separated list, or an one column text file containing sample
list, or a two column text file containing sample pairs.

=item -[-]ignore [=] <string>

Regex pattern for duplicated sample IDs (default: '_Re\d*$').

=for Euclid:
	string.default: '_Re\d*$'

=item -[-]twins [=] <list>

A list of twin pairs. Should be a file with two columns.

=for Euclid:
	list.type: readable

=item -[-]wrkdir [=] <dir>

Working directory for idcoefs.

=item -[-]strict

Under strict mode, both parents for non-founders must be specified.

=item -[-]verbose

Verbose mode.

=back

