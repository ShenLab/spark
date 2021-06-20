use strict;
use warnings;
use FindBin qw|$Bin|;
use List::Util qw|product first|;
use List::MoreUtils qw|any all none|;
use Data::Dumper;
use Storable qw|dclone|;
use Getopt::Euclid;
use FaSlice;
use Utils::List qw|parse_fields|;
use Utils::File::Iter qw|iter_file|;
use Genet::Var qw|normalize|;

use lib "$Bin/../../lib";
use Shared qw|expand_dat parse_tabfile|;
use Variants qw|var_type var_dist cluster_vars|;

# Need to sort table files by IID, then by OrigVarID, then by Distance
# IID,Chrom,Position,Ref,Alt,OrigVarID,Distance,PGT
my @keyfields = parse_fields($ARGV{'--fields'});
unless(@keyfields == 8) {
	die "Must provide key fields: IID,Chrom,Position,Ref,Alt,OrigVarID,Distance,PGT!";
}

my @colns;
{
	my ($it, $fnames) = iter_file($ARGV{'--tab'}, { fsep => qr/\t/ });
	foreach my $field (@keyfields) {
		my $num = first { $fnames->[$_-1] eq $field } 1..@{$fnames};
		unless(defined $num) {
			die "Cannot find required column $field";
		}
		push @colns, $num;
	}
}

open my $fin, "cat $ARGV{'--tab'} | body sort -k$colns[0],$colns[0] -k$colns[5],$colns[5] -k$colns[6],$colns[6]n |"
	or die "Cannot open body sort pipe";

#my ($it, $fnames, $keyfields) = parse_tabfile($ARGV{'--tab'}, $ARGV{'--fields'}, 8, 8);
my ($it, $fnames) = iter_file($fin, { fsep => qr/\t/ });

my @outfields = @$fnames;
if ($ARGV{'--fields-extra'}) {
	foreach my $field (split(',', $ARGV{'--fields-extra'})) {
		unless(grep { $field eq $_ } @$fnames) {
			die "Cannot find field $field in the input file";
		}
		push @outfields, $field;
	}
}

my $fa = FaSlice->new(file => $ARGV{'--fasta'});

my ($fout, $ferr);
if (defined $ARGV{'--output'}) {
	$ARGV{'--output'} =~ s/\.txt$//;
	open $fout, ">$ARGV{'--output'}.txt" or die "Cannot write to $ARGV{'--output'}.txt";
	open $ferr, "| uniq -c > $ARGV{'--output'}.missing.txt" or die "Cannot write to $ARGV{'--output'}.missing.txt"
} 
else {
	$fout = \*STDOUT;
	$ferr = \*STDERR;
}
print $fout join("\t", @outfields), "\n";


my %clsted;
my (@cluster);
my ($prev_varid, $prev_iid) = ("") x 2;
my $ncands = 0;
while(my $dat = $it->()) {
	# $varid here is the original VarID for variant clustering
	my ($iid, $chrom, $pos, $ref, $alt, $varid, $dist, $pgt) = @{$dat}{@keyfields};
	next if defined $clsted{$iid,$varid};

	# The start of a new MNV, either start a new varid or dist = 0
	# In case when varid is new but dist > 0, it should be skipped and recorded in err log
	if ( $varid ne $prev_varid || $iid ne $prev_iid || $dist == 0) {
		if (($varid ne $prev_varid || $iid ne $prev_iid) && $dist != 0) {
			print $ferr join("\t", @{$dat}{@keyfields}), "\n";
			next;
		}
		if (@cluster > 1) {
			output_cluster(@cluster);
		}
		# The first variant in the cluster must have dist==0
		unless($dist == 0) {
			die "The start of new variant cluster must have dist==0!";
		}
		@cluster = (dclone $dat);
	}
	else {
		# Need to check phase before adding variant to the cluster
		if (any { var_type($_->{$keyfields[3]}, $_->{$keyfields[4]}) eq 'Indel' } (@cluster, $dat)) {
			if (any { var_dist($pos, $ref, $_->{$keyfields[2]}, $_->{$keyfields[3]}) 
						<= $ARGV{'--dist-indel'}  } @cluster) {
				push @cluster, dclone $dat if coupled(@cluster, $dat);
			}
		}
		else {
			if (any { var_dist($pos, $ref, $_->{$keyfields[2]}, $_->{$keyfields[3]}) 
						<= $ARGV{'--dist-snv'}  } @cluster) {
				push @cluster, dclone $dat if coupled(@cluster, $dat);
			}
		}
	}
	$prev_varid = $varid;
	$prev_iid = $iid;
}

if (@cluster > 1) {
	output_cluster(@cluster);
}
close $fout;
close $ferr;

if ($ncands > 0) {
	print STDERR "Found a total of $ncands MNV candidates\n";
}
else {
	print STDERR "No MNV candidates was found\n";
	unlink "$ARGV{'--output'}.txt";
}

if (-z "$ARGV{'--output'}.missing.txt") {
	unlink "$ARGV{'--output'}.missing.txt";
}

# Test if variants in the cluster are coupled
sub coupled {
	my (@vars) = @_;
	# All variants must be phased or homozygotes alt (multi-allelic not supported)
	# The necessary condition is suitable for GATK RBP
	unless (all { $_->{$keyfields[7]} =~ /^[01]\|[01]$/ || $_->{$keyfields[7]} eq '1/1' } @vars) {
		return 0;
	}
	my (@hap1, @hap2);
	foreach my $var (@vars) {
		if ($var->{$keyfields[7]} eq '1/1') {
			push @hap1, 1;
			push @hap2, 1;
		}
		else {
			my @als = split(q{\|}, $var->{$keyfields[7]});
			die "$var->{$keyfields[7]} is not biallelic" unless @als == 2;
			unless(all { $_ == 0 || $_ == 1 } @als) {
				die "$var->{$keyfields[7]} is not a 0/1 genotype string";
			}
			push @hap1, $als[0];
			push @hap2, $als[1];
		}	
	}
	if (product(@hap1) == 0 && product(@hap2) == 0) {
		return 0;
	}
	else {
		return 1;
	}
}

sub output_cluster {
	my @cluster = @_;
	my $iid =  $cluster[0]{$keyfields[0]};
	my $varid = $cluster[0]{$keyfields[5]};

	# All variants in the cluster should have the same VarID and IID
	unless(all { $_->{$keyfields[5]} eq $varid } @cluster) {
		print STDERR join(", ", map { $_->{$keyfields[5]} } @cluster), "\n";
		die "Not all variants in the cluster have the same OrigVarID";
	}
	unless(all { $_->{$keyfields[0]} eq $iid } @cluster) {
		print Dumper \@cluster;
		die "Not all variants in the cluster have the same IID";
	}

	my @clstord = sort { $a->{Position} <=> $b->{Position} } @cluster;
	
	if ($ARGV{'--enum'}) {
		# Ouptut all possible combination of consecutive SNVs
		# 1,2,3,4: (1,2)(2,3)(3,4)(1,2,3)(2,3,4)(1,2,3,4) length: 2~tot
		for(my $len = 2; $len <= @clstord; $len ++) {
			for(my $ii = 0; $ii < scalar(@clstord)-$len+1; $ii ++) {
				my @mnv = cluster_vars($fa, map { [ @{$_}{@{keyfields}[1..4]} ] } @clstord[$ii..$ii+$len-1]);
				unshift @mnv, $iid;
				push @mnv, $varid, ".", ".";
				my %merged;
				@merged{@keyfields} = @mnv;
				print $fout join("\t", map { $_ // "." } @merged{@outfields}), "\n";
			}
		}
	}
	else {
		my @mnv = cluster_vars($fa, map { [ @{$_}{@{keyfields}[1..4]} ] } @clstord);
		unshift @mnv, $iid;
		push @mnv, $varid, ".", ".";
		my %merged;
		@merged{@keyfields} = @mnv;
		print $fout join("\t", map { $_ // "." } @merged{@outfields}), "\n";
	}
	foreach my $clsvar (@cluster) {
		my $vid = join(":", @{$clsvar}{@{keyfields}[1..4]});
		$clsted{$iid,$vid} = 1;
		print $fout join("\t", @{$clsvar}{@outfields}), "\n";
	}
	$ncands ++;
}


__END__

=head1 NAME

list_mnv_cands.pl -- Extract candidate MNVs for re-annotation.

=head1 NOTE

MNV candidates are generated from variants that are on the same haplotye.

We will collect nearby variants that are on the same haplotypes.

The input file are prepared by tabulating nearby variants which are collected from
dv_eval output by list_nearvars.pl.

=head1 REQUIRED ARGUMENTS
 
=over

=item -[-]tab[le] [=] <input>

The original input varinat table.

=item -[-]fa[sta] [=] <file>

Path to genome reference sequence file in fasta format.

=back

=head1 OPTIONS

=over

=item -[-]out[put] [=] <file>

The output file name.

=item -[-]fields [=] <string>

Standard fields for sample variant and nearby variants

=for Euclid:
	string.default: 'IID,Chrom,Position,Ref,Alt,OrigVarID,Distance,PGT'

=item -[-]fields-extra [=] <string>

Extra fields from the input file to be added to the output.

=item -[-]dist-snv [=] <distance>

Distance (bp) cutoff for nearby SNVs that forms part of candidate MNVs, default: 2.

=for Euclid:
	distance.default: 2

=item -[-]dist-indel [=] <distance>

Distance cutoff to indels that form complex indels/MNVs, default: 10.

=for Euclid:
	distance.default: 10

=item -[-]enum

Enumerate all possible combination of nearby variants.

=back