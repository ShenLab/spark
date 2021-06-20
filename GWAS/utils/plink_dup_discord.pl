#!!/usr/bin/perl

use strict;
use warnings;
use Carp;
use Data::Dumper;
use Graph::Undirected;
use IO::File;
use File::Path qw|make_path|;
use File::Temp qw|tempdir|;
use File::Copy qw|move copy|;
use List::MoreUtils qw|uniq all|;
use Perl6::Slurp;
use Utils::List qw|all_pairs|;
use Utils::File::Iter qw|iter_file|;
use Getopt::Euclid;


# Check input file
unless(all { -f "$ARGV{'--input'}.$_" } qw|fam bed bim|) {
	croak "Cannot find all necessary input files";
}


my @dup_edges;
if ($ARGV{'--nofid'}) {
	my %fids = map { (split)[1,0] } slurp "$ARGV{'--input'}.fam";
	open my $fin, $ARGV{'--dups'} or die "Cannot open dups file: $ARGV{'--dups'}";
	while(<$fin>) {
		my @iids = split;
		unless(all { defined $fids{$_} } @iids) {
			croak "Cannot find FIDs for $_";
		}
		push @dup_edges, [ map { "$fids{$_} $_" } @iids ];	
	}
}
else {
	open my $fin, $ARGV{'--dups'} or die "Cannot open dups file: $ARGV{'--dups'}";
	while(<$fin>) {
		my @fiids = split;
		if(@fiids % 2) {
			croak "FID and IID are not paired: $_";
		}
		my $nsamp = @fiids/2;
		push @dup_edges, [ map { $fiids[2*$_-2]." ".$fiids[2*$_-1] } 1..$nsamp ];
	}
}


# create graph to represent relationship between duplicated samples
my $gdup = Graph::Undirected->new;
for my $dups (@dup_edges) {
	$gdup->add_edge(@$_) foreach (all_pairs(@$dups));
}

# Find sub graphs of size 2,3,4...
my @cc = $gdup->connected_components();
my %dupcc; 
push @{$dupcc{scalar @$_}}, $_ foreach (@cc);
my @all_pairs = pairwise(\%dupcc);


my $tmpdir;
if ($ARGV{'--wrkdir'}) {
	$tmpdir = $ARGV{'--wrkdir'};
	make_path $tmpdir unless -d $tmpdir;
}
else {
	$tmpdir = tempdir(CLEANUP => 1);
}


my %discord;
foreach my $pair (@all_pairs) {
	my $prefix = $tmpdir."/".join("_", @$pair);
	$prefix =~ s/\s+/-/g;
  	# retrieve duplicates
  	foreach my $ii (0,1) {
  		my $fh = IO::File->new("$prefix.$ii.keep","w");
  		print $fh $pair->[$ii], "\n";
  		$fh->close;
  		system("plink --bfile $ARGV{'--input'} --keep $prefix.$ii.keep --make-bed --out $prefix.$ii");
  	}
  	# update individual IDs on the second file
  	move("$prefix.1.fam", "$prefix.1.fam.bak");
  	copy("$prefix.0.fam", "$prefix.1.fam");

  	# calculate discordant genotype calls using plink merge mode 6/7
  	my $mode = 7;
  	if($ARGV{'--allmis'}) {
  		$mode = 6;
  	}
  	system("plink --bfile $prefix.0 --bmerge $prefix.1.bed $prefix.1.bim $prefix.1.fam --merge-mode $mode --out $prefix");
  	my $it = iter_file("$prefix.diff", {chomp => 3});
  	while(my $dat = $it->()) {
  		$discord{$dat->{SNP}} ++;
  	}
}

my $fout;
if ($ARGV{'--output'}) {
	open $fout, ">$ARGV{'--output'}" or die "Cannot open output file";
}
else {
	open $fout, ">$ARGV{'--input'}.duperr" or die "Cannot write to duperr file";
}
print $fout "SNP\tDUPERR\n";
while ( my ($rsid, $num) = each %discord ) {
  print $fout "$rsid\t$num\n";
}


# pairing of duplicated individuals
sub pairwise {
	my $dupcc = shift;
	my @all_pairs;
	while ( my ($size, $subcc) = each %$dupcc ) {
		for my $ii (0..$size-1) {
			for my $jj ($ii+1..$size-1) {
				my @one_pair;
				foreach my $cc (@$subcc) {
					push @one_pair, [$cc->[$ii], $cc->[$jj]];
				}
				push @all_pairs, @one_pair;
			}
		}
	}
	return @all_pairs;
}


__END__

=head1 NAME

plink_dup_discord.pl -- Calculate number of discordant genotypes for duplicated samples.

=head1 DESCRIPTION

This script made use of plink's merge function to calculate discordance rate of genotypes.

=head1 REQUIRED ARGUMENTS

=over

=item -[-]in[put] [=] <input_prefix>

Input data file, should be in plink bped format.

=item -[-]dup[s] [=] <duplist>

A list of duplicated sample, each line contains FID IID of a pair of duplicated samples.
Or it can also be contain more than two individuals. If FID is not provided in this file,
it should be retrived from fam file.

=back

=head1 OPTIONS

=over

=item -[-]out[put] [=] <output_file>

Default output file will be input.duperr. Can specify a different name using this option.

=item -[-]wrkdir [=] <wrkdir>

Keep the intermediate files in the given directory.

=item -[-]allmis

Report all mismatching call, regardless of missingness

=item -[-]nofid

Indicating that the no FID can be found in the dups file.

=back

=cut
