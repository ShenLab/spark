#!/usr/bin/env perl

use strict;
use warnings;
use Data::Dumper;
use Perl6::Slurp;

my ($geno, $label, $popstr, $output) = @ARGV;

unless(@ARGV == 4) {
	print STDERR "$0 GENO LABEL POPSTR OUTPUT\n";
	exit 1;
}
else {
	unless(-f $geno.".bim" && -f $geno.".bed" && -f $geno.".fam") {
		die "Cannot find all genotype related files $geno.{bim,bed,fam}";
	}
	unless(-f $label) {
		die "Cannot find label file $label";
	}

}

my %popu = map { (split)[0,1] } slurp $label;

my %keep = parse_group($popstr);
# Check to make sure all population subset has individual samples

my %fids = map { (split)[1,0] } slurp $geno.".fam";

my %known;
open my $fin, "$geno.fam" or die "Cannot open $geno.fam";
open my $fout, ">$output.keep" or die "Cannot write to $output.keep";
open my $flab, ">$output.label" or die "Cannot write to $output.labels";
while(<$fin>) {
	my ($fid, $iid, $dad, $mom, $sex, $phe) = split;
	die "Cannot find population label for $iid" unless defined $popu{$iid};
	if ( ($dad eq '0' || !defined $fids{$dad}) &&
		 ($mom eq '0' || !defined $fids{$mom}) ) {
		if ($keep{$popu{$iid}}) {
			print $fout $fids{$iid}, "\t", $iid, "\n";
			print $flab $iid, "\t", $keep{$popu{$iid}}, "\n";
			$known{$popu{$iid}} = 1;
		}
	}
}
close $fout; close $flab;

foreach my $popu (keys %keep) {
	unless (defined $known{$popu}) {
		die "No individuals from population $popu was found";
	}
}

system("plink --bfile $geno --keep $output.keep --make-bed --out $output");


sub parse_group {
	my ($fstr) = @_; 
	my %rename;
	my @fields = split(',', $fstr);
	foreach my $fstr (@fields) {
		if ($fstr =~ /^(\w+)\[(\w+)\]$/) {
			my ($oldid, $newid) = ($1, $2);
			$rename{$oldid} = $newid;
		}
		else {
			$rename{$fstr} = $fstr;
		}
	}
	return %rename;
}

