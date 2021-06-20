#!/usr/bin/env perl

use strict;
use warnings;
use Carp;
use Getopt::Euclid;
use Data::Dumper;
use List::Util qw|first|;
use List::MoreUtils qw|all|;
use Utils::File::Iter qw|iter_file|;
use Utils::List qw|replace_all|;
use Genome::UCSC::BinKeeper;

use lib "$ENV{HOME}/Pipelines/lib/";
use Shared qw|fam_rels|;


# Parse pedigree file
# Relationships between samples that appear in the iuput table will be updated 
# to keep consistent with the PED file.
# All samples must be found in PED file!
my ($famsamp, $famrels) = fam_rels($ARGV{'--ped'}, 
				{ twins => $ARGV{'--twins'}, shorten => 1, strict => 0,
				  ignore => $ARGV{'--ignore'}, verbose => 0 });
my %pedinfo;
{
	open my $fin, $ARGV{'--ped'} or die "Cannot open PED file $ARGV{'--ped'}";
	while(<$fin>) {
		my ($famid, $iid, $dad, $mom, $sex, $phe) = split;
		$pedinfo{$iid} = { FamID => $famid,
			Sex => $sex eq '1' ? 'Male' : $sex eq '2' ? 'Female' : 'NA',
			Pheno => $phe eq '1' ? 'Unaffected' : $phe eq '2' ? 'Affected' : 'NA' };
	}
}

# Get input file fields
my ($it, $fnames) = iter_file($ARGV{'--input'}, { fsep => qr/\t/ });
foreach my $stdfd (qw|Chrom Position Ref Alt|) {
	unless(grep { $stdfd eq $_ } @$fnames) {
		die "Cannot find standard field $stdfd from input!";
	}
}

my @outfields = @$fnames;

# Primary fields for the variant carrier, family members, and family carriers
# For Variant carrier, should be FamID IID/Gender/Pheno or alias provided in --fields option
# For family members, should be XXX.FamMembers (FamMembers is the default suffix)
my @keys;
{
	@keys = split(',', $ARGV{'--fields'});
	unless(@keys == 3) {
		die "Incorrect number of primary fields for the CNV carrier!";
	}
	unless(grep { $_ eq $keys[0] } @$fnames) {
		die "Cannot find CNV carrier ID field $keys[0]!"
	}
	for(my $ii = 1; $ii < @keys; $ii ++) {
		unless(grep { $_ eq $keys[$ii] } @$fnames) {
			warn "Primary field $keys[$ii] will be added";
			insert_after(\@outfields, $keys[$ii-1], $keys[$ii]);
		}
	}
}


my @infofields;
my $suffix = $ARGV{'--suffix'};
my $geno = $ARGV{'--geno'};
# Info fields will be shared by carrier and family members with similar field names.
# XXX and XXX.FamMembers
{
	foreach my $stdfd (qw|FamMembers Relations Phenotypes|) {
		unless(grep { $stdfd eq $_ } @$fnames) {
			die "Cannot find standard field $stdfd";
		}
	}
	foreach my $memfd (grep { /\.$suffix$/ } @$fnames) {
		(my $origfd = $memfd) =~ s/\.$suffix$//;
		if (grep { $origfd eq $_ } @$fnames) {
			push @infofields, $origfd;
		}
	}
	unless(grep { $_ eq $geno } @infofields) {
		die "Cannot find genotype field $geno";
	}
	#if (@infofields > 0) {
		print STDERR "The following fields will be updated for both carrier and family members:\n";
		print STDERR join(",", @infofields), "\n";
	#}
}

my $fout;
if ($ARGV{'--output'}) {
	open $fout, ">$ARGV{'--output'}";
}
else {
	$fout = \*STDOUT;
}
print $fout join("\t", @outfields), "\n";

my $pattern;
if ($ARGV{'--relation'}) {
	$pattern = qr/$ARGV{'--relation'}/;
} 
else {
	$pattern = qr/^(\d*1st|\d*2nd|\d*3rd|\d*[^123]th)?(Daughter|Son)$/;
}
my $carrier;
if ($ARGV{'--carrier'}) {
	$carrier = qr/$ARGV{'--carrier'}/;
}
else {
	$carrier = qr/Het(?!O)|Alt/;
}

my %known;
while(my $dat = $it->()) {
	next if $dat->{FamMembers} eq '.';
	my $iid = $dat->{$keys[0]};
	my $varid = join(":", @{$dat}{qw|Chrom Position Ref Alt|});
	my @memberIDs = split(',', $dat->{FamMembers});
	my @memberGenos = split(',', $dat->{$geno.".FamMembers"});

	for(my $ii = 0; $ii < @memberIDs; $ii ++) {
		my $memID = $memberIDs[$ii];
		next if $known{$memID,$varid};

		my $relation = $famrels->{$iid}{$memID};
		unless(defined $relation) {
			die "Cannot find relationship between $iid and carrier $memID";
		}
		next unless $relation =~ /$pattern/;

		if ($memberGenos[$ii] =~ /$carrier/) {
			my $memInfo = $pedinfo{$memID};
			unless(defined $memInfo) {
				die "Cannot find pedigree info for carrier $memID";
			}

			# Update primary fields for carrier
			$dat->{$keys[0]} = $memID;
			$dat->{$keys[1]} = $memInfo->{Sex};
			$dat->{$keys[2]} = $memInfo->{Pheno};

			# Update family member fields
			my @newMemberIDs = @memberIDs;
			my $nrep = replace_all(\@newMemberIDs, $memID, $iid);
			unless($nrep == 1) {
				die "Cannot correctly replace family members for carrier $memID";
			}
			$dat->{FamMembers} = join(',', @newMemberIDs);
			$dat->{Relations} = join(',', map { $famrels->{$memID}{$_} // 
									do { "Cannot find relationship between $memID and $_" } } @newMemberIDs);
			$dat->{Phenotypes} = join(',', map { $pedinfo{$_}{Pheno} } @newMemberIDs);

			# Update(swap) info fields
			my $jj = first { $memberIDs[$_] eq $memID } 0..$#memberIDs;
			foreach my $field (@infofields) {
				my @memVals = split(',', $dat->{$field.".FamMembers"});
				unless(@memVals == @memberIDs) {
					die "Incorrect number of values in $field.FamMembers for family members";
				}
				($dat->{$field}, $memVals[$jj]) = ($memVals[$jj], $dat->{$field});
				$dat->{$field.".FamMembers"} = join(',', @memVals);
			}

			print $fout join("\t", @{$dat}{@outfields}), "\n";
			$known{$memID,$varid} = 1;
		}
	}
}
		

__END__

=head1 NAME

extract_transvar.pl -- Extract transmitted in parents or inherited CNVs in cases,

=head1 DESCRIPTION

Extract transmitted variants in parents and reformat them as case variants, or extract inherited variants
in cases and reformat them as parent variants. 

Required (updated) fields: IID,(Gender),(Pheno),(FamMembers, Relations, Phenotypes)

All other family info fields "XXX.FamMembers" will be updated accordingly.

=head1 REQUIRED ARGUMENTS

=over

=item -[-]in[put] [=] <file>

Input variant table. Must be tab separated with required fields.

=for Euclid:
	file.type: readable

=item -[-]ped [=] <pedfile>

Pedigree file. It will be used to look for family members.
It must include all samples and their family members appearing in the variant table.

=for Euclid:
	pedfile.type: readable

=back

=head1 OPTIONS

=over

=item -[-]out[put] [=] <file>

Output file name. If not provided, will write to STDOUT.

=item -[-]fields [=] <string>

Fields names for the primary CNV carrier, default: IID,Gender,Pheno.
The corresponding member fields will be FamMembers/Relations/Phenotypes.

=for Euclid:
	string.default: "IID,Gender,Pheno"

=item -[-]suffix [=] <string>

Field name suffix for family member information.

=for Euclid:
	string.default: "FamMembers"

=item -[-]geno [=] <field>

Genotype field.

=for Euclid:
	field.default: "GENO"

=item -[-]ignore [=] <string>

Regex pattern for duplicated sample IDs.

=for Euclid:
	string.default: '_Re\d*$'

=item -[-]twins [=] <list>

A list of twin pairs. Should be a file with two columns.

=item -[-]relation [=] <regex>

Pattern of relationships to be extracted, default is to extract offspring.

=item -[-]carrier [=] <regex>

Pattern of genotypes of a variant carrier, default is to match Het or Hom of alt allele.

=back



