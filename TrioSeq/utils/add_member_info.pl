use strict;
use warnings;
use Data::Dumper;
use FindBin qw|$Bin|; 
use Getopt::Euclid;
use Utils::List qw|insert_after|;
use Utils::File::Iter qw|iter_file slurp_file|;


use lib "$Bin/../../lib";
use Shared qw|fam_rels parse_fstr parse_tabfile slurp_xref|;


my ($iter, $infields) = iter_file($ARGV{'--input'}, { fsep => qr/\t/ });

my $idfield;
{
	my $xtrafds = parse_fstr($ARGV{'--xref-fields'}, 1);
	my @xtrafds = values %$xtrafds;
	unless(@xtrafds > 1) {
		die "The external reference field does not have more than one fields!";
	}
	$idfield = $xtrafds[0];
}

my ($xtrainfo, $xtrafields) = slurp_xref($ARGV{'--xref'}, $ARGV{'--xref-fields'});

my $cnvflag;
# Check the input file 
unless(grep { $_ eq $idfield } @$infields) {
	die "Cannot find primary sample ID field $idfield from input!";
}
unless (grep { $_ eq 'FamMembers' } @$infields) {
	my $prefix = $ARGV{'--prefix'};
	if (grep { $_ eq $prefix."_IDs" } @$infields) {
		$cnvflag = $prefix;
		print STDERR "Input file is in filtered CNV format\n";
	}
	else {
		die "Cannot find sample IDs for family members from input";
	}
}

# Determine the output fields
my (@outself, @outmembers);
if ($cnvflag) {
	for(my $ii = 0; $ii < @$xtrafields; $ii ++) {
		if (exists $ARGV{'--self'}) {
			push @outself, $ARGV{'--self'}."_".$xtrafields->[$ii];
		}
		push @outmembers, $cnvflag."_".$xtrafields->[$ii];
	}
}
else {
	for(my $ii = 0; $ii < @$xtrafields; $ii ++) {
		if (exists $ARGV{'--self'}) {
			push @outself, $xtrafields->[$ii]; 
		}
		push @outmembers, $xtrafields->[$ii].".FamMembers";
	}
}

foreach my $field (@outself, @outmembers) {
	if (grep { $field eq $_ } @$infields) {
		unless ($ARGV{'--over-write'}) {
			die "Output field $field already exist in the input file!";
		}

	}
}

# Slurp pedigree info
if ($ARGV{'--xref'} =~ /\.(ped|fam)$/) {
	print STDERR "Sample external reference is a PED file!";
	if ($ARGV{'--ped'}) {
		warn "The PED file provided by --ped will be ignored";
	}
	$ARGV{'--ped'} = $ARGV{'--xref'};
}


my ($famsamp, $famrels, %pedinfo);
if ($ARGV{'--ped'}) {
	($famsamp, $famrels) = fam_rels($ARGV{'--ped'}, 
		{ twins => $ARGV{'--twins'}, shorten => 1, strict => 0,
		  ignore => $ARGV{'--ignore'}, verbose => 0 });
	open my $fin, $ARGV{'--ped'} or die "Cannot open PED file $ARGV{'--ped'}";
	while(<$fin>) {
		my ($famid, $iid, $dad, $mom, $sex, $phe) = split;
		$pedinfo{$iid} = { FamID => $famid,
			Gender => $sex eq '1' ? 'Male' : $sex eq '2' ? 'Female' : 'NA',
			Pheno => $phe eq '1' ? 'Unaffected' : $phe eq '2' ? 'Affected' : 'NA' };
	}
	# FamID/Gender/Pheno for self, FamMembers/Relations/Phenotypes for family members.
	if ($cnvflag) {
		unshift @outmembers, "FamMember_Relations", "FamMember_Phenotypes";
	}
	else {
		unshift @outmembers, "Relations", "Phenotypes";
	}
	if ($ARGV{'--self'}) {
		unshift @outself, "FamID", "Gender", "Pheno";
	}
}

# Determine output fields
my @outfields = @$infields;
foreach my $field (reverse @outself) {
	unless(grep { $field eq $_ } @outfields) {
		insert_after(\@outfields, $idfield, $field);
	}
}
foreach my $field (reverse @outmembers) {
	unless(grep { $field eq $_ } @outfields) {
		if ($cnvflag) {
			insert_after(\@outfields, $cnvflag.'_IDs', $field);
		}
		else {
			insert_after(\@outfields, 'FamMembers', $field);
		}
	}
}
if ($ARGV{'--ped'}) {
	shift @outmembers for 1..2;
	if ($ARGV{'--self'}) {
		shift @outmembers for 1..3;
	}
}

my $fout;
if ($ARGV{'--output'}) {
	open $fout, ">$ARGV{'--output'}" or die "Cannot write $ARGV{'--output'}";
}
else {
	$fout = \*STDOUT;
}
print $fout join("\t", @outfields), "\n";
while(my $dat = $iter->()) {
	# Determine the ID of primary sample and family members
	my $iid = $dat->{$idfield};
	my @members;
	if ($cnvflag) {
		@members = split(',', $dat->{$cnvflag."_IDs"});
	}
	else {
		@members = split(',', $dat->{FamMembers});
	}
	my %known;
	# First update pedigree info if PED is provided
	if ($ARGV{'--ped'}) {
		if ($cnvflag) {
			$dat->{$cnvflag."_Relations"} = join(",", map { $famrels->{$iid}{$_} // "." } @members);
			$dat->{$cnvflag."_Phenotypes"} = join(",", map { $pedinfo{$_}{Pheno} // "." } @members);
			$known{$cnvflag."_Relations"} = 1; $known{$cnvflag."_Phenotypes"} = 1;
		}
		else {
			$dat->{Relations} = join(",", map { $famrels->{$iid}{$_} // "." } @members);
			$dat->{Phenotypes} = join(",", map { $pedinfo{$_}{Pheno} // "." } @members);
			$known{Relations} = 1; $known{Pheontypes} = 1;
		}
		if (exists $ARGV{'--self'}) {
			foreach my $fd (qw|FamID Gender Pheno|) {
				$dat->{$fd} = $pedinfo{$iid}{$fd};	
				$known{$fd} = 1;
			}
		}
	}
	for (my $ii = 0; $ii < @$xtrafields; $ii ++) {
		#$dat->{$outmembers[$ii]} = join(",", map { my $a = $xtrainfo->{$_}{$xtrafields->[$ii]};
		#						if (defined $a && $a !~ /^\s*$/) { $a } else { "." } } @members); 
		$dat->{$outmembers[$ii]} = join(",", map { $xtrainfo->{$_}{$xtrafields->[$ii]} // "." } @members);
		if (@outself) {
			next if defined $known{$outself[$ii]};
			#my $a = $xtrainfo->{$_}{$xtrafields->[$ii]};
			#if (defined $a && $a !~ /^\s*$/) {
			#	$dat->{$outself[$ii]} = $a;
			#}
			#else {
			#	$dat->{$outself[$ii]} = ".";
			#}
			$dat->{$outself[$ii]} = $xtrainfo->{$_}{$xtrafields->[$ii]} // ".";
		}
	}
	print $fout join("\t", map { $dat->{$_} // "." } @outfields), "\n";
}


__END__

=head1 NAME

add_member_info.pl -- Add packed family member info in filtered variant table.

=head1 NOTE

In rare variant filtering pipeline, all family members of the primary sample are packed in 'FamMembers'
as a comma separated list. Their relationships with the primary sample are given in 'Relations', and
phenotypes are in 'Phenotypes'. Other types of member-specific information are packed up in various 
XXX.FamMembers fields. rarevar_burden make use those member-specific information to perform family-based 
filter and tally transmission events to different groups of family members.

This utility adds additional sample-specific information into the table. Note that these information
are only related to sample and not variant-dependent. It supports the CNV table, which requies a prefix 
for family member related fields (default: FamMembers). And family member info are given FamMembers_XXX. 
Information for the primary sample can also be added. It can also be used to update family relationship 
if a PED file is provided. 

Convetions: Input variant should be tab separated. IDs of family members are listed in FamMembers,
Missing data will be ".", member specific info are comma separated. So the individual value 
should not contain comma or it will cause issues in downstream parsing. 

=head1 REQUIRED ARGUMENTS

=over

=item -[-]in[put] [=] <file>

Input variant or CNV table. It can be '-' to read from STDIN.
Input must be tab separated with a header.

=item -[-]xref [=] <file>

Lookup table for external sample info to be added. It can be tab separated text or comma separated csv file.
It can also be a PED file is file name ends with .ped or .fam (must also be tab separated!).

=for Euclid:
	file.type: readable

=item -[-]xref-fields [=] <fields>

Key fields for the external sample info, format: Samp:IID,Info1:Anno1,Info2
The first column will be taken as sample ID and must be remaed to match the sample ID column from input.
The remaining fields will be new annotations to be added to the input.

=back

=head1 OPTIONS

=over

=item -[-]out[put] [=] <file>

Output file name. If it is not provided, will dump the output to STDOUT.

=item -[-]prefix [=] <string>

For CNV table, the prefix for family member related fields.

=for Euclid:
	string.default: "FamMembers"

=item -[-]self [=] [<string>]

Also add the info for the primary him/her-self. 
For CNV table only: An optional string can be provided and used as prefix for field name.

=over -[-]over[-]write

Over-write the content if the field already exist in the input.

=item -[-]ped [=] <pedfile>

Pedigree file. It will be used to determine family relationships.
PED is not needed, if family member related fields exist in the input.

=for Euclid:
	pedfile.type: readable

=item -[-]ignore [=] <string>

Regex pattern for techical duplicate sample IDs in PED file.

=for Euclid:
	string.default: '_Re\d*$'

=item -[-]twins [=] <list>

A list of twin pairs. Should be a file with two columns.

=cut




