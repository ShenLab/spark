use strict;
use warnings;
use FindBin qw|$Bin|;
use Utils::List qw|insert_after|;
use Utils::File::Iter qw|iter_file slurp_file|;
use Getopt::Euclid;

use lib "$Bin/../../lib";
use Shared qw|fam_rels reltype slurp_xref|;


# Input and output
my ($prefix, $output);
if ($ARGV{'--input'} =~ /\.kin0?$/) {
	$prefix = $ARGV{'--input'}; $prefix =~ s/\.kin0?//;
}
else {
	$prefix = $ARGV{'--input'};
}
unless(-f "$prefix.kin") {
	die "Cannot find input file $prefix.kin!";
}

if (defined $ARGV{'--output'}) {
	$output = $ARGV{'--output'};
}
else {
	$output = $prefix;
}

my $pattern = qr/$ARGV{'--ignore'}/;

# Collect samp/ped info from PED
my (%sampinfo, %famsamp);
{
	open my $fin, $ARGV{'--ped'} or die "Cannot open $ARGV{'--ped'}";
	while(<$fin>) {
		my ($fid, $iid, $dad, $mom, $sex, $phe) = split;
		my $bipar;
		if ($dad ne '0' && $mom ne '0') {
			$bipar = "b";
		}
		elsif ($dad eq '0' && $mom ne '0') {
			$bipar = "m";
		}
		elsif ($dad ne '0' && $mom eq '0') {
			$bipar = "p";
		}
		else {
			$bipar = ".";
		}
		$sampinfo{$iid} = { FID => $fid, Bipar => $bipar,
			Sex => $sex eq '1' ? 'Male' : $sex eq '2' ? 'Female' : '.',
			Pheno => $phe eq '1' ? 'Unaffected' : $phe eq '2' ? 'Affected' : '.' };
		next if $iid =~ $pattern;
		push @{$famsamp{$fid}}, $iid;
	}
}

# Collect pairwise ped relationship using fam_rels 
# Expected degree of relatedness will be inferred based on familial relationship
# Currently only up to 2nd degree
my ($famsamp, $famrels) = fam_rels($ARGV{'--ped'},
	{ twins => $ARGV{'--twins'}, ignore => $ARGV{'--ignore'}, shorten => 1, fullsib => 1, 
	  strict => 0, verbose => 0 });

# Collect external samp info
my ($xsampinfo, $xsampfields);
if (defined $ARGV{'--sxref'}) {
	unless(defined $ARGV{'--sxref-fields'}) {
		die "Must define fields in external sample reference file!";
	}
	($xsampinfo, $xsampfields) = slurp_xref($ARGV{'--sxref'}, $ARGV{'--sxref-fields'});
}

my ($it, $fnames) = iter_file("$prefix.kin");
my @outfields = qw|FID ID1 ID2 IBS0 Kinship IBD2Seg IBD1Seg PedRel PedType InfType|;
foreach my $field (@outfields) {
	unless(grep { $_ eq $field } @$fnames) {
		next if $field =~ /^Ped/;
		die "Cannot find $field in $prefix.kin!";
	}
	if (grep { $_ eq $field } @$xsampfields) {
		die "Field $field in external sample ref already exists in $prefix.kin";
	}
}
insert_after(\@outfields, "ID1", "Sex1", "Pheno1", "Bipar1", 
		defined $xsampfields ? map { "${_}1" } @$xsampfields : ());
insert_after(\@outfields, "ID2", "Sex2", "Pheno2", "Bipar2", 
		defined $xsampfields ? map { "${_}2" } @$xsampfields : ());

my ($it2, $fnames2);
my @outfields2 = qw|FID1 ID1 FID2 ID2 IBS0 Kinship IBD2Seg IBD1Seg InfType|;
my $kin0flag;
if (-f "$prefix.kin0") {
	($it2, $fnames2) = iter_file("$prefix.kin0");
	foreach my $field (@outfields2) {
		unless(grep { $_ eq $field } @$fnames2) {
			die "Cannot find $field in $prefix.kin0!";
		}
		if (grep { $_ eq $field } @$xsampfields) {
			die "Field $field in external sample ref already exists in $prefix.kin0";
		}
	}
}
insert_after(\@outfields2, "FID1", "NSamp1");
insert_after(\@outfields2, "ID1", "Sex1", "Pheno1", "Bipar1", 
		defined $xsampfields ? map { "${_}1" } @$xsampfields : ());
insert_after(\@outfields2, "FID2", "NSamp2");
insert_after(\@outfields2, "ID2", "Sex2", "Pheno2", "Bipar2", 
		defined $xsampfields ? map { "${_}2" } @$xsampfields : ());


open my $fout, ">$output.kin.txt" or die "Cannot write to $output.kin.txt";
print $fout join("\t", @outfields), "\n";
open my $fout2, ">$output.kin0.txt" or die "Cannot write to $output.kin0.txt";
print $fout2 join("\t", @outfields2), "\n";

while(my $dat = $it->()) {
	for my $ii (1,2) {
		my $iid = $dat->{"ID$ii"};
		if (defined $sampinfo{$iid}) {
			foreach my $fd (qw|FID Sex Pheno Bipar|) {
				$dat->{"${fd}${ii}"} = $sampinfo{$iid}{$fd};
			}
		}
		else {
			warn "Cannot find sample $iid in PED file!";
		}
		if (defined $xsampinfo && defined $xsampinfo->{$iid}) {
			foreach my $fd (@$xsampfields) {
				my $val = $xsampinfo->{$iid}{$fd};
				if (!defined $val || $val =~ /^s*$/) {
					$val = ".";
				}
				$dat->{"${fd}${ii}"} = $val;
			}
		}
	}
	if (defined $dat->{FID1} && defined $dat->{FID2} && $dat->{FID1} ne $dat->{FID2}) {
		$kin0flag = 1;
		foreach my $ii (1,2) {
			my $fid = $dat->{"FID$ii"};
			if (defined $famsamp{$fid}) {
				$dat->{"NSamp$ii"} = scalar(@{$famsamp{$fid}});
			}
		}
		print $fout2 join("\t", map { $dat->{$_} // "." } @outfields2), "\n";
	}
	else {
		$dat->{PedRel} = $famrels->{$dat->{ID1}}{$dat->{ID2}} // $famrels->{$dat->{ID2}}{$dat->{ID1}};
		if (defined $dat->{PedRel}) {
			$dat->{PedType} = reltype($dat->{PedRel});
		}
		print $fout join("\t", map { $dat->{$_} // "." } @outfields), "\n";
	}
}

if (-f "$prefix.kin0") {
	while(my $dat = $it2->()) {
		for my $ii (1,2) {
			my $iid = $dat->{"ID$ii"};
			if (defined $sampinfo{$iid}) {
				foreach my $fd (qw|FID Sex Pheno Bipar|) {
					$dat->{"${fd}${ii}"} = $sampinfo{$iid}{$fd};
				}
			}
			else {
				warn "Cannot find sample $iid in PED file!";
			}
			if (defined $xsampinfo && defined $xsampinfo->{$iid}) {
				foreach my $fd (@$xsampfields) {
					my $val = $xsampinfo->{$iid}{$fd};
					if (!defined $val || $val =~ /^s*$/) {
						$val = ".";
					}
					$dat->{"${fd}${ii}"} = $val;
				}
			}
		}
		if (defined $dat->{FID1} && defined $dat->{FID2} && $dat->{FID1} ne $dat->{FID2}) {
			$kin0flag = 1;
			foreach my $ii (1,2) {
				my $fid = $dat->{"FID$ii"};
				if (defined $famsamp{$fid}) {
					$dat->{"NSamp$ii"} = scalar(@{$famsamp{$fid}});
				}
			}
			print $fout2 join("\t", map { $dat->{$_} // "." } @outfields2), "\n";
		}
		else {
			$dat->{PedRel} = $famrels->{$dat->{ID1}}{$dat->{ID2}} // $famrels->{$dat->{ID2}}{$dat->{ID1}};
			if (defined $dat->{PedRel}) {
				$dat->{PedType} = reltype($dat->{PedRel});
			}
			print $fout join("\t", map { $dat->{$_} // "." } @outfields), "\n";
		}
	}
}


__END__

=head1 NAME

anno_king_out.pl -- Annotate the output from king.

=head1 NOTE

The script requires a PED file to get basic sample info and to determine the pairwise 
pedigree relationships between samples. The sample level info like parents, sex, phenotype 
will be added after ID1 and ID2. And known pedigree relationships will be annotated.

King results can be found in two files: .kin for within-family relationships, kin0 for between
family cryptic relatedness. The output will also have two files .kin.txt and .kin0.txt for
reformatted within and between family relationships respecitvely. Inidividdal's family ID
may be different between king's output and PED file, so whether the realtionship is within or 
between will be determined by the specified PED file.
 
=head1 REQUIRED ARGUMENTS

=over 

=item -[-]in[put] [=] <prefix>

Input file or input file prefix to the output of king.
Output will be written to prefix.kin.txt and prefix.kin0.txt unless otherwise specified.

=item -[-]ped [=] <pedfile>

Pedigree file.

=for Euclid:
	pedfile.type: readable

=back

=head1 OPTIONS

=over

=item -[-]twins [=] <list>

A list of twin pairs. Should be a file with two columns.

=for Euclid:
	list.type: readable

=item -[-]ignore [=] <string>

Regex pattern for duplicated sample IDs (default: '_Re\d*$').

=for Euclid:
	string.default: '_Re\d*$'

=item -[-]sxref [=] <file>...

External reference file for additional sample level information.
Must be tab separated text file or comma separated csv format.
Can be more than one file

=item -[-]sxref-fields [=] <fields>...

Key fields to in the external sample info table, starting with sample ID followed by
the fields to be extracted. Fields can be renamed to avoild conflict with existing fields.
Example: ID,Age,Diagnosis:Dx.

=item -[-]output [=] <prefix>

Output file prefix, default is the same as input.

=back

