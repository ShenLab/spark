#!/usr/bin/env perl

use strict;
use warnings;
use List::MoreUtils qw|any all uniq|;
use Sort::Versions;
use Perl6::Slurp;
use FindBin qw|$Bin|;
use Getopt::Euclid;
use Genet::Ped;
use Genet::File::VCF;
use Utils::Stat qw|mean|;
use Utils::File::Iter qw|iter_file|;

use lib "$Bin/../../lib"; 
use Shared qw|merge_cols|;

# List file format: *VCF, OUTPUT, LABEL*.

# Find all samples in each group
my (%samps, %priority, %cols);
{
	my $pattern = qr/$ARGV{'--ignore'}/;
	my $ped = Genet::Ped->new($ARGV{'--ped'}, { ignore => $pattern, verbose => 0 });

	my @groups;
	open my $fin, $ARGV{'--list'} or die "Cannot read list $ARGV{'--list'}";
	while(<$fin>) {
		my @dat = split;
		unless (@dat == 3) {
			die "Incorrect number of columns in the list: $_";
		}
		my ($vcffile, $vartab, $label) = @dat;
		push @groups, $label;
				
		# Slurp in column names
		my ($it, $fnames) = iter_file($vartab, { fsep => qr/\t/ });
		if (defined $cols{$label}) {
			unless(join("\t", @{$cols{$label}}) eq join("\t", @$fnames)) {
				die "DNV table $vartab has different columns from other tables of the same group";
			}
		}
		else {
			$cols{$label} = $fnames;
		}

		my %vcfsamp;
		if ($vcffile =~ /\.vcf$/ || $vcffile =~ /\.vcf\.gz$/) {
			my $vcf = Genet::File::VCF->new($vcffile);
			%vcfsamp = map { $_ => 1 } $vcf->get_sampids();
		}
		else {
			%vcfsamp = map { (split)[0] => 1 } slurp $vcffile;
		}
		
		if ($ARGV{'--trio'}) {
			foreach my $famid ($ped->get_famids()) {
				foreach my $trio ($ped->get_trios($famid)) {
					if (all { defined $vcfsamp{$_} } @$trio) {
						$samps{$trio->[0]}{$label} = 1;
					}
				}
			}
		}
		else {
			foreach my $famid ($ped->get_famids()) {
				foreach my $iid ($ped->get_members($famid)) {
					if (defined $vcfsamp{$iid}) {
						$samps{$iid}{$label} = 1;
					}
				}
			}
		}	
	}
	# Determine priority list for each group based on the order of input list
	for(my $ii = 0; $ii < @groups; $ii ++) {
		next if defined $priority{$groups[$ii]};
		$priority{$groups[$ii]} = $ii;
	}
}

# Merge columns
my @mergecols = merge_cols(\%cols, \%priority);
print STDERR "Merged columns:\n";
print STDERR join(" ", @mergecols), "\n";

my @newcols = split(',', $ARGV{'--fields-add'}); 
unless(@newcols == 2) {
	die "Must provide names for two new columns corrsponding to sample and variant groups";
}
else {
	foreach my $field (@newcols) {
		if(grep { $field eq $_ } @mergecols) {
			die "New column name $field already exist in the merged input!";
		}
	}
}

# Go through each DNV table and slurp data from file
# Variants will be represented by Chrom,Pos,Ref,Alt
# {VarID}{Group} => VarData
# If the same variant appear more than once in the same group
# only the first one will be stored.
my (%known, %vars, %chrpos);
{
	my @varcols = split(',', $ARGV{'--fields'});
	unless (@varcols == 5) {
		die "Incorrect number of fields for IID,Chrom,Pos,Ref,Alt: $ARGV{'--fields'}";
	}
	open my $fin, $ARGV{'--list'} or die "Cannot read list $ARGV{'--list'}";
	while(<$fin>) {
		my ($vcffile, $vartab, $label) = split;
		next if $known{$vartab};
		my ($it, $fnames) = iter_file($vartab, { fsep => qr/\t/ });
		$known{$vartab} = 1;
		while(my $dat = $it->()) {
			unless(all { defined $dat->{$_} } @varcols) {
				die "Cannot find IID,Chrom,Pos,Ref,Alt for a variant in $vartab";
			}
			my $sampid = $dat->{$varcols[0]};
			# Skip samples that does not appear in PED file
			unless(defined $samps{$sampid}) {
				warn "Cannot find sample $sampid in PED file";
				next;
			}
			my $sampvar = join("\t", @{$dat}{@varcols});
			# Skip known sample-var combination in the same group 
			if (exists $vars{$sampvar}{$label}) {
				warn "Variants $sampvar has been observed in group $label";
			}
			else {
				$vars{$sampvar}{$label} = $dat;
			}
			if (!defined $chrpos{$sampvar}) {
				$chrpos{$sampvar} = [ @{$dat}{@varcols[0..2]} ];
			}
		}
	}
}

# For each variants, pick the data from the highest priority group
my $fout;
if (defined $ARGV{'--output'}) {
	open $fout, ">$ARGV{'--output'}" or die "Cannot write to output $ARGV{'--output'}";
}
else {
	$fout = \*STDOUT;
}

print $fout join("\t", @newcols, @mergecols), "\n";
foreach my $sampvar (sort { versioncmp($chrpos{$a}[1],$chrpos{$b}[1]) ||
						  	$chrpos{$a}[2] <=> $chrpos{$b}[2] ||
						  	$chrpos{$a}[0] cmp $chrpos{$b}[0]
						  	 } keys %vars) {
	my ($iid, $chr, $pos, $ref, $alt) = split(/\t/, $sampvar);

	my @groups = sort { $priority{$a} <=> $priority{$b} } keys %{$vars{$sampvar}};
	my $varGrp = join(',', @groups);
	my $sampGrp = join(',', sort { $priority{$a} <=> $priority{$b} } keys %{$samps{$iid}});
 
	my $vardat = $vars{$sampvar}{$groups[0]};
	print $fout join("\t", $sampGrp, $varGrp, map { $_ // "." } @{$vardat}{@mergecols}), "\n";
}


__END__

=head1 NAME

merge_vars.pl -- Merge variants from different batches or callers.

=head1 NOTES

This script merge variants by taking the union of both variants and table columns.

We assume variant table files are all tab separated and the files from the same group 
should share the same columns.

Duplicated entries will be automatically removed to keep only the first appearance.

The script is only intended to merge small number of variants.

=head1 REQUIRED ARGUMENTS
 
=over
 
=item -[-]ped [=] <file>   
 
The PED file containing all samples in the variant output.
Only samples present in both PED and VCF will be used.
It will also be used to find trios when --trio option is enabled.

=for Euclid:
	file.type: readable
 
=item -[-]list [=] <file>

A list of variant table files and associated original VCFs or sample list.
The list file should contain following three columns: *VCF VARTAB LABEL*. 
VCF is only used to determine if a sample or trio is present in the VCF, it can also be replaced by a sample list file.
LABEL is the sample group label associated with the VARTAB.
When a variant can be found in multiple variant table, the information from the first table that appears in this list
will be used in output.

=for Euclid:
	file.type: readable

=back
 
=head1 OPTIONS
 
=over

=item -[-]out[put] [=] <filename>

The output file name include the union of variants from all input files.
If not provided, we will dump the results to stdout.

=item -[-]ignore [=] <pattern>

Pattern for duplicated samples, they will be ignored when parsing PED file.

=for Euclid:
	pattern.default: '_Re\d*$'

=item -[-]trio

Switch on the trio mode, require the presence of a complete trio for a sample to be labeled in each group.

=item -[-]fields-add [=] <string>

Two additional columns will be added to the combined output indicating the groups that
sample appears and variant was called. Default: SampGroup,VarGroup
The column names should not exist in the input.

=for Euclid:
	string.default: "SampGroup,VarGroup"

=item -[-]fields [=] <colnames>

Column names for IID,chrom,pos,ref,and alt alleles. These columns are presummed have the 
same name across all input files. We also assume variants has already been normalized. 
So the data from the last four columns will be used in uniquely represent a variant. 
Chromosome and position will also be used in sorting the final output.

=for Euclid:
	colnames.default: "IID,Chrom,Position,Ref,Alt"

=back
