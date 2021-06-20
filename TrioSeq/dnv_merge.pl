#!/usr/bin/env perl

use strict;
use warnings;
use Carp;
use IO::Dir;
use IO::File;
use IO::Prompt;
use Data::Dumper;
use FindBin qw|$Bin|;
use List::Util qw|first|;
use List::MoreUtils qw|all any none uniq|;
use Getopt::Lucid qw|:all|; 
use Config::Std;
use Perl6::Slurp;
use Utils::Parser qw|sql_query|;
use Utils::Hash qw|merge_conf chk_default|;
use Utils::File::Iter qw|iter_file|;


use lib "$Bin/../lib";
use Shared qw|read_dir read_list expand_dat parse_tabfile|;

#############################
## Command line interface  ##
#############################

my @spec =  (
	Param("conf|c")->valid(sub { -r }),
	Param("in|i")->valid(sub { -r }),
	Param("out|o"),
	Keypair("param|par"),
	Switch("help|h"),
	);

my $opt = Getopt::Lucid->getopt(\@spec);


if ($opt->get_help) {
	print STDERR <<EOF;
Purpose:
	This is a stand alone script to merge DNV candidates while removing those with contradictory
	evidence from different input sources. 

Usage:
	dnv_merge.pl --conf Config [--in Input] --out Prefix

Options:
	--in: Input DNV table. It overrides the input variant table defined in the config.
	--out: (Required) Output file name prefix.

Notes:
	Merging DNVs is more than taking the union from all callers or across all batches. Multiple input 
	sources also provide information to re-evaluate the evidence of de novo, so that DNVs with contradicting
	contradicting evidence among multiple input sources can be removed.

	To do this, we first merge DNV candidates identified from multiple sources using utils/merge_vars.pl. 
	Then tabulate merged DNV candidates from multiple sources using dnv_table.pl or vartab_pipeline.pl
	to collect site or genotype level information from all sources. Here input sources can be genotypes from 
	different variant callers, the same caller but different genotyping modes, same or overlapping samples 
	from different sequencing platforms or even different genotyping batches. Input sources may or may not be 
	genotype calls that used for initial DNV filtering. The only purpose for them is to provide information
	for re-evaluate the evidence of a DNV if the variant can be found in the source. 

	The process of re-evaluation of evidence are represented by a series of filter expression specified in
	the config file. There are several types of contradicting evidences if a variant can be found, they include 
	but not limited to:
	 * Non-reference parental genotypes
	 * High cohort allele frequencies
	 * High parental alternative allele fraction
	 * Low confidence in parent reference genotype calls
 	Absence of a variant from the input source are not considered as contradictory evidence.

 	The contradictory filters can be customized to different input sources. For example, GATK provide a number of 
 	commonly used QC metrics but not available from other callers. Sites with low mapping quality (MQ) or high 
 	strand bias (SOR) are enriched for false positives and can be considered evidnece against a real variant.

Input/Output:
	The input to the script include merged DNV candidates, and one or more variant tables including site or genotype
	level information for DNVs that can be found from multiple input sources. The merged DNV candidates table should
	include a VarGroup field. It is a comma separated list of group names indicating the groups of input sources
	from which the DNV was initially identified. The config file specifies detailed filtering expressions and results 
	from their combinations to be considered as contradictory evidence. 

	Two files will be in the output: variants that are not contradicted by any input sources will be written to Prefix.txt;
	variants that failed contradictory filter will be written to Prefix.failed.txt. Group names in the VarGroup field from
	the input will be updated in the following way:

	    Group:  Group is DNV's initial call source, and no contradicting evidence was found.
	   +Group:  Group is not the DNV's intial call source, but the DNV can be found in Group's 
	    		variant table without contradicting evidence.
	   -Group:  Group is not DNV's intial call set, the DNV can be found in Group's var table
	   			but contradicting the evidence of DNV. 
	   			Note: It is possible to have some Group in the tabulated output file but not
	   			included in the initial call set. They can be some joint VCF files generated
	   			at later stage used specifically for calibrating genotype qualities.
	   ?Group:  Group is in the DNV's intial call set, but contradicting evidence was found in the 
	   			Group based on the current filter.
	   If a group label does not appear, then the variant cannot be found in that group of input source.

	 For DNVs that failed filter, an additional field will appear in Prefix.failed.txt to flag which filter is failed
	 when evaluating evidence from which Group of input source. 

EOF
	exit 1;
}

$opt->validate({ requires => [qw|conf out|] });

my %conf = merge_conf($opt->get_conf, $opt->get_param); 

if ($opt->get_in) {
	$conf{Input}{File} = $opt->get_in;
}
unless(-f $conf{Input}{File}) {
	die "Cannot find input file: $conf{Input}{File}";
}
unless(defined $conf{Vartab}{Dir} && -d $conf{Vartab}{Dir} || 
	   defined $conf{Vartab}{List} && -f $conf{Vartab}{List}) {
	die "Cannot find variant table directory or list";
}

# First identify group names in vartab directory
#my @groups;
#{
#	foreach my $file (sort grep { /\.txt$/ } IO::Dir->new($conf{Vartab}{Dir})->read()) {
#		$file =~ s/\.txt$//;
#		push @groups, $file;
#	}
#}

my %sources;
if (defined $conf{Vartab}{Dir}) {
	my $sources = read_dir($conf{Vartab}{Dir}, { suffix => "txt" });
	%sources = %$sources;
}
elsif (defined $conf{Vartab}{List}) {
	my $sources = read_list($conf{Vartab}{List}, { suffix => "txt", multi => 1 });
	%sources = %$sources;
}
else {
	die "Cannot find variant table directory or list";
}


my @groups = sort keys %sources;
if (@groups > 0) {
	print STDERR "The following var table groups were found: \n";
	print STDERR join(", ", @groups), "\n";
	# Check if any custom field or filters exist
	foreach my $field (qw|VarFields GenoFields Filter FilterName|) {
		my @customgrps = map { s/^${field}_//; $_ } grep { /^${field}_\w+/ } keys %{$conf{Vartab}};
		unless (all { defined $sources{$_} } @customgrps) {
			die "Not all custom groups for $field can be found in input sources";
		}
	}
}
else {
	print STDERR "No var table group was found\n";
	exit 1;
}


# Go through each de novo candidates in the vartab output
# keep track of their appearance and test if they meet the criteria
# filtflag is the note label for variants failed filter
# varflag indicate if only variant level information will be used for filter
# strict flag indicate if geno fields will can be filled in when their length differ
my (%filtflag, $varflag, $strictflag);
foreach my $group (@groups) {
	my (@filterstrs, @filternames);
	# Adding general filter applicable to all groups
	if (exists $conf{Vartab}{Filter}) {
		print STDERR "Parsing standard filters\n";
		if (ref $conf{Vartab}{Filter} eq 'ARRAY') {
			push @filterstrs => @{$conf{Vartab}{Filter}};
			if (exists $conf{Vartab}{FilterName}) {
				unless( ref $conf{Vartab}{FilterName} eq 'ARRAY' &&
					 @{$conf{Vartab}{Filter}} == @{$conf{Vartab}{FilterName}} ) {
					die "Numbers of filters and filter names do not match";
				}
				push @filternames => @{$conf{Vartab}{FilterName}};
			}
		}
		else {
			push @filterstrs => $conf{Vartab}{Filter};
			if (exists $conf{Vartab}{FilterName}) {
				if (ref $conf{Vartab}{FilterName} eq 'ARRAY') {
					die "Numbers of filters and filter names do not match";
				}
				push @filternames => $conf{Vartab}{FilterName};
			}
		}
	}
	if (exists $conf{Vartab}{"Filter_$group"}) {
		print STDERR "Parsing custom filter for $group\n";
		if (ref $conf{Vartab}{"Filter_$group"} eq 'ARRAY') {
			push @filterstrs => @{$conf{Vartab}{"Filter_$group"}};
			if (exists $conf{Vartab}{"FilterName_$group"}) {
				unless( ref $conf{Vartab}{"FilterName_$group"} eq 'ARRAY' &&
					 @{$conf{Vartab}{"Filter_$group"}} == @{$conf{Vartab}{"FilterName_$group"}} ) {
					die "Number of filters and filter names does not match for $group";
				}
				push @filternames => @{$conf{Vartab}{"FilterName_$group"}};
			}	
		}
		else {
			push @filterstrs => $conf{Vartab}{"Filter_$group"};
			if (exists $conf{Vartab}{"FilterName_$group"}) {
				if (ref $conf{Vartab}{"FilterName_$group"} eq 'ARRAY') {
					die "Number of filters and filter names does not match for $group";
				}
				push @filternames => $conf{Vartab}{"FilterName_$group"};
			}
		}
	}
	

	# Prepare callbacks 
	# Fields from all filters will be extracted
	my (@callbacks, @labels, @fields);		
	for(my $ii = 0; $ii < @filterstrs; $ii ++) {
		my $filterstr = $filterstrs[$ii];
		$filterstr =~ s/^["']//; $filterstr =~ s/["']$//; 
		if (@filternames) {
			print $filternames[$ii], ": ", $filterstr, "\n";
		}
		else {
			print $filterstr, "\n";
		}
		my ($callback, $tokens) = sql_query($filterstr, 1);
		push @callbacks, $callback;
		foreach my $tok (@$tokens) {
			push @fields, $tok->[1] if $tok->[0] eq 'FIELD';
		}
		if (@filternames) {
			push @labels, $filternames[$ii];
		}
	}

	if (@callbacks > 1) {
		unless (defined $conf{Vartab}{MultiAction}) {
			die "Must provide MultiAction for multiple filters";
		}
	}
	
	# Parse variant and genotype fields
	# When multiple input files are associated with a group, all of them should be checked
	my (@vfields, @gfields);
	if (exists $conf{Vartab}{"VarFields_$group"}) {
		@vfields = split(',', $conf{Vartab}{"VarFields_$group"});
	}
	else {
		@vfields = split(',', $conf{Vartab}{VarFields});
	}
	if (@vfields == 4) {
		print STDERR "Filter for $group is based on variant information only, ignoring individual genotypes\n";
		$varflag = 1;
	}
	else {
		die "Must provide IID,Chrom,Pos,Ref,Alt fields for Vartab" unless @vfields == 5;
		$strictflag = 1;
		if (exists $conf{Vartab}{GenoFields}) {
			push @gfields => split(',', $conf{Vartab}{GenoFields});
			if (exists $conf{Vartab}{GenoFillin} && $conf{Vartab}{GenoFillin} =~ /^Y/i) {
				$strictflag = 0;
			}
		}
		if (exists $conf{Vartab}{"GenoFields_$group"}) {
			push @gfields => split(',', $conf{Vartab}{"GenoFields_$group"});
			if (exists $conf{Vartab}{"GenoFillin_$group"} && $conf{Vartab}{"GenoFillin_$group"} =~ /^Y/i) {
				$strictflag = 0;
			}
		}
		
	}


	my @tabfiles;
	if (ref $sources{$group} eq 'ARRAY') {
		@tabfiles = @{$sources{$group}};
	}
	else {
		@tabfiles = ($sources{$group});
	}

	foreach my $tabfile (@tabfiles) {
		my ($it, $fnames) = iter_file($tabfile, {fsep => qr/\t/});
		foreach my $field (@vfields, @gfields, @fields) {
			unless(grep { $_ eq $field } @$fnames) {
				die "Cannot find field $field in the tabulated variant file: $tabfile";
			}
		}
		while(my $dat = $it->()) {
			my $viid = join(":", @{$dat}{@vfields});
			my @vdat;
			if (@gfields > 0) {
				@vdat = expand_dat($dat, {optional => \@gfields, strict => $strictflag, sep => ','});
			}
			else {
				@vdat = ($dat);
			}
			for(my $ii = 0; $ii < @callbacks; $ii ++) {
				my $callback = $callbacks[$ii];
				my $label;
				if (@labels > 0) {
					$label = $labels[$ii];
				}
				else {
					$label = "Filter".($ii+1);
				}
				foreach my $vdat (@vdat) {
					my $flag;
					if ($callback->($vdat)) {
						$flag = 'Pass';
					}
					else {
						$flag = 'Fail';
					}
					push @{$filtflag{$viid}{$group}{$label}} => $flag;
				}
			}
		}
	}
}

# Consolidate filter flags
my %dnv;
foreach my $viid (keys %filtflag) {
	foreach my $group (keys %{$filtflag{$viid}}) {
		foreach my $label (keys %{$filtflag{$viid}{$group}}) {
			if (lc($conf{Vartab}{DupAction}) eq 'any') {
				if (any { $_ eq 'Pass' } @{$filtflag{$viid}{$group}{$label}}) {
					$dnv{$viid}{$group}{$label} = 'Pass';
				}
				else {
					$dnv{$viid}{$group}{$label} = 'Fail';
				}
			}
			elsif (lc($conf{Vartab}{DupAction}) eq 'all') {
				if (all { $_ eq 'Pass' } @{$filtflag{$viid}{$group}{$label}}) {
					$dnv{$viid}{$group}{$label} = 'Pass';
				}
				else {
					$dnv{$viid}{$group}{$label} = 'Fail';
				}
			}
			elsif (lc($conf{Vartab}{DupAction}) eq 'ignore' || 
				   lc($conf{Vartab}{DupAction}) eq 'first') {
				if ($filtflag{$viid}{$group}{$label}[0] eq 'Pass') {
					$dnv{$viid}{$group}{$label} = 'Pass';
				}
				else {
					$dnv{$viid}{$group}{$label} = 'Fail';
				}
			}
			elsif (lc($conf{Vartab}{DupAction}) eq 'last') {
				if ($filtflag{$viid}{$group}{$label}[-1] eq 'Pass') {
					$dnv{$viid}{$group}{$label} = 'Pass';
				}
				else {
					$dnv{$viid}{$group}{$label} = 'Fail';
				}
			}
			else {
				die "Cannot recognize DupAction: $conf{Vartab}{DupAction}";
			}
		}
	}
}




# Writing results
my $output = $opt->get_out; $output =~ s/\.txt$//;

my ($it, $fnames, $keyfields) = parse_tabfile($conf{Input}{File}, 
		$conf{Input}{GroupFields}.",".$conf{Input}{VarFileds}, 6, 7);

my @keyfields = @$keyfields;
shift @keyfields if @keyfields == 7;


open my $fout, ">$output.txt" or die "Cannot write to $output.txt";
print $fout join("\t", @$fnames), "\n";
open my $ferr, ">$output.failed.txt" or die "Cannot write to failed variants $output.failed.txt";
if (exists $conf{Output}{ErrField}) {
	print $ferr join("\t", $conf{Output}{ErrField}, @$fnames), "\n";
}
else {
	print $ferr join("\t", @$fnames), "\n";
}

while(my $dat = $it->()) {
	my ($vargrp, $iid, $chrom, $pos, $ref, $alt) = @{$dat}{@keyfields};
	# Note: sampgrps are not used in the current version
	# VarGroup is the initial group in which variant was called
	# Need to strip off '+/-/?' prefix if the variant table was from previous dnv_merge output.
	# Although '-/?' prefix should not exist, because those variants were already filtered out.
	# They can exist if we are looking into dropped variant. So we also support special cases.
	# The following rules will apply:
	# 1. '+' Groups will be treated as a new support group. If filter support the variant, the prefix 
	# will be preserved; otherwise, prefix will change to '-'.
	# 2. '?/-' Groups will be treated as non-support group. Even if current filter support the variant
	# They will still be considered as failed group, with no changes made to the existing prefix.
	# But they will not be output to failed.txt.
	# Only if current filter contradict the variant, it will be written to failed.txt.
	my (@vargrps, %prefix, %reason);
	foreach my $group (split(',', $vargrp)) {
		if ($group =~ /^([\+\-\?])/) {
			my $pref = $1;
			$group =~ s/^[\+\-\?]//;
			push @vargrps, $group;
			$prefix{$group} = $pref;
		}
		else {
			push @vargrps, $group;
			$prefix{$group} = "";
		}	
	}
	# my @vargrps = map { $_ =~ s/^[\+\-\?]//; $_ } split(',', $vargrp);

	my $viid;
	if ($varflag) {
		$viid = join(":", $chrom, $pos, $ref, $alt);
	}
	else {
		$viid = join(":", $iid, $chrom, $pos, $ref, $alt);
	}

	my @failflag;
	if ($dnv{$viid}) {
		foreach my $group (sort keys %{$dnv{$viid}}) {
			unless(grep { $_ eq $group } @vargrps) {
				push @vargrps, $group;
			}

			my $fail;
			if (defined $conf{Vartab}{MultiAction}) {
				if ($conf{Vartab}{MultiAction} eq 'Any') {
					unless (any { $_ eq 'Pass' } values %{$dnv{$viid}{$group}}) {
						$fail = 1;
					}
				}
				elsif ($conf{Vartab}{MultiAction} eq 'All') {
					unless (all { $_ eq 'Pass' } values %{$dnv{$viid}{$group}}) {
						$fail = 1;
					}
				}
				#elsif ($conf{Vartab}{MultiAction} eq 'None') {
				#	unless (none { $_ eq 'Pass' } values %{$dnv{$viid}{$group}}) {
				#		$fail = 1;
				#	}
				#}
				else {
					die "Cannot recognize Vartab.MultiAction: $conf{Vartab}{MultiAction}";
				}
			}
			else {
				my @filts = keys %{$dnv{$viid}{$group}};
				unless (@filts == 1) {
					die "When MultiAction is not defined, should contain one and only one filter!";
				}
				unless ($dnv{$viid}{$group}{$filts[0]} eq 'Pass') {
					$fail = 1;
				}
			}

			# if Failed, need to change existing prefix, listing fail reasons
			if ($fail) {
				push @failflag, $fail;
				if (defined $prefix{$group}) {
					if ($prefix{$group} eq "") {
						$prefix{$group} = "?";
					}
					elsif ($prefix{$group} eq "+")  {
						$prefix{$group} = "-";
					}
				}
				else {
					$prefix{$group} = '-';
				}
				if (exists $conf{Output}{ErrField}) {
					$reason{$group} = join(',', grep { $dnv{$viid}{$group}{$_} ne 'Pass' } sort keys %{$dnv{$viid}{$group}});
				}
			}
			else {
				# If success, adding support if the group was not listed before
				unless(defined $prefix{$group}) {
					$prefix{$group} = "+";
				}
				# Otherwise, do not make changes to existing prefix
			}
		}
	}

	$dat->{$keyfields[0]} = join(',', map { $prefix{$_}.$_ } @vargrps);
	if (@failflag) {
		if (exists $conf{Output}{ErrField}) {
			my $note = join(";", map { $_.":".$reason{$_} } sort keys %reason);
			print $ferr join("\t", $note, @{$dat}{@$fnames}), "\n";
		}
		else {
			print $ferr join("\t", @{$dat}{@$fnames}), "\n";
		}
	}
	else {
		print $fout join("\t", @{$dat}{@$fnames}), "\n";
	}
}


