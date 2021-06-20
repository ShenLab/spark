#!/usr/bin/env perl

use strict;
use warnings;
use Carp;
use IO::Dir;
use IO::File;
use Perl6::Slurp;
use Data::Dumper;
use File::Temp qw|tempdir|;
use File::Path qw|make_path|;
use FindBin qw|$Bin|;
use Hash::Util qw|lock_hash_recurse|;
use List::Util qw|min max first|;
use List::MoreUtils qw|all any uniq none|;
use Getopt::Lucid qw|:all|; 
use Config::Std;
use Statistics::R;
use Utils::List qw|all_pairs insert_after|;
use Utils::Hash qw|chk_default|;
use Utils::Parser qw|sql_query|;
use Utils::File::Iter qw|iter_file|;
use Genome::Ranges::IntSet;
use Genome::UCSC::GeneTable qw|iter_geneTab iter_genePred gdat_list_coding_exons|;


use lib "$Bin/../lib";
use Shared qw|parse_filters parse_tabfile read_geneset parse_fstr slurp_xref expand_dat read_twinsibs|;
use Variants qw|expand_fam var_type|;


my @spec =  (
	Param("conf")->valid(sub { -r }),
	Param("input|in|i")->valid(sub { -r }),
	Param("wrkdir|wrk"),
	Param("output|out|o"),
	Switch("force"), # Under force mode, previously generated files will not be over-written
	Switch("help|h")
	);

my $opt = Getopt::Lucid->getopt(\@spec);


if ($opt->get_help) {
	print STDERR "rg_caf.pl --conf CONFIG -in INPUT -out OUTPUT\n";
	print STDERR <<EOF;
Purpose:
	This is a standalone script to evaluate and compare the burden of recessive genotypes.

Usage:
	rarevar_burden.pl [--in VarTable --wrk WorkDir] --conf Config --out Prefix

Options:
	--in: Input variant table. It overrides the variant table defined in the config.
	--wrk: Working directory. Intermediate filter-sorted variants are stored in this directory. If not provided
 		   a temp directory will be used. If working directory already exists and contain intermediate table
		   from previous run, it will be used directly to skip the time-consuming filtering and sorting steps. 
	--out: Output file name prefix. If not provided, it will use the prefix of config file.
	
Notes:
	The recessive genotype (RG) for a variant class are defined as either homozygous genotypes (hom) or compound heterozygotes
	(cHet) in a gene. To determine cHet, we currently use familial transmission to phase multiple heterozygotes in the
	same gene. Additional information like linked reads or distant relatives can be exploited by has not been implemented.
	So for example if an individual carries both a de novo and inherited heterzygous variants in the same gene, use familial 
	transmission alone cannot infer if they comprise a cHet. 

	The burden of RGs is summarized by the observed number of hom or cHet RGs of a variant class in a gene set from a group
	of samples. If multiple related individuals (MZ twins or sibs) exist in the same group, only RGs from the first one will 
	be kept. Under random mating, the expected number of RGs is sample size times the squared CAF; under recent common ancestry 
	(consangunity), it is sample size times CAF. To account for recent common ancestry, long runs of homozygosity (a.k.a. 
	homozygosity by descent (HBD) segments) can be inferred for each individual and incoporated in computing expected RGs in
	a cohort. We usually calculated CAF for different variant classes in each gene in parental genotypes that are filtered 
	in the same way as offspring used in evaluating the RG burden. When CAF is not available, burden can also be evaluated by 
	comparing RG rates between different sample groups. One-sided Poisson test is used in these comparisons. Currently, the 
	RG analysis is gene-centric. Variants that cannot be annotated to genes are not included.  

Output:
	The input files to the RG burden analysis are defined in the config. They should include: one or more variant tables,
	a gene-specific CAF table different classes of variants (genereated by rg_burden.pl), and one or more sample tables including
	phenotype, sex and other sample level information. Variant classes, sample groups, and gene sets are defined in the config file.

	The following files will appear in the output:

 	 * prefix.burden.txt -- Summary of the observed RG burden w.r.t expectations predicted by CAF of different variant class 
 	 						in different sample groups and gene sets.
	 * prefix.contra.txt -- Comparison of RG burdens between different sample groups and by gene set and variant class 
	 						if contrasting groups are defined.
	 * prefix.genetab.txt -- Table of RG counts and CAF for each gene for selected sample groups and variant classes 
	 * prefix.samptab.txt -- Table of RG counts in different gene sets and variant classes per sample.
	 						 This table includes all samples passed filter and inclusion/exclusion list. 
	 * prefix.RGs.txt -- list of RGs on the selected sample group, variant class and gene set.

	 Filter-sorted variants will be stored in an intermediate file vartab.txt.gz in the working directory. 

EOF
	exit 1;
}

$opt->validate({ requires => [qw|conf|] });


read_config $opt->get_conf => my %conf;
# If input is provided, we will over-ride the VarTab in config
if ($opt->get_input) {
	$conf{Variant}{Table} = $opt->get_input;
	unless(-f $conf{Variant}{Table}) {
		die "Cannot find variant table: $conf{Variant}{Table}";
	}
}
lock_hash_recurse(%conf);


my $output;
if($opt->get_output) {
	$output = $opt->get_output;
}
else {
	$output = $opt->get_conf;
	$output =~ s/\.conf$//;
}
if ($output =~ /\.txt$/) {
	$output =~ s/\.txt$//;
}


my $wrkdir = $opt->get_wrkdir;
if ($wrkdir) {
	mkdir $wrkdir unless -d $wrkdir;
} else {
	$wrkdir = tempdir(CLEANUP => 1);
}
make_path "$wrkdir/tmp";

#
# Sample Information
#
my (@sampgroup, %sampgroup, %sampsize, %sampincl, %samprm, %samprmbygrp, %famrm, @gpairs);
# @sampgroup is an array of sample groups defined in the config
# %sampgroup is a hash of hash to hold group membership of each sample $sampgroup{IID}{Group} = 1
# %sampsize store total sample size for each sample group $sampsize{Group}=Count
# -- note: we only analyze autosomes for RG, so sample sex is not modeled
# %sampincl/%samprm: hash of sample IDs to be included or removed
# If sample table is not provided, only one group "All" will be analyzed, and sample size must
# be explicitly specified.
# %samprmbygrp: $samprmbygrp{IID}{Group} = 1
# Group-specific removal list, this is derived from twin/sib pairs to keep samples unrelated 
# in the each sample groups used for summary 
{
	# Sample inclusion/exclusion list
	if (exists $conf{Sample}{Include}) {
		unless (-f $conf{Sample}{Include}) {
			die "Cannot find good sample list $conf{Sample}{Include}";
		}
		print STDERR "Reading sample inclusion list from $conf{Sample}{Include}\n";
		%sampincl = map { (split)[0] => 1 } slurp $conf{Sample}{Include};
	}
	if (exists $conf{Sample}{Exclude}) {
		unless (-f $conf{Sample}{Exclude}) {
			die "Cannot find bad sample list $conf{Sample}{Exclude}";
		}
		print STDERR "Reading sample exclusion list from $conf{Sample}{Exclude}\n";
		%samprm = map { (split)[0] => 1 } slurp $conf{Sample}{Exclude};
	}
	if (exists $conf{Sample}{FamExclude}) {
		unless (-f $conf{Sample}{FamExclude}) {
			die "Cannot find bad sample list $conf{Sample}{FamExclude}";
		}
		print STDERR "Reading family member exclusion list from $conf{Sample}{FamExclude}\n";
		%famrm = map { (split)[0] => 1 } slurp $conf{Sample}{FamExclude};
	}

	# Parse sample tables (or variant table when sample table is not available)
	my (@samptabs, $f_iid);
	if (exists $conf{Sample}{Table}) {
		if (ref $conf{Sample}{Table} eq 'ARRAY') {
			@samptabs = @{$conf{Sample}{Table}};
		}
		else {
			@samptabs = ($conf{Sample}{Table});
		}
		$f_iid = $conf{Sample}{SampID};
	}
	else {
		print STDERR "Gather sample info from variant table\n";
		if (ref $conf{Variant}{Table} eq 'ARRAY') {
			@samptabs = @{$conf{Variant}{Table}};
		}
		else {
			@samptabs = ($conf{Variant}{Table});
		}
		$f_iid = $conf{Variant}{SampID};
	}
	unless(all { -f $_ } @samptabs) {
		die "Not all sample table files exist";
	}

	my $filter;
	if (exists $conf{Sample}{Filter}) {
		my $filtexpr = $conf{Sample}{Filter};
		$filtexpr =~ s/^["']//; $filtexpr =~ s/["']$//;
		($filter, my $tokens) = sql_query($filtexpr, 1);
		foreach my $samptab (@samptabs) {
			my ($it, $fnames) = iter_file($samptab, { fsep => qr/\t/ });
			foreach my $tok (@$tokens) {
				if ($tok->[0] eq 'FIELD') {
					unless(grep { $tok->[1] eq $_} @$fnames) {
						die "Cannot find Filter field $tok->[1] in sample table $samptab";
					}
				}
			}
		}
	}

	# Add sample group membership to %sampgroup and tally group-specific %sampsize 
	my $callbacks;
	if (exists $conf{Sample}{Group}) {
		(my $samptypes, $callbacks, my $fields) = 
			parse_filters($conf{Sample}{Group}, $conf{Sample}{Group_Filter});
		foreach my $samptab (@samptabs) {
			my ($it, $fnames) = iter_file($samptab, { fsep => qr/\t/ });
			unless(grep { $f_iid eq $_ } @$fnames) {
				die "Cannot find sample ID field $f_iid from sample table";
			}
			foreach my $field (keys %$fields) {
				unless(grep { $field eq $_ } @$fnames) {
					die "Cannot find field $field from sample table";
				}
			}
		}
		@sampgroup = @$samptypes;
		# Init sample size
		foreach my $sgroup (@sampgroup) {
			$sampsize{$sgroup} = 0;
		}
	}
	else {
		unless(exists $conf{Sample}{Size}) {
			die "Must provide sample size when sample group is not defined";
		}
		unless($conf{Sample}{Size} =~ /^\d+$/) {
			die "Sample size for all must be a single integer!";
		}
		@sampgroup = qw|All|;
		# Set sample size to the pre-specified 
		$sampsize{All} = $conf{Sample}{Size};
	}

	# Parse sample table to determine samp size for each subgroup
	# Also need to check no sample overlaps in these sample tables
	my ($sampct, %known, %knowntabs);
	foreach my $samptab (@samptabs) {
		my ($it, $fnames) = iter_file($samptab, { fsep => qr/\t/ });
		while(my $dat = $it->()) {
			my $iid = $dat->{$f_iid};
			if (defined $filter) {
				# Sample failed sample-filter will be added to removal list
				unless ($filter->($dat)) {
					$samprm{$iid} = 1;
					next;
				}
			}
			if (%sampincl) {
				next unless defined $sampincl{$iid};
			}
			next if defined $samprm{$iid};
			if (defined $known{$iid}) {
				if (exists $conf{Sample}{Table}) {
					die "Sample $iid has appeard in the sample table previously";
				}
				else {
					if (defined $knowntabs{$iid} && $knowntabs{$iid} ne $samptab) {
						warn "Sample $iid has appeared in other variant table";
					}
					next;
				}
			}
			if (defined $callbacks) {
				for(my $ii = 0; $ii < @sampgroup; $ii ++) {
					if ($callbacks->[$ii]->($dat)) {
						$sampgroup{$iid}{$sampgroup[$ii]} = 1;
						$sampsize{$sampgroup[$ii]} ++;
					}
				}
			}
			else {
				$sampgroup{$iid}{All} = 1;
			}	
			$sampct ++;
			$known{$iid} = 1;
			$knowntabs{$iid} = $samptab;
		}
	}
	foreach my $sgroup (@sampgroup) {
		unless(defined $sampsize{$sgroup} && $sampsize{$sgroup} > 0) {
			die "No sample in group $sgroup was found";
		}
	}

	print STDERR "A total of $sampct samples are found in sample table.\n";
	
	# Check sample size if sample size is specified in the config
	if (exists $conf{Sample}{Size}) {
		my $sz;
		if ($conf{Sample}{Size} =~ /^\d+$/) {
			$sz = $conf{Sample}{Size};
			if (exists $conf{Sample}{Table}) {
				unless ($sz == $sampct) {
					die "Inconsistent sample size, $sampct found in sample table <> $sz specified in config";
				}
			}
			else {
				if($sampct > $sz) {
					die "Num of samples ($sampct) found in variant table > $sz specified in config";
				}
				elsif ($sampct < $sz) {
					warn "Num of samples ($sampct) found in variant table < $sz specified in config";
					if (@sampgroup > 1) {
						warn "You need to check consistency of sample sizes for other sample groups!";
					}
				}
			}
		}
		else {
			$sz = str2hash($conf{Sample}{Size}, { psep => ',', kvsep => ':' });
			foreach my $sgroup (@sampgroup) {
				unless(defined $sz->{$sgroup}) {
					die "Cannot find sample size for group $sgroup from sample size spec from config!";
				}
				if (defined $conf{Sample}{Table}) {
					unless($sz->{$sgroup} == $sampsize{$sgroup}) {
						die "Inconsistent sample size for group $sgroup, $sampsize{$sgroup} <> $sz->{$sgroup} specified in config";
					}
				}
				else {
					if ($sampsize{$sgroup} > $sz->{$sgroup}) {
						die "Num of samples for group $sgroup ($sampsize{$sgroup}) found in variant table > $sz->{$sgroup} specified in config";
					}
					elsif ($sampsize{$sgroup} < $sz->{$sgroup}) {
						warn "Num of samples for group $sgroup ($sampsize{$sgroup}) found in variant table < $sz->{$sgroup} specified in config";
						$sampsize{$sgroup} = $sz->{$sgroup};
						if (@sampgroup > 1) {
							warn "You need to check consistency of sample size for other sample groups!";
						}
					}
				}
			}
		}
	}

	# If twin/sibs are specified
	# When multiple related individuals exist in the same sample group, only first one will be kept.
	# This is achieved by adding samples to a group-specific removal list, so they will be skipped when 
	# tallying RG counts in that group.
	# We do not distinguish twins and sibs for the current implmentation
	if (defined $conf{Sample}{Twins} || defined $conf{Sample}{Sibs}) {
		my %twinsibs = read_twinsibs($conf{Sample}{Twins}, $conf{Sample}{Sibs});
		foreach my $iid (sort keys %sampgroup) {
			foreach my $sgroup (sort keys %{$sampgroup{$iid}}) {
				next if defined $samprmbygrp{$iid} && defined $samprmbygrp{$iid}{$sgroup};
				# For each sample and group, look for other sibs from the same group that have not been removed
				my @samegrp = grep { defined $sampgroup{$_} && defined $sampgroup{$_}{$sgroup} && 
									!(defined $samprmbygrp{$_} && defined $samprmbygrp{$_}{$sgroup}) } @{$twinsibs{$iid}};
				# They will be added to the group-specific removal list
				# Note that we will not remove their sample group membership
				# We also adjust group-specifc sample sizes accordingly
				foreach my $iid (@samegrp) {
					$samprmbygrp{$iid}{$sgroup} = 1;
					$sampsize{$sgroup} --;
				}
			}
		}
	}

	# Group pairs for comparison
	if (exists $conf{Sample}{Compare}) {
		if (ref $conf{Sample}{Compare} eq 'ARRAY') {
			foreach my $comp (@{$conf{Sample}{Compare}}) {
				push @gpairs => map { [ split(':') ] } split(',', $comp);
			}
		}
		else {
			@gpairs = map { [ split(':') ] } split(',', $conf{Sample}{Compare});
		}
	}
	foreach my $gpair (@gpairs) {
		my @groups = @$gpair;
		unless(@groups == 2 && (grep { $groups[0] eq $_ } @sampgroup) && 
							   (grep { $groups[1] eq $_ } @sampgroup)) {
			die "Cannot find comparison pair @{[ join(q|,|, @groups) ]} from sample groups";
		}
	}
	# Check that there is no sample overlap in contrasting group pairs
	foreach my $iid (keys %sampgroup) {
		foreach my $gpair (@gpairs) {
			if (defined $sampgroup{$iid}{$gpair->[0]} &&
		 		!(defined $samprmbygrp{$iid} && defined $samprmbygrp{$iid}{$gpair->[0]}) &&
			 	defined $sampgroup{$iid}{$gpair->[1]} &&
			 	!(defined $samprmbygrp{$iid} && defined $samprmbygrp{$iid}{$gpair->[1]})) {
				die "Sample $iid belongs to both constrasting groups: ".join(",", @$gpair);
			}
		}
	}
}

#
# ROH binkeeper for each individual
#
# Need to be sure that chromsome nomenclauture are the same as variant table!
# Format (4 columns): IID, CHROM, START, END
my (%roh, $chr_roh);
if (exists $conf{Sample}{ROH}) {
	open my $fin, $conf{Sample}{ROH} or die "Cannot read sample ROH";
	while(<$fin>) {
		my @dat = split;
		#unless(@dat == 4) {
		#	warn "Incorrect number of columns in ROH file: $_";
		#}
		#if (@dat > 4) {
		#	warn "More than four columns in ROH file: will use the first 4";
		#}
		if (@dat < 4) {
			die "Less than four columns in ROH file!";
		}
		my ($iid, $chrom, $start, $end) = @dat[0,1,2,3];
		next if defined $samprm{$iid};
		if (%sampincl) {
			next unless defined $sampincl{$iid};
		}
		unless(defined $chr_roh) {
			$chr_roh = $chrom =~ /^chr/ ? 1 : 0;
		}
		unless(defined $roh{$iid}) {
			$roh{$iid} = Genome::Ranges::IntSet->new();
		}
		$roh{$iid}->add($chrom, $start, $end);
	}
}

#
# Variant table and phase filter
#
my (@vartabs, @vartabfields, $f_iid, @varfilt, @varfiltfields, %varexcl, $hetfilt, $homfilt, @genofiltfields, @repulsefilt);
# @vartabs: One or more variant tables
# @vartabfields: common variant table fields found in all tables, will be used in dumped vartab
# 			It will be later updated to include additional gene level and sample level fields.
# @varfilt: an array of filters for parsing variant table
#			All filters in the array must be evaluate to true for variant to be included
# 			The first filter will be main site filter, all remaining are optional family-based filter
# @varfiltfields: an array of fields in the *first* filter defined above (used for exapnsion)			
# 
# $hetfilt : determine if a genotype is heterozygotes or not
# @repulsefilt: filters applies on expanded family data of two het sites to determine if they are in repulse phase
#			If any of them evaluates to be true, then two sites will be in repulse phase status
sub parse_genofilt {
	my ($filtexpr, $dupaction) = @_;
	$filtexpr=~ s/^["']//; $filtexpr =~ s/["']$//;
	my ($cb, $tokens) = sql_query($filtexpr, 1);
	foreach my $tok (@$tokens) {
		if ($tok->[0] eq 'FIELD') {
			unless(grep { $tok->[1] eq $_ } @vartabfields) {
				die "Cannot find Filter field $tok->[1] in variant table";
			}
			unless(grep { $tok->[1] eq $_ } @genofiltfields) {
				push @genofiltfields, $tok->[1];
			}
		}
	}
	my $callback;
	# the same variant level dupaction will be applied to het filter
	if ($dupaction eq 'Any') {
		$callback = sub { any { $cb->($_) } @_ };
	}
	elsif ($dupaction eq 'All') {
		$callback = sub { all { $cb->($_) } @_ };
	}
	elsif ($dupaction eq 'First') {
		$callback = sub {  $cb->($_[0])  };
	}
	else {
		die "Cannot recognize DupAction: $dupaction";
	}
	return $callback;
}

{
	if (ref $conf{Variant}{Table} eq 'ARRAY') {
		@vartabs = @{$conf{Variant}{Table}};
	}
	else {
		@vartabs = ($conf{Variant}{Table});
	}

	$f_iid = $conf{Variant}{SampID};


	# Overall site level filter and optional customization for SNVs and indels
	foreach my $label (qw|Filter SNV_Filter Indel_Filter MNV_Filter|) {
		if (exists $conf{Variant}{$label}) {
			my $filtexpr = $conf{Variant}{$label};
			$filtexpr =~ s/^["']//; $filtexpr =~ s/["']$//;
			my ($cb, $tokens) = sql_query($filtexpr, 1);
			if ($conf{Variant}{DupAction} eq 'Any') {
				$varfilt[0]{$label} = sub { any { $cb->($_) } @_ };
			}
			elsif ($conf{Variant}{DupAction} eq 'All') {
				$varfilt[0]{$label} = sub { all { $cb->($_) } @_ };
			}
			elsif ($conf{Variant}{DupAction} eq 'First') {
				$varfilt[0]{$label} = sub {  $cb->($_[0])  };
			}
			else {
				die "Cannot recognize DupAction: $conf{Variant}{DupAction}";
			}
			foreach my $tok (@$tokens) {
				if ($tok->[0] eq 'FIELD') {
					unless(grep { $tok->[1] eq $_ } @varfiltfields) {
						push @varfiltfields, $tok->[1];
					}
				}
			}
		}
		else {
			$varfilt[0]{$label} = undef;
		}
	}

	# Slurp variant exclusion list if any
	if (exists $conf{Variant}{Exclude}) {
		open my $fin, $conf{Variant}{Exclude} or die "Cannot open $conf{Variant}{Exclude}";
		while(<$fin>) {
			my ($varid, $iid) = split;
			if (defined $iid) {
				$varexcl{$iid,$varid} = 1;
			}
			else {
				$varexcl{$varid} = 1;
			}
		}
	}	

	my @allfnames;
	#my %knownids; # IID => VarTab
	foreach my $vartab (@vartabs) {
		unless(-f $vartab) {
			die "Cannot find variant table $vartab!";
		}
		print STDERR "Validating variant table $vartab\n";
		my ($it, $fnames) = iter_file($vartab, { fsep => qr/\t/ });
		push @allfnames, $fnames;
		# check standard fields
		foreach my $stdfd (qw|Chrom Position Ref Alt GeneID|, $f_iid) {
			unless (grep { $_ eq $stdfd } @$fnames) {
				die "Cannot find standard field $stdfd in variant table $vartab";
			}
		}
		# Check filtering fields for the main site level filter
		foreach my $fd (@varfiltfields) {
			unless(grep { $fd eq $_ } @$fnames) {
				die "Cannot find Filter field $fd in variant table $vartab";
			}
		}

		# Skip validating variant table if re-use existing intermediate results
		# The only validation here is to ensure samples from different variant tables have no overlaps
		#next if $reuse;
		#if (@vartabs > 1) {
		#	my %sampids;
		#	while(my $dat = $it->()) {
		#		my $iid = $dat->{$f_iid};
		#		if (defined $knownids{$iid}) {
		#			die "Sample $iid in $vartab appear in table $knownids{$iid}";
		#		}
		#		$sampids{$iid} = 1;
		#	}
		#	foreach my $iid (keys %sampids) {
		#		$knownids{$iid} = $vartab;
		#	}
		#}
	}
	
	# Merge variant table field names in common
	if (@allfnames > 1) {
		foreach my $field (@{$allfnames[0]}) {
			my $ntimes;
			for(my $ii = 1; $ii < @allfnames; $ii ++) {
				if (grep { $field eq $_ } @{$allfnames[$ii]}) {
					$ntimes ++;
				}
			}
			if ($ntimes == @allfnames-1) {
				push @vartabfields, $field;
			}
		}
	}
	else {
		@vartabfields = @{$allfnames[0]};
	}
	unless(@vartabfields) {
		die "No fields are in common by variant tables!";
	}

	# If any FamFilter is defined, we will parse the filter and check if related fields exist
	if (exists $conf{Variant}{FamFilter}) {
		my ($conditions, $callbacks, $fields) = 
			parse_filters($conf{Variant}{FamFilter_Cond}, $conf{Variant}{FamFilter});
		foreach my $cb (@$callbacks) {
			my $cond = shift @$conditions;
			if ($cond eq 'Any') {
				push @varfilt => sub { any { $cb->($_) } @_ };
			}
			elsif ($cond eq 'All') {
				push @varfilt => sub { all { $cb->($_) } @_ };
			}
			elsif ($cond eq 'None') {
				push @varfilt => sub { none { $cb->($_) } @_ };
			}
			elsif ($cond eq 'Notall') {
				push @varfilt => sub { notall { $cb->($_) } @_ };
			}
			elsif ($cond eq 'None') {
				push @varfilt => sub { none { $cb->($_) } @_ };
			}
			else {
				die "Cannot recognize condition: $cond";
			}
		}
		
		foreach my $reqfd (grep { $_ ne 'Relation' && $_ ne 'Pheno' && $_ ne 'IID' } keys %$fields) {
			$reqfd =~ s/^_// if $reqfd =~ /^_/;
			unless(grep { $reqfd.".FamMembers" eq $_ || $reqfd eq $_  } @vartabfields) {
				die "Cannot find FamFilter field $reqfd in variant table";
			}
		}
	}

	# Now parse repulse filter, filter fields will have .1/.2 suffix
	if (exists $conf{Phase}) {
		$hetfilt = parse_genofilt($conf{Phase}{HetFilter}, $conf{Phase}{HetFilter_DupAction});
		$homfilt = parse_genofilt($conf{Phase}{HomFilter}, $conf{Phase}{HomFilter_DupAction});

		my @repulse;
		if (ref $conf{Phase}{Repulse} eq 'ARRAY') {
			@repulse = @{$conf{Phase}{Repulse}};
		}
		else {
			@repulse = ($conf{Phase}{Repulse});
		}
		foreach my $repulse (@repulse) {
			$repulse =~ s/^["']//; $repulse =~ s/["']$//;
			my ($cb, $tokens) = sql_query($repulse, 1);
			foreach my $tok (@$tokens) {
				if ($tok->[0] eq 'FIELD') {
					my $field = $tok->[1]; $field =~ s/\.[12]$//;
					next if $field =~ /^(Relation|Pheno|IID)$/;
					unless(grep { $field.".FamMembers" eq $_ } @vartabfields) {
						die "Cannot find Repulse filter field $field in variant table";
					}
				}
			}
			push @repulsefilt, $cb;
		}
	}
	else {
		die "Must provide filter for phasing!"
	}
} 

#
# Variant classes parsers
#
# Note 1: All variant classes used in rg_burden are associated with genes!
# Note 2: we do not filter nearby variants as in rarevar_burden, but use family inheritance
# to resolve significance.
# Note 3: We also do not need a effect modifier that are used in dnv_burden or rarevar_burden
# but will treat overlapping genes by keep track of all RGs identified
my (@varclass, @varclsfilt, @varclsfiltfields);
{
	my ($varcls, $filters, $fields) = parse_filters($conf{Variant}{Class}, $conf{Variant}{Class_Filter});
	@varclsfiltfields = keys %$fields;

	@varclass = @$varcls;
	@varclsfilt = @$filters;

	foreach my $reqfd (keys %$fields) {
		unless(grep { $reqfd eq $_ } @vartabfields) {
			die "Cannot find variant class filter field $reqfd in variant table";
		}
	}
}

#
# Gene list (CAF) and gene sets
#
my (%caf, @varclswcaf, %geneset, @geneset, %generm, %genechr, %generng);
# We will collect the full gene list from CAF table, anything not in this table will not be
# included in the burden analysis.
# 
# %caf is a has of hash to store per-gene mutation rate $mutrate{GeneID}{VarClass} = rate
# Not all vairant class defined will have CAF
# @varclswcaf is the list of variant class with CAF, it can be empty, in which case we will not 
# compare counts with expecations predicted by CAF. But variant class wo CAF can still be used in
# between group comparison!
#
# %generm: will be the blacklist genes to be removed from analysis
# %genesets: is a hash of hash to hold geneset membership for each gene $geneset{Set}{GeneID} = 1
# @geneset: is an array of all gene sets starting with the defalt "All".
# %genechr: will store chrom for each gene
# %generng: will store coding regions for each gene, this is used to overlap ROH regions
{
	if (exists $conf{Gene}{Exclude}) {
		%generm = map { (split)[0] => 1 } slurp $conf{Gene}{Exclude};
	}
	my ($it, $fnames) = iter_file($conf{CAF}{Table}, { fsep => qr/\t/ });
	my $f_gid = $conf{CAF}{GeneID};
	unless(grep { $f_gid eq $_ } @$fnames) {
		die "Cannot find gene ID field $f_gid from mutation rate table";
	}

	my $varcls;
	if (exists $conf{CAF}{VarClass}) {
		$varcls = parse_fstr($conf{CAF}{VarClass}, 1);
		foreach my $class (keys %$varcls) {
			unless(grep { $_ eq $class } @varclass) {
				die "Cannot find definition of variant class: $class!";
			}
			my $colname = $varcls->{$class};
			unless (grep { $colname eq $_ } @$fnames) {
				die "Cannot find field $colname in mutation rate table!";
			}
			push @varclswcaf, $class;
		}
	}

	# Collect CAF for each variant class
	while(my $dat = $it->()) {
		my $gid = $dat->{$f_gid};
		next if defined $generm{$gid};
		if (defined $caf{$gid}) {
			warn "Mutation rate for gene $gid has already been defined";
			next;
		}
		# By default, genes that appear in the CAF table will be in "All" set
		$geneset{All}{$gid} = 1;
		unless(defined $varcls) {
			$caf{$gid} = {};
			next;
		}
		foreach my $class (sort keys %$varcls) {
			my $colname = $varcls->{$class};
			if ($dat->{$colname} > 0) {
				$caf{$gid}{$class} = $dat->{$colname};
			}
			elsif ($dat->{$colname} == 0) {
				# For CAF rate of zero, assign a value that is abitarily small
				$caf{$gid}{$class} = $conf{CAF}{Zero};
			}
			else {
				die "Cumulative allele freq cannot be less than 0: $colname = $dat->{$colname}";
			}	
		}
	}

	# Gene sets other than all
	@geneset = qw|All|;
	if (exists $conf{Gene}{Set}) {
		my $gs = read_geneset($conf{Gene}{Set}, $conf{Gene}{Set_Fields});
		while(my ($set, $genes) = each %$gs) {
			if ($set eq 'All') {
				die "Gene set All is reserved!";
			}
			foreach my $gid (@$genes) {
				# Restrict to genes with mutation rate!
				next unless defined $caf{$gid};
				$geneset{$set}{$gid} = 1;
			}
		}
		push @geneset, sort grep { $_ ne 'All' } keys %geneset;
	}

	my $git;
	if (-f $conf{Gene}{DBTable}) {
		$git = iter_genePred($conf{Gene}{DBTable});

	}
	else {
		$git = iter_geneTab($conf{Gene}{HgBuild}, $conf{Gene}{DBTable});
	}
	# Consider only autosomal coding genes!
	while(my $dat = $git->()) {
		next unless $dat->{chrom} =~ /^chr[0-9]+$/;
		next unless $dat->{cdsStart} < $dat->{cdsEnd};
		my $chr = $dat->{chrom};
		unless ($chr_roh) {
			$chr =~ s/^chr//;
		}
		chk_default(\%genechr, $dat->{name2}, $chr);
		foreach my $intv (gdat_list_coding_exons($dat)) {
			push @{$generng{$dat->{name2}}}, [$intv->[0], $intv->[1]];
		}
	}
}


# 
# Expected number of RGs for different class of variants
# 
my (%expctbygene, %expct);
if (@varclswcaf) {
	# For each samples in the group:
	# if he/she does not have ROH overlapping this gene, we will add CAF^2
	# if he/she has ROH with x% overlap of CDS of the gene, we will add CAF^2*(1-x%) + x%*CAF
	# Going through all variant classes and genes defined above 
	foreach my $vclass (@varclswcaf) {
		foreach my $gid (sort keys %caf) {
			my $gchr = $genechr{$gid} // do { die "Cannot find chromosome for gene $gid" };
			my $gcds = $generng{$gid} // do { die "Cannot find coding regions for gene $gid" };
			foreach my $iid (sort keys %sampgroup) {
				my $overlaproh;
				if (defined $roh{$iid} && defined $roh{$iid}{$gchr}) {
					my $iroh = $roh{$iid}{$gchr};
					$overlaproh = overlap_rohgene($iroh, $gcds);
				}
				else {
					$overlaproh = 0;
				}
				foreach my $sgroup (sort keys %{$sampgroup{$iid}}) {
					next if defined $samprmbygrp{$iid} && defined $samprmbygrp{$iid}{$sgroup};
					if ($overlaproh > 0) {
						$expctbygene{$sgroup}{$gid}{$vclass} = ($overlaproh * $caf{$gid}{$vclass} + (1-$overlaproh) * $caf{$gid}{$vclass} ** 2) * $sampsize{$sgroup};
					}
					else {
						$expctbygene{$sgroup}{$gid}{$vclass} = $caf{$gid}{$vclass} ** 2 * $sampsize{$sgroup};
					}
				}
			}
		}
		# Then we can sum up the gene-specific expectations to get gene-set exp counts
		foreach my $gset (@geneset) {
			foreach my $gid (sort keys %{$geneset{$gset}}) {
				foreach my $sgroup (@sampgroup) {
					$expct{$sgroup}{$gset}{$vclass} += $expctbygene{$sgroup}{$gid}{$vclass};
				}
			}
		}
	}
}

sub overlap_rohgene {
	my ($roh, $cds) = @_;
	my $query = Set::IntSpan->new($roh);
	my $target = Set::IntSpan->new($cds);
	my $overlap = $query * $target;
	my $iLen = $overlap->size();
	my $tLen = $target->size();
	if ($tLen == 0) {
		die "Coding gene has zero CDS length!";
	}
	return $iLen/$tLen;
}

#
# Initial filtering and variant table sorting
#
# A temp filtered variant table will be generated. 
# The site level and family based filter will be applied to input variant table.
# If duplicated variant-IID combination exist, only the first one will be kept.
# In case a variant was annotated to multiple genes, they will be *EXPANDED* to one row per gene.
# Then we sort the table by GeneID then by IID, so all variants within an individual will appear in
# consecutive rows.
# Note: in rg_burden, we split the overlapping genes and tally RGs on each of them.
# The expected number of RGs will be calculated in the same way. 
# It may slighly inflate the statistics, but given the rarity of event, it will be a reasonable approximation.
my $vartab = "$wrkdir/vartab.txt.gz";
if (-f $vartab) {
	print STDERR "Re-use the filter-sorted variant table in existing working directory $wrkdir/vartab.txt.gz\n";
}
else {
	my $bk;
	if (exists $conf{Variant}{Region}) {
		unless(-f $conf{Variant}{Region}) {
			die "Cannot find region file: $conf{Variant}{Region}";
		}
		$bk = Genome::UCSC::BinKeeper->new($conf{Variant}{Region}, { bed => 1});
	}

	# Go through each vartab, first apply filters and create temp file to sort variant tab
	my $c_gid = first { $vartabfields[$_-1] eq 'GeneID' } 1..scalar(@vartabfields);
	my $c_iid = first { $vartabfields[$_-1] eq $f_iid } 1..scalar(@vartabfields);
	my $c_chr = first { $vartabfields[$_-1] eq 'Chrom' } 1..scalar(@vartabfields);
	my $c_pos = first { $vartabfields[$_-1] eq 'Position' } 1..scalar(@vartabfields);
	unless(defined $c_gid && defined $c_iid && defined $c_chr && defined $c_pos) {
		die "Cannot find Gene ID, sample ID, Chrom, Position from variant table";
	}
	open my $fout, "| body sort -k $c_gid,$c_gid -k $c_iid,$c_iid -k $c_chr,$c_chr -k $c_pos,${c_pos}n -T $wrkdir/tmp | bgzip -c > $wrkdir/vartab.txt.gz"
			or die "Cannot open body sort pipe to write to $wrkdir/vartab.txt.gz";
	print $fout join("\t", @vartabfields), "\n";

	for(my $ii = 0; $ii < @vartabs; $ii ++) {
		my $vartab = $vartabs[$ii];
		my $it = iter_file($vartab, { fsep => qr/\t/ });
		print "Initial filtering on variant table $vartab\n"; 
		
		my %known;
		while(my $dat = $it->()) {
			my $iid = $dat->{$f_iid};
			if (%sampincl) {
				next unless defined $sampincl{$iid};
			}
			next if defined $samprm{$iid};

			my $varid = join(":", @{$dat}{qw|Chrom Position Ref Alt|});
			if (defined $known{$iid,$varid}) {
				warn "Duplicated variant found in the variant table $vartab: $iid,$varid";
				next;
			}
			# Variant exclusion
			if (%varexcl) {
				next if defined $varexcl{$iid,$varid} || defined $varexcl{$varid};
			}
			# Region-based inclusion
			if ($bk) {
				next unless $bk->any_overlap($dat->{Chrom}, $dat->{Position}, $dat->{Position}+length($dat->{Ref})-1);
			}

			my @data4filt = expand_dat($dat, { sep => ',', optional => \@varfiltfields });
			die "Data cannot be expanded at $dat->{Chrom}:$dat->{Position}!" unless @data4filt;
			
			if (defined $varfilt[0]{Filter}) {
				next unless $varfilt[0]{Filter}->(@data4filt);
			}
			my $VT = var_type($dat->{Ref}, $dat->{Alt});
			if (defined $varfilt[0]{"${VT}_Filter"}) {
				next unless $varfilt[0]{"${VT}_Filter"}->(@data4filt);
			}

			# If more than one filter is specified, the remaining will be FamMember-based
			my $failfilt;
			if (@varfilt > 1) {
				my @famdata = expand_fam($dat, \%famrm);
				# Only use the family-based filter when family data are available!
				if (@famdata) {
					for(my $jj = 1; $jj < @varfilt; $jj ++) {
						unless ($varfilt[$jj]->(@famdata)) {
							$failfilt = 1;
						}
					}
				}
			}
			next if $failfilt;

			# Write to output if the variant passed both basic and fam-based filters
			# We will split data if multiple gene IDs are present
			my @data = expand_dat($dat, { sep => ';', optional => ["GeneID", @varclsfiltfields] });
			foreach my $data (grep { $_->{GeneID} ne '.' } @data) {
				print $fout join("\t", @{$data}{@vartabfields}), "\n";
			}		
			$known{$iid,$varid} = 1;
		}
	}
	close $fout;
}


# 
# Tally recessive genotypes
# 
my %rg;
# %rg will store RGs found for all variant classes in all individuals
# $rg{GeneID}{VarClass}{IID} = "Hom"/"cHet"

# We will select some variant class / sample group to dump to output
my ($fdump, @ginfofields, $dumpxtraginfo, @sinfofields, $dumpxtrasinfo);
if (exists $conf{Dump}) {
	if (exists $conf{Dump}{GXref}) {
		my ($xgdat, $xgfds) = slurp_xref($conf{Dump}{GXref}, $conf{Dump}{GXref_Fields});
		$dumpxtraginfo = $xgdat;
		@ginfofields = @$xgfds;
		foreach my $gfield (@ginfofields) {
			if (grep { $gfield eq $_ } @vartabfields) {
				warn "Gene info field $gfield already appear in the original variant table!";
			}
		}
	}
	if (@ginfofields) {
		insert_after(\@vartabfields, 'GeneID', @ginfofields);
	}
	if (exists $conf{Dump}{SXref}) {
		my ($xsdat, $xsfds) = slurp_xref($conf{Dump}{SXref}, $conf{Dump}{SXref_Fields});
		$dumpxtrasinfo = $xsdat;
		@sinfofields = @$xsfds;
		foreach my $sfield (@sinfofields) {
			if (grep { $sfield eq $_ } @vartabfields) {
				warn "Sample info field $sfield already appear in the original variant table!";
			}
		}
	}
	if (@sinfofields) {
		insert_after(\@vartabfields, $f_iid, @sinfofields);
	}
	
	open $fdump, ">$output.RGs.txt" or die "Cannot write to $output.RGs.txt";
	print $fdump join("\t", @vartabfields), "\n";
}

my ($prev_gid, $prev_iid, @iHetVars, @iHomVars);

my $it = iter_file($vartab, { fsep => qr/\t/ });
while(my $dat = $it->()) {
	# find all variants of the same individual
	if (defined $prev_gid && $prev_gid ne $dat->{GeneID} ||
		defined $prev_iid && $prev_iid ne $dat->{$f_iid}) {
		# Go-through each variants of this individual, find out variants belonging to each class
		for(my $ii = 0; $ii < @varclass; $ii ++) {
			my @clsHomVars = grep { $varclsfilt[$ii]->($_) } @iHomVars;
			my @clsHetVars = grep { $varclsfilt[$ii]->($_) } @iHetVars;
			next if @clsHomVars == 0 && @clsHetVars == 0;
			if (@clsHomVars > 0) {
				$rg{$prev_gid}{$varclass[$ii]}{$prev_iid} = "Hom";
				if (defined $fdump && exists $conf{Dump}{SampGroup} && defined $sampgroup{$prev_iid}{$conf{Dump}{SampGroup}} &&
					exists $conf{Dump}{GeneSet}  && defined $geneset{$conf{Dump}{GeneSet}}{$prev_gid} &&
					exists $conf{Dump}{VarClass} && $varclass[$ii] eq $conf{Dump}{VarClass}) {
					foreach my $vardat (@clsHomVars) {	
						if (defined $dumpxtraginfo) {
							foreach my $gfd (@ginfofields) {
								$vardat->{$gfd} = $dumpxtraginfo->{$prev_gid}{$gfd} // ".";
							}			
						}
						if (defined $dumpxtrasinfo) {
							foreach my $sfd (@sinfofields) {
								$vardat->{$sfd} = $dumpxtrasinfo->{$prev_iid}{$sfd} // ".";
							}
						}
						print $fdump join("\t", @{$vardat}{@vartabfields}), "\n";	
					}
				}
			}
			elsif (@clsHetVars > 1) {
				ALLPAIRS:foreach my $pair (map { my @fam1 = expand_fam($_->[0], \%famrm); 
												 my @fam2 = expand_fam($_->[1], \%famrm);
												 ([\@fam1, \@fam2], [\@fam2, \@fam1]) } 
												 	all_pairs(@clsHetVars)) {
					unless(scalar(@{$pair->[0]}) == scalar(@{$pair->[1]})) {
						die "Incorrect number of elements in expanded family data at $prev_gid:$prev_iid";
					}

					for(my $kk = 0; $kk < scalar(@{$pair->[0]}); $kk ++) {
						unless($pair->[0][$kk]{IID} eq $pair->[1][$kk]{IID}) {
							die "Samples from two family expansions do not align: $pair->[0]{$kk}{IID} ne $pair->[1]{$kk}{IID}";
						}
						my %ipair;
						foreach my $ii (0,1) {
							my $jj = $ii + 1;
							foreach my $key (keys %{$pair->[$ii][$kk]}) {
								$ipair{"$key.$jj"} = $pair->[$ii][$kk]{$key};
							}
						}
						if (any { $_->(\%ipair) } @repulsefilt) {
							$rg{$prev_gid}{$varclass[$ii]}{$prev_iid} = "cHet";
							if (defined $fdump &&  exists $conf{Dump}{SampGroup} && defined $sampgroup{$prev_iid}{$conf{Dump}{SampGroup}} &&
								exists $conf{Dump}{GeneSet} && defined $geneset{$conf{Dump}{GeneSet}}{$prev_gid} &&
								exists $conf{Dump}{VarClass} && $varclass[$ii] eq $conf{Dump}{VarClass}) {
								foreach my $vardat (@clsHetVars) {
									if (defined $dumpxtraginfo) {
										foreach my $gfd (@ginfofields) {
											$vardat->{$gfd} = $dumpxtraginfo->{$prev_gid}{$gfd} // ".";
										}
									}
									if (defined $dumpxtrasinfo) {
										foreach my $sfd (@sinfofields) {
											$vardat->{$sfd} = $dumpxtrasinfo->{$prev_iid}{$sfd} // ".";
										}
									}
									print $fdump join("\t", @{$vardat}{@vartabfields}), "\n";
								}
							}
							last ALLPAIRS;
						}
					}
				}
			}
		}
		@iHetVars = ();
		@iHomVars = ();
	}
	if ($hetfilt->($dat)) {
		push @iHetVars, $dat;
	}
	elsif ($homfilt->($dat)) {
		push @iHomVars, $dat;
	}
	else {
		next;
	}
	$prev_gid = $dat->{GeneID};
	$prev_iid = $dat->{$f_iid};
}
# Last gene
if (@iHetVars > 0 || @iHomVars > 0) {
	for(my $ii = 0; $ii < @varclass; $ii ++) {
		my @clsHomVars = grep { $varclsfilt[$ii]->($_) } @iHomVars;
		my @clsHetVars = grep { $varclsfilt[$ii]->($_) } @iHetVars;
		next if @clsHomVars == 0 && @clsHetVars == 0;
		if (@clsHomVars > 0) {
			$rg{$prev_gid}{$varclass[$ii]}{$prev_iid} = "Hom";
			if (defined $fdump && 
				exists $conf{Dump}{SampGroup} && defined $sampgroup{$prev_iid}{$conf{Dump}{SampGroup}} &&
				exists $conf{Dump}{GeneSet} && defined $geneset{$conf{Dump}{GeneSet}}{$prev_gid} &&
				exists $conf{Dump}{VarClass} && $varclass[$ii] eq $conf{Dump}{VarClass}) {
				foreach my $vardat (@clsHomVars) {	
					if (defined $dumpxtraginfo) {
						foreach my $gfd (@ginfofields) {
							$vardat->{$gfd} = $dumpxtraginfo->{$prev_gid}{$gfd} // ".";
						}			
					}
					if (defined $dumpxtrasinfo) {
						foreach my $sfd (@sinfofields) {
							$vardat->{$sfd} = $dumpxtrasinfo->{$prev_iid}{$sfd} // ".";
						}
					}
					print $fdump join("\t", @{$vardat}{@vartabfields}), "\n";	
				}
			}
		}
		elsif (@clsHetVars > 1) {
			ALLPAIRS:foreach my $pair (map { my @fam1 = expand_fam($_->[0], \%famrm); 
											 my @fam2 = expand_fam($_->[1], \%famrm);
											 ([\@fam1, \@fam2], [\@fam2, \@fam1]) } 
											 	all_pairs(@clsHetVars)) {
				unless(scalar(@{$pair->[0]}) == scalar(@{$pair->[1]})) {
					die "Incorrect number of family data at $prev_gid:$prev_iid";
				}

				for(my $kk = 0; $kk < scalar(@{$pair->[0]}); $kk ++) {
					unless($pair->[0][$kk]{IID} eq $pair->[1][$kk]{IID}) {
						die "Samples from two family expansions do not align: $pair->[0]{$kk}{IID} ne $pair->[1]{$kk}{IID}";
					}
					my %ipair;
					foreach my $ii (0,1) {
						my $jj = $ii + 1;
						foreach my $key (keys %{$pair->[$ii][$kk]}) {
							$ipair{"$key.$jj"} = $pair->[$ii][$kk]{$key};
						}
					}
					if (any { $_->(\%ipair) } @repulsefilt) {
						$rg{$prev_gid}{$varclass[$ii]}{$prev_iid} = "cHet";
						if (defined $fdump && 
							exists $conf{Dump}{SampGroup} && defined $sampgroup{$prev_iid}{$conf{Dump}{SampGroup}} &&
							exists $conf{Dump}{GeneSet} && defined $geneset{$conf{Dump}{GeneSet}}{$prev_gid} &&
							exists $conf{Dump}{VarClass} && $varclass[$ii] eq $conf{Dump}{VarClass}) {
							foreach my $vardat (@clsHetVars) {
								if (defined $dumpxtraginfo) {
									foreach my $gfd (@ginfofields) {
										$vardat->{$gfd} = $dumpxtraginfo->{$prev_gid}{$gfd} // ".";
									}
								}
								if (defined $dumpxtrasinfo) {
									foreach my $sfd (@sinfofields) {
										$vardat->{$sfd} = $dumpxtrasinfo->{$prev_iid}{$sfd} // ".";
									}
								}
								print $fdump join("\t", @{$vardat}{@vartabfields}), "\n";
							}
						}
						last ALLPAIRS;
					}
				}
			}
		}
	}
}

# Perform aggregation %rg => %obsctbygene, %obsct
# $obsct{$sgroup}{$gset}{$vclass} = {tot/cHet/Hom count},
# $obsctbygene{$sgroup}{$gid}{$vclass} = {tot/cHet/Hom count}, 
# obsctbysamp{$iid}{$gset}{$vclass} = {tot/cHet/Hom count}.
my (%obsct, %obsctbygene, %obsctbysamp);
# CAF defines the full gene set
foreach my $gid (sort keys %caf) {
	foreach my $vclass (@varclass) {
		if (defined $rg{$gid}{$vclass}) {
			my @rgsamps = keys %{$rg{$gid}{$vclass}}; 
			foreach my $iid (@rgsamps) {
				foreach my $gset (@geneset) {
					if (defined $geneset{$gset}{$gid}) {
						$obsctbysamp{$iid}{$gset}{$vclass}{Tot} += 1;
						$obsctbysamp{$iid}{$gset}{$vclass}{Hom} += 1;
						$obsctbysamp{$iid}{$gset}{$vclass}{cHet} += 1;
					}
				}
			}
			foreach my $sgroup (@sampgroup) {
				my @sampingrp = grep { defined $sampgroup{$_}{$sgroup} &&
					!(defined $samprmbygrp{$_} && defined $samprmbygrp{$_}{$sgroup}) } @rgsamps;
				my $nTot = scalar(@sampingrp);
				my $nHom = scalar(grep { $rg{$gid}{$vclass}{$_} eq 'Hom' } @sampingrp);
				my $ncHet = scalar(grep { $rg{$gid}{$vclass}{$_} eq 'cHet' } @sampingrp);
				$obsctbygene{$sgroup}{$gid}{$vclass} = { Tot => $nTot, Hom => $nHom, cHet => $ncHet };
				# Also add gene set counts
				foreach my $gset (@geneset) {
					if (defined $geneset{$gset}{$gid}) {
						$obsct{$sgroup}{$gset}{$vclass}{Tot} += $nTot;
						$obsct{$sgroup}{$gset}{$vclass}{Hom} += $nHom;
						$obsct{$sgroup}{$gset}{$vclass}{cHet} += $ncHet;
					}
				}
			}
		}
		else {
			foreach my $sgroup (@sampgroup) {
				$obsctbygene{$sgroup}{$gid}{$vclass} = { Tot => 0, Hom => 0, cHet => 0 };
			}
		}
	}
}

# 
# Statistical test and summary
# 
my $R = Statistics::R->new();

# Summarize observed count and rate in each gene set and variant class
# and compare with expected count predicted by CAF if available
open my $fout, ">$output.burden.txt" or die "Cannot write $output.burden.txt";
print $fout join("\t", qw|SampleGroup N_Samps GeneSet N_Genes VarClass ObsCount ObsCntHom ObsCntcHet ObsRate
						  ExpCount ExpRate FoldOvsE P-value|), "\n";
foreach my $sgroup (@sampgroup) {
	foreach my $gset (@geneset) {
		foreach my $vclass (@varclass) {
			my $n_obs = $obsct{$sgroup}{$gset}{$vclass}{Tot} // 0;
			my $r_obs = $n_obs / $sampsize{$sgroup};
			my ($n_exp, $r_exp, $fold, $pval);
			if (grep { $vclass eq $_ } @varclswcaf) {
				$n_exp = $expct{$sgroup}{$gset}{$vclass};
				$r_exp = $n_exp / $sampsize{$sgroup};
				if ($r_exp > 0) {
					$fold = sprintf("%.3f", $n_obs/$n_exp);
					$R->set('obs', $n_obs);
					$R->set('exp', $n_exp);
					$R->run('pval<-poisson.test(obs, exp, alternative="two.sided")$p.value');
					$pval = sprintf("%.2E", $R->get('pval'));
				}
				else {
					$fold = ".";
					$pval = ".";
				}
			}
			else {
				$r_exp = ".";
				$n_exp = ".";
				$fold = ".";
				$pval = ".";
			}
			print $fout join("\t",  $sgroup, $sampsize{$sgroup}, $gset, scalar(keys %{$geneset{$gset}}), 
					$vclass, $n_obs, $obsct{$sgroup}{$gset}{$vclass}{Hom} // 0, $obsct{$sgroup}{$gset}{$vclass}{cHet} // 0,
					sprintf("%.4f", $r_obs), $n_exp eq '.' ? $n_exp : sprintf("%.1f", $n_exp), 
					$r_exp eq '.' ? $r_exp : sprintf("%.4f", $r_exp), $fold, $pval), "\n";
		}
	}
}

# Between group comparison (compare rate of all RGs)
if (@gpairs > 0) {
	open my $fcomp, ">$output.contra.txt" or die "Cannot write to $output.contra.txt";
	print $fcomp join("\t", qw|SampGroup1 N_Samps1 SampGroup2 N_Samps2 GeneSet N_Genes VarClass 
							   Count1 CntHom1 CntcHet1 Rate1 Count2 CntHom2 CntcHet2 Rate2 RateRatio P-value|), "\n";
	foreach my $gpair (@gpairs) {
		my @groups = @$gpair;
		# Compare observed rate per group
		foreach my $gset (@geneset) {
			foreach my $vclass (@varclass) {
				my $n_obs1 = $obsct{$groups[0]}{$gset}{$vclass}{Tot} // 0;
				my $r_obs1 = $n_obs1 / $sampsize{$groups[0]};
				my $n_obs2 = $obsct{$groups[1]}{$gset}{$vclass}{Tot} // 0;
				my $r_obs2 = $n_obs2 / $sampsize{$groups[1]};
				my ($ratio, $pval);
				if ($r_obs2 > 0) {
					$ratio = sprintf("%.3f", $r_obs1/$r_obs2);
					$R->set('x', $n_obs1);
					$R->set('n', $n_obs1+$n_obs2);
					$R->set('p', $sampsize{$groups[0]}/($sampsize{$groups[0]} + $sampsize{$groups[1]}));
					$R->run('pval<-binom.test(x, n, p, alternative="two.sided")$p.value');
					$pval = sprintf("%.2E", $R->get('pval'));
				}
				else {
					$ratio = ".";
					$pval = ".";
				}
				print $fcomp join("\t", $groups[0], $sampsize{$groups[0]}, $groups[1], $sampsize{$groups[1]}, $gset, scalar(keys %{$geneset{$gset}}), $vclass, 
								$n_obs1, $obsct{$groups[0]}{$gset}{$vclass}{Hom} // 0, $obsct{$groups[0]}{$gset}{$vclass}{cHet} // 0, sprintf("%.4f", $r_obs1),
								$n_obs2, $obsct{$groups[1]}{$gset}{$vclass}{Hom} // 0, $obsct{$groups[1]}{$gset}{$vclass}{cHet} // 0, sprintf("%.4f", $r_obs2), $ratio, $pval), "\n";
			}
		}
	}
}

# List gene-specific counts for selected "sample group - gene class" combination
if (exists $conf{GeneTab}) {
	my (@groupxclass, @header);
	foreach my $sgroup (split(',', $conf{GeneTab}{SampGroup})) {
		unless(grep { $sgroup eq $_ } @sampgroup) {
			die "Cannot find $sgroup in predefined sample groups";
		}
		foreach my $vclass (split(',', $conf{GeneTab}{VarClass})) {
			unless(grep { $vclass eq $_ } @varclass) {
				die "Cannot find $vclass in predefined variant class";
			}
			push @groupxclass, [$sgroup, $vclass];
			if (grep { $vclass eq $_ } @varclswcaf) {
				push @header, "${sgroup}_${vclass}_Tot", "${sgroup}_${vclass}_Hom", "${sgroup}_${vclass}_cHet",  "${sgroup}_${vclass}_Exp";
			}
			else {
				push @header, "${sgroup}_${vclass}_Tot", "${sgroup}_${vclass}_Hom", "${sgroup}_${vclass}_cHet";
			}
		}
	}
	# Additional gene level information
	my ($xtraginfo, @gfields);
	if (exists $conf{GeneTab}{GXref}) {
		my ($xgdat, $xgfds) = slurp_xref($conf{GeneTab}{GXref}, $conf{GeneTab}{GXref_Fields});
		$xtraginfo = $xgdat;
		@gfields = @$xgfds;
		foreach my $gfield (@gfields) {
			if (grep { $gfield eq $_ } @header) {
				warn "Field $gfield already appear in the genetab output!";
			}
		}
	}
	open my $fgen, ">$output.genetab.txt" or die "Cannot write to $output.genetab.txt";
	print $fgen join("\t", "GeneID", @gfields, @header), "\n";
	foreach my $gid (sort keys %caf) {
		my @output;
		foreach my $sv (@groupxclass) {
			my ($sgroup, $vclass) = @$sv;
			push @output, @{$obsctbygene{$sgroup}{$gid}{$vclass}}{qw|Tot Hom cHet|};
			if (grep { $vclass eq $_ } @varclswcaf) {
				push @output, sprintf("%.2E", $expctbygene{$sgroup}{$gid}{$vclass});
			}
		}
		my @ginfo;
		if (defined $xtraginfo) {
			@ginfo = map { $xtraginfo->{$gid}{$_} // "." } @gfields;  
		}
		print $fgen join("\t", $gid, @ginfo, @output), "\n";
	}
}

if (exists $conf{SampTab}) {
	my (@setxclass, @header);
	foreach my $gset (split(',', $conf{SampTab}{GeneSet})) {
		unless (grep { $gset eq $_ } (@geneset, ".")) {
			die "Cannot find $gset in predefined gene sets";
		}
		foreach my $vclass (split(',', $conf{GeneTab}{VarClass})) {
			unless(grep { $vclass eq $_ } @varclass) {
				die "Cannot find $vclass in predefined variant class";
			}
			push @setxclass, [$gset, $vclass];
			push @header, "${gset}_${vclass}_Tot", "${gset}_${vclass}_Hom", "${gset}_${vclass}_cHet"; 
		}
	}
	# Additional sample level information
	my ($xtrasinfo, @sfields);
	if (exists $conf{SampTab}{SXref}) {
		my ($xsdat, $xsfds) = slurp_xref($conf{SampTab}{SXref}, $conf{SampTab}{SXref_Fields});
		$xtrasinfo = $xsdat;
		@sfields = @$xsfds;
		foreach my $sfield (@sfields) {
			if (grep { $sfield eq $_ } @header) {
				warn "Field $sfield already appear in the genetab output!";
			}
		}
	}
	open my $fsamp, ">$output.samptab.txt" or die "Cannot write to $output.samptab.txt";
	print $fsamp join("\t", "IID", @sfields, @header), "\n";
	foreach my $iid (sort keys %sampgroup) {
		my @output;
		foreach my $sv (@setxclass) {
			my ($gset, $vclass) = @$sv;
			push @output, map { $obsctbysamp{$iid}{$gset}{$vclass}{$_} // 0 } qw|Tot Hom cHet|;
		}
		my @sinfo;
		if (defined $xtrasinfo) {
			@sinfo = map { $xtrasinfo->{$iid}{$_} // "." } @sfields;  
		}
		print $fsamp join("\t", $iid, @sinfo, @output), "\n";
	}
}



