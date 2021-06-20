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
use List::Util qw|min max first|;
use List::MoreUtils qw|uniq all any none notall one|;
use Getopt::Lucid qw|:all|; 
use Config::Std;
use Hash::Util qw|lock_hash_recurse|;
use Statistics::R;
use Iterator::Simple qw|ichain igrep|;
use Utils::Parser qw|sql_query|;
use Utils::Hash qw|chk_default str2hash|;
use Utils::List qw|all_pairs insert_after|;
use Utils::File::Iter qw|iter_file slurp_file|;


use lib "$Bin/../lib";
use Shared qw|parse_filters parse_tabfile parse_calibr read_geneset read_twinsibs read_twindups parse_fstr slurp_xref expand_dat|;
use Variants qw|expand_fam var_type|;


my @spec =  (
	Param("input|in|i")->valid(sub { -r }),
	Param("output|out|o"),
	Param("wrkdir|wrk"),
	Param("conf")->valid(sub { -r }),
	Switch("help|h")
	);

my $opt = Getopt::Lucid->getopt(\@spec);

if ($opt->get_help) {
	print STDERR "rarevar_burden.pl --conf CONFIG -out OUTPUT_PREFIX\n";
	print STDERR <<EOF;
Purpose:
	This is a standalone script to evaluate the burden of rare variants in different gene sets and 
	sub-groups of samples.

Usage:
	rarevar_burden.pl [--in VarTable --wrk WorkDir] --conf Config --out Prefix

Options:
	--in: Input variant table. It overrides the variant table defined in the config.
	--wrk: Working directory. Intemediate results including filtered and sorted variants, removed variants from
		   nearby variant clusters are stored in working directory. If not provided, a temp directory will be used.
		   If working directory already exists and contain intermediate results from previous run, then they will be
		   used directly to skip the time-consuming filtering-sorting steps.  
	--out: Output file name prefix. If not provided, it will use the prefix of config file.

Notes:
	Rare variants burdens are summarized as the total counts and rates per individual in each sample group for different
	classes of variants in different gene sets. If family members are available, we can further tally "familial events" for 
	each class of variants based on the carrier status in different family members (e.g. transmitted to affected or 
	unaffected). If group pairs for comparison are specified, the per-individual rate between groups will be compared 
	using Poisson test. 

	We treat related samples in the following way. Known twins and sib pairs lists can be specified in config file. 
	Each twin/sib pair adds an edge to a relationship (undirected) graph in which each individual from twin/sib pair 
	is a node. Groups of related individuals can be represented as connected components in this relationship graph. 
	When summarizing burdens, if multiple related individuals in the same connected component are in the same sample 
	group, only the first one (based on the lexical order of ID) will be kept. When tallying familial events, events 
	in MZ twins/dups should not be double counted. And if multiple twin or dups sample are listed in the "FamMembers" 
	column, only the first one will be used. Note that for comparisons between groups, related individuals in different 
	groups will not be de-duped by the above treatment. For comparing rare variant burdens between groups, samples from 
	different groups should also match for ancestry, depth of coverage, etc. It is user's responsibility to prepare 
	appropriate input data.   

	We also specially treated multiple nearby variants (defined by a distance cutoff in config) in the same individual.
	Each variant in can be represented as a node in a gene-specific graph, each nearby variant pair in the same gene adds
	an edge to the graph. Then for each connected component with more than one node in the gene-specific graph, we select
	one node with the most severe effect. The severity rank of variant classes is taken from the order they are defined 
	in the config. The rationale is if multiple near variants are located on the same haplotype of a gene, then the 
	combined effect is usually the most severe effect from all variants. When information of transmission among family 
	members can be used to phase nearby variants (see Repluse filters in the config), if two nearby variants are inferred 
	to be on different haplotypes, the edge between this pair of variants in the graph will be deleted.

	Rare variants analyses usually require sophisticated filtering. The script allows one main site/geno level filter and
	one or more family-based filters. The main filter can also include customization for different types of variants.
	To fine tune QC parameters so that summary statistics of certain class of "neutral" variants matches the null 
	expectation, the script can also be used to "calibrate" QC metrics. It requires the definition of a filter template
	and a grid of parameters. The calibration filter is applied after main and family-based filters. The calibration 
	mode will be run in parallel with each thread using one set of QC parameters. The summary from different runs 
	will be merged.    

Input/output:
	The input files to rare variant burden analyses are defined in the config. They should include one or more variant
	tables that should typically include both variant and sample level information, and one or more sample tables that 
	can be used to define sample groups. When multiple variant tables are provided as input, variants from one sample
	must exist in only one variant table. Variant class, sample groups and gene sets are defined in the config file. 

	The following files will be appear in the output:

 	 * prefix.summary.txt -- Summary of count/rate of different variant classes and the number of different familial events
 	 						 (if defined) in different sample groups and gene sets.
 	 * prefix.calibrate.txt -- Only available under calibration mode, combined summary of variant burdens and familial event
 	 						 (if defined) for different sets of QC parameters.
	 * prefix.contra.txt -- Comparison of burden between different sample groups that are stratified by gene set and 
	 						variant class if contrasting groups are defined.
	 * prefix.genetab.txt -- Table of variant counts and number of familial events (if defined) in each gene for selected
	 						 sample groups and variant classes.
	 * prefix.samptab.txt -- Table of variant counts and number of familial events (if defined) in each sample for selected 
	 						 gene sets and variant classes.
	 * prefix.rarevars.txt -- List of QC passed rare variants found in selected sample group, variant class and gene set.  
	 						  If multiple input files are provided, only fields common to all files will be kept.

	Filtered variant table will appear in the workding directory as vartab.#.txt.gz. And list of variants removed from 
	the analysis due to variant cluster will be written to varrm.txt

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
else {
	unless(defined $conf{Variant}{Table}) {
		die "Cannot find variant table spec in config file!";
	}
}


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

#### Under calibrate model, we can use parallel to run mulltiple instances 
#### A working directory will be set up to store intermediate results

if (exists $conf{Calibrate} && exists $conf{Parallel}) {
	my $option = "";
	if (defined $conf{Parallel}{jobs}) {
		$option = join(" ", map { "--$_ $conf{Parallel}{$_}" } keys %{$conf{Parallel}});
	}
	else {
		die "Must provide number of threads in config file for parallel";
	}

	my @paramgrid = parse_calibr($conf{Calibrate});
	# update keys to $conf{Calibrate}, and write to working directory
	my $ct = 1;
	foreach my $paramset (@paramgrid) {
		foreach my $key (keys %$paramset) {
			$conf{Calibrate}{$key} = $paramset->{$key};
		}
		delete $conf{Parallel};
		# Write updated config file
		write_config %conf, "$wrkdir/param_set.$ct.conf";
		$ct ++;
	}
	my $ntot = scalar(@paramgrid);
	my $command = qq[seq 1 $ntot | parallel --eta $option "perl $0 --wrkdir $wrkdir/param_set.{} --conf $wrkdir/param_set.{}.conf --out $wrkdir/param_set.{} 1>$wrkdir/param_set.{}.out 2>$wrkdir/param_set.{}.err"];
	print $command, "\n";
	system($command);

	# Then collect output
	open my $fout, ">$output.calibrate.txt" or die "Cannot write to $output.calibrate.txt";
	for(my $ii = 1; $ii <= @paramgrid; $ii ++) {
		open my $fin, "$wrkdir/param_set.$ii.calibrate.txt" or die "Cannot open $wrkdir/param_set.$ii.calibrate.txt";
		my $header = <$fin>;
		if ($ii == 1) {
			print $fout $header;
		}
		while(<$fin>) {
			my @a = split;
			$a[0] = $ii;
			print $fout join("\t", @a), "\n";
		}
	}
	exit 1;
}

lock_hash_recurse(%conf);

#
# Sample Information
#
# Keep track of sample group membership
my (@sampgroup, %sampgroup, %sampsize, %sampincl, %samprm, %samprmbygrp, %famrm, %twins, @gpairs);
# %sampgroup is hash of hash to hold group membership of each sample appeard in the table $sampgroup{IID}{Group} = 1
# %sampsize stores sample size for each sample group $sampsize{Group} = Count
# %sampincl: is hash of IDs of samples to be included
# %samprm is hash of IDs of samples to be removed from both the main table and family members
# %famrm is has of IDs of samples to be removed from family-based analysis
# %samprmbygrp: $samprmbygrp{IID}{Group} = 1
# For group comparison, we should have two groups composed of unrelated individuals
# So we parse twins and sibs list, and create group-specific removal list so that only one of related samples will be 
# counted in group comparison. This is different from de novo analysis (dnv_burden.pl).
# %twins: For family events, we should not double count events with MZ twins. So we will keep track of 
# all twin pairs, and apply family event callback to pairs and consider PASS if any one of the pair PASS.
{
	# Slurp inclusion/exclusion list
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
		$f_iid = exists $conf{Sample}{SampID} ? $conf{Sample}{SampID} : $conf{Variant}{SampID};
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
			die "Must provide full sample size when sample group is not defined";
		}
		unless($conf{Sample}{Size} =~ /^\d+$/) {
			die "Sample size for all must be a single integer!";
		}
		# Default group is All when callbacks are not specified
		@sampgroup = qw|All|;
		# Set sample size to the pre-specified 
		$sampsize{All} = $conf{Sample}{Size};
	}
		
	# Parse sample table to determine samp size for each subgroup
	# Also need to check no sample overlaps in these sample tables
	# Note: if samples are parsed from variant table, only the first appearance will be used
	# and we will check samples from different sample tables are disjoint
	# Checking large variant tables can be time-consuming, provide sample table whenever possible
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
						warn "Sample $iid has appeared in some previous variant table";
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
				if (exists $conf{Sample}{Table}) {
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

	
	# If twin/sibs pairs are specified
	# For each sample in each group, find out if other sibs are in the same group
	# Add them to the group-specific sample removal list, and adjust sample sizes for each group
	# Note: this treatment may not be valid for inherited variants if a pair of sibs are in 
	# contracting groups, because they may both inherit same variants from parents.
	if (exists $conf{Sample}{Twins} || exists $conf{Sample}{Sibs}) {
		my %twinsibs = read_twinsibs(exists $conf{Sample}{Twins} ? $conf{Sample}{Twins} : undef, 
									 exists $conf{Sample}{Sibs}  ? $conf{Sample}{Sibs}  : undef);
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

	# We will also keep track of MZ twin or dup pairs, and will use only one of them in counting family events
	# A graph is created to keep track of pairwise connections
	if (exists $conf{Sample}{Twins}) {
		#open my $fin, $conf{Sample}{Twins} or die "Cannot open twin pairs list: $conf{Sample}{Twins}";
		#while(<$fin>) {
		#	my @pair = split;
		#	unless(@pair == 2) {
		#		die "Incorrect line in twin pairs list: $_";
		#	}
		#	$twins{$pair[0]} = $pair[1];
		#}
		%twins = read_twindups($conf{Sample}{Twins});
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
}

#
# Variant table and family events parsers
#
my (@vartabs, @vartabfields, $f_iid, @varfilt, @varfiltfields, %varexcl, %famevents);
# @vartabs: an array of all variant tables defined in the config
# @vartabfiels: common variant table fields found in all tables, will be used in dumped vartab
#			It will be updated later to include additional gene level or sample level fields.
# @varfilt: an array of filters for parsing variant table(s)
# 			All filters in the array must evaluate to be true for the variants to be included.
#			The first filter will be the main site filter, all remaining are family based filter
#			All filters are optional. The main filter can have variant type specific customization.
# @varfiltfields: an array of filter fields in the *first* filter defined above (used for expansion).
# %famevents: a hash of filters to define familial events from genotypes of family members
tie %famevents, 'Tie::IxHash';
{
	if (ref $conf{Variant}{Table} eq 'ARRAY') {
		@vartabs = @{$conf{Variant}{Table}};
	}
	else {
		@vartabs = ($conf{Variant}{Table});
	}

	# Parse variant table filter and check related fields
	# We also to check there is no overlapping sample between different variant table
	# This is because when we need to sort the table and look for nearby variants in the same person.
	# The same requirement applies to rg analysis, because we need to sort table by sample ID
	# and variant positions for phasing. Note that we do not have this requirement for CNV data.
	$f_iid = $conf{Variant}{SampID};
	
	# Overall site level filter may be defined separately for SNVs and indels
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
			my ($varid, $iid) = (split)[0,1];
			if (defined $iid) {
				$varexcl{$iid,$varid} = 1;
			}
			else {
				$varexcl{$varid} = 1;
			}
		}
	}

	my @allfnames;
	#my %knownids;
	# IID => VarTable, IDs that were appeared in previous parsed variant tables
	foreach my $vartab (@vartabs) {
		unless(-f $vartab) {
			die "Cannot find variant table $vartab!";
		}
		print STDERR "Validating variant table $vartab\n";
		my ($it, $fnames) = iter_file($vartab, { fsep => qr/\t/ });
		push @allfnames, $fnames;
		# First check the existance of SampID field
		unless (grep { $_ eq $f_iid } @$fnames) {
			die "Cannot find sample ID field in variant table $vartab";
		}
		foreach my $stdfd (qw|Chrom Position Ref Alt GeneID|) {
			unless (grep { $_ eq $stdfd } @$fnames) {
				die "Cannot find standard field $stdfd in variant table $vartab";
			}
		}
		foreach my $fd (@varfiltfields) {
			unless(grep { $fd eq $_ } @$fnames) {
				die "Cannot find Filter field $fd in variant table $vartab";
			}
		}
		# Skip validating variant table if re-use existing intermediate files
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
			my $ntimes = 0;
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

	# If FamFilter or FamEvent is defined, we will parse the filter and check if related fields are available
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
			elsif ($cond eq 'One') {
				push @varfilt => sub { one { $cb->($_) } @_ };
			}
			else {
				die "Cannot recognize condition: $cond";
			}
		}

		foreach my $reqfd (grep { $_ ne 'Relation' && $_ ne 'Pheno' && $_ ne 'IID' } keys %$fields) {
			$reqfd =~ s/^_// if $reqfd =~ /^_/;
			unless(grep { $reqfd.".FamMembers" eq $_ || $reqfd eq $_ } @vartabfields) {
				die "Cannot find FamFilter field $reqfd in variant table";
			}
		}
	}

	if (exists $conf{Variant}{FamEvent}) {
		my ($eventnames, $callbacks, $fields) = 
			parse_filters($conf{Variant}{FamEvent}, $conf{Variant}{FamEvent_Filter});
		for (my $ii = 0; $ii < @$callbacks; $ii ++) {
			if ($eventnames->[$ii] eq 'Count') {
				die "Count is a reserved word and cannot be used as event name";
			}
			$famevents{$eventnames->[$ii]} = $callbacks->[$ii];
		}

		foreach my $reqfd (grep { $_ ne 'Relation' && $_ ne 'Pheno' && $_ ne 'IID' } keys %$fields) {
			$reqfd =~ s/^_// if $reqfd =~ /^_/;
			unless(grep { $reqfd.".FamMembers" eq $_ || $reqfd eq $_ } @vartabfields) {
				die "Cannot find FamFilter field $reqfd in variant table";
			}
		}
	}
}

#
# Gene lists and gene sets
#
my (%geneset, @geneset, %generm);
# %genesets is hash of hash to hold gene set membership for each gene: $geneset{Set}{GeneID} = 1
# by default each gene (not including ".") will also be assigned to an "All" set.
# %genekeep will be an optional full inclusion list of genes
# If such list is not provided, then we will only remove blacklist genes from the variant table
# %genekeep will keep track of All genes to be included in the analysis
{
	if (exists $conf{Gene}{Exclude}) {
		my @lists;
		if (ref $conf{Gene}{Exclude} eq 'ARAY') {
			@lists = @{$conf{Gene}{Exclude}};
		}
		else {
			@lists = ($conf{Gene}{Exclude});
		}
		foreach my $excllist (@lists) {
			unless(-f $excllist) {
				die "Cannot find gene exclusion list: $excllist";
			}
			open my $fin, $excllist or die "Cannot open $excllist";
			while(<$fin>) {
				my $gid = (split)[0];
				$generm{$gid} = 1;
			}
		}
	}
	push @geneset, "All";
	if (exists $conf{Gene}{Include}) {
		die "Cannot find gene inclusion list: $conf{Gene}{Include}" unless -f $conf{Gene}{Include};
		print STDERR "Reading gene inclusion list from $conf{Gene}{Include}\n";
		foreach my $gene (map { (split)[0] } slurp $conf{Gene}{Include}) {
			next if defined $generm{$gene};
			$geneset{All}{$gene} = 1;
		}			
	}
	else {
		print STDERR "Reading the full gene list from variant tables\n";
		foreach my $vartab (@vartabs) {
			my $it = iter_file($vartab, { fsep => qr/\t/ });
			while(my $dat = $it->()) {
				next if $dat->{GeneID} eq '.'; # <- '.' should not belong to any gene set.
				foreach my $gid (split(';', $dat->{GeneID})) {
					next if defined $generm{$gid};
					$geneset{All}{$gid} = 1;
				}
			}
		}
	}
	if (exists $conf{Gene}{Set}) {
		my $gs = read_geneset($conf{Gene}{Set}, $conf{Gene}{Set_Fields});
		while(my ($set, $genes) = each %$gs) {
			if ($set eq 'All') {
				die "Gene set All is reserved!";
			}
			foreach my $gid (@$genes) {
				# Restrict to genes in All set first
				# Then assigning gene set membership
				next unless defined $geneset{All}{$gid};
				$geneset{$set}{$gid} = 1;
			}
		}
		push @geneset, sort grep { $_ ne 'All' } keys %geneset;
	}
}


#
# Variant class
#
# we have separated variant class related filter from the rest of variant level filter
my (@varclass, @varclsfilt, %varclswogene, @varclsfiltfields, %varclsrank);
# @varclass is an array of variant class defined in the config in the order of appearance
# @varclsfilt is an array of callbacks corresponding to each variant class
# %varclswogene is a hash to keep variant classes that are not associated with gene
# @varclsfiltfields is an array of field names used by varclass filter, used to expand data of variant tab
# %varclsrank will be the ranks (priority) for each defined variant class, used in nearby variant cluster
{
	# Parse variant class filter and check related fields
	my ($varcls, $filters, $fields) = parse_filters($conf{Variant}{Class}, $conf{Variant}{Class_Filter});
	@varclsfiltfields = keys %$fields;

	my %known = (All => 1, Any => 1, First => 1);

	foreach my $vclass (@$varcls) {
		my ($class, $modifier) = split(':', $vclass);
		unless(defined $modifier && defined $known{$modifier}) {
			die "Cannot find or recognize modifier: $vclass";
		}
		push @varclass, $class;
		$varclsrank{$class} = scalar(@varclass);

		my $filter = shift @$filters;
		if ($modifier eq 'All') {
			push @varclsfilt, sub { all { $filter->($_) } @_ };
		}
		elsif ($modifier eq 'Any') {
			push @varclsfilt, sub { any { $filter->($_) } @_ };
		}
		elsif ($modifier eq 'First') {
			push @varclsfilt, sub {  $filter->($_[0]) };
		}
		else {
			die "Cannot recognize modifier: $modifier";
		}
	}

	foreach my $reqfd (keys %$fields) {
		unless(grep { $reqfd eq $_ } @vartabfields) {
			die "Cannot find variant class filter field $reqfd in variant table";
		}
	}

	if (exists $conf{Variant}{ClassWoGene}) {
		foreach my $vclass (split(',', $conf{Variant}{ClassWoGene})) {
			unless(grep { $vclass eq $_ } @varclass) {
				die "Cannot find class without gene $vclass in variant class definition";
			}
			$varclswogene{$vclass} = 1;
		}
	}
}

#
# Repulse filter for nearby variants
#
my @repulsefilt;
if (exists $conf{Variant}{Repulse}) {
	my @repulse;
	if (ref $conf{Variant}{Repulse} eq 'ARRAY') {
		@repulse = @{$conf{Variant}{Repulse}};
	}
	else {
		@repulse = ($conf{Variant}{Repulse});
	}
	foreach my $repulse (@repulse) {
		$repulse =~ s/^["']//; $repulse =~ s/["']$//;
		my ($cb, $tokens) = sql_query($repulse, 1);
		foreach my $tok (@$tokens) {
			if ($tok->[0] eq 'FIELD') {
				my $field = $tok->[1]; $field =~ s/\.[12]$//; $field =~ s/^_//;
				next if $field =~ /^(Relation|Pheno|IID)$/;
				unless(grep { $field.".FamMembers" eq $_ } @vartabfields) {
					die "Cannot find Repulse filter field $field in variant table";
				}
			}
		}
		push @repulsefilt, $cb;
	}
}

#
# Initial filtering and nearby variants clustering
#
# Results of initial filter will be written to a temp file sorted by chr/pos
# The purpose of using temp file is to store large data file and reduce memory usage
# We then search nearby variants from the same individual and store them in a graph
# Nearby variant filter will be disabled under Calibrate mode.
my %varrm;
# %varrm: The nearby variants that are not prioritized will be removed from analysis: 
#		  The removal list will be sample and gene-specific: $varrm{IID,VarID}{GeneID} = 1;
# Even if calibration was not done, the basic site level filter are applied here 
# and intermediate files will be stored in temp working directory $wrkdir for further use.
#my $wrkdir = tempdir( CLEANUP => 1 );
if (all { -f "$wrkdir/vartab.$_.txt.gz" } 0..$#vartabs) {
	print STDERR "Re-use filtered variant tables in existing working directory $wrkdir\n";
	if (-f "$wrkdir/varrm.txt") {
		open my $fin, "$wrkdir/varrm.txt" or die "Cannot open $wrkdir/varrm.txt";
		while(<$fin>) {
			my ($iid, $varid, $gid) = split;
			$varrm{$iid,$varid}{$gid} = 1;
		}
	}
}
else {
	my $bk;
	if (exists $conf{Variant}{Region}) {
		unless(-f $conf{Variant}{Region}) {
			die "Cannot find region file: $conf{Variant}{Region}";
		}
		$bk = Genome::UCSC::BinKeeper->new($conf{Variant}{Region}, { bed => 1});
	}

	my %known;
	# Go through each vartab, first apply filters and create temp file that sort the variant table
	for(my $ii = 0; $ii < @vartabs; $ii ++) {
		my $vartab = $vartabs[$ii];
		my ($it, $fnames) = iter_file($vartab, { fsep => qr/\t/ });
		
		my $c_iid = first { $fnames->[$_-1] eq $conf{Variant}{SampID} } 1..scalar(@$fnames);
		my $c_chr = first { $fnames->[$_-1] eq 'Chrom' } 1..scalar(@$fnames);
		my $c_pos = first { $fnames->[$_-1] eq 'Position' } 1..scalar(@$fnames);
		unless(defined $c_iid && defined $c_chr && defined $c_pos) {
			die "Cannot find sample ID, Chrom, Position from variane table $ii: $vartab";
		}

		print STDERR "Initial filtering on variant table $vartab\n";
		open my $fout, "| body sort -k $c_iid,$c_iid -k $c_chr,$c_chr -k $c_pos,${c_pos}n -T $wrkdir/tmp | gzip -c > $wrkdir/vartab.$ii.txt.gz" 
			or die "Cannot open body sort pipe to write to $wrkdir/vartab.$ii.txt.gz";
		print $fout join("\t", @$fnames), "\n";

		#my %known;
		while(my $dat = $it->()) {
			my $iid = $dat->{$conf{Variant}{SampID}};
			my $varid = join(":", @{$dat}{qw|Chrom Position Ref Alt|});
			if (defined $known{$iid,$varid}) {
				#warn "Duplicated variant found in the variant table $vartab: $iid,$varid";
				next;
			}

			# Sample inclusion/exclusion
			if (%sampincl) {
				next unless defined $sampincl{$iid};
			}
			next if defined $samprm{$iid};
			# Variant exclusion
			if (%varexcl) {
				next if defined $varexcl{$iid,$varid} || defined $varexcl{$varid};
			}
			# Region-based inclusion
			if (defined $bk) {
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
			print $fout join("\t", @{$dat}{@$fnames}), "\n";
			$known{$iid,$varid} = 1;
		}
		close $fout;
	}

	# Then store nearby variants into a graph, and keep variant information for determining variant class
	my (%genevars, %nearvars);
	# %genevars stores graphs of nearby variants for each sample
	# %nearvars stores the raw variant level information for all nearby varianrts
	if (exists $conf{Variant}{DistCutoff} && !exists $conf{Calibrate}) {
		for(my $ii = 0; $ii < @vartabs; $ii ++) {
			print STDERR "Identify and prioritize nearby variants in filtered variant table $ii\n";
			my ($it, $fnames) = iter_file("$wrkdir/vartab.$ii.txt.gz", { fsep => qr/\t/ });
			my ($prev_iid, $prev_varid, $prev_genes, $prev_dat) = ("", "", [], {});
			while(my $dat = $it->()) {
				my $iid = $dat->{$conf{Variant}{SampID}};
				# Bad samps have already been removed in previous step
				# next if defined $samprm{$iid};
				my $varid = join(":", @{$dat}{qw|Chrom Position Ref Alt|});

				my @genes = grep { $_ ne '.' } split(';', $dat->{GeneID});
				push @genes, "." if %varclswogene; # Also add a non-gene slot for variants without gene
				
				if ($iid eq $prev_iid && $dat->{Chrom} eq $prev_dat->{Chrom} &&
					$dat->{Position}-$prev_dat->{Position} <= $conf{Variant}{DistCutoff}) {
					foreach my $gid (@genes) {
						next unless grep { $_ eq $gid } @$prev_genes;
						unless (defined $genevars{$iid,$gid}) {
							$genevars{$iid,$gid} = Graph::Undirected->new();
						}
						$genevars{$iid,$gid}->add_edge($prev_varid, $varid);

						unless(defined $nearvars{$iid,$varid}) {
							$nearvars{$iid,$varid} = $dat;
						}
						unless(defined $nearvars{$iid,$prev_varid}) {
							$nearvars{$iid,$prev_varid} = $prev_dat;
						}
					}
				}
				($prev_iid, $prev_genes, $prev_varid, $prev_dat) = ($iid, \@genes, $varid, $dat);
			}
		}
		# If repulse filters are provided, we will delete edges in the graph for variant pairs that are on different haplotypes
		if (@repulsefilt) {
			foreach my $giid (keys %genevars) {
				my ($iid, $gid) = split($;, $giid);
				my $graph = $genevars{$giid};
				my $newgraph = Graph::Undirected->new();
				foreach my $cc ($graph->connected_components()) {
					ALLPAIRS:foreach my $pair ( all_pairs(@$cc) ) {
						my @fam1 = expand_fam($nearvars{$iid,$pair->[0]}, \%famrm); 
						my @fam2 = expand_fam($nearvars{$iid,$pair->[1]}, \%famrm);
						#unless(scalar(@fam1) == scalar(@fam2)) {
						#	die "Incorrect number of family data at $iid: $pair->[0] and $pair->[1]"
						#}
						my $flag;
						for(my $kk = 0; $kk < scalar(@fam1); $kk ++) {
							#unless($fam1[$kk]{IID} eq $fam2[$kk]{IID}) {
							#	die "Samples from two family expansions do not align: $fam1[$kk]{IID} ne $fam2[$kk]{IID}!"
							#}
							my $jj = first { $fam2[$_]{IID} eq  $fam1[$kk]{IID} } 0..$#fam2;
							next unless defined $jj;

							my (%ipair_1, %ipair_2);
							foreach my $key (keys %{$fam1[$kk]}) {
								$ipair_1{"$key.1"} = $fam1[$kk]{$key};
								$ipair_1{"$key.2"} = $fam2[$jj]{$key};
								$ipair_2{"$key.1"} = $fam2[$jj]{$key};
								$ipair_2{"$key.2"} = $fam1[$kk]{$key};
							}
							if ( (any { $_->(\%ipair_1) } @repulsefilt) || (any { $_->(\%ipair_2) } @repulsefilt) ) {
								#$graph->delete_edge($pair->[0], $pair->[1]);
								$flag = 1;	
							}
						}
						# Only add nearby variant pairs when no evidence suggest they are different haplotypes
						unless($flag) {
							$newgraph->add_edge($pair->[0], $pair->[1]);
						}
					}
				}
				$genevars{$giid} = $newgraph;
			}
		}

		# Now rank nearby variant clusters to dertermine variants to be removed
		while(my ($giid, $graph) = each %genevars) {
			my ($iid, $gid) = split($;, $giid);
			foreach my $cc ($graph->connected_components()) {
				my %rank;
				foreach my $varid (@$cc) {
					my $vardat = $nearvars{$iid,$varid} // do { die "Cannot find variant data for $iid,$varid" };
					my @vardata = grep { $_->{GeneID} eq $gid } 
									expand_dat($vardat, { sep => ';', optional => [@varclsfiltfields] });
					for(my $ii = 0; $ii < @varclsfilt; $ii ++) {
						if ($gid eq '.') {
							next unless defined $varclswogene{$varclass[$ii]};
						}
						else {
							next if defined $varclswogene{$varclass[$ii]};
						}
						if ($varclsfilt[$ii]->(@vardata)) {
							$rank{$varid} = $ii;
							last;
						}
					}
					# If rank not defined? Should be removed.
					unless(defined $rank{$varid}) {
						$varrm{$iid,$varid}{$gid} = 1; 
					}
				}
				my @ccsort = sort { $rank{$a} <=> $rank{$b} } grep { defined $rank{$_} } @$cc;
				if (@ccsort > 1) {
					foreach my $othervar (@ccsort[1..$#ccsort]) {
						$varrm{$iid,$othervar}{$gid} = 1; 
					}
				}
			}
		}
	}

	if (%varrm) {
		open my $fvrm, ">$wrkdir/varrm.txt" or die "Cannot write to $wrkdir/varrm.txt";
		foreach my $viid (sort keys %varrm) {
			my ($iid, $varid) = split($;, $viid);
			foreach my $gid (sort keys %{$varrm{$viid}}) {
				print $fvrm join("\t", $iid, $varid, $gid), "\n";
			}
		}
	}
}




# Tally variants given a variant table iterator, assuming variant table has been cleaned by QC filter
# and stored in a temp file
sub tallyvars {
	my ($vit, $dump) = @_;

	# If dump is required
	my ($fdump, @ginfofields, $dumpxtraginfo, @sinfofields, $dumpxtrasinfo, @outfields);
	if ($dump) {
		@outfields = @vartabfields;
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
			insert_after(\@outfields, 'GeneID', @ginfofields);
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
			insert_after(\@outfields, $f_iid, @sinfofields);
		}

		open $fdump, ">$output.rarevars.txt" or die "Cannot write to $output.rarevars.txt";
		print $fdump join("\t", @outfields), "\n";
	}

	my (%ct, %ctbygene, %ctbysamp);
	# $ct{SampGroup}{GeneSet}{VarClass} = { Count => N, Event => M, ... }
	# $ctbygene{SampGroup}{GeneID}{VarClass} = { Count => N, Event => M, ... }
	# $ctbysamp{IID}{GeneSet}{VarClass} = ...
	while(my $dat = $vit->()) {
		my $iid = $dat->{$f_iid};
		# Find out all associated sample groups
		my @sgroup = grep { !(defined $samprmbygrp{$iid} && defined $samprmbygrp{$iid}{$_}) } 
							sort keys %{$sampgroup{$iid}};
		#next unless @sgroup;

		my $varid = join(":", @{$dat}{qw|Chrom Position Ref Alt|});

		# First calculate family events based on family data.
		# Event counts will be added to %ct later
		# We will keep one of the MZ twin pairs in counting family events.
		my %events;
		if (%famevents) {
			my @famdata = expand_fam($dat, \%famrm);
			if (@famdata) {
				my %known;
				foreach my $famdat (@famdata) {
					if (defined $twins{$famdat->{IID}}) {
						next if defined $known{$twins{$famdat->{IID}}};
					}
					while(my ($name, $cb) = each %famevents) {
						$events{$name} ++ if $cb->($famdat);
					}
					if (defined $twins{$famdat->{IID}}) {
						#$known{$famdat->{IID}} = 1;
						$known{$twins{$famdat->{IID}}} = 1;
					}
					else {
						$known{$famdat->{IID}} = 1;
					}		
				}
			}
		}

		my @data = expand_dat($dat, { sep => ';', optional => ["GeneID", @varclsfiltfields] });

		# Gene-specific counts, only when ctbygene is needed in the output
		if (wantarray() && @sgroup > 0) {
			foreach my $data (grep { $_->{GeneID} ne '.' } @data) {
				next if defined $varrm{$iid,$varid}{$data->{GeneID}};
				my @vclass;
				for(my $ii = 0; $ii < @varclass; $ii ++) {
					next if $varclswogene{$varclass[$ii]};
					if ($varclsfilt[$ii]->($data)) {
						push @vclass => $varclass[$ii];
					}
				}
				next unless @vclass > 0;
				foreach my $sgroup (@sgroup) {
					foreach my $vclass (@vclass) {
						$ctbygene{$sgroup}{$data->{GeneID}}{$vclass}{Count} ++;
						foreach my $ename (sort keys %events) {
							$ctbygene{$sgroup}{$data->{GeneID}}{$vclass}{$ename} += $events{$ename};
						}
					}
				}
			}
		}

		# Before tallying gene-set specific counts, first determine variant class not associated with genes
		if (%varclswogene && !defined $varrm{$iid,$varid}{"."}) {
			for(my $ii = 0; $ii < @varclass; $ii ++) {
				next unless $varclswogene{$varclass[$ii]};
				if ($varclsfilt[$ii]->(@data)) {
					$ctbysamp{$iid}{"."}{$varclass[$ii]}{Count} ++;
					foreach my $ename (sort keys %events) {
						$ctbysamp{$iid}{"."}{$varclass[$ii]}{$ename} += $events{$ename};
					}
					foreach my $sgroup (@sgroup) {
						$ct{$sgroup}{"."}{$varclass[$ii]}{Count} ++;
						foreach my $ename (sort keys %events) {
							$ct{$sgroup}{"."}{$varclass[$ii]}{$ename} += $events{$ename};
						}
						if (defined $fdump && 
							exists $conf{Dump}{VarClass} && $varclass[$ii] eq $conf{Dump}{VarClass} &&
							exists $conf{Dump}{SampGroup} && $sgroup eq $conf{Dump}{SampGroup} &&
							exists $conf{Dump}{GeneSet} && $conf{Dump}{GeneSet} eq '.') {
							if (defined $dumpxtraginfo) {
								foreach my $gfd (@ginfofields) {
									$dat->{$gfd} = $dumpxtraginfo->{$dat->{GeneID}}{$gfd} // ".";
								}			
							}
							if (defined $dumpxtrasinfo) {
								foreach my $sfd (@sinfofields) {
									$dat->{$sfd} = $dumpxtrasinfo->{$dat->{$f_iid}}{$sfd} // ".";
								}
							}
							print $fdump join("\t", @{$dat}{@outfields}), "\n";	
						}
					}
				}
			}
		}

		# Then, we iterate through gene sets to determine set-specific variant class
		foreach my $gset (@geneset) {
			my @vclass;
			my @datainset = grep { defined $geneset{$gset}{$_->{GeneID}} &&
								  !defined $varrm{$iid,$varid}{$_->{GeneID}} } @data;
			next unless @datainset;
			for (my $ii = 0; $ii < @varclass; $ii ++) {
				next if $varclswogene{$varclass[$ii]};
				if ($varclsfilt[$ii]->(@datainset)) {
					push @vclass => $varclass[$ii];
				}
			}
			next unless @vclass > 0;

			foreach my $vclass (@vclass) {
				$ctbysamp{$iid}{$gset}{$vclass}{Count} ++;
				foreach my $ename (sort keys %events) {
					$ctbysamp{$iid}{$gset}{$vclass}{$ename} += $events{$ename};
				}
				foreach my $sgroup (@sgroup) {	
					$ct{$sgroup}{$gset}{$vclass}{Count} ++;
					foreach my $ename (sort keys %events) {
						$ct{$sgroup}{$gset}{$vclass}{$ename} += $events{$ename};
					}
					if (defined $fdump &&
						exists $conf{Dump}{VarClass} && $vclass eq $conf{Dump}{VarClass} &&
						exists $conf{Dump}{SampGroup} && $sgroup eq $conf{Dump}{SampGroup} &&
						exists $conf{Dump}{GeneSet} && $gset eq $conf{Dump}{GeneSet}) {
						if (defined $dumpxtraginfo) {
							foreach my $gfd (@ginfofields) {
								$dat->{$gfd} = $dumpxtraginfo->{$dat->{GeneID}}{$gfd} // ".";
							}			
						}
						if (defined $dumpxtrasinfo) {
							foreach my $sfd (@sinfofields) {
								$dat->{$sfd} = $dumpxtrasinfo->{$dat->{$f_iid}}{$sfd} // ".";
							}
						}
						print $fdump join("\t", @{$dat}{@outfields}), "\n";	
					}
				}
			}
		}
	}
	if (wantarray()) {
		return (\%ct, \%ctbygene, \%ctbysamp);
	}
	else {
		return \%ct;
	}
}

#
#  Tally observed variants
#
# Under calibration mode, we will store counts under different filtering parameter combinations.
# Otherwise, we will store tallied counts as well as individual genes
my (@paramgrid, @obsctbyfilt);
# @paramgrid is an array of hash of all tunable parameters
# @obsctbyfilt is an array of observed counts under different paramter combinations
my ($obsct, $obsctbygene, $obsctbysamp);
# $obsct will store the number of variants observed in different sample groups and gene sets
# and family events if defined: $obsct->{SampGroup}{GeneSet}{VarClass} = { Count => N, Event => NE, ... } 
if (exists $conf{Calibrate}) {
	print STDERR "Calibrating filtering parameters\n";
	@paramgrid = parse_calibr($conf{Calibrate});
	# Now collect observed counts under different filtering parameters
	foreach my $paramset (@paramgrid) {
		my $filter = $conf{Calibrate}{Filter};
		while(my ($vname, $value) = each %$paramset) {
			my @match = ($filter =~ /($vname)(?=\W)/g);
			unless(@match) {
				die "Cannot find variable $vname in filter expression: $filter";
			}
			else {
				if (@match > 1) {
					warn "Variable $vname appear multiple times in the filter: $filter";
				}
			}
			$filter =~ s/$vname(?=\W)/$value/g;
		}
		# Now instantiated the filter
		$filter =~ s/^['"]//; $filter =~ s/['"]$//;
		my ($filt, $tokens) = sql_query($filter, 1);
		# First check filter is valid for all vartabs
		my @iters;
		for(my $ii = 0; $ii < @vartabs; $ii ++) {
			my ($it, $fnames) = iter_file("$wrkdir/vartab.$ii.txt.gz", { fsep => qr/\t/ });
			push @iters, $it;
			foreach my $tok (@$tokens) {
				if ($tok->[0] eq 'FIELD') {
					unless(grep { $tok->[1] eq $_ } @$fnames) {
						die "Cannot find Filter field $tok->[1] in filtered vartab.$ii.txt.gz";
					}
					unless(grep { $tok->[1] eq $_ } @varfiltfields) {
						push @varfiltfields, $tok->[1];
					}
				}
			}
		}
		my $vit = igrep { $filt->($_) } ichain(@iters);
		my $ct = tallyvars($vit);		
		push @obsctbyfilt, $ct;
	}
}
else {
	print STDERR "Tallying observed variants by function, sample group and gene set\n";
	my @iters;
	for(my $ii = 0; $ii < @vartabs; $ii ++) {
		my $it = iter_file("$wrkdir/vartab.$ii.txt.gz", { fsep => qr/\t/ });
		push @iters, $it;
	}
	my $vit = ichain(@iters);
	my $dump = exists $conf{Dump} ? 1 : 0;
	if (exists $conf{GeneTab} || exists $conf{SampTab}) {
		($obsct, $obsctbygene, $obsctbysamp) = tallyvars($vit, $dump);
	}
	else {
		$obsct = tallyvars($vit, $dump);
	}
}


#
# Statistics and summary
#
if (exists $conf{Calibrate}) {
	open my $fout, ">$output.calibrate.txt" or die "Cannot write to $output.caribrate.txt";
	my @vnames = sort grep { $_ ne 'Filter' } keys %{$conf{Calibrate}};
	my @fields = qw|SampleGroup N_Samps GeneSet N_Genes VarClass Count Rate|;
	if (%famevents) {
		push @fields, keys %famevents;
	}
	print $fout join("\t", "QCParSet", @vnames, @fields), "\n";
	for(my $ii = 0; $ii < @paramgrid; $ii ++) {
		my $param = $paramgrid[$ii];
		foreach my $sgroup (@sampgroup) {
			foreach my $gset (@geneset, ".") {
				foreach my $vclass (@varclass) {
					if ($gset eq '.') {
						next unless defined $varclswogene{$vclass};
					}
					else {
						next if defined $varclswogene{$vclass};
					}
					my $ct =  $obsctbyfilt[$ii]{$sgroup}{$gset}{$vclass};
					my $rate = sprintf("%.4f", $ct->{Count}/$sampsize{$sgroup});
					my @output = (@{$param}{@vnames}, $sgroup, $sampsize{$sgroup}, $gset, 
								  $gset eq '.' ? "." : scalar(keys %{$geneset{$gset}}), $vclass, 
								  $ct->{Count}, $rate);
					if (%famevents) {
						push @output, map { $ct->{$_} // 0 } keys %famevents;
					}
					print $fout join("\t", $ii+1, @output), "\n";
				}
			}	
		}
	}
}
else {
	open my $fout, ">$output.summary.txt" or die "Cannot write to $output.summary.txt";
	my @fields = qw|SampleGroup N_Samps GeneSet N_Genes VarClass Count Rate|;
	if (%famevents) {
		push @fields, keys %famevents;
	}
	print $fout join("\t", @fields), "\n";
	foreach my $sgroup (@sampgroup) {
		foreach my $gset (@geneset, ".") {
			foreach my $vclass (@varclass) {
				if ($gset eq '.') {
					next unless defined $varclswogene{$vclass};
				}
				else {
					next if defined $varclswogene{$vclass};
				}
				my $ct = $obsct->{$sgroup}{$gset}{$vclass};
				my ($count, $rate);
				if (defined $ct) {
					$count = $ct->{Count};
					$rate = sprintf("%.4f", $count/$sampsize{$sgroup});				
				}
				else {
					$count = 0;
					$rate = 0;
				}
				my @output = ($sgroup, $sampsize{$sgroup}, $gset, $gset eq "." ? "." : scalar(keys %{$geneset{$gset}}), 
							  $vclass, $count, $rate);
				if (%famevents) {
					if (defined $ct) {
						push @output, map { $ct->{$_} // 0 } keys %famevents;
					}
					else {
						push @output, (0) x scalar(keys %famevents);
					}
				}
				print $fout join("\t", @output), "\n";
			}
		}
	}
}


if (@gpairs > 0 && !exists $conf{Calibrate}) {
	my $R = Statistics::R->new();
	open my $fout, ">$output.contra.txt" or die "Cannot write to $output.burden.txt";
	print $fout join("\t", qw|SampGroup1 N_Samps1 SampGroup2 N_Samps2 GeneSet N_Genes VarClass Count1 Rate1 Count2 Rate2 RateRatio P-value|), "\n";
	foreach my $gpair (@gpairs) {
		my @groups = @$gpair;
		foreach my $gset (@geneset, ".") {
			foreach my $vclass (@varclass) {
				if ($gset eq '.') {
					next unless $varclswogene{$vclass};
				}
				else {
					next if $varclswogene{$vclass};
				}
				my $n_obs1 = $obsct->{$groups[0]}{$gset}{$vclass}{Count} // 0;
				my $r_obs1 = $n_obs1 / $sampsize{$groups[0]};
				my $n_obs2 = $obsct->{$groups[1]}{$gset}{$vclass}{Count} // 0;
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
				print $fout join("\t", $groups[0], $sampsize{$groups[0]}, $groups[1], $sampsize{$groups[1]}, 
									$gset, $gset eq "." ? "." : scalar(keys %{$geneset{$gset}}),
									$vclass, $n_obs1, sprintf("%.4f", $r_obs1), $n_obs2, sprintf("%.4f", $r_obs2), $ratio, $pval), "\n";
			}
		}
	}
}

# Gene table will be omitted under calibration mode
if (exists $conf{GeneTab} && !exists $conf{Calibrate}) {
	my (@groupxclass, @header);
	foreach my $sgroup (split(',', $conf{GeneTab}{SampGroup})) {
		unless(grep { $sgroup eq $_ } @sampgroup) {
			die "Cannot find $sgroup in predefined sample groups";
		}
		foreach my $vclass (split(',', $conf{GeneTab}{VarClass})) {
			unless(grep { $vclass eq $_ } @varclass) {
				die "Cannot find $vclass in predefined variant class";
			}
			if (defined $varclswogene{$vclass}) {
				die "Variant class $vclass is not associated with gene";
			}
			push @groupxclass, [$sgroup, $vclass];
			push @header, "${sgroup}_${vclass}_Count";
			if (%famevents) {
				push @header, map { "${sgroup}_${vclass}_$_" } keys %famevents;
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
				warn "Gene info field $gfield already appear in the genetab output!";
			}
		}
	}
	open my $fgen, ">$output.genetab.txt" or die "Cannot write to $output.genetab.txt";
	print $fgen join("\t", "GeneID", @gfields, @header), "\n";
	foreach my $gid (sort keys %{$geneset{All}}) {
		my @output;
		foreach my $sv (@groupxclass) {
			my ($sgroup, $vclass) = @$sv;
			push @output, $obsctbygene->{$sgroup}{$gid}{$vclass}{Count} // 0;
			if (%famevents) {
				foreach my $ename (keys %famevents) {
					push @output, $obsctbygene->{$sgroup}{$gid}{$vclass}{$ename} // 0;
				}	
			}
		}
		my @ginfo;
		if (defined $xtraginfo) {
			@ginfo = map { $xtraginfo->{$gid}{$_} // "." } @gfields;  
		}
		print $fgen join("\t", $gid, @ginfo, @output), "\n";
	}
}


if (exists $conf{SampTab} && !exists $conf{Calibrate}) {
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
			push @header, "${gset}_${vclass}_Count";
			if (%famevents) {
				push @header, map { "${gset}_${vclass}_$_" } keys %famevents;
			}
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
				warn "Field $sfield already appear in the samptab output!";
			}
		}
	}
	open my $fsamp, ">$output.samptab.txt" or die "Cannot write to $output.samptab.txt";
	print $fsamp join("\t", "IID", @sfields, @header), "\n";
	foreach my $iid (sort keys %sampgroup) {
		my @output;
		foreach my $sv (@setxclass) {
			my ($gset, $vclass) = @$sv;
			push @output, $obsctbysamp->{$iid}{$gset}{$vclass}{Count} // 0;
			if (%famevents) {
				foreach my $ename (keys %famevents) {
					push @output, $obsctbysamp->{$iid}{$gset}{$vclass}{$ename} // 0;
				}
			}
		}
		my @sinfo;
		if (defined $xtrasinfo) {
			@sinfo = map { $xtrasinfo->{$iid}{$_} // "." } @sfields;  
		}
		print $fsamp join("\t", $iid, @sinfo, @output), "\n";
	}
}



