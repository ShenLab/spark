#!/usr/bin/env perl

use strict;
use warnings;
use Carp;
use IO::Dir;
use IO::File;
use Perl6::Slurp;
use Data::Dumper;
use File::Temp qw|tempdir|;
use FindBin qw|$Bin|;
use List::Util qw|min max first|;
use List::MoreUtils qw|all any uniq none|;
use Getopt::Lucid qw|:all|; 
use Config::Std;
use Hash::Util qw|lock_hash_recurse|;
use Utils::List qw|all_pairs|;
use Utils::Parser qw|sql_query|;
use Utils::File::Iter qw|iter_file|;


use lib "$Bin/../lib";
use Shared qw|parse_filters parse_tabfile read_geneset parse_fstr slurp_xref expand_dat|;
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
	This is a standalone script to calculate cumulative allele frequencies per gene for recessive burden analysis.

Usage:
	rg_caf.pl [--in VarTable] --conf Config --out Output 

Options:
	--in: Input variant table with inidividual genotypes. It overrides the variant table defined in config.
	--wrkdir: Working directory. Intermediate filter-sorted variants are stored in this directory. If not provided
			  a temp directory will be used. If working directory already exists and contain intermediate table
			  from previous run, it will be used directly to skip the time-consuming filtering-sorting steps. 

Notes:
	Cumulative allele frequency (CAF) for a class of variants in each gene is defined as the proportion of 
	haplotypes in the population that carry that class of variants. CAFs are typically calculated on the parents
	cohort of affected individuals. Parents are presumed to be unrelated and should have the same ancestry 
	background. The resulting CAFs will be used by rg_burden.pl to calculate expected number of RGs observed in
	offspring for each gene. The output file will have one row per gene, and one column for CAF of each varinat 
	class.

EOF
	exit 1;
}

$opt->validate({ requires => [qw|conf output|] });


my $infile = $opt->get_input;
my $outfile = $opt->get_output;
my $config = $opt->get_conf;

if (-f $outfile) {
	unless($opt->get_force) {
		print STDERR "Output file $outfile already exists!\n";
		exit 1;
	}
}

read_config $opt->get_conf => my %conf;

if (defined $infile) {
	unless (-f $infile || -d $infile) {
		die "Cannot find input file or directory: $infile";
	}
	if (-f $infile) {
		$conf{Variant}{Table} = $infile;
	}
}
else {
	unless(defined $conf{Variant}{Table} && -f $conf{Variant}{Table}) {
		die "Must provide input variant table!";
	}
}

lock_hash_recurse(%conf);

my $wrkdir = $opt->get_wrkdir;
if ($wrkdir) {
	mkdir $wrkdir unless -d $wrkdir;
} else {
	$wrkdir = tempdir(CLEANUP => 1);
}


#
# Sample inclusion/exclusion
# 
my ($sampsize, %sampincl, %samprm, %famrm);
{
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
	# Determine sample size from sample table, if sample table is provided, we will
	# use sample IDs that appear in the variant table
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

	my (%known, %knowntabs);
	# When sample table is provided, we will check to ensure that each sample appear in only one table.
	# Note: for RG analysis, we also require that one sample can only exist in one table
	#		although in our implementation, we have merged sorted all data tables into one
	foreach my $samptab (@samptabs) {
		my ($it, $fnames) = iter_file($samptab, { fsep => qr/\t/ });
		unless(grep { $f_iid eq $_ } @$fnames) {
			die "Cannot find sample ID field $f_iid from sample table";
		}
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
						die "Sample $iid has appeared in some previous variant table";
					}
					next;
				}
			}
			$known{$iid} = 1;
			$knowntabs{$iid} = $samptab;
			$sampsize ++;
		}
	}

	print STDERR "A total of $sampsize samples are found in sample table.\n";

	if (exists $conf{Sample}{Size}) {
		if (exists $conf{Sample}{Table}) {
			unless ($conf{Sample}{Size} == $sampsize) {
				die "Inconsistent sample size, $sampsize found in table <> $conf{Sample}{Size} specified in config";
			}
		}
		else {
			if($conf{Sample}{Size} < $sampsize) {
				die "Num of samples ($sampsize) found in variant table < $conf{Sample}{Size} specified in config";
			}
			elsif ($conf{Sample}{Size} > $sampsize) {
				warn "Num of samples ($sampsize) found in variant table > $conf{Sample}{Size} specified in config";	
			}
			else {
				print "Num of samples ($sampsize) found in variant table = $conf{Sample}{Size} specified in config";
			}
		}
	}
}

#
# Variant table and phase filters
#
# If over-written by the command line option, only one variant table is allowed
# If specified in config, it is possible to have multiple vartabs (no parallelization)
# Alternatiely, we can provide a directory of multiple vartabs and do parallelization
my (@vartabs, @vartabfields, $f_iid, @varfilt, @varfiltfields, %varexcl, $hetfilt, $homfilt, @genofiltfields, @repulsefilt);
# @vartabs: One or more variant tables
# @vartabfields: common variant table fields found in all tables
#
# @varfilt: an array of filters for parsing variant table
#			All filters in the array must be evaluate to true for variant to be included
# 			The first filter will be main filter, all remaining are optional family-based filter
# 			The main filter can have variant type specific customization
# @varfiltfields: an array of fields in the main filter defined above (used for exapnsion)			
# 
# $hetfilt : determine if a genotype is heterozygotes or not
# $homfilt : similarly for homozygotes
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
	my %knownids; # IID => VarTab
	foreach my $vartab (@vartabs) {
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

		if (@vartabs > 1) {
			my %sampids;
			while(my $dat = $it->()) {
				my $iid = $dat->{$f_iid};
				if (defined $knownids{$iid}) {
					die "Sample $iid in $vartab appear in table $knownids{$iid}";
				}
				$sampids{$iid} = 1;
			}
			foreach my $iid (keys %sampids) {
				$knownids{$iid} = $vartab;
			}
		}
	}

	# Keep variant table fields that are common to all input tables.
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
			unless(grep { $reqfd.".FamMembers" eq $_ } @vartabfields) {
				die "Cannot find FamFilter field $reqfd in variant table(s)";
			}
		}
	}

	# Now parse repulse filter, filter fields will have .1/.2 suffix added to existing family member fields
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
						die "Cannot find Repulse filter field $field in variant table(s)";
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
# Note: we do not filter nearby variants as in rarevar_burden, they will be collapsed into CAF
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
# Initial filtering and variant table sorting
#
# A temp table of filtered variants will be generated. 
# The site level and family based filter will be applied to input variant table(s).
# If duplicated variant-IID combination exists, only the first one will be kept.
# In case a variant was annotated to multiple genes, they will be *EXPANDED* to one row per gene.
# Then we sort the table by GeneID then by IID, so all variants within an individual will appear in
# consecutive rows.
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
	open my $fout, "| body sort -k $c_gid,$c_gid -k $c_iid,$c_iid -k $c_chr,$c_chr -k $c_pos,${c_pos}n | gzip -c > $wrkdir/vartab.txt.gz"
			or die "Cannot write open body sort pipe to write to $wrkdir/vartab.txt.gz";
	print $fout join("\t", @vartabfields), "\n";

	for(my $ii = 0; $ii < @vartabs; $ii ++) {
		my $vartab = $vartabs[$ii];
		my $it = iter_file($vartab, { fsep => qr/\t/ });
		print "Initial filtering on variant table $vartab\n"; 
		
		my %known;
		while(my $dat = $it->()) {
			my $iid = $dat->{$f_iid};
			# Autosome variants filter will be enabled by default!
			next unless $dat->{Chrom} =~ /^\d+$/ || $dat->{Chrom} =~ /^chr\d+$/; 

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
			# We will expand the data if multiple gene IDs are present
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
# Go-through filter-merged variant table to calculate CAFs
#
my %caf;
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
				$caf{$prev_gid}{$varclass[$ii]} += 2;
			}
			else {
				if (@clsHetVars == 1) {
					$caf{$prev_gid}{$varclass[$ii]} += 1;
				} 
				elsif (@clsHetVars > 1) {
					# If more than one heterozygotes are observed, we need resolve their phases.
					# If any pair of them are on different haplotypes (cHet), it contribute 2 copies to CAF.
					my $cHetFlag;
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
								$caf{$prev_gid}{$varclass[$ii]} += 2;
								$cHetFlag = 1;
								last ALLPAIRS;
							}
						}
					}
					unless ($cHetFlag) {
						$caf{$prev_gid}{$varclass[$ii]} += 1;
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
			$caf{$prev_gid}{$varclass[$ii]} += 2;
		}
		else {
			if (@clsHetVars == 1) {
				$caf{$prev_gid}{$varclass[$ii]} += 1;
			} 
			elsif (@clsHetVars > 1) {
				my $cHetFlag;
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
							$caf{$prev_gid}{$varclass[$ii]} += 2;
							$cHetFlag = 1;
							last ALLPAIRS;
						}
					}
				}
				unless ($cHetFlag) {
					$caf{$prev_gid}{$varclass[$ii]} += 1;
				}
			}
		}
	}
}


my @header = map { exists $conf{Output}{FieldPref} ? $conf{Output}{FieldPref}.$_ : $_ } @varclass;

my ($xtraginfo, @gfields);
if (defined $conf{Output}{GXref}) {
	my ($xgdat, $xgfds) = slurp_xref($conf{Output}{GXref}, $conf{Output}{GXref_Fields});
	$xtraginfo = $xgdat;
	@gfields = @$xgfds;
	foreach my $gfield (@gfields) {
		if (grep { $gfield eq $_ } @header) {
			die "Field $gfield already appear in the genetab output!";
		}
	}
}

open my $fout, ">$outfile" or die "Cannot write to $outfile";
print $fout join("\t", "GeneID", @gfields, @header), "\n";

foreach my $gid (sort keys %caf) {
	my @ginfo;
	if (defined $xtraginfo) {
		@ginfo = map { $xtraginfo->{$gid}{$_} // "." } @gfields;  
	}
	my @allcafs;
	foreach my $vclass (@varclass) {
		if (defined $caf{$gid}{$vclass}) {
			push @allcafs, $caf{$gid}{$vclass}/(2*$sampsize);
		}
		else {
			push @allcafs, 0;
		}
	}
	print $fout join("\t", $gid, @ginfo, @allcafs), "\n";
}




