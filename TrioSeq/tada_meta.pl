use strict;
use warnings;
use FindBin qw|$Bin|;
use Config::Std;
use Data::Dumper;
use Perl6::Slurp;
use String::ShellQuote;
use File::Temp qw|tempdir|;
use Hash::Util qw|lock_keys|;
use Getopt::Lucid qw|:all|;
use List::MoreUtils qw|uniq all any|;
use Utils::Parser qw|sql_query|;
use Utils::Hash qw|str2hash chk_default|;
use Utils::File::Iter qw|iter_file|;

use lib "$Bin/../lib";
use Shared qw|slurp_xref read_geneset parse_fstr|;


my @spec =  (
	Param("conf|c")->valid(sub { -r }),
	Param("wrkdir|wrk|w"),
	Param("output|out|o"),
	Switch("force"), # Under force mode, previously generated files will not be over-written
	Switch("help|h")
	);


my $opt = Getopt::Lucid->getopt(\@spec);

if ($opt->get_help) {
	print STDERR <<EOF;
Purpose:
	This is standalone script to run TADA.

Usage:
	tada_meta.pl --conf CONFIG -out OUTPUT

Options:
	--wrkdir: Working directory. It will store all intermediate result files.

Notes:
	The script implements the TADA framework for meta-analyzing de novo and case-control rare variants.
	The input should include one or more tables that summarizes the counts of different classes of variants
	in each gene. They can either be de novo variants (DNVs) from trios or ultra-rare variants from case-control
	comparisons. The DNV table should also provide mutation rate for each variant class.   

	The TADA framework also requires an estimation of prior fraction of disease genes and prior mean and
	dispersion of relative risk (assumed to follow a gamma distribution) for each class of variants.
	To integrate with gene level metrics, the current implementation allows gene-set specific or even gene specific 
	priors. The priors can be fixed parameters by MLE or samples from MCMC (extTADA). If full sampling distribution 
	of priors is used, they should be stored in one external file for each gene set, and we need to be sure that all 
	parameters including prior fraction of disease genes were jonitly sampled from the same run of MCMC. 
	
	The original TADA assumes that disease genes are targeted by all classes of variants. The assumption may not
	be valid in real data, for example, if some disease genes harbor missense variants only. To allow multiple classes 
	of disease genes that are affected by different classes of variants, the two-class mixture model of TADA can be 
	extended to multi-class mixture. For multi-class TADA model, BF calculation requires more than one sets of priors
	defined for each class of variants, and prior estimation can be done by modified extTADA.

	The output will include one column of BF for each variant class, final BF synthesized from all classes of variants
	(log10BF_All), posterior probability of disease gene (PP_All), and estimated FDR (Qvalue) for each gene.


EOF

	exit 1;
}


$opt->validate({ requires => [qw|conf|] });

my $config = $opt->get_conf;

my $outfile;
if($opt->get_output) {
	$outfile = $opt->get_output;
}
else {
	$outfile = $config;
	$outfile =~ s/\.conf$//;
	$outfile .= ".txt";
}


my $module = shell_quote("$Bin/module");

read_config $config => my %conf;
lock_keys(%conf);


my $wrkdir = $opt->get_wrkdir;
if ($wrkdir) {
	mkdir $wrkdir unless -d $wrkdir;
} else {
	$wrkdir = tempdir(CLEANUP => 1);
}

#
# First slurp variant tables: GeneID, Count1, MutRate1, Count2, MutRate2, Case1, Control1...
#
my ($varcnts, $varcols) = slurp_xref($conf{Input}{Table}, $conf{Input}{Fields});
# Validate the input tables fields
if (@$varcols % 2) {
	die "Incorrect number columns for variant counts and mutation rates!";
}
if (any { /\./ } @$varcols) {
	die "Incorrect column names for variant counts and mutation rates (should not contain dot)!"
}
if (scalar(@$varcols) != scalar(uniq sort @$varcols)) {
	die "Column name for variant counts and mutation rates are not unique!";
}

my (@tables, @varfields);
if (ref $conf{Input}{Table} eq 'ARRAY') {
	@tables = @{$conf{Input}{Table}};
	foreach my $fstr (@{$conf{Input}{Fields}}) {
		my $vfields = parse_fstr($fstr, 1);
		push @varfields, $vfields;
	}
}
else {
	@tables = ($conf{Input}{Table});
	my $vfields = parse_fstr($conf{Input}{Fields}, 1);
	@varfields = ($vfields);
}

# Optional: associate each variant table with a label for grouping
# It will be used later to sort priors into table
my @labs;
my %ncolbylab;
if (exists $conf{Input}{Label}) {
	if (ref $conf{Input}{Label} eq 'ARRAY') {
		@labs = @{$conf{Input}{Label}};
	}
	else {
		@labs = ($conf{Input}{Label});
	}
	unless(@labs == @tables) {
		die "If label is provided, it should be one for each table!";
	}
	# Check each label should be associatefd with the same number of variant columns
	for(my $ii = 0; $ii < @labs; $ii ++) {
		my $ncol = scalar(keys %{$varfields[$ii]})-1;
		chk_default(\%ncolbylab, $labs[$ii], $ncol);
	}
}


# Determine sample sizes for each variant class 
# For de novo data, it's number of families; for case-control, it's the sampe sizes for cases and controls 
my %sampsize;
{
	my @sz;
	if (ref $conf{Input}{SampSize} eq 'ARRAY') {
		@sz = @{$conf{Input}{SampSize}};
	}
	else {
		@sz = ($conf{Input}{SampSize});
	}
	unless(@sz == @tables) {
		die "Sample size should be specified for each table!";
	}

	for(my $ii = 0; $ii < @sz; $ii ++) {
		#my $vfields = parse_fstr($fields[$ii], 1);
		my $vfields = $varfields[$ii];
		my @vcols = values(%$vfields);
		unless (@vcols % 2) {
			die "Incorrect number of columns for variant table $ii";
		}
		if ($sz[$ii] =~ /,/) {
			# Case-control
			$sz[$ii] = [split(',', $sz[$ii])];
			unless(@{$sz[$ii]} == 2) {
				die "Incorrect sample size for variant table $ii";
			}
		}
		for(my $jj = 1; $jj <= (@vcols-1)/2; $jj ++) {
			if (ref $sz[$ii] eq 'ARRAY') {
				$sampsize{$vcols[$jj*2-1]} = $sz[$ii][0];
				$sampsize{$vcols[$jj*2]}   = $sz[$ii][1];
			}
			else {
				$sampsize{$vcols[$jj*2-1]} = $sz[$ii];
			}
		}
	}
}

# See if we need to remove some gene
my %generm;
if (defined $conf{Gene}{Exclude}) {
	%generm = map { (split)[0] => 1 } slurp $conf{Gene}{Exclude};
}
# The final set of "All" genes include genes appear in the variant counts table 
# removing those in the exclusion list
my @allgenes = grep { !defined $generm{$_} } sort keys %{$varcnts};
print scalar(@allgenes), "\n";

#
# Then slurp gene info: GeneID, GeneCol1, GeneCol2, ...
#
my ($geneinfo, $genecols);
if (exists $conf{Gene} && defined $conf{Gene}{GXref}) {
	($geneinfo, $genecols) = slurp_xref($conf{Gene}{GXref}, $conf{Gene}{GXref_Fields});
	unless(all { defined $geneinfo->{$_} } @allgenes) {
		my @gmiss = grep { !defined $geneinfo->{$_} } @allgenes;
		#print STDERR "Cannot find info for the following genes\n";
		#print STDERR join("\n", @gmiss), "\n";
		if ($conf{Gene}{Strict} =~ /^Y|T/i) {
			die "Gene info cannot be found for all genes!";
			exit 1;
		}
	}
	# Fill in missing data for gene info
	foreach my $gid (@allgenes) {
		foreach my $gcol (@$genecols) {
			$geneinfo->{$gid}{$gcol} = '.' unless defined $geneinfo->{$gid}{$gcol};
		}
	}
}

# Output fields
# GeneID, Gene level information, VarCount/MutRate, BF1, Case/Control, BF2, ... BF_all PP_all
my @outfields = qw|GeneID|;
push @outfields, @$genecols if defined $genecols;
for(my $ii = 0; $ii < scalar(@$varcols)/2; $ii ++) {
	push @outfields, $varcols->[2*$ii];
	push @outfields, $varcols->[2*$ii+1];
	push @outfields, "log10BF_$varcols->[2*$ii]";
}
push @outfields, "log10BF_All", "PP_All", "Qvalue";
unless(@outfields == scalar(uniq sort @outfields)) {
	print STDERR join(",", @outfields), "\n";
	die "Output file contains non-unique field names";
}


#
# Determine gene sets if gene set filter is provided. Otherwise, all genes belongs to "All" set.
#
# Each gene should be uniquely associated with one set!
my %gs;
if (exists $conf{Gene} && defined $conf{Gene}{SetFilter}) {
	my %callbacks;
	foreach my $filter (split(',', $conf{Gene}{SetFilter})) {
		my ($expr, $label) = split(':', $filter);
		unless(defined $label) {
			die "Cannot find gene set label: $filter";
		}
		else {
			if ($label eq 'All') {
				die "Gene set label 'All' is reserved!";
			}
			if(defined $callbacks{$label}) {
				die "Gene set label has been used before: $label"	
			}
		}
		$expr =~ s/^["']//; $expr =~ s/["']$//;
		my ($cb, $tokens) = sql_query($expr, 1);
		foreach my $tok (@$tokens) {
			if ($tok->[0] eq 'FIELD') {
				unless(grep { $_ eq  $tok->[1] } @$genecols) {
					die "Cannot find field from gene info files: $tok->[1]";
				}
			}
		}
		$callbacks{$label} = $cb;
	}
	# Classify genes into gene set
	foreach my $gid (@allgenes) {
		my @set;
		while(my ($label, $cb) = each %callbacks) {
			if ($cb->($geneinfo->{$gid})) {
				push @set, $label;
			}
		}
		unless(@set == 1) {
			die "Gene $gid cannot be uniquely classified to one gene set";
		}
		push @{$gs{$set[0]}}, $gid;
	}
}
else {
	push @{$gs{All}}, @allgenes;
}

#
# Slurp piror if priors are provided, priors will be written to gene-set specific tables
#
# Prior are specified for each gene class, and optionally for each sample group
my %priors;
if (exists $conf{Prior}) {
	unless(all { defined $conf{Prior}{$_} } sort keys %gs) {
		die "Cannot find priors for all gene sets!";
	}
	foreach my $geneset (sort keys %gs) {
		if (ref $conf{Prior}{$geneset} eq 'ARRAY') {
			# When multiple prior sets are specified, they should be provided as files 
			# or they should be packed together.
			if (any { -f $_ } @{$conf{Prior}{$geneset}}) {
				unless(all { -f $_ } @{$conf{Prior}{$geneset}}) {
					die "Not all prior table files can be found!";
				}
				$priors{$geneset} = $conf{Prior}{$geneset};
			}
			else {
				my @params = map { my $par = str2hash($_, { psep => ',', kvsep => ':', order => 1 }); $par  } @{$conf{Prior}{$geneset}};
				unless(all { scalar(keys %$_) == @$varcols + 1 } @params) {
					die "Incorrect number of prior parameters for $geneset";
				}
				for(my $ii = 0; $ii < @params; $ii ++) {
					my $jj = $ii + 1;
					my $param = $params[$ii];
					open my $fpar, ">$wrkdir/Prior_$geneset.$jj.txt" or die "Cannot write to $wrkdir/Prior_$geneset.$jj.txt";
					print $fpar join("\t", keys %$param), "\n";
					print $fpar join("\t", values %$param), "\n";
					push @{$priors{$geneset}} => "$wrkdir/Prior_$geneset.$jj.txt";
				}
			}
		}
		else {
			# If priors are provided as a file
			if (-f $conf{Prior}{$geneset}) {
				# Validate its column numbers
				my ($it, $fnames) = iter_file($conf{Prior}{$geneset});
				unless (@$fnames == @$varcols + 1 || @$fnames == @$varcols + 2) {
					die "Incorrect number of columns for $conf{Prior}{$geneset}";
				}
				$priors{$geneset} = $conf{Prior}{$geneset};
			}
			else {
				my $param = str2hash($conf{Prior}{$geneset}, { psep => ',', kvsep => ':', order => 1 });
				if (scalar(keys %$param) == 1) {
					my @keys = keys %$param;
					my @vals = values %$param;
					unless (any { exists $conf{Prior}{"$_.$geneset"} } @labs) {
						die "Cannot find prior for any table label";
					} 
					for(my $ii = 0; $ii < @labs; $ii ++) {
						if (exists $conf{Prior}{"$labs[$ii].$geneset"}) {
							my $setpar = str2hash($conf{Prior}{"$labs[$ii].$geneset"}, { psep => ',', kvsep => ':', order => 1 });
							unless(scalar(keys %$setpar) == $ncolbylab{$labs[$ii]}) {
							#print STDERR scalar(keys %$setpar), " != ", $ncolbylab{$labs[$ii]}, "\n";
							#print Dumper $setpar; 
								die "Number of prior params does not match variant class in $labs[$ii]";
							}
							push @keys, keys %$setpar;
							push @vals, values %$setpar;
						}
						else {
							warn "No prior was specified for $labs[$ii].$geneset";
							push @keys, ('.') x $ncolbylab{$labs[$ii]};
							push @vals, ('.') x $ncolbylab{$labs[$ii]};
						}
					}
					open my $fpar, ">$wrkdir/Prior_$geneset.txt" or die "Cannot write to $wrkdir/Prior_$geneset.txt";
					print $fpar join("\t", @keys), "\n";
					print $fpar join("\t", @vals), "\n";	
				}
				else {
					unless(scalar(keys %$param) == @$varcols + 1) {
						die "Incorrect number of prior parameters for $geneset";
					}
					open my $fpar, ">$wrkdir/Prior_$geneset.txt" or die "Cannot write to $wrkdir/Prior_$geneset.txt";
					print $fpar join("\t", keys %$param), "\n";
					print $fpar join("\t", values %$param), "\n";	
				}
				$priors{$geneset} = "$wrkdir/Prior_$geneset.txt";
			}
		}
	}
}

#
# Calculate BF/PP
#
# BF/PP will be calculated for each gene set defined, and results are then collected and combined.
my %BFPP;
foreach my $geneset (sort keys %gs) {
	open my $fvar, ">$wrkdir/Input_$geneset.txt" or die "Cannot write to $wrkdir/Input_$geneset.txt";
	print $fvar join("\t", "GeneID", map { defined $sampsize{$_} ? "$_.$sampsize{$_}" : $_ } @$varcols), "\n";
	foreach my $gid (@{$gs{$geneset}}) {
		print $fvar join("\t", $gid, map { $varcnts->{$gid}{$_} // '.' } @$varcols), "\n";
	}
	close $fvar;
	next unless defined $priors{$geneset};
	my $prioropt;
	if (ref $priors{$geneset} eq 'ARRAY') {
		$prioropt = join(" ", @{$priors{$geneset}});
	}
	else {
		$prioropt = $priors{$geneset};
	}
	print(qq|Rscript $module/calc_TADA_BFPP.R $wrkdir/Input_$geneset.txt $prioropt $wrkdir/Output_$geneset.txt|, "\n");
	system(qq|Rscript $module/calc_TADA_BFPP.R $wrkdir/Input_$geneset.txt $prioropt $wrkdir/Output_$geneset.txt|);
	my $it = iter_file("$wrkdir/Output_$geneset.txt", { fsep => qr/\t/ });
	while(my $dat = $it->()) {
		if (defined $BFPP{$dat->{GeneID}}) {
			die "Results for gene $dat->{GeneID} has been collected!";
		}
		$BFPP{$dat->{GeneID}} = $dat;
	}
}

unless(exists $conf{Prior}) {
	print STDERR "No prior is provided for TADA BF/PP calculation, exit!\n";
	exit 1;
}

# 
# Reformat final output
# 
open my $fout, ">$outfile" or die "Cannot write to $outfile";
print $fout join("\t", @outfields), "\n";

# Sort genes based on PP then by BF
# Calculate FDR on the fly
my ($cumsumq0, $nn) = (0, 0);
foreach my $gid (sort { $BFPP{$b}{PP_All} <=> $BFPP{$a}{PP_All} ||
				  	    $BFPP{$b}{log10BF_All} <=> $BFPP{$a}{log10BF_All} } keys %BFPP) {
	my @outvals = ($gid);
	if (defined $genecols) {
		push @outvals, map { $geneinfo->{$gid}{$_} // "." } @$genecols;
	}
	for(my $ii = 0; $ii < scalar(@$varcols)/2; $ii ++) {
		push @outvals, $varcnts->{$gid}{$varcols->[2*$ii]} // ".";
		
		if (defined $sampsize{$varcols->[2*$ii+1]}) {
			push @outvals, $varcnts->{$gid}{$varcols->[2*$ii+1]} // ".";
		}
		else {
			my $mutrate = $varcnts->{$gid}{$varcols->[2*$ii+1]};
			if (defined $mutrate) {
				push @outvals, sprintf("%.2E", $mutrate);
			}
			else {
				push @outvals, ".";
			}
		}

		my $BF = $BFPP{$gid}{"log10BF_$varcols->[2*$ii]"};
		if ($BF =~ /^\d+$/) {
			push @outvals, $BF;
		}
		else {
			push @outvals, sprintf("%.2f", $BF);
		}
	}

	push @outvals, sprintf("%.2f", $BFPP{$gid}{"log10BF_All"});
	push @outvals, sprintf("%.3f", $BFPP{$gid}{"PP_All"});

	$cumsumq0 += 1-$BFPP{$gid}{"PP_All"};
	$nn ++;

	my $qval = $cumsumq0/$nn;
	print $fout join("\t", @outvals, $qval < 0.01 ? sprintf("%.2E", $qval) : sprintf("%.2f", $qval)), "\n";
}


