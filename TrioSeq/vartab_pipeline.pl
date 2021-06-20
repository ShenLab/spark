#!/usr/bin/env perl

use strict;
use warnings;
use Carp;
use FindBin qw|$Bin|;
use POSIX qw|ceil|;
use Digest::MD5;
use Getopt::Lucid qw|:all|;
use IO::Prompt;
use Data::Dumper;
use Config::Std;
use File::Basename;
use Cwd qw|abs_path|;
use List::MoreUtils qw|all any uniq|;
use String::ShellQuote;
use Digest::MD5 qw|md5_hex|;
use File::Path qw|make_path|;
use Utils::Workflow;
use Utils::File qw|count_line|;
use Utils::Hash qw|merge_conf|;
use Genome::Ranges::IntSet;

use lib "$Bin/../lib";
use Shared qw|parse_fstr|;

#############################
## Command line interface  ##
#############################

my @spec =  (
	Param("conf|c")->valid(sub { -r }),
	Param("list|l")->valid(sub { -r }),
	Param("ped")->valid(sub { -r }),
	Param("tab")->valid(sub { -r }),
	Param("outdir|out"),
	Param("engine|eng")->valid(sub { $_ eq 'SGE' || $_ eq 'BASH' }),
	Keypair("param|par"),
	Switch("help|h"),
	Switch("dryrun|dry"),
	Switch("force")
	);


my $opt = Getopt::Lucid->getopt(\@spec);

if ($opt->get_help) {
	print STDERR <<EOF;
Purpose:
	This is a pipeline script to run variant tabulation and annotation in parallel.

Usage: 
	vartab_pipeline.pl --conf Config --list VCFList [--tab VarTab --ped Cohort.ped] --out OutDir 

Options:
	--list: VCF file list (required). It should contain 2~3 columns described below.
	--tab: Variant table file. Each line represent a variant in an individual. 
			It overrides the variant table in config.  
	--ped: A standard 6-col PED file describing pedigree relationship of samples in VCF file.
			It overrides the PED file in config.

Notes:
	This is wrapper of variant tabulation from VCF and annotation on tabulated variants.
	It currently supports one of dnv_table and var_table for variant tabulation. The configs 
	for tabulation and annotation should be specified in the main config. Due to the differences 
	in input VCFs, specified fields that do not appear in input will have missing values in 
	output. All tabulated fields will automatically included in annotation output.  

	The input list is a list of VCF or VCF-region combination.
		List file format: VCF [REGION] LABEL
	An optional REGION can be either a chromosome name, a region specification ("Chr:Start-End"), 
	or a BED file. When REGION is specified, so only variants within regions will be tabulated.
	It is preferrable to have multiple disjoint regions for a large VCF file. LABEL is used to 
	group differnet VCF-region combinations and used as the file name prefix in output. 

Output:
	Tabulation and annotation will be performed in parallel for each VCF or VCF-region combination.
	The output from each split will be merged by group. Merged variants of each group will
	appear in the out directory named after LABEL.

EOF
	exit 1;
}


$opt->validate({requires => [qw|list conf outdir|]});


my %conf = merge_conf($opt->get_conf, $opt->get_param); 


unless($conf{PATH}{PROG} eq 'dnv_table.pl' || $conf{PATH}{PROG} eq 'var_table.pl') {
	die "Currently only support pipelining dnv_table.pl or var_table.pl";
}
$conf{PATH}{BINDIR} = shell_quote($Bin);

if ($opt->get_ped) {
	$conf{Pedigree}{File} = $opt->get_ped;
}
unless(-f $conf{Pedigree}{File}) {
	die "Cannot find pedigree file: $conf{Pedigree}{File}";
}
if ($opt->get_tab) {
	$conf{Input}{File} = $opt->get_tab;
}
unless(-f $conf{Input}{File}) {
	die "Cannot find input variants table: $conf{Input}{File}";
}


#############################
## Input files preparation ##
#############################

my $rootdir = $opt->get_outdir;

my $wkf = Utils::Workflow->new($rootdir,
	{ engine => $opt->get_engine, force => $opt->get_force });

my (@vcfbeds, %groups);
{
	tie my %inbygrp, 'Tie::IxHash';
	open my $fin, $opt->get_list or die "Cannot open input list";
	while(<$fin>) {
		my @line = split;
		my ($vcf, $region, $label);
		if (@line == 2) {
			($vcf, $label) = @line;
			push @{$inbygrp{$label}}, [$vcf];
		}
		elsif (@line == 3) {
			($vcf, $region, $label) = @line;
			push @{$inbygrp{$label}}, [$vcf, $region];
		}
		else {
			die "Incorrect number of columns: $_";
		}
		
	}
	my $ct = 0;
	my @labels = keys %inbygrp;
	foreach my $label (@labels) {
		foreach my $dat (@{$inbygrp{$label}}) {
			push @vcfbeds, $dat;
			$ct ++;
			push @{$groups{$label}}, $ct;
		}
	}
	# Write DATA files
	for(my $ii = 0; $ii < @vcfbeds; $ii ++) {
		my $jj = $ii + 1;
		open my $fout, ">$rootdir/par/DATA.$jj" or die "Cannot write to DATA.$jj";
		print $fout join("\t", @{$vcfbeds[$ii]}), "\n";
	}
	# Write GROUP files
	open my $fgrp, ">$rootdir/par/GROUPS" or die "Cannot write to GROUPS";
	foreach my $label (@labels) {
		my @index = @{$groups{$label}};
		print $fgrp join("\t", $label, $index[0], $index[$#index]), "\n";
	}
}

# Customize and write the format and anno config file
{
	read_config $conf{PATH}{FORMAT} => my %fmt;
	$fmt{Pedigree} = $conf{Pedigree};
	$fmt{Input} = $conf{Input};
	write_config %fmt, "$rootdir/par/Format.conf";
	$conf{PATH}{FORMAT} = "$rootdir/par/Format.conf";

	my @fields;
	my %stfd = ('Chrom' => 1, 'Position' => 1, 'Ref' => 1, 'Alt' => 1);
	my $site_f = parse_fstr($fmt{Output}{Site}, 1);
	push @fields, grep { not defined $stfd{$_} } values %$site_f;

	if (defined $fmt{Input}{XtraFields}) {
		my $xtra_f = parse_fstr($fmt{Input}{XtraFields}, 1);
		push @fields, values %$xtra_f;
	}

	my $geno_f = parse_fstr($fmt{Output}{Geno}, 1);
	if ($conf{PATH}{PROG} =~ /^dnv_table/) {
		foreach my $samp (qw|Offspring Father Mother|) {
			push @fields, map { "$_.$samp" } values %$geno_f;
		}
	}
	elsif ($conf{PATH}{PROG} =~ /^var_table/) {
		push @fields, values %$geno_f;
		push @fields, qw|FamMembers Relations Phenotypes|;	
		push @fields, map { "$_.FamMembers" } values %$geno_f;
	}
	else {
		die "Cannot determine genotype fields: $conf{PATH}{PROG}";
	}

	if (exists $conf{PATH}{ANNO}) {
		read_config $conf{PATH}{ANNO} => my %anno;
		$anno{Input}{Fields} = "Chrom,Position,Ref,Alt";
		$anno{Input}{HG} = $fmt{Genome}{Build};
		$anno{Input}{XtraFields} = join(',', @fields);

		if ($conf{PATH}{PROG} =~ /^dnv_table/) {
			$anno{Sample}{IDField} = 'IID,FamID,Father,Mother,Gender,Pheno';
		}
		elsif ($conf{PATH}{PROG} =~ /^var_table/) {
			$anno{Sample}{IDField} = 'IID,FamID,Gender,Pheno';
		}
		else {
			delete $anno{Sample}{IDField};
		}
		write_config %anno, "$rootdir/par/Anno.conf";
		$conf{PATH}{ANNO} = "$rootdir/par/Anno.conf",
	}
}
write_config %conf, "$rootdir/par/run.conf" unless $opt->get_dryrun;


#################################
##  Workflow initialization &  ##
##  Working directory setup    ##
#################################

$wkf->add(var_tab(), { name => 'VarFormat', nslots => scalar(@vcfbeds), 
				expect => [ map { ["wrk/input.$_", "wrk/table.$_.txt", "wrk/table.$_.missing.txt"] } 1..scalar(@vcfbeds) ] ,
				callback => \&check_tab });
my %deparg;
if (exists $conf{PATH}{ANNO}) {
	$wkf->add(var_anno(), { name => 'VarAnno', nslots => scalar(@vcfbeds), 
				depend => 'VarFormat', deparray => 1,
				expect => [ map { "wrk/anno.$_.txt" } 1..scalar(@vcfbeds) ] });
	$deparg{depend} = 'VarAnno';
}
else {
	$deparg{depend} = 'VarFormat';
}
$wkf->add(var_merge(), { name => 'VarMerge',  %deparg,
			expect => [ map { "out/$_.txt" } sort keys %groups ] });


################################
## Kickstart workflow engine  ##
################################

$wkf->inst(\%conf);
$wkf->run({ conf => $conf{$wkf->{engine}}, dryrun => $opt->get_dryrun });


############################
## Workflow components    ##
############################

# Now split input variants 
# prepare one table for each VCF-Region combination
sub var_tab {
	my $script = <<'EOF';

read VCF REGION < _PARDIR_/DATA._INDEX_

perl _PATH.BINDIR_/module/subset_var_tab.pl _PARDIR_/run.conf $VCF $REGION \
	> _WRKDIR_/input._INDEX_

perl _PATH.BINDIR_/_PATH.PROG_ --conf _PATH.FORMAT_ --vcf $VCF \
	--tab _WRKDIR_/input._INDEX_ --out _WRKDIR_/table._INDEX_.txt

EOF
}

sub check_tab {
	my @exp = @_;
	if (all { -f "$rootdir/$_" } @exp) {
		my $n_in  = count_line("$rootdir/$exp[0]");
		my $n_out = count_line("$rootdir/$exp[1]");
		my $n_mis = count_line("$rootdir/$exp[2]");
		if ($n_in <= $n_out + $n_mis - 1) {
			return 1;
		}
		else {
			return 0;
		}
	}
	else {
		return 0;
	}
}


sub var_anno {
	my $script = <<'EOF';

NL=$(cat _WRKDIR_/table._INDEX_.txt | wc -l | awk '{print $1}')

if [[ $NL == 1 ]]; then
	touch _WRKDIR_/anno._INDEX_.txt
else 
	perl _PATH.BINDIR_/../VarAnno/anno_seqvars.pl --in _WRKDIR_/table._INDEX_.txt \
		--out _WRKDIR_/anno._INDEX_.txt --conf _PATH.ANNO_ --wrkdir _TMPDIR_/anno._INDEX_
fi

EOF
}

sub var_merge {
	my $script;
	if (exists $conf{PATH}{ANNO}) {
		$script = <<'EOF';

perl _PATH.BINDIR_/../utils/merge_var_splits.pl --input _WRKDIR_/anno --outdir _OUTDIR_ --group _PARDIR_/GROUPS

EOF
	}
	else {
		$script = <<'EOF';

perl _PATH.BINDIR_/../utils/merge_var_splits.pl --input _WRKDIR_/table --outdir _OUTDIR_ --group _PARDIR_/GROUPS

EOF
	}
	return $script;
}

