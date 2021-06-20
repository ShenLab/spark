#!/usr/bin/env perl

use strict;
use warnings;
use Carp;
use Data::Dumper;
use FindBin qw|$Bin|;
use Getopt::Lucid qw|:all|; 
use String::ShellQuote;
use Utils::Hash qw|merge_conf|;

my @spec =  (
	Param("conf")->valid(sub { -r }),
	Param("vcf|v")->valid(sub { -r $_ && -r "$_.tbi" }),
	Param("bed|b")->valid(sub { -r || /^(\w+):([\d,]+)\-([\d,]+)$/ || /^\w+$/ }),
	Param("out|o"),
	Keypair("param|par"),
	Switch("help|h")
	);

my $opt = Getopt::Lucid->getopt(\@spec);


if ($opt->get_help) {
	print STDERR <<EOF;
Purpose:
	This is a standalone utility to filter site-only VCF file.

Usages:
	var_filter.pl --conf Config --vcf Cohort.gvcf.gz --out Output [--bed Region.bed]

Notes:
	The script has similar interface as var_filter but only works on site-only VCFs. It is implemented as 
	a wrapper of tab_vcfsites.pl utility. Pedigree file is not required and can be skipped if they are 
	specified in the config file. Only site level filter is accepted. 

	Note: some public data sets like gnomAD use non-standard INFO field names that are not meet the 
	standard requirement for variable names used by the parser for site level filter. This can be resolved
	by using a site level customization module to patch the INFO field names.


EOF
	exit 1;
}

$opt->validate({ requires => [qw|conf vcf out|] });


#######################
# Parsing input files #
#######################
my %conf = merge_conf($opt->get_conf, $opt->get_param); 

if ($opt->get_bed) {
	$conf{Genome}{Region} = $opt->get_bed;
}

my $options = "--output " . $opt->get_out . " ";
if (defined $conf{Genome}{Fasta}) {
	$options .= "--seq $conf{Genome}{Fasta} ";
}
if (defined $conf{Genome}{Region}) {
	$options .= "--bed $conf{Genome}{Region} ";
}
if (defined $conf{Module}{File}) {
	$options .= "--module $conf{Module}{File} ";
}
if (defined $conf{Module}{Site}) {
	$options .= "--custom $conf{Module}{Site} ";
}
if (defined $conf{Site}{Filter}) {
	$options .= "--filter $conf{Site}{Filter} ";
}
if (defined $conf{Site}{SNV_Filter}) {
	$options .= "--filter-snv $conf{Site}{SNV_Filter} ";
}
if (defined $conf{Site}{Indel_Filter}) {
	$options .= "--filter-indel $conf{Site}{Indel_Filter} ";
}
if (defined $conf{Output}{Site}) {
	$options .= "--select $conf{Output}{Site} ";
}

my $utildir = shell_quote("$Bin/../utils");

my $cmd = qq|perl $utildir/tab_vcfsites.pl $options --vcf |.$opt->get_vcf;
print $cmd, "\n";
system($cmd);



