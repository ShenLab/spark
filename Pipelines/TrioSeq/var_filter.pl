#!/usr/bin/env perl

use strict;
use warnings;
use Carp;
use IO::File;
use Data::Dumper;
use FindBin qw|$Bin|;
use List::Util qw|sum|;
use List::MoreUtils qw|all any uniq|;
use Getopt::Lucid qw|:all|; 
use Config::Std;
use Perl6::Slurp;
use Module::Load;
use FaSlice;
use Genet::Var qw|normalize|;
use Genet::File::VCF;
use Genome::UCSC qw|hg_par|;
use Genome::UCSC::TwoBit;
use Utils::Parser qw|sql_query|;
use Utils::Hash qw|merge_conf chk_default|;
use Utils::File::Iter qw|iter_file|;

use lib "$Bin/../lib";
use Shared qw|parse_fstr parse_tabfile fam_rels|;
use Variants qw|var_type get_vcfped_samps get_fieldsinfo expand_site flatten_geno arrange_genodat|;


my @spec =  (
	Param("conf")->valid(sub { -r }),
	Param("vcf|v")->valid(sub { -r $_ && -r "$_.tbi" }),
	Param("ped|p")->valid(sub { -r }),
	Param("bed|b")->valid(sub { -r || /^(\w+):([\d,]+)\-([\d,]+)$/ || /^\w+$/ }),
	Param("out|o"),
	Switch("strict"),
	Keypair("param|par"),
	Switch("help|h")
	);

my $opt = Getopt::Lucid->getopt(\@spec);


if ($opt->get_help) {
	print STDERR <<EOF;
Purpose:
	This is standalone script to filter and tabulate variants from joint genotype calls.

Usage:
	var_filter.pl --conf Config --vcf Cohort.gvcf.gz [--ped Cohort.ped --bed Region.bed] --out Output

Options:
	--ped: A standard 6-col PED file describing pedigree relationship of samples in VCF file.
			It overrides the PED file in config.
	--bed: Genomic regions in 3-col BED file format. Only variants in these regions will be filtered.
			It override the BED file in config.
	--strict: Under strict mode, all fields used in filtering or output must have a defined type. 

Note:
	This is a general purpose utility for variant filtering from VCF. It makes use of existing and 
	customized fields from VCF file in site and genotype level filters. PED file must be provided
	even if all samples are unrelated. Only samples appearing in both VCF and PED will be processed
	and used in customizing VCF fields. 

	For family-based sequencing data, it is recommended to provide the most complete PED file such 
	that relationship between each pair of family members can be inferred from PED file. To tabulate 
	sample genotypes, a sample inclusion list must be provided in the config. Otherwise, only site 
	level filter will be applied and site level information will appear in the output.

Output:
	The output is a table of variants passed filter with columns of site and genotype level information
	specified in the config file. When family members are available, genotype related information for 
	available family members will be packed in one column for each field with FamMembers suffix. 

	Information from duplicated samples or MZ twins is used in filtering and output. 

EOF
	exit 1;
}

$opt->validate({ requires => [qw|conf vcf out|] });


#######################
# Parsing input files #
#######################
my %conf = merge_conf($opt->get_conf, $opt->get_param); 
if ($opt->get_ped) {
	$conf{Pedigree}{File} = $opt->get_ped;
}
unless(-f $conf{Pedigree}{File}) {
	croak "Cannot find pedigree file: $conf{Pedigree}{File}";
}

if ($opt->get_bed) {
	$conf{Genome}{Region} = $opt->get_bed;
}

my $seq;
if ($conf{Genome}{Fasta} =~ /\.2bit$/) {
	$seq = Genome::UCSC::TwoBit->new($conf{Genome}{Fasta});
}
else {
	$seq = FaSlice->new(file => $conf{Genome}{Fasta});
}
my $exclchr = {};
if (defined $conf{Genome}{IgnoreChr}) {
	$exclchr = parse_fstr($conf{Genome}{IgnoreChr});
}


my $vcf;
{
	# Get a VCF handler for creating VCF object
	my $fvcf;
	if (exists $conf{Genome}{Region}) {
		my $region = $conf{Genome}{Region};
		if (-f $region) {
			open $fvcf, "tabix -p vcf -h -R " . $region ." ". $opt->get_vcf . " |" 
				or die "Cannot open tabix pipe";
		}
		elsif ($region =~ /^(\w+):([\d,]+)\-([\d,]+)$/ || $region =~ /^(\w+)$/) {
			$region =~ s/,//g;
			open $fvcf, "tabix -p vcf -h " . $opt->get_vcf . " $region |"
				or die "Cannot open tabix pipe";
		}
		else {
			die "Cannot find region file or recognize region spec: $region";
		}
	}
	else {
		$fvcf = $opt->get_vcf;
	}
	$vcf = Genet::File::VCF->new($fvcf);
}

# Determine if we need to dynamically load custom module
if (exists $conf{Module}{File}) {
	unless (-f $conf{Module}{File}) {
		die "Cannot find module file: $conf{Module}{File}";
	}
	else {
		load $conf{Module}{File};
	}
}

my $samps = get_vcfped_samps($opt->get_vcf, $conf{Pedigree}{File}, 
	{ ignore => $conf{Pedigree}{Ignore}, trios => 0, twins => $conf{Pedigree}{Twins} });
unless(keys %$samps) {
	print STDERR "No sample was found in VCF/PED files\n";
	exit 1;
}
else {
	print STDERR "Found a total of ", scalar(keys %$samps), " samples in VCF/PED files\n";
}

# Determine how many samples will be in the inclusion list?
# If no inclusion is provided, we will output sites only by default.
my (@samps, %famembers, %famrels);
if (defined $conf{Pedigree}{Include}) {
	if (lc($conf{Pedigree}{Include}) eq 'all') {
		@samps = map { $samps->{$_} } sort keys %$samps;
	}
	else {
		open my $fin, $conf{Pedigree}{Include} or die "Cannot open inclusion list file: $conf{Pedigree}{Include}";
		while(<$fin>) {
			my $iid = (split)[0];
			if (defined $samps->{$iid}) {
				push @samps, $samps->{$iid};
			}
		}
	}
	if (@samps > 0) {
		print STDERR "The genotypes of ", scalar(keys @samps), " samples will be written to output\n";
		my ($famsamp, $famrels) = fam_rels($conf{Pedigree}{File},
	 		{ twins => $conf{Pedigree}{Twins}, shorten => 1, ignore => $conf{Pedigree}{Ignore},
	 		  strict => 0, verbose => 0 });
		foreach my $samp (@samps) {
			my @members = grep { defined $samps->{$_} && $samp->{IID} ne $_ } @{$famsamp->{$samp->{FID}}};
			#my @members = grep { defined $samps->{$_}  } @{$famsamp->{$samp->{FID}}};
			if (@members > 0) {
				$famembers{$samp->{IID}} = \@members;
				$famrels{$samp->{IID}} = [ map { $famrels->{$samp->{IID}}{$_} } @members ];
			}
		}
	}
}
unless (@samps) {
	print STDERR "Only site-level information will be written to the output\n";
}


###############################################################
## Parse the config file to parse filters and output fields  ##
###############################################################

# Create parsers for Site and Geno levels
# Geno level parse is now optional, if it not provided, we will only filter by sites
my (%filters, %fields);
foreach my $level (qw|Site Geno|) {
	foreach my $name (grep { /Filter$/ } keys %{$conf{$level}}) {
	#foreach my $name (qw|Filter SNV_Filter Indel_Filter|) {
		next unless defined $conf{$level}{$name};
		$conf{$level}{$name} =~ s/^["']//;
		$conf{$level}{$name} =~ s/["']$//;
		my ($cb, $tokens) = sql_query($conf{$level}{$name}, 1);
		$filters{$level}{$name} = $cb;
		foreach my $tok (@$tokens) {
			if ($tok->[0] eq 'FIELD') {
				$fields{$level}{$tok->[1]} ++;
			}
		}
	}
}

# Site inclusion listf
my (%include, %exclude);
if (defined $conf{Site}{Include}) {
	my ($it, $fnames, $key_fields) =
		parse_tabfile($conf{Site}{Include}, $conf{Site}{Include_Fields}, 4, 4);
	while(my $dat = $it->()) {
		my $varid = join(":", @{$dat}{@$key_fields});
		$include{$varid} = 1;
	}
}
if (defined $conf{Site}{Exclude}) {
	my ($it, $fnames, $key_fields) =
		parse_tabfile($conf{Site}{Exclude}, $conf{Site}{Exclude_Fields}, 4, 4);
	while(my $dat = $it->()) {
		my $varid = join(":", @{$dat}{@$key_fields});
		$exclude{$varid} = 1;
	}
}

# fields and field types
my (%sfields, %gfields); # original field names used in expand_site and flatten_geno
%sfields = (CHROM => 1, POS => 1, REF => 1, ALT => 'A');
my ($site_f, $geno_f);   # expanded field names used in the filters and outputs
{
	my ($site_type, $geno_type) = get_fieldsinfo($opt->get_vcf, $conf{Module}{Site}, $conf{Module}{Geno});

	$site_f = parse_fstr($conf{Output}{Site}, 1);	
	foreach my $field (keys %$site_f, keys %{$fields{Site}}) {
		my $origfield = (split(q|\.|, $field))[0];
		unless(defined $site_type->{$origfield}) {
			if ($opt->get_strict) {
				die "Cannot find type for site field: $field";
			}
			else {
				warn "Cannot find type for site field: $field";
			}
		}
		chk_default(\%sfields, $origfield, $site_type->{$origfield});
	}

	$geno_f = parse_fstr($conf{Output}{Geno}, 1);
	foreach my $field (keys %$geno_f, keys %{$fields{Geno}}) {
		my $origfield = (split(q|\.|, $field))[0];
		unless(defined $geno_type->{$origfield}) {
			if ($opt->get_strict) {
				die "Cannot find type for geno field: $field";
			}
			else {
				warn "Cannot find type for geno field: $field";
			}
		}
		chk_default(\%gfields, $origfield, $geno_type->{$origfield});
	}
}

#######################################################
## Now iterate through the VCF file filter variants  ##
#######################################################
my $fout = IO::File->new($opt->get_out, "w");
unless (@samps) {
	print $fout join("\t", values %$site_f), "\n";
}
else {
	print $fout join("\t", qw|FamID IID Gender Pheno|, (values %$site_f), (values %$geno_f),	
						qw|FamMembers Relations Phenotypes|, 
						(map { "$_.FamMembers" } values %$geno_f)), "\n";
}

# Haplo flag indicate we are processing haploid portion of the genomes now 
my $pattern;
if (exists $conf{Pedigree}{Ignore}) {
	my $str = $conf{Pedigree}{Ignore};
	$str =~ s/^['"]//; $str =~ s/['"]$//;
	$pattern = qr/$str/;
}

my $it = $vcf->iter({strict => 0});
while(my ($site, $geno) = $it->()) {
	next if defined $exclchr->{$site->{CHROM}};
	
	# Determine in in male hap regions
	if ($site->{CHROM} eq 'chrX' || $site->{CHROM} eq 'X' ||
		$site->{CHROM} eq 'chrY' || $site->{CHROM} eq 'Y') {
		unless (hg_par($site->{CHROM}, $site->{POS}, $conf{Genome}{Build})) {
			$site->{_MALEHAP} = 1;
		}
	}
	# Customize sites and geno
	if (exists $conf{Module}{Site} || exists $conf{Module}{Geno}) {
		Genet::File::VCF::custom_site_geno($site, $geno, $samps, $pattern);
	}

	my @sites = expand_site($site, \%sfields, { nastring => $conf{Output}{NAChar} });
	# For each alternative allele
	for(my $ii = 0; $ii < @sites; $ii ++) {
		next if $sites[$ii]{ALT} eq '*';
		# Normalize site
		normalize($sites[$ii], $seq);
		# Determine variant type (after normalization)
		my $VT = var_type($sites[$ii]{REF}, $sites[$ii]{ALT});

		# Site level filter, with possible type-dependent customization
		next unless $filters{Site}{Filter}->($sites[$ii]);
		if (defined $filters{Site}{"${VT}_Filter"}) {
			next unless $filters{Site}{"${VT}_Filter"}->($sites[$ii]);
		}

		# Site level inclusion
		if (defined $conf{Site}{Include} || defined $conf{Site}{Exclude}) {
			my $varid = join(":", @{$sites[$ii]}{qw|CHROM POS REF ALT|});
			if (defined $conf{Site}{Include}) {
				next unless $include{$varid};
			}
			if (defined $conf{Site}{Exclude}) {
				next if $exclude{$varid};
			}
		}

		# If sample list is not provided, only output site level information if the site level filter is passed
		unless(@samps) {
			print $fout join("\t", 
				map { ref $_ eq 'ARRAY' ? join($conf{Output}{SepChar}, @$_) : $_ }
				map { defined $_ ? $_ : $conf{Output}{NAChar} }	@{$sites[$ii]}{keys %$site_f}), "\n";
			next;
		}

		# Genotype level filter
		# Loop through each sample
		foreach my $samp (@samps) {
			my @iids = ($samp->{IID});
			if (defined $samps->{$iids[0]}{DUP}) {
				push @iids, @{$samps->{$iids[0]}{DUP}};
			}
			my @sampdat = map { flatten_geno($geno->{$_}, \%gfields, $ii+1, scalar(@sites),
									{ nastring => $conf{Output}{NAChar} }) } @iids;
			
			my $filter = $filters{Geno}{Filter};
			if ($site->{_MALEHAP}) {
				if ($samp->{GENDER} eq 'Male' && defined $filters{Geno}{Haploid_Filter} ) {
					$filter = $filters{Geno}{Haploid_Filter};		
				}
			}
			
			my $passflag;
			if ($conf{Pedigree}{DupAction} eq 'Any') {
				if (any { $filter->($_) } @sampdat) {
					$passflag = 1;
				}
			}
			elsif  ($conf{Pedigree}{DupAction} eq 'Ignore') {
				if ($filter->($sampdat[0])) {
					$passflag = 1;
				}
			}
			elsif ($conf{Pedigree}{DupAction} eq 'All') {
				if (all { $filter->($_) } @sampdat) {
					$passflag = 1;
				}
			}
			else {
				die "Cannot recognize the way to process duplicated samples: $conf{Pedigree}{DupAction}";
			}

			if ($passflag) {
				# Gather geno level data for all family members
				my %mdat;
				if (defined $famembers{$samp->{IID}}) {
					foreach my $member (@{$famembers{$samp->{IID}}}) {
						my $gdat = flatten_geno($geno->{$member}, \%gfields, $ii+1, scalar(@sites),
										{ nastring => $conf{Output}{NAChar} });
						foreach my $gf (keys %$geno_f) {
							push @{$mdat{$gf}}, $gdat->{$gf} // $conf{Output}{NAChar};
						}
					}
				}
				print $fout join("\t", @{$samp}{qw|FID IID GENDER PHENO|},
					map { ref $_ eq 'ARRAY' ? join($conf{Output}{SepChar}, @$_) : $_ }
					map { defined $_ ? $_ : $conf{Output}{NAChar} }	@{$sites[$ii]}{keys %$site_f},
					arrange_genodat(\@sampdat, $geno_f, $conf{Output}{NAChar}),
					$famembers{$samp->{IID}},  $famrels{$samp->{IID}},
					[map { $samps->{$_}{PHENO} } @{$famembers{$samp->{IID}}}], @mdat{keys %$geno_f}), "\n";
			}
		}
	}
}











