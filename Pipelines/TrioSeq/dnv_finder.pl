#!/usr/bin/env perl

use strict;
use warnings;
use Carp;
use IO::File;
use Data::Dumper;
use FindBin qw|$Bin|;
use File::Which;
use List::Util qw|sum|;
use List::MoreUtils qw|all any uniq|;
use Getopt::Lucid qw|:all|; 
use Perl6::Slurp;
use Module::Load;
use FaSlice;
use Genet::Ped;
use Genet::Var qw|normalize|;
use Genet::File::VCF;
use Genome::UCSC qw|hg_par|;
use Utils::Parser qw|sql_query|;
use Utils::Hash qw|merge_conf chk_default|;

use lib "$Bin/../lib";
use Shared qw|parse_fstr|;
use Variants qw|var_type get_vcfped_samps get_fieldsinfo expand_site flatten_geno arrange_genodat|;


my @spec =  (
	Param("conf|c")->valid(sub { -r }),
	Param("vcf|v")->valid(sub { -r $_ && -r "$_.tbi" }),
	Param("ped|p")->valid(sub { -r }),
	Param("bed|b")->valid(sub { -r $_  || /^(\w+):([\d,]+)\-([\d,]+)$/ || /^\w+$/ }),
	Param("out|o"),
	Switch("strict"),
	Keypair("param|par"),
	Switch("help|h")
	);

my $opt = Getopt::Lucid->getopt(\@spec);

if ($opt->get_help) {
	print STDERR <<EOF;

Purpose:
	This is a standalone script to filter and tabulate candidate de novo sequence variants (DNVs) from 
	joint genotype calls of family-based sequencing data.

Usage:
	dnv_finder.pl --conf Config --vcf Cohort.vcf.gz [--ped Cohort.ped] --out Output

Options:
	--ped: A standard 6-col PED file describing pedigree relationship of samples in VCF file.
			It overrides the PED file in config.
	--bed: Genomic regions in 3-col BED file format. Only variants in these regions will be filtered.
			It override the BED file in config.
	--strict: Under strict mode, all fields used in filtering or output must have a defined type. 

Notes:
	We assume that the input genotypes are in multi-sample VCF format that strictly conforms to the format 
	specification (https://samtools.github.io/hts-specs/VCFv4.2.pdf), and relationships of all individuals 
	are described in a 6-column PED file (http://zzz.bwh.harvard.edu/plink/data.shtml#ped). 
	Familial relationships should have been verified by genotypes. This script applies customizable filters 
	at both variant and genotype levels to look for candidate DNVs in each parents-offspring trio. 

	It is possible to a have non-familial samples in the PED file, as long as unrelated samples have different
	family IDs. It is also possible that VCF and PED file have only partial overlaps of samples. And only samples
	appear in both file will be used for filtering.

	Detailed DNV filter expressions are given in the config file. Site and genotype level filters make use of 
	fields from VCF. It is also possible to define additional customized fields from VCF using a perl module.
	The customized fields can be used the same way as original fields in filtering and output. When customizing
	fields from VCFs, only samples that present in both PED and VCF files will be used.

Output:
	The output is a table of DNV candidates passed filters with columns of site and genotype level information 
	specified in the config file. Genotype related fields will be repeated for offspring and both parents.  

	Information from duplicated samples or MZ twins can be used in filtering. When this is the case, genotype level 
	metrics of dup or MZ twin will be packed together with the orignial sample in the output.


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

my $seq = FaSlice->new(file => $conf{Genome}{Fasta});
my $exclchr = {};
if (defined $conf{Genome}{IgnoreChr}) {
	$exclchr = parse_fstr($conf{Genome}{IgnoreChr});
}

my $tabix = which("tabix");
my $vcf;
{
	# Get a VCF handler for creating VCF object
	my $fvcf;
	if (exists $conf{Genome}{Region}) {
		my $region = $conf{Genome}{Region};
		if (-f $region) {
			die "Cannot find tabix utility" unless $tabix;
			open $fvcf, "$tabix -p vcf -h -R " . $region ." ". $opt->get_vcf . " |" 
				or die "Cannot open tabix pipe";
		}
		elsif ($region =~ /^(\w+):([\d,]+)\-([\d,]+)$/ || $region =~ /^(\w+)$/) {
			die "Cannot find tabix utility" unless $tabix;
			$region =~ s/,//g;
			open $fvcf, "$tabix -p vcf -h " . $opt->get_vcf . " $region |"
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

my $trios = get_vcfped_samps($opt->get_vcf, $conf{Pedigree}{File}, 
		{ ignore => $conf{Pedigree}{Ignore}, trios => 1, twins => $conf{Pedigree}{Twins} });
my $samps = get_vcfped_samps($opt->get_vcf, $conf{Pedigree}{File}, 
		{ ignore => $conf{Pedigree}{Ignore}, trios => 0, twins => $conf{Pedigree}{Twins} });
my @trios = sort { $a->{FID} cmp $b->{FID} || $a->{IID} cmp $b->{IID} } values %$trios;
unless(@trios) {
	print STDERR "No trio was found in both VCF/PED files\n";
	exit 1;
} else {
	print STDERR "Found a total of ", scalar(@trios), " trios in VCF/PED files\n";
}

###############################################################
## Parse the config file to parse filters and output fields  ##
###############################################################
# Create parsers
my (%filters, %fields);
foreach my $level (qw|Site Geno|) {
	foreach my $name (keys %{$conf{$level}}) {
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


####################################################
## Now iterate through the VCF file to find DNVs  ##
####################################################
my $fout = IO::File->new($opt->get_out, "w");
print $fout join("\t", qw|FamID IID Father Mother Gender Pheno|, (values %$site_f), 
						(map { "$_.Offspring" } values %$geno_f),	(map { "$_.Father" } values %$geno_f),
						(map { "$_.Mother"  } values %$geno_f)), "\n";

my $pattern;
if (exists $conf{Pedigree}{Ignore}) {
	my $str = $conf{Pedigree}{Ignore};
	$str =~ s/^['"]//; $str =~ s/['"]$//;
	$pattern = qr/$str/;
}

my $it = $vcf->iter({strict => 0});
while(my $dat = $it->()) {
	my ($site, $geno) = @$dat;
	next if defined $exclchr->{$site->{CHROM}};
	
	# Determine in in male hap regions, this is also site cusomization
	# but considered as internal
	if ($site->{CHROM} eq 'chrX' || $site->{CHROM} eq 'X' ||
		$site->{CHROM} eq 'chrY' || $site->{CHROM} eq 'Y' ) {
		unless (hg_par($site->{CHROM}, $site->{POS}, $conf{Genome}{Build})) {
			$site->{_MALEHAP} = 1; # to indicate under the male-haploid mode now
		}
	}

	# Customize sites and geno
	if (exists $conf{Module}{Site} || exists $conf{Module}{Geno}) {
		Genet::File::VCF::custom_site_geno($site, $geno, $samps, $pattern);
	}

	my @sites = expand_site($site, \%sfields, { nastring => $conf{Output}{NAChar} });
	# For each alternative allele
	for(my $ii = 0; $ii < @sites; $ii ++) {
		# Skip the star allele
		next if $sites[$ii]{ALT} eq '*';
		# Normalize site
		normalize($sites[$ii], $seq);
		# Determine variant type (after normalization)
		my $VT = var_type($sites[$ii]{REF}, $sites[$ii]{ALT});

		# Site level filter, with possible type dependent customization
		next unless $filters{Site}{Filter}->($sites[$ii]);
		if (defined $filters{Site}{"${VT}_Filter"}) {
			next unless $filters{Site}{"${VT}_Filter"}->($sites[$ii]);
		}

		# Genotype level filter
		# For each IID-Father-Mother trio
		foreach my $trio (@trios) {
			my %sampdat;
			foreach my $role (qw|IID DAD MOM|) {
				my @iids = ($trio->{$role});
				if (defined $samps->{$iids[0]}{DUP}) {
					push @iids, @{$samps->{$iids[0]}{DUP}};
				}
				my @gdat = map { flatten_geno($geno->{$_}, \%gfields, $ii+1, scalar(@sites),
									{ nastring => $conf{Output}{NAChar}, 
									   extra => $_.":$site->{CHROM}:$site->{POS}:$site->{REF}:$sites[$ii]{ALT}" }) } @iids;
				$sampdat{$role} = \@gdat;
			}
			
			# Filter names
			my ($offspring, $father, $mother) = qw|Offspring Parent Parent|;
			if ($site->{_MALEHAP} ) {
				if ($trio->{GENDER} eq 'Male') {
					$offspring = defined $filters{Geno}{Haploid_Offspring} ? "Haploid_Offspring" : "Offspring";		
				}
				$father = defined $filters{Geno}{Haploid_Parent} ? "Haploid_Parent" : "Parent";
			}

			my $passflag;
			if ($conf{Pedigree}{DupAction} eq 'Any') {
				if( (any { $filters{Geno}{$offspring}->($_) } @{$sampdat{IID}}) &&
					(any { $filters{Geno}{$father}->($_) } @{$sampdat{DAD}}) &&
					(any { $filters{Geno}{$mother}->($_) } @{$sampdat{MOM}}) ) {
					$passflag = 1;
				}
			}
			elsif ($conf{Pedigree}{DupAction} eq 'Ignore') {
				if ( $filters{Geno}{$offspring}->($sampdat{IID}[0]) &&
					 $filters{Geno}{$father}->($sampdat{DAD}[0]) &&
					 $filters{Geno}{$mother}->($sampdat{MOM}[0]) ) {
					$passflag = 1;
				}
			}
			elsif ($conf{Pedigree}{DupAction} eq 'All') {
				if( (all { $filters{Geno}{$offspring}->($_) } @{$sampdat{IID}}) &&
					(all { $filters{Geno}{$father}->($_) } @{$sampdat{DAD}}) &&
					(all { $filters{Geno}{$mother}->($_) } @{$sampdat{MOM}}) ) {
					$passflag = 1;
				}
			}
			else {
				die "Cannot recognize the way to process duplicated samples: $conf{Pedigree}{DupAction}";
			}

			if ($passflag) {
				# Output if passing filters above, QC metrics from all duplicated sample will be in the output
				print $fout join("\t", @{$trio}{qw|FID IID DAD MOM GENDER PHENO|},
					map { ref $_ eq 'ARRAY' ? join($conf{Output}{SepChar}, @$_) : $_ }
					map { defined $_ ? $_ : $conf{Output}{NAChar} }	@{$sites[$ii]}{keys %$site_f},
					map { arrange_genodat($sampdat{$_}, $geno_f, $conf{Output}{NAChar}) } qw|IID DAD MOM|), "\n";
			}
		}
	}
}



