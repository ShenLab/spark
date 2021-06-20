#!/usr/bin/env perl

use strict;
use warnings;
use Carp;
use IO::File;
use IO::Prompt;
use Data::Dumper;
use FindBin qw|$Bin|;
use List::Util qw|sum|;
use List::MoreUtils qw|all any uniq|;
use Getopt::Lucid qw|:all|; 
use Config::Std;
use Perl6::Slurp;
use Module::Load;
use FaSlice;
use Genome::UCSC qw|hg_par|;
use Genet::Var qw|normalize|;
use Genet::File::VCF;
use Utils::Hash qw|merge_conf chk_default|;
use Utils::File::Iter qw|iter_file|;

use lib "$Bin/../lib";
use Shared qw|parse_fstr parse_tabfile fam_rels|;
use Variants qw|get_vcfped_samps get_fieldsinfo expand_site flatten_geno arrange_genodat|;


my @spec =  (
	Param("vcf|v")->valid(sub { -r $_ && -r "$_.tbi" }),
	Param("ped|p")->valid(sub { -r }),
	Param("tab|t")->valid(sub { -r }),
	Param("conf")->valid(sub { -r }),
	Param("out|o"),
	Keypair("param|par"),
	Switch("help|h")
	);

my $opt = Getopt::Lucid->getopt(\@spec);

if ($opt->get_help) {
	print STDERR <<EOF;
Purpose:
	This is a standalone script to tabulate a list of variants from VCFs with custmized output fields.

Usages:
	var_table.pl --conf Config --vcf Cohort.gvcf.gz [--tab Table --ped Cohort.ped] --out Prefix

Options:
	--out: (Required) Output file name prefix.
	--tab: A table of input DNVs. It overrides the input variant table in config.
	--ped: A standard 6-col PED file describing pedigree relationship of samples in VCF file.
		   Pedigree file is required even if all samples are unrelated. It overrides the PED file in config. 

Notes:
	The script extends dnv_table to tabulate any type of variants from VCFs. Site and genotype level
	information is extracted from known or customized fields from VCF. PED file must be provided
	even if all samples are unrelated. Only samples appearing in both VCF and PED wll be processed
	and used in customizing VCF fields.

Output:
	For variants in the input table, if the sample or variant cannot be found in the VCF, then it will 
	appear in Prefix.missing.txt. If both the sample and variant can be found in the VCF, the detailed
	site and genotype related fields specified in the config will be shown in Prefix.txt. Genotype related
	info for family members will be packed in one column for each fields.  

EOF
	exit 1;
}

$opt->validate({ requires => [qw|vcf out conf|] });

#######################
# Parsing input files #
#######################

my %conf = merge_conf($opt->get_conf, $opt->get_param); 
if ($opt->get_tab) {
	$conf{Input}{File} = $opt->get_tab;
}
if ($opt->get_ped) {
	$conf{Pedigree}{File} = $opt->get_ped;
}
unless(-f $conf{Pedigree}{File}) {
	croak "Cannot find pedigree file: $conf{Pedigree}{File}";
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

# For each sample, prepare a list of family members with data (sorted by ID)
my ($famsamp, $famrels) = fam_rels($conf{Pedigree}{File},
	 	{ twins => $conf{Pedigree}{Twins}, shorten => 1, ignore => $conf{Pedigree}{Ignore}, verbose => 0 });
my %famembers;
foreach my $iid (keys %$samps) {
	my $fid = $samps->{$iid}{FID} // do { die "Cannot find family ID for $iid" };
	my @members = grep { defined $samps->{$_} && $_ ne $iid } @{$famsamp->{$fid}};
	if (@members > 0) {
		$famembers{$iid} = \@members;
	}
}

# Sequence file needed when performing left alignment
my $vcf = Genet::File::VCF->new($opt->get_vcf);
my $seq;
if ($conf{VCF}{Norm} =~ /^T/i || $conf{VCF}{Norm} =~ /^Y/i) {
	$seq = FaSlice->new(file => $conf{Genome}{Fasta});
}


################################################
# Parse config file to find output fields and  #
# determine field types                        #
################################################
my (%sfields, %gfields); # original field names used in expand_site and flatten_geno
my ($site_f, $geno_f);   # field names used in the filters and output
{
	my ($site_type, $geno_type) = get_fieldsinfo($opt->get_vcf, $conf{Module}{Site}, $conf{Module}{Geno});

	$site_f = parse_fstr($conf{Output}{Site}, 1);	
	foreach my $field (keys %$site_f) {
		my $origfield = (split(q|\.|, $field))[0];
		chk_default(\%sfields, $origfield, 
					$site_type->{$origfield} // do { carp "Cannot find type for site field: $field" });
	}

	$geno_f = parse_fstr($conf{Output}{Geno}, 1);
	foreach my $field (keys %$geno_f) {
		my $origfield = (split(q|\.|, $field))[0];
		chk_default(\%gfields, $origfield, 
					$geno_type->{$origfield} // do { carp "Cannot find type for geno field: $field" });
	}
}


############################################
# Parse known variants table and keep the  #
# column names for the five key fields     #
############################################
my $pattern;
if (exists $conf{Pedigree}{Ignore}) {
	my $str = $conf{Pedigree}{Ignore};
	$str =~ s/^['"]//; $str =~ s/['"]$//;
	$pattern = qr/$str/;
}

my ($it, $fnames, $input_fields) =
	parse_tabfile($conf{Input}{File}, $conf{Input}{Fields}, 5, 5);
my ($xtrafds, @xfdsout);
if (exists $conf{Input}{XtraFields}) {
	$xtrafds = parse_fstr($conf{Input}{XtraFields}, 1);
	@xfdsout = keys %$xtrafds;
	foreach my $field (@xfdsout) {
		unless (grep { $field eq $_ } @$fnames) {
			croak "Cannot find field $field in the input table";
		}
	}
}
# Slurp all variants that should be tabulated, IID and extra fields will be stored
my %vars;
while(my $dat = $it->()) {
	my @site = @{$dat}{@$input_fields};
	my ($iid, $chr, $pos, $ref, $alt) = @site;
	my @xtrainfo = @{$dat}{@xfdsout};
	my $varid = join("\t", $chr, $pos, $ref, $alt);
	push @{$vars{$varid}}, [$iid, \@xtrainfo];
}

# Go through each line in the variant table
my $outprefix = $opt->get_out;
if ($outprefix =~ /\.txt/) {
	$outprefix =~ s/\.txt//;
}
my $fout = IO::File->new($outprefix.".txt", "w");
my $fmis = IO::File->new($outprefix.".missing.txt", "w");
print $fout join("\t", qw|FamID IID Gender Pheno|, (values %$site_f), values %$xtrafds, (values %$geno_f),
				qw|FamMembers Relations Phenotypes|, (map { "$_.FamMembers" } values %$geno_f)), "\n";
print $fmis join("\t", @$input_fields, @xfdsout, $conf{Output}{Miss}), "\n";

foreach my $varid (sort keys %vars) {
	my ($chr, $pos, $ref, $alt) = split(/\t/, $varid);
	my @samps = map { $_->[0] } @{$vars{$varid}};

	if (any { defined $samps->{$_} } @samps) {
		my $iter;
		if (exists $conf{Module}{Site}) {
			$iter = $vcf->iter({strict => 0, region => sprintf("%s:%d-%d", $chr, $pos-1-$conf{VCF}{Extend}, 
																$pos+length($ref)+$conf{VCF}{Extend})});	
		}
		else {
			# If customize site info is not needed, we can restrict to all family samples to save time
			my @famsamp;
			foreach my $iid (grep { defined $samps->{$_} } @samps) {
				# Dup samples will appear in famsamps
				push @famsamp, $iid;
				if (defined $famembers{$iid}) {
					push @famsamp, @{$famembers{$iid}};
				}
			}
			@famsamp = uniq sort @famsamp;
			$iter = $vcf->iter({strict => 0, samp => \@famsamp,
					region => sprintf("%s:%d-%d", $chr, $pos-1-$conf{VCF}{Extend}, 
										$pos+length($ref)+$conf{VCF}{Extend})});	
		}

		# Hit flag
		my $flag;
		while(my $vdat = $iter->()) {
			my ($site, $geno) = @$vdat;
			# Adding customization fields
			if ($site->{CHROM} eq 'chrX' || $site->{CHROM} eq 'X' ||
				$site->{CHROM} eq 'chrY' || $site->{CHROM} eq 'Y') {
				unless (hg_par($site->{CHROM}, $site->{POS}, $conf{Genome}{Build})) {
					$site->{_MALEHAP} = 1;
				}
			}
			if (exists $conf{Module}{Site} || exists $conf{Module}{Geno}) {
				Genet::File::VCF::custom_site_geno($site, $geno, $samps, $pattern);
			}
			
			my @sites = expand_site($site, \%sfields, { nastring => $conf{Output}{NAChar} });
			@sites = map { normalize($_, $seq) } @sites if defined $seq; 
			for(my $ii = 0; $ii < @sites; $ii ++) {
				next if $sites[$ii]{ALT} eq '*';
				my $var = $sites[$ii];
				my $jj = $ii + 1;
				if ($var->{CHROM} eq $chr && $var->{POS} == $pos && 
					$var->{REF} eq $ref && $var->{ALT} eq $alt) {
					$flag = 1;
					foreach my $sampinfo (@{$vars{$varid}}) {
						my ($iid, $xtrainfo) = @$sampinfo;
						if (defined $samps->{$iid}) {
							# If a variant was found, we will collect related genotype information for all family members
							my $idat = flatten_geno($geno->{$iid}, \%gfields, $jj, scalar(@sites),
														{ nastring => $conf{Output}{NAChar} });
							my %mdat;
							if (defined $famembers{$iid}) {
								foreach my $member (@{$famembers{$iid}}) {
									my $gdat = flatten_geno($geno->{$member}, \%gfields, $jj, scalar(@sites),
																{ nastring => $conf{Output}{NAChar} });
									foreach my $gf (keys %$geno_f) {
										push @{$mdat{$gf}}, $gdat->{$gf} // $conf{Output}{NAChar};
									}
								}
							}
							print $fout join("\t", @{$samps->{$iid}}{qw|FID IID GENDER PHENO|}, 
								map { ref $_ eq 'ARRAY' ? join($conf{Output}{SepChar}, @$_) : $_ }
								map { defined $_ ? $_ : $conf{Output}{NAChar} }
								@{$var}{keys %$site_f}, @$xtrainfo, @{$idat}{keys %$geno_f},
								$famembers{$iid}, [map { $famrels->{$iid}{$_} } @{$famembers{$iid}}],
								[map { $samps->{$_}{PHENO} } @{$famembers{$iid}}],	
								@mdat{keys %$geno_f}), "\n";
						}
						else {
							print $fmis join("\t", $iid, $chr, $pos, $ref, $alt, @$xtrainfo, "SampNotFound"), "\n";
						}
					}
					last;
				}
			}
		}
		unless ($flag) {
			# Site not found
			foreach my $sampinfo (@{$vars{$varid}}) {
				my ($iid, $xtrainfo) = @$sampinfo;
				print $fmis join("\t", $iid, $chr, $pos, $ref, $alt, @$xtrainfo, "VarNotFound"), "\n";
			}
		}
	}
	else {
		# None of the samples can be found in the VCF/PED
		foreach my $sampinfo (@{$vars{$varid}}) {
			my ($iid, $xtrainfo) = @$sampinfo;
			print $fmis join("\t", $iid, $chr, $pos, $ref, $alt, @$xtrainfo, "SampNotFound"), "\n";
		}
	}
}


