#!/usr/bin/env perl

use strict;
use warnings;
use Carp;
use IO::File;
use IO::Prompt;
use Data::Dumper;
use FindBin qw|$Bin|;
use List::Util qw|sum|;
use List::MoreUtils qw|all uniq|;
#use Hash::Util qw|lock_hash_recurse|;
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
use Shared qw|parse_fstr parse_tabfile|;
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
	This is a standalone script to tabulate a list of DNVs from VCFs with customized output fields.  

Usages:
	dnv_table.pl --conf Config --vcf Cohort.vcf.gz [--tab Table --ped Cohort.ped] --out Prefix

Options:
	--out: (Required) Output file name prefix.
	--tab: A table of input DNVs. It overrides the input variant table in config.
	--ped: A standard 6-col PED file describing pedigree relationship of samples in VCF file.
			It overrides the PED file in config.

Notes:
	The input multi-sample VCF must strictly conforms to the format specification. We extract existing INFO
	and GENO fields from VCF, or customize additional fields using a perl module. Only samples present in 
	both PED and VCF files will be used in variant tabulation and VCF fields customization.  

	The script is intended for tabulating DNVs from genotypes of trios. For tabulating variants from any type 
	of family structure use var_table.pl.

Output:
	For variants in the input DNV table, if a trio or variant cannot be found in the VCF, then it will appear
	in Prefix.missing.txt. If both the trio and variant can be found in the VCF, the detailed site and genotype 
	related fields specified in the config will be shown in Prefix.txt. Genotype related fields will be repeated
	for offspring and both parents. When duplicated sample or MZ twin pairs are specified, genotype level metrics
	of dup or MZ twin will be packed together with original sample in the output.

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

my $trios = get_vcfped_samps($opt->get_vcf, $conf{Pedigree}{File}, 
		{ ignore => $conf{Pedigree}{Ignore}, trios => 1, twins => $conf{Pedigree}{Twins} });
my $samps = get_vcfped_samps($opt->get_vcf, $conf{Pedigree}{File}, 
		{ ignore => $conf{Pedigree}{Ignore}, trios => 0, twins => $conf{Pedigree}{Twins} });
unless(keys %$trios) {
	print STDERR "No trio was found in VCF/PED files\n";
	exit 1;
}
else {
	print STDERR "Found a total of ", scalar(keys %$trios), " trios in VCF/PED files\n";
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

# Go through each line in the variant table
my $outprefix = $opt->get_out;
if ($outprefix =~ /\.txt/) {
	$outprefix =~ s/\.txt//;
}
my $fout = IO::File->new($outprefix.".txt", "w");
my $fmis = IO::File->new($outprefix.".missing.txt", "w");
print $fout join("\t", qw|FamID IID Father Mother Gender Pheno|, (values %$site_f), (values %$xtrafds),
						(map { "$_.Offspring" } values %$geno_f),	(map { "$_.Father" } values %$geno_f),
						(map { "$_.Mother"  } values %$geno_f) ), "\n";
print $fmis join("\t", @$input_fields, @xfdsout, $conf{Output}{Miss}), "\n";

while(my $dat = $it->()) {
	my @site = @{$dat}{@$input_fields};
	my ($iid, $chr, $pos, $ref, $alt) = @site;

	# Hit flag
	my $flag;
	if (defined $trios->{$iid}) {
		my $trio = $trios->{$iid};
		my @triosamp;
		foreach (qw|IID DAD MOM|) {
			push @triosamp, $trio->{$_};
			if (defined $samps->{$trio->{$_}}{DUP}) {
				push @triosamp, @{$samps->{$trio->{$_}}{DUP}};
			}
		}
		my $iter;
		if (exists $conf{Module}{Site}) {
			$iter = $vcf->iter({strict => 0, region => sprintf("%s:%d-%d", $chr, $pos-1-$conf{VCF}{Extend}, 
															 		$pos+length($ref)+$conf{VCF}{Extend})});	
		}
		else {
			# If customizing site info is not needed, we can restrict to trio samples 
			$iter = $vcf->iter({strict => 0, samp => [@triosamp],
				region => sprintf("%s:%d-%d", $chr, $pos-1-$conf{VCF}{Extend}, 
												$pos+length($ref)+$conf{VCF}{Extend})});	
		}

		while(my $vdat = $iter->()) {
			my ($site, $geno) = @$vdat;

			if ($site->{CHROM} eq 'chrX' || $site->{CHROM} eq 'X' ||
				$site->{CHROM} eq 'chrY' || $site->{CHROM} eq 'Y') {
				unless (hg_par($site->{CHROM}, $site->{POS}, $conf{Genome}{Build})) {
					$site->{_MALEHAP} = 1;
				}
			}

			# Adding customization fields
			if (exists $conf{Module}{Site} || exists $conf{Module}{Geno}) {
				Genet::File::VCF::custom_site_geno($site, $geno, $samps, $pattern);
			}
			my @sites = expand_site($site, \%sfields, { nastring => $conf{Output}{NAChar} });
			@sites = map { normalize($_, $seq) } @sites if defined $seq; 
			for(my $ii = 0; $ii < @sites; $ii ++) {
				next if $sites[$ii]{ALT} eq '*';
				my $var = $sites[$ii];
				if ($var->{CHROM} eq $chr && $var->{POS} == $pos && 
					$var->{REF} eq $ref && $var->{ALT} eq $alt) {
					$flag = 1;
					
					my @triodat;
					for my $role (qw|IID DAD MOM|) {
						my @iids = ($trio->{$role});
						if (defined $samps->{$iids[0]}{DUP}) {
							push @iids, @{$samps->{$iids[0]}{DUP}};
						}
						my @gdat = map { flatten_geno($geno->{$_}, \%gfields, $ii+1, scalar(@sites),
											{ nastring => $conf{Output}{NAChar} }) } @iids;
						push @triodat => \@gdat;
					}
			
					print $fout join("\t", @{$trio}{qw|FID IID DAD MOM GENDER PHENO|}, 
						map { ref $_ eq 'ARRAY' ? join($conf{Output}{SepChar}, @$_) : $_ }
						map { defined $_ ? $_ : $conf{Output}{NAChar} }
						@{$var}{keys %$site_f}, @{$dat}{@xfdsout},
						map { arrange_genodat($_, $geno_f, $conf{Output}{NAChar}) } @triodat), "\n";
					last;
				}
			}
		}
		unless ($flag) {
			# Site not found
			print $fmis join("\t", $iid, $chr, $pos, $ref, $alt, @{$dat}{@xfdsout}, "VarNotFound"), "\n";
		}
	}
	else {
		# Trio not found
		print $fmis join("\t", $iid, $chr, $pos, $ref, $alt, @{$dat}{@xfdsout}, "TrioNotFound"), "\n";
	}
}


