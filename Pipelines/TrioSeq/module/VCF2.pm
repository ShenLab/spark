package Genet::File::VCF;

use strict;
use warnings;
use Carp;
use Data::Dumper;
use List::Util qw|sum|;
use List::MoreUtils qw|any all|;

=head1 NOTES

Subroutines in this module will be useful for rare variants analysis that are robust
to low level sample contamination.

Assuming following GENO fields should be present in VCF:
  * AD, DP   Allelic and total depth
  * GT       The most likely genotype from VCF

The following custom GENO fields were calculated:
  * AB     Allele balance, calculated as AD/Total DP, it should have the same type as AD.
  * GENO   Genotype defined *for each alt allele*. Possible values of GENO are listed below
          Assuming we are currently looking at at allele 1:
          0/1 or 1/2: Het (on male X: Het?), 0/2 or 2/3: HetOther
          0/0: HomRef (on male X: HemiRef), 1/1: HomAlt (on male X: HemiAlt), 2/2: HomOther (on male X: HemiOther)
          ./1: Het(HalfMiss) (on male X: Het?(HalfMiss))
          ./0 or ./2: HetOther(HalfMiss) (on male X: HetOther?(HalfMiss))

The following custom INFO fields were calculated from sample genotypes:
  * NONMIS/FNONMIS   Number and percentage of samples with non-missing genotypes.
                     Non-missingness will be defined as anything other than './.' or '.' 
  * FAM_AC/AN/AF   Number and fraction of families with member(s) carrying the alt allele.
          This is similar to AC/AF/AN fields but use family as the unit, and carrier freq
          not allele freq is calculated. Only samples that appear in both PED and VCF
          files will be considered. Families in which *NONE* of the family members has genotype
          calls will be taken as missing.                

The following are newly added in this module
  * NHDP_DAGT/NHQ_DAGT/FHQ_DAGT  Number of high DP (defined as at least 15x) 
                                 diploid alt genotypes (DAGT, including Het or Alt), 
          And number/fraction of high-quality diploid non-ref genotypes.
          High-quality diploid alt genotypes are defined as not on male chrX and 
          having AB>0.25 & AB<0.75 for Het or Het(HalfMiss), AB>0.9 for HomAlt
          They will have one value for each alt allele (similar to AC/AF).

  * POP_AC/AN/AF   Population allele counts and frequencies.
          Population allele frequency is calculated based on nonmissing genotypes in founders.
          We will also account for the ploidy on chrX.

=cut


sub custom_site_geno {
	my ($site, $geno, $samps, $pattern) = @_;

	# Calculate AB and re-classify genotypes
	while(my ($iid, $gdat) = each %$geno) {
		#my $sex = defined $samps->{$iid} ? $samps->{$iid}{GENDER} : "NA";
		my $sex;
		if (defined $samps->{$iid}) {
			$sex = $samps->{$iid}{GENDER};	
		}
		elsif (defined $pattern && $iid =~ /$pattern/) {
			(my $origid = $iid) =~ s/$pattern//;
			$sex = defined $samps->{$origid} ? $samps->{$origid}{GENDER} : "NA";
		}
		else {
			$sex = "NA";
		}

		unless (defined $gdat->{AD} && defined $gdat->{DP}) {
			next;
		}

		$gdat->{AB} = [] if defined $gdat->{AB}; # customize AB, removing existing AB
		for(my $ii = 0; $ii < @{$gdat->{AD}}; $ii ++) {
			if ($gdat->{AD}[$ii] eq '.') {
				$gdat->{AB}[$ii] = '.';
			}
			else {
				if ($gdat->{DP} > 0) {
					$gdat->{AB}[$ii] = $gdat->{AD}[$ii]/$gdat->{DP};
					$gdat->{AB}[$ii] = sprintf("%.4f", $gdat->{AB}[$ii]);
				}
				else {
					$gdat->{AB}[$ii] = '.';
				}
			}
		}

		my @als = split(/\||\//, $gdat->{GT});
		for(my $ii = 1; $ii < @{$gdat->{AD}}; $ii ++) {
			if (any {  $_ eq '.' } @als) {
				if (all {  $_ eq '.' } @als) {
					$gdat->{GENO}[$ii-1] = '.';
				}
				else {
					unless(@als ==2) {
						print STDERR Dumper $site;
						die "Not a diploid genotype: $iid $gdat->{GT}";
					}
					if ($als[0] eq $ii || $als[1] eq $ii) {
						$gdat->{GENO}[$ii-1] = "Het(HalfMiss)";
					}
					else {
						$gdat->{GENO}[$ii-1] = "HetOther(HalfMiss)";
					}
					if ($site->{_MALEHAP} && $sex eq 'Male') {
						$gdat->{GENO}[$ii-1] =~ s/\(/?(/;
					}
				}
			}
			else {
				unless(@als ==2) {
					print STDERR Dumper $site;
					die "Not a diploid genotype: $iid $gdat->{GT}";
				}

				if ($als[0] eq $als[1] ) {
					my $atype;
					if ($als[0] eq '0') {
						$atype = 'Ref';
					} 
					elsif ($als[0] eq $ii) {
						$atype = 'Alt';
					} else {
						$atype = 'Other';
					}
					$gdat->{GENO}[$ii-1] = "Hom$atype";
					if ($site->{_MALEHAP} && $sex eq 'Male') {
						$gdat->{GENO}[$ii-1] = "Hemi$atype";
					}
				}
				else {
					if ($als[0] eq $ii || $als[1] eq $ii) {
						$gdat->{GENO}[$ii-1] = "Het";
					}
					else {
						$gdat->{GENO}[$ii-1] = "HetOther";
					}
					if ($site->{_MALEHAP} && $sex eq 'Male') {
						$gdat->{GENO}[$ii-1] .= "?";
					}
				}
			}
		}

		$geno->{$iid} = $gdat;
	}


	# Refill AC/AN if not present, keep in mind that are defined on ALL samples in VCF
	# and they are not population allele freq
	if (!defined $site->{INFO}{AC} || !defined $site->{INFO}{AN}) {
		allele_count($site, $geno);
	}

	# Find founders
	my %founders = map { $_ => $samps->{$_}{GENDER} }
		grep { $samps->{$_}{DAD} eq '0' && $samps->{$_}{MOM} eq '0' } keys %$samps;

	my $nalt = @{$site->{ALT}};

	my ($totsamp, $totnonmis, $fdnonmischr) = (0) x 3; 
	# Total number of samples, with nonmissing genotypes, number of founder chromsomes
	my @fdobschr = (0) x $nalt; # Number of observed chromosomes for each alt allele
	my (%famswal, %famstot); # Families with alt allele and with nonmissing genotypes
	my @nHiDpDaGT = (0) x $nalt; # Number of high-depth diploid alt-genotypes
	my @nHiQDaGT = (0) x $nalt; # Number of high-quality diploid alt-genotypes 


	foreach my $iid (keys %$samps) {
		# Sample must be present in both PED and VCF files
		next unless defined $geno->{$iid};
		# Count the total number of samples
		$totsamp ++;
		my $sex = $samps->{$iid}{GENDER};
		my $fid = $samps->{$iid}{FID};
		my $gdat = $geno->{$iid};
		# GENO should be defined, otherwise will be treated as missing
		next unless defined $gdat->{GENO};
		# First count the number of samples with non-missing genotypes
		# and number of founder chromosomes 
		if (defined $gdat->{GT} && $gdat->{GT} ne '.' && $gdat->{GT} ne './.' && $gdat->{GT} ne '.|.') {
			$totnonmis ++;
			$famstot{$fid} ++;
			if ($founders{$iid}) {
				if ($site->{_MALEHAP} && $sex eq 'Male') {
					$fdnonmischr += 1;
				}
				else {
					$fdnonmischr += 2;
				}
			}
		}
		# Tally the families carrying each alt allele
		# And count the number of founder chromsomes with alt allele
		for(my $jj = 1; $jj <= $nalt; $jj ++) {
			next if $gdat->{GENO}[$jj-1] eq '.';
			if ($gdat->{GENO}[$jj-1] =~ /^Het(?!O)/ || $gdat->{GENO}[$jj-1] =~ /Alt/) {
				$famswal{$jj}{$fid} ++;
				if ($gdat->{DP} >= 15 && !($site->{_MALEHAP} && $sex eq 'Male')) {
					$nHiDpDaGT[$jj-1] ++;
					if ($gdat->{GENO}[$jj-1] =~ /^Het(?!O)/ && $gdat->{AB}[$jj]>0.25 && $gdat->{AB}[$jj]<0.75 ||
						$gdat->{GENO}[$jj-1] =~ /Alt$/ && $gdat->{AB}[$jj]>0.9) {
						$nHiQDaGT[$jj-1] ++;
					}
				}
			}
			if ($founders{$iid}) {
				if ($gdat->{GENO}[$jj-1] =~ /^Het(?!O)/) {
					$fdobschr[$jj-1] += 1;
				}
				elsif ($gdat->{GENO}[$jj-1] eq 'HomAlt') {
					$fdobschr[$jj-1] += 2;
				}
				elsif ($gdat->{GENO}[$jj-1] eq 'HemiAlt') {
					$fdobschr[$jj-1] += 1;
				}		
			}	
		}
	}

	# Total number of samples with nonmissing genotype assignment
	my %info;
	$info{NONMIS} = $totnonmis;
	$info{FNONMIS} = $totsamp > 0 ? sprintf("%.4f", $totnonmis/$totsamp) : '.';
	$info{POP_AN} = $fdnonmischr;
	$info{FAM_AN} =  scalar(keys %famstot);

	for(my $jj = 1; $jj <= $nalt; $jj ++) {
		$info{POP_AC}[$jj-1] = $fdobschr[$jj-1];	
		if ($info{POP_AN} > 0) {
			my $af = $info{POP_AC}[$jj-1]/$info{POP_AN};
			if ($af >= 0.001) {
				$info{POP_AF}[$jj-1] = sprintf("%.4f", $af);
			}
			elsif ($af == 0) {
				$info{POP_AF}[$jj-1] = 0;
			}
			else {
				$info{POP_AF}[$jj-1] = sprintf("%.2e", $af);
			}
		}
		else {
			$info{POP_AF}[$jj-1] = '.';
		}

		$info{FAM_AC}[$jj-1] = scalar(keys %{$famswal{$jj}});
		if ($info{FAM_AN} > 0) {
			my $famaf = $info{FAM_AC}[$jj-1]/$info{FAM_AN};
			if ($famaf > 0.001) {
				$info{FAM_AF}[$jj-1] = sprintf("%.4f", $famaf);
			}
			elsif ($famaf == 0) {
				$info{FAM_AF}[$jj-1] = 0;
			}
			else {
				$info{FAM_AF}[$jj-1] = sprintf("%.2e", $famaf);
			}
		}
		else {
			$info{FAM_AF}[$jj-1] = ".";
		}

		$info{NHDP_DAGT}[$jj-1] = $nHiDpDaGT[$jj-1];
		$info{NHQ_DAGT}[$jj-1] = $nHiQDaGT[$jj-1];
		if ($nHiDpDaGT[$jj-1] > 0) {
			my $fHiQDaGT = $nHiQDaGT[$jj-1]/$nHiDpDaGT[$jj-1];
			if ($fHiQDaGT > 0.001) {
				$info{FHQ_DAGT}[$jj-1] = sprintf("%.4f", $fHiQDaGT);
			}
			elsif ($fHiQDaGT == 0) {
				$info{FHQ_DAGT}[$jj-1] = 0;
			}
			else {
				$info{FHQ_DAGT}[$jj-1] = sprintf("%.2e", $fHiQDaGT);
			}
		}
		else {
			$info{NHQ_DAGT}[$jj-1] = '.';
			$info{FHQ_DAGT}[$jj-1] = '.';
		}
	}
	
	while(my ($field, $dat) = each %info) {
		$site->{INFO}{$field} = $dat;
	}
	return 1;
}

1;