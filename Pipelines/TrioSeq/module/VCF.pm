package Genet::File::VCF;

use strict;
use warnings;
use Carp;
use Data::Dumper;
use List::Util qw|sum|;
use List::MoreUtils qw|any all|;

=head1 NOTES

Subroutines in this module will be useful for de novo variants and ultra-rare variants analysis.

Assuming following GENO fields should be present in VCF:
  * AD, DP   Allelic and total depth
  * GT       The most likely genotype from VCF

The following custom GENO fields were calculated:
  * AB    Allele balance, calculated as AD/Total DP, it should have the same type as AD.
  * GENO  Genotype defined *for each alt allele*. Possible values of GENO are listed below
          Assuming we are currently looking at at allele 1:
          0/1 or 1/2: Het (on male X: Het?), 0/2 or 2/3: HetOther
          0/0: HomRef (on male X: HemiRef), 1/1: HomAlt (on male X: HemiAlt), 2/2: HomOther (on male X: HemiOther)
          ./1: Het(HalfMiss) (on male X: Het?(HalfMiss))
          ./0 or ./2: HetOther(HalfMiss) (on male X: HetOther?(HalfMiss))

The following custom INFO fields were calculated from sample genotypes:
  * NONMIS/FNONMIS   Number and percentage of samples with non-missing genotypes.
                     Non-missing genotype will be defined as anything other than './.' or '.'  
  * FAM_AC/AN/AF   Number and fraction of families with member(s) carrying the alt allele.
          This is similar to AC/AF/AN fields but use family as the unit, and carrier freq
          not allele freq is calculated. Only samples that appear in both PED and VCF
          files will be considered. Families in which *NONE* of the family members has genotype
          calls will be taken as missing.

=cut

sub custom_site_geno {
	my ($site, $geno, $samps, $pattern) = @_;

	#  Customize genodata
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
					$gdat->{AB}[$ii] = sprintf("%.3f", $gdat->{AB}[$ii]);
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
		

	# In case AC/AN does not appear in the original INFO, we will fix them
	# AC/AN are defined on *ALL* samples in VCF 
	if (!defined $site->{INFO}{AC} || !defined $site->{INFO}{AN}) {
		allele_count($site, $geno);
	}

	my (%famswal, %famstot); # Families with alt allele and with nonmissing genotypes
	my ($nonmis, $ntot) = (0, 0); # Samples with nonmissing genotypes and total samples
	my $nalt = @{$site->{ALT}};
	
	while(my ($iid, $gdat) = each %$geno) {
		next unless defined $samps->{$iid};
		my $fid = $samps->{$iid}{FID};
		$ntot ++;
		if (defined $gdat->{GT} && $gdat->{GT} ne '.' && $gdat->{GT} ne './.' && $gdat->{GT} ne '.|.') {
			$nonmis ++;
			$famstot{$fid} ++;
			my @als = split(/\||\//, $gdat->{GT});
			foreach my $al (@als) {
				next if $al eq '.'; # Skip half missing genotypes
				$famswal{$al}{$fid} ++;
			}
		}
	}

	my %info;
	$info{NONMIS} = $nonmis;
	if ($ntot > 0) {
		my $p_nonmis =  $nonmis/$ntot;
		if ($p_nonmis > 0.001) {
			$info{FNONMIS} = sprintf("%.3f", $p_nonmis);
		}
		elsif ($p_nonmis = 0) {
			$info{FNONMIS} = 0;
		}
		else {
			$info{FNONMIS} = sprintf("%.2e", $p_nonmis);
		}
	}
	else {
		$info{FNONMIS} = ".";
	}
	
	$info{FAM_AN} = scalar(keys %famstot);
	for (my $jj = 1; $jj <= $nalt; $jj ++) {
		$info{FAM_AC}[$jj-1] = scalar(keys %{$famswal{$jj}});
		if ($info{FAM_AN} > 0) {
			my $famaf = $info{FAM_AC}[$jj-1]/$info{FAM_AN};
			if ($famaf > 0.001) {
				$info{FAM_AF}[$jj-1] = sprintf("%.3f", $famaf);
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
	}
	
	while(my ($field, $dat) = each %info) {
		$site->{INFO}{$field} = $dat;
	}
	return 1;
}

1;