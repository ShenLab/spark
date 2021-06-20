=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2018] EMBL-European Bioinformatics Institute

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

   http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

=head1 CONTACT

 Ensembl <dev@ensembl.org>
  
=cut

=head1 NAME

 SpliceRegion -- an inhouse VEP plugin for annotating splice regions.

=head1 SYNOPSIS

Provide further details on splice_region_variants. We annotate those variants
using the notation from Zhang et al. 2018, who found that :

1. In addition to the four “essential splice” nucleotides, positions D+3, D+4, D+5, D+6, D-1 and A+1
   are very significantly intolerant of mutations
2. Reference base T at D+6, G at D+5, A at D+4 and D+3, G at D-1, and G at A+1 are significantly less 
  tolerant of mutational alteration than are the other three reference bases at those same positions.

Reference: 

Zhang, S., Samocha, K. E., Rivas, M. A., Karczewski, K. J., Daly, E., Schmandt, B., … Daly, M. J. (2018). 
Base-specific mutational intolerance near splice sites clarifies the role of nonessential splice nucleotides. 
Genome Research. doi:10.1101/gr.231902.117

=cut

package SpliceRegion;

use strict;
use warnings;
use Data::Dumper;

use Bio::EnsEMBL::Variation::Utils::VariationEffect qw(overlap);
use Bio::EnsEMBL::Variation::Utils::Constants qw(%OVERLAP_CONSEQUENCES);

use base qw(Bio::EnsEMBL::Variation::Utils::BaseVepPlugin);


my %TERM_RANK;
for(my $ii = -3; $ii <= 8; $ii ++) {
  next unless $ii;
  my $dkey = $ii < 0 ? "D$ii" : "D+$ii";
  $TERM_RANK{$dkey} = 10+$ii;
}
for(my $ii = -8; $ii <= 3; $ii ++) {
  next unless $ii;
  my $akey = $ii < 0 ? "A$ii" : "A+$ii";
  $TERM_RANK{$akey} = $ii;
}


sub feature_types {
  return ['Transcript'];
}

sub get_header_info {
  return {
    SpliceRegion => "Splice region annotations",
  };
}

sub run {
  my ($self, $tva) = @_;

  my $vf = $tva->variation_feature;
  my ($vf_start, $vf_end) = ($vf->{start}, $vf->{end});

  my $is_insertion = 0;
  if($vf_start > $vf_end) {
    ($vf_start, $vf_end) = ($vf_end, $vf_start);
    $is_insertion = 1;
  }

  my $tv = $tva->transcript_variation;
  my $tr = $tv->transcript;
  my $vf_tr_seq = $tva->feature_seq;

  # define some variables depending on transcript strand
  my ($strand_mod, $donor_coord, $acc_coord);
  if($tr->strand > 0) {
    $strand_mod = 1;
    $donor_coord = 'start';
    $acc_coord = 'end';
  }
  else {
    $strand_mod = -1;
    $donor_coord = 'end';
    $acc_coord = 'start';
  }

  my %results;

  my @introns = @{$tv->_overlapped_introns($vf_start, $vf_end)};
  my @exons = @{$tv->_overlapped_exons($vf_start, $vf_end)};
 
  if (@introns > 1 || @exons > 1) {
    return { SpliceRegion => ['?']};
  }

  my @terms;
  if (@introns == 1) {
    my $intron = $introns[0];
    foreach my $ii (1..8) {
      push @terms, ["D+".$ii, $intron->{$donor_coord}+(($ii-1)*$strand_mod)];
      push @terms, ["A-".$ii, $intron->{$acc_coord}-(($ii-1)*$strand_mod)];
    }
  }
  if (@exons == 1) {
    my $exon = $exons[0];
    foreach my $ii (1..3) {
      push @terms, ["D-".$ii, $exon->{$acc_coord}-(($ii-1)*$strand_mod)];
      push @terms, ["A+".$ii, $exon->{$donor_coord}+(($ii-1)*$strand_mod)];
    }
  }


  foreach my $term (@terms) {
    my $pass = overlap($vf_start, $vf_end, $term->[1], $term->[1]);
    if ($pass) {
       $results{$term->[0]}++;
    }
  }

  return {} unless %results;

  return { SpliceRegion => [sort {$TERM_RANK{$a} <=> $TERM_RANK{$b}} keys %results] };
  
}

1;

