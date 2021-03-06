=head1 CONTACT                                                                                                       

 Konrad Karczewski <konradjkarczewski@gmail.com>
 
=cut

=head1 NAME

 Context -- VEP plugin to determine sequence context of SNVs.

=head1 SYNOPSIS

 mv context.pm ~/.vep/Plugins
 perl variant_effect_predictor.pl -i variations.vcf --plugin context
 perl variant_effect_predictor.pl -i variations.vcf --plugin context,N

=head1 DESCRIPTION

 A VEP plugin that writes N base (default = 1) context around the variant.

=cut

package Context;

use strict;
use warnings;

use Bio::EnsEMBL::Variation::Utils::BaseVepPlugin;
use base qw(Bio::EnsEMBL::Variation::Utils::BaseVepPlugin);
use Bio::Perl;

sub get_header_info {
    return {
        Context => "base context around the variant"
    };
}

sub feature_types {
    return ['Transcript', 'Intergenic'];
}

sub new {
    my $class = shift;

    my $self = $class->SUPER::new(@_);
    
    $self->{window} = $self->params->[0] || 1;
    
    return $self;
}
sub run {
    my ($self, $transcript_variation_allele) = @_;
    my $seq = $transcript_variation_allele->variation_feature->feature_Slice->expand($self->{window}, $self->{window})->seq;

    return { Context => $seq };
}

1;