package Genet::File::VCF;

use strict;
use warnings;


sub custom_site {
	my ($site) = @_;

	my @infofields = sort keys %{$site->{INFO}};
	foreach my $field (@infofields) {
		if ($field =~ /\-/) {
			(my $fieldrename = $field) =~ s/\-/_/g;
			if (defined $site->{INFO}{$fieldrename}) {
				warn "Renamed field name $fieldrename already exist!";
			}
			else {
				$site->{INFO}{$fieldrename} = $site->{INFO}{$field};
			}
		}
	}

	return 1;
}



1;
