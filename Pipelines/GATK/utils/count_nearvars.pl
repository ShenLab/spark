use strict;
use warnings;
use Data::Dumper;
use Utils::File::Iter qw|iter_file|;
use Getopt::Euclid;

# Assuming Chrom,Position,Ref,Alt are present in the table
# only need to specify the input table and nearbyvar field
# Also numerical code the DeepvarFilter flag
# PASS =. 1, RefCall or . => 0

# Format: 1_861295_C_T[Q57],...
my $nearvar = $ARGV{'--nearvar'};
my ($it, $fnames) = iter_file($ARGV{'--input'}, { fsep => qr/\t/ });
unless(grep { $nearvar eq $_ } @$fnames) {
    die "Cannot find nearbyvar field $nearvar from input table $ARGV{'--input'}";
}

my $fout;
if ($ARGV{'--output'}) {
    open $fout, ">$ARGV{'--output'}" or die "Cannot write to $ARGV{'--output'}";
}
else {
    $fout = \*STDOUT;
}
print $fout join("\t", @$fnames), "\n";

while(my $dat = $it->()) {
    if ($dat->{DeepvarFilter} =~ /PASS/) {
	   $dat->{DeepvarFilter} = 1;
    }
    #elsif ($dat->{DeepvarFilter} =~ /RefCall/) {
    #$dat->{DeepvarFilter} = 0.5;
    #}
    else {
	   $dat->{DeepvarFilter} = 0;
    }
    # Count number of nearby variants
    if ($dat->{$nearvar} eq '.') {
	   $dat->{$nearvar} = 0;
    }
    else {
       	my @nearvars = grep {  $_->[1] > $dat->{Position}-$ARGV{'--cutoff'} &&
			                 $_->[1] < $dat->{Position}+length($dat->{Ref})-1+$ARGV{'--cutoff'}  } 
                             parse_nearvars($dat->{$nearvar});
	   $dat->{$nearvar} = scalar(@nearvars);
    }
    print $fout join("\t", @{$dat}{@$fnames}), "\n";
}


sub parse_nearvars {
    my ($nearvar) = @_;
    my @vars;
    foreach my $var (split(',', $nearvar)) {
	# Format: 1_69511_A_G[Q58]
	$var =~ s/\[Q\d+\]$//;
	my ($chr, $pos, $ref, $alt) = split('_', $var);
	push @vars, [$chr, $pos, $ref, $alt];
    }
    return @vars;
}

__END__

=head1 NAME

count_nearvars.pl -- Count number of nearby variants.

=head1 NOTE

The script will change the NearbyVar field into the number of nearby variants of given distance.

=head1 REQUIRED ARGUMENTS

=over

=item -[-]in[put] [=] <infile>

The input table name.

=back

=head1 OPTIONS

=over

=item -[-]out[put] [=] <outfile>

The output file name (optional).

=item -[-]nearvar [=] <field>

=for Euclid:
    field.default: "NearbyVars"

=item -[-]cutoff [=] <distance>

=for Euclid:
    distance.default: 50

=back
