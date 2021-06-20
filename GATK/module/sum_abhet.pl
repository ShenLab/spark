use strict;
use warnings;
use Carp;
use IO::Dir;
use File::Temp;
use List::Util qw|sum|;
use Getopt::Euclid;
use Utils::Stat qw|std|;
use Utils::File::Iter qw|iter_file|;


my $indir = $ARGV{'--indir'}; $indir =~ s/\/$//;
my @abfiles = grep { /\.ab\.txt$/ } IO::Dir->new($indir)->read();

my %varcts;
my $iter = iter_file("$indir/varcounts.txt");
while(my $dat = $iter->()) {
    $varcts{$dat->{IID}} = { NHETS => $dat->{HET}, NHOMS => $dat->{HOM},
    PHOMS => sprintf("%.3f", $dat->{HOM}/($dat->{HET}+$dat->{HOM})) };
}


my @headfields = qw|IID NumHets NumHoms PropHoms MeanAB MedianAB StdevAB.Obs StdevAB.Exp StdevAB.OvE
                    ABPlus ABProp.Lower ABProp.Upper NumComps MixProps MixMeans MixVars|;

my $fout = IO::File->new($ARGV{'--output'}, "w");
print $fout join("\t", @headfields), "\n";
foreach my $abfile (@abfiles) {
    (my $iid = $abfile) =~ s/\.ab\.txt$//;
    my $stats = parse_stats_file("$indir/$iid.stats.summary");
    my ($ncomp, $params) = parse_model_file("$indir/$iid.best_model.summary");
    my ($exp_var_m, $exp_std_m) = calc_exp_var("$indir/$iid.ab.txt", $stats->{MEDIAN});

    my ($abplus, $ablo, $abhi) = (0, undef, undef);
    my @ablos = map { $_->{PROP} } grep { $_->{MEAN} <= $ARGV{'--lower'} } @$params;
    my @abhis = map { $_->{PROP} } grep { $_->{MEAN} >= $ARGV{'--upper'} } @$params;
    my ($lotail, $hitail);
    if (@ablos > 0 || @abhis > 0) {
        $abplus = 1;
        if (@ablos) {
            $ablo = sprintf("%.3f", sum(@ablos));
        }
        if (@abhis) {
            $abhi = sprintf("%.3f", sum(@abhis));
        }
    }

    print $fout join("\t", $iid, @{$varcts{$iid}}{qw|NHETS NHOMS PHOMS|},
                    (map { sprintf("%.3f", $_) }  (@{$stats}{qw|MEAN MEDIAN STDEV|}, 
                                                 $exp_std_m, $stats->{STDEV}/$exp_std_m)),
				    $abplus, $ablo // ".", $abhi // ".", $ncomp,
				    join(",", map { sprintf("%.3f", $_->{PROP}) } @$params),
				    join(",", map { sprintf("%.3f", $_->{MEAN}) } @$params),
				    join(",", map { sprintf("%.3f", $_->{VAR})  } @$params)
				  ), "\n";
}

# Calculate the variance assuming fixed reference bias
sub calc_exp_var {
    my ($file, $pnull) = @_;
    my (%depths, %dpfrac);
    my $it = iter_file($file);
    while(my $dat = $it->()){
       $depths{$dat->{DEPTH}} ++;
   }
   my $ntot = sum(values %depths);
   foreach my $dp (keys %depths) {
       $dpfrac{$dp} = $depths{$dp} / $ntot;
   }
   my $var_exp;
   foreach my $dp (keys %dpfrac) {
       $var_exp += $dpfrac{$dp} * ($pnull*(1-$pnull)/$dp);
   }
   my $std_exp = sqrt($var_exp);
   return ($var_exp, $std_exp);
}

sub parse_model_file {
    my ($file) = @_;
    my $fin = IO::File->new($file) or croak "Cannot open $file";
    my $header = <$fin>;
    my ($ncomp) = ($header =~ /NCOMP=(\d+)$/);
    <$fin>;
    my @params;
    while(<$fin>) {
        my ($prop, $mean, $var) = split;
	   push @params, { PROP => $prop, MEAN => $mean, VAR => $var };
    }
    return ($ncomp, \@params);
}

sub parse_stats_file {
    my ($file) = @_;
    my $fin = IO::File->new($file) or croak "Cannot open $file";
    my $header = <$fin>; my @fields = split(/\s+/, $header);
    my $line = <$fin>;   my @values = split(/\s+/, $line);
    my %stats;
    @stats{@fields} = @values;
    return \%stats;
}


__END__

=head1 NAME

sum_abhet.pl -- Summarize the results of AB test.

=head1 REQUIRED ARGUMENTS

=over

=item -[-]in[dir] [=] <indir>

Input directory. Output of ab_het.pl

=item -[-]out[put] [=] <outfile>

Output file name.

=back

=head1 OPTIONS

=over

=item -[-]lower [=] <threshold>

The lower boundary for contamination peak.

=for Euclid:
  threshold.default: 0.25

=item -[-]upper [=] <threshold>

The upper boundary for contamination peak.

=for Euclid:
  threshold.default: 0.8

=back



