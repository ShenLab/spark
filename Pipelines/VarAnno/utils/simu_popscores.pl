use strict;
use warnings;
use List::MoreUtils qw|all|;
use Config::Std;
use FindBin qw|$Bin|;
use IO::Dir;
use Data::Dumper;
use List::Util qw|min max|;
use Hash::Util qw|lock_keys|;
use Getopt::Euclid;
use Math::Random qw|random_set_seed random_binomial|;
use Utils::Parser qw|sql_query|;
use Utils::File::Iter qw|iter_file|;

use lib "$Bin/../../lib";
use Shared qw|parse_tabfile slurp_xref|;


my $infile = $ARGV{'--input'};
my $outfile = $ARGV{'--output'};
my $config = $ARGV{'--conf'};

read_config $config => my %conf;
lock_keys(%conf);


# Input variant table
my ($it, $fnames, $keyfields) = parse_tabfile($infile, $conf{Input}{Fields}, 2, 2);

# Variant table filter
my $filter;
if (defined $conf{Input}{Filter}) {
	$conf{Input}{Filter} =~ s/^['"]//;
	$conf{Input}{Filter} =~ s/['"]$//;
	($filter, my $tokens) = sql_query($conf{Input}{Filter}, 1);
	foreach my $tok (@$tokens) {
		if ($tok->[0] eq 'FIELD') {
			unless(grep { $tok->[1] eq $_ } @$fnames) {
				die "The field $tok->[1] in filter expression cannot be found in the input";
			}
		}
	}
}
# Also check standard fields
foreach my $field (qw|Chrom Position Ref Alt GeneID|) {
	unless(grep { $_ eq $field } @$fnames) {
		die "Cannot find required field $field in the input file";
	}
}

if (defined $conf{Param}{Seed}) {
	random_set_seed($conf{Param}{Seed}, $conf{Param}{Seed});
}

# Observed patho-scores, allele freq
# Only support one gene in a file!
my ($geneid, @pathoscores, @allelefreqs);
my %known;
while(my $dat = $it->()) {
	unless(defined $geneid) {
		$geneid = $dat->{GeneID};
	}
	else {
		die "Can only have one gene in a table file!" unless $geneid eq $dat->{GeneID};
	}
	my $varid = join(":", @{$dat}{qw|Chrom Position Ref Alt|});
	next if defined $known{$varid};
	$known{$varid} = 1;
	my ($pathoscore, $allelefreq) = @{$dat}{@$keyfields};
	next if $pathoscore eq '.';
	next unless $filter->($dat);
	push @pathoscores => $pathoscore;
	if ($allelefreq ne '.') {
		push @allelefreqs, $allelefreq;
	}
	else {
		push @allelefreqs, 0;
	}
}
# Now adding imputed AF
for(my $ii = 0; $ii < @allelefreqs; $ii ++) {
	if ($allelefreqs[$ii] == 0) {
		$allelefreqs[$ii] = $conf{Param}{ImpAF};	
	}
}

my @cutoffs;
foreach my $score (split(',', $conf{Param}{CutOffs})) {
	if ($score =~ /:/) {
		my ($start, $end, $step) = split(':', $score);
		unless(defined $start && defined $end && defined $step) {
			die "Cannot parse Start:End:Step from $score";
		}
		my $nsteps = int(($end-$start)/$step);
		for(my $ii = 0; $ii <= $nsteps; $ii ++) {
			push @cutoffs, $start+$ii*$step;
		}
	}
	else {
		push @cutoffs, $score;
	}
}
for(my $ii = 0; $ii < @cutoffs-1; $ii ++) {
	unless($cutoffs[$ii+1] > $cutoffs[$ii]) {
		die "Cutoffs are not ordered!";
	}
}
#print Dumper \@cutoffs;

my @hetmaxscores = ($cutoffs[0]) x $conf{Param}{NumSim};
my @het2ndmaxscores = ($cutoffs[0]) x $conf{Param}{NumSim};
my @hommaxscores = ($cutoffs[0]) x $conf{Param}{NumSim};

for(my $ii = 0; $ii < @allelefreqs; $ii ++) {
	print "Similating genotypes at site $ii\n";
	if ($pathoscores[$ii] < $cutoffs[0] || $pathoscores[$ii] >= $cutoffs[-1]) {
		die "Pathogenic score $pathoscores[$ii] is out of allowable range!";
	}
	my @genos = random_binomial($conf{Param}{NumSim}, 2, $allelefreqs[$ii]);
	for(my $jj = 0; $jj < @genos; $jj ++) {
		if ($genos[$jj] == 1) {
			if ($pathoscores[$ii] > $hetmaxscores[$jj]) {
				$het2ndmaxscores[$jj] = $hetmaxscores[$jj];
				$hetmaxscores[$jj] = $pathoscores[$ii];
			}
			elsif ($pathoscores[$ii] > $het2ndmaxscores[$jj]) {
				$het2ndmaxscores[$jj] = $pathoscores[$ii];
			}
		}
		elsif ($genos[$jj] == 2) {
			$hommaxscores[$jj] = $pathoscores[$ii];
		}
	}
}


my ($xtraginfo, @ginfo, @gfields);
if (defined $conf{Output}{GXref}) {
	my ($xgdat, $xgfds) = slurp_xref($conf{Output}{GXref}, $conf{Output}{GXref_Fields});
	$xtraginfo = $xgdat;
	@gfields = @$xgfds;
	foreach my $gfield (@gfields) {
		if (grep { $gfield eq $_ } qw|GeneID PopScore_Het PopScore_cHet PopScore_Hom|) {
			die "Field $gfield already appear in the genetab output!";
		}
	}
	@ginfo = map { $xtraginfo->{$geneid}{$_} // "." } @gfields;  
}

open my $fout, ">$outfile" or die "Cannot write to $outfile";
print $fout join("\t", "GeneID", @gfields, $keyfields->[0], qw|PopScore_Het PopScore_cHet PopScore_Hom|), "\n";
# Calculate pop-scores for each cut-off under three models
my %popscore;
foreach my $cutoff (@cutoffs) {
	my $n_extreme = grep { $_ >= $cutoff } @hetmaxscores;
	push @{$popscore{Het}} => $n_extreme/$conf{Param}{NumSim};
	$n_extreme = grep { $_ >= $cutoff } @het2ndmaxscores;
	push @{$popscore{cHet}} => $n_extreme/$conf{Param}{NumSim};
	$n_extreme = grep { $_ >= $cutoff } @hommaxscores;
	push @{$popscore{Hom}} => $n_extreme/$conf{Param}{NumSim};
}
#print Dumper \%popscore;
print $fout join("\t", $geneid, @ginfo, $conf{Param}{CutOffs}, map { join(",", @{$popscore{$_}}) } qw|Het cHet Hom|), "\n";


__END__

=head1 NAME

simu_popscore.pl -- Simulate genotypes to calculate PopScore for one gene.  

=head1 REQUIRED ARGUMENTS

=over 

=item -[-]in[put] [=] <table>

Input variant table.

=item -[-]conf [=] <config>

Config file.

=item -[-]out[put] [=] <outfile>

Output file.

=back

=back
