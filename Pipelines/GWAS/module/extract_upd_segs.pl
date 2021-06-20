use strict;
use warnings;
use List::Util qw|first|;
use List::MoreUtils qw|uniq|;
use Utils::File::Iter qw|iter_file|;
use Genome::UCSC qw|hg_chrom|;
use Utils::Hash qw|peek_hash|;
use Perl6::Slurp;
use Getopt::Euclid;


# Threshold of the fraction of informative events that support UPD event
$ARGV{'--thres'} = 0.99 unless defined $ARGV{'--thres'};

my $it = iter_file($ARGV{'--input'});

my (%updchr, $u_chr);
while(my $dat = $it->()) {
	push @{$updchr{$dat->{IID}}}, $dat->{Chr};
	unless (defined $u_chr) {
		if ($dat->{Chr} =~ /^chr/) {
			$u_chr = 1;
		}
		else {
			$u_chr = 0;
		}
	}
}

my $prefix = $ARGV{'--output'}; $prefix =~ s/\.txt$//;

# Find a region of consecutive UI/UA events
open my $fmd, ">$prefix.README.txt" or die "Cannot write to $prefix.README.txt";
print $fmd <<EOF;
IID 	Sample ID 
Chrom	Chromosome
Start 	Position for the first informative marker in the UPD segment
End		Position for the last informative marker in the UPD segment
POO		Parent of origin for the UPD event (Mat or Pat)
NumMarkers	Total number of informative markers
PropHom		Proportion of homozygotes among informative markers in the UPD segment
FracEvents	Fraction of informative markers in the chromosome supporting the events in the UPD segment
EOF

my (%chrsize, $s_chr);
if ($ARGV{'--chromsize'}) {
	print $fmd "FracChrom	Fraction of chromosome length spanned by the UPD segment\n";
	%chrsize = map { (split)[0,1] } slurp $ARGV{'--chromsize'};
	my $firstchr = peek_hash(\%chrsize);
	if ($firstchr =~ /^chr/) {
		$s_chr = 1;
	}
	else {
		$s_chr = 0;
	}
}

open my $fout, ">$ARGV{'--output'}" or die "Cannot write to $ARGV{'--output'}";
print $fout join("\t", qw|IID Chrom Start End POO NumMarkers PropHom FracEvents|);
if ($ARGV{'--chromsize'}) {
	print $fout "\tFracChrom\n";
}
else {
	print $fout "\n";
}

foreach my $iid (sort keys %updchr) {
	foreach my $chr (uniq sort @{$updchr{$iid}}) {
		#print $iid, "\t", $chr, "\n"; next;
		my ($ct_m, $ct_p) = (0, 0);
		my (@positions, @events, @genos);
		open my $fin, "awk '\$1==$chr && \$8!=\"uninformative\" \{print \$2,\$5,\$8\}' $ARGV{'--wrkdir'}/$iid/$iid.events_list |"
		#open my $fin, "$ARGV{'--wrkdir'}/$iid/$iid.events_list"
			or die "Cannot open awk pipe for reading events_list of $iid";
		while(<$fin>) {
			# my ($pos, $gt, $event) =  (split)[1,4,7];
			my ($pos, $gt, $event) = split;
			#next if $event eq 'uninformative';
			if ($event =~ /^U[AI]_([MP])/) {
				if ($1 eq 'M') {
					$ct_m ++;
					push @events, 'U_M';
				}
				else {
					$ct_p ++;
					push @events, 'U_P';
				}
			}
			else {
				push @events, 'B';
			}
			push @positions, $pos;
			if ($gt eq '01' || $gt eq '10') {
				push @genos, 'Het';
			}
			elsif ($gt eq '00' || $gt eq '11') {
				push @genos, 'Hom';
			}
			else {
				die "Cannot recognize genotype $gt";
			}
		}
		# Start with the first and last U events
		my $st_ii = first { $events[$_] =~ /^U/ } 0..$#events;
		my $ed_ii = first { $events[$_] =~ /^U/ } reverse 0..$#events;
		# Then narrow the range from both ends until it meets the threshold
		while($ed_ii > $st_ii && calc_uperc(\@events, $st_ii, $ed_ii) < $ARGV{'--thres'}) {
			# one step forward from start
			my $new_stii = first { $events[$_] =~ /^U/ } ($st_ii+1)..$ed_ii;
			my $new_len1;
			if ($new_stii < $ed_ii && calc_uperc(\@events, $st_ii, $ed_ii) >= $ARGV{'--thres'}) {
				$new_len1 = $positions[$ed_ii] - $positions[$new_stii];
			}
			my $new_edii = first { $events[$_] =~ /^U/ } reverse $st_ii..($ed_ii-1);
			my $new_len2;
			if ($new_edii > $st_ii && calc_uperc(\@events, $st_ii, $ed_ii) >= $ARGV{'--thres'}) {
				$new_len2 = $positions[$new_edii] - $positions[$st_ii];
			}
			if (defined $new_len1 || defined $new_len2) {
				if ($new_len1 > $new_len2) {
					$st_ii = $new_stii;
				}
				else {
					$ed_ii = $new_edii;
				}
				last;
			}
			else {
				$st_ii = $new_stii;
				$ed_ii = $new_edii;
			}
		}
		if (calc_uperc(\@events, $st_ii, $ed_ii) >= $ARGV{'--thres'}) {
			my $nhom = grep { $_ eq 'Hom' } @genos[$st_ii..$ed_ii];
			my $totevents = grep { $_ =~ /^U/ } @events;
			my $segevents;
			if ($ct_m > $ct_p) {
			 	$segevents = grep { $_ eq 'U_M' } @events[$st_ii..$ed_ii];
			}
			else {
				$segevents = grep { $_ eq 'U_P' } @events[$st_ii..$ed_ii];	
			}
			my @outputs = ($iid, $chr, $positions[$st_ii], $positions[$ed_ii],
							$ct_m > $ct_p ? 'Mat' : 'Pat', $ed_ii-$st_ii+1,
							sprintf("%.3f", $nhom/($ed_ii-$st_ii+1)), sprintf("%.3f", $segevents/$totevents));
			if ($ARGV{'--chromsize'}) {
				my $chrom = $chr;
				if ($u_chr == 0 && $s_chr == 1) {
					$chrom = hg_chrom($chr); 
				}
				elsif ($u_chr == 1 && $s_chr == 0) {
					$chrom =~ s/^chr//;
				}
				my $chrlen = $chrsize{$chrom};
				if (defined $chrlen && $chrlen > 0) {
					push @outputs, sprintf("%.3f", ($positions[$ed_ii]-$positions[$st_ii]+1)/$chrlen);
				}
				else {
					push @outputs, ".";
				}
			}
			print $fout join("\t", @outputs), "\n";
		}
		else {
			print STDERR join("\t", $st_ii, $ed_ii), "\n";
			warn "Cannot find UPD segment on chr$chr for $iid\n";
		}
	}
}

# Calculate percentage of U events
sub calc_uperc {
	my ($events, $st_ii, $ed_ii) = @_;
	my $num_u = grep { $_ =~ /^U/ } @$events;
	return $num_u/($ed_ii-$st_ii+1);
}


__END__

=head1 NAME

extract_upd_segs.pl -- Extract UPD segments from UPDio output.

=head1 NOTES

UPDio reports UPD event for trio offspring chromsome, but it does not give the segments.
Given the final UPD events table, we assume each chromsome contain one UPD event, then search 
for the max segment using greedy algorithm.

=head1 REQUIRED ARGUMENTS
 
=over
 
=item -[-]in[put] [=] <intab>

Table of UPD events.

=item -[-]wrk[dir] [=] <wrkdir>

Working directory of intermediate UPDio results.

=item -[-]out[put] [=] <outfile>

Output file name.

=back
 
=head1 OPTIONS
 
=over

=item -[-]frac [=] <fraction>

Fraction of informative markers consistent with UPD in the final segment.

=for Euclid:
	fraction.default: 0.99

=item -[-]chromsize [=] <chromsize>

Chromosome size table, should have two columns: chromosome name and length. 
If provided, the fraction of chromosome length for each segment will be calculated.

=back



