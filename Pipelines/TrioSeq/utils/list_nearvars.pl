use strict;
use warnings;
use FindBin qw|$Bin|;
use Getopt::Euclid;
use Genome::UCSC qw|hg_chrom|;
use List::MoreUtils qw|all|;
use Utils::File::Iter qw|iter_file|;

use lib "$Bin/../../lib";
use Shared qw|expand_dat parse_tabfile|;
use Variants qw|parse_nearby var_dist|;

my ($it, $fnames, $keyfields) = parse_tabfile($ARGV{'--tab'}, $ARGV{'--fields'}, 6, 6);

my @f_extra;
if (defined $ARGV{'--fields-extra'}) {
	@f_extra = split(',', $ARGV{'--fields-extra'});
	foreach my $field (@f_extra) {
		unless(grep { $field eq $_ } @$fnames) {
			die "Cannot find field $field from input file";
		}
	}
}

my (%varclst, %varextra);
while(my $dat = $it->()) {
	my ($iid, $chrom, $pos, $ref, $alt, $nearby) = @{$dat}{@$keyfields};
	next if $nearby eq '.';
	my @nearvars = parse_nearby($nearby, $ARGV{'--qual'});
	next unless @nearvars > 0;

	if ($ARGV{'--chr'} && $chrom !~ /^chr/) {
		$chrom = hg_chrom($chrom);
	}
	
	my @nearout;
	foreach my $nearvar (@nearvars) {
		my $dist = var_dist($nearvar->[1], $nearvar->[2], $pos, $ref);
		if ($dist <= $ARGV{'--dist'} && $dist > 0) {
			#print $fout join("\t", $iid, @$nearvar, $varid, $dist), "\n";
			push @nearout, [$iid, @$nearvar, join(":", @$nearvar), $dist];
		}
	}
	my $varid = join(":", $chrom, $pos, $ref, $alt);
	if (@f_extra) {
		$varextra{$iid,$varid} = join("\t", map { $dat->{$_} // "." } @f_extra);
	}
	if (@nearout) {
		$varclst{$iid}{$varid} = \@nearout;
	}
}

my %clusted;
my $fout;
if ($ARGV{'--output'}) {
	open $fout, ">$ARGV{'--output'}" or die "Cannot write to $ARGV{'--output'}";
}
else {
	$fout = \*STDOUT;
}
print $fout join("\t", qw|IID Chrom Position Ref Alt OrigVarID Distance|, @f_extra), "\n";
foreach my $iid (keys %varclst) {
	foreach my $varid (sort { scalar(@{$varclst{$iid}{$b}}) <=> 
							  scalar(@{$varclst{$iid}{$a}}) } keys %{$varclst{$iid}}) {
		my $cluster = $varclst{$iid}{$varid};
		unless ( (all { defined $clusted{$_->[0],$_->[5]} } @$cluster) && defined $clusted{$iid,$varid} ) {
			print $fout join("\t", $iid, split(':', $varid), $varid, 0);
			if (@f_extra) {
				print $fout $varextra{$iid,$varid} // join("\t", (".") x scalar(@f_extra)), "\n";
			}
			else {
				print $fout "\n";
			}					
			$clusted{$iid,$varid} = 1;
			foreach my $nearvar (sort { $a->[2] <=> $b->[2] } @$cluster) {
				print $fout join("\t", @{$nearvar}[0..4], $varid, $nearvar->[6]);
				if (@f_extra) {
					print $fout $varextra{$nearvar->[0],$nearvar->[5]} // join("\t", (".") x scalar(@f_extra)), "\n";
				}
				else {
					print $fout "\n";
				}	
				$clusted{$nearvar->[0],$nearvar->[5]} = 1;
			}
		}
	}
}



__END__

=head1 NAME

list_nearvars.pl -- Extract nearby variants from DeepVariant evaluation.

=head1 NOTE

The output will contain standard fields for sample variant plus original variant ID (OrigVarID)
of and physical distance to the original variant (Distance). OrigVarID will be in 
Chrom:Position:Ref:Alt format. Original variant will appear first followed by nearby variants
sorted by distance.

The phase information can be extracted from original VCF files (PGT) using var_table. And the phase of 
nearby variants can be used to infer parent of origin of DNVs or to identify MNVs.

=head1 REQUIRED ARGUMENTS
 
=over

=item -[-]tab[le] [=] <input>

The input varinat table.

=back

=head1 OPTIONS

=over


=item -[-]out[put] [=] <file>

The output file name.

=item -[-]fields [=] <string>

Standard fields for sample variant and nearby variants

=for Euclid:
	string.default: 'IID,Chrom,Position,Ref,Alt,NearbyVars'

=item -[-]fields-extra [=] <string>

Extra fields from the original input to be kept in the output.

=item -[-]dist [=] <distance>

Distance (bp) cutoff for nearby variants.

=for Euclid:
	distance.default: 1000

=item -[-]qual [=] <thres>

Qual score threshold for nearby variants.

=for Euclid:
	thres.default: 0

=item -[-]chr

Add chr prefix to chromosome names.

=back