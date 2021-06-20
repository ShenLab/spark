use strict;
use warnings;
use Data::Dumper;
use List::MoreUtils qw|any uniq|;
use Utils::File::Iter qw|iter_file|;
use FindBin qw|$Bin|;
use FaSlice;
use Genome::Ranges qw|range_overlap|;
use Genome::UCSC::BinKeeper;
use Getopt::Euclid;

use lib "$Bin/../../lib";
use Variants qw|parse_nearby cluster_vars|;
use Shared qw|parse_tabfile|;


#my ($it, $fnames) = iter_file($ARGV{'--dvout'}, { fsep => qr/\t/ });
my ($it, $fnames, $infields) = parse_tabfile($ARGV{'--dvout'}, $ARGV{'--fields'}, 7, 7);


# Look for adjacent indels that missed by deepvariant
# They are likely new MNV candidates
my %known;
my $bk = Genome::UCSC::BinKeeper->new();
while(my $dat = $it->()) {
	my ($iid, $chrom, $pos, $ref, $alt, $dvfilter, $nearstr) = @{$dat}{@$infields};
	next if $nearstr eq '.';
	if ($dvfilter eq '.' && (length($ref) > 1 || length($alt) > 1) ) {
		foreach my $var (parse_nearby($nearstr)) {
			my $varid = join(":", @{$var}[0..3]);
			next if $known{$iid,$varid};
			$bk->add($var->[0], $var->[1]-$ARGV{'--pad'}, $var->[1]+length($var->[2])-1+$ARGV{'--pad'},
					{ IID => $iid, Chrom => $var->[0], Position => $var->[1], 
					  Ref => $var->[2], Alt => $var->[3], DeepvarQual => $var->[4], 
					  VarID => $varid  });
			$known{$iid,$varid} = 1;
		}
	}
}

# Second round iter, looking for candidates
($it, $fnames, $infields) = parse_tabfile($ARGV{'--dvout'}, $ARGV{'--fields'}, 7, 7);

my @cands;
$it = iter_file($ARGV{'--dvout'}, { fsep => qr/\t/ });
print join("\t", @$fnames), "\n";
while(my $dat = $it->()) {
	my ($iid, $chrom, $pos, $ref, $alt, $dvfilter, $nearstr) = @{$dat}{@$infields};
	next unless $dvfilter eq '.';
	unless(defined $dat->{VarID}) {
		$dat->{VarID} = join(":", $chrom, $pos, $ref, $alt);
	}
	my @cluster = grep { $_->{VarID} ne $dat->{VarID} && $_->{IID} eq $iid } 
			map { $_->[2] } $bk->find_range($chrom, $pos, $pos+length($ref)-1);
	my @types = uniq sort map { var_type($_) } @cluster;
	if (@cluster >= 2 && @types >= $ARGV{'--ntype'} && 
		length($ref) <= $ARGV{'--length'} && length($alt) <= $ARGV{'--length'}) {
		push @cands, \@cluster;
		print join("\t", @{$dat}{@$fnames}), "\n";
	}
}

unless (@cands) {
	print STDERR "No MNV candidates were found\n";
	exit 1;
}

if ($ARGV{'--output'}) {
	unless(defined $ARGV{'--fasta'}) {
		die "Must provide genome ref seq file when output MNV candidates";
	}
	print STDERR "Writing MNV candidates to $ARGV{'--output'}\n";
	open my $fout, ">$ARGV{'--output'}" or die "Cannot write to $ARGV{'--output'}";
	my $refseq = FaSlice->new(file => $ARGV{'--fasta'});
	my @fields = qw|IID VarID Chrom Position Ref Alt DeepvarQual|;
	print $fout join("\t", @fields), "\n";
	foreach my $cluster (@cands) {
		my @clstord = sort { $a->{Chrom} cmp $b->{Chrom} || $a->{Position} <=> $b->{Position} } @$cluster;
		my $iid = $clstord[0]{IID};
		foreach my $var (@clstord) {
			print $fout join("\t", @{$var}{@fields}), "\n";
		}
		for(my $len = 2; $len <= @clstord; $len ++) {
			for(my $ii = 0; $ii < scalar(@clstord)-$len+1; $ii ++) {
				my @mnv = cluster_vars($refseq, map { [ @{$_}{qw|Chrom Position Ref Alt|} ] } @clstord[$ii..$ii+$len-1]);
				unshift @mnv, $iid, join(":", @mnv);
				push @mnv, ".";
				print $fout join("\t", @mnv), "\n";
			}
		}
	}
}

sub var_type {
	my ($var) = @_;
	if (length($var->{Ref})>1 && length($var->{Alt})==1) {
		return "Del";
	}
	elsif (length($var->{Ref})==1 && length($var->{Alt})>1) {
		return "Dup";
	}
	elsif (length($var->{Ref})==1 && length($var->{Alt})==1) {
		return "SNV";
	}
	else {
		return "MNV";
	}
}

__END__

=head1 NAME

find_extra_mnvcands.pl -- Find additional MNVs from DeepVar output.

=head1 NOTE

We look for indels that are not supported DeepVar but has multiple overlapping variants.

=head1 REQUIRED ARGUMENTS
 
=over

=item -[-]dv[out] [=] <input>

The DeepVar evaluation output.


=back

=head1 OPTIONS

=over

=item -[-]fields [=] <string>

Required fields for IID,Chrom,Position,Ref,Alt,DeepvarFilter,NearbyVars

=for Euclid:
	string.default: "IID,Chrom,Position,Ref,Alt,DeepvarFilter,NearbyVars"

=item -[-]out[put] [=] <file>

Output clustered variants for annotation.

=item -[-]fa[sta] [=] <file>

Path to genome reference sequence file in fasta format.
Must be provided when --output is present.

=item -[-]pad [=] <dist>

Distance padded to variants when testing for overlap, default: 1.

=for Euclid:
	dist.default: 1

=item -[-]ntype [=] <num>

Min number of overlapping variant types, default: 2.

=for Euclid:
	num.default: 2

=item -[-]len[gth] [=] <bp>

Max length of the original variant.

=for Euclid:
	bp.default: 10

=cut
