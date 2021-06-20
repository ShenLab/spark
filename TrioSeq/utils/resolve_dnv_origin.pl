use strict;
use warnings;
use FindBin qw|$Bin|;
use Data::Dumper;
use Getopt::Euclid;
use Utils::File::Iter qw|iter_file|; 

use lib "$Bin/../../lib";
use Shared qw|struct_dat parse_tabfile|;


my ($it, $fnames, $keyfields) = parse_tabfile($ARGV{'--tab'}, $ARGV{'--fields'}, 8, 8);
my $PGT = $keyfields->[6];
my $PID = $keyfields->[7];
foreach my $field (qw|FamMembers Relations|, $PGT.".FamMembers", $PID.".FamMembers") {
	unless (grep { $_ eq $field } @$fnames) {
		die "Cannot find $field in input file";
	}
}

# Store informative nearby markers
# Informative: het inherited from one parent, the other is ref geno
my (%dnvs, %nearvars);
while(my $dat = $it->()) {
	my ($iid, $chrom, $pos, $ref, $alt, $varid, $pgt, $pid) = @{$dat}{@$keyfields};
	if (join(":", $chrom, $pos, $ref, $alt) eq $varid) {
		$dnvs{$iid,$varid} = { PGT => $pgt, PID => $pid, MatOrigin => 0, PatOrigin => 0 };
	}	
}

($it, $fnames, $keyfields) = parse_tabfile($ARGV{'--tab'}, $ARGV{'--fields'}, 8, 8);

while(my $dat = $it->()) {
	my $info = struct_dat($dat, { key => "Relations", value => "$PGT.FamMembers" });
	my $geno = $info->{"$PGT.FamMembers"};
	next unless defined $geno->{Mother} && defined $geno->{Father};

	my ($iid, $chrom, $pos, $ref, $alt, $varid, $pgt, $pid) = @{$dat}{@$keyfields};
	next if join(":", $chrom, $pos, $ref, $alt) eq $varid;

	
	unless(defined $dnvs{$iid,$varid}) {
		warn "Cannot find original DNV $varid in the table";
		next;
	}
	# For nearby sites, first determine informativeness
	my $inherit;
	my @altal = grep { $_ ne '0' } split(/[\/\|]/, $pgt);
	if ($geno->{Mother} eq '0/0' && $geno->{Father} ne '0/0') {
		my @alts = grep { $_ ne '0' } split(/[\/\|]/, $geno->{Father});
		if ($altal[0] eq $alts[0]) {
			$inherit = 'Father';
		}
		else {
			warn "Cannot determine inheritance for $chrom:$pos:$ref:$alt";
		}
	}
	elsif ($geno->{Father} eq '0/0' && $geno->{Mother} ne '0/0') {
		my @alts = grep { $_ ne '0' } split(/[\/\|]/, $geno->{Mother});
		if ($altal[0] eq $alts[0]) {
			$inherit = 'Mother';
		}
		else {
			warn "Cannot determine inheritance for $chrom:$pos:$ref:$alt";
		}
	}
	# If marker is informative
	# And if marker is phased on the same haplotype as DNV
	if ($inherit) {
		if ($pid ne '.' && $pid eq $dnvs{$iid,$varid}{PID}) {
			if ($pgt =~ /^0\|\d$/ && $dnvs{$iid,$varid}{PGT} =~ /^0\|\d$/ ||
				$pgt =~ /^\d\|0$/ && $dnvs{$iid,$varid}{PGT} =~ /^\d\|0$/) {
				if ($inherit eq 'Father') {
					$dnvs{$iid,$varid}{PatOrigin} ++
				}
				else {
					$dnvs{$iid,$varid}{MatOrigin} ++
				}
			}
			elsif ($pgt =~ /^0\|\d$/ && $dnvs{$iid,$varid}{PGT} =~ /^\d\|0$/ ||
				$pgt =~ /^\d\|0$/ && $dnvs{$iid,$varid}{PGT} =~ /^0\|\d$/) {
				if ($inherit eq 'Mother') {
					$dnvs{$iid,$varid}{PatOrigin} ++
				}
				else {
					$dnvs{$iid,$varid}{MatOrigin} ++
				}
			}
		}
		#elsif ($dat->{PID} eq '.' && $dnvs{$iid,$varid}{PID} eq '.') {
		#		if ($inherit eq 'Mother') {
		#			$dnvs{$iid,$varid}{PatOrigin} ++
		#		}
		#		else {
		#			$dnvs{$iid,$varid}{MatOrigin} ++
		#		}
		#	}
		#}
	}	
}


my $fout;
if (defined $ARGV{'--output'}) {
	open $fout, ">$ARGV{'--output'}" or die "Cannot write to $ARGV{'--output'}";
} 
else {
	$fout = \*STDOUT;
}
print $fout join("\t", qw|IID VarID Chrom Position Ref Alt ParentOrigin NumInfoVars|), "\n";
foreach my $viid (keys %dnvs) {
	my $dat = $dnvs{$viid};
	my $poo;
	if ($dat->{PatOrigin} > 0 && $dat->{MatOrigin} == 0) {
		$poo = "Father";
	}
	elsif ($dat->{PatOrigin} == 0 && $dat->{MatOrigin} > 0) {
		$poo = "Mother";
	}
	elsif ($dat->{PatOrigin} > 0 && $dat->{MatOrigin} > 0) {
		$poo = "Ambiguous";
	}
	next unless defined $poo;
	my ($iid, $varid) = split($;, $viid);
	my @vardat = split(':', $varid);
	my $nmk = sprintf("%d(%d,%d)", $dat->{PatOrigin}+$dat->{MatOrigin}, $dat->{PatOrigin}, $dat->{MatOrigin});
	print $fout join("\t", $iid, $varid, @vardat, $poo, $nmk), "\n";
}




__END__

=head1 NAME

resolve_dnv_origin.pl -- Infer parent-of-origin for DNVs

=head1 NOTE

The input should contain standard fields for sample nearby variants (see --fields option)
and standard fields for phasing information: PGT PID and those in the family members from 
read backed phasing using WhatsHap.

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

Standard fields for tabulated nearby variant 

=for Euclid:
	string.default: 'IID,Chrom,Position,Ref,Alt,OrigVarID,PGT,PID'


=back