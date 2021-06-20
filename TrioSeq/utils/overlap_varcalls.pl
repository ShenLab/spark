use strict;
use warnings;
use Data::Dumper;
use FindBin qw|$Bin|;
use File::Basename;
use Getopt::Euclid;
use Genome::UCSC qw|hg_chrom|;
use Genome::UCSC::BinKeeper;

use lib "$Bin/../../lib";
use Shared qw|parse_tabfile|;

# Store all variants found from variant list
my %vars;
my %nearbk;

my @varlist;
if ($ARGV{'--list-fields'}) {
	@varlist = ([ $ARGV{'--list'},  $ARGV{'--list-label'} // (split(/\./, basename($ARGV{'--list'})))[0], $ARGV{'--list-fields'} ])
}
else {
	open my $fin, $ARGV{'--list'} or die "Cannot open $ARGV{'--list'}";
	while(<$fin>) {
		my ($vartab, $label, $fields) = split;
		push @varlist, [$vartab, $label, $fields];
	}
}


my %dups;
if ($ARGV{'--dups'}) {
	open my $fin, $ARGV{'--dups'} or die "Cannot open $ARGV{'--dups'}";
	while(<$fin>) {
		my @a = split;
		$dups{$a[0]}{$a[1]} = 1;
		$dups{$a[1]}{$a[0]} = 1;
	}
}

my %known;
my ($it, $fnames, $keyfields) = parse_tabfile($ARGV{'--input'}, $ARGV{'--fields'}, 5);
my $t_chr;
{
	my $dat = $it->();
	my ($iid, $chrom, $pos, $ref, $alt) = @{$dat}{@$keyfields}; 
	if ($chrom =~ /^chr/) {
		$t_chr = 1;
	}
	else {
		$t_chr = 0;
	}
}

my $ct = 1;
foreach my $varlst (@varlist) {
	my ($vartab, $label, $fields) = @$varlst;
	unless(-f $vartab) {
		die "Cannot find variant table $vartab";
	}
	unless(defined $label) {
		$label = (split(/\./, basename($vartab)))[0];
	}
	if (!defined $fields || $fields eq "") {
		$fields = $ARGV{'--fields'};
	}
	$known{$label} = $ct ++ unless defined $known{$label};
	#$known{$label} = 1;
	if (grep { $label eq $_ } @$fnames) {
		die "Label $label already exists in input file!";
	}
	my %knownvar;
	print STDERR "Slurping variants from $vartab\n";
	my ($it, $fnames, $keyfields) = parse_tabfile($vartab, $fields, 5, 5);
	my $l_chr;
	while(my $dat = $it->()) {
		my ($iid, $chrom, $pos, $ref, $alt) = @{$dat}{@$keyfields};
		unless(defined $l_chr) {
			$l_chr = $chrom =~ /^chr/ ? 1 : 0; 
			if ($l_chr != $t_chr) {
				warn "Chromosome nomenclature in $vartab is different from input";
			}
		}
		# Change chr code to match the input
		if ($l_chr == 1 && $t_chr == 0) {
			$chrom =~ s/^chr//; $chrom = 'MT' if $chrom eq 'M';
		}
		elsif ($l_chr == 0 && $t_chr == 1) {
			$chrom = hg_chrom($chrom);	
		}

		my $varid = join(":", $chrom, $pos, $ref, $alt);
		next if $knownvar{$iid,$varid};
		$vars{$iid}{$varid}{$label} = 1;
		unless(defined $nearbk{$iid}) {
			$nearbk{$iid} = Genome::UCSC::BinKeeper->new(); 
		}
		$nearbk{$iid}->add($chrom, $pos-$ARGV{'--nearby'}, $pos+length($ref)-1+$ARGV{'--nearby'}, $varid, $label);
		$knownvar{$iid,$varid} = 1;
	}
}

my $fout;
if ($ARGV{'--output'}) {
	open $fout, ">$ARGV{'--output'}" or die "Cannot write to $ARGV{'--output'}";
}
else {
	$fout = \*STDOUT;
}

my @labs = sort { $known{$a} <=> $known{$b} } keys %known;
#my @labs = sort keys %known;
my @outfields = @$fnames;
unshift @outfields, @labs;
print $fout join("\t", @outfields), "\n";

($it, $fnames, $keyfields) = parse_tabfile($ARGV{'--input'} , $ARGV{'--fields'}, 5);
while(my $dat = $it->()) {
	my ($iid, $chrom, $pos, $ref, $alt) = @{$dat}{@$keyfields}; 
	my $varid = join(":", $chrom, $pos, $ref, $alt);
	
	my %nearby;
	if (defined $nearbk{$iid}) {
		$nearby{$iid} = [ $nearbk{$iid}->find_range($chrom, $pos, $pos+length($ref)-1) ];
	}
	if (defined $dups{$iid}) {
		foreach my $alias (keys %{$dups{$iid}}) {
			if (defined $nearbk{$alias}) {
				$nearby{$alias} = [ $nearbk{$alias}->find_range($chrom, $pos, $pos+length($ref)-1) ];
			}
		}
	}
	foreach my $lab (@labs) {
		if (defined $vars{$iid}{$varid} && defined $vars{$iid}{$varid}{$lab}) {
			$dat->{$lab} = "x";
		}
		else {
			my @nearvars = grep { $_->[3] eq $lab } @{$nearby{$iid}};			
			if (@nearvars > 0) {
				$dat->{$lab} = join(',', map { $_->[2] } @nearvars);
			}
			else {
				$dat->{$lab} = ".";
			}
		}
		if ($dat->{$lab} eq '.' && defined $dups{$iid}) {
			foreach my $alias (keys %{$dups{$iid}}) {
				if (defined $vars{$alias}{$varid} && defined $vars{$alias}{$varid}{$lab}) {
					$dat->{$lab} = "[x]";
					last;
				}
				else {
					my @nearvars = grep { $_->[3] eq $lab } @{$nearby{$alias}};			
					if (@nearvars > 0) {
						$dat->{$lab} = "[".join(',', map { $_->[2] } @nearvars)."]";
						last;
					}
				}
			}
		}
	}
	print $fout join("\t", @{$dat}{@outfields}), "\n";
}




__END__

=head1 NAME

overlap_varcalls.pl -- Markup the overlapping variant calls. 

=head1 NOTE

For SNVs, exact match is required. For MNVs, we also look for overlapping variants.
And for indels, we also look for nearby or overlapping variants from this distance cutoff.

M: exact match of Chrom,Pos,Ref,Alt in the same individual, 
O: overlap variants (for indels & MNVs), N: nearby variants (for indels) 

We will store the entire variant from the list to binKeeper, so this script can only
work for small number of variants.

=head1 REQUIRED ARGUMENTS

=over

=item -[-]in[put] [=] <vartab>

The input variant table.

=item -[-]list [=] <list>

A list of known variant calls, format: VarTab, Label, [Fields]
Label will be used as additional column names in the output.
Optional fields list can be used to specify the key columns in those variant tables.

The list can also be a single variant table. In such case, the keys fields must be provided
by --list-fields option and the base file name as the output column name.

=back

=head1 OPTIONS

=over

=item -[-]out[put] <file>

Output file with additional columns.

=item -[-]fields [=] <fstr>

Key fields in the input table. Default: IID,Chrom,Position,Ref,Alt

=for Euclid:
	fstr.default: "IID,Chrom,Position,Ref,Alt"

=item -[-]list-fields [=] <fstr>

Key fields from the variant list file.

=item -[-]list-label [=] <string>

Label for the variant list file.

=item -[-]dups [=] <list>

A list of duplicates. So the variant matching will also be done between duplicates. 
Priority will be given to variants found in the same individual.

=item -[-]nearby [=] <distance>

Distance cutoff for nearby indels.

=for Euclid:
	distance.default: 10

=back

