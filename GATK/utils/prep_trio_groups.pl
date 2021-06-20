use strict;
use warnings;
use Perl6::Slurp;
use Getopt::Euclid;
use List::MoreUtils qw|any all|;
use Genet::Ped;

my $dupsamp = qr/$ARGV{'--ignore'}/;

my $ped = Genet::Ped->new($ARGV{'--ped'}, { ignore => $dupsamp, verbose => 0 });

# Also identify all Dup samples
my %dups;
{
	open my $fin, $ARGV{'--ped'} or die "Cannot open PED file";
	while(<$fin>) {
		my $iid = (split)[1];
		if ($iid =~ /$dupsamp/) {
			(my $origid = $iid) =~ s/$dupsamp//;
			push @{$dups{$origid}}, $iid;
		}
	}
}

my %alltrios = $ped->get_trios();

my (%keep, %remove);
if ($ARGV{'--keep'}) {
	if ($ARGV{'--strict'}) {
		my %incl = map { (split)[0] => 1 } slurp $ARGV{'--keep'};
		%keep = map { $_ => 1 } grep { defined $incl{$_} } map { (split)[1] } slurp $ARGV{'--ped'};
		unless(%keep) {
			print STDERR "No sample left for inclusion!";
			exit 1;
		}
	}
	else {
		%keep = map { (split)[0] => 1 } slurp $ARGV{'--keep'};
	}
}
else {
	if ($ARGV{'--strict'}) {
		%keep = map { (split)[1] => 1 } slurp $ARGV{'--ped'};
	}
}
if ($ARGV{'--remove'}) {
	%remove = map { (split)[0] => 1 } slurp $ARGV{'--remove'};
}

my $fout;
if ($ARGV{'--output'}) {
	open $fout, ">$ARGV{'--output'}" or die "Cannot write to group list";
}
else {
	$fout = \*STDOUT;
}

my %known;
while(my ($famid, $trios) = each %alltrios) {
	foreach my $trio (@$trios) {
		if (defined $known{$trio->[0]}) {
			die "Offspring $trio->[0] has already been included in a trio group"; 
		}
		if ($ARGV{'--caseonly'}) {
			my $phe = $ped->get_aff($famid, $trio->[0]);
			next unless $phe == 2;
		}
		if (%keep) {
			unless (all { defined $keep{$_} } @$trio) {
				next;
			}
		}
		if (%remove) {
			if (any { defined $remove{$_} } @$trio) {
				next;
			}
		}

		for(my $ii = 0; $ii <= 2; $ii ++) {
			print $fout $trio->[$ii], "\t", $trio->[0], "\n";
			if (defined $dups{$trio->[$ii]}) {
				foreach my $iid (@{$dups{$trio->[$ii]}}) {
					print $fout $iid, "\t", $trio->[0], "\n";
				}
			}
		}
	}
}


__END__

=head1 NAME

prep_trio_groups.pl -- Prepare a sample group file that associate each trio with offspring ID.

=head1 NOTE

Trio group will use child ID for group ID.

=head1 REQUIRED ARGUMENTS
 
=over
 
=item -[-]ped [=] <pedfile>

Pedigree file.

=back

=head1 OPTIONS

=over

=item -[-]remove [=] <samplist>

A list of bad sample to be removed.

=item -[-]keep [=] <samplist>

A list of good sample to be kept.

=item -[-]out[put] [=] <outfile>

Output group list file.

=item -[-]caseonly

Only output case trios (affected child).

=item -[-]strict

Under strict mode, only output a trio when all of its sampes are present in the PED file. 

=item -[-]ignore [=] <pattern>

Regex pattern for duplicated sample. Duplicates will be renamed and added
to the callable region of original sample.

=for Euclid:
	pattern.default: '_Re\d*$'

=back
