use strict;
use warnings;
use Perl6::Slurp;
use Getopt::Euclid;

my $ignore = qr/$ARGV{'--ignore'}/; 

my %known = map { $_->[0].$;.$_->[1] => 1 } grep { $_->[1] !~ /$ignore/ } 
	map { my @a = split; [@a] } slurp $ARGV{'--ped'};
#print STDERR scalar(keys %known), "\n";

my %keep;
if ($ARGV{'--keep'}) {
	my $fin;
	if ($ARGV{'--keep'} eq '-') {
		$fin = \*STDIN;
	} 
	else {
		open $fin, $ARGV{'--keep'} or die "Cannot open $ARGV{'--keep'}";
	}
	while(<$fin>) {
		my @a = split;
		if (@a == 1) {
			$keep{$a[0]} = 1;
		}
		elsif (@a == 2) {
			$keep{$a[0],$a[1]} = 1;
		}
		else {
			warn "Incorrect number of columns in sample list: $_";
		}
	}
}

my $fout;
if ($ARGV{'--output'}) {
	open $fout, ">$ARGV{'--output'}" or die "Cannot write to $ARGV{'--output'}";
}
else {
	$fout = \*STDOUT;
}
open my $fin, $ARGV{'--ped'} or die "Cannot open $ARGV{'--ped'} for reading";
print $fout join("\t", qw|FamID Child Father Mother Sex Pheno|), "\n";
while(<$fin>) {
	my ($fid, $iid, $father, $mother, $sex, $aff) = split;
	next unless $father ne '0' && $mother ne '0';
	next unless defined $known{$fid,$father} && defined $known{$fid,$mother};
	if ($ARGV{'--keep'}) {
		next unless (defined $keep{$iid} || defined $keep{$fid,$iid}) &&
				(defined $keep{$father} || defined $keep{$fid,$father}) && 
				(defined $keep{$mother} || defined $keep{$fid,$mother});
	}
	my $gender;
	if ($sex eq '1') {
		$gender = "Male";
	}
	elsif ($sex eq '2') {
		$gender = "Female";
	}
	else {
		warn "Cannot determine gender of $fid,$iid";
		$gender = "Unknown";
	}
	my $pheno;
	if ($aff eq '1') {
		$pheno = "Unaffected";
	}
	elsif ($aff eq '2') {
		$pheno = "Affected";
	}
	else {
		$pheno = "Unknown";
	}
	print $fout join("\t", $fid, $iid, $father, $mother, $gender, $pheno), "\n";
}

__END__

=head1 NAME

list_trios.pl -- List all trios from PED file.

=head1 NOTE

Default is to extract all complete trios, which require that sample in each trio 
to appear in 2nd column of PED file.

=head1 REQUIRED ARGUMENTS

=over 

=item -[-]ped [=] <pedfile>

Pedigree file.

=for Euclid:
	pedfile.type: readable

=back

=head1 OPTIONS

=over

=item -[-]out[put] <outfile>

The output file name, write to STDOUT if not providd.

=item -[-]keep [=] <list>

A list of samples to be kept. The intersection of this list with PED file will be used to 
determine the final list of samples from which trios are extracted.
The list can be one or two columns, and can also be redirected from STDIN.

=item -[-]ignore [=] <string>

Regex pattern for duplicated sample IDs.

=for Euclid:
	string.default: '_Re\d*$'

=back


