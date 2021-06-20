use strict;
use warnings;
use FindBin qw|$Bin|;
use List::Util qw|sum|;
use Getopt::Euclid;
use Utils::Stat qw|mean|;
use Utils::List qw|parse_fields|;
use Utils::File::Iter;

use lib "$Bin/../../lib";
use Shared qw|parse_fstr|;


my ($it, $fnames) = iter_file($ARGV{'--input'});

# Check standard fields
foreach my $fd (qw|GeneID GeneEff prob|) {
	unless(grep { $fd eq $_ } @$fnames) {
		die "Cannot find standard field $fd from input!";
	}
}

# Determine the fields to be filled in
my @fillin = qw|Chrom|;
foreach my $fd (parse_fields($ARGV{'--fill'})) {
	if (grep { $fd eq $_ } qw|GeneID GeneEff prob|) {
		die "Fill in field $fd is one of standard fields";
	}
	if (grep { $fd eq $_ } @$fnames) {
		push @fillin, $fd;
	}
}
unless (@fillin) {
	warn "No fields to be filled in";
}

my $fout;
if ($ARGV{'--output'}) {
	open $fout, ">$ARGV{'--output'}" or die "Cannot write to $ARGV{'--output'}";
}
else {
	$fout = \*STDOUT;
}
print $fout join("\t", @$fnames), "\n";

# Tally the mutation rate by class stratified by fill in fields 
my (%genert, %genepos);
while(my $dat = $it->()) {
	my $class;
	if (@fillin) {
		$class = join("\t", map { $dat->{$_} } @fillin);
		$genert{$dat->{GeneID}}{$class}{$dat->{GeneEff}} += $dat->{prob};
	}
	else {
		$genert{$dat->{GeneID}}{$dat->{GeneEff}} += $dat->{prob};
	}
	push @{$genepos{$dat->{GeneID}}} => $dat->{Position};
	print $fout join("\t", @{$dat}{@$fnames}), "\n";
}

# Determine how the new class will be derived
my $newcls = parse_fstr($ARGV{'--fix'}, 1);
my %newvarcls;
foreach my $clslab (sort keys %$newcls) {
	my $clsdef = $newcls->{$clslab};
	my (@coeffs, @fields);
	foreach my $term (split('\+', $clsdef)) {
		my @spterm = split('\*', $term);
		if (@spterm == 1) {
			push @coeffs, 1;
			push @fields, $spterm[0];
		}
		elsif (@spterm == 2) {
			unless($spterm[0] =~ /^[0-9\.]+$/) {
				die "Incorrect format of coefficient in combination term: $term";
			}
			push @coeffs, $spterm[0];
			push @fields, $spterm[1];
		}
		else {
			die "Cannot split $term to find coefficient and field";
		}
	}
	$newvarcls{$clslab} = { Coeff => \@coeffs, Fields => \@fields }
}

# Output results for added variant class 
foreach my $gid (sort keys %genert) {
	foreach my $label (sort keys %newvarcls) {
		my @coeffs = @{$newvarcls{$label}{Coeff}};
		my @fields = @{$newvarcls{$label}{Fields}};

		my %data = (GeneID => $gid, GeneEff => $label, Position => int(mean($genepos{$gid})) );
		if (@fillin) {
			foreach my $class (keys %{$genert{$gid}}) {
				my @fillvals = split(/\t/, $class);
				@data{@fillin} = @fillvals;
				$data{prob} = sum(map { $coeffs[$_]*$genert{$gid}{$class}{$fields[$_]} // 0 } 0..$#fields);
				print $fout join("\t", map { $data{$_} // $ARGV{'--nastr'} } @$fnames), "\n";
			}
		}
		else {
			$data{prob} = sum(map { $coeffs[$_]*$genert{$gid}{$fields[$_]} // 0 } 0..$#fields);
			print $fout join("\t", map { $data{$_} // $ARGV{'--nastr'} } @$fnames), "\n";
		}
	}
}




__END__

=NOTE

This script fix the indel rate for per-bp mutation rate table used in DenovoWEST input.

Assuming per-bp allelic SNV prob has been annotated, then:

We add customized cumulative mut rate for non-SNV classes derived from existing SNV class, and fill in other fields 
for determining weights. Other fields that exist in the table will be tallied from original SNV class for the same gene.
Score and other field will all be set to ".".

=head1 REQUIRED ARGUMENTS

=over

=item -[-]in[put] [=] <table>

Input file, must tab-separated and include GeneID,GeneEff,prob!

=item -[-]fix [=] <fstr>

Define new variant class from exising class based on GeneEff and scaling factor.
Example: frameshift:1.3*stop_gained,inframe_indel:0.03*missense

=item -[-]fill [=] <fields>

Fields to be filled in for the new variant class from existing SNV class. Fields that do not exist in the input will be skipped.

=back

=head1 OPTIONS

=over

=item -[-]out[put] [=] <file>

Output file.

=item -[-]nastr [=] <string>

NA string in the output.

=for Euclid:
	string.default: "."

=back





