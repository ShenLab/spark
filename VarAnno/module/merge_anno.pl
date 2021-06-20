use strict;
use warnings;
use Utils::File::Iter;

# Merge one or more anno_seqvars.pl output:
# Take the first as the primary table, then adding new fields from subsequent tables.
# All input/output will be tab separated.

my $output = pop @ARGV;

unless(@ARGV > 1) {
	die "Merging requires at least two input tables!";
}


my ($iter, @outfields, @extrafields);

# Slurp extra fields from secondary tables
my (%known, %extrainfo);
for(my $ii = 0; $ii < @ARGV; $ii ++) {
	my $table = $ARGV[$ii];
	my ($it, $fnames) = iter_file($table, { fsep => qr/\t/ });
	foreach my $stdfd (qw|Chrom Position Ref Alt|) {
		unless(grep { $_ eq $stdfd } @$fnames) {
			die "Cannot find standard field $stdfd from input!";
		}
	}
	my @newfields;
	foreach my $field (@$fnames) {
		unless(defined $known{$field}) {
			push @newfields => $field;
		}
		$known{$field} = 1;
	}
	if ($ii == 0) {
		$iter = $it;
		@outfields = @$fnames;
	}
	else {
		push @extrafields, @newfields;
		while(my $dat = $it->()) {
			my $varid = join(":", @{$dat}{qw|Chrom Position Ref Alt|});
			foreach my $field (@newfields) {
				if (defined $extrainfo{$varid} && defined $extrainfo{$varid}{$field}) {
					if ($extrainfo{$varid}{$field} ne $dat->{$field}) {
						warn "Field $field for variant $varid was already defined by a different value!";
					}
				}
				else {
					$extrainfo{$varid}{$field} = $dat->{$field};
				}
			}
		}
	}
}

open my $fout, ">$output" or die "Cannot write to $output";
print $fout join("\t", @outfields, @extrafields), "\n";
while (my $dat = $iter->()) {
	my $varid = join(":", @{$dat}{qw|Chrom Position Ref Alt|});
	my $info = $extrainfo{$varid};
	if (defined $info) {
		print $fout join("\t", @{$dat}{@outfields}, map { $_ // "." } @{$info}{@extrafields}), "\n";
	}
	else {
		print $fout join("\t", @{$dat}{@outfields}, ('.') x @extrafields), "\n";
	}
}



