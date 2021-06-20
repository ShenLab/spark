use strict;
use warnings;
use Getopt::Euclid;
use Utils::Hash qw|chk_default|;
use Utils::File::Iter qw|iter_file|;

my (%remove, %update);

unless (defined $ARGV{'--remove'} || defined $ARGV{'--update'}) {
	die "Must provide at least one remove or update table";
}

my @fields = split(',', $ARGV{'--fields'});

if ($ARGV{'--remove'}) {
	my ($it, $fnames) = iter_file($ARGV{'--remove'}, { fsep => qr/\t/, sheet => $ARGV{'--remove-sheet'} });
	while(my $dat = $it->()) {
		my ($iid, $chrom, $pos, $ref, $alt) = @{$dat}{@fields};
		my $varid = join(":", $chrom, $pos, $ref, $alt);
		$remove{$iid,$varid} = 1;
	}
}

my @newcols;
if ($ARGV{'--update'}) {
	my ($it, $fnames) = iter_file($ARGV{'--update'}, { fsep => qr/\t/, sheet => $ARGV{'--update-sheet'} });
	@newcols = @$fnames;
	while(my $dat = $it->()) {
		my ($iid, $chrom, $pos, $ref, $alt) = @{$dat}{@fields};
		my $varid = join(":", $chrom, $pos, $ref, $alt);
		$update{$iid,$varid} = $dat;
	}
}

my (%lookup, $f_key);
if (defined $ARGV{'--lookup'}) {
	my @cols = split(',', $ARGV{'--lookup'});
	my ($it, $fnames) = iter_file($ARGV{'--input'} eq '-' ? \*STDIN : $ARGV{'--input'}, 
		{ fsep => qr/\t/, sheet => $ARGV{'--input-sheet'} });
	foreach my $col (@cols) {
		unless(grep { $_ eq $col } @$fnames) {
			die "Cannot find column $col in input table";
		}
		$lookup{$col} = {};
	}
	$f_key = shift @cols;
	while(my $dat = $it->()) {
		foreach my $f_val (@cols) {
			chk_default($lookup{$f_val}, $dat->{$f_key}, $dat->{$f_val});
		}
	}
}

my ($it, $fnames) = iter_file($ARGV{'--input'} eq '-' ? \*STDIN : $ARGV{'--input'},
		{ fsep => qr/\t/, sheet => $ARGV{'--input-sheet'} });

my @addcols;
if ($ARGV{'--add-col'}) {
	foreach my $col (@newcols) {
		unless (grep { $_ eq $col } @$fnames) {
			push @addcols, $col;
		}
	}
}

my @outfields = @$fnames;
push @outfields, @addcols;

my $fout;
if ($ARGV{'--output'}) {
	open $fout, ">$ARGV{'--output'}" or die "Cannot write to $ARGV{'--output'}";
}
else {
	$fout = \*STDOUT;
}

print $fout join("\t", @outfields), "\n";
while(my $dat = $it->()) {
	my ($iid, $chrom, $pos, $ref, $alt) = @{$dat}{@fields};
	my $varid = join(":", $chrom, $pos, $ref, $alt);
	if ($remove{$iid,$varid}) {
		delete $remove{$iid,$varid};
		next;
	}
	elsif ($update{$iid,$varid}) {
		foreach my $field (@outfields) {

			$dat->{$field} = $update{$iid,$varid}{$field};
		}
		delete $update{$iid,$varid};
	}
	print $fout join("\t", @{$dat}{@outfields}), "\n";
}

foreach my $viid (sort keys %update) {
	my $vdat = $update{$viid};
	foreach my $field (keys %lookup) {
		if (defined $lookup{$field}{$vdat->{$f_key}}) {
			$vdat->{$field} = $lookup{$field}{$vdat->{$f_key}};
		} 
	}
	print $fout join("\t", map { $_ // "." } @{$vdat}{@outfields}), "\n";
}

if (keys %remove) {
	print STDERR "The following variants in the removal list cannot be found in input\n";
	foreach my $viid (sort keys %remove) {
		my ($iid, $vid) = split($;, $viid);
		print STDERR $iid, "\t", $vid, "\n";
	}
}

__END__

=head1 NAME

update_vars.pl -- Update variants table

=head1 REQUIRED ARGUMENTS
 
=over
 
=item -[-]in[put] [=] <table>

The input variant table.

=back
 
=head1 OPTIONS
 
=over

=item -[-]out[put] [=] <prefix>

The updated output variant table.

=item -[-]fields [=] <string>

Default fields for sample ID and variant. Should be the same across all input files.

=for Euclid:
	string.default: 'IID,Chrom,Position,Ref,Alt'

=item -[-]input-sheet [=] <num>

If input list is in excel file, provude sheet number or name.

=item -[-]remove [=] <table>

A table of variants that should be removed from original input.
If the input is text file, it must be tab separated.

=item -[-]remove-sheet [=] <num>

If remove list is in excel file, provide sheet number or name.

=item -[-]update [=] <table>

A table of new variants or variants with new annotations.

=item -[-]update-sheet [=] <num>

If update list in excel file, provide sheet number of name.

=item -[-]lookup [=] <string>

For columns listed by this option, fetch values from original table by looking up data from original table.
e.g.: "IID,Father,Mother", will lookup parental IDs that are known to sample's IID.

=item -[-]add-col

Adding new columns to the original table. The default is to ignore.
All newly added columns will appear at the end of existing columns.

=back