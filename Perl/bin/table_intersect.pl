#!/usr/bin/env perl

use strict;
use warnings;
use FindBin qw|$Bin|;
use Data::Dumper;
use Getopt::Euclid;
use List::MoreUtils qw|all any uniq zip6|;
use Utils::Hash qw|str2hash|;
use Utils::List qw|parse_fields|;
use Utils::File::Iter qw|iter_file|;



# Parse file headers to determine command fields.
if ($ARGV{'--afile'} eq $ARGV{'--bfile'}) {
	die "A and B file must be different";
}
my ($a_it, $a_fnames) = prep_iter("a");
my ($b_it, $b_fnames) = prep_iter("b");

# Find fields for matching ("_keys", with remaining as "_vals")
# If matching fields is not provided, we will look for share field names
my (@a_keys, @a_vals, @b_keys, @b_vals);
if ($ARGV{'--match'}) {
	foreach my $fpair (split(',', $ARGV{'--match'})) {
		my @fds = split(':', $fpair);
		unless(grep { $fds[0] eq $_ } @$a_fnames) {
			die "Cannot find matching field $fds[0] in file A";
		}
		unless(grep { $fds[1] eq $_ } @$b_fnames) {
			die "Cannot find matching field $fds[1] in file B";
		}
		push @a_keys, $fds[0];
		push @b_keys, $fds[1];
	}
	foreach my $field (@$a_fnames) {
		unless(grep { $field eq $_ } @a_keys) {
			push @a_vals, $field;
		}
	}
	foreach my $field (@$b_fnames) {
		unless(grep { $field eq $_ } @b_keys) {
			push @b_vals, $field;
		}
	}
}
else {
	foreach my $field (@$a_fnames) {
		if (grep { $_ eq $field } @{$b_fnames}) {
			push @a_keys, $field;
			push @b_keys, $field;
		}
		else {
			push @a_vals, $field;
		}
	}
	foreach my $field (@$b_fnames) {
		unless(grep { $field eq $_ } @b_keys) {
			push @b_vals, $field;
		}
	}
}

unless(@a_keys > 0 && @b_keys > 0) {
	die "Cannot find key fields share by A and B files";
}
else {
	print STDERR "Two files will be joined with following fields: ", 
		join(',', map { "$_->[0]:$_->[1]" } zip6 @a_keys, @b_keys), "\n";
}

# Slurp B file data to xref
my %xref;
while(my $dat = $b_it->()) {
	my $key = join("\t", @{$dat}{@b_keys});
	if (defined $xref{$key}) {
		warn "Data for [$key] in B file has already been found!";
		next;
	}
	if (@b_vals > 0 || exists $ARGV{'--with-b'}) {
		$xref{$key} = $dat;
	}
	else {
		$xref{$key} = 1;
	}
}

# Determine output fields
my (@allfields, @outfields);
if (defined $ARGV{'--negative'}) {
	if (exists $ARGV{'--with-b'} || defined $ARGV{'--with-loj'}) {
		warn "Under negative mode, B file fields will not appear in the output!";
	}
	@allfields = @$a_fnames;
}
else {
	@allfields = @$a_fnames;
	if (@b_vals > 0 || exists $ARGV{'--with-b'}) {
		unless(@b_vals > 0) {
			die "B file has no value field for output --with-b!"
		}
		foreach my $bfield (map { defined $ARGV{'--with-b'} ? $ARGV{'--with-b'}.$_ : $_ } @b_vals) {
			push @allfields, $bfield;
		}
	}
	else {
		if (defined $ARGV{'--with-loj'}) {
			die "Option --with-loj must be used together with --with-b!";
		}
	}
}

unless ($ARGV{'--select'} || $ARGV{'--select-pattern'}) {
	@outfields = @allfields;
}
else {
	my %known;
	if ($ARGV{'--select'}) {
		foreach my $field (parse_fields($ARGV{'--select'})) {
			unless (grep { $field eq $_ } @allfields) {
				die "Cannot find field $field in output file";
			}
			$known{$field} = 1;
		}
	}
	if ($ARGV{'--select-pattern'}) {
		my $pattern = qr/$ARGV{'--select-pattern'}/;
		@outfields = grep { /$pattern/ || defined $known{$_} } @allfields;
	}
	else {
		@outfields = grep { defined $known{$_} } @allfields;
	}
}

my ($fout, $ofsep);
if (defined $ARGV{'--ofsep'}) {
	$ofsep = $ARGV{'--ofsep'};
}
else {
	if ($ARGV{'--output'} && $ARGV{'--output'} =~ /\.csv$/) {
		$ofsep = ",";
	}
	else {
		$ofsep = "\t";
	}
}
if ($ARGV{'--output'}) {
	open $fout, ">$ARGV{'--output'}";
}
else {
	$fout = \*STDOUT;
}

unless(all { $_ =~ /^\d+$/ } @outfields) {
	print $fout join("\t", @outfields), "\n";
}

while(my $dat = $a_it->()) {
	my $key = join("\t", @{$dat}{@a_keys});

	unless(defined $xref{$key}) {
		if ($ARGV{'--negative'}) {
			print $fout join("\t", @{$dat}{@outfields}), "\n";
		}
		elsif ($ARGV{'--with-loj'}) {
			print $fout join("\t", map { $_ // $ARGV{'--na-str'} } @{$dat}{@outfields}), "\n";
		}
		next;
	}
	else {
		next if $ARGV{'--negative'};
	}

	if (@b_vals > 0 || exists $ARGV{'--with-b'}) {
		my $bprefix = $ARGV{'--with-b'} // "";
		foreach my $bfield (@b_vals) {
			$dat->{$bprefix.$bfield} = $xref{$key}{$bfield};
		}
	}
	print $fout join("\t", @{$dat}{@outfields}), "\n";
}


sub prep_iter {
	my ($prefix) = @_;
	my %opt;
	if (defined $ARGV{"--$prefix-fsep"}) {
		$opt{fsep} = qr/$ARGV{"--$prefix-fsep"}/;
	}
	my $fields;
	if (defined $ARGV{"--$prefix-select"}) {
		if (defined $ARGV{"--$prefix-alias"}) {
			warn "--$prefix-alias cannot be used together with --$prefix-select!";
		}
		$fields = str2hash($ARGV{"--$prefix-select"}, { psep => ',', kvsep => ':' });
		$opt{select} = $fields;
	}
	elsif (defined $ARGV{"--$prefix-alias"}) {
		$fields = str2hash($ARGV{"--$prefix-alias"}, { psep => ',', kvsep => ':' });
		$opt{alias} = $fields;
	}
	# If force noheader is provided
	if (defined $ARGV{"--$prefix-noheader"}) {
		$opt{header} = 0;
	}
	else {
		if (defined $fields) {
			if (all { /^\d+$/ } keys %$fields) {
				$opt{header} = 0;
			}
			else {
				$opt{header} = 1;
			}
		}
	}
	if ($ARGV{'--match'}) {
		my @fpairs = map { [split(':', $_)] } split(',', $ARGV{'--match'});
		unless(all { @$_ == 2 } @fpairs) {
			die "Incorrect field pairs: $ARGV{'--match'}";
		}
		my %suf = (a => 0, b => 1);
		unless(defined $opt{header}) {
			if (all { $_->[$suf{$prefix}] =~ /^\d+$/ } @fpairs) {
				$opt{header} = 0;
			}
			else {
				$opt{header} = 1;
			}
		}
	}

	my ($it, $fnames) = iter_file($ARGV{"--${prefix}file"} eq '-' ? \*STDIN : $ARGV{"--${prefix}file"},
		{ skip => $ARGV{"--$prefix-skip"}, sheet => $ARGV{"--$prefix-sheet"}, 
		  maxcol => $ARGV{"--$prefix-maxcol"}, maxrow => $ARGV{"--$prefix-maxrow"}, %opt });

	return ($it, $fnames);
}

__END__

=head1 NAME

table_intersect.pl -- Find overlapping entris in two tables.

=head1 NOTES

This utility is similar to range_intersect.pl. 

The purpose is to select entries from A files that overlap with those in B file based
on the common fields they share. Common fields are determined by the same field names.
There are options to rename fields from both files.

=head1 REQUIRED ARGUMENTS

=over

=item -[-]a[file] [=] <file>

File "A", use "-" if passing A from stdin.

=for Euclid:
	file.type: readable

=item -[-]b[file] [=] <file>

File "B", use "-" if passing B from stdin. But only one of A or B can be read from stdin.

=for Euclid:
	file.type: readable

=back

=head1 OPTIONS

=over

=item -[-]out[put] [=] <file>

The output file name. If not provided, will write to STDOUT.

=item -[-]with-b [=] [<prefix>]

Also output fields in B that are not used for matching. 
Optionally adding prefix to column names from B to avoid name conflict.

=item -[-]with-loj

Output all entries in A regardless of their overlap with B. "Left outer join".
It can only be used together with --with-b.

=item -[-]na-str [=] <string>

Under loj mode, string for missing data.

=for Euclid:
	string.default: "."

=item -[-]match [=] <fpairs>

Matching values of specified field pairs in A and B files. Field pairs should be comma 
separated, field in A and B should be joined by colon. For example: IID:SampID,FID:FamID
will match values of IID and IID from A file with values of SampID,FamID in B file.
The field names used for matching are renamed field names if file fields are renamed.

=item -[-]ofsep [=] <string>

The output field separator, default is tab for text file or comma for csv file.

=item -[-]neg[ative]

Reverse the operation to output A file entries that cannot be matched in B file.

=item -[-]select [=] <string>

Selected fields to appear in the output.

=item -[-]select-pattern [=] <regex>

Select fields that match the provided regex.

=item -[-]a[-]select [=] <fstr>

Select fields to be retained for A file, optionally with alias.
If this option is turned on, only specified fields will appear in the output.
If a file has no header, then either "*-select" or "*-alias" must be specified.

=item -[-]a-alias [=] <fstr>

Fields to be renamed for A file. All fields will appear in the output.
Only use this option to rename fields for matching with B file.

=item -[-]a-skip [=] <nrow>

=item -[-]a-fsep [=] <regex>

=item -[-]a-sheet [=] <number>

=item -[-]a-maxrow [=] <number>

=item -[-]a-maxcol [=] <number>

=item -[-]a-noheader

Options for paring A file.

=item -[-]b[-]select [=] <fstr>

Select fields to be retained for B file.

=item -[-]b-alias [=] <fstr>

Fields to to renamed for B file. 

=item -[-]b-skip [=] <nrow>

=item -[-]b-fsep [=] <regex>

=item -[-]b-sheet [=] <number>

=item -[-]b-maxrow [=] <number>

=item -[-]b-maxcol [=] <number>

=item -[-]b-noheader

Options for paring B file.

=back




