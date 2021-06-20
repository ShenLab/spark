#!/usr/bin/env perl

use strict;
use warnings;
use Data::Dumper;
use Getopt::Euclid;
use List::MoreUtils qw|all any uniq notall none one|;
use Utils::File::Iter qw|iter_file|;
use Utils::Parser qw|sql_query|;
use Utils::List qw|parse_fields|;
use Utils::Hash qw|str2hash|;

my %opt;
if (defined $ARGV{'--fsep'}) {
	$opt{fsep} = qr/$ARGV{'--fsep'}/;
}
if (defined $ARGV{'--alias'}) {
	$opt{alias} = str2hash($ARGV{'--alias'}, { psep => ',', kvsep => ':' });
	# If all alias are numbers, will assume the file has no header
	if (all { /^\d+$/ } keys %{$opt{alias}}) {
		$opt{header} = 0;
	}
}

my ($it, $fnames) = iter_file($ARGV{'--input'} eq '-' ? \*STDIN : $ARGV{'--input'},
		{ skip => $ARGV{'--skip'}, %opt, 
		  sheet => $ARGV{'--sheet'}, maxcol => $ARGV{'--maxcol'}, maxrow => $ARGV{'--maxrow'} });

my @ofields;
unless ($ARGV{'--select'} || $ARGV{'--select-pattern'}) {
	@ofields = @$fnames;
}
else {
	my %known;
	if ($ARGV{'--select'}) {
		foreach my $field (parse_fields($ARGV{'--select'})) {
			unless (grep { $field eq $_ } @$fnames) {
				die "Cannot find field $field in the input file";
			}
			$known{$field} = 1;
		}
	}
	if ($ARGV{'--select-pattern'}) {
		my $pattern = qr/$ARGV{'--select-pattern'}/;
		@ofields = grep { /$pattern/ || defined $known{$_} } @$fnames;
	}
	else {
		@ofields = grep { defined $known{$_} } @$fnames;
	}
}


my (@efields, $callback);
if ($ARGV{'--filter'}) {
	my ($cb, $tokens) = sql_query($ARGV{'--filter'}, 1);
	foreach my $tok (@$tokens) {
		if ($tok->[0] eq 'FIELD') {
			unless(grep { $tok->[1] eq $_ } @$fnames) {
				die "The field $tok->[1] in filter expression cannot be found in the input";
			}
			push @efields, $tok->[1];
		}
	}
	$callback = $cb;
}
else {
	$callback = sub { 1 };
}


my (@mvfields, $mvcallback);
if ($ARGV{'--unpack'}) {
	@mvfields = split(',', $ARGV{'--unpack'});
	foreach my $field (@mvfields) {
		unless (grep { $field eq $_ } @efields) {
			die "Cannot find multi-value field $field in the expression";
		}
	}
	if ($ARGV{'--dupaction'} eq 'any') {
		$mvcallback = sub { any { $callback->($_) } @_ };
	}
	elsif ($ARGV{'--dupaction'} eq 'all') {
		$mvcallback = sub { all { $callback->($_) } @_ };
	}
	elsif ($ARGV{'--dupaction'} eq 'none') {
		$mvcallback = sub { none { $callback->($_) } @_ };
	}
	elsif ($ARGV{'--dupaction'} eq 'notall') {
		$mvcallback = sub { notall { $callback->($_) } @_ };
	}
	elsif ($ARGV{'--dupaction'} eq 'one') {
		$mvcallback = sub { one { $callback->($_) } @_ };
	}
	elsif ($ARGV{'--dupaction'} eq 'first') {
		$mvcallback = sub { $callback->($_[0]) };
	}
	else {
		die "Cannot recognize the dupaction $ARGV{'--dupaction'}";
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

print $fout join($ofsep, @ofields), "\n";

while(my $dat = $it->()) {
	if (@mvfields) {
		# Fix MV fields to avoid zero expansion
		foreach my $field (@mvfields) {
			if ($dat->{$field} eq "") {
				$dat->{$field} = ".";
			}
		}
		my @nelem = uniq map {  scalar(split($ARGV{'--vsep'}, $dat->{$_})) } @mvfields;
		unless(@nelem == 1) {
			print Dumper $dat;
			print join(",", @nelem), "\n";
			die "Not all multi-var fields has the same number of elements";
		}
		my @data;
		foreach my $field (@efields) {
			if (grep { $field eq $_ } @mvfields) {
				my @val = split($ARGV{'--vsep'}, $dat->{$field});
				foreach my $ii (0..($nelem[0]-1)) {
					$data[$ii]{$field} = $val[$ii];
				}
			}
			else {
				foreach my $ii (0..($nelem[0]-1)) {
					$data[$ii]{$field} = $dat->{$field};
				}
			}
		}
		if ($ARGV{'--negative'}) {
			unless ($mvcallback->(@data)) {
				print $fout join($ofsep, @{$dat}{@ofields}), "\n";
			}
		}
		else {
			if ($mvcallback->(@data)) {
				print $fout join($ofsep, @{$dat}{@ofields}), "\n";
			}
		}
	}
	else {
		if ($ARGV{'--negative'}) {
			unless ($callback->($dat)) {
				print $fout join($ofsep, @{$dat}{@ofields}), "\n";
			}
		}	
		else {
			if ($callback->($dat)) {
				print $fout join($ofsep, @{$dat}{@ofields}), "\n";
			}
		}
	}
}


__END__

=head1 NAME

flatdb.pl -- Query a flat file database.

=head1 USAGE

flatdb [options] -in INPUT -f FILTER -out OUTPUT 

=head1 NOTE

The input file should be database like data table stored in a text file (including csv format)
or in an Excel spreadsheet. The input file must contain a header with appropriate field 
names assigned to each column. The field name should match regex: qr/[A-Z_][A-Z0-9_\.]*/i

The output will be tab separated text file, unless csv suffix is given.

It is possible that multiple values are packed in the same value. We can unpack such fields
before applying the filtering expression. In such case, we need to specify the way we evaluate
data with multiple entries. Is it true when any, all, or first expression to be true?

=head1 REQUIRED ARGUMENTS

=over

=item -[-]in[put] [=] <file>

Input data file to be checked. Use '-' to indicate the input is from STDIN.
When input is from STDIN, header is required.

=for Euclid:
  file.type: readable

=back

=head1 OPTIONS

=over

=item -[-]f[ilter] [=] <expr>

Filtering expression.

=item -[-]neg[ative]

Negate the filtering expression.

=item -[-]out[put] [=] <file>

Output file name. If not provided, will dump to STDOUT.

=for Euclid:
  file.type: writable

=item -[-]alias [=] <fstr>

Fields to be renamed. When alias is used, field names in the filtering expression and
in the selected output should all be aliases.

=item -[-]select [=] <string>

Selected fields to appear in the output.

=item -[-]select-pattern [=] <regex>

Select fields that match the provided regex.

=item -[-]skip [=] <nrow>

Number of rows to skip at the begnining.

=item -[-]fsep [=] <regex>

For text file: Field separator, the default is comma for csv file
and white space for other types of text file.

=item -[-]ofsep [=] <string>

The output field separator, default is tab for text file or comma for csv file.

=item -[-]sheet [=] <num>

If the input is an Excel workbook, specify the sheet name or number.

=item -[-]maxrow [=] <num>

=item -[-]maxcol [=] <num>

For spreadsheet: the max range of columns and rows.
maxcol can also be used when input is from STDIN and does not have a header.

=item -[-]unpack [=] <fields>

Unpack the values in the given list of fields. When this list is given
it is required that number of unpacked values should be the same in all 
fields in this list.

=item -[-]vsep [=] <sepchar>

The separator for multiple values in the field, default is ",".

=for Euclid:
    sepchar.default: ','

=item -[-]dup[action] [=] <action>

How to evaluate with multiple values. Evaluate to be true when "any", "all", "none", "notall", "one",
or "first" value is true. Default: any.

=for Euclid:
    action.default: 'any'

=back 

