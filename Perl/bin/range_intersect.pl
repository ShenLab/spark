#!/usr/bin/env perl

use strict;
use warnings;
use Data::Dumper;
use Getopt::Euclid;
use List::MoreUtils qw|all any uniq|;
use Utils::List qw|parse_fields|;
use Utils::Hash qw|str2hash|;
use Utils::File::Iter qw|iter_file|;
use Genome::Ranges qw|validate_elem|;
use Genome::UCSC qw|hg_chrom is_hgchr|;
use Genome::UCSC::BinKeeper;

# Sanity check
if ($ARGV{'--afile'} eq $ARGV{'--bfile'}) {
	die "A and B file must be different";
}

# If field matching is required, store matched fields
my (@amatch, @bmatch);
if ($ARGV{'--match'}) {
	foreach my $fpair (split(',', $ARGV{'--match'})) {
		my @fds = split(':', $fpair);
		unless (@fds == 2) {
			die "Incorrect field pair: $fpair";
		}
		push @amatch, $fds[0];
		push @bmatch, $fds[1];
	}
	print STDERR "The following fields in A: ", join(",", @amatch), 
		" will be matched to the following fields in B: ", join(',', @bmatch), "\n";
}

# Extract database entry and store in BinKeeper
my $bk = Genome::UCSC::BinKeeper->new();
my ($bchr, @bfields);
{
	my ($it, $fnames, $keyfields) = prep_iter("b");
	if (defined $ARGV{'--with-b'}) {
		@bfields = map { $ARGV{'--with-b'}.$_ } @$fnames;
	}
	else {
		@bfields = @$fnames;
	}
	foreach my $bmat (@bmatch) {
		unless(grep { $bmat eq $_ } @$fnames) {
			die "Cannot find matching field $bmat in B file";
		}
	}
	while(my $dat = $it->()) {
		my ($chrom, $start, $end) = @{$dat}{@$keyfields};
		validate_elem($chrom, $start, $end);
		unless(defined $bchr) {
			if ($chrom =~ /^chr/) {
				$bchr = 1;
			}
			else {
				$bchr = 0;
			}
		}
		$start += 1 if $ARGV{'--b-bed'};
		if (exists $ARGV{'--with-b'} || defined $ARGV{'--match'}) {
			$bk->add($chrom, $start, $end, $dat);
		}
		else {
			$bk->add($chrom, $start, $end);
		}
	}
}

# Iterator over a file and find overlaps
my ($it, $fnames, $keyfields) = prep_iter("a");
my $achr;
my @afields = @$fnames;
foreach my $amat (@amatch) {
	unless(grep { $amat eq $_ } @afields) {
		die "Cannot find matching field $amat in A file";
	}
}

# Determine the fields to appear in the output file
my (@allfields, @outfields);
if (defined $ARGV{'--negative'}) {
	if (exists $ARGV{'--with-b'} || exists $ARGV{'--with-calc'} || defined $ARGV{'--with-loj'}) {
		warn "Under negative mode, B file fields or calculated fields will not appear in the output."
	}
	@allfields = @afields;
}
else {
	if (exists $ARGV{'--with-b'} || exists $ARGV{'--with-calc'} || defined $ARGV{'--with-loj'}) {
		@allfields = @afields;
		if (exists $ARGV{'--with-b'}) {
			foreach my $bfield (@bfields) {
				if (grep { $bfield eq $_ } @allfields) {
					die "Output field from B file $bfield was already used";
				}
				push @allfields, $bfield;
			}
			if (exists $ARGV{'--with-calc'}) {
				my @cfields;
				if (defined $ARGV{'--with-calc'}) {
					@cfields = map { $ARGV{'--with-calc'}.$_ } qw|iLen fA fB fAxfB|;
				}
				else {	
					@cfields = qw|iLen fA fB fAxfB|;
				}
				foreach my $cfield (@cfields) {
					if (grep { $cfield eq $_ } @allfields) {
						die "Calculated field $cfield was already used";
					}
					push @allfields, $cfield;
				}
			}
		}
		else {
			if (exists $ARGV{'--with-calc'}) {
				die "Option --with-calc must be used together with --with-b";
			}
			if (defined $ARGV{'--with-loj'}) {
				die "Option --with-loj must be used together with --with-b";
			}
			if (defined $ARGV{'--with-first'}) {
				die "Option --with-first must be used together with --with-b";
			}
		}
	}
	else {
		@allfields = @afields;
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

my $fout;
if ($ARGV{'--output'}) {
	open $fout, ">$ARGV{'--output'}";
}
else {
	$fout = \*STDOUT;
}
my %bkopt;
if (exists $ARGV{'--with-calc'} || defined $ARGV{'--with-first'}) {
	$bkopt{calc} = 1;
}

unless(all { $_ =~ /^\d+$/ } @outfields) {
	print $fout join("\t", @outfields), "\n";
}
while(my $dat = $it->()) {
	my ($chrom, $start, $end) = @{$dat}{@$keyfields};
	validate_elem($chrom, $start, $end);
	unless(defined $achr) {
		if ($chrom =~ /^chr/) {
			$achr = 1;
		}
		else {
			$achr = 0;
		}
		unless($achr == $bchr) {
			warn "Chromosome nomenclatures in A and B files are different";
		}
	}

	if ($achr == 0 && $bchr == 1) {
		$chrom = hg_chrom($chrom);
	}
	elsif ($achr == 1 && $bchr == 0) {
		$chrom =~ s/^chr//; $chrom = 'MT' if $chrom eq 'M';
	}
	$start += 1 if $ARGV{'--a-bed'};

	
	my @overlaps;
	if ($ARGV{'--match'}) {
		my $matchval = join("\t", map { $dat->{$_} } @amatch);
		@overlaps = grep { join("\t", @{$_->[2]}{@bmatch}) eq $matchval }
					 $bk->find_range($chrom, $start, $end, 
									{ qCover => $ARGV{'--a-overlap'}, 
									  tCover => $ARGV{'--b-overlap'}, %bkopt });
	}
	else {
		@overlaps = $bk->find_range($chrom, $start, $end, 
									{ qCover => $ARGV{'--a-overlap'}, 
									  tCover => $ARGV{'--b-overlap'}, %bkopt });
	}

	unless(@overlaps) {
		if ($ARGV{'--negative'}) {
			print $fout join("\t", @{$dat}{@outfields}), "\n";
		}
		elsif ($ARGV{'--with-loj'}) {
			print $fout join("\t", map { $_ // "." } @{$dat}{@outfields}), "\n";
		}
		next;
	}
	else {
		next if $ARGV{'--negative'};
	}

	# Sort based on TxQ
	if ($bkopt{calc}) {
		@overlaps = sort { $b->[3]{TxQ} <=> $a->[3]{TxQ} } @overlaps;
	}

	if (exists $ARGV{'--with-b'}) {
		my $bprefix = $ARGV{'--with-b'} // "";
		foreach my $overlap (@overlaps) {
			my $info = $overlap->[2];
			foreach my $key (keys %$info) {
				$dat->{$bprefix.$key} = $info->{$key};
			}
			if (exists $ARGV{'--with-calc'}) {
				my $cprefix = $ARGV{'--with-calc'} // "";
				my $stat = $overlap->[3];
				$dat->{$cprefix."iLen"} = $stat->{iLen};
				$dat->{$cprefix."fA"} = sprintf("%.$ARGV{'--ndigit'}f", $stat->{Q});
				$dat->{$cprefix."fB"} = sprintf("%.$ARGV{'--ndigit'}f", $stat->{T});
				$dat->{$cprefix."fAxfB"} =  sprintf("%.$ARGV{'--ndigit'}f", $stat->{TxQ});
			}
			print $fout join("\t", @{$dat}{@outfields}), "\n";
			last if $ARGV{'--with-first'};
		}
	}
	else {
		print $fout join("\t", @{$dat}{@outfields}), "\n";
	}
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
	if (defined $fields) {
		if (all { /^\d+$/ } keys %$fields) {
			$opt{header} = 0;
		}
		else {
			$opt{header} = 1;
		}
	}
	my @fields;
	if (defined $ARGV{"--$prefix-fields"}) {
		@fields = parse_fields($ARGV{"--$prefix-fields"});
		unless(defined $opt{header}) {
			if (all { /^\d+$/ } @fields) {
				$opt{header} = 0;
			}
			else {
				$opt{header} = 1;
			}
		}
	}
	else {
		unless($ARGV{"--${prefix}file"} =~ /\.bed$/ || $ARGV{"--${prefix}file"} =~ /\.bed\.gz$/) {
			die "Must provide field names for Chrom,Start,End for non-BED file: " . $ARGV{"--${prefix}file"};
		}
		if (defined $opt{alias}) {
			@fields = map { $opt{alias}{$_} // $_ } qw|1 2 3|;
		}
		elsif (defined $opt{select}) {
			@fields = map { $opt{select}{$_} // $_ } qw|1 2 3|;
		}
		else {
			$opt{header} = 0;
			@fields = qw|1 2 3|;
		}
	}
	my ($it, $fnames) = iter_file($ARGV{"--${prefix}file"} eq '-' ? \*STDIN : $ARGV{"--${prefix}file"},
		{ skip => $ARGV{"--$prefix-skip"}, sheet => $ARGV{"--$prefix-sheet"}, 
		  maxcol => $ARGV{"--$prefix-maxcol"}, maxrow => $ARGV{"--$prefix-maxrow"}, %opt });
	if ($ARGV{"--${prefix}file"} =~ /\.bed$/ || $ARGV{"--${prefix}file"} =~ /\.bed\.gz$/) {
		print STDERR "$prefix file is in BED format by suffix\n";
		$ARGV{'--${prefix}-bed'} = 1;
	}
	# Check field names first
	foreach my $field (@fields) {
		unless(grep { $field eq $_ } @$fnames) {
			die "Cannot find field $field in $prefix-file: ".$ARGV{"--${prefix}file"};
		}
	}
	return ($it, $fnames, \@fields);
}


__END__

=head1 NAME

range_intersect.pl -- Finding overlap genomic features.

=head1 NOTES

This is the inhouse implementation of BEDtools' intersect utility with that allows more flexible 
input file formats. 

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

Also output fields in file B. Optionally adding prefix to column names from file B. 

=item -[-]with-loj

Output all entries in A regardless of their overlap with B, "Left outer join".
Can only be used together with --with-b.

=item -[-]with-calc [=] [<prefix>]

Also output calculated fields of overlapping length and fractions.
Currently, there are four calculated fields: iLen, fA, fB, fAxfB. An optional prefix can be added to
these calculated fields. It must be used together with --with-b. But Both --with-b and --with-calc 
will be inactivated if --negative is switchd on.

=item -[-]with-first

Only output first matched entry in file B. Must be used together with --with-b.
The entry in B that have the largest overlap (pT*pQ) will be shown in the output.

=item -[-]ndigit [=] <num>

Number of digits for calculated numerical values.

=for Euclid
	num.default: 3

=item -[-]match [=] <fpairs>

Matching values of specified field pairs in A and B files. Field pairs should be comma separated,
fields in A and B should be joined by colon. For example: IID:IID,Type:CNType will match values of 
IID and Type from A file with values of IID and CNType in B file. Field names used for matching
are renamed field names if files fields have been renamed.

=item -[-]neg[ative]

Instead of extracting overlapping regions, reverse the operation to output non-overlapping regions.
with-b, with-loj and with-calc will be disabled under negative mode. If fields matching is enabled,
negative mode will also output overlapping regions that not matching specified fields.

=item -[-]select [=] <string>

Selected fields to appear in the output.

=item -[-]select-pattern [=] <regex>

Select fields that match the provided regex.

=item -[-]a[-]fields [=] <fstr>

Specify key fields (Chrom,Start,End) in A file. When alias is used, they must be renamed field names.
If files have .bed suffix and fields are not specified, then we will assume they are bed format.

=item -[-]a-select [=] <fstr>

Select and rename fields in A file.

=item -[-]a-alias [=] <fstr>

Rename selected fields in A file. Example: 1:Chrom,2:Start,3:End.
If original field names are all digit, file will be assumed to have no header.

=item -[-]a-bed

Indicate that file A has BED style genomic intervals.

=item -[-]a-skip [=] <nrow>

=item -[-]a-fsep [=] <regex>

=item -[-]a-sheet [=] <number>

=item -[-]a-maxrow [=] <number>

=item -[-]a-maxcol [=] <number>

Options for paring A file.

=item -[-]b[-]fields [=] <fstr>

Specify key fields in B file.

=item -[-]b-select [=] <fstr>

Select and rename fiels in B file.

=item -[-]b-alias [=] <fstr>

Rename selected fields in B file.

=item -[-]b-bed

Indicate that file B has BED style genomic intervals.

=item -[-]b-skip [=] <nrow>

=item -[-]b-fsep [=] <regex>

=item -[-]b-sheet [=] <number>

=item -[-]b-maxrow [=] <number>

=item -[-]b-maxcol [=] <number>

Options for parsing B file.

=item -[-]a-over[lap] [=] <perc>

Minimial overlap requried as fraction of A.

=for Euclid:
	perc.default: 1e-10

=item -[-]b-over[lap] [=] <perc>

Minimal overlap required as fraction of B.

=for Euclid:
	perc.default: 1e-10

=back



