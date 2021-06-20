use strict;
use warnings;
use File::Copy qw|move|;
use Genome::UCSC qw|hg_chrom|;
use Utils::File qw|swap_files|;
use Utils::File::Iter qw|iter_file|;


unless(@ARGV == 2 || @ARGV == 3) {
	print STDERR "Fix chromsome code in variant table\nUsage:\n";
	print STDERR "$0 Input ChrCol [Output]\n";
	exit 1;
}

my ($input, $f_chr, $output) = @ARGV;

my $fout;
if (defined $output) {
	if ($input eq $output) {
		die "If output file is specified, it cannot be the same as input!";
	}
	open $fout, ">$output" or die "Cannot write to $output";
}
else {
	$output = "$input.bak";
	open $fout, ">$output" or die "Cannot write to $output";
}

my %opt;
if ($f_chr =~ /^\d+$/) {
	$opt{header} = 0;
}
else {
	$opt{header} = 1;
}
if ($input =~ /\.txt$/ || $input =~ /\.tsv/) {
	$opt{fsep} = qr/\t/;
}
elsif ($input =~ /\.xlsx$/) {
	$opt{sheet} = 1;
}

my ($it, $fnames) = iter_file($input, \%opt);
unless(grep { $f_chr eq $_ } @$fnames) {
	die "Cannot find chromosome column $f_chr from input!";
}

if ($opt{header}) {
	print $fout join("\t", @$fnames), "\n";
}

my $flag;
while(my $dat = $it->()) {
	my $chrom = $dat->{$f_chr};
	unless(defined $flag) {
		$flag = $chrom =~ /^chr/ ? 1 : 0;
	}
	if ($flag) {
		$chrom =~ s/^chr//; $chrom = 'MT' if $chrom eq 'M';
	}
	else {
		$chrom = hg_chrom($chrom);
	}
	$dat->{$f_chr} = $chrom;
	print $fout join("\t", @{$dat}{@$fnames}), "\n";
}
close $fout;

if ($output eq "$input.bak") {
	swap_files($input, "$input.bak");
}

