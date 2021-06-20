use strict;
use warnings;
use Data::Dumper;
use FindBin qw|$Bin|;
use Config::Std;
use Genome::UCSC::BinKeeper;

use lib "$Bin/../../lib";
use Shared qw|parse_tabfile|;
use Variants qw|get_vcfped_samps|;

# Subset the input variant table to keep sample appearing in VCF
# and variants overlapping the given regions
my ($conf, $vcffile, $region);
if (@ARGV == 3) {
	($conf, $vcffile, $region) = @ARGV;
}
elsif (@ARGV == 2) {
	($conf, $vcffile) = @ARGV;
}
else {
	print STDERR "$0 CONF VCFFILE [REGION]\n";
}


# Read configs to find how to parse PED file and INPUT
read_config $conf => my %config;

my %arg;
if ($config{PATH}{PROG} =~ /dnv_table/) {
	$arg{trios} = 1;
}

my $samps = get_vcfped_samps($vcffile, $config{Pedigree}{File}, { ignore => $config{Pedigree}{Ignore}, %arg });


my $callback;
if (!defined $region) {
	$callback = sub { 1 }
}
elsif (-f $region) {
	my $bk = Genome::UCSC::BinKeeper->new($region);
	$callback = sub {
		my ($chr, $pos, $ref, $alt) = @_;
		my $start = $pos - 1;
		my $end = $pos + length($ref) + 1;
		return $bk->any_overlap($chr, $start, $end);
	}
}
elsif ($region =~ /^(\w+):([\d,]+)\-([\d,]+)$/ || $region =~ /^(\w+)$/) {
	my ($chrom, $start, $end) = ($1, $2, $3);
	if (defined $start && defined $end) {	
		$start =~ s/,//g; $end =~ s/,//g;
	}
	$callback = sub {
		my ($chr, $pos, $ref, $alt) = @_;
		my $st = $pos - 1;
		my $ed = $pos + length($ref) + 1;
		if ($chr eq $chrom) {
			unless(defined $start && defined $end) {
				return 1;
			}
			else {
				if (range_overlap($start, $end, $st, $ed)) {
					return 1;
				}
				else {
					return 0;
				}
			}
		}
		else {
			return 0;
		}
	}
}
else {
	die "Cannot find region file or recognize region spec: $region";
}

my ($it, $fnames, $fields) = parse_tabfile($config{Input}{File}, $config{Input}{Fields}, 5, 5);
print join("\t", @$fnames), "\n";
while(my $dat = $it->()) {
	my ($iid, $chr, $pos, $ref, $alt) = @{$dat}{@$fields};
	if (defined $samps->{$iid} && $callback->($chr, $pos, $ref, $alt)) {
		print join("\t", @{$dat}{@$fnames}), "\n";
	}	
}

