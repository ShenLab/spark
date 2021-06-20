use strict;
use warnings;
use File::Path qw|make_path|;
use File::Temp qw|tempdir|;
use List::Util qw|first|;
use Perl6::Slurp;
use List::MoreUtils qw|none any all|;
use Data::Dumper;
use Getopt::Euclid;
use Utils::List qw|all_combs|;
use Utils::Hash qw|peek_hash|;
use Utils::Stat qw|mean|;
use Utils::File qw|count_line|;
use Utils::File::Iter qw|iter_file|;
use FaSlice;
use Genome::UCSC qw|hg_chrom|;
use Genome::UCSC::TwoBit;


my ($it, $fnames) = iter_file($ARGV{'--input'});
my @addfields = split(',', $ARGV{'--field-add'});
foreach my $field (@addfields) {
	if (grep { $field eq $_ } @$fnames) {
		die "Field to be added $field already exists in the input";
	}
}

my @required;
if (grep { /PL\.(HOMREF|HET|HOMALT)/ } @$fnames) {
	@required = qw|DP GQ PL.HOMREF PL.HET PL.HOMALT|;
}
else {
	@required = qw|DP GQ|;
}
my @reqbysamp;
foreach my $samp (qw|Offspring Father Mother|) {
	push @reqbysamp, map { $_.".$samp" } @required;
}
foreach my $field (qw|IID Gender Chrom Position Ref Alt|, @reqbysamp) {
	unless(grep { $field eq $_ } @$fnames) {
		die "Cannot find required field $field in the input file";
	}
}

my $wrkdir;
if ($ARGV{'--wrkdir'}) {
	$wrkdir = $ARGV{'--wrkdir'};
	make_path $wrkdir unless -d $wrkdir;
}
else {
	$wrkdir = tempdir(CLEANUP => 1);
}

if ($ARGV{'--mutrate'}) {
	unless($ARGV{'--seq'}) {
		die "Must also provide genome reference sequence to look for sequence context!";
	}
}

open my $fvcf, ">$wrkdir/input.vcf" or die "Cannot write to $wrkdir/input.vcf";
open my $fvcf2, ">$wrkdir/input_maleX.vcf" or die "Cannot write to $wrkdir/input_maleX.vcf";
print $fvcf join("\t", qw|#CHROM POS ID REF ALT QUAL FILTER INFO FORMAT Offspring Father Mother|), "\n";
print $fvcf2 join("\t", qw|#CHROM POS ID REF ALT QUAL FILTER INFO FORMAT Offspring Father Mother|), "\n";

my $inchr;
while(my $dat = $it->()) {
	unless($inchr) {
		if ($dat->{Chrom} =~ /^chr/) {
			$inchr = 1;
		}
		else {
			$inchr = 0;
		}
	}

	my %var;
	foreach my $field (qw|Chrom Position Ref Alt|) {
		$var{$field} = $dat->{$field};
	}
	$var{ID} = join(":", @{$dat}{qw|IID Chrom Position Ref Alt|});
	$var{Format} = "GT:DP:PL";

	# Deal with dup/twin case, when multiple data appear in the same cell
	# We will take the first one
	foreach my $samp (qw|Offspring Father Mother|) {
		if (any { $dat->{$_.".$samp"} =~ /,/ } @required) {
			unless(all { $dat->{$_.".$samp"} =~ /,/ } @required) {
				foreach my $fd (sort @required) {
					print STDERR $fd, "\t", $dat->{$fd}, "\n";
				}
				die "For twin/dup pairs, all required fields should have multiple values";
			}
			foreach my $field (@required) {
				$dat->{$field.".$samp"} = (split(',', $dat->{$field.".$samp"}))[0];
			}
		}
	}

	# Fix GLnexus's missing value in PL for homo ref genotype
	foreach my $samp (qw|Offspring Father Mother|) {
		if ((any { !defined $dat->{"PL.$_.$samp"} } qw|HOMREF HET HOMALT|) ||
			(all { $dat->{"PL.$_.$samp"} eq '.' } qw|HOMREF HET HOMALT|)) {
			if ($samp eq 'Offspring') {
				$dat->{"PL.HET.$samp"} = 0;
				$dat->{"PL.HOMREF.$samp"} = $dat->{"GQ.$samp"};
				$dat->{"PL.HOMALT.$samp"} = 10*$dat->{"GQ.$samp"};
			}
			else {
				$dat->{"PL.HOMREF.$samp"} = 0;
				$dat->{"PL.HET.$samp"} = $dat->{"GQ.$samp"};
				$dat->{"PL.HOMALT.$samp"} = 10*$dat->{"GQ.$samp"};
			}
		}
	}

	my $flag;
	foreach my $samp (qw|Offspring Father Mother|) {
		if (none { $dat->{"PL.$_.$samp"} == 0 } qw|HOMREF HET HOMALT|) {
			warn "None of ${samp}'s PL reached max 0: $var{ID}";
			$flag = 1;
		}
	}
	# Ignore problematic variants
	next if $flag;

	my %gcode = (HOMREF => '0/0', HET => '0/1', HOMALT => '1/1');
	my %gt;
	foreach my $samp (qw|Offspring Father Mother|) {
		my $geno = first { $dat->{"PL.$_.$samp"} == 0 } qw|HOMREF HET HOMALT|;
		if ($samp =~ /Father|Mother/) {
			unless($geno eq 'HOMREF') {
				die "Parents do not have HOMREF genotype!";
			}
		}
		else {
			if ($geno eq 'HOMREF') {
				die "Offspring have HOMREF genotype!";
			}
			elsif ($geno eq 'HOMALT') {
				# HOMALT can occasionally happen when there is a deletion, we will convert it to HET
				if ($dat->{Chrom} eq 'X' && $dat->{Gender} eq 'Female' || $dat->{Chrom} ne 'X') {
					warn "HOMALT DNV on female X or autosome";
					$geno = "HET";
					($dat->{"PL.HOMALT.Offspring"}, $dat->{"PL.HET.Offspring"}) = 
					($dat->{"PL.HET.Offspring"}, $dat->{"PL.HOMALT.Offspring"});
				} 
			}
		}
		$gt{$samp} = $geno;
		$var{$samp} = "$gcode{$geno}:".$dat->{"DP.$samp"}.":".
			join(",", map { $dat->{"$_.$samp"} } qw|PL.HOMREF PL.HET PL.HOMALT|);
	}
	my $line = join("\t", @var{qw|Chrom Position ID Ref Alt|}, (".") x 3,  @var{qw|Format Offspring Father Mother|});
	if ($gt{Offspring} eq 'HOMALT') {
		# Male X het will be treated as autosomal variants, but should be ideally removed from evaluation
		if ($dat->{Chrom} eq 'X' && $dat->{Gender} eq 'Male') {
			print $fvcf2 $line, "\n";
		}
		else {
			if ($dat->{Chrom} eq 'X' && $dat->{Gender} eq 'Female') {
				warn "HOMALT DNV on female chrX: $var{ID}";
			}
			else {
				warn "HOMALT DNV on non-chrX: $var{ID}"
			}
		}
	}
	else {
		print $fvcf $line, "\n";
	}
}
close $fvcf;
close $fvcf2;

my %dq;
open my $fped, ">$wrkdir/input.ped" or die "Cannot write to $wrkdir/input.ped";
print $fped <<EOF;
Trio Offspring Father Mother 2
Trio Father	0	0	1
Trio Mother	0	0	2
EOF
close $fped;
system(qq|triodenovo --minDQ 0 --minDepth 0 --ped $wrkdir/input.ped --in_vcf $wrkdir/input.vcf --out_vcf $wrkdir/output.vcf|);
collect_dq(\%dq, "$wrkdir/output.vcf");

if (count_line("$wrkdir/input_maleX.vcf")>1) {
	open my $fped2, ">$wrkdir/input_maleX.ped" or die "Cannot write to $wrkdir/input_maleX.ped";
	print $fped2 <<EOF;
Trio Offspring Father Mother 1
Trio Father	0	0	1
Trio Mother	0	0	2
EOF
	close $fped2;
	system(qq|triodenovo --minDQ 0 --minDepth 0 --ped $wrkdir/input_maleX.ped --in_vcf $wrkdir/input_maleX.vcf --out_vcf $wrkdir/output_maleX.vcf|);
	collect_dq(\%dq, "$wrkdir/output_maleX.vcf")
}

($it, $fnames) = iter_file($ARGV{'--input'});
my @fields = @$fnames;


my (%mutrate, $sq, $refchr, $motiflen, $adjlen);
if ($ARGV{'--mutrate'}) {
	%mutrate = map { (split)[0,1] } slurp $ARGV{'--mutrate'};
	my $context = peek_hash(\%mutrate); 
	$motiflen = length($context);
	if ($motiflen % 2 != 1) {
		die "Incorrect context motif length: $motiflen";
	}
	$adjlen = ($motiflen-1)/2;

	if ($ARGV{'--seq'} =~ /\.2bit$/) {
        $sq = Genome::UCSC::TwoBit->new($ARGV{'--seq'});
        my $firstchr =  peek_hash($sq->{SEQLEN});
		if ($firstchr =~ /^chr/) {
			$refchr = 1;
		}
		else {
			$refchr = 0;
		}
	}
	else {
        #croak "Currently only support .2bit file";
        $sq = FaSlice->new(file => $ARGV{'--seq'});
        open my $fin, $ARGV{'--seq'} or die "Cannot open $ARGV{'--seq'}";
		my $firstline = <$fin>;
		if ($firstline =~ /^>chr/) {
			$refchr = 1;
		}
		else {
			$refchr = 0;
		}
	}
	unless($refchr != $inchr) {
		warn "Chromosoome nomenclauture in variant table and refseq are different!";
	}
	unless(@addfields == 2) {
		die "Must provide two additional fields when --mutrate and --seq are provided."
	}
	push @fields, @addfields;
}
else {
	push @fields, $addfields[0];
}


# Appending DenovoBF,DenovoPP at the end of input file
open my $fout, ">$ARGV{'--output'}" or die "Cannot write to $ARGV{'--output'}";
print $fout join("\t", @fields), "\n";
while(my $dat = $it->()) {my $chrom = $dat->{Chrom};
	my $varid = join(":", @{$dat}{qw|IID Chrom Position Ref Alt|});
	my $dnBF = $dq{$varid} // do { die "Cannot find dnBF for $varid" };
	$dat->{$addfields[0]} = $dnBF;

	if (defined $ARGV{'--mutrate'}) {
		# find prior mut rate
		if ($inchr == 0 && $refchr == 1) {
			$chrom = hg_chrom($inchr);
		}
		elsif ($inchr == 1 && $refchr == 0) {
			$chrom =~ s/^chr//; $chrom = 'MT' if $chrom eq 'M'; 
		}
		my $prior;
		if ($dat->{Ref} =~ /^[ACGT]$/ && $dat->{Alt} =~ /^[ACGT]$/) {
			my $context = $sq->get_slice($chrom, $dat->{Position}-$adjlen, $dat->{Position}+$adjlen);
			if ($context =~ /^[ACGT]+$/) {
				unless(defined $mutrate{$context}) {
					die "Cannot find prior mutation rate for $context";
				}
				$prior = $mutrate{$context};
			}
			elsif ($context =~ /^[ACGTN]+$/) {
				my @contexts;
				for(my $ii = 0; $ii < length($context); $ii ++) {
					my $currbp = substr($context, $ii, 1);
					if ($currbp eq 'N') {
						push @contexts, [qw|A C G T|];
					}
					else {
						push @contexts, [$currbp];
					}
				}
				my @comb_contexts = all_combs(@contexts);
				unless(all { defined $mutrate{$_} } @comb_contexts) {
					die "Cannot find prior mutation rates for possible combinations of $context";
				}
				$prior = mean(map { $mutrate{$_} } @comb_contexts);
			}
			else {
				die "Cannot recognize sequence context: $context";
			}
		}
		else {
			$prior = $ARGV{'--indelrate'};
		}
		$dat->{$addfields[1]} = $prior * 10**$dnBF / (1-$prior+ $prior * 10**$dnBF);
	}

	print $fout join("\t", @{$dat}{@fields}), "\n";
}


sub collect_dq {
	my ($dq, $outvcf) = @_;
	my $it = iter_file($outvcf, { ignore => qr/^##/, fsep => qr/\t/ });
	while(my $dat = $it->()) {
		my @keys = split(':', $dat->{FORMAT});
		my @vals = split(':', $dat->{Offspring});
		my %gdat;
		@gdat{@keys} = @vals;
		$dq->{$dat->{ID}} = $gdat{DQ};
	}
	return $dq;
}


__END__

=head1 NAME

eval_dnvs -- Evaluate the confidence for DNVs using triodenovo

=head1 DESCRIPTION

triodenovo calculate a Bayes factor (BF) for each DNV based on the genotype likelihood (PLs) for offspring
and parents. The latest version of denovo (0.06) no longer use AF information, but instead assums prior for 
variants in parents is 1e-3 for SNVs and 1e-4 for indels. So DNV candidates should be filtered to remove high 
freq variants present in the population. 

We also support the calculation of posterior probability of de novo from prior mutation rates under different
sequence contexts.

=head1 REQUIRED ARGUMENTS

=over

=item -[-]input [=] <dnvtab>

Input DNV table, should be DNV finder's output.

Require the following fields present in the input: DP, GQ, PL.HOMREF, PL.HET, PL.HOMALT for proband and parents.

=item -[-]output [=] <output>

The output file, appending additional fields to the input.

=back

=head1 OPTIONS

=over

=item -[-]mutrate [=] <table>

The context dependent SNV mutation rate table, with two columns: Context MutRate. 
It should be absolute rate averaged over alt alleles.

=for Euclid:
	table.default: "$ENV{HOME}/Dropbox/Data/Genetics/AvgRate_3merDenovoNear.txt"

=item -[-]indelrate [=] <rate>

Specify default rates for insertions/deletions/MNVs, default: 1.3e-9.

=for Euclid:
	rate.default: 1.3e-9

=item -[-]seq [=] <fasta>

Genome reference sequence in FASTA or 2bit format. It muse be provided if we need to lookup 
sequence context to determine mutation rate.

=item -[-]wrkdir [=] <wrkdir>

Temp working directory for storing intermediate files.

=item -[-]field-add [=] <string>

The additional field name, default: DenovoBF,DenovoPP.

=for Euclid:
	string.default: "DenovoBF,DenovoPP"

=back




