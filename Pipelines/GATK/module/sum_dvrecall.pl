#!/usr/bin/env perl

use strict;
use warnings;
use Carp;
use IO::Dir;
use FindBin qw|$Bin|;
use List::MoreUtils qw|all uniq|;
use Utils::File::Iter qw|iter_file|;
use Getopt::Euclid;
use Perl6::Slurp;
use FaSlice;
use Genet::Var qw|normalize|;
use Genome::UCSC::BinKeeper;

use lib "$Bin/../../lib/";
use Shared qw|parse_fstr fam_rels|;

# Parse input table file.
my $finput = parse_fstr($ARGV{'--fields'}, 1);
my @fvar = keys %$finput;
my $input_header = 1;
if (all { /^\d+$/ } @fvar) {
	$input_header = 0;
}

my ($it, $fnames) = iter_file($ARGV{'--input'}, { header => $input_header, fsep => qr/\t/ });
foreach my $field (@fvar) {
	unless(grep { $field eq $_ } @$fnames) {
		croak "Cannot find field $field in the input file";
	}
}

if (defined $ARGV{'--oldvar'}) {
	unless(grep { $_ eq $ARGV{'--oldvar'} }  @$fnames) {
		croak "Cannot find oldvar field $ARGV{'--oldvar'} in the input file";
	}
	#@fields = grep { $_ ne $ARGV{'--oldvar'} } @$fnames;
}
my @fields = @$fnames;

my $refseq = FaSlice->new(file => $ARGV{'--fasta'});

# Collect DV results
my @dv_vcfs = grep { /\.vcf\.gz$/ } IO::Dir->new($ARGV{'--dvout'})->read();
my (%vars, %bks);
foreach my $dv_vcf (@dv_vcfs) {
	(my $iid = $dv_vcf) =~ s/\.vcf\.gz$//;
	my ($varcall, $bk) = collect_var("$ARGV{'--dvout'}/$dv_vcf", $refseq);
	$vars{$iid} = $varcall;
	$bks{$iid} = $bk;
}

my @xtrafields = split(',', $ARGV{'--fields-add'});
unless(@xtrafields == 3) {
	die "The length of extra fields must be 3: Filter,Qual,Nearby";
}

my (%fid, $famsamps, $famrels);
if (defined $ARGV{'--ped'}) {
	%fid = map { (split)[1,0] } slurp $ARGV{'--ped'};
	($famsamps, $famrels) = fam_rels($ARGV{'--ped'}, 
		{ ignore => $ARGV{'--ped-ignore'}, twins => $ARGV{'--ped-twins'},
		  shorten => 1,  strict => 0, verbose => 0 }); 
	push @xtrafields, "FamMembers.$xtrafields[0]PASS", "FamMembers.$xtrafields[0]RefCall", "FamMembers.$xtrafields[1]";
}

my $fout;
if ($ARGV{'--output'}) {
	$fout = IO::File->new($ARGV{'--output'}, "w");
}
else {
	$fout = \*STDOUT;
}

if ($input_header) {
	print $fout join("\t", @fields, @xtrafields), "\n";
}

while(my $dat = $it->()) {
	# Restore genome coordinates for the original variant
	# Then store the current variant to the field
	# typically used when liftover was enabled 
	if ($ARGV{'--oldvar'}) {
		my $currvar = join(":", @{$dat}{@fvar[1..4]});
		@{$dat}{@fvar[1..4]} = split(':', $dat->{$ARGV{'--oldvar'}});
		$dat->{$ARGV{'--oldvar'}} = $currvar;
	}

	my ($iid, $chr, $pos, $ref, $alt) = @{$dat}{@fvar};
	my $varid = join("_", $chr, $pos, $ref, $alt);

	my $var = { CHROM => $chr, POS => $pos, REF => $ref, ALT => $alt };
	my $varnorm = normalize($var, $refseq);
	my $varid2 = join("_", @{$varnorm}{qw|CHROM POS REF ALT|});
	if ($varid2 ne $varid) {
		print STDERR "Normalizing input: $varid => $varid2\n";
	}

	# First look for nearby variants
	my @nearvars;
	if (defined $bks{$iid}) {
	 	@nearvars = map { $_->[2]."[Q".sprintf("%d",$_->[3]+0.5)."]" } grep { $_->[2] ne $varid2 } 
					$bks{$iid}->find_range($varnorm->{CHROM}, $varnorm->{POS}-$ARGV{'--nearby'},
											$varnorm->{POS}+length($varnorm->{REF})-1+$ARGV{'--nearby'});
	}

	my $nearvar;
	if (@nearvars) {
		$nearvar = join(q|,|, @nearvars);
	}
	else {
		$nearvar = ".";
	}

	# Then test if the variant was found in Deepvar output by matching standard IDs
	if (defined $vars{$iid}{$varid2}) {
		my ($filter, $qual) = @{$vars{$iid}{$varid2}};
		if (defined $ARGV{'--ped'}) {
			if (defined $fid{$iid}) {
				my (@carriers, @lkcarriers, @quals);
				foreach my $sampid (@{$famsamps->{$fid{$iid}}}) {
					next if $sampid eq $iid;
					if (defined $vars{$sampid}{$varid2}) {
						my ($filt, $q) = @{$vars{$sampid}{$varid2}};
						if ($filt =~ /PASS/) {
							push @carriers, $famrels->{$iid}{$sampid};
							push @quals, $q;
						}
						elsif ($filt =~ /RefCall/) {
							push @lkcarriers, $famrels->{$iid}{$sampid};
						}
					}
				}
				if (@carriers || @lkcarriers) {
					print $fout join("\t", @{$dat}{@fields}, $filter, $qual, $nearvar, 
									@carriers > 0 ? join(',', @carriers) : $ARGV{'--nastring'}, 
									@lkcarriers > 0 ? join(',', @lkcarriers) : $ARGV{'--nastring'}, 
									@quals > 0 ? join(',', @quals) : $ARGV{'--nastring'}), "\n";
				}
				else {
					print $fout join("\t", @{$dat}{@fields}, $filter, $qual, $nearvar, ($ARGV{'--nastring'}) x 3), "\n";
				}
			}
			else {
				print $fout join("\t", @{$dat}{@fields}, $filter, $qual, $nearvar, ($ARGV{'--nastring'}) x 3), "\n";
			}
		}
		else {
			print $fout join("\t", @{$dat}{@fields}, $filter, $qual, $nearvar), "\n";
		}
	}
	else {
		if (defined $ARGV{'--ped'}) {
			print $fout join("\t", @{$dat}{@fields}, ($ARGV{'--nastring'}) x 2, $nearvar, ($ARGV{'--nastring'}) x 3), "\n";
		}
		else {
			print $fout join("\t", @{$dat}{@fields}, ($ARGV{'--nastring'}) x 2, $nearvar), "\n";
		}
	}
}


# All variants found in VCF will be normalized
sub collect_var {
	my ($vcf, $refseq) = @_;
	
	my $bk = Genome::UCSC::BinKeeper->new();
	my %vars;

	my ($it, $fnames) = iter_file($vcf, { fsep => qr/\t/, ignore => qr/^##/ });
	while(my $dat = $it->()) {
		#next unless $dat->{FILTER} qw 'PASS';
		my @alts = grep { /^[ACGTacgt]+$/ } split(',', $dat->{ALT});
		if (@alts > 1) {
			$dat->{FILTER} .= "_MultiAllelic";	
		}
	
		foreach my $alt (@alts) {
			my $varnorm = normalize({CHROM => $dat->{'#CHROM'}, POS => $dat->{POS},
								 	 REF => $dat->{REF}, ALT => $alt}, $refseq);
			my $vid = join("_", @{$varnorm}{qw|CHROM POS REF ALT|});
			# Only variant PASS filter will be kept in nearby vars
			if ($dat->{FILTER} =~ /PASS/) {
				$bk->add($varnorm->{CHROM}, $varnorm->{POS}, $varnorm->{POS}+length($varnorm->{REF})-1, 
						 $vid, $dat->{QUAL});
			}
			$vars{$vid} = [$dat->{FILTER}, $dat->{QUAL}];	
		}
	}

	return (\%vars, $bk);
}


__END__

=head1 NAME

sum_dvrecall.pl -- Summarize deepvariant recall results.

=head1 REQUIRED ARGUMENTS
 
=over

=item -[-]in[put] [=] <file>

Input variant table file.

=for Euclid:
	file.type: readable

=item -[-]dv[out] [=] <dir>

Deep variant output directory. Looking for variants.tfrecord.gz and vcf.gz

=item -[-]fasta [=] <seqfile>

Reference genome sequence file. Used for normalizing indels.

=for Euclid:
	seqfile.type: readable

=back

=head1 OPTIONS

=over

=item -[-]fields [=] <string>

Field string for parsing the variant table file. Should be "Chrom,Position,Ref,Alt"

=for Euclid:
	string.default: "IID,Chrom,Position,Ref,Alt"

=item -[-]ped [=] <pedfile>

Pedigree file. If provided, will also output DV results on all family members.

=item -[-]ped-ignore [=] <pattern>

Pattern for duplicated samples.

=for Euclid:
	pattern.default: '_Re\d*$'

=item -[-]ped-twins [=] <list>

A list of twin pairs.

=item -[-]out[put] [=] <file>

Output table file. If not provided, results will be written to stdout.

=item -[-]fields-add [=] <string>

Output field strings for filter, quality, and neraby variants (comma separated).
Filter and quality will also be given for family members: FamMembers.DeepvarFilterPASS, 
FamMembers.DeepvarQual.

=for Euclid:
	string.default: "DeepvarFilter,DeepvarQual,NearbyVars"

=item -[-]nastring [=] <string>

=for Euclid:
	string.default: "."

=item -[-]oldvar [=] <field>

Restore genome coordinates to the original variant given by this field. Represented by "Chrom:Pos:Ref:Alt".
The replace the current coordinates into this field.

=item -[-]nearby [=] <distance>

Flag nearby variants within this distance cutoff, default: 2bp.
Flag will be merged with the filter field. 
Even for variants failed DeepVar, nearby variants will also be written to the output.

=for Euclid:
	distance.default: 2

=back

