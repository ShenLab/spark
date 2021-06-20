#!/usr/bin/env perl

use strict;
use warnings;
use Carp;
use FindBin qw|$Bin|;
use List::Util qw|first|;
use Perl6::Slurp;
use Data::Dumper;
use Getopt::Euclid;
use Utils::Seq qw|iub3to1|;
use List::MoreUtils qw|any all uniq|;
use Utils::Hash qw|str2hash chk_default|;
use Utils::File::Iter qw|iter_file|;

use lib "$Bin/../../lib";
use Shared qw|parse_fstr|;
use Variants qw|%eff|;


# Functional consequences and their ranks are given in this web page
# https://useast.ensembl.org/info/genome/variation/prediction/predicted_data.html#consequences
# In the ranks below, we modified default ranking so the synonymous has higher priority than splice_region.
# For details see lib/Variants.pm

# It is possible for a variant to have multiple effects even for one transcript.
# E.g. splice_region can appear together with missense or synonymous or intronic.
# We also noted that functional consequences below ignored transcript biotype, e.g. splice site variant of a 
# noncoding transcript will be ranked as high as that of coding transcript. But when it affects noncoding
# transcript, non_coding_transcript will also shown.
# => So when comparing the priority of annotation terms, we need take into account all effects or biotypes

# To resolve these issues, we implemented following strategies in the current script:
# 1. When "--coding" option is switched on, only protein coding transcripts will be used in annotation
# 2. When multiple functional effects are associated with a variant-transcript combination, the default is
# to select one with higher rank. But when NMD_transcript or non_coding_transcript appear, they will be 
# given a higher priority regarless of other effects.
# 3. When biotype can be found in the output, we will show biotype for each transcript.
# 4. By default, we select one functional effect for each variant-transcript combination (with priority given by
# the effect ranks pr by Rule 3 above), we can also output full effects by '--full-effect'.

my ($flag_symbol, $flag_biotype, $flag_canon, $flag_hgvs, $flag_ancestral, $flag_context, $flag_splicereg);
my %known;
{
	# Parse the header from VEP output
	open my $fin, $ARGV{'--input'} or die "Cannot open input file $ARGV{'--input'}";
	while(<$fin>) {
		last unless /^##/;
		my $field = (split)[1];
		$known{$field} = 1;
		if ($field eq 'SYMBOL') {
			$flag_symbol = 1;
		} elsif ($field eq 'BIOTYPE') {
			$flag_biotype = 1;
		} elsif ($field eq 'CANONICAL') {
			$flag_canon = 1;
		} elsif ($field eq 'HGVSc' || $field eq 'HGVSp') {
			$flag_hgvs = 1;
		} elsif ($field eq 'AncestralAllele') {
			$flag_ancestral = 1;
		} elsif ($field eq 'Context') {
			$flag_context = 1;
		} elsif ($field eq 'SpliceRegion') {
			$flag_splicereg = 1;
		} 
	}
}


my $it = iter_file($ARGV{'--input'}, { fsep => qr/\t/, ignore => qr/^##/});

my (%inclgene, %exclgene, %incltrans, %excltrans);
if ($ARGV{'--include-gene'}) {
	%inclgene = map { (split)[0] => 1 } slurp $ARGV{'--include-gene'};
}
if ($ARGV{'--exclude-gene'}) {
	%exclgene = map { (split)[0] => 1 } slurp $ARGV{'--exclude-gene'};
}
if ($ARGV{'--include-trans'}) {
	%incltrans = map { (split)[0] => 1 } slurp $ARGV{'--include-trans'};
}
if ($ARGV{'--exclude-trans'}) {
	%excltrans = map { (split)[0] => 1 } slurp $ARGV{'--exclude-trans'};
}

my (%extragene, %extratrans, %extravar);
# Extra variant level information (mostly population allele freq) should have the same value for each chrom-pos-ref-alt combination
# We also support extra gene-variant annotation (e.g. PolyPhen,SIFT), that have one value for each chrom-pos-ref-alt plus gene ID combination
if (defined $ARGV{'--extra-gene'}) {
	my $fields = parse_fstr($ARGV{"--extra-gene"}, 1);
	unless(all { defined $known{$_} } keys %$fields) {
		die "Not all extra-gene fields can be found in the input";
	}
	%extragene = ( infields => [keys %$fields], outfields => [values %$fields] );
}
if (defined $ARGV{'--extra-trans'}) {
	my $fields = parse_fstr($ARGV{"--extra-trans"}, 1);
	unless(all { defined $known{$_} } keys %$fields) {
		die "Not all extra-trans fields can be found in the input";
	}
	%extratrans = ( infields => [keys %$fields], outfields => [values %$fields] );
}
if (defined $ARGV{'--extra-var'}) {
	my $fields = parse_fstr($ARGV{"--extra-var"}, 1);
	unless(all { defined $known{$_} } keys %$fields) {
		die "Not all extra-var fields can be found in the input";
	}
	%extravar = ( infields => [keys %$fields], outfields => [values %$fields] );
}

# Standard and optional output fields:
# VarID Chrom Position Ref Alt (Anc) (Context) (Symbol) GeneID GeneEff TransCount (TransCanon) TransIDs (TransBiotypes) TransEffs (SpliceReg) CodonChg (cDNAChg) AAChg
my @outfields = qw|VarID Chrom Position Ref Alt GeneID GeneEff TransCount TransIDs TransEffs CodonChg AAChg|;
if ($flag_symbol) {
	my $ii = first { $outfields[$_] eq 'GeneID' } 0..$#outfields;
	splice @outfields, $ii, 0, 'Symbol';
}
if ($flag_canon) {
	my $ii = first { $outfields[$_] eq 'TransIDs' } 0..$#outfields;
	splice @outfields, $ii, 0, 'TransCanon';
}
if ($flag_context) {
	my $ii = first { $outfields[$_] eq 'Alt' } 0..$#outfields;
	splice @outfields, $ii+1, 0, 'Context';
}
if ($flag_ancestral) {
	my $ii = first { $outfields[$_] eq 'Alt' } 0..$#outfields;
	splice @outfields, $ii+1, 0, 'Anc';
}
if ($flag_biotype) {
	my $ii = first { $outfields[$_] eq 'TransIDs' } 0..$#outfields;
	splice @outfields, $ii+1, 0, 'TransBiotypes';
}
if ($flag_splicereg) {
	my $ii = first { $outfields[$_] eq 'TransEffs' } 0..$#outfields;
	splice @outfields, $ii+1, 0, 'SpliceReg';
}
if ($flag_hgvs) {
	# When HGVS is available, cDNA changes in HGVS format will be added after codon changes
	# And AA changes will be in HGVS format
	my $ii = first { $outfields[$_] eq 'CodonChg' } 0..$#outfields;
	splice @outfields, $ii, 0, 'cDNAChg';
}
if (defined $extragene{outfields}) {
	my $ii = first { $outfields[$_] eq 'TransIDs' } 0..$#outfields;
	splice @outfields, $ii, 0, @{$extragene{outfields}};
}
# Extra transcript and variant level annotations will be appended to the end of standard/optional annotations
if (defined $extratrans{outfields}) {
	push @outfields, @{$extratrans{outfields}};
}
if (defined $extravar{outfields}) {
	push @outfields, @{$extravar{outfields}};
}



if (-f "$ARGV{'--output'}.done") {
	unlink("$ARGV{'--output'}.done");
}

my $fout = IO::File->new($ARGV{'--output'}, "w");
print $fout join("\t", @outfields), "\n";

my $prevar;
# ID of previous coding variants
my %anno; 
while(my $dat = $it->()) {
	# For gene-based annotation, feature type must be transcript
	# Currently we do not use VEP for another types of annotation
	next unless $dat->{Feature_type} eq 'Transcript';

	my %info = str2hash($dat->{Extra});
	
	# Under coding mode, it will restrict to coding genes and polymorphic_pseudogenes
	# otherwise, all types of transcripts will be kept
	if ($ARGV{'--coding'}) {
		die "Must have biotypes in the VEP output to select coding transcript" unless $flag_biotype;
		next unless $info{BIOTYPE} eq 'protein_coding' || $info{BIOTYPE} eq 'polymorphic_pseudogene';
	}
	
	# Skip excluded genes or transcript
	if ($ARGV{'--include-gene'}) {
		next unless defined $inclgene{$dat->{Gene}};
	}
	if ($ARGV{'--exclude-gene'}) {
		next if defined $exclgene{$dat->{Gene}};
	}
	if ($ARGV{'--include-trans'}) {
		next unless defined $incltrans{$dat->{Feature}};
	}
	if ($ARGV{'--exclude-trans'}) {
		next if defined $excltrans{$dat->{Feature}};
	}
	
	# When a new variants is found, output previous var
	if (defined $prevar && $dat->{'#Uploaded_variation'} ne $prevar) {
		# output one line of variant
		print $fout $prevar, "\t", refmt_varid($prevar), "\t", refmt_outline(\%anno), "\n";
		%anno = ();
	}

	###########################
	### Process annotation  ###
	###########################
	if ($flag_ancestral && !defined $anno{anc}) {
		$anno{anc} = $info{AA};
	}

	# Context will be fetched for SNVs only
	# Variant ID must be properly formatted to get ref and alt alleles
	if ($flag_context && !defined $anno{context}) {
		my ($ref, $alt) = (split(':', $dat->{'#Uploaded_variation'}))[2,3];
		if ($ref =~ /^[ACGT]$/ && $alt =~ /^[ACGT]$/) {
			my $length = length($info{Context});
			my $midlen = int($length/2);
			croak "Incorrect length of context: $length" unless $length % 2 == 1; 
			croak "Incorrect seq of context: $info{Context} for variant $dat->{'#Uploaded_variation'}" 
				unless substr($info{Context}, $midlen, 1) eq $ref;
			my $context_alt = $info{Context};
			substr($context_alt, $midlen, 1, $alt);
			$anno{context} = $info{Context}.'>'.$context_alt;
		}
		else {
			$anno{context} = '.';
		}
	}

	push @{$anno{gene}} => $dat->{Gene};
	if ($flag_symbol) {
		# $anno{symbol}{$dat->{Gene}} = $info{SYMBOL};
		$anno{symbol} = {} unless defined $anno{symbol};
		chk_default($anno{symbol}, $dat->{Gene}, $info{SYMBOL})
	}
	
	push @{$anno{trans}{$dat->{Gene}}} => $dat->{Feature};
	if ($flag_canon) {
		if (defined $info{CANONICAL} && $info{CANONICAL} eq 'YES') {
			#$anno{canon}{$dat->{Gene}} = $dat->{Feature};
			$anno{canon} = {} unless defined $anno{canon};
			chk_default($anno{canon}, $dat->{Gene}, $dat->{Feature});
		}
	}
	if ($flag_biotype) {
		push @{$anno{biotype}{$dat->{Gene}}} => $info{BIOTYPE};
	}

	# Whem multiple functional effects are present, we will pickup the *first* one
	# based on the rank shown above. For example, if a variant was annoteted as 
	# both missense and splice_region, then only missense will be shown in the output.
	# But when the effects was on non-coding transcript, we will give priority to
	# noncoding transcript to avoid confoundings in intrepretation.

 	my @csqs = sort { $eff{$a}[1] <=> $eff{$b}[1] } 
 		map { s/_variant$//; $_ } split(',', $dat->{Consequence});

	if (@csqs == 1) {
		push @{$anno{csq}{$dat->{Gene}}} => $csqs[0];
	} 
	elsif (@csqs > 1) {
		if (any { /^(non_coding|NMD)_transcript/ } @csqs) {
 			my $ii = first { $csqs[$_] =~ /^(non_coding|NMD)_transcript/ } 0..$#csqs;
 			if ($ii != 0) {
 				unshift @csqs, $csqs[$ii];
 				splice @csqs, $ii+1, 1;
 			}
 		}
 		if ($ARGV{'--full-effect'}) {
 			push @{$anno{csq}{$dat->{Gene}}} => $csqs[0].'('.join(",",@csqs[1..$#csqs]).')';
 		}
 		else {
 			push @{$anno{csq}{$dat->{Gene}}} => $csqs[0];
 		}

	}
	else {
		die "Incorrect number of consequences";
	}
	push @{$anno{effrank}{$dat->{Gene}}} => $eff{$csqs[0]}[1];

	if ($flag_splicereg) {
		if (defined $info{SpliceRegion}) {
			my @splicepos = split(',', $info{SpliceRegion});
			if (@splicepos > 1) {
				push @{$anno{splice}{$dat->{Gene}}} => $splicepos[0].'~'.$splicepos[-1];
			}
			else  {
				push @{$anno{splice}{$dat->{Gene}}} => $info{SpliceRegion};
			}	
		}
		else {
			push @{$anno{splice}{$dat->{Gene}}} => '.';
		}
	}
	else {
		if (any { $_ =~ /splice/ } @csqs) {
			if (any { $_ =~ /splice_donor/ } @csqs) {
				push @{$anno{splice}{$dat->{Gene}}} => 'D';
			}
			elsif (any { $_ =~ /splice_acceptor/ } @csqs) {
				push @{$anno{splice}{$dat->{Gene}}} => 'A';
			}
			else {
				push @{$anno{splice}{$dat->{Gene}}} => 'R';
			}	
		} 
		else {
			push @{$anno{splice}{$dat->{Gene}}} => '.';
		}
	}
	
	if ($flag_hgvs) {
		if (defined $info{HGVSc}) {
			my @hgvsc = split(':', $info{HGVSc});
			croak "Incorrect HGVSc nomenclature" unless @hgvsc == 2;
			push @{$anno{hgvsc}{$dat->{Gene}}} => $hgvsc[1];	
		}
		else {
			push @{$anno{hgvsc}{$dat->{Gene}}} => '.';
		}
		if ($dat->{CDS_position} ne '-' && $dat->{Codons} ne '-') {
			push @{$anno{codon}{$dat->{Gene}}} => $dat->{CDS_position}.":".$dat->{Codons};
		}
		else {
			push @{$anno{codon}{$dat->{Gene}}} => '.';
		}
		if (defined $info{HGVSp}) {
			my @hgvsp = split(':', $info{HGVSp});
			croak "Incorrect HGVSp nomenclature" unless @hgvsp == 2;
			push @{$anno{hgvsp}{$dat->{Gene}}} => abbrev_hgvsp($hgvsp[1]);
		}
		else {
			push @{$anno{hgvsp}{$dat->{Gene}}} => '.';
		}
	}
	else {
		if ($dat->{CDS_position} ne '-' && $dat->{Codons} ne '-') {
			#$dat->{Codons} =~ s/\//>/;
			push @{$anno{codon}{$dat->{Gene}}} => $dat->{CDS_position}.":".$dat->{Codons};
		}
		else {
			push @{$anno{codon}{$dat->{Gene}}} => '.';
		}
		if ($dat->{Protein_position} ne '-' && $dat->{Amino_acids} ne '-') {
			#$dat->{Amino_acids} =~ s/\//>/;
			push @{$anno{hgvsp}{$dat->{Gene}}} => $dat->{Protein_position}.":".$dat->{Amino_acids};
		}
		else {
			push @{$anno{hgvsp}{$dat->{Gene}}} => '.';
		}
	}

	if (defined $extravar{infields}) {
		$anno{extra_var} = {} unless defined $anno{extra_var};
		foreach my $key (@{$extravar{infields}}) {
			chk_default($anno{extra_var}, $key, $info{$key} // ".");
		}
	}
	if (defined $extragene{infields}) {
		$anno{extra_gene}{$dat->{Gene}} = {} unless defined $anno{extra_gene}{$dat->{Gene}};
		foreach my $key (@{$extragene{infields}}) {
			chk_default($anno{extra_gene}{$dat->{Gene}}, $key, $info{$key} // ".");
		}
	}
	if (defined $extratrans{infields}) {
		foreach my $key (@{$extratrans{infields}}) {
			push @{$anno{extra_trans}{$dat->{Gene}}{$key}}, $info{$key} // ".";
		}
	}

	$prevar = $dat->{'#Uploaded_variation'};
}

# Output the last variant
if (defined $prevar) {
	print $fout $prevar, "\t", refmt_varid($prevar), "\t", refmt_outline(\%anno), "\n";
	$fout->close;
}

open my $flog, ">$ARGV{'--output'}.done" or die "Cannot write to flag file";
print $flog "OK\n";
close $flog;

sub refmt_varid {
	my ($varid) = @_;
	return join("\t", split(':', $varid));
}

# Reformat the output, VarID will be concat later
# gene symbols, genes IDs, gene level effect, number of transcripts with annotations,
# transcript IDs, conequences, cDNA changes, AA changes
sub refmt_outline {
	my ($anno) = @_;
	my @line;

	if ($flag_ancestral) {
		push @line, uc($anno->{anc});
	}
	if ($flag_context) {
		push @line, $anno->{context};
	}
	
	# Decide which genes to be included in the final list
	my @genes;
	foreach my $gene (uniq sort @{$anno->{gene}}) {
		if ( any { $_ <= $ARGV{'--cutoff'} } @{$anno->{effrank}{$gene}} ) {
			push @genes, $gene;
		}
	}
	if (@genes) {
		if ($flag_symbol) {
			my @uqsym = uniq @{$anno->{symbol}}{@genes};
			if (@uqsym == 1) {
				push @line, $uqsym[0];
			}
			else {
				push @line, join(";", @{$anno->{symbol}}{@genes});
			}
		}
		push @line, join(";", @genes);

		# summarize gene level effects, pick most severe for each gene
		push @line, gene_eff([@{$anno->{csq}}{@genes}], [@{$anno->{effrank}}{@genes}]);
		# the number of annotated transcripts (and optional canonical transcripts)
		# Check that transcripts for each genes are not redundant
		foreach my $gene (@genes) {
			my @transids = @{$anno->{trans}{$gene}};
			my @uqtids = uniq sort @transids;
			unless(@transids == @uqtids) {
				die "Redundant transcript IDs found for $gene";
			}
		}
		push @line, join(";", map { scalar(@$_) } @{$anno->{trans}}{@genes});
		if ($flag_canon) {
			push @line, join(";", map { $_ // "." } @{$anno->{canon}}{@genes});
		}
		if (defined $extragene{infields}) {
			foreach my $key (@{$extragene{infields}}) {
				push @line, join(";", map { $_->{$key} =~ s/;/|/g; $_->{$key} } @{$anno->{extra_gene}}{@genes});
			}
		}

		foreach my $key (qw|trans biotype csq splice hgvsc codon hgvsp|) {
			if ($key eq 'biotype') {
				next unless $flag_biotype;
			}
			if ($key eq 'splice') {
				next unless $flag_splicereg;
			}
			if ($key eq 'hgvsc') {
				next unless $flag_hgvs;
			}
			push @line, join(";", map { join(',', @$_) } @{$anno->{$key}}{@genes});
		}
		if (defined $extratrans{infields}) {
			foreach my $key (@{$extratrans{infields}}) {
				push @line, join(";", map { join(',', map { $_ =~ s/,/|/g; $_ } @{$_->{$key}}) } @{$anno->{extra_trans}}{@genes});
			}
		}
	}
	else {
		my $nxtra = grep { defined $_ } ($flag_symbol, $flag_canon, $flag_biotype, $flag_splicereg, $flag_hgvs);
		push @line, ('.') x (7 + $nxtra);
		if (defined $extragene{infields}) {
			push @line, ('.') x scalar(@{$extragene{infields}});
		}
		if (defined $extratrans{infields}) {
			push @line, ('.') x scalar(@{$extratrans{infields}});
		}
	}
	
	# Extra variant level annotations independent of gene or transcript
	foreach my $key (@{$extravar{infields}}) {
		push @line, $anno->{extra_var}{$key};
	}

	return join("\t", map { $_ // "." } @line);
}

# Summarize gene level effect, by taking the most severe consequnce based on customized ranking
sub gene_eff {
	my ($csqs, $effranks) = @_;
	unless(@$csqs == @$effranks) {
		die "Incorrect length of csqs and effranks";
	}
	# will pick most servere consequence for each gene
	my @csqout;
	for(my $ii = 0; $ii < @$csqs; $ii ++) {
		my $csq = $csqs->[$ii];
		my $effrank = $effranks->[$ii];
		my @csqrnk = sort { $effrank->[$a] <=> $effrank->[$b] } 0..(@$csq-1);
		push @csqout, $csq->[$csqrnk[0]];
	}	
	if (@csqout == 1) {
		return $csqout[0];	
	}
	else {
		return join(";", @csqout);
	}
}

# Abbreviate HGVSp representation
sub abbrev_hgvsp {
	my ($hgvsp) = @_;
	$hgvsp =~ s/Ter/*/;
	$hgvsp =~ s/%3D/=/;
	my $aa3regex = qr/Ala|Val|Leu|Ile|Pro|Trp|Phe|Met|Gly|Ser|Thr|Tyr|Cys|Asn|Gln|Lys|Arg|His|Asp|Glu|Sec/;
	while($hgvsp =~ /($aa3regex)/) {
		my $aa3 = $1;
		my $aa1 = iub3to1(uc($aa3), 1);
		$hgvsp =~ s/$aa3/$aa1/;
	}
	return $hgvsp;
}

__END__

=head1 NAME

collect_vepout_coding.pl -- Collect and reformat VEP annotation for each sequence variant.

=head1 DESCRIPTION

We mainly use VEP for annotating functional consequences of coding sequence variants.
Wehn parsing the VEP output, we will only extract those with "Feature_type == Transcript",
(RegulatoryFeature, MotifFeature will be ignored.

The script will stack all annotation for each variant in one line. The Standard and optional 
output fields are listed below. When fields for extra variant information were provided, they 
will be inserted into the appropritate:
	VarID Chrom Position Ref Alt (Anc) (Context) (Symbol) GeneID GeneEff [extra-gene...]
	TransCount (TransCanon) TransIDs (TransBiotypes) TransEffs (SpliceReg) CodonChg (cDNAChg) AAChg 
	[extra-trans...] [extra-var...]
Where
	Anc is ancestral allele by AncestralAllele plugin
	Context is tri-nucleotide sequence context for SNV by Context plugin
	Symbol is gene symbol when --symbol is turned on, when overlapping genes are associated with the same symbol, it will be collapsed into one
	TransCanon is when --canonical is turned on
	SpliceReg: by inhouse SpliceRegion plugin
	cDNAChg: is HGVSc when --hgvs is turned on (and AAChg will be HGVSp in this case)

	 
Note: extra-gene/transcript level information is typically associated with each variant.
If they are variant-independent, they can also beetched from a lookup table. 

To get get and alt allele, uploaded variant ID MUST be in 'Chr:Pos:Ref:Alt' format.

=head1 REQUIRED ARGUMENTS

=over

=item -[-]in[put] [=] <file>

The standard VEP output file.

=for Euclid:
	file.type: readable 

=item -[-]out[put] [=] <file>

The reformatted output file.

=back

=head1 OPTIONS

=over

=item -[-]coding

Restrict annotations to only coding transcripts (biotype=="protein_coding" or "polymorphic_pseudogene").

=item -[-]full-effect

Output full functional effect. The default is to select one for each variant-transcript combination.
Under full-effect mode, the un-selected effect(s) will be shown in bracket.

=item -[-]cutoff [=] <threshold>

Max. rank for functional effects to be appear in the GeneEff, default: 18.
Note: all possible functional effects for that gene will still appear in TransEffs column regardless of cutoff.

=for Euclid:
	threshold.default: 18

=item -[-]include-gene [=] <genelist>

A gene IDs list, to be included. Can be used together with --exclude-gene.

=item -[-]exclude-gene [=] <genelist>

A gene IDs list, to be excluded. 
To avoid cluttered output, some genes or transcripts can be included or excluded from annotations.

=for Euclid:
	genelist.type: readable

=item -[-]include-trans [=] <translist>

A transcript ID list, to be included.

=item -[-]exclude-trans [=] <translist>

A transcript IDs list, to be excluded.

=item -[-]extra-gene [=] <fstr>

Extra gene level information per variant. Similar as above, except that it should have same value for each
unique chr:pos:ref:alt:gene combination

=item -[-]extra-trans [=] <fstr>

Extra transcript level information per variant.

=item -[-]extra-var [=] <fstr>

Extra variant level fields to be collected from VEP output's Extra field. It should have the same value
for each unique chr:pos:ref:alt combination. "fstr" can list extracted field name and (optionally) alias 
used in output. This is mainly for extracting extra fields annotated by custom tracks.

=back






