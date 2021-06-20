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
#use Utils::Number qw|pretty_bps|;
use Utils::Hash qw|str2hash chk_default|;
use Utils::File::Iter qw|iter_file|;
use Genome::Ranges qw|spec_to_range|;

use lib "$Bin/../../lib";
use Shared qw|parse_fstr|;
use Variants qw|%eff cnv_id cnv_type|;

# This is a modified version of collect_vepout utility for CNV annotations.

my ($flag_symbol, $flag_biotype, $flag_number);
{
	# Parse the header from VEP output
	open my $fin, $ARGV{'--input'} or die "Cannot open input file $ARGV{'--input'}";
	while(<$fin>) {
		last unless /^##/;
		my $field = (split)[1];
		if ($field eq 'SYMBOL') {
			$flag_symbol = 1;
		} 
		elsif ($field eq 'BIOTYPE') {
			$flag_biotype = 1;
		}
		elsif ($field eq 'EXON') {
			unless($ARGV{'--no-exon'}) {
				$flag_number = 1;
			}
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


# Standard and optional output fields:
# VarID Chrom Start End Type GeneCount (Symbols) GeneIDs GeneEffs TransCounts TransIDs (TransBiotypes) TransEffs (TransExons)
my @outfields = qw|VarID Chrom Start End Type GeneCount GeneIDs GeneEffs TransCounts TransIDs TransEffs|;
my @alltransfds = qw|csq|;
unless($ARGV{'--no-overlap'}) {
	push @outfields, 'TransOverlaps';
	push @alltransfds, 'overlap';
}
if ($flag_symbol) {
	my $ii = first { $outfields[$_] eq 'GeneIDs' } 0..$#outfields;
	splice @outfields, $ii, 0, 'Symbols';
}
if ($flag_biotype) {
	my $ii = first { $outfields[$_] eq 'TransIDs' } 0..$#outfields;
	splice @outfields, $ii+1, 0, 'TransBiotypes';
	unshift @alltransfds, 'biotype';
}
if ($flag_number) {
	push @outfields, 'TransExons';
	push @alltransfds, 'exon';
}

if (-f "$ARGV{'--output'}.done") {
	unlink("$ARGV{'--output'}.done");
}

my $fout = IO::File->new($ARGV{'--output'}, "w");
print $fout join("\t", @outfields), "\n";

my ($prevar, $currvar);
# IDs of previous and currents coding variants
my %anno; 
while(my $dat = $it->()) {
	# For gene-based annotation, feature type must be transcript
	# Currently we do not use VEP for another types of annotation
	next unless $dat->{Feature_type} eq 'Transcript';

	my %info = str2hash($dat->{Extra});
	
	# Under coding mode, it will restrict to coding genes and polymorphic_pseudogenes
	# otherwise, all types of transcripts will be kept
	if ($ARGV{'--coding'}) {
		die "Must have biotypes in the VEP output for selecting coding transcripts" unless $flag_biotype;
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
	
	# For CNV, we will reformat the variant ID
	my ($chrom, $start, $end) = spec_to_range($dat->{Location});
	my $currvar = join("\t", cnv_id(spec_to_range($dat->{Location}),$dat->{Allele}), 
								$chrom, $start, $end, cnv_type($dat->{Allele}));
	if (defined $prevar && $currvar ne $prevar) {
		# output one line of variant
		print $fout $prevar, "\t", refmt_outline(\%anno), "\n";
		%anno = ();
	}

	###########################
	### Process annotation  ###
	###########################
	
	push @{$anno{gene}} => $dat->{Gene};
	if ($flag_symbol) {
		# $anno{symbol}{$dat->{Gene}} = $info{SYMBOL};
		$anno{symbol} = {} unless defined $anno{symbol};
		chk_default($anno{symbol}, $dat->{Gene}, $info{SYMBOL})
	}
	
	push @{$anno{trans}{$dat->{Gene}}} => $dat->{Feature};
	if ($flag_biotype) {
		push @{$anno{biotype}{$dat->{Gene}}} => $info{BIOTYPE};
	}
	push @{$anno{overlap}{$dat->{Gene}}} => $info{OverlapPC} // ".";
	if ($flag_number) {
		foreach my $part (qw|exon intron|) {
			next if $part eq 'intron'; # skip intron
			my $PART = uc($part);
			if (defined $info{$PART}) {
				if ($info{$PART} eq '1/1' || $info{$PART} =~ /^1-(\d+)\/\1$/ ) {
					push @{$anno{$part}{$dat->{Gene}}} => 'All';
				}
				else {
					push @{$anno{$part}{$dat->{Gene}}}, $info{$PART};
				}
			}
			else {
				push @{$anno{$part}{$dat->{Gene}}} => '.';
			}
		}
	}

	# Whem multiple functional effects are present, we will pickup the *first* one
	# based on the rank.
 	my @csqs = sort { $eff{$a}[1] <=> $eff{$b}[1] } 
 		map { s/_variant$//; $_ } split(',', $dat->{Consequence});

	if (@csqs == 1) {
		push @{$anno{csq}{$dat->{Gene}}} => $csqs[0];
	} 
	elsif (@csqs > 1) {
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
	# The effect rank will be stored and used later when applying rank cutoff
	push @{$anno{effrank}{$dat->{Gene}}} => $eff{$csqs[0]}[1];

	$prevar = $currvar;
}

# Output the last variant
if (defined $prevar) {
	print $fout $prevar, "\t", refmt_outline(\%anno), "\n";
	$fout->close;
}

open my $flog, ">$ARGV{'--output'}.done" or die "Cannot write to flag file";
print $flog "OK\n";
close $flog;


sub refmt_outline {
	my ($anno) = @_;
	my @line;

	# Decide which genes to be included in the final list
	# Keep the gene as long as effect on any of its transcript fall below the cutoff
	my @genes;
	foreach my $gene (uniq sort @{$anno->{gene}}) {
		if ( any { $_ <= $ARGV{'--cutoff'} } @{$anno->{effrank}{$gene}} ) {
			push @genes, $gene;
		}
	}
	if (@genes) {		
		# Transcript specific annotations
		# First reorder transcript IDs based on the frequency of annotations
		# and also determining the fields that can be collapsed
		my %collapse;
		foreach my $gene (@genes) {
			my @transvals;
			foreach my $transfd (@alltransfds) {
				unless (@{$anno->{trans}{$gene}} == @{$anno->{$transfd}{$gene}}) {
					die "Incorrect length of trans and $transfd for gene $gene";
				}
				for(my $ii = 0; $ii < @{$anno->{trans}{$gene}}; $ii ++) {
					push @{$transvals[$ii]}, $anno->{$transfd}{$gene}[$ii];
				}
				my @vals = @{$anno->{$transfd}{$gene}};
				my @uqvals = uniq @vals;
				if (@uqvals == 1) {
					#if ($transfd eq 'exon') {
					#	$collapse{$gene}{$transfd} = 1 if $uqvals[0] eq 'All';
					#}
					#else {
						$collapse{$gene}{$transfd} = 1;
					#}
				}
			}
			my $transii = order_trans($anno->{trans}{$gene}, \@transvals);
			foreach my $transfd ("trans", "effrank", @alltransfds) {
				$anno->{$transfd}{$gene} = [ @{$anno->{$transfd}{$gene}}[@$transii] ];
			}
		}

		# Sort gene based on the number of identical gene effs and number collapsed fields
		my ($gene_eff, $genes_sorted) = gene_effs(\@genes, $anno, \%collapse);
		push @line, scalar(@$genes_sorted);
		if ($flag_symbol) {
			push @line, join(";", @{$anno->{symbol}}{@$genes_sorted});
		}
		push @line, join(";", @$genes_sorted);

		# Summarize the most severe gene level effects
		push @line, join(";", trunc_end(map { $gene_eff->{$_} } @$genes_sorted));
		# The counts of annotated transcripts 
		push @line, join(";", map { scalar(@$_) } @{$anno->{trans}}{@$genes_sorted});
		
		# Output transcript level fields
		my @transout;
		foreach my $transfd (@alltransfds) {
			my @output;	
			foreach my $gene (@$genes_sorted) {
				if ($collapse{$gene}{$transfd}) {
					push @output, $anno->{$transfd}{$gene}[0];
				}
				else {
					push @output, join(',', trunc_end(@{$anno->{$transfd}{$gene}}));
				}
			}
			push @transout, join(";", trunc_end(@output));
			#push @transout, join(";", @output);
		}

		# If all transcript field for a gene are collapsed, then we can omit the tranIDs for this gene?
		my @transids;
		foreach my $gene (@$genes_sorted) {
			if ( (all { defined $collapse{$gene}{$_} } @alltransfds) ||  @{$anno{trans}{$gene}} == 1 ) {
				unless(defined $ARGV{'--full-trans'}) {
					push @transids, "...";
				}
				else {
					push @transids, join(',', @{$anno{trans}{$gene}});
				}
			}
			else {
				push @transids, join(',', @{$anno{trans}{$gene}});
			}
		}
		push @line, join(';', trunc_end(@transids)); 
		push @line, @transout;
	}
	else {
		# No genes were affected, determine the number of empty fields
		# [VarID Chrom Start End Type] GeneCount (Symbols) GeneIDs GeneEffs TransCounts TransIDs (TransBiotype) TransEffs (TransExons)
		my $nxtra = grep { defined $_ } ($flag_symbol, $flag_biotype, $flag_number);
		if ($ARGV{'--no-overlap'}) {
			push @line, ('.') x (6 + $nxtra);
		}
		else {
			push @line, ('.') x (7 + $nxtra);
		}
	}
	return join("\t", map { $_ // "." } @line);
}

# Summarize gene level effect, by taking the most severe consequnce based on customized ranking
# Return two variables, 1. gene_id => most severe eff, and 2. sorted gene IDs based on freq of eff, 
# genes with less frequent effs (typically genes at CNV boundaries) will appear first 
sub gene_effs {
	my ($genes, $anno, $collapse) = @_;
	my %gene_eff; 
	foreach my $gene (@$genes) {
		my $csq = $anno->{csq}{$gene};
		my $effrank = $anno->{effrank}{$gene};
		unless(@$csq == @$effrank) {
			die "Incorrect length of csq and effrank for $gene";
		}
		my @csqrnk = sort { $effrank->[$a] <=> $effrank->[$b] } 0..(@$csq-1);
		$gene_eff{$gene} = $csq->[$csqrnk[0]];
	}
	my $genes_sorted = sort_genes(\%gene_eff, $collapse);
	return (\%gene_eff, $genes_sorted);
}


# Sort genes based on the frequency of effect, then by the number of collapsed fields
sub sort_genes {
	my ($gene_eff, $collapse) = @_;
	my %eff_count;
	while(my ($gene, $eff) = each %$gene_eff) {
		$eff_count{$eff} ++;
	}
	my @genes_sorted = map { $_->[0] } 
		sort { $eff_count{$a->[1]} <=> $eff_count{$b->[1]} ||
			   scalar(keys %{$collapse->{$a->[0]}}) <=> scalar(keys %{$collapse->{$b->[0]}}) } 
			map { [ $_, $gene_eff->{$_} ] } keys %$gene_eff;
	return \@genes_sorted;
}

# Sort transcript based on the frequency of combined transcript fields
sub order_trans {
	my ($transids, $transvals) = @_;
	unless(@$transids == @$transvals) {
		die "Incorrect of length of transcript IDs and effects";
	}
	my @uqtids = uniq sort @$transids;
	unless(@$transids == @uqtids) {
		print join("\n", @$transids), "\n";
		die "Redundant transcript IDs";
	}
	my %trans_val;
	@trans_val{@$transids} = map { join("\t", @$_) } @$transvals;
	my %val_count;
	while(my ($trans, $val) = each %trans_val) {
		$val_count{$val} ++;
	}
	my @transii_order = sort { $val_count{$trans_val{$transids->[$a]}} <=> 
							   $val_count{$trans_val{$transids->[$b]}} } 0..$#$transids;
	return \@transii_order;
}


# Truncate the array at the end by stop repeating the same content
sub trunc_end {
	my (@dat) = @_;
	return @dat if @dat == 1;
	my $ii;
	for($ii = $#dat; $ii > 0; $ii --) {
		if ($dat[$ii-1] ne $dat[$#dat]) {
			last;
		}
	}
	return @dat[0..$ii];
}



__END__

=head1 NAME

collect_vepout_coding.pl -- Collect and reformat VEP annotation for CNVs.

=head1 DESCRIPTION

This is the modified version of VEP collection scripts used for CNVs. We have changed default and extra 
columns to those most relevant for interpreting CNVs.

We also the following strategies to reformat the output.

1. GeneIDs, gene symbols if available, transcript counts will be displayed in full by this script.
But may be further simplified in further processing.

2. We want to simplify gene level functional effects. For each gene, we have selected the most severe
functional consequence based on effect ranking. When CNV affect many genes, majority of them are likely
have identical functional effects and can be unified. But consequences in other genes may be equally
or even more important. For example, a gene whose coding sequence was partially deleted by a deletion
should have the same consequence as entire gene deletion. But we may need to inspect more closely the
exons or introns that are affected.  To do this, we first sort genes based on the frequency of their 
functional effects. In case of multiple effects, the most frequent effects appearing at the end of 
the list will be collapsed into one. In this way, we remove redundant annotations but the results can
still be parsed programmably.

3. We used the same strategy to clean up trancript level information within each gene. In addition to
transcript effects, we may also have exon/intron number and transcript biotype associated with each
transcript. They will be processed as additional transcript level information in addition to TransEff.
 
4. For transcript level information, if all transcripts for a gene has the same functional effects,
then they can be collapsed into a single gene-level annotation. And if all transcript related fields can be 
collapsed, we may omit the transcript IDs (replace by ...) for that gene to further cleanup the output.

=head1 REQUIRED ARGUMENTS

=over

=item -[-]in[put] [=] <file>

The standard VEP output file for CNVs.

=for Euclid:
	file.type: readable 

=item -[-]out[put] [=] <file>

The re-formatted output file.

=back

=head1 OPTIONS

=over

=item -[-]coding

Restrict annotations to only coding transcripts (biotype=="protein_coding" or "polymorphic_pseudogene").
Recommended.

=item -[-]full-effect

Output full functional effect. The default is to select one for each variant-transcript combination.
Under full-effect mode, the un-selected effect(s) will be shown in bracket.

=item -[-]full-trans

Output full list of transcripts. The default is to replace transcript IDs with ... if all annotated transcripts
can be collaposed into the same value.

=item -[-]cutoff [=] <threshold>

Max. rank for functional effects to be appear in the GeneEff, default: 18 (coding_sequence).
Note: all possible functional effects for that gene will still appear in TransEffs column regardless of cutoff.

=for Euclid:
	threshold.default: 18

=item -[-]include-gene [=] <genelist>

A list of gene IDs, to be included. Can be used together with --exclude-gene.

=item -[-]exclude-gene [=] <genelist>

A list of gene IDs, to be excluded. 

=for Euclid:
	genelist.type: readable

=item -[-]include-trans [=] <translist>

A list of transcript ID, to be included.

=item -[-]exclude-trans [=] <translist>

A list of transcript IDs, to be excluded.

=item -[-]no-overlap

Do not reformat overlap PC.

=item -[-]no-exon

Do not reformat exon numbers.

=back






