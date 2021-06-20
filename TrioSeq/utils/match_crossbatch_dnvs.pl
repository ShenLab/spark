use strict;
use warnings;
use List::Util qw|sum|;
use Perl6::Slurp;
use Data::Dumper;
use Getopt::Euclid;
use Utils::Parser qw|sql_query|;
use Utils::Hash qw|str2hash chk_default peek_hash|;
use Utils::File::Iter qw|iter_file|;
use FaSlice; 
use Genet::Var qw|is_cpg|;
use Genome::UCSC qw|hg_chrom|;
use Genome::UCSC::TwoBit;
use Genome::UCSC::BinKeeper;

my (%hotspots, %varinfo, %sampsex, %samppheno);

if ($ARGV{'--hotspots'}) {
	%hotspots = map { (split)[0] => 1 } slurp $ARGV{'--hotspots'};
}

my $seq;
if (defined $ARGV{'--seq'}) {
	if ($ARGV{'--seq'} =~ /\.2bit$/) {
		$seq = Genome::UCSC::TwoBit->new($ARGV{'--seq'});
	}
	else {
		$seq = FaSlice->new(file => $ARGV{'--seq'});
	}
}

my ($bk, $bk_chr);
if (defined $ARGV{'--bed'}) {
	unless(-f $ARGV{'--bed'}) {
		die "Cannot find bedfile $ARGV{'--bed'}";
	}
	$bk = Genome::UCSC::BinKeeper->new($ARGV{'--bed'}, { bed => 1 });
	my $chrom = peek_hash($bk);
	$bk_chr = $chrom =~ /^chr/ ? 1 : 0;
}


my $vars1 = slurp_vars($ARGV{'--b1'}, $ARGV{'--b1-alias'}, $ARGV{'--b1-remove'}, $ARGV{'--b1-filter'});

my $vars2;
unless (defined $ARGV{'--b2'}) {
	print STDERR "We will examine within batch overlaps\n";
	$vars2 = $vars1;
}
else {
	if ($ARGV{'--b1'} eq $ARGV{'--b2'}) {
		print STDERR "b1 and b2 are the same: we will examine within batch overlaps\n";	
	}
	else {
		$vars2 = slurp_vars($ARGV{'--b2'}, $ARGV{'--b2-alias'}, $ARGV{'--b2-remove'}, $ARGV{'--b2-filter'}); 
	}
}

# IID_1, IID_2, NVars_1, NVars_2, N_SharedVars, AdjN_SharedVars SharedVars
my $fout;
if (defined $ARGV{'--output'}) {
	open $fout, ">$ARGV{'--output'}" or die "Cannot write to $ARGV{'--output'}";
}
else {
	$fout = \*STDOUT;
}

my $suf_1 = $ARGV{'--t1'} // "1";
my $suf_2 = $ARGV{'--t2'} // "2";
print $fout join("\t", "FamID_$suf_1", "IID_$suf_1", "Sex_$suf_1", "Pheno_$suf_1", "NVars_$suf_1",
					   "FamID_$suf_2", "IID_$suf_2", "Sex_$suf_2", "Pheno_$suf_2", "NVars_$suf_2",
					   "N_SharedVars",  "Nadj_SharedVars", "SharedVars"), "\n";
my %known;
foreach my $iid1 (sort keys %$vars1) {
	foreach my $iid2 (sort keys %$vars2) {
		if (!defined $ARGV{'--b2'} || $ARGV{'--b1'} eq $ARGV{'--b2'}) {
			next if $iid1 eq $iid2;
			next if $known{$iid2,$iid1};
		}
		my @sharedvars = shared_vars($vars1->{$iid1}, $vars2->{$iid2});
		my $nvars1 = scalar(keys %{$vars1->{$iid1}});
		my $nvars2 = scalar(keys %{$vars2->{$iid2}});
		my $nadj = sum(map { $hotspots{$_} ? $ARGV{'--weight'} : 1 } @sharedvars);
 		if (@sharedvars >= $ARGV{'--cutoff'}) {
 			my $varlist;
 			if (@sharedvars <= $ARGV{'--maxlist'}) {
 				if (defined $ARGV{'--details'}) {
 					$varlist = join(" ", map { defined $hotspots{$_} ? "[$_($varinfo{$_})]" : "$_($varinfo{$_})" } @sharedvars);
 				}
 				else {
 					$varlist = join(" ", map { defined $hotspots{$_} ? "[$_]" : $_ } @sharedvars);
 				}
 			}
 			print $fout join("\t", $iid1, $sampsex{$ARGV{'--b1'}}{$iid1} // ".",  $samppheno{$ARGV{'--b1'}}{$iid1} // ".", $nvars1, $iid2, 
 							defined $ARGV{'--b2'} ? $sampsex{$ARGV{'--b2'}}{$iid2} // "." : $sampsex{$ARGV{'--b1'}}{$iid2} // ".",
 							defined $ARGV{'--b2'} ? $samppheno{$ARGV{'--b2'}}{$iid2} // "." : $samppheno{$ARGV{'--b1'}}{$iid2} // ".",
 							$nvars2, scalar(@sharedvars), $nadj, $varlist // "."), "\n";
		}
		else {
			# If IDs are the same in two batches, we can also optinally output
			if ($ARGV{'--idmatch'} && $iid1 eq $iid2) {
				print $fout join("\t", $iid1, $sampsex{$ARGV{'--b1'}}{$iid1} // ".", $samppheno{$ARGV{'--b1'}}{$iid1} // ".", $nvars1, $iid2, 
								defined $ARGV{'--b2'} ? $sampsex{$ARGV{'--b2'}}{$iid2} // "." : $sampsex{$ARGV{'--b1'}}{$iid2} // ".",
								defined $ARGV{'--b2'} ? $samppheno{$ARGV{'--b2'}}{$iid2} // "." : $samppheno{$ARGV{'--b1'}}{$iid2} // ".",
								$nvars2, scalar(@sharedvars), $nadj // 0, "."), "\n";
			}
		}
		$known{$iid1,$iid2} = 1;
	}
}


sub slurp_vars {
	my ($file, $fields, $remove, $filter) = @_;
	my $alias;
	if (defined $fields) {
		$alias = str2hash($fields, {psep => ',', kvsep => ':'});
	}
	my %remove;
	if (defined $remove) {
		%remove = map { (split)[0] => 1 } slurp $remove;
	}
	my ($it, $fnames) = iter_file($file eq '-' ? \*STDIN : $file, { fsep => qr/\t/, alias => $alias });
	foreach my $field (qw|IID Chrom Position Ref Alt|) {
		unless(grep { $field eq $_ } @$fnames) {
			die "Cannot find standard field $field from $file";
		}
	}
	unless(grep { $_ eq 'Context'} @$fnames) {
		unless(defined $seq) {
			warn "Cannot determine sequence context for variants in input $file";
		}
	}
	my @fdet;
	if ($ARGV{'--details'}) {
		foreach my $field (split(',', $ARGV{'--details'})) {
			unless(grep { $field eq $_ } @$fnames) {
				warn "Cannot find variant info field $field from input $file";
			}
			else {
				push @fdet, $field;
			}
		}
	}
	my $fid;
	if (grep { $_ eq 'FamID' } @$fnames) {
		$fid = 'FamID';
	}
	else {
		warn "Use IID for FamID for $file";
		$fid = 'IID';
	}
	my $callback;
	if (defined $filter) {
		($callback, my $tokens) = sql_query($filter, 1);
		foreach my $tok (@$tokens) {
			if ($tok->[0] eq 'FIELD') {
				unless(grep { $tok->[1] eq $_ } @$fnames) {
					die "Cannot find Filter field $tok->[1] in variant table $file";
				}
			}
		}
	}

	my (%vars, %recur, %sex, %pheno, $tab_chr);
	while(my $dat = $it->()) {
		unless(defined $tab_chr) {
			$tab_chr = $dat->{Chrom} =~ /^chr/ ? 1 : 0;
		}
		if (defined $filter) {
			next unless $callback->($dat);
		}
		next if defined $remove{$dat->{IID}};
		if (defined $bk) {
			my $chrom = $dat->{Chrom};
			if ($bk_chr == 0 && $tab_chr == 1) {
				$chrom =~ s/^chr//; $chrom = 'MT' if $chrom eq 'M';
			}
			elsif ($bk_chr == 1 && $tab_chr == 0) {
				$chrom = hg_chrom($chrom);
			}
			unless ($bk->any_overlap($chrom, $dat->{Position}, $dat->{Position}+length($dat->{Ref})-1)) {
				next;
			}
		}

		my $varid = join(":", @{$dat}{qw|Chrom Position Ref Alt|});
		$varid =~ s/^chr//;
		my $fiid = $dat->{$fid}."\t".$dat->{IID};

		if (defined $dat->{Sex}) {
			chk_default(\%sex, $fiid, $dat->{Sex});
		}
		if (defined $dat->{Pheno}) {
			chk_default(\%pheno, $fiid, $dat->{Pheno});
		}
		$vars{$fiid}{$varid} = 1; 
		$recur{$varid}{$dat->{$fid}} = 1;
		if (@fdet > 0 && (!defined $varinfo{$varid} || $varinfo{$varid} =~ /^\.(\,\.)*$/)) {
			$varinfo{$varid} = join(',', @{$dat}{@fdet});
		}
		next if defined $hotspots{$varid};
		if (defined $dat->{Context}) {
			if ($dat->{Context} =~ /^([ACGT])CG>\1TG$/ || 
				$dat->{Context} =~ /^CG([ACGT])>CA\1$/) {
				$hotspots{$varid} = 1;
			}
		}
		elsif (defined $seq) {
			if ($dat->{Ref} =~ /^[ACGT]$/ && $dat->{Alt} =~ /^[ACGT]$/ && 
				is_cpg({ CHROM => $dat->{Chrom}, POS => $dat->{Position},
						 REF => $dat->{Ref}, ALT => $dat->{Alt} }, $seq)) {
				$hotspots{$varid} = 1;
			}
		}
	}
	# Adding recurrent DNVs to hotspots
	if (defined $ARGV{'--b2'} && $ARGV{'--b1'} ne $ARGV{'--b2'}) {
		foreach my $varid (sort keys %recur) {
			if (scalar(keys %{$recur{$varid}}) > 1) {
				$hotspots{$varid} = 1;
			}
		}
	}
	$sampsex{$file} = \%sex;
	$samppheno{$file} = \%pheno;
	return \%vars;
}

sub shared_vars {
	my ($muts1, $muts2) = @_;
	my @allmuts = sort keys %$muts1;
	push @allmuts => sort keys %$muts2;
	my %varct;
	foreach my $mutid (@allmuts) {
		$varct{$mutid} ++;
	}
	my @sharedmuts;
	foreach my $mutid (sort keys %varct) {
		if ($varct{$mutid} > 1) {
			push @sharedmuts, $mutid;
		}
	}
	return @sharedmuts;
}



__END__

=head1 NAME 

match_crossbatch_dnvs.pl -- Calculate the fraction of overlapping variants between samples across cohorts.

=head1 NOTES

We can use matched de novo variants to find likely overlapping samples between different cohorts.
For this purpose, we need to account for DNVs that are likely match by chance hotpots mutations).
Currently, we use the following three types of hotspots:
1. Recurrent DNVs between families within each batch (if FamID is avialable, or between individuals)
2. DNVs within CpG context (use Context if available, otherwise to query neighboring sequence)
3. Custom list of hotspot mutations (In VarID format)
Sharing a hotspot DNV will be given a customized weight (default 0.5)

Note: Because we take recurrent DNVs as hotspots, within-batch duplicated families must be removed
before this analysis. 

When the same file is provided for two batches, we will look for within-cohort duplicated samples. 
Under this setting, recurrent DNVs within cohort will not be counted as hotspots.

=head1 REQUIRED ARGUMENTS

=over

=item -[-]b1 [=] <file>

DNV table for batch 1

=back

=head1 OPTIONS

=over

=item -[-]b2 [=] <file>

DNV table for batch 2. 

=item -[-]out[put] [=] <outfile>

Output file name.

=item -[-]seq [=] <file>

Sequence file, will be used to look for sequence context when Context is not available in variant table. 

=item -[-]hotspots [=] <list>

An additional list of hotspots variants.

=item -[-]weight [=] <score>

Customized weight for hotspots variants, default: 0.5.

=for Euclid:
	score.default: 0.5

=item -[-]bed [=] <bedfile>

A bed file of genomic regions. Only variants in this region will be used to test for overlaps.

=item -[-]t1 [=] <string>

Tag for Batch 1, will be used as suffix to IID/FID (default: "_1").

=item -[-]t2 [=] <string>

Tag for Batch 2, will be used as suffix to IID/FID (default: "_2").

=item -[-]b1-alias [=] <rename>

Rename fields for batch 1 DNV table to standard names from anno_seqvars's output.
The following fields are used (optional fields are in square bracket):
IID,[Sex],[FamID],Chrom,Position,Ref,Alt,[Context],[HGNC],[GeneEff] 

=item -[-]b2-alias [=] <rename>

Key fields for batch 2 DNV table, must be renamed to standard names.

=item -[-]b1-filter [=] <expr>

Filtering expression to Batch 1 variant table.

=item -[-]b2-filter [=] <expr>

Filtering expression to Batch 2 variant table.

=item -[-]b1-remove [=] <list>

A list of sample to be removed from Batch 1.

=item -[-]b2-remove [=] <list>

A lisf of sample to be removed from Batch 2.

=item -[-]details [=] <fields>

Additional info fields for shared variants.

=item -[-]id[-]match

Also output sample pair with 0 shared var if they have the same ID.

=item -[-]cutoff [=] <number>

Min. number of shared variants to output a pair.

=for Euclid:
	number.default: 1

=item -[-]max[-]list [=] <number>

Max. number of variants to be listed in details column.

=for Euclid:
	number.default: 10

=cut


