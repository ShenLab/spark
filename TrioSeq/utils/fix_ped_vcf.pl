use strict;
use warnings;
use Carp;
use IO::Prompt;
use Graph::Undirected;
use Data::Dumper;
use File::Which;
use File::Copy qw|move|;
use Storable qw|dclone|;
use Utils::Hash qw|peek_hash array2hash|;
use Utils::List qw|all_pairs|;
use List::MoreUtils qw|any uniq all natatime|;
use Utils::File::Iter qw|slurp_file|;
use Genet::Ped;
use Genet::File::VCF qw|slurp_vcf_header|;
use Getopt::Euclid;


my %samps = map { $_->{FID}.$;.$_->{IID} => $_ } 
	slurp_file($ARGV{'--ped'}, { header => 0, fields => [qw(FID IID DAD MOM SEX AFF)] });

if ($ARGV{'--id-delim'}) {
	$ARGV{'--famid'} = $ARGV{'--id-delim'};
}

# Initial QC, basic ped file check and fix after slurp
my %FIDs;
unless ($ARGV{'--famid'}) {
	foreach my $samp (values %samps) {
		croak "Check PED input: Sample $samp->{IID} has already been observed!" if defined $FIDs{$samp->{IID}};
		$FIDs{$samp->{IID}} = $samp->{FID};
	}
}

my $SUFF;
unless ($ARGV{'--remove-dup'}) {
	$SUFF = $ARGV{'--rename-dup'};
}
my $SEP;
if ($ARGV{'--famid'}) {
	$SEP = $ARGV{'--famid'};
}

my %stdsuffs = (badsamp => 'badsamps.txt', swap => 'swapped.txt', sex => 'sex_update.txt',
	nonpar => 'nonpar.txt', dup => 'dups.txt', cryprel => 'rels.txt');

if ($ARGV{'--patch'}) {
	print STDERR "Option --patch will look for files with default suffix and override other options\n";
	while(my ($opt, $suf) = each %stdsuffs) {
		if (-f "$ARGV{'--patch'}.$suf") {
			print STDERR "Found $ARGV{'--patch'}.$suf\n";
			$ARGV{"--$opt"} = "$ARGV{'--patch'}.$suf";
		}
		else {
			print STDERR ".$suf is not found\n";
		}
	}
	prompt("Please confirm files are correctly specified (y/n): ", '-yn');
	if ($_ ne 'y') {
		print STDERR "exit\n";
		exit 1;
	}
}

my %PARSEX = (DAD => 1, MOM => 2);
foreach my $samp (values %samps) {
	foreach my $par (qw|DAD MOM|) {
		if ($samp->{$par} ne '0') {
			my $parid = $samp->{FID}.$;.$samp->{$par};
			my $parent = $samps{$parid};
			if (!defined $parent) {
				carp "Check PED input: In $samp->{FID}, $samp->{IID}'s $par was not observed";
				$samp->{$par} = '0';
			}
			else {
				# parents must have correct gender!
				if ($parent->{SEX} ne $PARSEX{$par}) {
					carp "Check PED input: In $samp->{FID}, $samp\'s $par has incorrect sex $parent->{SEX}";
					$parent->{SEX} = $PARSEX{$par};
				}
			}
		}
	}
}

# First pass, exclude bad samples, and zero out exclude parents
my %remove;
if (defined $ARGV{'--badsamp'}) {
	open my $fin, $ARGV{'--badsamp'} or croak "Cannot open sample exclude list";
	while(<$fin>) {
		next if /^\s*$/ || /^#/;
		my ($fid, $iid);
		if ($ARGV{'--famid'}) {
			my @row = split;
			croak "Badsamp: incorrect number of columns: $_" unless @row == 2;
			($fid, $iid) =  @row #(split)[0,1];
		}
		else {
			my @row = split;
			croak "Badsamp: incorrect number of columns: $_" unless @row == 1;
			$iid = $row[0]; #(split)[0];
			$fid = $FIDs{$iid} // do {croak "Cannot find FID for Sample $iid"};
		}
		if (!defined $samps{$fid,$iid}) {
			carp "Bad sample: Sample $iid of Family $fid in the exclusion list is not found in ped file";
		}
		else {
			#delete $samps{$fid,$iid};
			$remove{$fid,$iid} = 1;
		}
	}
}

# Second, fetch all SEX updates
my %sexupdate;
if (defined $ARGV{'--sex'}) {
	open my $fin, $ARGV{'--sex'} or croak "Cannot open sex update list";
	while(<$fin>) {
		next if /^\s*$/ || /^#/;
		my ($fid, $iid, $sex);
		if ($ARGV{'--famid'}) {
			my @row = split;
			croak "Sex update: incorrect number of columns: $_" unless @row == 3;
			($fid, $iid, $sex) = @row; # (split)[0,1,2];
		}
		else {
			my @row = split;
			croak "Sex update: incorrect number of columns: $_" unless @row == 2;
			($iid, $sex) =  @row;#(split)[0,1];
			$fid = $FIDs{$iid} // do {croak "Cannot find FID for Sample $iid"};
		}
		if (uc($sex) eq 'M' || uc($sex) eq 'MALE') {
			$sex = 1;
		}
		elsif (uc($sex) eq 'F' || uc($sex) eq 'FEMALE') {
			$sex = 2;
		}
		croak "Incorrect sex code update : $sex" unless $sex eq '1' || $sex eq '2';
		if (!defined $samps{$fid,$iid}) {
			carp "Sex update: Sample $iid of Family $fid in the sex update list is not found in ped file";
		}
		else {
			if ($samps{$fid,$iid}{SEX} eq $sex) {
				carp "Sex update: The sex of sample $iid of Family $fid does not change!";
			}
			else {
				$sexupdate{$fid,$iid} = $sex;
			}
		}
		if (defined $remove{$fid,$iid}) {
			carp "Sex update: sample $iid of Family $fid will be removed as bad sample";
		}
	}
	foreach my $samp (grep { $_->{SEX} !~ /^[12]$/ } values %samps) {
		unless ($sexupdate{$samp->{FID},$samp->{IID}} &&
			 	! defined $remove{$samp->{FID},$samp->{IID}} ) {
			carp "Sex update: Did not find update for $samp->{FID} - $samp->{IID} with missing gender";
		}
	}
}

# Third, fetch swapped samples
# We actually swap (rename) the samples in VCF.
my %swap;
if (defined $ARGV{'--swap'}) {
	open my $fin, $ARGV{'--swap'} or croak "Cannot open swap sample list";
	while(<$fin>) {
		next if /^\s*$/ || /^#/;
		my ($fid1, $iid1, $sex1, $fid2, $iid2, $sex2);
		if ($ARGV{'--famid'}) {
			my @row = split;
			croak "Swapped samples: incorrect number of columns: $_" unless @row == 4;
			($fid1, $iid1, $fid2, $iid2) = @row; #(split)[0,1,2,3];
		}
		else {
			my @row = split;
			croak "Swapped samples: incorrect number of columns: $_" unless @row == 2;
			($iid1, $iid2) = @row; # (split)[0,1];
			$fid1 = $FIDs{$iid1} // do {croak "Cannot find FID for Sample $iid1"};
			$fid2 = $FIDs{$iid2} // do {croak "Cannot find FID for Sample $iid2"};
		}
		unless (defined $samps{$fid1,$iid1} && defined $samps{$fid2,$iid2}) {
			carp "Swapped samples: Sample $iid1 in Family $fid1 or $iid2 in Family $fid2 is not found in ped file";
			next;
		}
		if ($fid1.$;.$iid1 eq $fid2.$;.$iid2) {
			carp "Swapped samples: Swapped sample have identical IDs: $fid1-$iid1"; 
			next;
		}
		$sex1 = $samps{$fid1,$iid1}{SEX};
		if ($sex1 !~ /^[12]$/) {
			$sex1 = $sexupdate{$fid1,$iid1};
		}
		$sex2 = $samps{$fid2,$iid2}{SEX};
		if ($sex2 !~ /^[12]$/) {
			$sex2 = $sexupdate{$fid2,$iid2};
		}
		$swap{$fid1,$iid1} = $fid2.$;.$iid2;
		$swap{$fid2,$iid2} = $fid1.$;.$iid1;

		if (defined $remove{$fid1,$iid1} && defined $remove{$fid2,$iid2}) {
			carp "Swapped samples: Both of the swap pair $fid1 - $iid1 <=> $fid2 - $iid2 will be removed";
		}
		elsif (defined $remove{$fid1,$iid1} && ! defined $remove{$fid2,$iid2} ||
			! defined $remove{$fid1,$iid1} && defined $remove{$fid2,$iid2}) {
			carp "Swapped samples: One of the swap pair $fid1 - $iid1 <=> $fid2 - $iid2 will be removed";
			# NOTE: when processing PED and VCF file, we should remove sample first, then rename based on swap
		}

		# If the swapped samples has sex update, then check if updated sex is correct
		if ($sex1 ne $sex2) {
			foreach my $fiids ([$fid1, $iid1, $sex2], [$fid2, $iid2, $sex1]) {
				my ($fid, $iid, $sex) = @$fiids;
				if (defined $sexupdate{$fid,$iid}) {
					# if sex update is defined
					if ($sexupdate{$fid,$iid} eq $sex)  {
						delete $sexupdate{$fid,$iid};
					}
					else {
						croak "Swapped samples: Sex update for $iid in Family $fid ($sexupdate{$fid,$iid}) is incorrect, should be $sex!";
					}
				}
				else {
					croak "Swapped samples: The swap pair $fid1 - $iid1 <=> $fid2 - $iid2 have different sex, but not found in sex update."
				}
			}
		}
	}
}


# Fourth, rename error duplicates
my %dups; # hold ped info of primary duplicated sample 
my %rename; # in case the sample to be kept is different from primary duplicated sample, rename kept samp
my %duplcats; # indicate if a samp is dup of another
my %sampdups; 

# Duplicated samples shuold be automatically added to relpairs unless they are removed.
my %samprels;
if ($ARGV{'--dup'}) {
	my $nn = $ARGV{'--famid'} ? 2 : 1;
	open my $fin, $ARGV{'--dup'} or croak "Cannot open duplication list";
	while(<$fin>) {
		next if /^\s*$/ || /^#/;
		my @row = split;
		my $it = natatime $nn, @row;
		my (@dups, $keep, $flag);
		while(my @id = $it->()) {
			croak "Duplication: Incorrect number columns in duplication list"  unless scalar(@id) == $nn;
			if (any { /^\[\S+\]$/ } @id) {
				foreach (@id) {
					s/^\[//; s/\]$//;
				}
				$flag = 1;
			}
			else {
				$flag = 0;
			}
			if (@id == 1) {
				my $fid = $FIDs{$id[0]} // do { croak "Cannot find FID for Sample $id[0]" };
				unshift @id, $fid;
			}
			my $fiid = join($;, @id);
			if (!defined $samps{$fiid}) {
				carp "Duplication: Cannot find $id[1] in Family $id[0] in original PED!";
				next;
			}
			$keep = $fiid if $flag;
			push @dups, $fiid;
		}
		unless (scalar(uniq sort @dups) > 1) {
			carp "Duplication: Incorrect number of duplicates";
			next;
		}
		croak "Duplication: Cannot find samples to be kept from duplicates" unless defined $keep;
		croak "Duplication: Different duplication list overlap with each other, put them in the same line" if (any { $dups{$_} } @dups);
		# Check consistency with sex update.
		# Note if one of dup samp is part of swapped pair, sexupdate has been deleted.
		my @dupsex = map { defined $swap{$_} ? $samps{$swap{$_}}{SEX} : $sexupdate{$_} // $samps{$_}{SEX}  } @dups;
	 	unless (scalar(uniq sort @dupsex) == 1) {
	 		print STDERR join(",", map { join '-', split($;, $_)  } @dups), "\t", join(',', @dupsex), "\n";
	 		croak "Duplication: The duplicated sample have inconsistent sex";
	 	}

	 	# Check if primary dup'ed sample (appeared at beginning of dup line) is swapped
	 	# We only allow primary dup'ed sample to appear in swapped list
		my $primary;
		if ($swap{$dups[0]}) {
			$primary = $swap{$dups[0]};
			carp "Duplication: The primary sample of duplication @{[ join '-', split($;, $dups[0]) ]} ".
				"was swapped with @{[ join '-', split($;, $primary) ]}";
		}
		else {
			$primary = $dups[0];
		}

		# Check if kept sample is swapped

		# Keep record of duplicated samples
		# Also check with swapped sample.
		for(my $ii = 0; $ii < @dups; $ii ++) {
			$dups{$dups[$ii]} = $samps{$primary};
			if ($ii > 0) {
				$duplcats{$dups[$ii]} = 1;
				if ($swap{$dups[$ii]}) {
					# Check if the duplication list has overlap with swapped pairs
					# This is not allowed
					croak "Duplication: The duplicate $ii of @{[ join '-', split($;, $primary) ]} was also found in the swap list.";
				}
			}
		}

		# $keep is the best quality duplicated sample, and will be considered to be kept.
		if (defined $remove{$keep}) {
			my ($fid, $iid) = split($;, $primary);
			# In this case, all duplicated samples should be removed.
			if (all { $remove{$_} } @dups) {
				carp "Duplication: Sample $iid in Family $iid and all duplicates are removed";
			}
			else {
				# The kept sample should not be in removal list
				croak "Duplication: Sample $iid in Family $fid was suppose to be kept but found in the remove list";
			}
		}
		else {
			if ($keep ne $dups[0]) {
				$rename{$keep} = $primary;
			}
			if ($ARGV{'--remove-dup'}) {
				foreach my $badsamp (grep { $_ ne $keep } @dups) {
					$remove{$badsamp} = 1;
				}
			}
			else {
				# Keep all pairwise combination of duplicates in duppair and relpair
				my @othersamp = grep { $_ ne $keep  && !defined $remove{$_} } @dups;
				if (@othersamp > 0) {
					my ($fid, $iid) = split($;, $primary);
					$rename{$othersamp[0]} = $primary.$SUFF;
					push @{$sampdups{$fid}}, [$iid, $iid.$SUFF];
					# for relpairs
					my @iids = ($iid, $iid.$SUFF);
					if (@othersamp > 1) {
						for(my $ii = 1; $ii < @othersamp; $ii ++) {
							$rename{$othersamp[$ii]} = $primary.$SUFF.$ii;
							push @{$sampdups{$fid}}, [$iid, $iid.$SUFF.$ii];
							push @iids, $iid.$SUFF.$ii;
						}
					}
					push @{$samprels{$fid}}, all_pairs(@iids);
				}
			}
		}
	}
}

# Non-parantality
my %nonpar;
if ($ARGV{'--nonpar'}) {
	open my $fin, $ARGV{'--nonpar'} or croak "Cannot open duplication list";
	while(<$fin>) {
		next if /^\s*$/ || /^#/;
		my ($fid, $iid, $par, $fid2, $iid2);
		if ($ARGV{'--famid'}) {
			my @row = split;
			croak "Non-parantality: incorrect number of columns: $_" unless @row == 3 || @row == 5;
			($fid, $iid, $par, $fid2, $iid2) = @row; #(split)[0,1,2];
		}
		else {
			my @row = split;
			croak "Non-parantality: incorrect number of columns: $_" unless @row == 2 || @row == 3;
			($iid, $par, $iid2) = @row; #(split)[0,1];
			$fid = $FIDs{$iid} // do { croak "Cannot find FID for Sample $iid" };
			if (defined $iid2) {
				$fid2 = $FIDs{$iid2} // do { croak "Cannot find FID for Sample $iid2" };	
			}
		}
		croak "Non-parantality: Incorrect parent $par" unless $par eq 'DAD' || $par eq 'MOM';
		unless (defined $samps{$fid,$iid}) {
			carp "Non-parantality: Cannot find $iid in Family $fid";
		}
		else {
			# Do not consider non-parentality caused due to sample swap or duplication
			if (defined $swap{$fid,$iid} || defined $duplcats{$fid,$iid}) {
				carp "Non-parantality: $fid - $iid already appear in sample swap or duplicates list ".
					"and should be excluded from nonpar list";
			}
			else {
				if (!defined $fid2) {
					if ($samps{$fid,$iid}{$par} ne '0') {
						#push @{$nonpar{$fid,$iid}}, $par;
						$nonpar{$fid,$iid}{$par} = '0';
					}
					else {
						carp "Non-parantality: $iid in Family $fid has no $par in original PED";
					}
				}
				else {
					if ($fid2 eq $fid) {
						croak "Non-parantality: Cannot find new $par $iid2 in Family $fid" unless defined $samps{$fid,$iid2};
						croak "Non-parantality: New $par $fid - $iid2 alredy appear in sample duplicates list" if defined $duplcats{$fid,$iid2};
						$nonpar{$fid,$iid}{$par} = $iid2;
					}
					else {
						croak "Non-parantality: currently do not support new parents from different family";
					}
				}
			}
		}
	}
}

# Fix PED
# First update SEX and parents, then remove bad samples, then rename duplicated samples
# finally fix missing parents issues.
foreach my $fiid (keys %sexupdate) {
	#next unless $samps{$fiid};
	$samps{$fiid}{SEX} = $sexupdate{$fiid};
}
while(my ($fiid, $pars) = each %nonpar) {
	#next unless $samps{$fiid};
	while(my ($pa, $id2) = each %$pars) {
		$samps{$fiid}{$pa} = $id2;
	}
}
foreach my $fiid (keys %rename) {
	$samps{$fiid} = dclone $dups{$fiid};
	$samps{$fiid}{IID} = (split($;, $rename{$fiid}))[1];
}
foreach my $fiid (keys %remove) {
	delete $samps{$fiid};
}

# Then we consider cryptic relatedness that not due to sample swap or err duplication
#my (%fidrename, %samprels);
my %fidrename;
if ($ARGV{'--cryprel'}) {
	# if there is a between family relatedness, then FIDs will be changed
	my $famrel = Graph::Undirected->new();
	my @relpairs;
	my $nn = $ARGV{'--famid'} ? 2 : 1;
	open my $fin, $ARGV{'--cryprel'} or croak "Cannot open relatedness list";
	while(<$fin>) {
		next if /^\s*$/ || /^#/;
		my @row = split;
		my $it = natatime $nn, @row;
		my (@fids, @fiids);
		while(my @id = $it->()) {
			croak "Relatedness: Incorrect number of columns in duplication list" unless scalar(@id) == $nn;
			if (@id == 1) {
				my $fid = $FIDs{$id[0]} // do { croak "Cannot find FID for Sample $id[0]" };
				unshift @id, $fid;
			}
			my $fiid = join($;, @id);
			if (!defined $samps{$fiid}) {
				carp "Relatedness: Cannot find $id[1] in Family $id[0] in original PED";
				next;
			}
			if ($swap{$fiid}) {
				@id = split($;, $swap{$fiid});
				$fiid = join($;, @id);
				carp "Relatedness: Sample @{[ join ' - ', split($;, $fiid) ]} is swapped with ".
					"@{[ join('-', @id) ]}";
			}
			elsif ($duplcats{$fiid}) {
				carp "Relatedness: Sample @{[ join ' - ', split($;, $fiid) ]} is part of duplicated sample, ignore";
				next;
			}
			push @fids, $id[0];
			push @fiids, $fiid;
		}
		unless (scalar(uniq sort @fiids) > 1) {
			carp "Relatedness: Incorrect number of related samples" ;
			next;
		}
		my @uqfids = uniq sort @fids;
		# related pairs from different families
		if (@uqfids > 1) {
			foreach my $fpair (all_pairs(@uqfids)) {
				$famrel->add_edge(@$fpair);
			}
		}
		push @relpairs, [@fiids];
	}
	# Now update FIDs for merged families
	foreach my $cc ($famrel->connected_components) {
		my $newfid = join($ARGV{'--merge-fam'}, sort @$cc);
		print STDERR "Merged family: $newfid\n";
		# Check this newfid has not been observed before
		if (grep { $newfid eq $_ } values %FIDs) {
			croak "New FID $newfid exist in original PED";
		}
		# Check samples with newfid has unique names
		my %cc = array2hash(@$cc);
		my @iids = map { $_->{IID} } grep { defined $cc{$_->{FID}} } values %samps;
		my @uqiids = uniq sort @iids;
		if (@iids != @uqiids) {
			croak "Samples in merged family $newfid have non-unique IDs";
		}
		foreach (@$cc) {
			$fidrename{$_} = $newfid;
		}
	}
	# Now rename FIDs in samps
	foreach my $indiv (values %samps) {
		my $newfid = $fidrename{$indiv->{FID}};
		if ($newfid) {
			$indiv->{FID} = $newfid;
		}
	}
	# Now second pass, slurp all related samples
	foreach my $fiids (@relpairs) {
		my (@fids, @iids);
		foreach my $fiid (@$fiids) {
			my @id = split $;, $fiid;
			if ($fidrename{$id[0]}) {
				push @fids, $fidrename{$id[0]};
			}
			else {
				push @fids, $id[0];
			}
			push @iids, $id[1];
		}
		unless (scalar(uniq sort @fids) == 1) {
			croak "Relatedness: In second pass, related sample should have same FID" ;
		}
		if (@iids < 2) {
			carp "Relatedness: Incorrect number of related samples samples";
			next;
		}
		push @{$samprels{$fids[0]}}, all_pairs(@iids);	
	}
	# We also add sampdups into samprels
	while(my ($fid, $dpairs) = each %sampdups) {
		if ($fidrename{$fid}) {
			$fid = $fidrename{$fid};
		}
		push @{$samprels{$fid}}, @$dpairs;
	}
}


my $ped = Genet::Ped->new([values %samps]);

my $fout  = IO::File->new("$ARGV{'--output'}.ped", "w");
# Now output PED file
# NOTE: Disconnected sub-fams that are initially merged in the same family will be split, so must provide
# cryptic relatedness pair to keep them merged.
# NOTE: Disconnected components due to additional duplication will also be kept.
my @scc_fids = grep {  $ped->get_connected($_) == 1 } $ped->get_famids();
my $newped = $ped->clone(\@scc_fids);
$newped->write_ped($fout);

my @mcc_fids = grep {  $ped->get_connected($_) > 1 } $ped->get_famids();
if ($ARGV{'--no-split'}) {
	my $newped = $ped->clone(\@mcc_fids);
	$newped->write_ped($fout);
}
else {
	foreach my $fid (@mcc_fids) {
		if (defined $samprels{$fid}) {
			my $grel = dclone $ped->{$fid};
			my @relpairs = @{$samprels{$fid}};
			foreach my $rpair (@relpairs) {
				$grel->add_edge(@$rpair);
			}
			my @cc = sort { @$b <=> @$a } $grel->weakly_connected_components();
			if (@cc > 1) {
				my $splitped = $ped->sub_peds($fid, \@cc, { sep => $ARGV{'--sub-fam'}, verbose => 1 });
				$splitped->write_ped($fout);
			}
			else {
				my $singleped = $ped->clone($fid);
				$singleped->write_ped($fout);
			}
		}
		else {
			my $splitped = $ped->sub_peds($fid, undef, { sep => $ARGV{'--sub-fam'}, verbose => 1 });
			$splitped->write_ped($fout);
		}
	} 
}


## Fix VCF
if ($ARGV{'--vcf'}) {
	my $bcftools = which("bcftools");
	croak "Cannot find bcftools" unless $bcftools;

	# Remove bad samples (REMOVE BAD SAMPLE BEFORE RENAME)
	my $vcfclean;
	if (keys %remove) {
		$vcfclean = "$ARGV{'--output'}.clean.vcf.gz";
		my @badsamps;
		foreach my $fiid (keys %remove) {
			my ($fid, $iid) = split($;, $fiid);
			if ($ARGV{'--famid'}) {
				push @badsamps, $fid.$SEP.$iid;
			}
			else {
				push @badsamps, $iid;
			}
		}
		my $badsamps = join(',', @badsamps);
		print STDERR "Run bcftools to remove bad samples ...\n";
		if (@badsamps < 10) {
			system(qq|bcftools view -s ^$badsamps -a -c 1 -O z -o $vcfclean --threads $ARGV{'--threads'} $ARGV{'--vcf'}|);
		}
		else {
			my $tmpdir = tempdir();
			open my $fout, "$tmpdir/badsamps.txt" or die "Cannot write to badsamps.txt";
			print $fout join("\n", @badsamps), "\n";
			close $fout;
			system(qq|bcftools view -S ^$tmpdir/badsamps.txt -a -c 1 -O z -o $vcfclean --threads $ARGV{'--threads'} $ARGV{'--vcf'}|);
		}
		system(qq{ tabix -p vcf -f $vcfclean });
		print STDERR "Done\n";
	}
	else {
		print STDERR "No bad sample will be removed from VCF.\n";
		$vcfclean = $ARGV{'--vcf'};
	}
		
	# Swap and rename
	if (keys %swap || keys %rename) {
		my ($header, $sampids) = slurp_vcf_header($vcfclean); 
		my $swap = vcf_sampids(\%swap);
		if (defined $swap) {
			foreach my $iid (@$sampids) {
				if (defined $swap->{$iid}) {
					$iid = $swap->{$iid};
				}
			}
		}
		my $rename = vcf_sampids(\%rename);
		if (defined $rename) {
			foreach my $iid (@$sampids) {
				if (defined $rename->{$iid}) {
					$iid = $rename->{$iid};
				}
			}
		}
		open my $fout, ">$ARGV{'--output'}.vcf_header" or croak "Cannot write new vcf header";
		print $fout $header, join("\t", qw|#CHROM POS ID REF ALT QUAL FILTER INFO FORMAT|, @$sampids), "\n";
		$fout->close;

		print STDERR "Run bcftools to rename samples ...\n";
		system(qq|bcftools reheader -h $ARGV{'--output'}.vcf_header -o $ARGV{'--output'}.vcf.gz $vcfclean|);
		system(qq|tabix -p vcf -f $ARGV{'--output'}.vcf.gz|);
		print STDERR "Done\n";
	}
	else {
		print STDERR "No individual needs to be renamed.\n";
		if ($vcfclean ne $ARGV{'--vcf'}) {
			move($vcfclean,  "$ARGV{'--output'}.vcf.gz");
			move("$vcfclean.tbi",  "$ARGV{'--output'}.vcf.gz.tbi");
		}
	}	
}

# Write the rename list
write_samplst("$ARGV{'--output'}.vcf_remove", \%remove);
my %fullrenm;
while(my ($fiid1, $fiid2) = each %swap) {
	if (defined $rename{$fiid2}) {
		$fullrenm{$fiid1} = $rename{$fiid2};
		delete $rename{$fiid2};
	}
	else {
		$fullrenm{$fiid1} = $fiid2;
	}
}
while(my ($fiid1, $fiid2) = each %rename) {
	$fullrenm{$fiid1} = $fiid2;
}
write_samplst("$ARGV{'--output'}.vcf_rename", \%fullrenm);


# Remove intermediate files
if (-f "$ARGV{'--output'}.clean.vcf.gz") {
	unlink("$ARGV{'--output'}.clean.vcf.gz", "$ARGV{'--output'}.clean.vcf.gz.tbi");
}

sub write_samplst {
	my ($file, $samps) = @_;
	return unless keys %$samps > 0;
	open my $fout, ">$file" or croak "Cannot write $file";
	my $firstval = (peek_hash($samps))[1];
	if (index($firstval, $;) !=  -1) {
		while(my ($fiid1, $fiid2) = each %$samps) {
			my ($fid1, $iid1) = split($;, $fiid1);
			my ($fid2, $iid2) = split($;, $fiid2);
			if ($ARGV{'--famid'}) {
				print $fout $fid1.$SEP.$iid1, "\t", $fid2.$SEP.$iid2, "\n";
			}
			else {
				print $fout $iid1, "\t", $iid2, "\n";
			}
		}
	}
	else {
		foreach my $fiid (keys %$samps) {
			my ($fid, $iid) = split($;, $fiid);
			if ($ARGV{'--famid'}) {
				if ($samps->{$fiid} eq '1') {
					print $fout $fid.$SEP.$iid, "\n";
				}
				else {
					print $fout $fid.$SEP.$iid, "\t", $samps->{$fiid}, "\n";
				}	
			}
			else {
				if ($samps->{$fiid} eq '1') {
					print $fout $iid, "\n";
				}
				else {
					print $fout $iid, "\t", $samps->{$fiid}, "\n";
				}
			}
		}
	}
}

sub vcf_sampids {
	my ($samps) = @_;
	my %newsamps;
	return unless keys %$samps > 0;
	my $firstval = (peek_hash($samps))[1];
	if (index($firstval, $;) !=  -1) {
		while(my ($fiid1, $fiid2) = each %$samps) {
			my ($fid1, $iid1) = split($;, $fiid1);
			my ($fid2, $iid2) = split($;, $fiid2);
			if ($ARGV{'--famid'}) {
				$newsamps{$fid1.$SEP.$iid1} = $fid2.$SEP.$iid2;
			}
			else {
				$newsamps{$iid1} = $iid2;
			}
		}
	}
	else {
		foreach my $fiid (keys %$samps) {
			my ($fid, $iid) = split($;, $fiid);
			if ($ARGV{'--famid'}) {
				my $iid2 = $samps->{$fiid};
				my $fid2 = $FIDs{$iid2};
				$newsamps{$fid.$SEP.$iid} = $fid2.$SEP.$iid2;
			}
			else {
				$newsamps{$iid} = $samps->{$fiid};
			}
		}		
	}
	if (wantarray) {
		return %newsamps;
	}
	else {
		return \%newsamps;
	}
}


__END__

=head1 NAME

fix_ped_vcf -- Fix gender and relatedness issues in PED and (optionally) VCF.

=head1 DESCRIOTION

Gender and relatedness errors are two major sample quality issues in trio-based
sequencing studies. After sample level QC with C<samp_qc>, we get a list of samples
with gender errors and a pairwise kinship coefficient estimates. Then we can use
that information to infer possible errors like sample swap, error duplication, etc.
Given that, we can fix the PED and VCF file.

The input of this script include following four (optional) files:
[In all files below, FID is optinal]
1. Bad Samples. Samples considered to be removed because of low sequencing depth, 
   excessive missingness, etc.
   Format: FID IID
2. Gender Updates.
   Format: FID IID GENDER
3. Sample Swap. Pairs of swapped samples.
   Format: FID1 IID1 FID2 IID2
4. Duplications.
   Format: FID1 IID1 [FID2] [IID2] FID3 IID3... 
   The target sample that was duplicated should appear first, the sample that should 
   be kept as the main data is shown in square bracket.
5. Error Parent-offspring relationships (after accounting for swap and duplication).
   Format: FID1 IID1 DAD/MOM (FID2 IID2)
   Pairs of samples that did not show expected PO relationship as indicated in the 
   original ped. When FID2/IID2 is provided, it will be used as real parent.
6. Cryptic relatedness.
   Format: FID1 IID1 [FID2] [IID2] FID3 IID3... 
   All samples within the same line are related with each other, but do not reflected
   in the original PED file.

The four list above will be internally checked for consistency. For example,
if a pair of parents were swapped, they should also appear in gender update list.
Another example, if one of the swapped pair was erroneously duplicated, then
the one appear in the swapped pair must be the duplication target (primary) sample, 
otherwise, the sample in swapped pair would appear more than once. Further example, 
all samples in the duplication list must have the same SEX.

Processing order for PED: remove badsampes, update sex/nonpar, rename/remove duplicates,
finally patch relatedness

Processing order for VCF: remove badsampes and duplicates (if any), rename swapped samples, 
rename duplicates if not removed


=head1 REQUIRED ARGUMENTS

=over

=item -[-]ped [=] <file>

Input ped file.

=for Euclid:
  file.type: readable

=item -[-]out[put] [=] <prefix>

Output file prefix.

=back

=head1 OPTIONS

=over

=item -[-]vcf [=] <file>

Input VCF file.

=for Euclid:
  file.type: readable

=item -[-]patch [=] <prefix>

Instead of specifying patched files individually, they can be provided by a prefix
and using standard suffices as follows.

--badsamp  badsamps.txt
--swap      swapped.txt
--sex       sex_update.txt
--nonpar    nonpar.txt
--dup       dups.txt
--cryprel   rels.txt

=item -[-]bad[samp] [=] <badsamp>

A list of bad samples that should be excluded.

=item -[-]swap [=] <pairlist>

A list of samples that are swapped.
Format: FID1 IID1 FID2 IID2 "Comments"

=item -[-]sex [=] <sexupdate>

A list of sex update.

=item -[-]dup [=] <duplist>

A list of duplicated sample. 

Format: FID1 IID1 [FID2] [IID2] FID3 IID3
the one that was duplicated will appear first, the replicates whose data we wish
to keep will shown in square brackets.

=item -[-]remove-dup

Remove duplicated samples with less good quality.
The default is to keep all duplicated samples, and rename samples with less good quality,
unless they were also specified in the badsamp list.

=item -[-]no[n]par [=] <nonparent>

A list of non-parentality.

=item -[-][cryp]rel [=] <relpairs>

A list of cryptic related pairs.
By default, we will use the FID in the original PED. When a family was split due to 
the removal of family members, FID will also be changed to FID.1,2.... 

=item -[-]id-delim [=] <char>

=item -[-]fam[i]id [=] <char>

By default, we will use IID found in PED file to represent individuals in VCF file.
When IIDs are not unique across families, then we should use FID-IID combination. Use this
option to specify the FID-IID delimiter. This is similar to PLINK2's --id-delim option.
When this option is enabled, all input files should use both FID and IID to specify an 
individual.

=item -[-]rename-dup[s] [=] <suffix>

When --remove-dup is not enabled, duplicated samples will be kept and renamed.
The suffix used to rename duplicated sample, default is '_Re'.

=for Euclid:
  suffix.default: "_Re"

=item -[-]no-split

Do not split families with multiple disconnected components.

=item -[-]sub-fam [=] <string>

Provide a separator used to add numerical suffix disconnected sub-families, default is ".".

=for Euclid:
  string.default: "."

=item -[-]merge-fam [=] <string>

When where are cryptic relatedness between families or between different disconnected
components of families. Fam IDs will be automatically concated, default is "-".

=for Euclid:
  string.default: "-" 

=item -[-]threads [=] <num>

Number of threads when running bcftools, default: 4.

=for Euclid:
  num.default: 4

=back

=cut





