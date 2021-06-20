use strict;
use warnings;
use IO::Prompt;
use Data::Dumper;
use Getopt::Euclid;
use Utils::File::Iter qw|slurp_file|;

my @samps = slurp_file($ARGV{'--ped'}, { header => 0, fields => [qw(FID IID DAD MOM SEX AFF)] });

# Validity check for original input PED file
validate_ped(\@samps);

my %stdsuffs = ( 'keep' => 'keep.txt', 'remove' => 'remove.txt', 
				'update-fids' => 'update_fids.txt', 'update-iids' => 'update_iids.txt',
				'update-sex' => 'update_sex.txt', 'update-pars' => 'update_pars.txt', 
				'update-pheno' => 'update_pheno.txt');

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

my %keep;
# Step 0: Keep good samples
if (defined $ARGV{'--keep'}) {
	open my $fin, $ARGV{'--keep'} or die "Cannot open sample removal list: $ARGV{'--keep'}";
	while(<$fin>) {
		next if /^\s*$/ || /^#/;
		my ($fid, $iid);
		my @row = split;
		if (@row == 2) {
			($fid, $iid) = @row;
			$keep{$fid,$iid} = 1;
		}
		elsif (@row == 1) {
			$iid = $row[0];
			$keep{$iid} = 1;
		}
		else {
			die "Keep list: Incorrect number of columns : $_";
		}	
	}
	for(my $ii = $#samps; $ii >= 0; $ii --) {
		#print Dumper $ii, $samps[$ii]; exit 1; 
		unless (defined $keep{$samps[$ii]{FID},$samps[$ii]{IID}} || defined $keep{$samps[$ii]{IID}}) {
			splice @samps, $ii, 1;
		}
	}
}

my %remove;
# Step 1: Remove bad samples
if (defined $ARGV{'--remove'}) {
	open my $fin, $ARGV{'--remove'} or die "Cannot open sample removal list: $ARGV{'--remove'}";
	while(<$fin>) {
		next if /^\s*$/ || /^#/;
		my ($fid, $iid);
		my @row = split;
		if (@row == 2) {
			($fid, $iid) = @row;
			$remove{$fid,$iid} = 1;
		}
		elsif (@row == 1) {
			$iid = $row[0];
			$remove{$iid} = 1;
		}
		else {
			die "Remove list: Incorrect number of columns : $_";
		}	
	}
	for(my $ii = $#samps; $ii >= 0; $ii --) {
		if (defined $remove{$samps[$ii]{FID},$samps[$ii]{IID}} || defined $remove{$samps[$ii]{IID}}) {
			splice @samps, $ii, 1;
		}
	}
}

# Step 2: FamID update
my %renamefam;
if (defined $ARGV{'--update-fids'}) {
	open my $fin, $ARGV{'--update-fids'} or die "Cannot open FamID update list: $ARGV{'--update-fids'}";
	while(<$fin>) {
		next if /^\s*$/ || /^#/;
		my @row = split;
		if (@row == 2) {
			my ($oldfid, $newfid) = @row;
			$renamefam{$oldfid} = $newfid;
		}
		elsif (@row == 3) {
			my ($oldfid, $oldiid, $newfid) = @row;
			$renamefam{$oldfid,$oldiid} = $newfid;
		}
		else {
			die "FamIDs update: Incorrect number of columns : $_";
		}
	}
	foreach (@samps) {
		if (defined $renamefam{$_->{FID}} || defined $renamefam{$_->{FID},$_->{IID}}) {
			$_->{FID} = $renamefam{$_->{FID}} // $renamefam{$_->{FID},$_->{IID}};
		}
	}
}

# Step 3: IID update
my %renamesamp;
if (defined $ARGV{'--update-iids'}) {
	open my $fin, $ARGV{'--update-iids'} or die "Cannot open IID update list: $ARGV{'--update-iids'}";
	while(<$fin>) {
		next if /^\s*$/ || /^#/;
		my @row = split;
		my ($fid, $oldiid, $newiid);
		if (@row == 2) {
			($oldiid, $newiid) = @row;
			$renamesamp{$oldiid} = $newiid;
		}
		elsif (@row == 3) {
			($fid, $oldiid, $newiid) = @row;
			$renamesamp{$fid,$oldiid} = $newiid;
		}
		else {
			die "IID update: Incorrect number of columns : $_";
		}
	}
	foreach (@samps) {
		if (defined $renamesamp{$_->{IID}} || defined $renamesamp{$_->{FID},$_->{IID}}) {
			$_->{IID} = $renamesamp{$_->{IID}} // $renamesamp{$_->{FID},$_->{IID}};
		}
		foreach my $PAR (qw|DAD MOM|) {
			if ($_->{$PAR} ne '0' && (defined $renamesamp{$_->{$PAR}} || defined $renamesamp{$_->{FID},$_->{$PAR}})) {
				$_->{$PAR} = $renamesamp{$_->{$PAR}} // $renamesamp{$_->{FID},$_->{$PAR}};
			}
		}
	}
}

# Step 4: Update sex
if (defined $ARGV{'--update-sex'}) {
	my %updatesex;
	open my $fin, $ARGV{'--update-sex'} or die "Cannot open sex update list: $ARGV{'--update-sex'}";
	while(<$fin>) {
		next if /^\s*$/ || /^#/;
		my @row = split;
		my ($fid, $iid, $sex);
		if (@row == 2) {
			($iid, $sex) = @row;
		}
		elsif (@row == 3) {
			($fid, $iid, $sex) = @row;
		}
		else {
			die "Sex update: Incorrect number of columns : $_";
		}
		if ($sex eq '1' || $sex =~ /^m/i) {
			if (defined $fid) {
				$updatesex{$fid,$iid} = 1;
			}
			else {
				$updatesex{$iid} = 1;
			}	
		}
		elsif ($sex eq '2' || $sex =~ /^f/i) {
			if (defined $fid) {
				$updatesex{$fid,$iid} = 2;
			}
			else {
				$updatesex{$iid} = 2
			}
		}
		else {
			warn "Cannot recognize sex code: $sex";
			if (defined $fid) {
				$updatesex{$fid,$iid} = 0;
			}
			else {
				$updatesex{$iid} = 0;
			}			
		}
	}
	foreach (@samps) {
		if (defined $updatesex{$_->{IID}} || defined $updatesex{$_->{FID},$_->{IID}}) {
			$_->{SEX} = $updatesex{$_->{IID}} // $updatesex{$_->{FID},$_->{IID}};
		}

	}
}

# Step 5: update parants
if (defined $ARGV{'--update-pars'}) {
	my %updatepars;
	open my $fin, $ARGV{'--update-pars'} or die "Cannot open parents update list: $ARGV{'--update-pars'}";
	while(<$fin>) {
		next if /^\s*$/ || /^#/;
		my @row = split;
		my ($fid, $iid, $dad, $mom);
		if (@row == 3) {
			($iid, $dad, $mom) = @row;
			$updatepars{$iid} = [$dad, $mom];
		}
		elsif (@row == 4) {
			($fid, $iid, $dad, $mom) = @row;
			$updatepars{$fid,$iid} = [$dad, $mom];
		}
		else {
			die "Parents update: Incorrect number of columns : $_";
		}
	}
	foreach (@samps) {
		if (defined $updatepars{$_->{IID}}) {
			$_->{DAD} = $updatepars{$_->{IID}}[0];
			$_->{MOM} = $updatepars{$_->{IID}}[1];
		}
		elsif (defined $updatepars{$_->{FID},$_->{IID}}) {
			$_->{DAD} = $updatepars{$_->{FID},$_->{IID}}[0];
			$_->{MOM} = $updatepars{$_->{FID},$_->{IID}}[1];
		}
	}
}

# Step 6: update phenotype
if (defined $ARGV{'--update-pheno'}) {
	my %updatepheno;
	open my $fin, $ARGV{'--update-pheno'} or die "Cannot open pheno update list: $ARGV{'--update-pheno'}";
	while(<$fin>) {
		next if /^\s*$/ || /^#/;
		my @row = split;
		my ($fid, $iid, $pheno);
		if (@row == 2) {
			($iid, $pheno) = @row;
			$updatepheno{$iid} = $pheno;
		}
		elsif (@row == 3) {
			($fid, $iid, $pheno) = @row;
			$updatepheno{$fid,$iid} = $pheno;
		}
		else {
			die "Pheno update: Incorrect number of columns : $_";
		}
	}
	foreach (@samps) {
		if (defined $updatepheno{$_->{IID}} || defined $updatepheno{$_->{FID},$_->{IID}}) {
			$_->{AFF} = $updatepheno{$_->{IID}} // $updatepheno{$_->{FID},$_->{IID}};
		}
	}
}

# Before output, validate again with stringency
my $FIDs = validate_ped(\@samps, 1);

open my $fped, ">$ARGV{'--output'}" or die "Cannot write to $ARGV{'--output'}";
foreach my $samp (sort { $a->{FID} cmp $b->{FID} || $a->{IID} cmp $b->{FID} } @samps) {
	print $fped join("\t", @{$samp}{qw|FID IID DAD MOM SEX AFF|}), "\n";
}

# Validate sample info from PED file
sub validate_ped {
	my ($samps, $strict) = @_;
	my $type = $strict ? "output" : "input";

	# First round, unique FID-IID combination, keep sex info
	my %known;
	foreach my $samp (@$samps) {
		if (defined $known{$samp->{FID},$samp->{IID}}) {
			die "Check PED $type: Sample $samp->{FID},$samp->{IID} are duplicated!";
		}
		unless($samp->{SEX} eq '1' || $samp->{SEX} eq '2' || $samp->{SEX} eq '0') {
			if ($strict) {
				die "Check PED $type: Sample $samp->{FID},$samp->{IID} has incorrect sex : $samp->{SEX}";
			}
			else {	
				warn "Check PED $type: Sample $samp->{FID},$samp->{IID} has incorrect sex : $samp->{SEX}";
				$samp->{SEX} = 0;
			}
		}
		$known{$samp->{FID},$samp->{IID}} = $samp->{SEX};
	}
	# Second round, check parent sex
	my $flag;
	foreach my $samp (@$samps) {
		if ($samp->{DAD} ne '0') {
			if (defined $known{$samp->{FID},$samp->{DAD}} && $known{$samp->{FID},$samp->{DAD}} ne '1') {
				if ($strict) {
					die "Check PED $type: Sample $samp->{FID},$samp->{IID}'s father has incorrect sex : $known{$samp->{FID},$samp->{DAD}}";
				}
				else {
					warn "Check PED $type: Sample $samp->{FID},$samp->{IID}'s father has incorrect sex : $known{$samp->{FID},$samp->{DAD}}";
					$known{$samp->{FID},$samp->{DAD}} = 1; $flag = 1;
				}
			}
			else {
				$known{$samp->{FID},$samp->{DAD}} = 1;
			}
		}
		if ($samp->{MOM} ne '0') {
			if (defined $known{$samp->{FID},$samp->{MOM}} && $known{$samp->{FID},$samp->{MOM}} ne '2') {
				if ($strict) {
					die "Check PED $type: Sample $samp->{FID},$samp->{IID}'s mother has incorrect sex : $known{$samp->{FID},$samp->{MOM}}";
				}
				else {
					warn "Check PED $type: Sample $samp->{FID},$samp->{IID}'s mother has incorrect sex : $known{$samp->{FID},$samp->{MOM}}";
					$known{$samp->{FID},$samp->{MOM}} = 2; $flag = 1;
				}
			}
			else {
				$known{$samp->{FID},$samp->{MOM}} = 2;
			}
		}
	}
	if ($flag) {
		foreach (@$samps) {
			if ($_->{SEX} ne $known{$_->{FID},$_->{IID}}) {
				$_->{SEX} = $known{$_->{FID},$_->{IID}}
			}
		}
	}
	# Third round, further check if no --id-delim is provided IIDs should be unique across families
	# and parents should also be unique across families even if they do not appear in IID column
	# %FIDs also store FamIDs for parents who do not appear on IID column
	my %FIDs;
	unless ($ARGV{'--id-delim'}) {
		foreach my $samp (@$samps) {
			unless (defined $FIDs{$samp->{IID}}) {
				$FIDs{$samp->{IID}} = $samp->{FID};
			}
			else {
				if ($FIDs{$samp->{IID}} ne $samp->{FID}) {
					if ($strict) {
						die "Check PED $type: Sample ID $samp->{IID} is not unique across families (no --id-delim)."; 
					}
					else {
						warn "Check PED $type: Sample ID $samp->{IID} is not unique across families (no --id-delim), only first FamID will be kept."; 
					}	
				}	
			}
			foreach my $PAR (qw|DAD MOM|) {
				if ($samp->{$PAR} ne '0') {
					unless(defined $FIDs{$samp->{$PAR}}) {
						$FIDs{$samp->{$PAR}} = $samp->{FID};
					}
					else {
						if ($FIDs{$samp->{$PAR}} ne $samp->{FID}) {
							my $par = lc($PAR);
							if ($strict) {
								die "Check PED $type: $samp->{IID}'s $par is in a different family : $FIDs{$samp->{$PAR}} <> $samp->{FID}";
							}
							else {
								warn "Check PED $type: $samp->{IID}'s $par is in a different family : $FIDs{$samp->{$PAR}} <> $samp->{FID}, only first FamID will be kept";
							}
						}
					}
				}
			}
		}
	}
	return \%FIDs;
}


__END__

=head1 NAME

update_ped_vcf.pl -- Update sample info in PED file.

=head1 DESCRIPTION

The script implement plink-style procedures to update sample info with some extensions.
It can update family ID, sample ID, parents, sexs, and phenotypes.

Order of operations:
	1. Sample removal will always be the first, so the remove list use the original IDs.
	2. FamID update, all individuals or selected in the family will have their FamID updated.
	3. IID update, i.e. individual rename (we extend plink's update-ids to seperately update FamID and IID).
	4. Sex update, should use updated FamID/IID in sex update list.
	4. Parents update, parents may or may not exist in input PED file, if exist sex should be compatible with role.
	5. Pheno update, should use updated FamID/IID in pheno update list.

If VCF is provided and sample has been removed or renamed, VCF will be modified to reflect changes.

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

=item -[-]famid [=] <char>

The default require IID to be unique across families. Turn on this option to allow IDs to be unique only within family.

=item -[-]patch [=] <prefix>

Instead of specifying patched files individually, they can be provided by a prefix
and using standard suffices as follows.

--remove		remove.txt
--update-fids 	update_fids.txt
--update-iids 	update_iids.txt
--update-sex	update_sex.txt
--update-pars 	update_pars.txt
--update-pheno	update_pheno.txt

=item -[-]remove [=] <list>

Expects list with one column (IID) or two columns (FamID, IID).

=item -[-]update-fids [=] <list>

Expects list with two columns (old FamID, new FamID) or three columns (Old FamID, old IID, new FamID).

=item -[-]update-iids [=] <list>

Expects list with two columns (Old IID, new IID) or three columns (Old FamID, Old IID, new IID).
Note: IDs in the parent columns will also be updated! 

=item -[-]update-sex [=] <list>

Expects list with two columns (IID, Sex) or three columns (FamID, IID, Sex)

=item -[-]update-pars [=] <list>

Expects list with three columns (IID, Dad, Mom) or four columns (FamID, IID, Dad, Mom)

=item -[-]update-pheno [=] <list>

Expect list with two columns (IID, Pheno) or three columns (FamID, IID, Pheno)

=back

=cut





