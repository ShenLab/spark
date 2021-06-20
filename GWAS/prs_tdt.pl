use strict;
use warnings;
use Data::Dumper;
use Getopt::Euclid;
use Perl6::Slurp;
use Math::Gauss qw|cdf|;
use Utils::Stat qw|mean std|;
use Utils::File::Iter qw|iter_file|;


# Slurp PRS
my @thres;
if ($ARGV{'--p-vals'}) {
	@thres = split(',', $ARGV{'--p-vals'});
}

# Standard output from PRSice
my $f_id = "IID";

# Slurp PRS and calculate difference between child and mid-parents
my (%PRS, @pvalthres);
{
	my ($it, $fnames) = iter_file($ARGV{'--profile'}, { chomp => 3 });
	if (@thres) {
		foreach my $thres (@thres) {
			unless(grep { $_ =~ /^[0-9\.]+$/ && $_ == $thres } @$fnames) {
				warn "Cannot find p-value threshold $thres in the table";
			}
			else {
				push @pvalthres, $thres;
			}
		}
	}
	else {
		@pvalthres = grep { /^[0-9\.]+$/ } @$fnames;
	}
	unless (@pvalthres > 0) {
		warn "No p-value threshold left!";
		exit 1;
	}
	while(my $dat = $it->()) {
		$PRS{$dat->{$f_id}} = [@{$dat}{@pvalthres}];
	}
}

my (%keep, %remove);
if ($ARGV{'--remove'}) {
	if ($ARGV{'--no-fid'}) {
		%remove = map { (split)[0] => 1 } slurp $ARGV{'--remove'};
	}
	else {
		%remove = map { (split)[1] => 1 } slurp $ARGV{'--remove'};
	}
}
if ($ARGV{'--keep'}) {
	if ($ARGV{'--no-fid'}) {
		%keep = map { (split)[0] => 1 } slurp $ARGV{'--keep'};
	}
	else {
		%keep = map { (split)[1] => 1 } slurp $ARGV{'--keep'};
	}
}

# Find samples with both parents available.
my (%trios, %famid, %sex, %pheno);
{
	open my $fin, $ARGV{'--ped'} or die "Cannot open PED file";
	while(<$fin>) {
		my ($fid, $iid, $dad, $mom, $sex, $phe) = split;
		next if defined $remove{$iid};
		if ($ARGV{'--keep'}) {
			next unless defined $keep{$iid};
		}
		if (defined $PRS{$iid} && defined $PRS{$dad} && defined $PRS{$mom}) {
			$trios{$iid} = [$dad, $mom];
			$famid{$iid} = $fid;
			$sex{$iid} = $sex == 1 ? "Male" : $sex == 2 ? "Female" : "Unknown"; 
			$pheno{$iid} = $phe == 1 ? "Unaffected" : $phe == 2 ? "Affected" : "Unknown";
		}
	}
}

my %sampgroup;
if (defined $ARGV{'--group'}) {
	open my $fin, $ARGV{'--group'} or die "Cannot open $ARGV{'--group'}";
	while(<$fin>) {
		my ($iid, $label) = split;
		$sampgroup{$iid}{$label} = 1; 
	}
}
else {
	foreach my $iid (keys %pheno) {
		$sampgroup{$iid}{$pheno{$iid}} = 1;
	}
}

# Calculate PRS-TD, and output individual and summary level results at each threshold.
open my $fout, ">$ARGV{'--output'}.prstd.txt" or die "Cannot write to $ARGV{'--output'}.prstd.txt";
print $fout join("\t", qw|FamID IID Sex Pheno|, @pvalthres), "\n";

my (%PTDs, @MidPars);
foreach my $iid (sort keys %PRS) {
	next unless defined $trios{$iid} && defined $sampgroup{$iid};
	my @iPTDs;
	for(my $ii = 0; $ii < @pvalthres; $ii ++) {
		my $midpar = 0.5*($PRS{$trios{$iid}[0]}[$ii]+$PRS{$trios{$iid}[1]}[$ii]);
		#push @{$MidPars[$ii]}, $midpar;
		my $parents = join(",", @{$trios{$iid}});
		$MidPars[$ii]{$parents} = $midpar;
		my $PRSTD = $PRS{$iid}[$ii]-$midpar;
		push @iPTDs, $PRSTD;
		foreach my $group (keys %{$sampgroup{$iid}}) {
			push @{$PTDs{$group}[$ii]}, $PRSTD;

		}
	}
	print $fout join("\t", $famid{$iid}, $iid, $sex{$iid}, $pheno{$iid}, @iPTDs), "\n";
}

open my $fsum, ">$ARGV{'--output'}.summary.txt" or die "Cannot write to $ARGV{'--output'}.summary.txt";
print $fsum join("\t", qw|Group N Thres Mean SE Zscore Pvalue|), "\n";
foreach my $group (sort keys %PTDs) {
	for(my $ii = 0; $ii < @pvalthres; $ii ++) {
		my $N = scalar(@{$PTDs{$group}[$ii]});
		my $midparSE = std( values %{$MidPars[$ii]} );
		my $mean = mean(@{$PTDs{$group}[$ii]})/$midparSE;
		my $SE = std(@{$PTDs{$group}[$ii]})/sqrt($N)/$midparSE;
		my $Zsc = $mean/$SE;
		print $fsum join("\t", $group, $N, $pvalthres[$ii], $mean, $SE, $Zsc, 1-cdf($Zsc)), "\n";
	}
}

__END__

=head1 NAME

prs_tdt.pl -- Polygenic score transmission disequilibrium.

=head1 NOTES

=head1 REQUIRED ARGUMENTS

=over

=item -[-]profile [=] <table>

PRS profile scores. It should be a table from the output of PRsice.

=item -[-]ped [=] <pedfile>

Pedigree file. We will extract complete trios whose PRS are available for analysis.

=item -[-]out[put] [=] <prefix>

Output file prefix. Two files will be created: 
prefix.prstd.txt lists delta-PRS for each trio offspring (compared with mid-parent).
prefix.summary.txt summarize PRS-TD for each sample groups.

=back

=head1 OPTIONS

=over

=item -[-]keep [=] <list>

List of samples to be included. It should have two columns (FID,IID) or one column (IID) if --no-fid is turned on.
Only IID will be used.

=item -[-]remove [=] <list>

List of samples to be removed. It should have two columns (FID,IID) or one column (IID) if --no-fid is turned on.
Only IID will be used.

=item -[-]no-fid

Indicate that keep/remove file does not have FID.

=item -[-]group [=] <list>

Group membership for offspring sample, two columns: IID => Group.
One sample can exist in multiple groups. Default is to group sample based on phenotype in PED file.

=item -[-]p[-]val[s] [=] <thres>

Select p-value thresholds, comma separate list. Should match column names in PRsice output.

=back

=cut
