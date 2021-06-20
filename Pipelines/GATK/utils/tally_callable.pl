use strict;
use warnings;
use Perl6::Slurp;
use Getopt::Euclid;
use Utils::Number qw|commafy|;
use List::MoreUtils qw|all uniq|;
use Genome::Ranges::IntSet;
use Genome::UCSC::BinKeeper;

my $dupid = qr/$ARGV{'--ignore'}/;

my %remove;
if (defined $ARGV{'--remove'}) {
	%remove = map { (split)[0] => 1 } slurp $ARGV{'--remove'};
}

# Individual group membership: 
# $gsamp{$iid}{$gid} = 1
my %gsamp;
if (defined $ARGV{'--group'}) {
	open my $fin, $ARGV{'--group'} or die "Cannot open group file";
	while(<$fin>) {
		my ($iid, $gid) = split;
		unless(defined $iid && defined $gid) {
			die "Incorrect line of group file: $_"
		}
		my $origid = $iid;
		if ($iid =~ /$dupid$/) {
			$origid = s/$dupid$//;
		}
		$gsamp{$origid}{$gid} = 1;
	}
}

my $bk = Genome::UCSC::BinKeeper->new();

# Keep track of samples with callable region file in the list
my %known;
open my $fin, $ARGV{'--list'} or die "Cannot open list file";
while(<$fin>) {
	my ($gzbed, $iid) = (split)[0,1];
	unless(defined $gzbed && defined $iid) {
		die "Incorrect number of columns for list file: $_";
	}
	unless (-f $gzbed) {
		die "Cannot find gz'ed bed file $gzbed";
	}
	next if defined $remove{$iid};
	if ($iid =~ /$dupid$/) {
		$iid =~ s/$dupid$//;
	}
	$known{$iid} = 1;
	open my $fh, "tabix -p bed -R $ARGV{'--bed'} $gzbed |" or die "Cannot open pipe";
	while(<$fh>) {
		my ($chr, $start, $end, $value) = (split)[0,1,2,3];
		if ($value eq $ARGV{'--label'}) {
			$bk->add($chr, $start+1, $end, $iid);
		}
	}
}

my (@samps, $ntot, %group);
if ($ARGV{'--group'}) {
	@samps = grep { defined $gsamp{$_} } keys %known;
	foreach my $iid (@samps) {
		foreach my $gid (keys %{$gsamp{$iid}}) {
			push @{$group{$gid}} => $iid;
		}
	}
	#$ntot = scalar(uniq sort map { keys %{$gsamp{$_}} } @samps);
	$ntot = scalar(keys %group);
	print STDERR "A total of @{[ scalar(@samps) ]} samples found in $ntot groups\n";
	my $nsamprm = grep { !defined $gsamp{$_} } keys %known;
	if ($nsamprm > 0) {
		print STDERR "A total of $nsamprm samples are not associated with any group and will be ignored\n";
	}
}
else {
	@samps = sort keys %known;
	$ntot = scalar(@samps);
	print STDERR "Total number of samples found in the list: $ntot\n";
}
unless ($ntot) {
	exit 1;
}

my $set = Genome::Ranges::IntSet->new($ARGV{'--bed'}, {bed => 1});
if (defined $ARGV{'--mask'}) {
	my $mask = Genome::Ranges::IntSet->new($ARGV{'--mask'}, {bed => 1});
	$set *= $mask;
	print STDERR "Total base pair length of the targeted regions after mask: @{[ commafy($set->size()) ]}\n";
}
else {
	print STDERR "Total base pair length of the targeted regions: @{[ commafy($set->size()) ]}\n";
}


my $fout;
if (defined $ARGV{'--output'}) {
	open $fout, ">$ARGV{'--output'}" or die "Cannot write to $ARGV{'--output'}";
}
else {
	$fout = \*STDOUT;
}

# Iterate over each base pair of masked target region
my $it = $set->iter();
while(my ($chr, $pos) = $it->()) {
	#print $chr, "\t", $pos, "\n";
	# Find all overlapping samples
	my @overlaps = uniq sort map { $_->[2] } $bk->find_range($chr, $pos);
	my $ngroup = 0;
	if ($ARGV{'--group'}) {
		my %hits; @hits{@overlaps} = @overlaps;
		foreach my $gid (keys %group) {
			# all samples in the group
			if (all { defined $hits{$_} } @{$group{$gid}}) {
				$ngroup ++;
			}
		}
	}
	else {
		$ngroup = scalar(@overlaps);
	}
	my $perc = sprintf("%.4f", $ngroup/$ntot);
	print $fout  join("\t", $chr, $pos, $perc), "\n";
}



__END__

=head1 NAME

tally_callable -- Tally percentage of callable samples/groups at each position across target regions. 

=head1 REQUIRED ARGUMENTS

=over
 
=item -[-]list [=] <file>

Interval list. bgzip'ed BEDs file of DP-quantized regions and sample ID list.
DP-quantized regions should be bgzipped and tabix indexed (MosDepth output).

=for Euclid:
	file.type: readable

=item -[-]bed [=] <file>

Targeted regions BED file. In parallel enviroment, this is a subet of full target regions. 

=for Euclid:
	file.type: readable

=back
 
=head1 OPTIONS
 
=over

=item -[-]mask [=] <file>

Mask region BED file. Callable regions will be intersected with the mask regions.

=item -[-]group [=] <file>

Two columns: sample ID and group ID. Each sample is allowed to be present in multiple groups.
When group info is provied, sample not in the group will be ignored. A region is callable for
the group only if ALL samples are callable for this region.

=item -[-]out[put] [=] <file>

The output file. The output will have three or three columns: chr, pos, and perc.
Perc is the percentage of callable samples/groups at this position. When group is proivided
only group level statistics will be in the output

=item -[-]label [=] <label>

Label for callable intervals

=for Euclid:
	label.default: "CALLABLE"

=item -[-]ignore [=] <pattern>

Regex pattern for duplicated sample. Callable region of dups will be merged with original sample.

=for Euclid:
	pattern.default: '_Re\d*$'

=item -[-]remove [=] <list>

A list of bad samples. It should contain sample IDs.

=for Euclid:
	list.type: readable
 
=back





