use strict;
use warnings;
use Perl6::Slurp;
use Utils::Number qw|commafy|;
use Getopt::Euclid;
use Utils::File qw|open_file|;
use List::MoreUtils qw|uniq|;
use Genome::Ranges::IntSet;

my $dupid = qr/$ARGV{'--ignore'}/;

my %remove;
if (defined $ARGV{'--remove'}) {
	%remove = map { (split)[0] => 1 } slurp $ARGV{'--remove'};
}

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

# Keep gzbed file for each unique sample
my %list;
open my $fin, $ARGV{'--list'} or die "Cannot open list file";
while(<$fin>) {
	my ($gzbed, $iid) = (split)[0,1];
	unless(defined $gzbed && defined $iid) {
		die "Incorrect number of columns for list file: $_";
	}
	unless (-f $gzbed) {
		die "Cannot find gz'ed bed file: $gzbed";
	}
	next if defined $remove{$iid};
	if ($iid =~ /$dupid$/) {
		$iid =~ s/$dupid$//;
	}
	push @{$list{$iid}}, $gzbed;
}

my (@samps, $ntot, %group);
if ($ARGV{'--group'}) {
	@samps = grep { defined $gsamp{$_} } keys %list;
	foreach my $iid (@samps) {
		foreach my $gid (keys %{$gsamp{$iid}}) {
			push @{$group{$gid}} => $iid;
		}
	}
	#$ntot = scalar(uniq sort map { @{$group{$_}} } @samps);
	$ntot = scalar(keys %group);
	print STDERR "A total of @{[ scalar(@samps) ]} samples found in $ntot groups\n";
	my $nsamprm = grep { !defined $gsamp{$_} } keys %list;
	if ($nsamprm > 0) {
		print STDERR "A total of $nsamprm samples are not associated with any group and will be ignored\n";
	}
}
else {
	@samps = sort keys %list;
	$ntot = scalar(@samps);
	print STDERR "Total number of samples found in the list: $ntot\n";
}
unless($ntot) {
	exit 1;
}

my $target = bed_to_set($ARGV{'--bed'});
my $mask;
if (defined $ARGV{'--mask'}) {
	$mask = Genome::Ranges::IntSet->new($ARGV{'--mask'}, {bed => 1});
	#$target *= $mask;
}
my $tglen = $target->size();
print STDERR "Total base pair length of the targeted regions: @{[ commafy($tglen) ]}\n";

my $fout;
if (defined $ARGV{'--output'}) {
	open $fout, ">$ARGV{'--output'}" or die "Cannot write to $ARGV{'--output'}";
}
else {
	$fout = \*STDOUT;
}

my $outdir;
if ($ARGV{'--outdir'}) {
	$outdir = $ARGV{'--outdir'};
	mkdir $outdir unless -d $outdir;
}

foreach my $gid (sort keys %group) {
	my $gset;
	foreach my $iid (@{$group{$gid}}) {
		my $iset = bed_to_set($list{$iid}[0], $target, $ARGV{'--label'});
		if (@{$list{$iid}} > 1) {
			for(my $ii = 1; $ii < @{$list{$iid}}; $ii ++) {
				# For duplicated sample, we take union
				$iset += bed_to_set($list{$iid}[$ii], $target, $ARGV{'--label'});
			}
		}
		# For samples in the group, we take intersection
		if (defined $gset) {
			$gset *= $iset;
		}
		else {
			$gset = $iset;
		}
	}
	# Final intersection with target (mask should be applied before)
	$gset *= $target;
	$gset *= $mask if $ARGV{'--mask'};
	print $fout $gid, "\t", sprintf("%.4f", $gset->size()/$tglen), "\n";
	if ($outdir) {
		$gset->write("$outdir/$gid.callable.bed", { bed => 1 });
	}
}

sub bed_to_set {
	my ($bedfile, $targs, $label) = @_;
	my %ranges;
	my $fbed = open_file($bedfile);
	while(<$fbed>) {
		my ($chr, $start, $end, $value) = (split)[0,1,2,3];	
		if (defined $targs) {
			next unless defined $targs->{$chr};
		}
		if (defined $label) {
			next unless $value eq $label;
		}
		push @{$ranges{$chr}}, [$start+1, $end];
	}
	my $set = Genome::Ranges::IntSet->new(\%ranges);
	return $set;
}


__END__

=head1 NAME

collapse_callable -- "Collapse" callable regions per sample/group.

=head1 REQUIRED ARGUMENTS
 
=over
 
=item -[-]list [=] <file>

Interval list. bgzip'ed BED files of DP-quantized regions and sample ID list.
Typically, this is a subet of full region-sample list in parallel enviroment.
Note: if region-sample list was split, the duplicated sample or group sample
should be in the same split.

=for Euclid:
	file.type: readable

=item -[-]bed [=] <file>

Targeted regions BED file. We will calculate percetage of target regions that are callable.

=for Euclid:
	file.type: readable

=back

=head1 OPTIONS
 
=over

=item -[-]mask [=] <file>

Mask region BED file. We will only count callable based within mask file.

=item -[-]group [=] <file>

Two columns: sample ID and group ID. Each sample is allowed to be present in multiple groups.
When group info is provided, samples not in the group will be ignored.

=item -[-]out[put] [=] <file>

The output file. The output will have two columns: sample/group, perc.
Perc is the percentage of targeted regions that are callable.
When group is provided only group level statistics will be in the output.

=item -[-]outdir [=] <dir>

If an output directory is provided, then per-sample/group callable regions will be written
to BED files (one file per sample or group).

=item -[-]label [=] <label>

Label for callable intervals

=for Euclid:
	label.default: "CALLABLE"

=item -[-]ignore [=] <pattern>

Regex pattern for duplicated sample. Duplicates will be renamed and added
to the callable region of original sample.

=for Euclid:
	pattern.default: '_Re\d*$'

=item -[-]remove [=] <list>

A list of bad samples. It should contain sample IDs.

=for Euclid:
	list.type: readable
 
=back




