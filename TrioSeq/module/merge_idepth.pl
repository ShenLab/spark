use strict;
use warnings;
use Utils::File::Iter qw|iter_file|;
use List::Util qw|sum|;

# Merge idepth calculated on chrX, Y, and chr21, and autosomes

my ($prefix, $pedfile, $outfile) = @ARGV;
unless(defined $outfile) {
	$outfile = "${prefix}_merged.txt";
}

my %sex;
{
	unless(defined $pedfile) {
		die "Must provide pedfile!";
	}
	my %code2gender = ('1' => 'Male', 2 => 'Female');
	open my $fin, $pedfile or die "Cannot open $pedfile";
	while(<$fin>) {
		my ($iid, $sexcode) = (split)[1,4];
		$sex{$iid} = $code2gender{$sexcode} // "Unknown";
	}
}

my (%idepth, %dpbysamp);
foreach my $chrom (qw|chr21 chrX chrY auto|) {
	my @idpfiles = glob("${prefix}_$chrom.*.idepth");
	if (-f "${prefix}_$chrom.idepth") {
		push @idpfiles, "${prefix}_$chrom.idepth"
	}
	unless(@idpfiles) {
		if ($chrom eq 'chrY') {
			warn "Cannot find idepth files for $chrom!";
		}
		else {
			die "Cannot find idepth files for $chrom!";
		}
	}
	foreach my $idpfile (@idpfiles) {
		slurp_idepth($idpfile, $chrom);
	}
	# Take weighted average of splits
	foreach my $iid (sort keys %idepth) {
		$dpbysamp{$iid}{$chrom} = sum(map { $_->[0]*$_->[1] } @{$idepth{$iid}{$chrom}})/
									sum(map { $_->[0] } @{$idepth{$iid}{$chrom}});
	}
}


open my $fout, ">$outfile" or die "Cannot write to $outfile";
print $fout join("\t", qw|IID DP_Auto DP_Chr21 DP_ChrX DP_ChrY Gender|), "\n";
foreach my $iid (sort keys %dpbysamp) {
	print $fout join("\t", $iid,  (map { $dpbysamp{$iid}{$_} // "." } qw|auto chr21 chrX chrY|), $sex{$iid} // "Unknown"), "\n";
}


# INDV N_SITES MEAN_DEPTH
sub slurp_idepth {
	my ($file, $chrom) = @_;
	my $it = iter_file($file) or die "Cannot iter over $file";
	while(my $dat = $it->()) {
		next unless $dat->{N_SITES} > 0;
		unless($dat->{MEAN_DEPTH} =~ /^[0-9\.]+$/) {
			die "$file: Incorrect mean depth for $dat->{INDV}: $dat->{MEAN_DEPTH}";
		}
		push @{$idepth{$dat->{INDV}}{$chrom}}, [$dat->{N_SITES}, $dat->{MEAN_DEPTH}];
	}
	return 1;
}
