use strict;
use warnings;
use IO::Dir;
use IO::Compress::Gzip;
use FindBin qw|$Bin|;
use List::MoreUtils qw|all|;
use Utils::Hash qw|chk_default|;
use Utils::List qw|parse_fields|;
use Utils::File::Iter qw|iter_file|;
use Getopt::Euclid;

use lib "$Bin/../../lib";
use Shared qw|expand_dat|;


my @infiles;
if (defined $ARGV{'--nsplit'}) {
	unless(all { "$ARGV{'--prefix'}.$_.txt" } 1..$ARGV{'--nsplit'}) {
		die "Not all variant annotation tables can be found";
	}
	else {
		@infiles = map { "$ARGV{'--prefix'}.$_.txt" } 1..$ARGV{'--nsplit'};
	}
}
else {
	unless(-f $ARGV{'--prefix'} || $ARGV{'--prefix'} eq '-') {
		die "Cannot find input variant annotation table: $ARGV{'--prefix'}";
	}
	@infiles = ($ARGV{'--prefix'});
}

unless (-d $ARGV{'--outdir'}) {
	print STDERR "Creating output directory: $ARGV{'--outdir'}\n";
	mkdir $ARGV{'--outdir'};
}

# Max. dist to safely close gene-specific file handle
# We will assume input files are sorted by chrom,position
my $maxdist = $ARGV{'--max-dist'};

# Columns to include or exclude
my %exclude;
if ($ARGV{'--col-remove'}) {
	%exclude = map { $_ => 1 } parse_fields($ARGV{'--col-remove'});
}

my %varrm;
if (defined $ARGV{'--varrm'}) {
	open my $fin, $ARGV{'--varrm'} or die "Cannot open variant exclusion list";
	while(<$fin>) {
		my @dat = split;
		unless(@dat == 3) {
			die "Incorrect number of columns for variant removal list!";
		}
		my ($iid, $varid, $geneid) = @dat;
		unless($varid =~ /^\w+\:\d+\:[ACGT]+\:[ACGT]+$/) {
			die "Incorrect variant ID: $varid";
		}
		$varrm{$iid,$varid}{$geneid} = 1;
	}
}


# Expand packed annotations for overlapping genes
# Note that it is possible to have duplicated entry per gene (may be caused by VCF split from previous steps)
# It should be dealt by downstream analysis pipeline.
my (%genechr, %genepos, %fouts);

foreach my $infile (@infiles) {
	print STDERR "Processing $infile\n";
	my ($it, $fnames) = iter_file($infile eq '-' ? \*STDIN : $infile, { fsep => qr/\t/ });
	my @fields = grep { !defined $exclude{$_} } @$fnames;
	
	# Check standard fields
	foreach my $field (qw|Chrom Position Ref Alt GeneID GeneEff|) {
		unless(grep { $_ eq $field } @$fnames) {
			die "Cannot find all standard fields $field from $infile";
		}
	}
	# When varrm list is provided, should also contain IID
	if (%varrm) {
		unless(grep { $_ eq "IID" } @$fnames) {
			die "Cannot find sample ID field IID"
		}
	}

	while(my $dat = $it->()) {
		next if $dat->{GeneID} eq '.';
		my $varid = join(":", @{$dat}{qw|Chrom Position Ref Alt|});
		foreach my $gdat (expand_dat($dat, { sep => ';', optional => \@fields })) {
			my $gid = $gdat->{GeneID};
			if (%varrm) {
				next if defined $varrm{$dat->{IID},$varid}{$gid};
			}
			unless (defined $fouts{$gid}) {
				if(-f "$ARGV{'--outdir'}/$gid") {
					die "Gene table $ARGV{'--outdir'}/$gid already exists: at $dat->{Chrom}:$dat->{Position}!";
				}
				# Before open up new files, close tables for genes that are enough farther away
				foreach my $gid (keys %fouts) {
					if ($dat->{Chrom} ne $genechr{$gid} || $dat->{Position}-$genepos{$gid} > $maxdist) {
						delete $fouts{$gid};
					}
				}
				if ($ARGV{'--gzip'}) {
					if ($ARGV{'--suffix'}) {
						$fouts{$gid} = IO::Compress::Gzip->new("$ARGV{'--outdir'}/${gid}.$ARGV{'--suffix'}.gz")
							or die "Cannot write to compressed file $ARGV{'--outdir'}/${gid}.$ARGV{'--suffix'}.gz";
					}
					else {
						$fouts{$gid} = IO::Compress::Gzip->new("$ARGV{'--outdir'}/$gid.gz")
							or die "Cannot write to compressed file $ARGV{'--outdir'}/$gid.gz";
					}
				}
				else {
					if ($ARGV{'--suffix'}) {
						$fouts{$gid} = IO::File->new("$ARGV{'--outdir'}/${gid}.$ARGV{'--suffix'}", "w") 
							or die "Cannot write to output $ARGV{'--outdir'}/${gid}.$ARGV{'--suffix'}";
					}
					else {
						$fouts{$gid} = IO::File->new("$ARGV{'--outdir'}/$gid", "w") 
							or die "Cannot write to output $ARGV{'--outdir'}/$gid";
					}		
				}
				
				unless(defined $fouts{$gid}) {
					print scalar(keys %fouts), "\n";
					die "Cannot write $gid to $ARGV{'--outdir'}!";
				}
				$fouts{$gid}->print(join("\t", @$fnames), "\n");
			}

			# Add current location to gene table
			chk_default(\%genechr, $gid, $gdat->{Chrom});
			$genepos{$gid} = $gdat->{Position};

			my $varid = join(":", @{$gdat}{qw|Chrom Position Ref Alt|});
			$fouts{$gid}->print(join("\t", @{$gdat}{@$fnames}), "\n");
		}
	}
}

__END__

=head1 NAME

extract_anno_pergene.pl -- Extract annotated variants per-gene from a variant table.


=head1 REQUIRED ARGUMENTS

=over

=item -[-]prefix [=] <prefix>

Input table file or file prefix. Should be sorted by Chrom,Position.

=item -[-]outdir [=] <directory>

Output directory.

=back

=head1 OPTIONS

=over

=item -[-]suffix [=] <string>

Output file name suffix. Default is null or gz  when gzip is enabled.
Then suffix is provided, output will become GeneID.suffix

=item -[-]nsplit [=] <count>

Number of splits for the input. The input file will be named as "prefix.$_.txt".
Otherwise, input file should be named as prefix.

=item -[-]max[-]dist [=] <distance>

Max. distance for closing a gene file handle.

=for Euclid:
	distance.default: 1000000

=item -[-]gzip

Gzip output files.

=item -[-]varrm [=] <list>

Gene-specific variant removal list, format: VarID,IID,Gene.

=item -[-]col[-]remove [=] <columns>

Columns to be excluded from expansion.

=back



