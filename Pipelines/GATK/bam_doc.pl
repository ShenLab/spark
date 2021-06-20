#!/usr/bin/env perl

use strict;
use warnings;
use Carp;
use POSIX qw|ceil|;
use IO::File;
use Data::Dumper;
use FindBin qw|$Bin|;
use File::Copy;
use File::Path qw|make_path|;
use Cwd qw|abs_path|;
use Getopt::Lucid qw|:all|;
use List::Util qw|min max|;
use List::MoreUtils qw|all none uniq|;
use Config::Std;
use Utils::Workflow;
use Utils::Hash qw|merge_conf|;
use Utils::File qw|count_line|;

use lib "$Bin/../lib/";
use Shared qw|read_list read_chrlen|;

my @spec =  (
	Param("conf|c")->valid(sub { -r }),
	Param("list|l")->valid(sub { -r }),
	Param("rename")->valid(sub { -r }),
	Param("remove")->valid(sub { -r }),
	Param("outdir|out"),
	Param("engine|eng")->valid(sub { $_ eq 'SGE' || $_ eq 'BASH' }),
	Keypair("param|par"),
	Switch("help|h"),
	Switch("dryrun|dry"),
	Switch("force")
	);

my $opt = Getopt::Lucid->getopt(\@spec);

if ($opt->get_help) {
	print STDERR <<EOF;
Purpose:
	This is a pipeline script to calculate depth of coverage (DoC) from BAM/CRAMs.

Usage:
	bam_doc.pl --conf Conf --list BAMList --outdir OutDir 

Options:
	--list: BAM/CRAM file list or directory. The list should have 1 or 2 columns per line: path to 
		BAM/CRAM file (required) and sample ID. If only one columns is available, sample ID will 
		be taken from the basename after removing suffix. If a directory is provided, we will 
		look for BAM/CRAM files from the directory. Only one file type is allowed.
	--rename: Sample rename list.
	--remove: A list of bad samples to be removed.

Note:
	The script supports DoC calculation using mosdepth and/or GATK. There are rich options for DoC calculation
	by each tools. To allow flexibility in downstream analysis, it's possible to apply different tools 
	and options for calculating DoC for the same individual BAM/CRAM file. 

	In the config file, apart from common options, different tool and options are organized into multiple
	sections. Each section starting with 'MOSDEPTH' or 'GATK' will represent a combination of tool and
	associated options. All combinations of tools and options will be applied to each individual BAM/CRAMs. 
	Each combination will be associated with an output directory with the same name (in lower cases) to
	store per-individual DoC results. Output for each individual will have sample ID as basename. It's
	possible to have multiple output files per sample. The per-target DoC is the required output for each
	tool-options combination, and will be reformatted to tabix indexed BED4 format with .cov.bed.gz suffix.

	For both tools, output of per-base DoC was disabled by default to save time and space. 

Dependencies:
	mosdepth, gatk (v3, invoked by GenomeAnalysisTK)
	Note: At the time this script was developed, GATK4 did not support DepthOfCoverage function.
		  So it makes use of GATK v3 and must be invoked by calling GenomeAnalysisTK.

EOF
	exit 1;
}

$opt->validate({requires => [qw|conf list outdir|]});

my $rootdir = $opt->get_outdir;

my %conf = merge_conf($opt->get_conf, $opt->get_param); 

# Also determine the DoC sets to be calcualted
my @sets = grep { /^(MOSDEPTH|GATK)/ } sort keys %conf;
if (@sets > 0) {
	print STDERR "The following tools for calculating DoC can be found in config:\n";
	print STDERR join("\t", @sets), "\n"; 
}
else {
	print STDERR "No DoC set was found in config.";
	exit 1;
}


my $wkf = Utils::Workflow->new($rootdir,
	{ dir => [map { lc($_) } @sets], strict_var => 1, 
	  engine => $opt->get_engine, force => $opt->get_force });


while(my ($label, $file) = each %{$conf{PATH}}) {
	unless(-f $file || -d $file) {
		croak "Cannot find $label file $file";
	}
}

############################################
## Input files validation and collection  ##
############################################
# Skip testing for files on s3
my %opt = exists $conf{AWS} && exists $conf{AWS}{STAGEIN} ? (notest => 1) : ();

my ($bams, $filetype) = read_list($opt->get_list,
	{ suffix => { bam => qr/\.bam$/, cram => qr/\.cram$/ }, 
	  rename => $opt->get_rename, remove => $opt->get_remove, %opt });

my @samps = sort keys %$bams;
unless (@samps > 0) {
	print STDERR "No sample was found in the list\n";
	exit 1;
}

# Write to sample list
for(my $ii = 1; $ii <= @samps; $ii ++) {
	open my $fout, ">$rootdir/par/IID2BAM.$ii" or die "Cannot write to IID2BAM.$ii";
	print $fout $samps[$ii-1], "\t", $bams->{$samps[$ii-1]}, "\n";
}

#################################
##  Workflow initialization    ##
##  Working directory setup    ##
#################################

my $chrlen = read_chrlen($conf{PATH}{SEQDICT});

# Check targets files, and determine the number of targeted intervals (#lines in output)
my %ntargs;
foreach my $set (@sets) {
	unless (ref $conf{$set}{OPTION}) {
		my @dat = ($conf{$set}{OPTION});
		$conf{$set}{OPTION} = [@dat];
	}
	my $target;
	foreach my $option (@{$conf{$set}{OPTION}}) {
		if ($set =~ /^MOSDEPTH/) {
			if ($option =~ /\-\-by (\S+)/ || $option =~ /\-b (\S+)/) {
				$target = $1;
				if ($target =~ /^\d+$/) {
					$ntargs{$set} = num_of_targs($chrlen, $target);
				}
				elsif (-f $target) {
					$ntargs{$set} = count_line($target);
				}
				else {
					croak "Incorrect target specification: -b $target";
				}
				last;
			}
		}
		elsif ($set =~ /^GATK/) {
			if ($option =~ /\-\-intervals (\S+)/ || $option =~ /\-L (\S+)/) {
				$target = $1;
				if (-f $target) {
					$ntargs{$set} = count_line($target);
				}
				else {
					# Note we only support interval files followed by -L
					croak "Cannot find target interval list: $target";
				}
				last;
			}
		}
		else {
			croak "Cannot recognize set: $set";
		}
	}
	unless (defined $target) {
		unless (exists $conf{PATH}{TARGETS}) {
			croak "Cannot find target specification!";
		}
		else {
			$ntargs{$set} = count_line($conf{PATH}{TARGETS});
			if ($set =~ /^MOSDEPTH/) {
				push @{$conf{$set}{OPTION}}, "--by $conf{PATH}{TARGETS}";
			}
			elsif ($set =~ /^GATK/) {
				push @{$conf{$set}{OPTION}}, "-L $conf{PATH}{TARGETS}";
			}
		}
	}
}

write_config %conf, "$rootdir/par/run.conf" unless $opt->get_dryrun;

foreach my $set (@sets) {
	my $tag = ucfirst(lc($set));
	my ($script, $suffix) = calc_doc($set);
	$wkf->add($script, { name => $tag,
			expect => [ exp_bysamp(lc($set), @$suffix) ],
			nslots => scalar(@samps)  })
		->add(refmt_doc($set), { name => "Refmt".$tag, 
			depend => $tag, deparray => 1,
			expect => [ exp_bysamp(lc($set), "cov.bed.gz", "cov.bed.gz.tbi") ], 
			nslots => scalar(@samps) });
}


$wkf->inst(\%conf);
$wkf->run({ conf => $conf{$wkf->{engine}}, dryrun => $opt->get_dryrun  });


############################
## Workflow components    ##
############################

# Mosdepth
sub calc_doc {
	my ($docset) = @_;
	my $suffix;
	my $script = "\nread IID BAMFILE < _PARDIR_/IID2BAM._INDEX_\n";
	if (exists $conf{AWS}{STAGEIN}) {
		$script .= <<'EOF';
TMPDIR=$(stagein.pl --input $BAMFILE --output _TMPDIR_/BAMFILE._INDEX_ --profile _AWS.PROFILEIN_ --tmpdir _AWS.STAGEIN_ --all)

trap "echo Clean up $TMPDIR; rm -fR $TMPDIR; exit 1" EXIT SIGINT SIGTERM

read BAMFILE < _TMPDIR_/BAMFILE._INDEX_

EOF
	}

	if ($docset =~ /^MOSDEPTH/) {
		$script .= <<'EOF';

mosdepth -n --fasta _PATH.FASTA_ _DOCSET.OPTION[ \
	]_ _DOCSETDIR_/$IID $BAMFILE

EOF
		$suffix = [qw|regions.bed.gz regions.bed.gz.csi|];
	}
	elsif ($docset =~ /^GATK/) {
		$script .= <<'EOF';

# Note: GATK4 has not yet ported the DepthOfCoverage module
GenomeAnalysisTK _DOCSET.GATKOPT_ \
	-T DepthOfCoverage -l INFO \
	--omitDepthOutputAtEachBase \
	-R _PATH.FASTA _DOCSET.OPTION[ \
	]_ -o _DOCSETDIR_/$IID -I $BAMFILE

EOF
		$suffix = [qw|sample_interval_summary sample_interval_statistics|];
	}
	else {
		croak "Cannot recognize set: $docset";
	}
	$script =~ s/docset/${(\ lc($docset) )}/g;
	$script =~ s/DOCSET/${(\ uc($docset) )}/g;
	return ($script, $suffix);
}


sub refmt_doc {
	my ($docset) = @_;
	my $script = "\nread IID BAMFILE < _PARDIR_/IID2BAM._INDEX_\n";
	if ($docset =~ /^MOSDEPTH/) {
		$script .= <<'EOF';

NCOL=$(head -n1 <(zcat _DOCSETDIR_/$IID.regions.bed.gz) | awk '{print NF}')

if [[ $NCOL == 5 ]]; then
	zcat _DOCSETDIR_/$IID.regions.bed.gz | cut -f1-3,5 | bgzip -c > _DOCSETDIR_/$IID.cov.bed.gz
elif [[ $NCOL == 4 ]]; then
	echo "Creat symlink for $IID.cov.bed.gz"
	ln -s $IID.regions.bed.gz _DOCSETDIR_/$IID.cov.bed.gz
else 
	echo "Incorrect number of columns in $IID.regions.bed.gz" 1>&2
	exit 1
fi

tabix -f -p bed _DOCSETDIR_/$IID.cov.bed.gz

EOF
		if (exists $conf{AWS}{STAGEOUT}) {
			$script .= <<'EOF';
stageout.pl --files $IID.regions.bed.gz,$IID.regions.bed.gz.csi,$IID.cov.bed.gz,$IID.cov.bed.gz.tbi \
	--profile _AWS.PROFILEOUT_ --indir _DOCSETDIR_ --outdir _AWS.STAGEOUT_/docset

EOF
		}
	}
	elsif ($docset =~ /^GATK/) {
		$script .= <<'EOF';

perl -ape '@chrpos = ($F[0]=~/(\w+):(\d+)-(\d+)/); $chrpos[1]-=1; $_=join("\t", @chrpos, $F[2])."\n"' \
	_DOCSETDIR_/$IID.sample_interval_summary | awk 'NR>1' | bgzip -c > _DOCSETDIR_/$IID.cov.bed.gz

tabix -f -p bed _DOCSETDIR_/$IID.cov.bed.gz

EOF
		if (exists $conf{AWS}{STAGEOUT}) {
			$script .= <<'EOF';
stageout.pl --files $IID.sample_interval_summary,$IID.sample_interval_statistics,$IID.cov.bed.gz,$IID.cov.bed.gz.tbi \
	--profile _AWS.PROFILEOUT_ --indir _DOCSETDIR_ --outdir _AWS.STAGEOUT_/docset

EOF
		}
	}
	else {
		croak "Cannot recognize set: $docset";
	}
	$script =~ s/docset/${(\ lc($docset) )}/g;
	$script =~ s/DOCSET/${(\ uc($docset) )}/g;
	return $script;
}

sub exp_bysamp {
	my $dir = shift @_;
	my @exp;
	foreach my $iid (@samps) {
		push @exp, [map { "$dir/$iid.$_" } @_]
	} 
	return @exp;
}

sub exp_checker {
	my ($docset) = @_;
	my $nline = $ntargs{$docset};
	return sub {
		if (all { -f "$rootdir/$_" } @_) {
			my $docfile = shift @_;
			if ($docfile =~ /\.bed\.gz$/) {
				return count_line($rootdir."/".$docfile) == $nline ? 1 : 0;
			}
			else {
				return count_line($rootdir."/".$docfile) == $nline + 1 ? 1 : 0;
			}
		}
		else {
			return 0;
		}
	}
}

sub num_of_targs {
	my ($chrlen, $binsize) = @_;
	my $ntot;
	foreach my $len (values %$chrlen) {
		$ntot += ceil($len/$binsize);
	}
	return $ntot;
}

