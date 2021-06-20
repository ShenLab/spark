#!/usr/bin/env perl

use strict;
use warnings;
use Data::Dumper;
use FindBin qw|$Bin|;
use Cwd qw|abs_path|;
use Getopt::Lucid qw|:all|;
use Config::Std;
use Perl6::Slurp;
use Utils::Workflow;
use Utils::Hash qw|merge_conf|;

use lib "$Bin/../lib/";
use Shared qw|read_list|;

my @spec =  (
	Param("conf|c")->valid(sub { -r }),
	Param("list|l")->valid(sub { -r $_ }),
	Param("rename")->valid(sub { -r $_ }),
	Param("remove")->valid(sub { -r $_ }),
	Param("outdir|out"),
	Param("bed")->valid(sub { -r $_ || -d $_ }),
	Param("group")->valid(sub { -r $_ }),
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
	To extract subset of reads from BAM/CRAM files.

Usage:
	bam_subset.pl --conf Config --list BAMList --bed TARGET --outdir ROOTDIR

Options:
	--list: BAM or CRAM file list or directory. If a list file is provided, it should have 
			two columns per line: path to the BAM/CRAM file (cloud bucket supported) and sample ID. 
			If a directory is provided or the list only contain file name, sample ID will be taken
			from the basename of file after stripping of suffix.
	--rename: A rename list for BAM/CRAM files to match the IDs in the target file.
	--remove: A removal list for bad BAM/CRAM files to be removed. If rename is also provided,
				renaming will precede removal.
	--group: Group individuals into groups. It should have two columns: IID GID (group label). 
				If target directory is provided (below), then it will look for target file GID.bed for IID.
				If a sample does not belong to any group, it will fall back to look for IID.bed.
	--bed:  BED file or a directory containing BED files. If a file is provided, it contains the 
			region to subset all BAM/CRAM files. The format should be compatible with SAMtools or GATK. 
			If a directory is provided, it should contain individual or group specific target files in 
			BED format with name IID.bed or GID.bed. 
			And only BAM/CRAMs of individuals with target file will be written to output.
	
Notes: 
	The reads in selected target file will be extracted using either SAMtools (view) or GATK (PrintReads).
	SAMtools and GATK have similar functionality, except that SAMtools support AWS s3 and GATK support
	google cloud bucket. The output can be BAM or CRAM file, one file per individual.

	Sample ID in BAM/CRAM header may be different from the designated ID. As a convention of all pipeline
	tools, we do not modify BAM/CRAM header.

EOF
	exit 1;
}

$opt->validate({requires => [qw|conf list outdir|]});

my $rootdir = $opt->get_outdir;

my %conf = merge_conf($opt->get_conf, $opt->get_param); 


my $wkf = Utils::Workflow->new($rootdir, { engine => $opt->get_engine, force => $opt->get_force });

# Currently we skip testing for file if stage in from AWS
my %opt;
if (defined $conf{AWS}{STAGEIN}) {
	if (defined $conf{SSH}) {
		warn "SSH config will be ignored when AWS is specified!";
	}
	$opt{notest} = 1;
}
elsif (defined $conf{SSH}{HOST} && $conf{SSH}{HOST} !~ /localhost/i ) {
	if ($conf{SSH}{NOTEST} =~ /^Y|T/ ) {
		$opt{notest} = 1;
	}
}

my ($bam, $filetype) = read_list($opt->get_list,
	{ suffix => { bam => qr/\.bam$/, cram => qr/\.cram$/ }, 
	  rename => $opt->get_rename, remove => $opt->get_remove, %opt });

my %group;
if ($opt->get_group) {
	%group = map { (split)[0,1] } slurp $opt->get_group;
}

my ($target, %targets);
if ($opt->get_bed) {
	$conf{PATH}{TARGETS} = $opt->get_bed;
}
unless(defined $conf{PATH}{TARGETS}) {
	die "Must provide target files!";
}

my @samps;
if (-f $conf{PATH}{TARGETS}) {
	$target = abs_path($conf{PATH}{TARGETS});
	@samps = sort keys %$bam;
}
elsif (-d $conf{PATH}{TARGETS}) {
	foreach my $iid (sort keys %$bam) {
		# Fall back to IID if GID cannot be found
		my $gid = $group{$iid} // $iid;
		if (-f "$conf{PATH}{TARGETS}/$gid.bed") {
			$targets{$iid} = abs_path("$conf{PATH}{TARGETS}/$gid.bed");
		}
	}
	@samps = sort keys %targets;
}
else {
	die "Path to target is not a file or directory: $conf{PATH}{TARGETS}!";
}


printf STDERR "A total number of %d individual %s files are found, %d have associated target files\n",
	scalar(keys %$bam), $filetype, scalar(@samps);

exit 1 unless @samps > 0;

for(my $ii = 1; $ii <= @samps; $ii ++) {
	my $iid = $samps[$ii-1];
	open my $fout, ">$rootdir/par/IID2BAM.$ii" or die "Cannot write to IID2BAM.$ii";
	print $fout $iid, "\t", $bam->{$iid}, "\t", $targets{$iid} // $target, "\n";
}


$wkf->add(subset($conf{SUBSET}{TOOL}, $conf{SUBSET}{OUTFORMAT}), 
		 { name => "Subset".uc($conf{SUBSET}{OUTFORMAT}),  nslots => scalar(@samps),
		   expect => [ expect($conf{SUBSET}{OUTFORMAT}) ] });

$wkf->inst(\%conf);
$wkf->run({ conf => $conf{$wkf->{engine}}, dryrun => $opt->get_dryrun  });


sub subset {
	my ($tool, $outfmt) = @_;
	unless($outfmt eq 'bam' || $outfmt eq 'cram') {
		die "Cannot recognize output format: $outfmt";
	}
	unless($outfmt eq $conf{SUBSET}{OUTFORMAT}) {
		die "Incorrect output format: $outfmt <> $conf{SUBSET}{OUTFORMAT}";
	}
	my $script = "\nread IID BAMFILE TARGET < _PARDIR_/IID2BAM._INDEX_\n";
	if (lc($tool) eq 'samtools') {
		if ($conf{SSH}{HOST} && $conf{SSH}!~/localhost/i) {
			if ($outfmt eq 'cram') {
				$script .=<< 'EOF';

cat $TARGET | ssh _SSH.HOST_ \
	"_SSH.SAMTOOLS_ view _SUBSET.OPTION[ ]_ -M -h -T _PATH.FASTA_ -L - -C -o - $BAMFILE" > _OUTDIR_/$IID.cram 

samtools index _OUTDIR_/$IID.cram

EOF
			}
			elsif ($outfmt eq 'bam') {
				$script .= <<'EOF';

cat $TARGET | ssh _SSH.HOST_ \
	"_SSH.SAMTOOLS_ view _SUBSET.OPTION[ ]_ -M -h -T _PATH.FASTA_ -L - -b -o - $BAMFILE" > _OUTDIR_/$IID.bam

samtools index _OUTDIR_/$IID.bam

EOF
			}
		}
		else {
			if ($outfmt eq 'cram') {
				$script .=<< 'EOF';

samtools view _SUBSET.OPTION[ ]_ -M -h -T _PATH.FASTA_ -L $TARGET -C -o _OUTDIR_/$IID.cram $BAMFILE

samtools index _OUTDIR_/$IID.cram

EOF
			}
			elsif ($outfmt eq 'bam') {
				$script .= <<'EOF';

samtools view _SUBSET.OPTION[ ]_ -M -h -T _PATH.FASTA_ -L $TARGET -b -o _OUTDIR_/$IID.bam $BAMFILE

samtools index _OUTDIR_/$IID.bam

EOF
			}
		}	
	}
	elsif (lc($tool) eq 'gatk') {
		if ($conf{SSH}{HOST}) {
			die "We do not support extracting BAM/CRAMs via SSH using GATK!";
		}
		$script .= <<'EOF';

gatk --java-options "_RSRC.GATKOPT_" PrintReads \
	-I $BAMFILE -R _PATH.FASTA_ \
	-L $TARGET _SUBSET.OPTION[ \
	]_ -O _OUTDIR_/$IID._SUBSET.OUTFORMAT_

EOF
		if ($outfmt eq 'cram') {
			$script .= <<'EOF'

samtools index _OUTDIR_/$IID.cram

EOF
		}
		else {
			$script .= <<'EOF';

samtools index _OUTDIR_/$IID.bam

EOF
		}
	}

	if (exists $conf{AWS}{STAGEOUT}) {
		if ($outfmt eq 'cram') {
			$script .= <<'EOF';

stageout.pl --files $IID.cram,$IID.cram.crai --profile _AWS.PROFILEOUT_ --indir _OUTDIR_ --outdir _AWS.STAGEOUT_

EOF
		}
		elsif ($outfmt eq 'bam') {
			$script .= <<'EOF';

stageout.pl --files $IID.bam,$IID.bam.bai --profile _AWS.PROFILEOUT_ --indir _OUTDIR_ --outdir _AWS.STAGEOUT_

EOF
		}
	}
	return $script;
}

sub expect {
	my ($outfmt) = @_;
	my @exp;
	foreach my $iid (@samps) {
		if ($outfmt eq 'cram') {
			push @exp, ["out/$iid.cram", "out/$iid.cram.crai"];
		}
		else {
			push @exp, ["out/$iid.bam", "out/$iid.bam.bai"];
		}
	}
	return @exp;
}
