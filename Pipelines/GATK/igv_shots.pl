#!/bin/env/perl

use strict;
use warnings;
use Carp;
use IO::Dir;
use IO::File;
use FindBin qw|$Bin|;
use List::Util qw|max first|;
use List::MoreUtils qw|all|;
use Config::Std;
use Data::Dumper;
use Getopt::Lucid qw|:all|;
use Perl6::Slurp;
use String::ShellQuote;
use File::Path qw|make_path remove_tree|;
use Utils::Hash qw|merge_conf|;
use Utils::File::Iter;
use Utils::Workflow;
use Genome::Ranges;
use Genome::Ranges::IntSet;

use lib "$Bin/../lib";
use Shared qw|read_list|;


############################
## Command line interface ##
############################

my @spec = (
	Param("conf|c")->valid(sub { -r }),
	Param("list|l")->valid(sub { -r }),
	Param("tab")->valid(sub { -r }),
	Param("ped")->valid(sub { -r }),
	Param("rename")->valid(sub { -r }),
	Param("outdir|out"),
	Keypair("param|par"),
	Switch("local"),
	Switch("help|h"),
	Switch("dryrun|dry"),
	Switch("force")
	);

my $opt = Getopt::Lucid->getopt(\@spec);

if ($opt->get_help) {
	print STDERR <<EOF;
Purpose:
	This is a pipeline script that runs on Desktop to create IGV screenshots around input 
	variants.

Usage:
	igv_shots.pl --conf Config --list BAMDir --outdir RootDir

Options:
	--list: BAM/CRAM file list or directory. The list should have 1 or 2 columns per line: path to 
			BAM/CRAM file (required) and sample ID. If only one column is available or if a directory 
			is used, sample ID will be taken from file basename after removing suffix.
	--tab:  The input variant table. It overrides variant table in config.
	--ped:  The pedigree file used to look for familial relationships. It is recommended to provide
			most complete pedigree file so that all pairwise relationships can be determined.
			It overrides the PED file in config.
	--rename: A rename list. It can be used to rename samples from BAM/CRAM list to match samples in the variant table.
	--local: (Advanced) Extract and store BAM files around candidate variants locally before running IGV.
			This option only works for small number of variants per individual. When this option is swithed on
			The input must BAM/CRAM list with path on the remote server, and it cannot be directory.

Notes:
	To create IGV screeshots, the IGV GUI must be able to run. We typically run this script on Desktop.
	If BAM/CRAM files are stored on remote server, we can create a directory to store symlinks to orignal
	BAM/CRAM files and mount that directory to the desktop by ssh-fs. When PED file is provided, IGV 
	screenshots will also be generated for all available family members.

Dependencies:
	IGV, ImageMagick.
	The utility make_IGV_snapshots.py requires igv.jar from IGV binary distribution. This is available
	for all MacOS App. For linux, only v2.4 or before had binary distribution. Also need to ensure the 
	command line jave version is compatible with the IGV version.


EOF
	exit 1;
}

$opt->validate({requires => [qw|conf outdir|]});

my $rootdir = $opt->get_outdir;

my %conf = merge_conf($opt->get_conf, $opt->get_param); 
$conf{PATH}{UTILDIR} = shell_quote("$Bin/utils");
$conf{PATH}{MODULE} = shell_quote("$Bin/../utils");

if ($opt->get_tab) {
	$conf{Input}{File} = $opt->get_tab;
}
croak "Cannot find input table file: $conf{Input}{File}" unless -f $conf{Input}{File};

my %remote;
if ($opt->get_local) {
	$remote{remote} = $conf{REMOTE}{HOST};
	if ($conf{REMOTE}{NOTEST} =~ /^Y|T/ ) {
		$remote{notest} = 1;
	}
}

my ($bams, $filetype); 
if ($opt->get_list) {
	($bams, $filetype) = read_list($opt->get_list, 
		{ suffix => ['bam', 'cram'], %remote, rename => $opt->get_rename });
}
else {
	if (exists $conf{Input}{BAMLIST} && (-f $conf{Input}{BAMLIST} || -d $conf{Input}{BAMLIST})) {
		($bams, $filetype) = read_list($conf{Input}{BAMLIST}, 
			{ suffix => 'bam', %remote, rename => $opt->get_rename });
	} 
	elsif (exists $conf{Input}{CRAMLIST} && (-f $conf{Input}{CRAMLIST} || -d $conf{Input}{CRAMLIST})) {
		($bams, $filetype) = read_list($conf{Input}{CRAMLIST}, 
			{ suffix => 'cram', %remote, rename => $opt->get_rename });
	}
	else {
		die "Must provide BAMLIST or CRAMLIST in the config!";
	}
}
#my ($bams, $filetype) = read_list($opt->get_list, 
#	{ suffix => ['bam', 'cram'], %remote, rename => $opt->get_rename });
print STDERR "Finished reading BAM/CRAM file list\n";

my (%fids, %famsamps);
if ($opt->get_ped) {
	$conf{PED}{FILE} = $opt->get_ped;
}
if (exists $conf{PED}) {
	$conf{PED}{OPTION} = "";
	unless(-f $conf{PED}{FILE}) {
		die "Cannot find PED file $conf{PED}{FILE}";
	}
	%fids = map { (split)[1,0] } slurp $conf{PED}{FILE};
	#($famsamps, $famrels) = fam_rels($opt->get_ped)
	while(my ($iid, $fid) = each %fids) {
		push @{$famsamps{$fid}}, $iid;
	}
	if (defined $conf{PED}{IGNORE}) {
		$conf{PED}{OPTION} .= "--ped-ignore $conf{PED}{IGNORE} ";
	}
	if (defined $conf{PED}{TWINS}) {
		$conf{PED}{OPTION} .= "--ped-twins $conf{PED}{TWINS} ";
	}
}


#########################################################
# Parse input file to prepare sample-specific intervals #
#########################################################

my $wkf = Utils::Workflow->new($rootdir,
	{ engine => 'BASH', force => $opt->get_force, strict_var => 1 });

my (@samps, %sampvars);
{
	# parse the input file to determine the ranges to extract from cases
	# If pedigree file is givem, the same range will be added to all relatives
	my $input_header = 1;
	my @input_fields = split(',', $conf{Input}{Fields});
	croak "Must provide 5 fields IID,CHR,POS,REF,ALT" unless @input_fields == 5;
	if (all { /^\d+$/ } @input_fields) {
		$input_header = 0;
	}
	#croak "Input fields must be renamed to IID,Chrom,Position,Ref,Alt"
	#	unless join(',', values %$finput) eq 'IID,Chrom,Position,Ref,Alt';
	my ($it, $fnames) = iter_file($conf{Input}{File}, { header => $input_header, fsep => qr/\t/ });
	foreach my $field (@input_fields) {
		unless(grep { $field eq $_ } @$fnames) {
			croak "Cannot find field $field in the input file";
		}
	}

	my (%bamrngs);
	open my $flst, ">$rootdir/par/var_list.txt" or die "Cannot write to variant list";
	while(my $dat = $it->()) {
		my ($iid, $chr, $pos, $ref, $alt) = @{$dat}{@input_fields};
		$pos =~ s/,//g;
		my $varid = join("_", $chr, $pos, $ref, $alt);
		print $flst $varid, "\t", $iid, "\n";
		if (defined $fids{$iid}) {
			foreach my $sampid (@{$famsamps{$fids{$iid}}}) {
				push @{$sampvars{$sampid}}, [$chr, $pos-$conf{Param}{VarPad}, 
										  	 $pos+length($ref)-1+$conf{Param}{VarPad}, $varid];
				if ($opt->get_local) {
					push @{$bamrngs{$sampid}}, [$chr, $pos-$conf{Param}{BamPad}, 
										 		$pos+length($ref)-1+$conf{Param}{BamPad}];
				}
			}
		}
		else {
			push @{$sampvars{$iid}}, [$chr, $pos-$conf{Param}{VarPad}, 
								   		$pos+length($ref)-1+$conf{Param}{VarPad}, $varid];
			if ($opt->get_local) {
				push @{$bamrngs{$iid}}, [$chr, $pos-$conf{Param}{BamPad}, 
										 $pos+length($ref)-1+$conf{Param}{BamPad}];
			}
		}
	}

	my $nsamp = scalar(keys %sampvars);
	my $nincl = grep { defined $bams->{$_} } keys %sampvars; 
	print STDERR "A total of $nsamp samples (including relatives) found in the variants table, ",
		"$nincl of them were associated with $filetype files\n";

	unless ($nsamp > 0) {
		print STDERR "No sample found in the variant table with BAM/CRAM files\n";
		exit 1;
	}

	# Samples with both variants and BAM files
	@samps = sort grep { defined $bams->{$_} } keys %sampvars;

	my $it2 = iter_file($conf{Input}{File}, { header => $input_header, fsep => qr/\t/ });

	
	# Write the list of variants that needs to be processed for each sample
	for (my $ii = 0; $ii < @samps; $ii ++) {
		my $iid = $samps[$ii];
		make_path "$rootdir/wrk/$iid";

		my $jj = $ii + 1;
		open my $fout, ">$rootdir/par/IID2BAM.$jj" or die "Cannot write to IID2BAM";
		print $fout $iid, "\t", $bams->{$iid}, "\n";	

		open $fout, ">$rootdir/wrk/$iid/vars.bed" or die "Cannot write vars.bed for $iid";
		foreach my $var (@{$sampvars{$iid}}) {
			print $fout join("\t", @$var), "\n";
		}

		if ($opt->get_local) {
			my $unionrng = Genome::Ranges::IntSet->new($bamrngs{$iid})->to_ranges();
			$unionrng->write("$rootdir/wrk/$iid/bam_regions.bed", {bed => 1});
			#my @regions;
			#my $iter = $unionrng->iter();
			#while(my $dat = $iter->()) {
			#	push @regions, sprintf("%s:%d-%d", @$dat);
			#}
			#open my $frng, ">$rootdir/wrk/$iid/bam_regions.txt" or die "Cannot write to bam_regions.txt";
			#print $frng join(" ", @regions), "\n";
		}	
	}
}

write_config %conf, "$rootdir/par/run.conf" unless $opt->get_dryrun;

#####################################
##  Workflow initialization & Run  ##
#####################################

my %deparg;
if ($opt->get_local) {
	$wkf->add(extract_bams(), { name => "ExtractBAMs", 
			expect => [ map { ["out/$_.cram", "out/$_.cram.crai"] } @samps ],
			nslots => scalar(@samps) });
	$deparg{depend} = "ExtractBAMs";
}
$wkf->add(IGV_snapshots(),	{ name => "IGVSnapShots", %deparg, 
		expect => [ get_snapshot_exp() ], nslots => scalar(@samps)  })
	->add(organize_plots(), { name => "OrganizePlots", depend => "IGVSnapShots",
	 	expect => [ get_final_exp() ] });


$wkf->inst(\%conf);
$wkf->run({ conf => $conf{$wkf->{engine}}, dryrun => $opt->get_dryrun  });


############################
## Workflow components    ##
############################

##
# Extract BAM files ##
##
# We indeed only store CRAM files locally
sub extract_bams {
	my $script;
	if ($conf{REMOTE}{HOST} =~ /localhost/i) {
		$script = <<'EOF'; 

read IID BAM < _PARDIR_/IID2BAM._INDEX_

_REMOTE.SAMTOOLS_ view -T _REMOTE.FASTA_ -Ch -L _WRKDIR_/$IID/bam_regions.bed -M \
	-o _OUTDIR_/$IID.cram $BAM

samtools index _OUTDIR_/$IID.cram

EOF
	}
	else {
		$script = <<'EOF';

read IID BAM < _PARDIR_/IID2BAM._INDEX_

cat _WRKDIR_/$IID/bam_regions.bed | \
	ssh _REMOTE.HOST_ \
	"_REMOTE.SAMTOOLS_ view -T _REMOTE.FASTA_ -Ch -L - -M -o - $BAM" \
	> _OUTDIR_/$IID.cram

samtools index _OUTDIR_/$IID.cram

EOF
	}
	return $script;
}

# Cannot do parallel on this
sub IGV_snapshots {
	my $script = <<'EOF';

read IID BAM < _PARDIR_/IID2BAM._INDEX_

#if [[ -f _OUTDIR_/$IID.cram ]]; then 
#	BAM=_OUTDIR_/$IID.cram
#fi

# To process alls sample 
# find _PARDIR_ -name 'IID2BAM.*' | while read IID2BAMFILE; do
#	IID=$(cat $IID2BAMFILE | cut -f1)
#	BAM=$(cat $IID2BAMFILE | cut -f2)

cd _WRKDIR_/$IID

if [[ -f ../../out/$IID.cram ]]; then 
	BAM=../../out/$IID.cram
fi

python _PATH.UTILDIR_/make_IGV_snapshots.py $BAM \
	-g _IGV.GENOME_ \
	-ht _IGV.HEIGHT_ \
	-r vars.bed  -o . \
	-bin _IGV.PATH_ \
	-mem _IGV.MEM_ \
	-nf4

#done

EOF
	return $script;
}

sub get_snapshot_exp {
	my @exp = map { [ "wrk/$_/IGV_snapshots.bat" ] } @samps;
	for(my $ii = 0; $ii < @samps; $ii ++) {
		my $iid = $samps[$ii];
		push @{$exp[$ii]}, map { "wrk/$iid/$_->[3].png" } @{$sampvars{$iid}};
	}
	return @exp;
}

sub organize_plots {
	my $script;
	if (exists $conf{PED}) {
		$script = <<'EOF';

perl _PATH.MODULE_/organize_plots.pl --rootdir _WRKDIR_/.. --suffix png \
	--ped _PED.FILE_ _PED.OPTION_ --montage --caption

EOF
	}
	else {
		$script = <<'EOF';

perl _PATH.MODULE_/organize_plots.pl --rootdir _WRKDIR_/.. --suffix png

EOF
	}
	return $script;
}

sub get_final_exp {
	my @exp;
	open my $fin, "$rootdir/par/var_list.txt" or die "Cannot read varlist";
	while(<$fin>) {
		my ($varid, $iid) = split;	
		push @exp, "byiid/$iid-$varid.png";
		push @exp, "byvar/$varid-$iid.png";
		if (exists $conf{PED}) {
			push @exp, "wrk/$iid/$varid-montage.png";
		}
	}
	if (exists $conf{PED}) {
		foreach my $iid (@samps) {
			foreach my $var (@{$sampvars{$iid}}) {
				push @exp, "wrk/$iid/$var->[3]-caption.png";
			}
		}
	}
	return @exp;
}


