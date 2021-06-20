#!/usr/bin/env perl

use strict;
use warnings;
use Carp;
use IO::File;
use Data::Dumper;
use Config::Std;
use FindBin qw|$Bin|;
use File::Copy qw|copy|;
use File::Path qw|make_path|;
use File::Temp qw|tempdir|;
use List::MoreUtils qw|all uniq any none|;
use String::ShellQuote;
use Hash::Util qw|lock_keys|;
use Getopt::Lucid qw|:all|;
use Utils::Number qw|pretty_bps|;
use Utils::List qw|insert_after insert_before|;
use Utils::File qw|count_line|;
use Utils::Hash qw|str2hash peek_hash|;
use Utils::File::Iter qw|iter_file|;


use lib "$Bin/../lib/";
use Shared qw|parse_fstr merge_tabfiles slurp_xref read_geneset struct_dat|;
use Variants qw|cnv_type cnv_id|;

############################
## Command line interface ##
############################

my @spec =  (
	Param("conf|c")->valid(sub { -r }),
	Param("input|in|i")->valid(sub { -r }),
	Param("wrkdir|wrk|w"),
	Param("output|out|o"),
	Switch("force"),
	Switch("help|h")
	);

my $opt = Getopt::Lucid->getopt(\@spec);

if ($opt->get_help) {
	print STDERR <<EOF;
Purpose:
	Annotate copy number variations (CNVs) using VEP and custom tools.

Usage:
	anno_cnvs.pl --conf CONFIG [--in INPUT] --wrk WRKDIR -out OUTPUT

Options:
	--in: Input CNV table. It overrides the CNV table in the config file.
	--wrk: Workding directory for storing intermediate files.
	--force: Force running the program even if the output exist. The default is to exit if output file already exists.

Notes:
	The script collects the following types of annotations for each CNV:

	1. Functional effects to genes & transcripts
	   Gene-based annotation of deletion and duplications is performed by VEP. Note that VEP also support annotation of 
	   insertion(INS) or tandem duplication(TDUP), but they have not been incoporated in this script. The default gene 
	   source for annotation is GENCODE Basic set of Ensembl genes. Functional effects to genes and transcripts are 
	   simplified to collapse most frequent annotation terms (Details are given in config template).
	  
	2. Highlighting genes in selected gene sets
	   It is also possible to gather additional gene level information for each gene. But because large CNVs can overlap 
	   many genes, it is more informative to define gene sets from gene level information, and pickup important genes 
	   overlapping the CNV. There are multiple ways to define gene sets. GeneID is used to link genes with gene set. 

	3. Overlapping genomic intervals
	   Genomic interval-based annotations are performed by overlap_cnv.pl. It provides flexible ways to define overlapping
	   proportions and functions to summarize the overlapping features. In addition to genomic features, CNVs in population 
	   cohorts and morbidity map can also be used as annotation sources. When CNVs are used as an annotation source, it can 
	   further require that overlapping CNVs match the type. 

	Some standalone utilities can be used to add further annotations to the output:
	 	highlight_cnvgenes.pl - To highlight selected genes overlapping CNV region. 
	 							It extends the function of this script to allow use gene names or other gene IDs in gene set,
	 							display extra gene level information together with highlighted genes, and use overlapping genes
	 							from selected gene table of UCSC database. 
	 	overlap_cnvs.pl - This is also used by interval-based annotation of this script.

Input/output:
	Input CNV table should be tab separated text file or comma separated csv file. 
	Output will have one line for each CNV. Default fields are explained in config template. Additional output fields
	can be customized by different options.

EOF
	exit 1;
}

$opt->validate({ requires => [qw|conf output|] });

my $infile = $opt->get_input;
my $outfile = $opt->get_output;
my $config = $opt->get_conf;

if (-f $outfile) {
	unless($opt->get_force) {
		print STDERR "Output file $outfile already exists!\n";
		exit 1;
	}
}

read_config $config => my %config;


my $wrkdir = $opt->get_wrkdir;
if ($wrkdir) {
	make_path $wrkdir unless -d $wrkdir;
}
else {
	$wrkdir = tempdir(CLEANUP => 1);
}

# Annotate all inputs from a directory
if (defined $infile && -d $infile) {
	my @varlists = grep { !/^\./ && -f "$infile/$_" } IO::Dir->new($infile)->read();
	my $option = "";
	if (exists $config{Parallel} && defined $config{Parallel}{jobs}) {
		$option = join(" ", map { "--$_ $config{Parallel}{$_}" } keys %{$config{Parallel}});
	}
	else {
		die "Must provide number of threads in config file for parallel";
	}
	my $command = qq(find $infile/ -type f | parallel --eta $option "perl $0 --in {} --wrk $wrkdir/{/}.wrkdir --conf $config --out $wrkdir/{/} 1>$wrkdir/{/}.out 2>$wrkdir/{/}.err");
	#print $command, "\n"; exit 1;
	system($command);
	# Then merge the output
	open my $fout, ">$outfile" or die "Cannot write to $outfile";
	for(my $ii = 0; $ii < @varlists; $ii ++) {
		unless(-f "$wrkdir/$varlists[$ii]") {
			die "Cannot find split $ii output $varlists[$ii] in the working directory";
		}
		open my $fin, $wrkdir."/".$varlists[$ii] or die "Cannot open $varlists[$ii]";
		<$fin> if $ii > 0;
		while(<$fin>) {
			print $fout $_;
		}
	}
	exit 0;
}

######################################
### Parsing config & prepare inputs ##
######################################

if (defined $infile && -f $infile) {
	$config{Input}{File} = $infile;
	unless (-f $config{Input}{File}) {
		croak "Cannot find input file: $config{Input}{File}";
	}
}
lock_keys(%config);

while (my ($label, $file) = each %{$config{Path}}) {
	if ($label =~ /DB$/) {
		croak "Cannot find $label directory: $file" unless -d $file;
	}
	else {
		croak "Cannot find $label file: $file" unless -f $file;
	}
}

croak "Cannot find input file" unless -f $config{Input}{File}; 

my $input_header = 1;
# Parse input CNV file
my %known;
{
	my $finput = parse_fstr($config{Input}{Fields}, 1);
	my @fvar = keys %$finput;
	if (all { /^\d+$/ } @fvar) {
		$input_header = 0;
	}
	croak "Input fields must be renamed to Chrom,Start,End,Type" 
		unless join(',', values %$finput) eq 'Chrom,Start,End,Type';
	
	my $it = iter_file($config{Input}{File}, { header => $input_header, fsep => qr/\t/, select => $finput});
	open my $fout, ">$wrkdir/input.txt" or die "Cannot write to VEP input";
	while(my $dat = $it->()) {
		$dat->{Start} =~ s/,//g; $dat->{End} =~ s/,//g;
		unless($dat->{Start} =~ /^\d+$/ && $dat->{End} =~ /^\d+/ && $dat->{Start} < $dat->{End}) {
			warn "Invalid range $dat->{Chrom}:$dat->{Start}-$dat->{End}";
			next;
		}
		my $cntype = cnv_type($dat->{Type});
		if ($cntype eq 'Normal') {
			warn "Cannot recognize CNV type: $dat->{Type}";
			next;
		}
		my $inputline = join("\t", @{$dat}{qw|Chrom Start End|}, uc($cntype));
		next if defined $known{$inputline};
		print $fout $inputline, "\n";
		$known{$inputline} = 1;
	}
}

#########################################
### Run VEP for gene level annotation ###
#########################################
my $utildir = shell_quote("$Bin/utils");
{
	my ($opt_a, $opt_c, $ver_a);
	if ($config{Input}{HG} eq 'hg19' || $config{Input}{HG} eq 'b37') {
		$ver_a = "GRCh37";
		$opt_a = "-a $ver_a --port 3337";
		if ($config{Input}{HG} eq 'b37') {
			$config{Input}{HG} = "hg19";
		}
	}
	elsif ($config{Input}{HG} eq 'hg38' || $config{Input}{HG} eq 'b38' ) {
		$ver_a = "GRCh38";
		$opt_a = "-a $ver_a";
		if ($config{Input}{HG} eq 'b38') {
			$config{Input}{HG} = "hg38";
		}
	}
	else {
		croak "Unsupported hg build: $config{Input}{HG}";
	}
	
	if (exists $config{VEP} && exists $config{VEP}{DBVer}) {
		$opt_c = "--cache_version $config{VEP}{DBVer}";
		# Check VEP Cache for the specific gene version
		unless(-d "$config{Path}{VEPDB}/homo_sapiens/$config{VEP}{DBVer}_$ver_a") {
			croak "Cache version $config{VEP}{DBVer} does not exist!";
		}
	}
	else {
		# If version is not provided, will assume to use the default from cache
		# that matches the VEP version
		$opt_c = "";
	}
	print "Running VEP for gene level annotations\n";
	my $vepcmd = "$config{Path}{VEP} -i $wrkdir/input.txt -o $wrkdir/vep_out.txt --force --no_stats ".
		"$opt_a --cache --dir_cache $config{Path}{VEPDB} $opt_c ";
	if (defined $config{VEP}{Option} && $config{VEP}{Option} !~ /^\s*$/) {
		$vepcmd =  $config{VEP}{Option}." ";
		$vepcmd =~ s/^["'']//; $vepcmd =~ s/["'']$//;
	}
	else {
		$vepcmd .= "--gencode_basic --symbol --transcript_version --biotype --numbers ";
	}
	print $vepcmd, "\n";
	system($vepcmd);

	# Validate the output: should have same number of variants as input
	my %vars;
	open my $fin, "$wrkdir/vep_out.txt" or die "Cannot open $wrkdir/vep_out.txt";
	while(<$fin>) {
		next if /^#/;
		my ($region, $type) = (split)[1,2];
		$vars{$region,$type} = 1;
	}
	unless (scalar(keys %vars) == scalar(keys %known)) {
		print STDERR "Input=", scalar(keys %known), ", VEPOutput=", scalar(keys %vars), "\n";
		die "Incorrect number of variants in $wrkdir/vep_out.txt";
	}

	# Collect VEP results into a table
	my ($inclusion, $exclusion) = ("", "");
	if (exists $config{Gene} && exists $config{Gene}{Include}) {
		croak "Cannot find gene inclusion list file" unless -f $config{Gene}{Include};
		$inclusion = "--include-gene $config{Gene}{Include} ";
	}
	if (exists $config{Gene} && exists $config{Gene}{Exclude}) {
		croak "Cannot find gene exclusion list file" unless -f $config{Gene}{Exclude};
		$exclusion .= "--exclude-gene $config{Gene}{Exclude} ";
	}
	if (exists $config{Transcript} && exists $config{Transcript}{Include}) {
		croak "Cannot find transcript inclusion list file" unless -f $config{Transcript}{Include};
		$inclusion .= "--include-trans $config{Transcript}{Include} ";
	}
	if (exists $config{Transcript} && exists $config{Transcript}{Exclude}) {
		croak "Cannot find transcript exclusion list file" unless -f $config{Transcript}{Exclude};
		$exclusion .= "--exclude-trans $config{Transcript}{Exclude} ";
	}
	my $coding = "";
	if ($config{VEP}{Coding} =~ /^Y|T/i) {
		$coding = "--coding";
	}
	my $cutoff = "";
	if (defined $config{VEP}{RankCutoff}) {
		$cutoff = "--cutoff $config{VEP}{RankCutoff}";
	}

	print "Collect VEP output into tabular format\n";
	my $sumcmd = "perl $utildir/collect_vepout_cnv.pl -in $wrkdir/vep_out.txt ".
		"$exclusion $coding $cutoff -out $wrkdir/vep_out_tab.txt --no-overlap ";
	if (exists $config{Transcript} && exists $config{Transcript}{Full} && $config{Transcript}{Full} =~ /^Y|T/i) {
		$sumcmd .= "--full-trans"
	}
	print $sumcmd, "\n";
	system($sumcmd);

	unless(-f "$wrkdir/vep_out_tab.txt.done") {
		die "collect_vepout_cnv did not finish successfully";
	}
}

######################################################
### Add gene level information and highlight genes ###
######################################################

my $vepout;
if (exists $config{Gene} && exists $config{Gene}{GXref}) {
	my ($it, $fnames) = iter_file("$wrkdir/vep_out_tab.txt", { fsep => qr/\t/ });
	my @gfields = @$fnames;

	my ($xgdat, $xgfds);
	if (exists $config{Gene}{GXref}) {
		($xgdat, $xgfds) = slurp_xref($config{Gene}{GXref}, $config{Gene}{GXref_Fields});
		foreach my $xfd (@$xgfds) {
			if (any { $xfd eq $_ } @$fnames) {
				croak "Extra geneinfo field $xfd already exist in the vep_out_tab.txt!";
			}
		}
		# Identify GeneID field and splice in additional gene related fields
		my ($gid_ii) = grep { $gfields[$_] eq 'GeneIDs' } 0..$#gfields;
		die "Cannot find GeneID column" unless defined $gid_ii;
		splice @gfields, $gid_ii+1, 0, @$xgfds;
	}
	# Reformat outputs
	my $fout = IO::File->new("$wrkdir/vep_out_tab_ginfo.txt", "w");
	print $fout join("\t", @gfields), "\n";
	while(my $dat = $it->()) {
		my @gids = split(';', $dat->{GeneIDs});
		foreach my $xgfd (@$xgfds) {
			$dat->{$xgfd} = join(";", map { $xgdat->{$_}{$xgfd} // "." } @gids);
		}
		print $fout join("\t",  map { $_ // "." } @{$dat}{@gfields}), "\n";
	}
	$vepout = "$wrkdir/vep_out_tab_ginfo.txt";
}

# Highlight genes
if (exists $config{Gene}{Set}) {
	# First read gene sets
	# $geneset is a hash ref: $geneset{Set} = [Genes ...]
	# $gsmember of hash of hash ref: $gsmember{Gene}{Set} = 1
	my ($geneset, $gsmember) = read_geneset($config{Gene}{Set}, $config{Gene}{Set_Fields});

	# Determine output columns and ther positions
	my @hlfields;
	my $superset;
	my %setgroup;
	if (exists $config{Gene}{SuperSet}) {
		if (ref $config{Gene}{SuperSet} eq 'ARRAY') {
			$config{Gene}{SuperSet} = join(',', @{$config{Gene}{SuperSet}});
		}
		$superset = str2hash($config{Gene}{SuperSet}, {psep => ',', kvsep => ':', order => 1});

		while(my ($set, $field) = each %$superset) {
			push @{$setgroup{$field}}, $set;
			unless (grep { $field eq $_ } @hlfields) {
				push @hlfields, $field;
			}
		}
	}
	else {
		$superset = {};
		foreach my $set (keys %$geneset) {
			$superset->{$set} = $set;
			push @{$setgroup{$set}}, $set;
			push @hlfields, $set;
		}
	}

	my ($it, $fnames) = iter_file($vepout, { fsep => qr/\t/ });
	my @gfields = @$fnames;
	my $highlight;
	if (exists $config{Gene}{Highlight}) {
		$highlight = $config{Gene}{Highlight};
	}
	else {
		$highlight = "Symbol";
	}
	unless (grep { $_ eq $highlight } @$fnames) {
		die "Cannot find field $highlight in $vepout";
	}
	foreach my $field (@hlfields) {
		if(grep { $field eq $_  } @$fnames) {
			die "Gene higlight field $field already exist in $vepout";
		}
	}

	insert_after(\@gfields, $highlight, @hlfields);

	my $outfile = $vepout; $outfile =~ s/\.txt$/_highlight.txt/;
	open my $fout, ">$outfile" or die "Cannot write to $outfile";
	print $fout join("\t", @gfields), "\n";
	while(my $dat = $it->()) {
		my $info = struct_dat($dat, {key => 'GeneIDs', value => $highlight, sep => ';', clone => 1});
		foreach my $gid (split(';', $dat->{GeneIDs})) {
			while(my ($set, $field) = each %$superset) {
				if (defined $gsmember->{$gid}{$set}) {
					if (@{$setgroup{$field}} > 1) {
						push @{$info->{$field}{$info->{$highlight}{$gid}}} => $set;
					}
					else {
						push @{$info->{$field}}, $info->{$highlight}{$gid};
					}
				}
			}
		}
		foreach my $field (@hlfields) {
			unless (defined $info->{$field}) {
				$dat->{$field} = ".";
			}
			elsif (ref $info->{$field} eq 'HASH') {
				$dat->{$field} = join(",", map {  $_.'('.join('|', @{$info->{$field}{$_}}).')' } 
										sort keys %{$info->{$field}});
			}
			elsif (ref $info->{$field} eq 'ARRAY') {
				$dat->{$field} = join(",", sort @{$info->{$field}});
			}
			else {
				die "Cannot recognize the structure of $field";
			}
		}
		print $fout join("\t", @{$dat}{@gfields}), "\n";
	}
	$vepout = $outfile;
}

##########################################
### Find overlapping features or CNVs  ###
###########################################

if (exists $config{Overlap}{DBFile}) {
	if (all { ref $config{Overlap}{"DBFile$_"} eq '' } ("", qw|_Fields _Operation _Overlap|)) {
		foreach ("", qw|_Fields _Operation _Overlap|) {
			$config{Overlap}{"DBFile$_"} = [$config{Overlap}{"DBFile$_"}];
		}
	}
	else {
		if (any { ref $config{Overlap}{"DBFile$_"} eq '' } ("", qw|_Fields _Operation _Overlap|)) {
			die "Not all DBFile specifications can be found";
		}
		else {
			unless(all { @{$config{Overlap}{"DBFile$_"}} == @{$config{Overlap}{DBFile}} } qw|_Fields _Operation _Overlap|) {
				die "Not all DBFile specifications can be found";
			}
		}
	}
	my $prefix = $vepout; $prefix =~ s/\.txt$//;
	for(my $ii = 0; $ii < @{$config{Overlap}{DBFile}}; $ii ++) {
		my $dbfile = $config{Overlap}{DBFile}[$ii];
		unless(-f $dbfile) {
			die "Cannot find DBFile: $dbfile";
		}
		my $dbfields = $config{Overlap}{DBFile_Fields}[$ii];
		$dbfields =~ s/([)(])/\\$1/g;
		my $dbops = $config{Overlap}{DBFile_Operation}[$ii];
		$dbops =~ s/([)(])/\\$1/g;
		my %overlap = str2hash($config{Overlap}{DBFile_Overlap}[$ii], { psep => qr/[,;]/, kvsep => qr/[:=]/ });
		unless(join(",", sort keys %overlap) eq 'Probes,Query,Target') {
			die "Incorrect specification for overlap: $config{Overlap}{DBFile_Overlap}[$ii]";
		}
		my $opt = "";
		if ($overlap{Probes} =~ /^Y/i) {
			unless(exists $config{Overlap}{ProbeFile}) {
				die "Must specify a probes file to enable Probe option in DBFile_Overlap";
			}
			unless(-f $config{Overlap}{ProbeFile}) {
				die "Cannot find probes file: $config{Overlap}{ProbeFile}";
			}
			$opt = "--probes $config{Overlap}{ProbeFile}";
		}

		my $input = $ii > 0 ? sprintf("%s.%d.txt", $prefix, $ii) : "$prefix.txt";
		$vepout = sprintf("%s.%d.txt", $prefix, $ii+1); 

		my $overlapcmd = "perl $utildir/overlap_cnvs.pl -in $input -out $vepout -db $dbfile -db-fields $dbfields ".
			"-op $dbops --query-overlap $overlap{Query} --targ-overlap $overlap{Target} $opt";
		print $overlapcmd, "\n";
		system($overlapcmd);
		unless(-f $vepout) {
			die "Cannot find CNV overlap annotation output: $vepout";
		}
		# Check input and output have the same number of lines
		my $inline = count_line($input);
		my $outline = count_line($vepout);
		unless($inline == $outline) {
			die "Input and output of CNV overlap annotation do not have the same number of lines: $input($inline), $vepout($outline)";
		}
	}
}

##################################
### Merge and reformat results ###
##################################

{
	my $infstr = $config{Input}{Fields};
	# Key to IID
	my $f_iid;
	if (exists $config{Input}{XtraFields}) {
		$infstr .= ",".$config{Input}{XtraFields};
	}
	if (exists $config{Sample} && exists $config{Sample}{IDField}) {
		my $ftmp = parse_fstr($config{Sample}{IDField}, 1);
		$f_iid = (values %$ftmp)[0];
		$infstr = $config{Sample}{IDField}.",".$infstr;
	}
	
	my $finput = parse_fstr($infstr, 1);
	my @ifields = values %$finput;
	# Splice in VarID column into appropriate positions
	#my $chrii = first { $ifields[$_] eq 'Chrom' } 0..$#ifields;
	#splice @ifields, $chrii, 0, 'VarID';
	insert_before(\@ifields, 'Chrom', 'VarID');
	insert_after(\@ifields, 'Type', 'Length', 'Size');

	my ($it, $fnames) = iter_file("$config{Input}{File}",
								{ header => $input_header, fsep => qr/\t/, 
								  select => $finput });

	my ($xdat, $xfds);
	if (exists $config{Sample} && exists $config{Sample}{IXref}) {
		($xdat, $xfds) = slurp_xref($config{Sample}{IXref}, $config{Sample}{IXref_Fields});		
		foreach my $xfd (@$xfds) {
			if (any { $xfd eq $_ } @$fnames) {
				croak "Field name $xfd already exist in input file"
			}
		}
		# Extra sample level info will appear at the beginning after IID
		splice @ifields, 1, 0, @$xfds;
	}

	my $fout = IO::File->new("$wrkdir/input_full.txt", "w");
	print $fout join("\t", @ifields), "\n";
	while(my $dat = $it->()) {
		$dat->{Start} =~ s/,//g; $dat->{End} =~ s/,//g;
		unless($dat->{Start} =~ /^\d+$/ && $dat->{End} =~ /^\d+/ && $dat->{Start} < $dat->{End}) {
			next;
		}
		$dat->{Type} = cnv_type($dat->{Type});
		next if $dat->{Type} eq 'Normal';

		$dat->{VarID} = cnv_id(@{$dat}{qw|Chrom Start End Type|});
		$dat->{Length} = $dat->{End}-$dat->{Start}+1;
		$dat->{Size} = pretty_bps($dat->{End}-$dat->{Start}+1);

		if (defined $xdat && defined $xdat->{$dat->{$f_iid}}) {
			@{$dat}{@$xfds} = @{$xdat->{$dat->{$f_iid}}}{@$xfds};
		}
		print $fout join("\t", map { $_ // "." } @{$dat}{@ifields}), "\n"
	}
	close $fout;

	my ($it2, $fnames2) = iter_file($vepout, {fsep => qr/\t/});
	my @igkeys;
	foreach my $gfield (@$fnames2) {
		if (grep { $gfield eq $_ } @ifields) {
			push @igkeys, $gfield;
		}
	}
	unless(join(',', @igkeys) eq 'VarID,Chrom,Start,End,Type') {
		croak "Incorrect keys fields for merging input and vep output";
	}
	print "Merge with selected fields from input file\n";
	my $igkeys = merge_tabfiles("$wrkdir/input_full.txt", join(',', @ifields),
									$vepout, join(',', @$fnames2), $outfile);
	unless($igkeys eq join(',', @igkeys)) {
		croak "Incorrect keys were used for merging input and vep output: $igkeys";
	}
}


