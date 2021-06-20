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
use List::MoreUtils qw|all uniq any|;
use String::ShellQuote;
use Hash::Util qw|lock_keys|;
use Getopt::Lucid qw|:all|;
use Utils::List qw|insert_before|;
use Utils::File qw|count_line|;
use Utils::File::Iter qw|iter_file|;


use lib "$Bin/../lib/";
use Shared qw|check_conf parse_fstr merge_tabfiles slurp_xref|;

############################
## Command line interface ##
############################

my @spec =  (
	Param("conf|c")->valid(sub { -r }),
	Param("input|in|i")->valid(sub { -r }),
	Param("wrkdir|wrk|w"),
	Param("output|out|o"),
	Switch("force"), # Under force mode, previously generated files will not be over-written
	Switch("help|h")
	);

my $opt = Getopt::Lucid->getopt(\@spec);

if ($opt->get_help) {
	print STDERR <<EOF;
Purpose:
	Annotate sequence variants using VEP/ANNOVAR and custom tools.

Usage:
	anno_seqvars.pl --conf CONFIG [--in INPUT] --wrk WRKDIR -out OUTPUT	

Options:
	--in: Input variant table. It overrides the variant table in config file.
	--wrk: Workding directory for storing intermediate files.
	--force: Force running the program even if the output exist. The default is to exit if output file already exists.

Notes:
	The script collects the following types of annotations for each variant: 

	1. Functional effects to genes & transcripts
	   VEP is used to annotate functional effects to genes and transcripts, including gene/transcript IDs, cDNA changes, 
	   codon and amino acid changes. The default is to use GENCODE Basic set of Ensembl genes, but it is customizable to 
	   use all Ensembl genes or any sources of gene models defined in GFF/GTF. The Ensembl genes used depend on the version
	   of VEP and cache database files. It is recommended to use the version that matches the genes used for cross-referencing 
	   additional gene/transcript level information. When a variant is annotated to multiple genes, each type of gene-specific 
	   annotations will be packed into the one field and joined by ";". It is typical for a gene to have multiple transcripts, 
	   each type of transcript-specific annotations also packed and joined by ',' for each gene.

	2. Gene-specific variant (pathogenicity) annotations
	   Most exome variant pathogenicity scores (including missense or splicing effect prediction scores) are typically 
	   available as lookup tables. They can be extracted from those lookup tables by querying variant position and alleles, 
	   and further matching gene ID or amino acid change (for missense prediction). This is done by query_genevars.pl utility
	   with options specified in Missense section of the config. 

	3. Position and allele-specific annotations
	   Annotations that only depend on genomic position and/or variant alleles includes population allele frequencies, CADD 
	   scores, and evolutionary conservations, etc. VEP can be used to extract selected INFO fields from VCF file by querying 
	   variant positions and matching alleles. It can also extract scores from bigWig file by querying positions. ANNOVAR
	   also has a similar function to query position and allele specific annotations from its database files. For large 
	   ANNOVAR database table, it is recommended to index them (create_annovar_idx.pl) to speed up the position-based query. 
	   Options can be specified in either VEP or ANNOVAR section of config file.

	4. Overlapping genomic features.
	   Genomic features are stored in UCSC BED format, with an optional fourth column for feature ID. Examples include 
	   evolutionary conserved elements, missense constrained regions, CpG islands, cytobands, etc. Genomic position-based
	   query of overlapping features can be done by VEP custom annotation which also require the BED to be bgzipped and 
	   tabix indexed, or by ANNOVAR region-based annotation which require the BED file as part of its database. Options 
	   can be specified in either VEP or ANNOVAR section of config file.

	5. Gene/transcript level information 
       This type of information include gene metrics like loss-of-function intolerance, ontology terms like associated 
       pathways or phenotyes, quantitative measures like expression level etc. They can be extracted from one or more
       lookup tables cross-referencing by gene or transcript IDs. It is recommended to standardize the available 
       data sources to the same set and versions of genes used by the annotation. Cross-reference by gene name is not 
       used by this script because gene names changes over time, and it is diffcult to keep track of versions. Options
       for gathering gene level information are specified in Gene section of config file.

	6. Sample level information
	   This type of information is extracted from one or more lookup table by cross-referencing sample IDs if sample ID field
	   is available in the input variant table. Options are specified in Sample section of config file.

	Several additional standalone utilities can be used to add further annotations to the output:

		query_vars.pl - Query position-based annotations from VCF or ANNOVAR database file.
						It extends the functions of VEP/ANNOVAR: it supports matching additional fields between	input variant 
						table and database tht are not limited to ref and alt alleles. For VCF file, it expands the field values
						by field type, performs variant normaliation internally before matching to the input variant. 
						When multiple matching records can be found	in the database, it supports collecting all values. 
		query_genevars.pl - Query gene-specific variant annotations. This is used as part of this pipeline to collect 
						missense or splicing effect prediction scores.
		overlap_vars.pl - Query overlapping genomic features from genomic intervals or scores from begraph or bigWig.
						It extends the functions of VE/ANNOVAR: it supports any file with chromosome, start and end columns to
						store genomic intervals. File format of BED, bedGraph, bigWig will be automatially recognized. 
						When multiple overlapping features are found, it is possible to apply collapsing functions to summarize
						the results.
		anno_genes.pl - Annotate additional gene level information from one external lookup table.
						It supports the use gene names or other IDs to lookup gene information, and accounts for multiple
						overlapping genes for each variant.
		lookup_varwt.pl - Annotated gene or category specific variant weight.
						  It maps a continuous score for each variant to a weight that is specific to its gene or other categorical
						  annotations of the gene. It is used to annotate PopScore from PSAP (Wilfert 2016) or variant weights for 
						  DenovoWEST.
		calc_pext.pl - Calculate pext measure for tissue-specific expressions (Cummings B et al. 2020).


Input/output:
	Input variant table can be specified in the config file, or overriden by the command line option -i.
	If input table contain indels but not generated by inhouse variant filtering pipeline, it is recommended to normalize
	the variants running annotations (norm_vars.pl). Input file and other lookup tables must be tab separated text or comma 
	separated csv files. Input can also be a directory of many input files, each of which should have the same format. 
	In such case, annotations to each input file will be run in parallel using GNU parallel (SGE is not supported), then 
	individual annotation outputs will be merged to create the final output.

	Ouptut will have one line for each variant, when one variant can be mapped to multiple genes, different gene-level 
	annotations will be packed in one field. Default fields in the output are explained in config template. Additional 
	fields can be customized by different options. 

Dependencies:
	VEP (>= Ensembl release 74) and ANNOVAR (modified, included in the pipeine directory)

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
check_conf(\%config);

my $wrkdir = $opt->get_wrkdir;
if ($wrkdir) {
	make_path $wrkdir unless -d $wrkdir;
}
else {
	$wrkdir = tempdir(CLEANUP => 1);
}

# Annotate all inputs from a directory
if (defined $infile && -d $infile) {
	my @varlists = sort grep { !/^\./ && -f "$infile/$_" } IO::Dir->new($infile)->read();
	my $findopt = "";
	if (exists $config{Input}{Suffix}) {
		my $suffix = $config{Input}{Suffix};
		@varlists = grep { /\.$suffix$/ } @varlists;
		$findopt = "-name '*.$suffix'";
	}
	my $option = "";
	if (exists $config{Parallel} && defined $config{Parallel}{jobs}) {
		$option = join(" ", map { "--$_ $config{Parallel}{$_}" } keys %{$config{Parallel}});
	}
	else {
		die "Must provide number of threads in config file for parallel";
	}
	my $args = "";
	if ($opt->get_force) {
		$args = "--force";
	}

	my $command = qq(find $infile/ -type f $findopt | parallel --eta $option "perl $0 --in {} --wrk $wrkdir/{/}.wrkdir --conf $config --out $wrkdir/{/} 1>$wrkdir/{/}.out $args 2>$wrkdir/{/}.err");
	print $command, "\n"; 
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
}
unless (-f $config{Input}{File}) {
	croak "Cannot find input file: $config{Input}{File}";
}
lock_keys(%config);

while (my ($label, $file) = each %{$config{Path}}) {
	if ($label =~ /DB$/ || $label =~ /^ANNO/) {
		croak "Cannot find $label directory: $file" unless -d $file;
	}
	else {
		croak "Cannot find $label file: $file" unless -f $file;
	}
}

croak "Cannot find input file" unless -f $config{Input}{File}; 

my $input_header = 1;
# Parse and prepare input variant file
my %known;
{
	my $finput = parse_fstr($config{Input}{Fields}, 1);
	my @fvar = keys %$finput;
	if (all { /^\d+$/ } @fvar) {
		$input_header = 0;
	}
	croak "Input fields must be renamed to Chrom,Position,Ref,Alt" 
		unless join(',', values %$finput) eq 'Chrom,Position,Ref,Alt';
	my ($it, $fnames) = iter_file($config{Input}{File}, { header => $input_header, fsep => qr/\t/ });
	foreach my $field (@fvar) {
		unless(grep { $field eq $_ } @$fnames) {
			croak "Cannot find key field $field in the input file";
		}
	}
	if (exists $config{Input}{XtraFields}) {
		my $fxtra = parse_fstr($config{Input}{XtraFields}, 1);
		foreach my $field (keys %$fxtra) {
			unless(grep { $field eq $_ } @$fnames) {
				croak "Cannot find extra field $field in the input file";
			}
		}
	}
	if (exists $config{Sample} && exists $config{Sample}{IDField}) {
		my $fsamp = parse_fstr($config{Sample}{IDField}, 1);
		foreach my $field (keys %$fsamp) {
			unless(grep { $field eq $_ } @$fnames) {
				croak "Cannot find sample field $field in the input file";
			}
		}
	}

	# Prepare VCF
	my $fvcf = IO::File->new("$wrkdir/input.vcf", "w");
	print $fvcf "##fileformat=VCFv4.0\n";
	print $fvcf '#'.join("\t", qw|CHROM POS ID REF ALT QUAL FILTER INFO|), "\n";

	while(my $dat = $it->()) {
		my ($chr, $pos, $ref, $alt) = @{$dat}{@fvar};
		unless ($ref =~ /^[ACGT]+$/i && $alt =~ /^[ACGT]+$/i) {
			warn "Incorrect allele : $chr:$pos:$ref:$alt";
			next;
  		}
		if (uc($ref) eq uc($alt)) {
			warn "Variant $chr:$pos:$ref:$alt has the same ref and alt allele, skipped";
			next;
		}
		# VEP is based on ensembl, chr prefix is not needed
  		$chr =~ s/^chr//;
  		$pos =~ s/,//g;
  		#my $vid = join(":", @{$dat}{@fvar});
  		my $vid = join(":", $chr, $pos, $ref, $alt);
  		next if defined $known{$vid};
  		print $fvcf join("\t", $chr, $pos, $vid, uc($ref), uc($alt), '100', 'PASS', '.'), "\n";
  		$known{$vid} = 1;
	}

	# If ANNOVAR configs present
	# Using ANNOVAR utility to convert to avinput format (keep VarID in the last column)
	if (exists $config{ANNOVAR}) {
		system("$config{Path}{ANNOVAR}/convert2annovar.pl -format vcf4 $wrkdir/input.vcf ".
		   "-includeinfo | cut -f1-5,8 > $wrkdir/input.avinput");
		# Validate ANNOVAR output
		my $nline = count_line("$wrkdir/input.avinput");
		unless(scalar(keys %known) == $nline) {
			print STDERR "Known=", scalar(keys %known), ", NLines=$nline\n";
			die "Incorrect number of variants in $wrkdir/input.avinput";
		}	
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
		my $subdir = "homo_sapiens";
		if (defined $config{VEP}{Option}) {
			if ($config{VEP}{Option} =~ /\-\-species (\w+)/) {
				$subdir = $1;
			}
			if ($config{VEP}{Option} =~ /\-\-merged/) {
				$subdir .= "_merged";
			}
			elsif ($config{VEP}{Option} =~ /\-\-refseq/) {
				$subdir .= "_refseq";
			}
		}
		unless(-d "$config{Path}{VEPDB}/$subdir/$config{VEP}{DBVer}_$ver_a") {
			croak "Cache version $config{VEP}{DBVer} for $subdir does not exist!";
		}
	}
	else {
		# If version is not provided, will assume to use the default from cache
		# that matches the VEP version
		$opt_c = "";
	}

	print "Running VEP for gene level annotations\n";
	my $vepcmd = "$config{Path}{VEP} -i $wrkdir/input.vcf -o $wrkdir/vep_out.txt --force --no_stats ".
		"$opt_a --cache --dir_cache $config{Path}{VEPDB} $opt_c --fasta $config{Path}{FASTA} ";
	if (defined $config{VEP}{Option} && $config{VEP}{Option} !~ /^\s*$/) {
		$config{VEP}{Option} =~ s/^["']//; $config{VEP}{Option} =~ s/["']$//;
		$vepcmd .=  $config{VEP}{Option}." ";
	}
	else {
		$vepcmd .= "--gencode_basic --hgvs --symbol --transcript_version --biotype --numbers --offline ";
	}
	my $vepxtra = "";
	# $config{VEP} must also exist
	if (exists $config{VEP}{Custom}) {
		my @vepxtra;
		if (ref $config{VEP}{Custom} eq 'ARRAY') {
			for(my $ii = 0; $ii < @{$config{VEP}{Custom}}; $ii ++) {
				$vepcmd .= "--custom $config{VEP}{Custom}[$ii] ";
				push @vepxtra, $config{VEP}{Custom_Fields}[$ii];
			}
			$vepxtra .= " --extra-var ".join(',', @vepxtra)." ";
		}
		else {
			$vepcmd  .= "--custom $config{VEP}{Custom} ";
			$vepxtra .= "--extra-var"." ".$config{VEP}{Custom_Fields}." ";
		}
	}
	if (exists $config{VEP}{Gene_Fields}) {
		$vepxtra .= "--extra-gene $config{VEP}{Gene_Fields} ";
	}
	if (exists $config{VEP}{Trans_Fields}) {
		$vepxtra .= "--extra-trans $config{VEP}{Trans_Fields} ";
	}
	if (exists $config{VEP}{Plugin}) {
		$vepcmd .= "--dir_plugin ".shell_quote("$Bin/module")." ";
		if (ref $config{VEP}{Plugin} eq 'ARRAY') {
			$vepcmd .= join(" ",  map { "--plugin $_" } @{$config{VEP}{Plugin}});
		}
		else {
			$vepcmd .= "--plugin $config{VEP}{Plugin}";
		}
	}
	print $vepcmd, "\n";
	system($vepcmd);
	# Validate the output: should have the same number of variants as input
	my %vars;
	open my $fin, "$wrkdir/vep_out.txt" or die "Cannot open $wrkdir/vep_out.txt";
	while(<$fin>) {
		next if /^#/;
		my $vid = (split)[0];
		$vars{$vid} = 1;
	}
	unless (scalar(keys %vars) == scalar(keys %known)) {
		print STDERR "Input=", scalar(keys %known), ", VEPOutput=", scalar(keys %vars), "\n";
		die "Incorrect number of variants in $wrkdir/vep_out.txt";
	}

	# Collect VEP results into table, one line per variant
	my ($inclusion, $exclusion) = ("", "");
	if (exists $config{Gene} && exists $config{Gene}{Exclude}) {
		croak "Cannot find gene exclusion list file" unless -f $config{Gene}{Exclude};
		$exclusion = "--exclude-gene $config{Gene}{Exclude} ";
	}
	if (exists $config{Gene} && exists $config{Gene}{Include}) {
		croak "Cannot find gene inclusion list file" unless -f $config{Gene}{Include};
		$inclusion = "--include-gene $config{Gene}{Include} ";
	}
	if (exists $config{Transcript} && exists $config{Transcript}{Exclude}) {
		croak "Cannot find transcript exclusion list file" unless -f $config{Transcript}{Exclude};
		$exclusion .= "--exclude-trans $config{Transcript}{Exclude} ";
	}
	if (exists $config{Transcript} && exists $config{Transcript}{Include}) {
		croak "Cannot find transcript inclusion list file" unless -f $config{Transcript}{Include};
		$inclusion .= "--include-trans $config{Transcript}{Include} ";
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
	my $sumcmd = "perl $utildir/collect_vepout.pl -in $wrkdir/vep_out.txt ".
		"$exclusion $inclusion $coding $cutoff -out $wrkdir/vep_out_tab_raw.txt $vepxtra";
	print $sumcmd, "\n";
	system($sumcmd);
	
	#my $nline = count_line("$wrkdir/vep_out_tab_raw.txt");
	#unless(scalar(keys %known) == $nline-1) {
	#	print STDERR "Known=", scalar(keys %known), ", NLines=$nline\n";
	#	die "Incorrect number of variants in the tabular output: $wrkdir/vep_out_tab_raw.txt";
	#}
	# tabulated vepout may not have the same number of variants as input
	# due to the removal of some genes
	unless(-f "$wrkdir/vep_out_tab_raw.txt.done") {
		die "collect_vepout did not finish successfully";
	}
}


#############################################
### Add gene/transcript level information ###
#############################################

my @gfields;
if (exists $config{Gene} && exists $config{Gene}{GXref} ||
	exists $config{Transcript} && exists $config{Transcript}{TXref}) {
	my ($it, $fnames) = iter_file("$wrkdir/vep_out_tab_raw.txt", { fsep => qr/\t/ });
	@gfields = @$fnames;

	my ($xgdat, $xgfds, $xtdat, $xtfds);
	if (exists $config{Gene}{GXref}) {
		($xgdat, $xgfds) = slurp_xref($config{Gene}{GXref}, $config{Gene}{GXref_Fields});
		foreach my $xfd (@$xgfds) {
			if (any { $xfd eq $_ } @$fnames) {
				croak "Field name $xfd already exist in the VEP output!";
			}
		}
		# Identify GeneID field and splice in additional gene related fields
		my ($gid_ii) = grep { $gfields[$_] eq 'GeneID' } 0..$#gfields;
		die "Cannot find GeneID column" unless defined $gid_ii;
		splice @gfields, $gid_ii+1, 0, @$xgfds;
	}
	# $config{Transcript} may not exist
	if (exists $config{Transcript} && exists $config{Transcript}{TXref}) {
		($xtdat, $xtfds) = slurp_xref($config{Transcript}{TXref}, $config{Transcript}{TXref_Fields});
		foreach my $xfd (@$xtfds) {
			if (any { $xfd eq $_ } @$fnames) {
				croak "Field name $xfd already exist in the VEP output!";
			}
		}
		my ($tid_ii) = grep { $gfields[$_] eq 'TransIDs' } 0..$#gfields; 
		die "Cannot find TransIDs column" unless defined $tid_ii;
		splice @gfields, $tid_ii+1, 0, @$xtfds;
	}

	# Reformat outputs
	my $fout = IO::File->new("$wrkdir/vep_out_tab.txt", "w");
	print $fout join("\t", @gfields), "\n";
	while(my $dat = $it->()) {
		if (defined $xgdat) {
			my @gids = split(';', $dat->{GeneID});
			foreach my $xgfd (@$xgfds) {
				$dat->{$xgfd} = join(";", map { $xgdat->{$_}{$xgfd} // "." } @gids);
			}
		}
		if (defined $xtdat) {
			# Need to first map gene to trans
			# Then add trans level information
			my @gids = split(';', $dat->{GeneID});
			my @tidgrps = split(';', $dat->{TransIDs});
			unless(@gids == @tidgrps) {
				print Dumper $dat;
				die "Gene IDs and trans ID groups do not have the same length";
			}
			foreach my $xtfd (@$xtfds) {
				my @tdat;
				foreach my $tgrp (@tidgrps) {
					push @tdat, join(',', map { $xtdat->{$_}{$xtfd} // "." } split(',', $tgrp));
				}
				$dat->{$xtfd} = join(";", @tdat);
			}
		}
		print $fout join("\t",  map { $_ // "." } @{$dat}{@gfields}), "\n";
	}
}
else {
	copy("$wrkdir/vep_out_tab_raw.txt", "$wrkdir/vep_out_tab.txt");
	my ($it, $fnames) = iter_file("$wrkdir/vep_out_tab.txt", { fsep => qr/\t/ });
	@gfields = @$fnames;
}

#############################################
### Run in-house D-mis annotation script  ###
##############################################

my $vepout = "$wrkdir/vep_out_tab.txt";
if (exists $config{Missense}) {
	my %stdfields = (Chrom => 1, Position => 1, Ref => 1, Alt => 1, GeneID => 1, AAChg => 1);
	if (ref $config{Missense}{DBFile} eq 'ARRAY') {
		unless(ref $config{Missense}{DBFile_Fields} eq 'ARRAY' && 
			   @{$config{Missense}{DBFile}} == @{$config{Missense}{DBFile_Fields}}) {
			die "Unequal length of DBFile files and Fields specs";
		}
	}
	else {
		unless (defined $config{Missense}{DBFile_Fields} || 
				ref $config{Missense}{DBFile_Fields} eq 'ARRAY') {
			die "Cannot find field spec for DBfile or more than one field specs for one DBFile";
		}
		$config{Missense}{DBFile} = [$config{Missense}{DBFile}];
		$config{Missense}{DBFile_Fields} = [$config{Missense}{DBFile_Fields}];	
	}
	for(my $ii = 0; $ii < @{$config{Missense}{DBFile}}; $ii ++) {
		my $dbfile = $config{Missense}{DBFile}[$ii];
		my $dbfields = $config{Missense}{DBFile_Fields}[$ii];
		my $f_mis = parse_fstr($dbfields, 1);
		my $misoption = "";
		foreach my $stdfd (qw|Chrom Position Ref Alt|) {
			unless(grep { $_ eq $stdfd } values %$f_mis) {
				die "Cannot find standard field $stdfd from database field specification: $dbfields";
			}
		}
		if (grep { $_ eq 'AAChg' } values %$f_mis) {
			print STDERR "Collecting missense pathogenicity predictions from $dbfile\n";
			$misoption = '--misonly';
		}
		else {
			unless(grep { $_ eq 'GeneID' } values %$f_mis) {
				die "Must specify GeneID field when AAChg is not available: $dbfields!";
			}
			print STDERR "Collecting gene-specific variant scores from $dbfile\n";
		}
	
		my @f_pred = grep { !defined $stdfields{$_} } values %$f_mis;
		unless(@f_pred) {
			die "Must specify one or more custom field(s) to extract from $dbfile";
		}
		push @gfields, @f_pred;
		my $input = $ii > 0 ? sprintf("%s/vep_out_tab.%d.txt", $wrkdir, $ii-1) : "$wrkdir/vep_out_tab.txt";
		$vepout = sprintf("%s/vep_out_tab.%d.txt", $wrkdir, $ii); 
		my $command = qq|perl $utildir/query_genevars.pl --in $input --output $vepout --db $dbfile --db-fields $dbfields $misoption|;
		print $command, "\n";
		system($command);
		unless(-f $vepout) {
			die "Cannot find missense annotation output: $vepout";
		}
		# Always check if the input and output have the same number of lines
		my $inline = count_line($input);
		my $outline = count_line($vepout);
		unless($inline == $outline) {
			die "Input and output of missense annotation do not have the same number of lines: $input($inline), $vepout($outline)";
		}
	}
}


######################################################
### Run ANNOVAR for filter/region based annotation ###
######################################################

my $afields;
if (exists $config{ANNOVAR}) {
	my (@prots, @ops, @fstr);

	if (exists $config{ANNOVAR}{Filter}) {
		if (ref $config{ANNOVAR}{Filter} eq 'ARRAY') {
			push @prots => @{$config{ANNOVAR}{Filter}};
			push @ops => ('f') x @{$config{ANNOVAR}{Filter}};
			push @fstr => @{$config{ANNOVAR}{Filter_Fields}};
		}
		else {
			push @prots => $config{ANNOVAR}{Filter};
			push @ops => 'f';
			push @fstr => $config{ANNOVAR}{Filter_Fields};
		}
	}
	
	if (exists $config{ANNOVAR}{Region}) {
		my $regtabs = parse_fstr($config{ANNOVAR}{Region}, 1);
		my @regs = keys %$regtabs;
		push @prots, @regs;
		push @ops => ('f') x @regs;
		unshift @fstr, $config{ANNOVAR}{Region};
	}

	for(my $ii = 0; $ii < @prots; $ii ++) {
		unless (-f "$config{Path}{ANNODB}/$config{Input}{HG}_$prots[$ii].txt") {
			croak "Cannot find table $config{Input}{HG}_$prots[$ii] in ANNOVAR DB directory"
		}
	}

	my $protocol = join(',', @prots);
	my $operation = join(',', @ops);
	$afields = join(',', @fstr);

  	system("$config{Path}{ANNOVAR}/table_annovar.pl -buildver $config{Input}{HG} ".
		"-protocol $protocol -operation $operation -nastring . -otherinfo ".
		"-outfile $wrkdir/annovar_out $wrkdir/input.avinput $config{Path}{ANNODB}");
  	my $nline = count_line("$wrkdir/annovar_out.$config{Input}{HG}_multianno.txt");
  	unless(scalar(keys %known) == $nline-1) {
  		print STDERR "Known=", scalar(keys %known), ", NLines=$nline\n";
  		die "Incorrect number of variants in $wrkdir/annovar_out.$config{Input}{HG}_multianno.txt";
  	}
}

##################################
### Merge and reformat results ###
##################################

# Select sample/variant information from the original input file
# Sample level information will appear at the beginning
# Extra variant level info will be appended at the end of the input
my @ifields;
{
	my $infstr = $config{Input}{Fields};
	my $fkey;

	if (exists $config{Input}{XtraFields}) {
		#my $ftmp = parse_fstr($config{Input}{XtraFields}, 1);
		$infstr .= ",".$config{Input}{XtraFields};
	}
	if (exists $config{Sample} && exists $config{Sample}{IDField}) {
		my $ftmp = parse_fstr($config{Sample}{IDField}, 1);
		$fkey = (values %$ftmp)[0];
		$infstr = $config{Sample}{IDField}.",".$infstr;
	}
	
	my $finput = parse_fstr($infstr, 1);
	@ifields = values %$finput;
	# Splice in VarID column into appropriate position
	#my $chrii = first { $ifields[$_] eq 'Chrom' } 0..$#ifields;
	#splice @ifields, $chrii, 0, 'VarID';
	insert_before(\@ifields, 'Chrom', 'VarID');

	my ($it, $fnames) = iter_file("$config{Input}{File}",
								{ header => $input_header, fsep => qr/\t/, 
								  select => $finput });
								#  alias => $finput, subset => [values %$finput] });

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

	my $fout = IO::File->new("$wrkdir/input.txt", "w");
	print $fout join("\t", @ifields), "\n";
	while(my $dat = $it->()) {
		$dat->{Chrom} =~ s/^chr//; $dat->{Chrom} = 'MT' if $dat->{Chrom} eq 'M';
		$dat->{Position} =~ s/,//g;
		$dat->{Ref} = uc($dat->{Ref});
		$dat->{Alt} = uc($dat->{Alt});
		if ($dat->{Ref} !~ /^[ACGT]+$/ || $dat->{Alt} !~ /^[ACGT]+$/ || $dat->{Ref} eq $dat->{Alt}) {
			next;
		}
		$dat->{VarID} = join(":", @{$dat}{qw|Chrom Position Ref Alt|});
		if (defined $xdat && defined $xdat->{$dat->{$fkey}}) {
			@{$dat}{@$xfds} = @{$xdat->{$dat->{$fkey}}}{@$xfds};
		}
		print $fout join("\t", map { $_ // "." } @{$dat}{@ifields}), "\n";
	}
}


# Then merge input, vep_output, and variant annotation
# input and vep_output will be merged based on the VarID/Chrom/Position/Ref/Alt
# the merged file will then be merged with ANNOVAR annotation based on VarID
# All operations will be "left join"
print STDERR "Merging VEP tabular output and selected fields from input table\n";
my @igfields = @ifields;
my @igkeys;
foreach my $gfield (@gfields) {
	if (grep { $gfield eq $_ } @ifields) {
		push @igkeys, $gfield;
	}
	else {
		push @igfields, $gfield;
	}
}
unless (join(',', @igkeys) eq 'VarID,Chrom,Position,Ref,Alt') {
	croak "Incorrect key fields for merging input and vep output";
}


# Compare ifields and gfields to find fields in common.
# Should have {Chrom,Position,Ref,Alt}
my $igkeys = merge_tabfiles("$wrkdir/input.txt",  join(',', @ifields),
							$vepout, join(',', @gfields),
							"$wrkdir/input_vepout_tab.txt");
unless ($igkeys eq 'VarID,Chrom,Position,Ref,Alt') {
	croak "Incorrect key fields were used to merge input and vep output: $igkeys";
}

# ANNOVAR may be optional
if (exists $config{ANNOVAR}) {
	print STDERR "Merging VEP and ANNOVAR results\n";
	my $fanno = parse_fstr($afields);
	my @gvkeys;
	foreach my $afield (values %$fanno) {
		if (grep { $afield eq $_ } @igfields) {
			push @gvkeys, $afield;
		}
	}
	croak "Incorrect key fields for merging input_vep and annovar_out" if @gvkeys > 0;

	my $gvkeys = merge_tabfiles("$wrkdir/input_vepout_tab.txt", join(',', @igfields),
							  	"$wrkdir/annovar_out.$config{Input}{HG}_multianno.txt", "Otherinfo:VarID,$afields",
							  	$outfile);
	unless ($gvkeys eq 'VarID') {
		croak "Incorrect key fields were used to merge input_vepout and annovar_out: $gvkeys";
	}
}
else {
	copy("$wrkdir/input_vepout_tab.txt", $outfile);
}
