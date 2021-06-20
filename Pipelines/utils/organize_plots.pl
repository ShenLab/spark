#!/usr/bin/env perl
use strict;
use warnings;
use Carp;
use IO::File;
use List::Util qw|max|;
use FindBin qw|$Bin|;
use Perl6::Slurp;
use POSIX qw|ceil|;
use File::Path qw|make_path remove_tree|;
use Genet::Ped;
use Utils::List qw|split_list|;
use Getopt::Euclid;
use String::ShellQuote;

use lib "$Bin/../lib";
use Shared qw|fam_rels|;


my $rootdir = $ARGV{'--rootdir'};
my $rootdirq = shell_quote($rootdir);
my $varlist;
unless(-f "$rootdir/par/var_list.txt") {
	if (defined $ARGV{'--prefix'}) {
		carp "Cannot find variant list $rootdir/par/var_list.txt";
	}
	else {
		croak "Cannot find variant list $rootdir/par/var_list.txt"
	}
}
else {
	$varlist = "$rootdir/par/var_list.txt";
}

my ($prefix, $wrkdir);
if (defined $ARGV{'--prefix'}) {
	$prefix = $ARGV{'--prefix'}."_";
	$wrkdir = $ARGV{'--prefix'};
	$varlist = "$rootdir/$ARGV{'--prefix'}/var_list.txt";
	unless(-f "$rootdir/$ARGV{'--prefix'}/var_list.txt") {
		croak "Can also not find variant list $varlist";
	}
}
else {
	$prefix = "";
	$wrkdir = "wrk";
}

unless(-d "$rootdir/$wrkdir") {
	croak "Cannot find $wrkdir directory under $rootdir";
}

#open my $flst, "$rootdir/par/exp_list.txt" or die "Cannot write to expected list";

#my $prefix = $ARGV{'--prefix'} ? $ARGV{'--prefix'}."_" : "";

foreach my $subdir (qw|byvar byiid byfiid|) {
	unless ($ARGV{'--ped'}) {
		next if $subdir eq "byfiid";
	}
	remove_tree("$rootdir/${prefix}$subdir") if -d "$rootdir/${prefix}$subdir";
	make_path "$rootdir/${prefix}$subdir";
	if ($ARGV{'--montage'}) {
		remove_tree("$rootdir/${prefix}${subdir}_montage") if -d "$rootdir/${prefix}${subdir}_montage";
		make_path "$rootdir/${prefix}${subdir}_montage";
	}
}

my (%fids, %pheno, %gender, $famsamps, $famrels, $famlayout);
if ($ARGV{'--ped'}) {
	#%fids = map { (split)[1,0] } slurp $ARGV{'--ped'};
	open my $fin, $ARGV{'--ped'} or die "Cannot open $ARGV{'--ped'}";
	while(<$fin>) {
		my ($fid, $iid, $sex, $pheno) = (split)[0,1,4,5];
		$fids{$iid} = $fid;
		$gender{$iid} = $sex == 2 ? 'Female' : $sex == 1 ? 'Male' : 'NA';
		$pheno{$iid} = $pheno == 2 ? 'Aff' : $pheno == 1 ? 'Unaff' : 'NA';
	}
	($famsamps, $famrels) = fam_rels($ARGV{'--ped'}, { shorten => 1, strict => 0,
		 twins => $ARGV{'--ped-twins'}, ignore => $ARGV{'--ped-ignore'}  });
	if ($ARGV{'--montage'}) {
		$famlayout = fam_layout($ARGV{'--ped'});
	}
}
else {
	if($ARGV{'--montage'}) {
		warn "The \"--montage\" option will be ignored because PED file is not provided";
	}
}

my $suff = $ARGV{'--suffix'};

my $fin = IO::File->new($varlist) or die "Cannot open $varlist";
while(<$fin>) {
	my ($varid, $iid, $carrier) = split;
	my $varidout = abbrev_varid($varid);
	unless (-d "$rootdir/$wrkdir/$iid" && -f "$rootdir/$wrkdir/$iid/$varid.$suff") {
		croak "Cannot find plot of $varid for $iid";
	}

	my @carriers;
	if (defined $carrier) {
		@carriers = split(',', $carrier);
	}

	my $width = `identify -format %w $rootdirq/$wrkdir/$iid/$varid.$suff`;
	my $height = `identify -format %h $rootdirq/$wrkdir/$iid/$varid.$suff`;

	my $size;
	unless ($ARGV{'--montage-size'}) {
		$size = sprintf("%dx%d", $width, $height);
	}
	else {
		$size = $ARGV{'--montage-size'};
	}

	my ($c_width, $c_height);
	if ($ARGV{'--caption'}) {
		$c_width = $width; # int($ARGV{'--caption-width'}*$width);
		$c_height = int($ARGV{'--caption-height'}*$height);
	}

	if ($ARGV{'--caption'}) {	
		my $iidout;
		if (defined $pheno{$iid}) {
			$iidout = "$iid ($gender{$iid},$pheno{$iid})";
		}
		else {
			$iidout = $iid;
		}
		system(qq|convert -background black -fill white -gravity center -size ${width}x${c_height} |.
			qq|caption:\"$iidout : $varidout\" $rootdirq/$wrkdir/$iid/$varid.$suff +swap -gravity $ARGV{'--caption-pos'} |.
			qq|-composite $rootdirq/$wrkdir/$iid/${varid}-caption.$suff|);
		symlink "../$wrkdir/$iid/$varid-caption.$suff", "$rootdir/${prefix}byvar/${varid}-${iid}.$suff";
		symlink "../$wrkdir/$iid/$varid-caption.$suff", "$rootdir/${prefix}byiid/${iid}-${varid}.$suff";
		#print $flst "$rootdir/$wrkdir/$iid/$varid-caption.$suff\n";
	}
	else {
		symlink "../$wrkdir/$iid/$varid.$suff", "$rootdir/${prefix}byvar/${varid}-${iid}.$suff";
		symlink "../$wrkdir/$iid/$varid.$suff", "$rootdir/${prefix}byiid/${iid}-${varid}.$suff";
	}
	#print $flst "$rootdir/$wrkdir/$iid/$varid-caption.$suff\n";
	#print $flst "$rootdir/${prefix}byiid/${iid}-${varid}.$suff\n";

	if (defined $fids{$iid}) {
		my $famid = $fids{$iid};
		if ($ARGV{'--caption'}) {
			symlink "../$wrkdir/$iid/$varid-caption.$suff", "$rootdir/${prefix}byfiid/${famid}-${iid}-$varid.$suff";
		}
		else {
			symlink "../$wrkdir/$iid/$varid.$suff", "$rootdir/${prefix}byfiid/${famid}-${iid}-$varid.$suff";
		}
		#print $flst "$rootdir/${prefix}byfiid/${famid}-${iid}-$varid.$suff\n";
		
		foreach my $sampid (@{$famsamps->{$famid}}) {
			next if $sampid eq $iid;
			my $relation = $famrels->{$iid}{$sampid} //
				do { croak "Cannot find relation between $iid and $sampid" };
			unless (-f "$rootdir/$wrkdir/$sampid/$varid.$suff") {
				carp "Cannot find plot ${iid}'s $relation: $rootdir/$wrkdir/$sampid/$varid.$suff, skipping ...";
				next;
			}
			
			if ($ARGV{'--caption'}) {
				(my $relationship = $relation) =~ s/\'s/'s /;
				if (grep { $sampid eq $_ || uc($relation) eq uc($_) } @carriers) {
					system(qq|convert -background black -fill white -gravity center -size ${c_width}x${c_height} |.
						qq|caption:\"$sampid ($relationship,$pheno{$sampid}) : $varidout\" $rootdirq/$wrkdir/$sampid/$varid.$suff +swap -gravity $ARGV{'--caption-pos'} |.
						qq|-composite $rootdirq/$wrkdir/$sampid/${varid}-caption.$suff|);
				}
				else {
					system(qq|convert -background black -fill white -gravity center -size ${c_width}x${c_height} |.
						qq|caption:\"$sampid ($relationship,$pheno{$sampid})\" $rootdirq/$wrkdir/$sampid/$varid.$suff +swap -gravity $ARGV{'--caption-pos'} |.
						qq|-composite $rootdirq/$wrkdir/$sampid/${varid}-caption.$suff|);
				}
				symlink "../$wrkdir/$sampid/$varid-caption.$suff", "$rootdir/${prefix}byvar/${varid}-${iid}'s${relation}(${sampid}).$suff";
				symlink "../$wrkdir/$sampid/$varid-caption.$suff", "$rootdir/${prefix}byiid/${iid}'s${relation}(${sampid})-${varid}.$suff";
				symlink "../$wrkdir/$sampid/$varid-caption.$suff", "$rootdir/${prefix}byfiid/${famid}-${iid}'s${relation}(${sampid})-$varid.$suff";
				#print $flst "$rootdir/$wrkdir/$sampid/${varid}-caption.$suff\n";
			}
			else {
				symlink "../$wrkdir/$sampid/$varid.$suff", "$rootdir/${prefix}byvar/${varid}-${iid}'s${relation}(${sampid}).$suff";
				symlink "../$wrkdir/$sampid/$varid.$suff", "$rootdir/${prefix}byiid/${iid}'s${relation}(${sampid})-${varid}.$suff";
				symlink "../$wrkdir/$sampid/$varid.$suff", "$rootdir/${prefix}byfiid/${famid}-${iid}'s${relation}(${sampid})-$varid.$suff";
			}
			#print $flst "$rootdir/${prefix}byvar/${varid}-${iid}'s${relation}(${sampid}).$suff\n";
			#print $flst "$rootdir/${prefix}byiid/${iid}'s${relation}(${sampid})-${varid}.$suff\n";
			#print $flst "$rootdir/${prefix}byfiid/${famid}-${iid}'s${relation}(${sampid})-$varid.$suff\n";
		}
	

		if ($ARGV{'--montage'}) {
			my $rows = $famlayout->{$iid} // do { croak "Cannot find layout for $iid" };
			my $nrow = @$rows;
			my $ncol = max(map { scalar(@{$rows->[$_-1]}) } 1..$nrow);
			# determine the files associated with each sample on the grid
			my @files;
			for(my $ii = 0; $ii < $nrow; $ii ++) {
				for(my $jj = 0; $jj < $ncol; $jj ++) {
					my $sampid = $rows->[$ii][$jj];
					if (defined $sampid && -f "$rootdir/$wrkdir/$sampid/$varid.$suff") {
						if ($ARGV{'--caption'}) {
							push @files, "\'$rootdir/$wrkdir/$sampid/$varid-caption.$suff\[$size\]\'";
						}
						else {
							push @files, "\'$rootdir/$wrkdir/$sampid/$varid.$suff\[$size\]\'";
						}
					}
					else {
						push @files, "\'xc:white[$size]\'";
					}
				}
			}
			my $cmd = "montage -mode concatenate -tile ${ncol}x${nrow} @{[ join(qq| |, @files) ]} $rootdirq/$wrkdir/$iid/${varid}-montage.$suff";
			#print $cmd, "\n";
			system($cmd);

			symlink "../$wrkdir/$iid/${varid}-montage.$suff", "$rootdir/${prefix}byvar_montage/${varid}-${iid}.$suff";
			symlink "../$wrkdir/$iid/${varid}-montage.$suff", "$rootdir/${prefix}byiid_montage/${iid}-${varid}.$suff";
			#print $flst "$rootdir/${prefix}byvar_montage/${varid}-${iid}.$suff\n";
			#print $flst "$rootdir/${prefix}byiid_montage/${iid}-${varid}.$suff\n";
			if (defined $fids{$iid}) {
				symlink "../$wrkdir/$iid/${varid}-montage.$suff", "$rootdir/${prefix}byfiid_montage/$fids{$iid}-${iid}-$varid.$suff";
			#print $flst "$rootdir/${prefix}byfiid_montage/$fids{$iid}-${iid}-$varid.$suff\n";
			}
		}
	}
	else {
		carp "Cannot find family ID for $iid\n";
	}
}


=head2 fam_layout

Create layout for plotting family members. Typically used for nuclear families.
Foreach sample: it gives the rows of sample alias that will be appear in the output.
In each family, founders will be in first row, other samples in second row.
In each row, samples are sorted by sex and then their name
Example: SP012 => [[SP012'sFather(SP010),SP012'sMother(SP011),undef],[SP012,SP12'sBrother(SP013),SP12'Brother2(SP014)]]

When the number of plots in two rows differ by 2 fold, additional rows will be inserted to balance
the plots arranged in each row.

=cut


sub fam_layout {
	my ($pedfile) = @_;
	my $ped = Genet::Ped->new($pedfile);

	my %layout;
	foreach my $fid ($ped->get_famids) {
		my @founders = sort { $ped->get_sex($fid, $a) cmp $ped->get_sex($fid, $b) || $a cmp $b } 
							$ped->get_founders($fid);
		my @nonfounders = sort { $ped->get_sex($fid, $a) cmp $ped->get_sex($fid, $b) || $a cmp $b } 
							$ped->get_nonfounders($fid);

		foreach my $iid ($ped->get_members($fid)) {
			my @rows;
			foreach my $founder (@founders) {
				if ($founder eq $iid) {
					push @{$rows[0]}, $iid;
				}
				else {
					push @{$rows[0]}, $founder;
				}
 			}
 			foreach my $nonfounder (@nonfounders) {
 				if ($nonfounder eq $iid) {
 					push @{$rows[1]}, $iid;
 				}
 				else {
 					push @{$rows[1]}, $nonfounder;
 				}
 			}
 			# Re-balance the rows if needed
 			if (defined $rows[0] && defined $rows[1]) {
 				my ($top, $bottom) = (scalar(@{$rows[0]}), scalar(@{$rows[1]}));
 				if ( $top >= 2 && $bottom >= $top * 2 || $top == 1 && $bottom >= 4  ) {
 					#print STDERR $top, "\t", $bottom, "\n";
 					my @segs = split_list($rows[1], calc_nseg($top, $bottom));
 					splice @rows, 1, 1, @segs;
 				}
 				elsif ( $bottom >= 2 && $top >= $bottom * 2 || $bottom == 1 && $top >= 4 ) {
 					#print STDERR $top, "\t", $bottom, "\n";
 					my @segs = split_list($rows[0], calc_nseg($bottom, $top));
 					splice @rows, 0, 1, @segs;
 				}
 			}
 			$layout{$iid} = \@rows;
		}
	}

	if (wantarray) {
		return %layout;
	}
	else {
		return \%layout;
	}
}

# Determine the number of splits (make the output as "square" as possible)
sub calc_nseg {
	my ($small, $large) = @_;
	my $ncol = ceil(sqrt($large));
	if ($ncol < $small) {
		$ncol = $small;
	}
	return $ncol;
}

sub abbrev_varid {
	my ($varid) = @_;
	my $varlen = length($varid);
	if ($varlen <= $ARGV{'--max-varlen'}) {
		return $varid;
	}
	else {
		my $halflen = int($ARGV{'--max-varlen'}/2)-1;
		return substr($varid, 0, $halflen)."..".substr($varid,$varlen-$halflen,$halflen);
	}
}



__END__

=head1 NAME

organize_plots.pl - Re-orgnize variant plots for easy manual curation.

=head1 REQUIRED ARGUMENTS
 
=over
 
=item -[-]root[dir] [=] <dir>

The workflow root directory. 
The script will look for images for each sample in wrk, and look for variant list in par.
Additional dirs will be created: byvar, byiid, and byfiid (if ped file is provided)
And if --montage is switched on, the three additional dirs with _montage suffix will 
be created.

=back
 
=head1 OPTIONS
 
=over

=item -[-]ped [=] <file>

Pedigree file. We will parse the PED file under non-strict mode.
i.e., even if a parent does not appear in the PED file they will
be considered parent.

=for Euclid:
	file.type: readable
 
=item -[-]ped-ignore [=] <pattern>

Pattern for duplicated samples in the PED file.

=for Euclid:
	pattern.default: '_Re(\d*)$'

=item -[-]ped-twins [=] <list>

List of twin pairs.

=item -[-]prefix [=] <string>

Directory name prefix. When prefix is provided, it will also be used as working directory name.
And the variant list can also be optionally found in this working directory.

=item -[-]suffix [=] <string>

The plot file name suffix.

=for Euclid:
	string.default: "JPG"

=item -[-]montage

Also put all family members in one image.

=item -[-]montage-size [=] <size>

Image size under montage mode, default is the actual size.

=item -[-]caption

Add caption to image.

=item -[-]caption-pos [=] <location>

The location of added captions.

=for Euclid:
	location.default: 'north'

=item -[-]caption-width [=] <width>

The relative width of added caption

=for Euclid:
	width.default: 0.90

=item -[-]caption-height [=] <heigth>

The relative height of added caption

=for Euclid:
	heigth.default: 0.05

=item -[-]max-varlen [=] <length>

MAx. length of VarID. VarID longer than this will be abbreviated.

=for Euclid:
	length.default: 36


=back

