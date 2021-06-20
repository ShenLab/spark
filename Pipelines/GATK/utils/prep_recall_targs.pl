use strict;
use warnings;
use File::Path qw|make_path|;
use FindBin qw|$Bin|;
use Utils::File::Iter qw|iter_file|;
use Genome::Ranges::IntSet;
use Getopt::Euclid;

use lib "$Bin/../../lib";
use Shared qw|parse_tabfile|;
use Variants qw|sort_vars|;

if (-d $ARGV{'--outdir'}) {
	unless($ARGV{'--force'}) {
		print STDERR "The output directory $ARGV{'--outdir'} already exists!\n";
		exit 1;
	}
}
else {
	make_path $ARGV{'--outdir'};
}

my (%famsamps, %famid);
if ($ARGV{'--ped'}) {
	open my $fin, $ARGV{'--ped'} or die "Cannot open PED file!";
	while(<$fin>) {
		my ($fid, $iid) = (split)[0,1];
		$famid{$iid} = $fid;
		push @{$famsamps{$fid}} => $iid;
	}
}

my ($it, $fnames, $keyfields)  = parse_tabfile($ARGV{'--tab'},  $ARGV{'--fields'}, 4, 5);

my (%targs, %vars, %allvars);
my $merged = Genome::Ranges::IntSet->new();

while(my $dat = $it->()) {
	my ($iid, $chrom, $pos, $ref, $alt) = @{$dat}{@$keyfields};
	my $grpid;
	if (defined $famid{$iid}) {
		$grpid = $famid{$iid};
	}
	else {
		$grpid = $iid;
	}
   unless (defined $targs{$grpid}) {
        $targs{$grpid} = Genome::Ranges::IntSet->new();
    }


    my $start = $pos-100 < 0 ? 1 : $pos - $ARGV{'--padding'};
    my $end = $pos + length($ref) + $ARGV{'--padding'};
    $chrom = "chr".$chrom if $ARGV{'--chr'} && $chrom !~ /^chr/;
    
    $targs{$grpid}->add($chrom, $start, $end);
    $merged->add($chrom, $start, $end);

    next unless $ARGV{'--alleles'};

   	my $varid = join(":", $chrom, $pos, $ref, $alt);

   	$vars{$grpid}{$varid} = 1;
   	$allvars{$varid} = 1;
}

$merged->write("$ARGV{'--outdir'}.bed", { bed => 1 });
write_vcf(\%allvars, "$ARGV{'--outdir'}.vcf") if $ARGV{'--alleles'};


my $fgrp;
if ($ARGV{'--ped'}) {
	open $fgrp, ">$ARGV{'--outdir'}.groups.txt" or die "Cannot write to $ARGV{'--outdir'}.groups.txt";
}
else {
	open $fgrp, ">$ARGV{'--outdir'}.samples.txt" or die "Cannot write to $ARGV{'--outdir'}.samples.txt";
}

foreach my $grpid (sort keys %targs) {
    if ($ARGV{'--bysamp'} && defined $famsamps{$grpid}) {
    	foreach my $sampid (@{$famsamps{$grpid}}) {
    		$targs{$grpid}->write("$ARGV{'--outdir'}/$sampid.bed", { bed => 1 });
    		write_vcf($vars{$grpid}, "$ARGV{'--outdir'}/$sampid.vcf") if $ARGV{'--alleles'};
    	}
    }
    else {
    	$targs{$grpid}->write("$ARGV{'--outdir'}/$grpid.bed", { bed => 1 });
    	write_vcf($vars{$grpid}, "$ARGV{'--outdir'}/$grpid.vcf") if $ARGV{'--alleles'};
    }
   	#next unless defined $fgrp;

   	if ($ARGV{'--ped'}) {
   		if (defined $famsamps{$grpid}) {
    		foreach my $iid (@{$famsamps{$grpid}}) {
    			print $fgrp $iid, "\t", $grpid, "\n";
    		}
   		}
	    else {
	    	print $fgrp $grpid, "\t", $grpid, "\n";
    	}
   	}
   	else {
   		print $fgrp $grpid, "\n";
   	}
    
}


sub write_vcf {
	my ($vars, $filename) = @_;
	open my $fvcf, "| bgzip -c > $filename.gz" or die "Cannot write to $filename.gz!";
	print $fvcf "##fileformat=VCFv4.0\n";
	print $fvcf '#'.join("\t", qw|CHROM POS ID REF ALT QUAL FILTER INFO|), "\n";
	foreach my $var (sort_vars(keys %$vars)) {
		my ($chrom, $pos, $ref, $alt) = @$var;
		my $varid = join(":", @$var);
		print $fvcf join("\t", $chrom, $pos, $varid, uc($ref), uc($alt), '100', 'PASS', '.'), "\n";
	}
	close $fvcf;
	system(qq|tabix -p vcf $filename.gz|);
}


__END__

=head1 NAME 

prep_recall_targs.pl -- Prepare targets for variant recall.

=head1 REQUIRED ARGUMENTS
 
=over
 
=item -[-]tab [=] <table>

Input variant table. Must be tab-separated.

=item -[-]out[dir] [=] <dir>

Output directory, containing sample or group-specific intervals.
Also used as output prefix for group file if PED file is provided.
 
=back
 
=head1 OPTIONS
 
=over

=item -[-]fields [=] <fields>

Fields for IID,Chrom,Position from the input variant table.

=for Euclid:
	fields.default: "IID,Chrom,Position,Ref,Alt"

=item -[-]ped [=] <pedfile>

If pedigree file is provided, variant will be re-called on all family members.

=item -[-]pad[ding] [=] <length>

Flanking length added to each of the variant (default: 100bp).

=for Euclid:
	length.default: 100

=item -[-]chr

Adding chr prefix to chromosome names to the output.

=item -[-]bysamp

Write one target file per-sample instead of per-family when PED file is provided.
The recall regions still include the union of all family members.

=item -[-]allele

Also output alleles in VCF format.

=item -[-]force

Over-write output directory.

=back

