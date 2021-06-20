use strict;
use warnings;
use FindBin qw|$Bin|;
use Getopt::Euclid;
use Data::Dumper;
use Module::Load;
use Genet::Var qw|normalize|;
use Genet::File::VCF;
use FaSlice;
use Genome::UCSC::TwoBit;
use Utils::Hash qw|chk_default|;
use Utils::Parser qw|sql_query|;
use Utils::File::Iter qw|iter_file|;

use lib "$Bin/../lib";
use Shared qw|parse_fstr|;
use Variants qw|var_type get_fieldsinfo expand_site|;

my $vcf;
{
	my $fvcf;
	my $region = $ARGV{'--bed'};
	if ($region) {
		if (-f $region) {
			open $fvcf, "tabix -p vcf -h -R $region $ARGV{'--vcf'} |" or die "Cannot open tabix pipe";
		}
		elsif ($region =~ /^(\w+):([\d,]+)\-([\d,]+)$/ || $region =~ /^(\w+)$/) {
			$fvcf = Genet::File::VCF->new($ARGV{'--vcf'});
		}
		else {
			die "Cannot find region file or recognize region spec: $region";
		}
	}
	else {
		$fvcf = $ARGV{'--vcf'};
	}
	$vcf = Genet::File::VCF->new($fvcf);
}

if (defined $ARGV{'--module'}) {
	load $ARGV{'--module'};
}

my $seq;
if ($ARGV{'--seq'}) {
	if ($ARGV{'--seq'} =~ /\.2bit$/) {
		$seq = Genome::UCSC::TwoBit->new($ARGV{'--seq'});
	}
	else {
		FaSlice->new(file => $ARGV{'--seq'});
	}
}

# Get fields type from VCF header or custom specification
my ($sitetype, undef) = get_fieldsinfo($ARGV{'--vcf'}, $ARGV{'--custom'});

# %sfields will be the original field names to be used in expand_site
# Under strict mode, all fields should have a type
# Under default mode, fields without a type is tolerated, and will not be expanded
my %filters;
my %sfields = (CHROM => 1, POS => 1, REF => 1, ALT => 'A');
foreach my $label (qw|filter filter-snv filter-indel|) {
	if (defined $ARGV{"--$label"}) {
		my ($callback, $tokens) = sql_query($ARGV{"--$label"}, 1);
		foreach my $tok (@$tokens) {
			if ($tok->[0] eq 'FIELD') {
				my $origfield = (split(q|\.|, $tok->[1]))[0];
				unless(defined $sitetype->{$origfield}) {
					die "Cannot find type for field used in filter: $origfield!";
				}
				chk_default(\%sfields, $origfield, $sitetype->{$origfield});
			}
		}
		$filters{$label} = $callback;
	}
}


# Parse all original site level fields 
my @sitefields;
foreach my $field (keys %$sitetype) {
	# R:.REF/.ALT, G:.HOMREF/.HET/.HOMALT, number: .1/.2...
	if ($sitetype->{$field} =~ /^(\d+)$/) {
		my $nel = $1;
		if ($nel > 1) {
			push @sitefields, map { $field.".$_" } 1..$nel;
		}
		else {
			push @sitefields, $field;
		}
	}
	elsif ($sitetype->{$field} eq 'G') {
		push @sitefields, $field.".HOMREF", $field.".HET", $field.".HOMALT";
	}
	elsif ($sitetype->{$field} eq 'R') {
		push @sitefields, $field.".REF", $field.".ALT";
	}
	else {
		push @sitefields, $field;
	}
}


my $alias;
if ($ARGV{'--select'}) {
	$alias = parse_fstr($ARGV{'--select'}, 1);
	foreach my $field (keys %$alias) {
		unless (grep { $field eq $_  } @sitefields) {
			die "Cannot find original field $field used in alias";
		}
		my $origfield = (split(q|\.|, $field))[0];
		unless(defined $sitetype->{$origfield}) {
			die "Cannot find type for field used for output: $origfield!";
		}
		chk_default(\%sfields, $origfield, $sitetype->{$origfield});
	}
}
else {
	foreach my $origfield (keys %$sitetype) {
		chk_default(\%sfields, $origfield, $sitetype->{$origfield});
	}
}

# Check output
my $fout;
if ($ARGV{'--output'}) {
	open $fout, ">$ARGV{'--output'}" or die "Cannot write to $ARGV{'--output'}";
}
else {
	$fout = \*STDOUT;
}

my @outfields;
if ($ARGV{'--select'}) {
	@outfields = grep { defined $alias->{$_} } @sitefields;
	print $fout join("\t", map { $alias->{$_} } @outfields), "\n";
}
else {
	@outfields = @sitefields;
	print $fout join("\t", @outfields), "\n";
}


my $it = $vcf->iter({samp => [], strict => 1});
while(my ($site, $geno) = $it->()) {
	# Only for some public available VCF that may not conform to standard format.
	next unless defined $site->{ALT}; 
	# Customize site level INFO if module is provided
	if (defined $ARGV{'--module'}) {
		Genet::File::VCF::custom_site($site);
	}

	my @sites = expand_site($site, \%sfields);
	for(my $ii = 0; $ii < @sites; $ii ++) {
		next if $sites[$ii]{ALT} eq '*';
		if (defined $seq && $sites[$ii]{REF} =~ /^[ACGT]+$/ && $sites[$ii]{ALT} =~ /^[ACGT]+$/) {
			normalize($sites[$ii], $seq);
		}

		# To account for the presence non-strict variants in public data set
		my $VT = var_type($sites[$ii]{REF}, $sites[$ii]{ALT}, 1);
		next if $VT eq 'Unknown';

		my $vt = lc($VT);
		if (defined $filters{filter}) {
			next unless $filters{filter}->($sites[$ii]);
		}
		if (defined $filters{"filter-$vt"}) {
			next unless $filters{"filter-$vt"}->($sites[$ii]);
		}

		print $fout join("\t", map { defined $_ ? (ref $_ eq 'ARRAY' ? join(",", @$_) : $_) : "." } 
							@{$sites[$ii]}{@outfields}), "\n"; 
	}
}


__END__

=head1 NAME

tab_vcfsites.pl -- Extract and reformat variants in site-only VCF file 

=head1 REQUIRED ARGUMENTS

=over

=item -[-]vcf [=] <input>

The input VCF file.

=back

=head1 OPTIONS

=over

=item -[-]out[put] [=] <table>

The output table file. Will output to STDOUT if not provided.

=item -[-]seq [=] <file>

Reference genome sequence file, if provided, will be used for normalization.

=for Euclid:
	file.type: readable

=item -[-]bed [=] <bed>

Region of interest BED file or a region spec.

=item -[-]module [=] <file>

Customizing the site level INFO field by a perl module.

=for Euclid:
	file.type: readable

=item -[-]custom [=] <fstr>

Customizing the type of selected INFO field. By default, we will find field types from VCF header.

=item -[-]filter [=] <expr>

Filter expression. Can only use original field names with expansion if needed.

=item -[-]filter-snv [=] <expr>

Filter expression for SNVs. Can only use original field names with expansion if needed.

=item -[-]filter-indel [=] <expr>

Filter expression for indels. Can only use original field names with expansion if needed.

=item -[-]select [=] <fstr>

Selected fields and alias. Original field names should be expanded if expansion is needed.

=back

=cut

