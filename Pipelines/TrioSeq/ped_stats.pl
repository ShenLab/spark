#!/usr/bin/env perl

use strict;
use warnings;
use Carp;
use Data::Table;
use Data::Dumper;
use List::Util qw|sum|;
use List::MoreUtils qw|all|;
use Text::Table;
use Data::Table::Excel qw|tables2xlsx|;
use Getopt::Euclid;
use Genet::Ped;
use Utils::Hash qw|array2hash array2hist|;


# Transform sex/aff coding into human readable strings
my %plural = ("Unknown Sex" => "Unknown Sex", "Unaffected" => "Unaffected",
	"Affected" => "Affected", "Unknown Pheno" => "Unknown Pheno", 
	"Child(?Sex)" => "Children(?Sex)", "Sib(?Sex)" => "Sibs(?Sex)",
	"Parent(?Sex)" => "Parents(?Sex)");

my %SEX = (1 => "Male", 2 => "Female", 0 => "Unknown Sex");
my %AFF = (1 => "Unaffected", 2 => "Affected", 0 => "Unknown Pheno");
my %AF = (1 => "Unaff", 2 => "Aff", 0 => "NA");
my %PARSEX = (1 => "Father", 2 => "Mother", 0 => 'Parent(?Sex)');
my %CHDSEX = (1 => "Son", 2 => "Daughter", 0 => 'Child(?Sex)');
my %SIBSEX = (1 => "Brother", 2 => "Sister", 0 => 'Sib(?Sex)');
my %POSEX = ('11' => 'Father-Son', '12' => 'Father-Daughter', '10' => 'Father-Child(?Sex)', 
	'22' => 'Mother-Daughter', '21' => 'Mother-Son', '20' => 'Mother-Child(?Sex)');

my $ignore_regex;
if (defined $ARGV{'--ignore'}) {
	$ignore_regex = qr/$ARGV{'--ignore'}/;
}

my $ped = Genet::Ped->new($ARGV{'--ped'}, { splitfam => $ARGV{'--splitfam'}, sep => $ARGV{'--famsep'},
 											ignore => $ignore_regex });

my $fout;
if ($ARGV{'--output'}) {
	$fout = IO::File->new($ARGV{'--output'}, "w");
}
else {
	$fout = \*STDOUT;
}

print  $fout "########################################\n";
printf $fout "A Summary of Different Types of Families\n";
print  $fout "########################################\n";

# Add sort to keep results consistent between different runs
my @fids = sort $ped->get_famids();
my $n_tot = $ped->get_sampsize();
my %n_aff = (affected => 0, unaffected => 0, unknown => 0);
foreach my $fid (@fids) {
	my %pheno = $ped->get_aff($fid);
	while(my ($iid, $phe) = each %pheno) {
		unless (defined $phe) {
			$n_aff{unknown} ++;
		}
		else {
			if ($phe eq '1') {
				$n_aff{unaffected} ++;
			}
			elsif ($phe eq '2') {
				$n_aff{affected} ++;
			}
			else {
				die "Incorrect phenotype code: $phe";
			}
		}
	}
}
my @phecode = sort grep { $n_aff{$_} > 0 } keys %n_aff;

printf $fout "\nFound a total of %d %s in %d %s,\n",  $n_tot, $n_tot > 1 ? "individuals" : "individual", 
	scalar(@fids), @fids > 1 ? "families" : "family";
print $fout "including ", join(", ", map { "$n_aff{$_} $_" } @phecode), "\n";

my %fids_bytyp;
foreach my $fid (@fids) {
	my $fam_type = fam_type($ped, $fid);
	push @{$fids_bytyp{$fam_type}}, $fid;
}

my @alltyps = qw|single duo single-parent trio quad multi-marriage extended multi-parent multi-generation multi-subfam|;
my %tabular = (single => \&tab_singles, duo => \&tab_duos, trio => \&tab_trios, quad => \&tab_quads);

my (@sheets, @labels);
for my $famtyp (@alltyps) {
	if ($fids_bytyp{$famtyp}) {
		my $fids = $fids_bytyp{$famtyp};
		printf $fout "\n* %s: N = %d\n", @$fids > 1 ? ucfirst($famtyp)."s" : ucfirst($famtyp), scalar(@$fids);
		if (@$fids > $ARGV{'--cutoff'} && $tabular{$famtyp}) {
			my $table = $tabular{$famtyp}->($ped, $fids);
			print_txt_tab($table);
			if ($ARGV{'--xlsx'}) {
				push @sheets, $table;
				push @labels, $famtyp;
			}
		}
		else {
			if ($famtyp ne 'multi-subfam') {
				list_fams($ped, $fids, $famtyp);
			}
			else {
				foreach my $fid (@$fids) {
					my @cc = $ped->get_connected($fid);
					my $sep = $ARGV{'--famsep'};
					my $subpeds = $ped->sub_peds($fid, \@cc, { sep => $sep });
					printf $fout " - $fid : %d sub-families\n", scalar(@cc);
					list_fams($subpeds, [map { "$fid$sep$_" } 1..@cc], undef, "  + ");
				}
			}	
		}
	}
}
if ($ARGV{'--xlsx'}) {
	my $xlsx = $ARGV{'--output'};
	$xlsx =~ s/\.\w+$/.xlsx/;
	tables2xlsx($xlsx, [@sheets], [@labels]);
}


# ---------------------------------------------- #
#                 SUBROUTINES                    #
# ---------------------------------------------- #

# Determine the type of family
sub fam_type {
	my ($ped, $famid) = @_;
	my $n_cc = $ped->get_connected($famid);
	if ($n_cc > 1) {
		return "multi-subfam";
	}
	my $f_size = $ped->get_famsize($famid);
	my $n_gen = $ped->get_numgen($famid);
	my $n_founder = $ped->get_founders($famid);
	my $n_offspng = $ped->get_nonfounders($famid);
	if ($f_size == 1) {
		return 'single';
	}
	else {
		if ($n_gen == 2) {
			if ($n_founder == 1) {
				if ($n_offspng == 1) {
					return "duo";
				}
				else {
					return "single-parent"
				}
			}
			elsif ($n_founder > 2) {
				return "multi-marriage";
			}
			else {
				# bi-parental nuclear families
				if ($n_offspng == 1) {
					return "trio";
				}
				elsif ($n_offspng == 2) {
					my %pars = $ped->get_parents($famid);
					if (all { @$_ == 2 } values %pars) {
						return "quad";
					}
					else {
						return "extended";
					}
				}
				else {
					return "extended";
				}
			}
		}
		else {
			return "multi-generation";
		}
	}
}


sub _hist2desc {
	my (%hist) = @_;
	my @list;
	foreach my $type (keys %hist) {
		my $count = $hist{$type};
		my $item = $type;
		if ($hist{$type} > 1) {
			my @elem = split(/\s+/, $type);
			if ($elem[0] =~ /^Parent/) {
				$elem[0] = "Parents(?Sex)";
			}
			elsif ($elem[0] =~ /^Child/) {
				$elem[0] = "Children(?Sex)";
			}
			elsif ($elem[0] =~ /^Unknown/) {
				$elem[1] = $elem[1]."s";
			}
			$item = join(" ", @elem);
		}
		if ($count > 1) {
			if (defined $plural{$item}) {
				push @list, "$count $plural{$item}";
			}
			else {
				push @list, "$count ${item}s";
			}	
		}
		else {
			push @list, "$count $item";
		}
	}
	return join(", ", @list);
}

sub _po_sex {
	my ($par_gender, $chd_gender) = @_;
	my %par_hist = array2hist(@$par_gender);
	my %chd_hist = array2hist(@$chd_gender);

	if (@$par_gender == 1) {
		if (@$chd_gender == 1) {
			return $par_gender->[0]."-".$chd_gender->[0];
		}
		else {
			return $par_gender->[0]."-"._hist2desc(%chd_hist);
		}
	}
	elsif (@$par_gender == 2) {
		if (@$chd_gender == 1) {
			return "Both Parents and 1 ".$chd_gender->[0];
		}
		else {
			return "Both Parents and "._hist2desc(%chd_hist);
		}
	}
	else {
		return _hist2desc(%par_hist)." and "._hist2desc(%chd_hist);
	}
}

sub _par_aff {
	my ($paraff, $parsex) = @_;
	if (all { $paraff->[$_] eq $paraff->[0] } 0..$#$parsex ) {
		if (@$parsex == 1) {
			return $PARSEX{$parsex->[0]}." ".lcfirst($AFF{$paraff->[0]});
		}
		elsif (@$parsex == 2) {
			return "Both parents ".lcfirst($AFF{$paraff->[0]});
		}
		else {
			return "All parents ".lcfirst($AFF{$paraff->[0]});
		}
		
	}
	else {
		return join(", ", map { $PARSEX{$parsex->[$_]}." ".lcfirst($AFF{$paraff->[$_]}) } 0..$#$parsex);
	}
}

sub _chd_aff {
	my ($chdaff, $chdsex) = @_;
	if (all { $chdaff->[$_] eq $chdaff->[0] } 0..$#$chdsex) {
		if (@$chdsex == 1) {
			return $CHDSEX{$chdsex->[0]}." ".lcfirst($AFF{$chdaff->[0]});
		}
		elsif (@$chdsex == 2) {
			return "Both children ".lcfirst($AFF{$chdaff->[0]});
		}
		else {
			return "All children ".lcfirst($AFF{$chdaff->[0]})
		}
	}
	else {
		my %hist;
		for(my $ii = 0; $ii < @$chdaff; $ii ++) {
			my $key = $CHDSEX{$chdsex->[$ii]}." ".lcfirst($AFF{$chdaff->[$ii]});
			$hist{$key} ++;
		}
		return _hist2desc(%hist);
	}
}

sub list_fams {
	my ($ped, $famids, $type, $pref) = @_;
	$pref = " - " unless defined $pref;

	foreach my $famid (@$famids) {
		my ($fid, $ftype);
		if (defined $type) {
			$fid = $famid;
			$ftype = $type;
		}
		else {
			$ftype = fam_type($ped, $famid);
			$fid = $famid." (".ucfirst($ftype).")";

		}

		my @founders = $ped->get_founders($famid);
		my @parent_sex = map { $ped->get_sex($famid, $_) // 0 } @founders;
		my @parent_gender = map { $PARSEX{$_} } @parent_sex;
		my @parent_aff = map { $ped->get_aff($famid, $_) // 0 } @founders;
		my @parent_pheno = map { $AFF{$_} } @parent_aff;
		
		my @offspring = $ped->get_nonfounders($famid);
		# We will order offspring by their phenotypes, affected come first
		my %child_pheno = map { my $phe = $ped->get_aff($famid, $_) // 0; ($_ => $phe) } @offspring;
		#@offspring = sort { $child_pheno{$b} <=> $child_pheno{$a} } @offspring;

		my @child_sex = map { $ped->get_sex($famid, $_) // 0 } @offspring;
		my @child_gender = map { $CHDSEX{$_} } @child_sex;
		my @child_aff = map { $child_pheno{$_} } @offspring;
		my @child_pheno = map { $AFF{$_} } @child_aff;

		if ($ftype eq 'single') {
			printf $fout $pref."$fid : %s; %s\n", $SEX{$parent_sex[0]}, $AFF{$parent_aff[0]};
		}
		elsif ($ftype eq 'duo' || $ftype eq 'single-parent' || $ftype eq 'trio' || 
			$ftype eq 'quad' || $ftype eq 'multi-marriage' || $ftype eq 'extended') {
			printf $fout $pref."$fid : %s; %s; %s\n", _po_sex(\@parent_gender, \@child_gender),
				_par_aff(\@parent_aff, \@parent_sex), _chd_aff(\@child_aff, \@child_sex);
		}
		elsif ($ftype eq 'multi-generation') {
			my %sex_hist = array2hist(map { $SEX{$_} } @parent_sex, @child_sex);
			my %phe_hist = array2hist(@parent_pheno, @child_pheno);
			printf $fout $pref."$fid : Family size = %d, Number of generation = %d; %s; %s\n",
				$ped->get_famsize($famid), $ped->get_numgen($famid),
				_hist2desc(%sex_hist), _hist2desc(%phe_hist);
		}
		else {
			croak "Unsupported family type: $type";
		}
	}
}

sub _sum_tab {
	my ($t) = @_;
	my @fields = qw|CHDAFF PARAFF POAFF|;
	my @g;
	$g[0] = $t->group(['TYPE'], ['FAMID'], [sub {scalar @_}], ['All'], 0);
	foreach my $field (@fields) {
		push @g, $t->group(['TYPE', $field], ['FAMID'], [sub {scalar @_}], ['Count'], 0)->
			cast(['TYPE'], $field, Data::Table::STRING, 'Count', sub {$_[0]});
	}
	# Join tables ?
	# Depending on the number of rows in the intermediate tables
	my $tab = $g[0]->join($g[1], 0, ['TYPE'], ['TYPE'])->join($g[2], 0, ['TYPE'], ['TYPE']);
	if ($g[1]->nofCol > 2 && $g[2]->nofCol > 2) {
		$tab = $tab->join($g[3], 0, ['TYPE'], ['TYPE']);
	}
	$tab->rename("TYPE", "Family Type");
	#$tab->col(2);
	#add last raw for total
	my @tot = "Total";
	foreach my $ii (1..$tab->lastCol) {
		push @tot, sum(grep { defined $_ } $tab->col($ii));
	}
	$tab->addRow(\@tot, undef, {addNewCol => 1});
	return $tab;

}

sub tab_singles {
	my ($ped, $fids) = @_;
	my $t = Data::Table->new();
	foreach my $fid (@$fids) { 
		my @samp = $ped->get_members($fid);
		my $sampsex = $ped->get_sex($fid, $samp[0])//0;
		my $sampaff = $ped->get_aff($fid, $samp[0])//0;
		$t->addRow({ FAMID => $fid, TYPE => $SEX{$sampsex}, AFF => $AFF{$sampaff} });
	}
	my @g;
	$g[0] = $t->group(['TYPE'], ['FAMID'], [sub {scalar @_}], ['All'], 0);
	$g[1] = $t->group(['TYPE', 'AFF'], ['FAMID'], [sub {scalar @_}], ['Count'], 0)->
		cast(['TYPE'], 'AFF', Data::Table::STRING, 'Count', sub {$_[0]});

	my $tab = $g[0]->join($g[1], 0, ['TYPE'], ['TYPE']);
	$tab->rename('TYPE', "Gender");
	my @tot = "Total";
	foreach my $ii (1..$tab->lastCol) {
		push @tot, sum(grep { defined $_ } $tab->col($ii));
	}
	$tab->addRow(\@tot, undef, {addNewCol => 1});
	return $tab;
}

sub _child_pheno {
	my ($chdphe, $chdphe2);
	if (@_ == 1) {
		my ($chdaff) = @_;
		$chdphe = "Child\n".$AFF{$chdaff};
		$chdphe2 = "Chld ".$AF{$chdaff};
	}
	elsif (@_ >= 2) {
		my @sibaff = @_;
		my $sibaff = join('', @sibaff);
		my $all = @_ == 2 ? "Both" : "All";
		$chdphe = $sibaff eq '0' x @_ ? "$all Sibs\nUnknown" :
			$sibaff eq '1' x @_ ? "$all Sibs\nUnaffected" :
			$sibaff eq '2' x @_ ? "$all Sibs\nAffected" :
			@_ == 2 ? "1st Sib ".$AFF{$sibaff[0]}."\n"."2nd Sib ".$AFF{$sibaff[1]} :
			sprintf("%d Affected,\n %d Unaffected", scalar(grep { $_ eq '2'} @sibaff),
				scalar(grep { $_ eq '1'} @sibaff));
		$chdphe2 = $sibaff eq '0' x @_ ? "$all Sib NA" :
			$sibaff eq '1' x @_ ? "$all Sib Unaff" :
			$sibaff eq '2' x @_ ? "$all Sib Aff" :
			@_ == 2 ? "Sib1 ".$AF{$sibaff[0]}." ,Sib2 ".$AF{$sibaff[1]} : 
			sprintf("Aff=%d, Unaff=%d", scalar(grep { $_ eq '2'} @sibaff),
				scalar(grep { $_ eq '1' } @sibaff));
	}
	return ($chdphe, $chdphe2);
}


sub _parent_pheno {
	my ($parphe, $parphe2);
	if (@_ == 1) {
		my ($paraff) = @_;
		$parphe = "Parent\n".$AFF{$paraff};
		$parphe2 = "Par ".$AF{$paraff};
	}
	elsif (@_ == 2) {
		my ($dadaff, $momaff) = @_;
		my $paraff = join('', @_);
		$parphe = $paraff eq '00' ? "Both Parents\nUnknown" :
			$paraff eq '11' ? "Both Parents\nUnaffected" :
			$paraff eq '22' ? "Both Parents\nAffected" :
			"Father ".$AFF{$dadaff}.",\nMother ".$AFF{$momaff};
		$parphe2 = $paraff eq '00' ? "Both Par NA" :
			$paraff eq '11' ? "Both Par Unaff" : $paraff eq '22' ? "Both Par Aff" :
			"Fa ".$AF{$dadaff}.", Mo ".$AF{$momaff};
	}
	return ($parphe, $parphe2);
}

sub tab_duos {
	my ($ped, $fids) = @_;
	my $t = Data::Table->new();
	foreach my $fid (@$fids) {
		my @parent = $ped->get_founders($fid); 
		my @child = $ped->get_nonfounders($fid);
		my $parsex = $ped->get_sex($fid, $parent[0])//0;
		my $chdsex = $ped->get_sex($fid, $child[0])//0;
		
		my ($chdphe, $chdphe2) = _child_pheno($ped->get_aff($fid, $child[0])//0);
		my ($parphe, $parphe2) = _parent_pheno($ped->get_aff($fid, $parent[0])//0);

		$t->addRow({ FAMID => $fid, TYPE => $POSEX{"$parsex$chdsex"}, 
			CHDAFF => $chdphe,  PARAFF => $parphe, POAFF =>  $chdphe2."\n".$parphe2 });
	}
	return _sum_tab($t);
}

sub tab_trios {
	my ($ped, $fids) = @_;
	my $t = Data::Table->new();
	foreach my $fid (@$fids) {
		my @child = $ped->get_nonfounders($fid);
		my $chdsex = $ped->get_sex($fid, $child[0])//0;
		my @parents = $ped->get_parents($fid, $child[0]);
		
		my ($chdphe, $chdphe2) = _child_pheno($ped->get_aff($fid, $child[0])//0);
		my ($parphe, $parphe2) = _parent_pheno($ped->get_aff($fid, $parents[0])//0,
											   $ped->get_aff($fid, $parents[1])//0);

		$t->addRow({ FAMID => $fid,	TYPE => "Parents-".$CHDSEX{$chdsex}, 
			CHDAFF => $chdphe,  PARAFF => $parphe, POAFF =>  $chdphe2."\n".$parphe2 });
	}
	return _sum_tab($t);
}

sub tab_quads {
	my ($ped, $fids) = @_;
	my %SIBSEX = ('11' => 'Brothers', '12' => 'Brother-Sister', '10' => 'Brother-Sib(?Sex)',
		'21' => 'Sister-Brother', '22' => 'Sisters', '20' => 'Sister-Sib(?Sex)', 
		'00' => 'Sib-Pair(?Sex)');

	my $t = Data::Table->new();
	foreach my $fid (@$fids) {
		my @sibs = sort { $ped->get_sex($fid, $b) // 0 <=> $ped->get_sex($fid, $a) // 0 }
			$ped->get_nonfounders($fid);
		my @sibaff = map { $ped->get_aff($fid, $_) // 0 } @sibs;	
		my $sibsex = join('', map { $ped->get_sex($fid, $_) // 0 } @sibs);
		my @parents = $ped->get_parents($fid, $sibs[0]);

		my ($chdphe, $chdphe2) = _child_pheno(@sibaff);
		my ($parphe, $parphe2) = _parent_pheno($ped->get_aff($fid, $parents[0]) // 0,
											   $ped->get_aff($fid, $parents[1]) // 0);

		$t->addRow({ FAMID => $fid, TYPE => $SIBSEX{$sibsex}, 
			CHDAFF => $chdphe, PARAFF => $parphe, POAFF =>  $chdphe2."\n".$parphe2 });
	}
	return _sum_tab($t);
}

sub print_txt_tab {
	my ($t) = @_;
	$fout = \*STDOUT unless defined $fout;
	my $tb = Text::Table->new($t->header);
	$tb->load(@{$t->{data}});
	my $nline = @{$t->{data}};
	print $fout $tb->rule('-');
	print $fout $tb->title;
	print $fout $tb->rule('-');
	print $fout $tb->body(0, $nline-1);
	print $fout $tb->rule('-');
	print $fout $tb->body($nline-1);
}

__END__

=head1 NAME

ped_stats - Summarize pedigrees defined in ped file.

=head1 DESCRIPTION

The scripts cross tabulates the counts of sex and affection status in different
types of nuclear families (duo, trio, and quad) and single individuals. Other types 
of families will be listed.

This script is mainly designed for trio or quad based sequencing studies. It has
limited function to summarize complex pedigrees. For the latter purpose, C<pedstat>
in C<merlin> software package is recommended.


=head1 REQUIRED ARGUMENTS

=over

=item -[-]p[ed] [=] <file> 

Input pedigree file.

=for Euclid:
  file.type: readable

=back

=head1 OPTIONS

=over

=item -[-]out[put] [=] <file>

Output file name. Default is to write summaries to STDOUT.

=for Euclid:
  file.type: writable

=item -[-]xlsx

Also output summary table to an excel file (prefix.xlsx).

=item -[-]split[-]fam

When one family has more than one connected components, switching on this option
will split them into different families when counting families.

=item -[-]famsep [=] <separator>

For families with multiple connected components, use this separator to connect FID and suffix.

=for Euclid:
  separator.default: '.'

=item -[-]ignore [=] <pattern>

Pattern for sample names to be ignored from analysis.

=item -[-]cutoff [=] <thres>

Sample size cutoff. For family types with counts no greater than this, output a brief list;
otherwise output a table of counts.

=for Euclid:
  thres.default: 5

=back

=cut




