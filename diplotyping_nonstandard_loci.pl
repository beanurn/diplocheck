#!/usr/bin/perl
use warnings; use strict;

my $FIRSTDATACOL = 3;
my $MINCOVER = 5;

die "Please, enter the dataset to be analysed (gp_f1... or fam?_set?)\n" unless (@ARGV == 1);
my ($set) = @ARGV;

my $varfile = "gp_variants.txt";
my $vf;
my $popsyncfile = "$set.sync";
my $psf;

my $keyfile = "$set.key";
my $outfile = "${set}_rescored.txt";
my $kf;
my $out;

my $line;
my @array;
my @basecount;
my @rawcounts;
my %markercount;
my %markertype;
my %markerstate;
my %typecounts;
my %sampledata;
my @samples;
my $sample;
my @string;
my @diplo;
my @missing;
my $i;
my $j;
my $k;
my $locus;
my $first;
my $current_locus;
my $pos;
my $count;
my $key;

my @states = ('A','T','C','G','N','D');

open ($vf,"< $varfile") or die "Can't open $varfile\n";
while ($line = <$vf>) {
	chomp($line);
	if ($line =~ /^\d/) {
		@array = split ('\t',$line);
		if (@array > 6 && length $array[6] > 0) {
			unless (exists $markertype{$array[0]}) {
				my %positions;
				$markertype{$array[0]} = \%positions;
				my %counts;
				$typecounts{$array[0]} = \%counts;
				my %mstate;
				$markerstate{$array[0]} = \%mstate;
			}
			$markertype{$array[0]}->{$array[1]} = $array[6];
			$typecounts{$array[0]}->{$array[6]} += 1;
			$markercount{$array[0]}++;
			$markerstate{$array[0]}->{$array[1]} = $array[4];
		}
	}
}
close ($vf);

if ($set =~ /gp_f1/) {
	push @samples, 'BvGP';
	push @samples, 'BbGP';
	push @samples, 'F1F6';
	push @samples, 'F1M';
	push @samples, 'F1F7';
	for ($i=0;$i<@samples;$i++) {
		push @string, '';
		push @diplo, '';
		push @missing, 0;
	}
}
else {
	open ($kf,"< $keyfile") or die "Can't open $keyfile\n";
	while ($line = <$kf>) {
		($sample) = $line =~ /(Fam[67]_\d+)/;
		push @samples, $sample;
		push @string, '';
		push @diplo, '';
		push @missing, 0;
	}
	close ($kf);
}

$first = 1;

open ($psf,"< $popsyncfile") or die "Can't open $popsyncfile\n";
open ($out,"> $outfile") or die "Can't open $outfile\n";
while ($line = <$psf>) {
	($locus) = $line =~ /contig(\d+)/;
	if (exists $markerstate{$locus}) {
		if ($first) { $current_locus = $locus; $first = 0 }
		if ($locus != $current_locus) {
			assign_diplotypes($current_locus);
			$current_locus = $locus;
			foreach $i (@string) { $i = '' }
			foreach $i (@missing) { $i = 0 }
		}
		chomp($line);
		@array = split ('\t',$line);
		if (exists $markerstate{$locus}->{$array[1]}) {
			$pos = $array[1];
			for ($i=$FIRSTDATACOL;$i<@array;$i++) {
				$j = $i - $FIRSTDATACOL;   # index of sample in @samples
				@rawcounts = split (':',$array[$i]);
				@basecount = ();
				for ($k=0;$k<@rawcounts;$k++) {
					push @basecount, { count => $rawcounts[$k], state => $states[$k] }
				}						
				@basecount = sort { $b->{count} <=> $a->{count} } @basecount;
#				if ($locus == 20496) { print "locus $locus, pos $pos, sample $samples[$j], 0: $basecount[0]->{count}, 1: $basecount[1]->{count}\n" }
				if ($basecount[1]->{count} > $MINCOVER) {
					if ($basecount[0]->{state} eq $markerstate{$locus}->{$pos} || $basecount[1]->{state} eq $markerstate{$locus}->{$pos}) {
#						print "locus $locus, position $pos, markerstate $markerstate{$locus}->{$pos}\n";
#						print "locus $locus, position $pos, markertype $markertype{$locus}->{$pos}\n";
						if ($string[$j] =~ /$markertype{$locus}->{$pos}/ && $string[$j] =~ /het/) {
							($count) = $string[$j] =~ /${markertype{$locus}->{$pos}}_het_(\d+)/;
							$count++;
							$string[$j] =~ s/${markertype{$locus}->{$pos}}_het_\d+/${markertype{$locus}->{$pos}}_het_$count/;
						}
						else { 	$string[$j] .= ":" . $markertype{$locus}->{$pos} . "_het_1" }
					}
				}
				elsif ($basecount[0]->{count} > $MINCOVER) {
#					if ($locus == 20496) { print "locus $locus, pos $pos, sample $samples[$j], 0th state > mincover, $markertype{$locus}->{$pos}\n" }
					if ($markertype{$locus}->{$pos} eq 'dom') { $string[$j] .= ':dom_pres' }
					else {
						if ($basecount[0]->{state} eq $markerstate{$locus}->{$pos}) {
							if ($string[$j] =~ /${markertype{$locus}->{$pos}}_hom/) {
								($count) = $string[$j] =~ /${markertype{$locus}->{$pos}}_hom_(\d+)/;
								$count++;
								$string[$j] =~ s/${markertype{$locus}->{$pos}}_hom_\d+/${markertype{$locus}->{$pos}}_hom_$count/;
							}
							else { 	$string[$j] .= ":" . $markertype{$locus}->{$pos} . "_hom_1" }
						}
					}
#					if ($locus == 22235) { print "locus $locus, pos $pos, sample $samples[$j], string $string[$j]\n" }
				}
				else { 
					if ($markertype{$locus}->{$pos} eq 'dom') { $string[$j] = 'dom_abs' }
					else                                      { $missing[$j]++ }
					# note that locus 22235 with Bb deletion and Bvv variant will have 1 missing score at v1 variant position when d/d
				}
			} # samples
		} # markerpositions
	} # markerloci
}

assign_diplotypes($current_locus);

close ($psf);
close ($out);

sub assign_diplotypes {

my ($loc) = @_;
my $i;
	
	for ($i=0;$i<@samples;$i++) {
		if    ($string[$i] eq 'dom_abs')     { $diplo[$i] = 'd/d' }
		elsif ($string[$i] =~ 'dom_pres') {                          # important not to use eq here, because of colon in :dom_pres
			if ($samples[$i] eq 'BvvGP') { $diplo[$i] = 'v/v' }
			elsif ($samples[$i] =~ /F1/) { $diplo[$i] = 'v/d' }
			else                         { $diplo[$i] = 'v/v|v/d' }
		} 
		elsif ($missing[$i] == $markercount{$loc}) { $diplo[$i] = 'NA' }
		else {
			if (length $string[$i] == 0) { $diplo[$i] = 'v/v' }
			elsif ($string[$i] =~ /Bb_hom/ && $string[$i] =~ /Bb_het/) { print "error - locus $loc, sample $samples[$i], string $string[$i]\n"; $diplo[$i] = 'error' }
			elsif ($string[$i] =~ /Bb_hom/) { 
				$diplo[$i] = 'b/b';
				if	($string[$i] =~ /b1_het/) { $diplo[$i] = 'b1/b2' }
				elsif	($string[$i] =~ /b1_hom/) { $diplo[$i] = 'b1/b1' }
				elsif   ($string[$i] =~ /b2_hom/) { $diplo[$i] = 'b2/b2' }
				elsif 	($string[$i] =~ /[zw]/) {
						$diplo[$i] = '';
						if ($string[$i] =~ /z/)	{ $diplo[$i] .= 'z' }
						if ($string[$i] =~ /w/)	{ $diplo[$i] .= 'w' }
						$diplo[$i] = $diplo[$i] . "/" . $diplo[$i];
				}
				elsif (exists $typecounts{$loc}->{b1}) { $diplo[$i] = 'b2/b2' }   # only b1 has distinguishing variant
			}
			elsif ($string[$i] =~ /Bb_het/) {
				$diplo[$i] = 'b/v';
				if	($string[$i] =~ /b1_het/) { $diplo[$i] = 'b1/v' }
				elsif	($string[$i] =~ /b2_het/) { $diplo[$i] = 'b2/v' }
				elsif 	($string[$i] =~ /[zw]/) {
						$diplo[$i] = '';
						if ($string[$i] =~ /z/)	{ $diplo[$i] .= 'z' }
						if ($string[$i] =~ /w/)	{ $diplo[$i] .= 'w' }
						$diplo[$i] .= '/v';
				}
				elsif (exists $typecounts{$loc}->{b1} && $string[$i] !~ /b1/) { $diplo[$i] = 'b2/v' }  # only b1 has variant, b2 differs from Bvv only by homdiffs
			}
			else {
				if	($string[$i] =~ /b1_hom/) { $diplo[$i] = 'b1/b1' }
				elsif   ($string[$i] =~ /b2_hom/) { $diplo[$i] = 'b2/b2' }
				elsif   ($string[$i] =~ /b1_het/ && $string[$i] =~ /b2_het/) { $diplo[$i] = 'b1/b2' }
				elsif   ($string[$i] =~ /b1_het/) { $diplo[$i] = 'b1/v' }
				elsif   ($string[$i] =~ /b2_het/) { $diplo[$i] = 'b2/v' }
			}
		}
	}
	if ($loc == 4146) { adjust_diplotypes_4146() }  # zw-locus,b_het_deletion
#	if ($loc == 5568) { adjust_diplotypes_5568() }  # 4 haplos
	if ($loc == 218180 || $loc == 6647) { adjust_diplotypes_218180_6647() }  # taxdiff_b_het_deletion
	if ($loc == 383456) { adjust_diplotypes_383456() }  # taxdiff_b_het_deletion, only F1M carries b haplotype
	if ($loc == 22235) { adjust_diplotypes_22235() }   # dominant_v1
	if ($loc == 24822 || $loc == 13428) { adjust_diplotypes_24822_13428() }   # mixed locus
	if ($loc == 840699) { adjust_diplotypes_840699() }   # mixed locus, v1
	if ($loc == 30421) { adjust_diplotypes_30421() }   # taxdiff_v_het_deletion
	if ($loc == 3143 || $loc == 126712) { adjust_diplotypes_3143() }   # taxdiff_v_het_deletion
	if ($loc == 4005) { adjust_diplotypes_4005() }   # taxdiff_v_het_deletion
	if ($loc == 4517 || $loc == 19758) { adjust_diplotypes_4517() }   # taxdiff_b_het_deletion
	if ($loc == 18751) { adjust_diplotypes_18751() }   # taxdiff_v1
	if ($loc == 3843) { adjust_diplotypes_3843() }   # taxdiff_v1
	if ($loc == 119758) { adjust_diplotypes_119758() }   # b1_v1 (v var fixed in Bb)
	if ($loc == 335939 || $loc == 21066) { adjust_diplotypes_335939($loc) }  # Bb het for deletion
	if ($loc == 211636 || $loc == 243058 || $loc == 403637 || $loc == 13968 || $loc == 1210 || $loc == 21292 || $loc == 163368 || $loc == 179810) { adjust_diplotypes_BvvGP_low_cover()  }
	if ($loc == 192445 || $loc == 1886999 || $loc == 23393 || $loc == 15684 || $loc == 844386 || $loc == 8087 || $loc == 20031) { adjust_diplotypes_BvvGP_low_cover()  }
	if ($loc == 26573 || $loc == 184052 || $loc == 3724211 || $loc == 1222927 || $loc == 524337 || $loc == 300778) { adjust_diplotypes_BvvGP_low_cover()  }

	for ($i=0;$i<@samples;$i++) {
		print $out "$loc\t$samples[$i]\t$diplo[$i]\t$missing[$i]\n";
	}
		
}


sub adjust_diplotypes_4146  {

my $s;

	for ($s=0;$s<@samples;$s++) {
		if ($samples[$s] =~ /Fam7/) {
			if ($diplo[$s] eq 'v/v')    { $diplo[$s] = 'v/v|v/d' }
			if ($diplo[$s] eq 'zw/zw')  { $diplo[$s] = 'zw/d' }
		}
		elsif ($samples[$s] eq 'F1F7')  { $diplo[$s] = 'v/d' }
		elsif ($samples[$s] eq 'BbGP')  { $diplo[$s] = 'zw/d' }
	}
}

sub adjust_diplotypes_5568  {

my $s;

	for ($s=0;$s<@samples;$s++) {
		if ($samples[$s] eq 'F1F6')          { $diplo[$s] = 'v2/b1' }
		elsif ($samples[$s] eq 'F1F7')       { $diplo[$s] = 'v2/b1' }
		elsif ($samples[$s] eq 'F1M')        { $diplo[$s] = 'v1/b2' }
		elsif ($string[$s] =~ /Bb_het/) {
			if ($string[$s] =~ /b1_het/) { $diplo[$s] = 'b1/v1' }
			else 			     { $diplo[$s] = 'b2/v2' }
		}
		elsif ($string[$s] =~ /v1_het/)      { $diplo[$s] = 'v1/v2' }
	}
}


sub adjust_diplotypes_218180_6647  {

my $i;

	for ($i=0;$i<@samples;$i++) {
		if 	($samples[$i] eq 'BvvGP') { $diplo[$i] = 'v/v' }
		elsif 	($samples[$i] eq 'BbGP')  { $diplo[$i] = 'b/d' }
		elsif 	($samples[$i] eq 'F1F6')   { $diplo[$i] = 'v/d' }
		elsif 	($samples[$i] eq 'F1F7')   { $diplo[$i] = 'v/d' }
		elsif 	($samples[$i] eq 'F1M')    { $diplo[$i] = 'b/v' }
		else {
			if ($diplo[$i] eq 'v/v') { $diplo[$i] = 'v/v|v/d' }
			if ($diplo[$i] eq 'b/b')  { $diplo[$i] = 'b/d' }
		}
	}
}


sub adjust_diplotypes_22235  {

my $i;
	# case of dom_abs is dealt with in assign diplotypes
	for ($i=0;$i<@samples;$i++) {
		if      ($string[$i] =~ /v1_het/)   { $diplo[$i] = 'v1/v2' }
		elsif   ($string[$i] =~ /v1_hom/)   { $diplo[$i] = 'v1/d' }
		elsif 	($string[$i] =~ /dom_pres/) { $diplo[$i] = 'v2/d' }
	}
}


sub adjust_diplotypes_24822_13428  {

my $i;
        # score by presence of absence of sex-linked w haplotype, w+Bb_het is ambiguous, because the v haplotype might originate from the second locus, z/v inds do occur
        for ($i=0;$i<@samples;$i++) {
                if ($string[$i] =~ /w/)  {
                        if ($string[$i] =~ /Bb_het/) { $diplo[$i] = 'b/v|b/b' }
                        else                         { $diplo[$i] = 'b/b' }
                }
                else                             { $diplo[$i] = 'v/v' }
        }
}


sub adjust_diplotypes_30421  {

my $i;
	# het deletion in BvvGP passed on to F1F7
	for ($i=0;$i<@samples;$i++) {
		if 	($samples[$i] eq 'BvvGP')        { $diplo[$i] = 'v/d' }
		elsif   ($samples[$i] eq 'F1F7')         { $diplo[$i] = 'b/d' }
		elsif   ($samples[$i] =~ /Fam7/) {
			if    ($diplo[$i] =~ /b\/b/)      { $diplo[$i] = 'b/b|b/d' }
			elsif ($diplo[$i] =~ /v\/v/)      { $diplo[$i] = 'v/d' }
		}
	}
}


sub adjust_diplotypes_335939  {

my ($locus) = @_;
my $i;
	# taxdiff_b_het_deletion: b haplotype passed on to F1F6 only
	for ($i=0;$i<@samples;$i++) {
		if 	($samples[$i] eq 'BbGP')        { $diplo[$i] = 'b/d' }
		elsif   ($samples[$i] eq 'F1F7')        { $diplo[$i] = 'v/d' }
		elsif   ($samples[$i] eq 'F1M')         { $diplo[$i] = 'v/d' }
		elsif   ($samples[$i] =~ /Fam6/) {
			if    ($diplo[$i] =~ /v\/v/)      { $diplo[$i] = 'v/v|v/d' }
			elsif ($diplo[$i] =~ /b\/b/)      { $diplo[$i] = 'b/d' }
		}
		elsif   ($samples[$i] =~ /Fam7/) {
			if    ($diplo[$i] =~ /v\/v/)      { $diplo[$i] = 'v/v|v/d' }
			elsif ($missing[$i] == $markercount{$locus}) { $diplo[$i] = 'd/d' }
		}
	}
}


sub adjust_diplotypes_3143  {

my $i;
	# taxdiff_v_het_deletion: F1M is b/-, f1 females are b/v. So, corrections apply to both families
	for ($i=0;$i<@samples;$i++) {
		if      ($samples[$i] eq 'BbGP')  { $diplo[$i] = 'b/b' }
		elsif   ($samples[$i] eq 'BvvGP') { $diplo[$i] = 'v/d' }
		elsif   ($samples[$i] eq 'F1M')   { $diplo[$i] = 'b/d' }
		elsif   ($diplo[$i] eq 'b/b')     { $diplo[$i] = 'b/b/|b/d' }
		elsif   ($diplo[$i] eq 'v/v')     { $diplo[$i] = 'v/d' }
	}
}

sub adjust_diplotypes_4005  {

my $i;
	# taxdiff_v_het_deletion: F1M is b/v, f1 females are b/-. So, corrections apply to both families
	for ($i=0;$i<@samples;$i++) {
		if      ($samples[$i] eq 'BvvGP') { $diplo[$i] = 'v/d' }
		elsif   ($samples[$i] eq 'F1F6')  { $diplo[$i] = 'b/d' }
		elsif   ($samples[$i] eq 'F1F7')  { $diplo[$i] = 'b/d' }
		elsif   ($diplo[$i] eq 'b/b')     { $diplo[$i] = 'b/b/|b/d' }
		elsif   ($diplo[$i] eq 'v/v')     { $diplo[$i] = 'v/d' }
	}
}


sub adjust_diplotypes_4517  {

my $i;
	# taxdiff_b_het_deletion: F1M is v/-, f1 females are b/v. So, corrections apply to both families
	for ($i=0;$i<@samples;$i++) {
		if      ($samples[$i] eq 'BbGP')  { $diplo[$i] = 'b/d' }
		elsif   ($samples[$i] eq 'F1M')   { $diplo[$i] = 'v/d' }
		elsif   ($diplo[$i] eq 'v/v')     { $diplo[$i] = 'v/v/|v/d' }
		elsif   ($diplo[$i] eq 'b/b')     { $diplo[$i] = 'b/d' }
	}
}


sub adjust_diplotypes_840699  {

my $i;
	# zw locus, v1 variant, scored on the basis of w/v1 in all F1s
	for ($i=0;$i<@samples;$i++) {
		if      ($samples[$i] eq 'BvvGP')  { $diplo[$i] = 'v1/v2' }
		elsif   ($samples[$i] eq 'F1F6')   { $diplo[$i] = 'v1/w' }
		elsif   ($samples[$i] eq 'F1F7')   { $diplo[$i] = 'v1/w' }
		elsif   ($samples[$i] eq 'F1M' )   { $diplo[$i] = 'v1/w' }
		elsif   ($string[$i] =~ /w/ && $string[$i] =~ /v/) { $diplo[$i] = 'v1/w' }
		elsif   ($string[$i] =~ /v/)                       { $diplo[$i] = 'v1/v1' }
		elsif   ($string[$i] =~ /w/)                       { $diplo[$i] = 'w/w' }
	}
}


sub adjust_diplotypes_383456  {

my $i;
	# taxdiff_b_het_deletion: only F1M carries b haplotype
	for ($i=0;$i<@samples;$i++) {
		if      ($samples[$i] eq 'BbGP')  { $diplo[$i] = 'b/d' }
		elsif   ($samples[$i] eq 'F1F6')  { $diplo[$i] = 'v/d' }
		elsif   ($samples[$i] eq 'F1F7')  { $diplo[$i] = 'v/d' }
		elsif   ($diplo[$i] eq 'v/v')     { $diplo[$i] = 'v/v/|v/d' }
		elsif   ($diplo[$i] eq 'b/b')     { $diplo[$i] = 'b/d' }
	}
}

sub adjust_diplotypes_119758  {

my $i;
	# b1 vars and a varpos that is het in both GP ('v1')
	for ($i=0;$i<@samples;$i++) {
		if      ($samples[$i] eq 'BbGP')  { $diplo[$i] = 'b1/b2' }
		elsif   ($samples[$i] eq 'BvvGP') { $diplo[$i] = 'v1/v2' }
		elsif   ($samples[$i] eq 'F1F6')  { $diplo[$i] = 'b1/v2' }
		elsif   ($samples[$i] eq 'F1F7')  { $diplo[$i] = 'b2/v2' }
		elsif   ($samples[$i] eq 'F1M')   { $diplo[$i] = 'b1/v1' }
		elsif   ($samples[$i] =~ /Fam6/) {
			if    ($string[$i] =~ /b1_hom/)                            { $diplo[$i] = 'b1/b1' }
			elsif ($string[$i] =~ /b1_het/ && $string[$i] =~ /v1_het/) { $diplo[$i] = 'b1/v1' }
			elsif ($string[$i] =~ /b1_het/ && $string[$i] !~ /v1/)     { $diplo[$i] = 'b1/v2' }
			elsif ($string[$i] !~ /b1/     && $string[$i] =~ /v1_het/) { $diplo[$i] = 'v1/v2' }
		}
		elsif   ($samples[$i] =~ /Fam7/) {
			if    ($string[$i] =~ /b1_het/ && $string[$i] =~ /v1_het/) { $diplo[$i] = 'b1/b2' }
			elsif ($string[$i] =~ /b1_het/ && $string[$i] !~ /v1/)     { $diplo[$i] = 'b1/v2' }
			elsif ($string[$i] =~ /v1_hom/)                            { $diplo[$i] = 'b2/v1' }
			elsif ($string[$i] !~ /b1/     && $string[$i] =~ /v1_het/) { $diplo[$i] = 'v1/v2' }
		}
	}
}


sub adjust_diplotypes_18751  {

my $i;
	# taxdiff_v1
	for ($i=0;$i<@samples;$i++) {
		if      ($samples[$i] eq 'BvvGP') { $diplo[$i] = 'v1/v2' }
		elsif   ($samples[$i] eq 'F1F6')  { $diplo[$i] = 'v1/b' }
		elsif   ($samples[$i] eq 'F1F7')  { $diplo[$i] = 'v2/b' }
		elsif   ($samples[$i] eq 'F1M')   { $diplo[$i] = 'v1/b' }
		elsif   ($string[$i] =~ /v1_hom/) { $diplo[$i] = 'v1/v1' }
		elsif   ($string[$i] =~ /Bb_het/ && $string[$i] =~ /v1_het/)     { $diplo[$i] = 'b/v1' }
		elsif   ($string[$i] =~ /Bb_het/)                                { $diplo[$i] = 'b/v2' }
		elsif   ($string[$i] =~ /v1_het/)                                { $diplo[$i] = 'v1/v2' }
	}
}

sub adjust_diplotypes_3843  {

my $i;
	# taxdiff_v1
	for ($i=0;$i<@samples;$i++) {
		if      ($samples[$i] eq 'BvvGP') { $diplo[$i] = 'v1/v2' }
		elsif   ($samples[$i] eq 'BbGP')  { $diplo[$i] = 'b/b' }
		elsif   ($samples[$i] eq 'F1F6')  { $diplo[$i] = 'b/v2' }
		elsif   ($samples[$i] eq 'F1F7')  { $diplo[$i] = 'b/v1' }
		elsif   ($samples[$i] eq 'F1M')   { $diplo[$i] = 'b/v1' }
		elsif   ($string[$i] =~ /v1_hom/) { $diplo[$i] = 'v1/v1' }
		elsif   ($string[$i] =~ /Bb_het/ && $string[$i] =~ /v1_hom/)     { $diplo[$i] = 'b/v1' }
		elsif   ($string[$i] =~ /Bb_het/ && $string[$i] =~ /v1_het/)     { $diplo[$i] = 'b/v2' }
		elsif   ($string[$i] =~ /v1_het/)                                { $diplo[$i] = 'v1/v2' }
		elsif   ($string[$i] =~ /v1_hom/)                                { $diplo[$i] = 'v1/v1' }
	}
}


sub adjust_diplotypes_BvvGP_low_cover {

my $i;

	for ($i=0;$i<@samples;$i++) {
		if ($samples[$i] eq 'BvvGP')  { $diplo[$i] = 'v/v' }
	}
}
	
