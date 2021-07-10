#!/usr/bin/perl
use strict; use warnings;

# Script produces a tabular data file with a list of variant (non-reference) states and their read counts in the two grandparents.
# Variants used for rescoring are manually recorded in this file.
# input: an annotation file for the GP and F1 samples, the popoolation2 mpileup2sync file for these samples (in the same order as in the annotation file), a list of loci to analyse

my $FIRSTDATACOL = 3;     # First sample data in in sync file in column with index 3
my $MINCOVER = 10;        # Only variants with read counts > $MINCOVER will be written to the output file


my $indfile = "gp_f1_inds_annotated.txt";
my $indf;
my $syncfile = "gp_f1_locus5568.sync";
my $syncf;
my $outfile = "gp_variants.txt";
my $out;
my $locusfile = "locus_list.txt";
my $lf;


# the six sequence state in the order in which they appear in the sync file, D = deletion
my @states = ('A','T','C','G','N','D');

my $line;
my @array;
my @rawcounts;
my @basecount;
my @samples;
my $i;
my $j;
my $k;
my $mincover;
my %selectedloci;
my $locus;
my $current_locus;
my $key;

open ($indf,"< $indfile") or die "Can't open $indfile\n";
while ($line = <$indf>) {
	chomp($line);	
	@array = split ('\t',$line);
	push @samples, $array[0];
}
close ($indf);

open ($lf,"< $locusfile") or die "Can't open $locusfile\n";
while ($line = <$lf>) {
	chomp($line);
	$selectedloci{$line} = 1;
}
close ($lf);

open ($syncf,"< $syncfile") or die "Can't open $syncfile\n";
open ($out,"> $outfile") or die "Can't open $outfile\n";
$current_locus = -1;
while ($line = <$syncf>) {
	($locus) = $line =~ /contig(\d+)_/;
	if (exists $selectedloci{$locus}) {
		if ($current_locus != $locus) {
			print $out "\n\nLOCUS $locus\n";
			print $out "locus\tposition\tsample\ttopdown_order\tstate\tcoverage\ttype\n";
			$current_locus = $locus;
		}
		chomp($line);	
		@array = split ('\t',$line);
		# here consider only the GP data
		for ($i=$FIRSTDATACOL;$i<$FIRSTDATACOL+2;$i++) {
			$j = $i - $FIRSTDATACOL;   # index of sample in @samples
			@rawcounts = split (':',$array[$i]);
			@basecount = ();
			for ($k=0;$k<@rawcounts;$k++) {
				push @basecount, { count => $rawcounts[$k], state => $states[$k] }
			}
			# write only non-ref states with read counts > $MINCOVER to the output file
			@basecount = sort { $b->{count} <=> $a->{count} } @basecount;
			if ($basecount[0]->{count} > $MINCOVER && $basecount[0]->{state} ne $array[2]) {
				print $out "$locus\t$array[1]\t$samples[$j]\t1\t$basecount[0]->{state}\t$basecount[0]->{count}\n";
			}
			if ($basecount[1]->{count} > $MINCOVER && $basecount[1]->{state} ne $array[2]) {           
				print $out "$locus\t$array[1]\t$samples[$j]\t2\t$basecount[1]->{state}\t$basecount[1]->{count}\n";
			}
			if ($basecount[2]->{count} > $MINCOVER && $basecount[2]->{state} ne $array[2]) { 
				print $out "$locus\t$array[1]\t$samples[$j]\t3\t$basecount[2]->{state}\t$basecount[2]->{count}\n";
			}
		}
	}
}

close ($syncf);
close ($out);
