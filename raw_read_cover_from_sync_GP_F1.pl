#!/usr/bin/perl
use strict; use warnings;

# Output from this script is used to plot the raw read count data for the two GP and the F1 individuals with R/ggplot2
# input: an annotation file for the five samples, the popoolation2 mpileup2sync file for these samples (in the same order as in the annotation file), a list of loci to analyse


my $FIRSTDATACOL = 3;     # first sample data in sync file is in column with index 3
my $MINCOVER1 = 2;        # small plot symbol (likely contamination) if minor allele read count < $MINCOVER1 and max allele read count  < $HIGHCOUNT
my $MINCOVER2 = 7; 	  # small plot symbol (likely contamination) if minor allele read count < $MINCOVER2 and max allele read count  >= $HIGHCOUNT
my $HIGHCOUNT = 100;
my $SMALLDOT = 0;
my $LARGEDOT = 1;


my $indfile = "gp_f1_inds_annotated.txt";        # file with sample annotation
my $indf;
my $syncfile = "gp_f1_locus5568.sync";           # popoolation2 mpileup2sync output file with data for the two grandparents and the three F1 individuals
my $syncf;
my $outfile = "gp_f1_for_plotting.txt";          # output file
my $out;
my $locusfile = "locus_list.txt";                # list of loci to analyse (contig numbers, one per line)
my $lf;

# the six sequence state in the order in which they appear in the sync file, D = deletion
my @states = ('A','T','C','G','N','D');

my $line;
my @array;
my @rawcounts;
my @basecount;
my @samples;
my @taxon;
my @sex;
my @generation;
my $i;
my $j;
my $k;
my $mincover;
my %selectedloci;
my $locus;

open ($indf,"< $indfile") or die "Can't open $indfile\n";
while ($line = <$indf>) {
	chomp($line);	
	@array = split ('\t',$line);
	push @samples, $array[0];
	push @taxon, $array[1];
	push @sex, $array[2];
	push @generation, $array[3];
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
print $out "locus\tposition\tsample\ttaxon\tsex\tgeneration\tcoverage\tstate\tsymbolsize\n";
while ($line = <$syncf>) {
	($locus) = $line =~ /contig(\d+)_/;
	if (exists $selectedloci{$locus}) {
		chomp($line);	
		@array = split ('\t',$line);
		for ($i=$FIRSTDATACOL;$i<@array;$i++) {
			$j = $i - $FIRSTDATACOL;   # index of sample in @samples
			@rawcounts = split (':',$array[$i]);
			@basecount = ();
			for ($k=0;$k<@rawcounts;$k++) {
				push @basecount, { count => $rawcounts[$k], state => $states[$k] }
			}
			@basecount = sort { $b->{count} <=> $a->{count} } @basecount;
			# output N (NULL) if the 'most common' state has read count of 0
			if ($basecount[0]->{count} == 0) {
				print $out "$locus\t$array[1]\t$samples[$j]\t$taxon[$j]\t$sex[$j]\t$generation[$j]\t0\tN\t$SMALLDOT\n";
				next;
			}
			# set threshold for read counts that are plotted with a small symbol (likely contamination)
			if ($basecount[0]->{count} > $HIGHCOUNT) { $mincover = $MINCOVER2 }
			else                                     { $mincover = $MINCOVER1 }
			if ($basecount[1]->{count} > $mincover) {            # heterozygote
				print $out "$locus\t$array[1]\t$samples[$j]\t$taxon[$j]\t$sex[$j]\t$generation[$j]\t$basecount[0]->{count}\t";
				# determine whether the state with highest read count is REF (plotted as 'R' with a small symbol)
				if ($basecount[0]->{state} eq $array[2]) {   print $out "R\t$SMALLDOT\n" }
				else                                     {   print $out "$basecount[0]->{state}\t$LARGEDOT\n" }
				print $out "$locus\t$array[1]\t$samples[$j]\t$taxon[$j]\t$sex[$j]\t$generation[$j]\t$basecount[1]->{count}\t";
				# determine whether the state with the 2nd highest read count is REF
				if ($basecount[1]->{state} eq $array[2]) {   print $out "R\t$SMALLDOT\n" }
				else                                     {   print $out "$basecount[1]->{state}\t$LARGEDOT\n" }
			}
			else {
				# state with highest read count
				print $out "$locus\t$array[1]\t$samples[$j]\t$taxon[$j]\t$sex[$j]\t$generation[$j]\t$basecount[0]->{count}\t";
				if ($basecount[0]->{state} eq $array[2]) {   print $out "R\t$SMALLDOT\n" }
				else                                     {   print $out "$basecount[0]->{state}\t$LARGEDOT\n" }
				# state with 2nd highest read count (smaller than $mincover) will be plotted with a small symbol
				if ($basecount[1]->{count} > 0) {
					print $out "$locus\t$array[1]\t$samples[$j]\t$taxon[$j]\t$sex[$j]\t$generation[$j]\t$basecount[1]->{count}\t$basecount[1]->{state}\t$SMALLDOT\n";
				}
			}
			# plot any third state with non-zero read count
			if ($basecount[2]->{count} > $mincover) {
				print $out "$locus\t$array[1]\t$samples[$j]\t$taxon[$j]\t$sex[$j]\t$generation[$j]\t$basecount[2]->{count}\t$basecount[2]->{state}\t$LARGEDOT\n";
			}
			elsif ($basecount[2]->{count} > 0) {
				print $out "$locus\t$array[1]\t$samples[$j]\t$taxon[$j]\t$sex[$j]\t$generation[$j]\t$basecount[2]->{count}\t$basecount[2]->{state}\t$SMALLDOT\n";
			}
		}
	}
}

close ($syncf);
close ($out);
