#!/usr/bin/perl
use strict; use warnings;

my $testlocus = -1;

die "Please, provide an input file\n" unless (@ARGV == 1);
my ($infile) = @ARGV;

my $in;
my $outfile = $infile;
$outfile =~ s/\.txt/_lepmap.txt/;
my $out;
my $errfile = $infile;
$errfile =~ s/\.txt/_lepmap_error.txt/;
my $ef;
my $logfile = $infile;
$logfile =~ s/\.txt/_lepmap_log.txt/;
my $log;

my $line;
my @array;
my @harray;
my %states;
my @ncodes;
my @haplos;
my @altstates;
my $alt;
my $h;
my $i;
my %onecodes;
my $locus;
my $code;
my $mflag;
my $key;

my @nucleos = ('A','C','G','T');
my @diplos = ('AA','AC','AG','AT','CC','CG','CT','GG','GT','TT');
my $missing = '1 1 1 1 1 1 1 1 1 1';

open ($in,"< $infile") or die "Can't open $infile\n";
while ($line = <$in>) {
	chomp($line);
	@array = split ('\t',$line);
	unless (exists $states{$array[0]}) {
		my %statelist;
		$states{$array[0]} = \%statelist;
	}
	unless ($array[2] eq 'error') {
		@altstates = split ('\|',$array[2]);
		if ($array[0] == $testlocus) { print "$array[2]\t" }
		foreach $alt (@altstates) {
			@haplos = split ('\/',$alt);
			foreach $h (@haplos) { $states{$array[0]}->{$h} = 1 }
		}
	}	
}
close ($in);

if ($testlocus != -1) { print "states: " }
foreach $key (keys %{$states{$testlocus}}) {  print "$key\t" }

open ($log,"> $logfile") or die "Can't open $logfile\n";
foreach $locus (keys %states) {
	print $log "$locus\n";
	@harray = ();
	foreach $h (keys %{$states{$locus}}) {
		unless (length $h == 0 || $h eq 'NA') {
			# discard the odd mis-scoring of a state as b or v at a locus where distinct b or v haplotypes have been identified
			if (($h eq 'b' && exists $states{$locus}->{b1}) || ($h eq 'v' && exists $states{$locus}->{v1})) { next }
			push @harray, $h;
			print $log "$h ";
		}
	}
	print $log "\n";
	@harray = sort { $a cmp $b } @harray;
	for ($i=0;$i<@harray;$i++) {
		$states{$locus}->{$harray[$i]} = $nucleos[$i];
		print $log "$harray[$i]-$nucleos[$i] ";
	}
	print $log "\n";
}
close ($log);

open ($in,"< $infile") or die "Can't open $infile\n";
open ($out,"> $outfile") or die "Can't open $outfile\n";
open ($ef,"> $errfile") or die "Can't open $errfile\n";
while ($line = <$in>) {
	chomp($line);
	@array = split ('\t',$line);
	if ($array[0] == $testlocus) { print "\n" }
	if (length $array[2] == 0 || $array[2] eq 'NA' || $array[2] eq 'error') { print $out "$line\t$missing\n" }
	else {
		# mflag is set when a haplotype is found for which no ACGT code has been assigned.
		# This happends when e.g. a 'b' has been scored for a locus for which distinct b1 and b2 haplotypes have been defined.
		# It should be rare and occur in low coverage samples. Such instances are reported in the error file.
		$mflag = 0;
		%onecodes = ();
#		print "locus $array[0]\tscore $array[2]\t";
		@altstates = split ('\|',$array[2]);
		foreach $alt (@altstates) {
#			print "alt $alt\t";
			@ncodes = ();
			if ($array[0] == $testlocus) { print "altstate: $alt\t" }
			@haplos = split ('\/',$alt);
			foreach $h (@haplos) {
				if ($array[0] == $testlocus) { print "haplo $h\t" }
#				print "haplo $h\tACGT-code = $states{$array[0]}->{$h}\n";
				# valid haplotypes have an ACGT code assigned to them
				if (!(defined $states{$array[0]}->{$h})) { print "loc $array[0], haplo $h undefined\n" }
				if ($states{$array[0]}->{$h} eq '1') { $mflag = 1; last }
				else                                 { push @ncodes, $states{$array[0]}->{$h} }
			 	if ($array[0] == $testlocus) { print "nucleo $states{$array[0]}->{$h}\t" }
			}
			unless ($mflag) {
				if (@ncodes == 0)  { print "locus $array[0]\n" }
				@ncodes = sort { $a cmp $b } @ncodes;
				$code = $ncodes[0] . $ncodes[1];
				$onecodes{$code} = 1;
				if ($array[0] == $testlocus) { print "code $code\t" }
			}
		}
		if ($array[0] == $testlocus) { print "\n" }
		if ($array[0] == $testlocus) { foreach $key (keys %onecodes) { print "$key has onecode\t" }}
		if ($mflag) {  
			print $out "$line\t$missing\n";
			print $ef "$line - diplotype set to missing\n";
		}
		else {
			print $out "$line\t";
			for ($i=0;$i<@diplos-1;$i++) {
				if (exists $onecodes{$diplos[$i]}) { print $out "1 " }
				else                               { print $out "0 " }
			}
			if (exists $onecodes{$diplos[$i]}) { print $out "1\n" }
			else                               { print $out "0\n" }
		}
	}
}

close ($in);
close ($out);
close ($ef);	
			
		
