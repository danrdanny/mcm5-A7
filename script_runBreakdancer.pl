#!/usr/bin/perl

use strict;
use Getopt::Std;

## Command-line options + usage statement
my %opts;
getopts('s:t:abcdghx', \%opts); # Values in %opts

#if ($opts{'h'} || !$opts{'t'} || !$opts{'s'}) {
#	print "
#
#	Required:
#	\n";
#	exit 0;
#}

#SV detection with Breakdancer
my $breakdancer 	= "breakdancer-max";
my $breakdancercfg 	= "bam2cfg.pl";
my $minReadPairs 	= 10; # breakdancer: minimum number of read pairs needed to call a SV

#Data
my %chr;
   $chr{'chrX'} = 23542271;
   $chr{'chr2L'} = 23513712;
   $chr{'chr2R'} = 25286936;
   $chr{'chr3L'} = 28110227;
   $chr{'chr3R'} = 32079331;
   $chr{'chr4'} = 1348131;

## subroutines used
sub executeCommand {
	print "[".localtime()."] CMD: $_[0]\n";
	my $output = `$_[0]`;
	return($output);
}

my @stocks = qw/ Df-Mcm5.realigned mcm5-01.realigned mcm5-02.realigned mcm5-03.realigned mcm5-04.realigned mcm5-05.realigned mcm5-06.realigned mcm5-07.realigned mcm5-08.realigned mcm5-09.realigned mcm5-10.realigned mcm5-11.realigned mcm5-12.realigned mcm5-13.realigned mcm5-14.realigned mcm5-15.realigned mcm5-16.realigned mcm5-17.realigned mcm5-18.realigned mcm5-19.realigned mcm5-20.realigned mcm5-21.realigned mcm5-22.realigned mcm5-23.realigned mcm5-24.realigned mcm5-25.realigned mcm5-26.realigned mcm5-27.realigned mcm5-28.realigned Mcm5-A7.realigned /;

## Loop over the stocks we should be aligning. 
foreach my $newStockName (@stocks) {
	if (!-e "$newStockName.breakDancer.out") {
		print "[".localtime(time)."] $newStockName | RUN: Running breakdancer on $newStockName\n";
		executeCommand("$breakdancercfg -m $newStockName.bam > $newStockName.cfg");
	
		executeCommand("$breakdancer -r $minReadPairs $newStockName.cfg > $newStockName.breakDancer.out");
		executeCommand("$breakdancer -r $minReadPairs -t $newStockName.cfg > $newStockName.breakDancer.transChrom.out");

		executeCommand("rm -f $newStockName.cfg");
	} else {
		print "[".localtime(time)."] $newStockName | OK:  Breakdancer file $newStockName.ctx exists, skipping.\n";
	}
}
