#!/usr/bin/perl

use strict;
use Getopt::Std;

## This script looks for novel indels seen in the vcf file

## Command-line options
my %opts;
getopts('c:h', \%opts); # Values in %opts

my $numberToReport = 1; # marks unique sites less than or equal to this number

## Usage Statement
if ($opts{'h'} || !$opts{'c'}) {
	print "

	-c	Class: wt, c3g-hom, c3g-het, or corolla, mcm5

	-s	Look at samtools output (.vcf.gz)
	-g	Look at GATK HaplotypeCaller output 

	-h	This helpful help.

	\n";
	exit 0;
}

my $repeatMasker = "/home/dem/projects/ref-genomes/dm6/rmsk.dm6.tsv";
my %repeats;
print "[".localtime(time)."] Getting repeats from $repeatMasker.\n";
open INF,"$repeatMasker";
while (<INF>) {
	next if $_ =~ /Simple_repeat/;
	my(@F) = split /\t/, $_;
	foreach my $id ($F[6]..$F[7]) {
		$repeats{$F[5]}{$id} = 1;
		my $chr = $F[5];
		$repeats{$chr}{$id} = 1;
	}
}
close INF;

my %chrSizes;
   $chrSizes{'chrX'}  = 23542271;
   $chrSizes{'chr2L'} = 23513712;
   $chrSizes{'chr2R'} = 25286936;
   $chrSizes{'chr3L'} = 28110227;
   $chrSizes{'chr3R'} = 32079331;
   $chrSizes{'chr4'}  = 1348131;

my $baseDir = "/n/projects/c3g-corolla";
my $dir;
if ($opts{'c'} =~ /c3g-hom/) {
	$dir = "$baseDir/c3g-homozygotes";
} elsif ($opts{'c'} =~ /c3g-het/) {
	$dir = "$baseDir/c3g-heterozygotes";
} elsif ($opts{'c'} =~ /corolla/) {
	$dir = "$baseDir/corolla-homozygotes";
} elsif ($opts{'c'} =~ /wt/) {
	$dir = "/n/projects/co-ncogc-wt-2016";
} elsif ($opts{'c'} =~ /mcm5/) {
	$dir = "/n/projects/mcm5a7/vcfFiles";
} else {
	die "unknown class\n";
}

my @dirlist = `ls -1 $dir/`;
my @stocks;
foreach my $stock (@dirlist) {
	chomp($stock);

	my $skip = 0;
	if ($opts{'c'} =~ /c3g-hom/) {
		$skip = 1 unless $stock =~ /c3g-hom-\d+\.\d+/;
	} elsif ($opts{'c'} =~ /c3g-het/) {
		$skip = 1 unless $stock =~ /c3g-het-\d+\.\d+/;
	} elsif ($opts{'c'} =~ /corolla/) {
		$skip = 1 unless $stock =~ /corolla-hom-\d+\.\d+/;
	} elsif ($opts{'c'} =~ /wt/) {
		$skip = 1 unless $stock =~ /(cs\d+\.\d+|w\d+\.\d+)/;
	} elsif ($opts{'c'} =~ /mcm/) {
		$skip = 1 unless $stock =~ /mcm5-\S+samtools.vcf.gz/;
	}

	next if $skip == 1;

	push(@stocks,"$stock");
}

my(%vcfData,%snpCount,$stockCount);

foreach my $stock (@stocks) {
	#my $file = "$stock.realigned.HaplotypeCaller.vcf.gz";
	my $file;
	if ($opts{'c'} =~ /mcm/) {
		$file = $stock; # for mcm only
	} else {
		$file = "$stock.samtools.vcf.gz"; # okay for c3g-hom, c3g-het, corolla-hom, wt
	}

	my $filename = "$dir/$stock/$file";
	   $filename = "$dir/$file" if $opts{'c'} =~ /mcm/;
	if (-e $filename) {
		++$stockCount;
		print "[".localtime()."] opening $filename\n";
       		open(INF,sprintf("zcat %s |", $filename)) or die "Can't open $filename: $!";
       		#open(INF, $filename) or die "Can't open $filename: $!";
       		while (<INF>) {
			my(@F) = split /\t/, $_;
			#next unless $F[0] eq "chrX";
			next unless $_ =~ /INDEL/;
			#next unless $_ =~ /0\/1/;
			#next unless $F[4] =~ /^(A|G|C|T)$/;

			my $refLen = length($F[3]);
			my $altLen = length($F[4]);

			next unless $altLen < $refLen;

			my $min = $F[1] - 20;
			my $max = $F[1] + 20;

			foreach my $loc ($min..$max) {
				$snpCount{$F[0]}{$loc}++;
			}

			#my($indv) = $_ =~ /\;IDV\=(\d+)\;/;
			#next unless $indv > 5;
			#next unless $F[5] > 100;

			$vcfData{$stock}{$F[0]}{$F[1]} = "$F[3]|$F[4]|$F[5]";
		}
		close INF;
	} else {
		print "[".localtime(time)."] WARN: $filename missing\n";
	}
}

print "Stocks analyzed: $stockCount\n";
my $commonSites = $stockCount - 1;

my(%singletons);
my $singletonOutput = "Stock\tChr\tPos\n";
foreach my $stock (@stocks) {
	#next if $stock =~ /(wt|net-cn)/;
	#print "$stock\n";

	foreach my $chr (keys %chrSizes) { 
	#foreach my $chr ("chrX") { #,"chr2L","chr2R") { 
		my $singletonCount = 0;
		foreach my $id (sort {$a<=>$b} keys %{$vcfData{$stock}{$chr}}) {
			next unless $snpCount{$chr}{$id} == $numberToReport;
			#next if $snpCount{$chr}{$id} == $commonSites;
			my($ref,$alt,$score) = $vcfData{$stock}{$chr}{$id} =~ /(\S+)\|(\S+)\|(\S+)/;
	
			next if $repeats{$chr}{$id};
			next unless $score > 200;
	
			++$singletonCount;
			$singletons{$stock}{$chr}{$id} = 1;
			print "$stock\t$chr\t$id\t$chr:$id\t$score\t$alt\t$ref\n";
			$singletonOutput .= "$stock\t$chr\t$id\n";
		}
		#print "$singletonCount\t";
	}
	print "\n";
}

open OUTF,">./out_novelDeletions.samtools.$opts{'c'}.tsv";
print OUTF $singletonOutput;
close OUTF;
