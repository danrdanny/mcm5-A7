#!/usr/bin/perl

use strict;

#my $coverage_with_X = 10;
#my $coverage_without_X = 10;
#my $readLength = 100;
#my $meanDNAFragSize = 450;
#my $stddevDNAFragSize = 50;

my $art = "art_illumina";

# you make one genome with no X and one with no Y to make a heterozygous male
# only make deletions in the genome with the X chromosome
#next if $chr =~ /chrY/ && $genome eq "genomeA"; # make a genome with no Y
#next if $chr =~ /chrX/ && $genome eq "genomeB"; # make a genome with no X

my $size = $ARGV[0];
next unless $size eq "1-10" || $size eq "10-1000";

# filename: fly##_1-10_bpDel.dm6.genome#.fa
my @files = `ls -1 *$size\_bpDel.dm6.genomeA.fa`;

sub runCommand {
	my $cmd = $_[0];
	print "$cmd\n";
	`$cmd`;

	return 0;
}

foreach my $genomeA (@files) {
	chomp($genomeA);

	my $genomeB = $genomeA;
	$genomeB =~ s/genomeA/genomeB/;

	runCommand("$art -p -na -ss HS20 -i ./$genomeA -m 450 -s 50 -l 100 -f 11 -o $genomeA");
	runCommand("$art -p -na -ss HS20 -i ./$genomeB -m 450 -s 50 -l 100 -f 9 -o $genomeB");

	my($out) = $genomeA =~ /(\S+)\.genomeA\.fa/;
	my $outFile1 = "$out\.1.fq";
	my $outFile2 = "$out\.2.fq";

	runCommand("mv ".$genomeA."1.fq $outFile1");
	runCommand("mv ".$genomeA."2.fq $outFile2");

	runCommand("cat ".$genomeB."1.fq >> $outFile1");
	runCommand("cat ".$genomeB."2.fq >> $outFile2");

	runCommand("rm -f ".$genomeB."1.fq ".$genomeB."2.fq");

	runCommand("gzip $outFile1 $outFile2");
}
