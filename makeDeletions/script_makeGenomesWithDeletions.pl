#!/usr/bin/perl

use strict;
use Getopt::Std;

# datasets
my %chrs;
   $chrs{1} = "chrX";
   $chrs{2} = "chr2L";
   $chrs{3} = "chr2R";
   $chrs{4} = "chr3L";
   $chrs{5} = "chr3R";
   $chrs{6} = "chr4";

my %chromosomeSizes;
   $chromosomeSizes{"chrX"} = 23542271;
   $chromosomeSizes{"chr3R"} = 32079331;
   $chromosomeSizes{"chr3L"} = 28110227;
   $chromosomeSizes{"chr2L"} = 23513712;
   $chromosomeSizes{"chr2R"} = 25286936;
   $chromosomeSizes{"chr4"} = 1348131;

my %chromosomeFraction;
my $genomeSize = 23542271 + 32079331 + 28110227 + 23513712 + 25286936 + 1348131;
my($low,$newHigh) = (0,0);
foreach my $chr (keys %chromosomeSizes) {
	my $fraction = sprintf("%0.3f", $chromosomeSizes{$chr} / $genomeSize);
	$fraction *= 1000;
	$newHigh += $fraction;
	foreach ($low..$newHigh) {
		$chromosomeFraction{$_} = $chr;
	}
	#print "$chr\t$low - $newHigh\n";
	$low = $newHigh + 1;
}

my $SNPDensity = 500; # a SNP every 1 in this number;

# open ref genome and put it in a hash
my $chr;
my %cnsData;
my $refGenome = "/Users/danny/projects/genomes/dmel/dm6.fa";

open INF,"$refGenome" or die "Can't open file $refGenome: $!\n";
while (<INF>) {
	chomp($_);

	my $length = length($_);

	if ($_ =~ /^\>(\w+)$/ && $length < 60) {
		$chr = $1;
		#print "Found chr: $chr\n";
	} elsif ($_ =~ /^(\+)$/) {
		$chr .= "|+";
	} else {
		$cnsData{$chr} .= $_;
	}
}

# make a series of reference genomes with a varying number of deletions on one of 5 chromosomes

my($deletionNumMin,$deletionNumMax) = (5,20);
foreach my $deletionSizes ("1-10","10-1000") {
	my($deletionSizeMin,$deletionSizeMax) = $deletionSizes =~ /(\d+)\-(\d+)/;
	my $deletionTracker = "fly,delCount,chr,pos,delSize\n";

	foreach my $fly (0..99) {
		$fly = "0$fly" if $fly < 10;
		my $filename = "fly$fly\_$deletionSizeMin-$deletionSizeMax\_bpDel.dm6";

		my %posToDelete;
		my $deletions = int(rand($deletionNumMax)) + $deletionNumMin;
		foreach my $deletionNum (1..$deletions) {
			## which chromatid? We only recover 1
			my $chromatid = int(rand(4)) + 1;
			next unless $chromatid == 1;

			my $chrNum = int(rand(1000)) + 1;
			my $chrName = $chromosomeFraction{$chrNum};

			my $maxbp = $chromosomeSizes{$chrName};
			my $pos = int(rand($maxbp)) + 1;

			my $deletionSize = int(rand($deletionSizeMax)) + 1;

			my $deletionEnd = $pos + $deletionSize - 1;
			foreach my $id ($pos..$deletionEnd) {
				$posToDelete{$chrName}{$id} = 1;
			}
			#$deletionSize--; # subtract one on purpouse

			print "$fly,$chrName,$pos,$chrName:$pos-$deletionEnd,$deletionSize\n";
			$deletionTracker .= "$fly,$chrName,$pos,$chrName:$pos-$deletionEnd,$deletionSize\n";
		}

		foreach my $genome ("genomeA","genomeB") {
			my $output;
			foreach my $chr (keys %cnsData) {
				next if $chr =~ /\|\+/;
				# you make one genome with no X and one with no Y to make a heterozygous male
				# only make deletions in the genome with the X chromosome
				next if $chr =~ /chrY/ && $genome eq "genomeA"; # make a genome with no Y
				next if $chr =~ /chrX/ && $genome eq "genomeB"; # make a genome with no X

				$output .= ">$chr\n";

				my($pos,$count) = (0,0);
				foreach my $foo (split //, $cnsData{$chr}) {
					next unless $foo =~ /\w/;
					++$pos;
					next if $posToDelete{$chr}{$pos} && $genome eq "genomeA";
					++$count;

					# ?SNP
					my $addSNP = int(rand($SNPDensity)) + 1;
					if ($addSNP == 150) {
						my $newsnp = int(rand(4)) + 1;
						$foo = "A" if $newsnp == 1;
						$foo = "G" if $newsnp == 2;
						$foo = "C" if $newsnp == 3;
						$foo = "T" if $newsnp == 4;
					}

					$output .= $foo;
					$output .= "\n" if $count == 100;
					$count = 0 if $count == 100;
				}
				$output .= "\n";
			}
	
			open OUTF,">$filename.$genome.fa";
			print OUTF $output;
			close OUTF;
		}
	}

	open OUTF,">out_deletionTracker_$deletionSizeMin-$deletionSizeMax\_bpDeletion";
	print OUTF $deletionTracker;
	close OUTF;
}
