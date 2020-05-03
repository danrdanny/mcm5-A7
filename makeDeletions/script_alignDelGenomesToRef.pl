#!/usr/bin/perl

use strict;
use Getopt::Std;

my $refGenome = "/home/dem/projects/ref-genomes/dm6/dm6.fa";

my %opts;
getopts('t:d:r:n:hsz', \%opts);

if ($opts{'h'} || !$opts{'t'}) {
	print "
	Required:
		-t	threads
		-d	deletion size, 1-10, 10-1000, or both, which is converted to both

	Optional:
		-h	this helpful help
		-s	call SNPs with samtools when aligning
		-z	only call SNPs (no alignment)
	\n";
	exit 0;
}

my $runForReal = 1;
my $threads = $opts{'t'};

my $flystart = 0;
my $flyend   = 99;

my @deletionSize;
push(@deletionSize,$opts{'d'});
@deletionSize = qw / 1-10 10-1000/ if $opts{'d'} eq "both";

sub runCommand {
	print "$_[0]\n";
	`$_[0]` if $runForReal == 1;

	return 0;
}

my %chr;
   $chr{'chrX'} = 23542271;
   $chr{'chr2L'} = 23513712;
   $chr{'chr2R'} = 25286936;
   $chr{'chr3L'} = 28110227;
   $chr{'chr3R'} = 32079331;
   $chr{'chr4'} = 1348131;

foreach my $fly ($flystart..$flyend) {
	$fly = "fly0$fly" if $fly < 10;
	$fly = "fly$fly" if $fly >= 10;
		
	foreach my $deletionSize (@deletionSize) {
		my $finalAlnName = "$fly\_$deletionSize\_bpDel";
		my $fastq_f = "$fly\_$deletionSize\_bpDel.dm6.1.fq.gz";
		my $fastq_r = "$fly\_$deletionSize\_bpDel.dm6.2.fq.gz";

		my $samtoolsCmd = "bcftools mpileup -Ou -f $refGenome $finalAlnName.bam | bcftools call -Ou -mv | bcftools filter -Oz -s LowQual -e '%QUAL<20' > $finalAlnName.vcf.gz &";

		if ($opts{'z'}) {
			print "Skipping samtools for $finalAlnName because .vcf.gz file exists!\n" if -e "$finalAlnName.vcf.gz";
			next if -e "$finalAlnName.vcf.gz";
			runCommand($samtoolsCmd);
		} else {
			if (!-e "$finalAlnName.bam") {
				runCommand("bwa mem -t $threads $refGenome $fastq_f $fastq_r > $finalAlnName.sam");
				runCommand("samtools view -bS $finalAlnName.sam > $finalAlnName.bam");
				runCommand("rm -f $finalAlnName.sam");
				runCommand("samtools sort $finalAlnName.bam -o $finalAlnName.sorted.bam");
				runCommand("mv $finalAlnName.sorted.bam $finalAlnName.bam");
				runCommand("samtools index $finalAlnName.bam");
			}

			# Samblaster

			my $targetDir = `pwd`;
			chomp($targetDir);

			my $genomeFile = "/home/dem/projects/ref-genomes/dm6/dm6.chrom.sizes";
			
        		# Calculate genome coverage
        		if (!-e "$targetDir/$finalAlnName.bw") {
				#sleep 120;

                		print "[".localtime(time)."] $finalAlnName | RUN: Calculating genome coverage for $finalAlnName (forked)\n";
                		#my $pid = fork();
                		#die "unable to fork: $!" unless defined($pid);
                		#if (!$pid) {
                        		runCommand("genomeCoverageBed -bga -trackline -trackopts 'name=\"$finalAlnName\" visibility=2' -ibam $finalAlnName.bam > $finalAlnName.coverage");
                        		runCommand("head -1 $targetDir/$finalAlnName.coverage > $targetDir/$finalAlnName.coverage.sorted");
                        		runCommand("cat $targetDir/$finalAlnName.coverage | grep -v track > $targetDir/$finalAlnName.coverage.tmp");
                        		runCommand("sort -k1,1 -k2,2n $targetDir/$finalAlnName.coverage.tmp >> $targetDir/$finalAlnName.coverage.sorted");
                        		runCommand("mv $targetDir/$finalAlnName.coverage.sorted $targetDir/$finalAlnName.coverage");
	
                        		runCommand("bedGraphToBigWig $targetDir/$finalAlnName.coverage $genomeFile $targetDir/$finalAlnName.bw");
                        		runCommand("rm -f $targetDir/$finalAlnName.coverage $targetDir/$finalAlnName.coverage.tmp");
                        		runCommand("bigWigToBedGraph $targetDir/$finalAlnName.bw $targetDir/$finalAlnName.bed");

                        		my(%depthOfCov,%totalDepth,%totalDepthCount);
                        		print "[".localtime(time)."] $finalAlnName | Opening $targetDir/$finalAlnName.bed.\n";
                        		open INF,"$targetDir/$finalAlnName.bed" or die "Can't open $finalAlnName.bed: $!";
                        		while (<INF>) {
                                		my(@F) = split /\t/, $_;
                                		next unless $chr{$F[0]} > 0;
                                		my $chr = $F[0];

                                		foreach my $pos ($F[1]..$F[2]) {
                                        		$depthOfCov{$chr}{$pos} = $F[3];
                                        		$totalDepth{$chr}+= $F[3];
                                        		$totalDepthCount{$chr}++;
                                		}
                        		}
                        		close INF;
                        		print "[".localtime(time)."] $finalAlnName | Data collected, analyzing.\n";

                        		runCommand("gzip $targetDir/$finalAlnName.bed");

                        		my $output = "Chr\tPos\tlog2depth\n";
                        		foreach my $chr (sort keys %chr) {
                                		my $aveDepth = sprintf("%0.1f", $totalDepth{$chr} / $totalDepthCount{$chr});
	
                                		$output .= "\#$chr\t$aveDepth\n";

                                		my $max         = 10000;
                                		my $step        = 5000;
                                		my $min         = 0;

                                		until ($min >= $chr{$chr}) {
                                        		my($totalDepth,$count);
                                        		foreach my $id ($min..$max) {
                                                		next unless $depthOfCov{$chr}{$id} > 0; #id < $min || $id >= $max;
                                                		$totalDepth += $depthOfCov{$chr}{$id};
                                                		$count++;
                                        		}

                                        		my($currAveDepth,$ratio,$logDepth);
                                        		if ($count == 0) {
                                                		$logDepth = -10;
                                        		} else {
                                                		$currAveDepth = $totalDepth / $count;
                                                		$ratio = $currAveDepth / $aveDepth;
	
                                                		$logDepth = log($ratio)/log(2);
                                                		$logDepth = sprintf("%0.3f", $logDepth);
                                                		$currAveDepth = sprintf("%0.1f", $currAveDepth);
                                        		}

                                        		$output .= "$chr\t$min\t$logDepth\n";
                                
                                        		$min += $step;
                                        		$max += $step;
                                		}
                        		}
                                	
                        		open OUTF,">$targetDir/$finalAlnName.log2depth.tsv";
                        		print OUTF $output;
                        		close OUTF;
	
                        		#exit 1;
                		#}
        		} else {
                		print "[".localtime(time)."] $finalAlnName | OK:  Genome coverage file exists for $finalAlnName\n";
        		}


			if ($opts{'s'} && !-e "$finalAlnName.vcf.gz") {
                		my $pid = fork();
                		die "unable to fork: $!" unless defined($pid);
                		if (!$pid) {
					runCommand($samtoolsCmd);
                        		exit 1;
				}
			} else {
				print "Skipping samtools for $finalAlnName because .vcf.gz file exists!\n";
			}
		}
	}
}
