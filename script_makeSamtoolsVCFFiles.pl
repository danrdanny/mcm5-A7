#!/usr/bin/perl

use strict;

sub doCommand {
	print "$_[0]\n";
	`$_[0]`;

	return 1;
}

my @files = `ls -1 mcm5-*.realigned.bam`;
foreach my $file (@files) {
	chomp($file);

	doCommand("samtools view -b $file chr3L > tmp.chr3L.bam");
	doCommand("samtools view -b $file chr3R > tmp.chr3R.bam");

	my $mergeBam = $file;
	$mergeBam =~ s/\.bam//;
	doCommand("samtools merge $mergeBam.chr3.bam tmp.chr3L.bam tmp.chr3R.bam");
	doCommand("rm -f tmp.chr3L.bam tmp.chr3R.bam");

#bcftools mpileup -Ou -f ~/projects/genomes/dmel/dm6.fa mcm5-01.realigned.chr3L.bam | bcftools call -Ou -mv | bcftools filter -s LowQual -e '%QUAL<20 || DP>100' > mcm5-01.chr3L.3-4.vcf

}
