#!/usr/bin/perl

foreach (63..70) {
	print "rm -f fly$_*bam fly$_*bai fly$_*vcf.gz fly$_*cnv fly$_*bw fly$_*bed.gz fly$_*tsv\n";
}
