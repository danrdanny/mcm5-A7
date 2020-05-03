#!/usr/bin/perl

use strict;

foreach my $chr ("chrX","chr2L","chr2R","chr3L","chr3R") {
	my $output = "#!/bin/sh\n\n";

	foreach (7..28) {
		my $fly = "$_";
	   	$fly = "0$_" if $_ < 10;

		my $config = "mcm5-$fly\.realigned.bam 200 mcm5-$fly";
		print "echo \"$config\" > pindel_config.$chr.txt\n";
		$output .= "echo \"$config\" > pindel_config.$chr.txt\n";
	
		#open OUTF,">pindel_config.txt";
		#print OUTF $config;
		#close OUTF;

		my $cmd = "pindel -f /home/dem/projects/ref-genomes/dm6/dm6.fa -i pindel_config.$chr.txt -c $chr -o pindel.out/mcm5-$fly\_$chr -T 2 -w 40"; 
		print "$cmd\n";
		$output .= "$cmd\n";
		#`$cmd`;
	}

	open OUTF,">runPindel.$chr.sh";
	print OUTF $output;
	close OUTF;
}
