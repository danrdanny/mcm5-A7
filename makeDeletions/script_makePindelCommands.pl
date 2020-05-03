#!/usr/bin/perl

use strict;

foreach my $chr ("chrX","chr2L","chr2R","chr3L","chr3R") {
        my $output = "#!/bin/sh\n\n";

        foreach (32..99) {

		my $fly = "$_";
	   	$fly = "0$_" if $_ < 10;
	
		my $config = "fly$fly\_10-1000_bpDel.bam 200 fly$fly\_10-1000_bpDel";
		print "echo \"$config\" > pindel_config.$chr.txt\n";
		$output .= "echo \"$config\" > pindel_config.$chr.txt\n";
	
		my $cmd = "pindel -f /home/dem/projects/ref-genomes/dm6/dm6.fa -i pindel_config.$chr.txt -c $chr -o fly$fly\_10-1000_bpDel.pindel/fly$fly\_$chr -T 1 -w 40";
		print "$cmd\n";
		#`$cmd`;
		$output .= "$cmd\n";

	}

	open OUTF,">runPindel.mcm5.testReads.$chr.sh";
	print OUTF $output;
	close OUTF;
}
