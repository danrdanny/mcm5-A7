#!/usr/bin/perl

use strict;

my $dirLoc = "/n/projects/mcm5a7";
my $pwd = `pwd`;
chomp($pwd);

my %chr;
$chr{'chrX'} = 23542271;
$chr{'chr2L'} = 23513712;
$chr{'chr2R'} = 25286936;
$chr{'chr3L'} = 28110227;
$chr{'chr3R'} = 32079331;
#$chr{'chr4'} = 1348131;

my @stocks;
foreach my $fly (1..5) {
	$fly = "mcm5-0$fly" if $fly < 10;
	$fly = "mcm5-$fly" if $fly >= 10;

	push(@stocks,"$fly");
}

#my %repeats;
#print "[".localtime(time)."] Getting repeats from $repeatMasker.\n";
#open INF,"~/projects/genomes/dmel/dm6.rmsk.txt";
#while (<INF>) {
	#next if $_ =~ /Simple_repeat/;
	#my(@F) = split /\t/, $_;
	#foreach my $id ($F[6]..$F[7]) {
		#$repeats{$F[5]}{$id} = 1;
		#my $chr = $F[5];
		##$chr =~ s/chr//;
		#$repeats{$chr}{$id} = 1;
	#}
#}
#close INF;
#print "[".localtime(time)."] Repeats gathered.\n";


foreach my $chr (keys %chr) {
	next unless $chr eq "chr2L";

	my(%candidateDeletions,%deletionCount,%deletionDetail);
	foreach my $stock (@stocks) {
		my $pindelData;
		#my($pindelName) = $stock =~ /^(fly\d\d)\_/;
		my $pindelFile = "$stock\_$chr\_D";
		my $fullPath = "$dirLoc/pindel.out/$pindelFile";
		print "opening $fullPath\n";
		open INF,"$fullPath" or die "can't open file: $fullPath $!\n"; 
		while (<INF>) {
			chomp($_);
			$pindelData .= "$_\n";
		}
		close INF;

		foreach my $chunk (split /\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#/, $pindelData) {
			$chunk =~ s/\s+/ /g;
			my($supports) = $chunk =~ /Supports\s(\d+)/;
			my($delStart,$delEnd) = $chunk =~ /BP\s(\d+)\s(\d+)\sBP_range/;
			my $gap = $delEnd - $delStart;

			#next unless $delStart >1000000 && $delEnd < 2000000;

			next unless $supports > 5;
			next unless $gap > 10;

			my $delStartWindow = $delStart - 50;
			my $delEndWindow = $delEnd + 50;

			foreach my $pos ($delStartWindow..$delEndWindow) {
				$deletionCount{$chr}{$pos}++;
			}

			foreach my $pos ($delStart..$delEnd) {
				$deletionDetail{$chr}{$pos} .= "$stock\,";
			}
			#push(@candidateDeletions,"$stock,$chr,$supports,$delStart,$delEnd");
			$candidateDeletions{$chr}{$delStart} = "$stock,$supports,$delEnd";
			#print "$stock\t$chr\t$delStart\t$delEnd\n";
		}
	}

	my %printCandidateDeletions;
	foreach my $delStart (keys %{$candidateDeletions{$chr}}) { 
		my($stock,$supports,$delEnd) = split /\,/, $candidateDeletions{$chr}{$delStart};

		my $skip = 0;
		foreach my $pos ($delStart..$delEnd) {
			$skip = 1 if $deletionCount{$chr}{$pos} > 1;
		}
		next if $skip == 1;

		my $gap = $delEnd - $delStart;
		next unless $gap > 10;
		next unless $supports > 10;
		## we consider the deletion real if we don't see it overlapping any other deletions
		$printCandidateDeletions{$stock}{$delStart} = "$stock\t$chr\t$supports\t$delStart\t$delEnd\t$gap\t$chr:$delStart..$delEnd";
		#print "$stock\t$chr\t$supports\t$delStart\t$delEnd\t$gap\t$chr:$delStart..$delEnd\n";
	}

	foreach my $stock (sort keys %printCandidateDeletions) {
		foreach my $delStart (sort keys %{$printCandidateDeletions{$stock}}) {
			print "$printCandidateDeletions{$stock}{$delStart}\n";
		}
	}
}
