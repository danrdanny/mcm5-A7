#!/usr/bin/perl

use strict;
use Getopt::Std;

my $minCNSBaseQuality	= 60;
my $minVCFScore		= 200; #only use SNPs with scores greater than or equal to this number
my $maxIndelLength	= 1; #ignore indels greater than or equal to this number
my $skipIndels		= 1; # 1 for yes, 0 for no
my $minParentalDepth 	= 20;
my $minChildDepth	= 8;

my $repeatMasker	= "rmsk.txt";
my $chromSizes		= "dm6.chrom.sizes";

my $parentA	= "Df-Mcm5";
my $parentB	= "Mcm5-A7";

my @chrList = qw/ chr3L chr3R /;

## Command-line options
my %opts;
getopts('ehls:c:', \%opts); # options as above. Values in %opts

## Usage Output
if ($opts{'h'}) {
        print "
	This script identifies CO and NCO events in child .vcf files from two known
	parental genotypes. 
	
	Required:

		None.

	Optional:

		-e If out_uniqueParentalVariants.tsv already exists you can use it 
		   and not re-make it every time.

		-l Populate or add to the cogcPositionsToSkip.tsv file.

		-s Specific stock you want to check.

		-c Chromosome: X, A(utosome), or B(oth). Default B(oth).

		-h This helpful help.
        \n";
        exit 0;
}

my %qualityScores;
foreach (0..90) {
        $_ += 33;
        my $ascii = chr($_);
        $_ -= 33;
        $qualityScores{$ascii} = $_;
}

## Gather system data
my $pwd = `pwd`;
chomp($pwd);

## Print out all variables for the user
print "[".localtime(time)."] \n";
print "[".localtime(time)."]       Script begin. Variables:\n";
printf "[".localtime(time)."] %30s %8d\n", "minimum cns base quality", $minCNSBaseQuality; 
printf "[".localtime(time)."] %30s %8d\n", "minimum VCF score", $minVCFScore; 
printf "[".localtime(time)."] %30s %8d\n", "min depth req for parents", $minParentalDepth; 
printf "[".localtime(time)."] %30s %8d\n", "min depth req for child", $minChildDepth; 
if ($skipIndels == 0) {
	printf "[".localtime(time)."] %30s %8d\n", "maximum INDEL length", $maxIndelLength; 
} else {
	printf "[".localtime(time)."] %30s\n", "Skipping all indels"; 
}
print "[".localtime(time)."] \n";

## Grab chromosome names and sizes from the chrom.sizes file
my %chromosomeSizes;
   #$chromosomeSizes{'chrX'} = 23542271;
   #$chromosomeSizes{'chr2L'} = 23513712;
   #$chromosomeSizes{'chr2R'} = 25286936;
   $chromosomeSizes{'chr3L'} = 28110227;
   $chromosomeSizes{'chr3R'} = 32079331;
   #$chromosomeSizes{'chr4'} = 1348131;

if (!$opts{'e'}) {  # -e flag is to use existing out_uniqueParentalVariants.tsv file
	## Open repeatmasker file
	my %repeats;
	print "[".localtime(time)."] Getting repeats from $repeatMasker.\n";
	open INF,"~/projects/genomes/dmel/dm6.rmsk.txt";
	while (<INF>) {
		next if $_ =~ /Simple_repeat/;
		my(@F) = split /\t/, $_;
		foreach my $id ($F[6]..$F[7]) {
			$repeats{$F[5]}{$id} = 1;
			my $chr = $F[5];
			#$chr =~ s/chr//;
			$repeats{$chr}{$id} = 1;
		}
	}
	close INF;
	print "[".localtime(time)."] Repeats gathered.\n";

	## Grab parental VCF data
	my(%parentalSNPs,%variantPositions);
        print "[".localtime(time)."] Getting parental data.\n";
        foreach my $parent ($parentA, $parentB) {
                my($SNPcount,$skipRepeat,$skipScore,$skipHet,$skipAlt,$skipINDEL,$skipDepth,$indelCount) = (0,0,0,0,0,0,0,0);
		#my $parentVCFFile = "$parent.realigned.HaplotypeCaller.vcf.gz";
		my $parentVCFFile = "$parent.chr3.samtools.vcf.gz";
                print "[".localtime(time)."] - Getting vcf data from $parentVCFFile.\n";
                open INF,"gunzip -c vcfFiles/$parentVCFFile |" or die "Can't open $parentVCFFile: $!";
                while (<INF>) {
                        my(@F) = split /\t/, $_;
                        my($chr,$id,$vcfScore,$ref,$alt) = ($F[0],$F[1],$F[5],$F[3],$F[4]);
			#$vcfScore =~ s/\.\d+$//g;

#next unless $id > 1990000 && $id < 1997000;

			my $lengthRef = length($ref);
			my $lengthAlt = length($alt);
			next if $lengthRef > 1;
			next if $lengthAlt > 1;

			next if $_ =~ /INDEL/;
                        next unless $chr =~ /^(chr3L|chr3R)$/; # Drosophila specific
#next unless $chr =~ /^(chr3L)$/; # Drosophila specific
                        $skipRepeat++ if $repeats{$chr}{$id} == 1;      # skip if in repetative region
                        next if $repeats{$chr}{$id} == 1;
                        $skipScore++ if $vcfScore < $minVCFScore;       # skip based on vcf score
                        next unless $vcfScore >= $minVCFScore;
                        #$skipHet++ if $F[9] =~ /0\/1/;                  # skip it if it's het when looking at parents
                        #next if $F[9] =~ /0\/1/;
                        $skipAlt++ if $alt =~ /\,/;                     # skip if the alt allele has two bases
                        next if $alt =~ /\,/;

#print "$parent\t$chr,$id,$vcfScore,$ref,$alt\n";
                        # calculate depth
                        #my($totalDepth) = $_ =~ /\;DP\=(\d+)\;Excess/;
                        my($totalDepth) = $_ =~ /DP\=(\d+)\;VDB\=/;
                        $skipDepth++ if $totalDepth <= $minParentalDepth;
                        next if $totalDepth <= $minParentalDepth;

			$parentalSNPs{$parent}{$chr}{$id} = "$ref|$alt|$vcfScore|$totalDepth";
			$variantPositions{$chr}{$id} = 1;
			++$SNPcount unless $_ =~ /INDEL/;
		}
		close INF;
	
		printf "[".localtime(time)."] %30s %8d\n", "Skip d/t repeat", $skipRepeat;
                printf "[".localtime(time)."] %30s %8d\n", "Skip d/t score < $minVCFScore", $skipScore;
                printf "[".localtime(time)."] %30s %8d\n", "Skip d/t het SNP", $skipHet;
                printf "[".localtime(time)."] %30s %8d\n", "Skip d/t two alt alleles", $skipAlt;
                #printf "[".localtime(time)."] %30s %8d\n", "Skip d/t INDEL >= $maxIndelLength", $skipINDEL;
                printf "[".localtime(time)."] %30s %8d\n", "Skip d/t depth <= $minParentalDepth", $skipDepth;
                #printf "[".localtime(time)."] %30s %8d\n", "Total INDELs", $indelCount;
                printf "[".localtime(time)."] %30s %8d\n", "Total SNPs", $SNPcount;
                #printf "[".localtime(time)."] %30s %8d\n", "Total variants", $total;
                print "[".localtime(time)."] \n";
	}

	print "[".localtime(time)."] Done gathering parental data.\n";
	print "[".localtime(time)."] \n";

	## Grab parental consensus sequence data
	my %cns;
	print "[".localtime(time)."] Opening consensus file for each parental line - can be slow.\n";
	foreach my $parent ($parentA, $parentB) {
		print "[".localtime(time)."] - Opening consensus files for $parent.\n";
		foreach my $chr (@chrList) {
			my $cns = "$parent.realigned.cns.$chr.gz";
			print "[".localtime(time)."]    $cns\n";
        		#open INF,"$cns" or die "Can't open file $cns: $!\n";
			open INF,"gunzip -c $cns |" or die "Can't open $cns: $!";
			while (<INF>) {
				chomp($_);
				my(@F) = split /\t/, $_;
			
				$F[1] =~ tr/[a-z]/[A-Z]/;
                		$cns{$parent}{$chr}{$F[0]} = "$F[1]|$F[2]";
                	}
			close INF;
       		}
	}
	print "[".localtime(time)."] Done collecting consensus sequence.\n";
	print "[".localtime(time)."] \n";

	## Identify positions that differ between $parentA and $parentB
	my %diffVariants;
	my %diffVariantCounts;
	my $output = "Chr\tID\tRef\t$parentA\_cns\tScore\t$parentA\_vcf\tScore\tDepth\t$parentB\_cns\tScore\t$parentB\_vcf\tScore\tDepth\n";
	foreach my $chr (keys %variantPositions) {
		my %skip;
		foreach my $id (sort {$a<=>$b} keys %{$variantPositions{$chr}}) {
			my($parentAcns,$parentAcnsScore) = $cns{$parentA}{$chr}{$id} =~ /(\w+)\|(\d+)/;
			my($parentBcns,$parentBcnsScore) = $cns{$parentB}{$chr}{$id} =~ /(\w+)\|(\d+)/;

#print "\nID: $id\ncnsA: $cns{$parentA}{$chr}{$id}\ncnsB: $cns{$parentB}{$chr}{$id}\nparentA: $parentalSNPs{$parentA}{$chr}{$id}\nparentB: $parentalSNPs{$parentB}{$chr}{$id}\n";

			$skip{'missingCNS'}++ if !$parentAcns || !$parentBcns;
			next if !$parentAcns || !$parentBcns;

			$skip{'LowCNSScore'}++ if $parentAcnsScore < $minCNSBaseQuality || $parentBcnsScore < $minCNSBaseQuality;
			next if $parentAcnsScore < $minCNSBaseQuality || $parentBcnsScore < $minCNSBaseQuality;

			my($refA,$parentAalt,$AvcfScore,$AvcfDepth) = $parentalSNPs{$parentA}{$chr}{$id} =~ /(\w+)\|(\w+)\|(\d+)\|(\d+)/;
			my($refB,$parentBalt,$BvcfScore,$BvcfDepth) = $parentalSNPs{$parentB}{$chr}{$id} =~ /(\w+)\|(\w+)\|(\d+)\|(\d+)/;

			my $ref = $refA;
		   	$ref = $refB if !$ref;

			$parentAalt 	= "."  if !$parentAalt; 
			$AvcfScore 	= "..." if $parentAalt eq '.'; 
			$AvcfDepth	= ".." if $parentAalt eq '.'; 
			$parentBalt 	= "."  if !$parentBalt; 
			$BvcfScore 	= "..." if $parentBalt eq '.';
			$BvcfDepth 	= ".." if $parentBalt eq '.';

#print "$parentAcns\t$parentAcnsScore\t$parentAalt\t$AvcfScore\t$AvcfDepth\t$parentBcns\t$parentBcnsScore\t$parentBalt\t$BvcfScore\t$BvcfDepth\n";

			$skip{'IdenticalAlt'}++ if $parentAalt eq $parentBalt;
			next if $parentAalt eq $parentBalt;

			$diffVariants{$chr}{$id} = "$parentAcns\t$parentAcnsScore\t$parentAalt\t$AvcfScore\t$AvcfDepth\t$parentBcns\t$parentBcnsScore\t$parentBalt\t$BvcfScore\t$BvcfDepth";
			$output .= "$chr\t$id\t$ref\t$diffVariants{$chr}{$id}\n";
			$diffVariantCounts{$chr}++;
	
			# need to add depth of coverage filter at the cns positions
		}

		print "[".localtime(time)."] \n";
		print "[".localtime(time)."] $chr skip missing CNS: $skip{'missingCNS'}\n";
		print "[".localtime(time)."] $chr skip low CNS score: $skip{'lowCNSScore'}\n";
		print "[".localtime(time)."] $chr skip identical alt alleles: $skip{'IdenticalAlt'}\n";
	}

	# deleting structures not needed. See this thread: http://www.perlmonks.org/?node_id=182343
	undef %parentalSNPs;
	undef %cns;

	open OUTF,">$pwd/out_uniqueParentalVariants.tsv";
	print OUTF $output;
	close OUTF;

	undef $output;

	printf "[".localtime(time)."] \n";
	printf "[".localtime(time)."] %39s\n", "Per chromosome counts";
	my $totalCount;
	foreach my $chr (keys %diffVariantCounts) {
		#my $tmpChr = "chr".$chr;
		my $snpPerBP = sprintf("%0.0f", $chromosomeSizes{$chr} / $diffVariantCounts{$chr});
		printf "[".localtime(time)."] %15s %8d - 1 snp / %4d bp\n", $chr, $diffVariantCounts{$chr}, $snpPerBP; 
		$totalCount += $diffVariantCounts{$chr};
	}
	printf "[".localtime(time)."] %30s %8d\n", "Total Count", $totalCount;
	print "[".localtime(time)."] \n";
	print "[".localtime(time)."] Unique variants saved in out_uniqueParentalVariants.tsv\n";
	print "[".localtime(time)."] \n";
} else {
	die "Error: out_uniqueParentalVariants.csv doesn't exist! Quitting.\n" if !-e "$pwd/out_uniqueParentalVariants.tsv";
	print "[".localtime(time)."] -e flag passed, out_uniqueParentalVariants.tsv exists.\n";
}

#exit 0;

# open out_uniqueParentalVariants.csv
my(%parental,%SNP,$countVariants,%countSNPs);
print "[".localtime(time)."] Opening out_uniqueParentalVariants.tsv.\n";
open INF,"out_uniqueParentalVariants.tsv" or die "Can't open out_uniqueParentalVariants.tsv: $!";
while (<INF>) {
	my(@F) = split /\t/, $_;
	next unless $F[1] =~ /[0-9]/;

	my($chr,$id,$ref,$wcns,$wvcf,$cscns,$csvcf) = ($F[0],$F[1],$F[2],$F[3],$F[5],$F[8],$F[10]);

	#next if $wvcf ne $wcns && $wvcf !~ /\./;
	#next if $csvcf ne $cscns && $csvcf !~ /\./;
	#next if $wcns !~ /[A|G|C|T]/ || $cscns !~ /[A|G|C|T]/;
	#next if $wcns eq $cscns;

	$SNP{$chr}{$id} = $ref;
	$parental{$parentA}{$chr}{$id} = $wcns;
	$parental{$parentB}{$chr}{$id} = $cscns;
	++$countVariants;
	$countSNPs{$chr}++;
}
print "[".localtime(time)."] $countVariants total variants identified.\n";
printf "[".localtime(time)."] \n";
printf "[".localtime(time)."] %39s\n", "Per chromosome counts";
foreach my $chr (@chrList) {
	my $snpPerBP = sprintf("%0.0f", $chromosomeSizes{$chr} / $countSNPs{$chr});
	printf "[".localtime(time)."] %17s %10d - 1 snp / %4d bp\n", $chr, $countSNPs{$chr}, $snpPerBP; 
}
print "[".localtime(time)."] \n";


######
open INF,"./cogcPositionsToSkip.tsv";
while (<INF>) {
	my(@F) = split /\t/, $_;
	my($chr,$count) = ($F[0],$F[2]);
	my($id1,$id2) = $F[1] =~ /(\d+)\-(\d+)/;

	next if $count == 1;
	print "Skipping: $chr\t$id1-$id2\n";
	
	foreach ($id1..$id2) {
		#print "Skipping: $chr\t$_\n";
		$parental{$parentA}{$chr}{$_} = undef;
		$parental{$parentB}{$chr}{$_} = undef;
	}
}
close INF;

my($stockCount,%eventCount,%eventStocks);

my $sep= "+-----+------++-----------+-----+----++-----------+-----+----+-----------+-----++-----------+-----+----++----------+----------++-------------+-------------+";
#print "$sep\n";
#print "| Stk |  Chr ||   Last ID | Bas | Pa ||        ID | Bas | Pa |        ID | Bas ||   Next ID | Bas | Pa ||     fGap |     bGap || Lst->CurSNP | SNPsInEvent | \n";
#print "Chr\tStock\tLastParent\tRange\t\n";

my @files = `ls -1 vcfFiles/ | grep mcm5 | grep samtools`;
foreach my $file (@files) {
	chomp($file);
	my($stock) = $file =~ /^(mcm5-\d\d)\.chr3/;
	print "[".localtime(time)."] $stock\n";
	$stockCount++;

	my(%vcfData,%allVCFData);

	# open VCF file, store data in %vcf
	#my $vcfFile = "vcfFiles/$stock.realigned.HaplotypeCaller.vcf.gz";
	my $vcfFile = "vcfFiles/$stock.chr3.samtools.vcf.gz";
	#next if !-e "$vcfFile";
	open INF,"gunzip -c $vcfFile |" or die "Can't open $vcfFile: $!";
	while (<INF>) {
		my(@F) = split /\t/, $_;
		my($chr,$id,$vcfScore,$ref,$alt) = ($F[0],$F[1],$F[5],$F[3],$F[4]);
		$allVCFData{$chr}{$id} = 1;
		#$vcfScore =~ s/\.\d+$//g;

		next if $_ =~ /INDEL/;
		next unless $chr =~ /^(chr3L|chr3R)$/; # Drosophila specific
		next unless $parental{$parentA}{$chr}{$id} && $parental{$parentB}{$chr}{$id};
		next unless $vcfScore >= $minVCFScore;	
               	next if $alt =~ /\,/; 		

		my $lenRef = length($ref);
		my $lenAlt = length($alt);
		next if $lenRef > 1;
		next if $lenAlt > 1;

		my $type = "het";
		   $type = "hom" if $F[9] =~ /1\/1/;

               	next if $type eq "het" && $chr eq "chrX"; 	# skip X chr calls that are het

		my $parent = "neither";
		
		if ($type eq "hom" ) {
			$parent = "$parentA" if $alt eq $parental{$parentA}{$chr}{$id};
			$parent = "$parentB" if $alt eq $parental{$parentB}{$chr}{$id};
		} elsif ($type eq "het") {
			$parent = "both" if $alt eq $parental{$parentA}{$chr}{$id} && $ref eq $parental{$parentB}{$chr}{$id};
			$parent = "both" if $ref eq $parental{$parentA}{$chr}{$id} && $alt eq $parental{$parentB}{$chr}{$id};
		}

		$vcfData{$chr}{$id} = "$parent|$type|$ref|$alt";
		#print "$vcfData{$chr}{$id}\n" if $id < 1000000 && $chr eq "chr2L";
	}

	foreach my $chr (@chrList) {
		#next if $chr eq "chrX";
		foreach my $id (sort {$a<=>$b} keys %{$SNP{$chr}}) {
			next if $vcfData{$chr}{$id};
			next if $allVCFData{$chr}{$id} == 1;

			my $parent = "neither";
			$parent = "$parentA"  if $SNP{$chr}{$id} eq $parental{$parentA}{$chr}{$id};
			$parent = "$parentB" if $SNP{$chr}{$id} eq $parental{$parentB}{$chr}{$id};

			$vcfData{$chr}{$id} = "$parent|ref|$SNP{$chr}{$id}|N";
		}

		#print "[".localtime(time)."] $stock | $chr \n";
		my($lastParent,$origLastParentSNPChange,$lastParentSNPChange,$lastSNP,$snpCount);
		foreach my $id (sort {$a<=>$b} keys %{$vcfData{$chr}}) {
			my($parent,$type,$ref,$alt) = split /\|/, $vcfData{$chr}{$id};
			#print "$vcfData{$chr}{$id}\n" ; #if $id == 7023293;

			next if $parent eq "neither";

			if (!$lastParent) {
				$lastParent = $parent;
				$origLastParentSNPChange = $id;
				$lastParentSNPChange = $id;
				$snpCount = 1;
				$lastSNP = $id;
			} elsif ($parent eq $lastParent) {
				$snpCount++;
				$lastSNP = $id;
			} elsif ($parent ne $lastParent) {
				my($depth,$allele,$alleles) = getAllele($stock,$chr,$id);

				my $tmpParent = "neither";
				   $tmpParent = "$parentA" if $allele eq $parental{$parentA}{$chr}{$id} && $alleles !~ /\|\w\|/;
				   $tmpParent = "$parentB" if $allele eq $parental{$parentB}{$chr}{$id} && $alleles !~ /\|\w\|/;

				if ($alleles =~ /(\w)\|(\w)\|/) {
					my($a,$b) = ($1,$2);
					$tmpParent = "both" if $a eq $parental{$parentA}{$chr}{$id} && $b eq $parental{$parentB}{$chr}{$id};
					$tmpParent = "both" if $b eq $parental{$parentA}{$chr}{$id} && $a eq $parental{$parentB}{$chr}{$id};
				}

				#print "$stock\t$chr\t$id\t$depth\t$allele\t$tmpParent\t$parent\n" if $id == 24865481;
				if ($depth > $minChildDepth && $tmpParent eq $parent) {
					# Three gaps: before event, event, after event
					my $beforeEventGap = $lastParentSNPChange - $origLastParentSNPChange;
					my $eventGap =       $lastSNP - $lastParentSNPChange;
					my $afterEventGap =  $id - $lastSNP;

					my $printID1 = $origLastParentSNPChange - 10;
					my $printID2 = $id + 10;

					print "$chr,$stock,$lastParent,$parent,$chr:$printID1-$printID2,$chr:$lastParentSNPChange-$lastSNP,->,$snpCount,$beforeEventGap,$eventGap,$afterEventGap\n";

					$eventCount{$chr}{"$lastParentSNPChange-$lastSNP"}++ if $snpCount < 100;
					$eventStocks{$chr}{"$lastParentSNPChange-$lastSNP"} .= "$stock," if $snpCount < 100;
				
					# update variables
					$lastParent = $parent;
					$origLastParentSNPChange = $lastSNP;
					$lastParentSNPChange = $id;
					$lastSNP = $id;
					$snpCount = 1;
				}
			}
		}

		my $beforeEventGap = $lastParentSNPChange - $origLastParentSNPChange;
		my $eventGap = $lastSNP - $lastParentSNPChange;
		print "$chr,$stock,$lastParent,END_OF_CHR,$chr:$origLastParentSNPChange-$lastSNP,$chr:$lastParentSNPChange-$lastSNP,->,$snpCount,$beforeEventGap,$eventGap,-\n";
	}
}

if ($opts{'l'}) {
	open OUTF,">cogcPositionsToSkip.tsv";
	print OUTF "Chr\tRange\tCount\tStocks Seen In (stocks checked: $stockCount)\n";
	foreach my $chr (keys %eventCount) {
		foreach my $id (sort {$a<=>$b} keys %{$eventCount{$chr}}) {
			my $count = $eventCount{$chr}{$id};
			#next unless $count > 3;
			my($id1,$id2) = $id =~ /(\d+)\-(\d+)/;
			my $gap = $id2 - $id1;
			print OUTF "$chr\t$id\t$count\t$gap\t$eventStocks{$chr}{$id}\n";
		}
	}
	close OUTF;
}

#####################

sub getAllele {
	my($stock,$chr,$id) = ($_[0],$_[1],$_[2]);
	$stock =~ s/^(mcm5-\d+)\.\S+/$1/;

	my %alleleCount;

	foreach my $chunk (split /\n/, `samtools view $stock.realigned.bam $chr:$id-$id 2>/dev/null`) {
		my(@F) = split /\t/, $chunk;
		next if $F[4] < 60;
		my($splitInfo,$originalSequence,$originalQuality) = ($F[5],$F[9],$F[10]);
	
       		my @seq;
       		foreach (split //, $originalSequence) {
       			push(@seq,$_);
       		}
       		$originalSequence = undef;
       		my $count = 0;
       		foreach (split //, $originalQuality) {
       			my $score = $qualityScores{$_};
       			my $base = $seq[$count];
       			$base = "n" if $qualityScores{$_} < 15;
       			$originalSequence .= $base;
       			++$count;
       		}

        	my $unknownAction = 0;
		my $seq;
        	foreach (split /([0-9]+[A-Z])/, $splitInfo) {
        		my($len,$action) = $_ =~ /([0-9]+)([A-Z])/;
        		next unless $len =~ /[0-9]/;

        		if ($action =~ /M/) {
        			$originalSequence =~ /^(.{$len})/;
        			$seq .= $1;
        			$originalSequence =~ s/^.{$len}(.+)/$1/;
        		} elsif ($action =~ /D/) {
        			foreach (1..$len) {
        				$seq .= ".";
        			}
        		} elsif ($action =~ /I/) {
        			$originalSequence =~ s/^.{$len}(.+)/$1/;
        		} elsif ($action =~ /S/) {
        			$originalSequence =~ s/^.{$len}(.+)/$1/;
        		} else {
				$unknownAction = 1;
			}
		}

		next if !$seq;
		my $length = length($seq);
		my $startID = $F[3];
		my $endID = $F[3] + $length - 1;

		$startID = $endID if $startID > $id;

		my $currID = $startID;
		foreach (split //, $seq) {
			#print "$currID\t$_\n" if $currID == 7023293;
			$alleleCount{$_}++ if $currID == $id;
			$currID++;
		}

	}

	my $readCount = $alleleCount{'A'} + $alleleCount{'G'} + $alleleCount{'C'} + $alleleCount{'T'};
	my($allele,$alleles);
	if ($readCount > 0) {
		my $tmpReadCount = $readCount / 2;
		$allele = "A" if $alleleCount{'A'} >= $tmpReadCount;
		$allele = "G" if $alleleCount{'G'} >= $tmpReadCount;
		$allele = "C" if $alleleCount{'C'} >= $tmpReadCount;
		$allele = "T" if $alleleCount{'T'} >= $tmpReadCount;

		my $alleleFreq = sprintf("%0.0f", (($alleleCount{$allele} / $readCount ) * 100));

		if ($alleleFreq >= 30 && $alleleFreq <= 70) {
			$alleles .= "A|" if $alleleCount{'A'} >= 2;
			$alleles .= "G|" if $alleleCount{'G'} >= 2;
			$alleles .= "C|" if $alleleCount{'C'} >= 2;
			$alleles .= "T|" if $alleleCount{'T'} >= 2;
		} elsif ($alleleFreq >= 5 && $alleleFreq <= 95) {
			$allele = "N";
		}

		#print "$alleleFreq\t$allele\t$alleles\t$alleleCount{'A'}, $alleleCount{'G'}, $alleleCount{'C'}, $alleleCount{'T'}\n" if $id == 24865481;
	}

	return($readCount,$allele,$alleles);
}
