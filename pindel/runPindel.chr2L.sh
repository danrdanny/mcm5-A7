#!/bin/sh

echo "Df-Mcm5.realigned.bam 200 Df-Mcm5" > pindel_config.chr2L.txt
pindel -f /home/dem/projects/ref-genomes/dm6/dm6.fa -i pindel_config.chr2L.txt -c chr2L -o pindel.out/Df-Mcm5_chr2L -T 1 -w 40
echo "Mcm5-A7.realigned.bam 200 Mcm5-A7" > pindel_config.chr2L.txt
pindel -f /home/dem/projects/ref-genomes/dm6/dm6.fa -i pindel_config.chr2L.txt -c chr2L -o pindel.out/Mcm5-A7_chr2L -T 1 -w 40
