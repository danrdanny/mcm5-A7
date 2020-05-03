#!/bin/sh

echo "fly08_10-1000_bpDel.bam 200 fly08_10-1000_bpDel" > pindel_config.chrX.txt
pindel -f /home/dem/projects/ref-genomes/dm6/dm6.fa -i pindel_config.chrX.txt -c chrX -o fly08_10-1000_bpDel.pindel/fly08_chrX -T 5 -w 40
