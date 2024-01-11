#!/bin/bash

set -eu -o pipefail

conda activate py37

# Reading required data 
var=2000
var1=Drosophila_melanogaster.BDGP6.32.dna.toplevel.fa
var2=Drosophila_melanogaster.BDGP6.32.107.gtf

# Extract information for genes from input gff file
# awk '/gbkey=Gene/' $var2 > genes.gff
cat $var2 | grep -v "^#" |  awk '$3 == "gene"' > genes.gff


# Run R script to convert gff to bed format
Rscript gff2bed.r

#Generate index for fasta file, details are here: http://www.htslib.org/doc/samtools.html
samtools faidx $var1

#Create a table (from index file) that contains contig/chromosome sizes
cut -f1-2 ${var1}.fai > sizes.chr

#Create bed file that contains locations of promoters
bedtools flank -i genes.bed -g sizes.chr -l $var -r 0 -s > promoters.bed

#Extract promoter regions from fasta file with contigs/chromosomes using bed file
bedtools getfasta -s -fi $var1 -bed promoters.bed -fo promoters.fa -name

echo Results are in promoters.fa file

#Count promoters
result=$(grep -c '^>' promoters.fa)

echo Number of obtained promoters: $result


# use seqkit to extract promoter and genebody sequence
seqkit subseq --gtf Drosophila_melanogaster.BDGP6.32.107.gtf --feature gene Drosophila_melanogaster.BDGP6.32.dna.toplevel.fa --up-stream 2000  > genebody-and-up2kbp.fa

# Get 2000 bp upstream sequence 
seqkit subseq --gtf Drosophila_melanogaster.BDGP6.32.107.gtf --feature gene Drosophila_melanogaster.BDGP6.32.dna.toplevel.fa --up-stream 2000 --only-flank > up2kbp.fa

cat promoters.fa | seqkit locate -p AGGGCGG > match-AGGGCGG.txt
cat promoters.fa | seqkit locate -p AGGGTGG > match-AGGGTGG.txt
cat promoters.fa | seqkit locate -p AGGATAA > match-AGGATAA.txt
cat promoters.fa | seqkit locate -p AGGACAA > match-AGGACAA.txt

cat promoters.fa | seqkit locate -p AAGCT > match-AAGCT.txt


cat up2kbp.fa | seqkit locate -p AGGGCGG > match-AGGGCGG-seqkitSubseq.txt 
# ``seqkit subseq`` result has minor difference with ``bedtools getfasta``
# 结果略有不同！

cut -f 9 genes.gff | grep protein_coding > genes.txt
