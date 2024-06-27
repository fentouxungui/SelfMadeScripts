ln -s ~/reference/Genome/Drosophlia/Flybase/r6.51/dmel-all-r6.51.gtf
ln -s ~/reference/Genome/Drosophlia/Flybase/r6.51/dmel-all-chromosome-r6.51.fasta
ln -s ~/reference/Genome/Drosophlia/Flybase/r6.51/dmel-all-chromosome-r6.51.fasta.fai


# extract all transcript from gtf file for promoter analysis
cat dmel-all-r6.51.gtf | awk '$3=="mRNA"' > transcripts.gtf
# extract all genes from gtf file for gene body analysis
cat dmel-all-r6.51.gtf |awk '$3 == "gene"' > genes.gtf
# get chrosome size
samtools faidx dmel-all-chromosome-r6.51.fasta 
cat dmel-all-chromosome-r6.51.fasta.fai | cut -f1,2 > chrom.sizes


# R code
# please run extract.R in rstudio


sort -k1,1 -k2,2n promoter.2kb.bed > promoter.2kb.sorted.bed
bedtools merge -i promoter.2kb.sorted.bed > promoter.2kb.sorted.merged.bed


