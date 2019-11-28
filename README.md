# Scripts

Scripts used in My lab!

##  AtacSeqPipeline-SE.sh
This script is a simple atacseq-pipeline for Single-end data of Drosophlia! put all your uncompressed fastq files in a directory named 'fastq'! put ``AtacSeqPipeline-SE.sh`` under the same directory with fastq dir.

**using**

``` sh
bash AtacSeqPipeline-SE.sh \
  -i /path/to/drosophlia/Bowtie2Index/genome \
  -b /path/to/dm6-blacklist.v2.bed \
  -o output-folder-name
```

## RNAseq-Hisat2-FeatureCounts-for-SE.sh

Pipeline: Hisat2 + FeatureCounts + RSeQC + MultiQC, more info please refer to [twbattaglia/RNAseq-workflow](https://github.com/twbattaglia/RNAseq-workflow)

Suitable for Single end data! put all your compressed fastq files in a directory named 'fastq'! put ``RNAseq-Hisat2-FeatureCounts-for-SE.sh`` under the same directory with fastq dir.

**using**

``` sh
sh RNAseq-Hisat2-FeatureCounts-for-SE.sh \
-i /data0/reference/hisat2_index/D.melanogaster.BDGP6.Genome/bdgp6/genome \
-b /data0/reference/ensembl/Drosophila_melanogaster/BDGP6/Annotation/Archives/archive-2015-07-23-16-41-33/Genes/genes.bed \
-g /data0/reference/ensembl/Drosophila_melanogaster/BDGP6/Annotation/Archives/archive-2015-07-23-16-41-33/Genes/genes.gtf \
```

## Chipseeker.r

For Drosophlia!

input: MACS2 opuput (narrowpeaks and summits.bed)

output: plots(feature distribution, DistributionRelativeToTSS, GO-BP, KEGG), Annotated peaks, extracted peaks for Homer(downstream Motif analysis)

**using**

Before using, pleas change the Rscript PATH is shebang: `` which Rscript``.

``` sh
./chipseeker.r -p xxx.narrowPeak \
               -s xxx.summits.bed \
               -n prefix
```      

## shinyapp-scRNAseq-DataSearch.R

This is the R code for building the Fly Gut EEs's scRNAseq Data Website: https://xilab.shinyapps.io/database/, where you can search the genes expression pattern. Thanks to https://www.shinyapps.io/ for holding this app.
