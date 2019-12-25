# Scripts

Scripts used in My lab!

## ATACseq/Cut&Tag
###  AtacSeqPipeline-SE.sh
This script is a simple atacseq-pipeline for Single-end data of Drosophlia! put all your uncompressed fastq files in a directory named 'fastq'! put ``AtacSeqPipeline-SE.sh`` under the same directory with fastq dir.

**using**:

``` sh
bash AtacSeqPipeline-SE.sh \
  -i /path/to/drosophlia/Bowtie2Index/genome \
  -b /path/to/dm6-blacklist.v2.bed \
  -o output-folder-name
```
###  ATACseq-PE.sh
**using**:

``` sh
bash AtacSeqPipeline-SE.sh \
  -i /path/to/drosophlia/Bowtie2Index/genome \
  -b /path/to/dm6-blacklist.v2.bed \
  -f fastqdir \
  -o output-folder-name
```
### Chipseeker.r

For Drosophlia!

input: MACS2 opuput (narrowpeaks and summits.bed)

output: plots(feature distribution, DistributionRelativeToTSS, GO-BP, KEGG), Annotated peaks, extracted peaks for Homer(downstream Motif analysis)

**using**:

Before using, pleas change the Rscript PATH in shebang: `` which Rscript``.

``` sh
./chipseeker.r -p xxx.narrowPeak \
               -s xxx.summits.bed \
               -n prefix
```      

## RNAseq
### RNAseq-Hisat2/STAR-FeatureCounts-for-SE.sh

Pipeline-1: Hisat2 + FeatureCounts + RSeQC + MultiQC

Pipeline-2: STAR + FeatureCounts + RSeQC + MultiQC

more info please refer to [twbattaglia/RNAseq-workflow](https://github.com/twbattaglia/RNAseq-workflow)

Suitable for Single end data! put all your compressed fastq files in a directory named 'fastq'! put ``RNAseq-Hisat2-FeatureCounts-for-SE.sh`` or ``RNAseq-STAR-FeatureCounts-for-SE.sh`` under the same directory with fastq dir.

**Attention**:

1. for Pipeline-1: Because Hisat2 do not support python3.7, I installed Hisat2 in Conda environment python2.7. So before using this script, Please change ``conda activate python2.7`` to your hisat2 conda environment!

2. Pipeline-2: less than 20 fastq files are used for a computer with 48 threads!

**using**:

``` sh
sh RNAseq-Hisat2-FeatureCounts-for-SE.sh \
-i /path/to/hisat2_index_or_STAR_index \
-b /path/to/genes.bed \
-g /path/to/genes.gtf \
```

## shinyapp-scRNAseq-DataSearch.R

This is the R code for building the Fly Gut EEs's scRNAseq Data Website: https://xilab.shinyapps.io/database/, where you can search the genes expression pattern. Thanks to https://www.shinyapps.io/ for holding this app.
