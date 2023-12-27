# Scripts

Scripts used in My lab!

## ATACseq/Cut&Tag
###  AtacSeqPipeline-SE.sh
This script is a simple atacseq-pipeline for Single-end data of Drosophlia! put all your uncompressed fastq files in a directory named 'fastq'! put ``AtacSeqPipeline-SE.sh`` under the same directory with fastq dir. for example:  ``sample1.fastq, sample2.fastq...`` .

**using**:

``bash AtacSeqPipeline-SE.sh -h`` will print the usage.

``` sh
bash AtacSeqPipeline-SE.sh \
  -i /path/to/Bowtie2Index/genome \
  -b /path/to/blacklist.bed \
  -j threads \
  -m minimal-reads-length-for-cutadapt \
  -o output-folder-name
```

###  ATACseq-PE.sh
all the paired reads should put into a same dir, make sure "R1" and "R2" is the only difference for a paired compressed fastqs.
such as ``sample1_R1.fastq.gz`` , ``sample1_R2.fastq.gz``,``sample2_R1.fastq.gz`` , ``sample2_R2.fastq.gz``...
outputs： in the ``multiqc`` dir, you can find the QC report for each step!

**using**:

``bash ATACseq-PE.sh -h`` will print the usage.

``` sh
bash ATACseq-PE.sh \
  -i /path/to/Bowtie2Index/genome \
  -b /path/to/blacklist.bed \
  -f fastq-folder-name \
  -j threads \
  -m minimal-reads-length-for-cutadapt \
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

more info please refer to [twbattaglia/RNAseq-workflow](https://github.com/twbattaglia/RNAseq-workflow)

### RNAseq-Hisat2/STAR-FeatureCounts-for-SE.sh

Pipeline-1: Hisat2 + FeatureCounts + RSeQC + MultiQC

Pipeline-2: STAR + FeatureCounts + RSeQC + MultiQC

Suitable for Single end data! put all your compressed fastq files in a directory named 'fastq'! put ``RNAseq-Hisat2-FeatureCounts-for-SE.sh`` or ``RNAseq-STAR-FeatureCounts-for-SE.sh`` under the same directory with fastq dir.

**Attention**:

1. for Pipeline-1: Because Hisat2 do not support python3.7, I installed Hisat2 in Conda environment python2.7. So before using this script, Please change ``conda activate python2.7`` to your hisat2 conda environment!

2. Pipeline-2: less than 20 fastq files are used for a computer with 48 threads!

**using**:

``` sh
sh RNAseq-Hisat2-FeatureCounts-for-SE.sh \
-i /path/to/hisat2_indes-or-STAR_index \
-b /path/to/genes.bed \
-g /path/to/genes.gtf
```

### RNAseq-STAR-FeatureCounts-PE.sh

put all your paired reads into same dir, file names should like: ``sample1_R1.fastq.gz`` , ``sample1_R2.fastq.gz``,``sample2_R1.fastq.gz`` , ``sample2_R2.fastq.gz``...

**using**:

``` sh
sh RNAseq-Hisat2-FeatureCounts-for-SE.sh \
-f /path/to/fastqDir \
-j threads \
-m minimal-reads-length-for-cutadapt \
-i /path/to/hisat2_indes-or-STAR_index \
-b /path/to/genes.bed \
-g /path/to/genes.gtf \
-o output-folder-name
```

## shinyapp-scRNAseq-DataSearch.R

This is the R code for building the Fly Gut EEs's scRNAseq Data Website: https://xilab.shinyapps.io/database/, where you can search the genes expression pattern. Thanks to https://www.shinyapps.io/ for holding this app.
