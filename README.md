# Scripts
Scripts used in My lab!

##  AtacSeqPipeline-SE.sh
This script is a simple atacseq-pipeline for Single-end data of Drosophlia! put all your uncompressed fastq files in a directory named 'fastq'! put AtacSeqPipeline-SE.sh under the same directory with fastq dir.

**using**

``` sh
bash AtacSeqPipeline-SE.sh \
  -i /path/to/drosophlia/Bowtie2Index/genome \
  -b /path/to/dm6-blacklist.v2.bed \
  -o: output folder name.
```

## Chipseeker.r
For Drosophlia!
input: MACS2 opuput (narrowpeaks and summits.bed)
output: plots(feature distribution, DistributionRelativeToTSS, GO-BP, KEGG), Annotated peaks, extracted peaks for Homer(downstream Motif analysis)

**using**

``` sh
./chipseeker.r -p xxx.narrowPeak \
               -s xxx.summits.bed \
               -n prefix
```      
