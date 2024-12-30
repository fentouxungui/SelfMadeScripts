#! /bin/bash
source ~/.bashrc
conda activate /home/xilab/software/miniconda-envs/bioinfo
annotatePeaks.pl df.peaks.txt dm6 > df.peaks.annotated.txt -gtf /home/xilab/reference/Genome/Drosophlia/UCSC/dm6/dm6.refGene.gtf