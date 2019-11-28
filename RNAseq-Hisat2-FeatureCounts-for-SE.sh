#!/bin/bash
# Using getopt

abort()
{
    echo >&2 '
***************
*** ABORTED ***
***************
'
    echo "An error occurred. Exiting..." >&2
    exit 1
}

trap 'abort' 0

set -e

##############################################################
## It assumes fastqc,trim_galore,hisat2,bedtools,samtools,RSeQC,
## featurecounts and multiQC are all available.
## We have tested this script on an Centos system.
##############################################################

threads=24
length=25

# Augument Parsing
print_usage_and_exit(){
    echo "Usage: $0 [-t <number of threads> =24 default] [-m <Galore minimum length> =$cutadapt_m] -i <hisat2 index> =path/to/hisat2/index -b <genes bed file> =path/to/genes.bed -g <genes gtf file> =path/to/genes.gtf "
    exit 1
}

while getopts ":t:m:i:b:g:h:" opt; do
    case $opt in
        j)
            threads="$OPTARG"
            echo "-j <threads used> = $threads"
            ;;
        m)
            length="$OPTARG"
            echo "-m <Galore minimum length> = $length"
            ;;
        i)
            hisat2_index="$OPTARG"
            echo "-i <hisat2 index> = $hisat2_index"
            ;;
        g)
            genes_gtf="$OPTARG"
            echo "-g <genes gtf file> = $genes_gtf"
            ;;
        b)
            genes_bed="$OPTARG"
            echo "-b <genes bed file> = $genes_bed"
            ;;
        h)
            print_usage_and_exit
            ;;
        \?)
            echo "Invalid option: -$OPTARG"
            print_usage_and_exit
            ;;
        :)
            echo "Option -$OPTARG requires an argument."
            print_usage_and_exit
            ;;
    esac
done

# Required arguments
if [ -z "$hisat2_index" -o -z "$genes_gtf" -o -z "$genes_bed" ]
then
    echo "Error: missing required argument(s)"
    print_usage_and_exit
fi

if [ -f $genes_bed ]
then
    echo "Genes bed file found! "
else
    echo "Genes bed file not exist! Please input correct Genes bed file! "
    print_usage_and_exit
fi

if [ -f $genes_gtf ]
then
    echo "Genes GTF file found! "
else
    echo "Genes GTF file not exist! Please input correct Genes GTF file! "
    print_usage_and_exit
fi

if [ -d ${hisat2_index%\/*} ]
then
    echo "Hisat2 index dir found! "
else
    echo "Hisat2 index dir not exist! Please input correct Hisat2 index dir! "
    print_usage_and_exit
fi




echo "####################################################################"
echo "RNA-Seq Pipline (Hisat2 + FeatureCounts) for SE data Starts At:"
time=$(date "+%Y-%m-%d %H:%M:%S")
echo "${time}"
echo "####################################################################"
echo ""
echo "####################################################################"
echo "########## Step 1. Analysing Sequence Quality with FastQC ##########"
echo "####################################################################"
mkdir results
cd results
mkdir 1_initial_qc
cd ../fastq

ls *.fastq.gz | while read id;
do
fastqc \
-o ../results/1_initial_qc/ \
--noextract \
${id} 
done
wait

echo "####################################################################"
echo "##### Step 2:  Removing Low Quality Sequences with Trim_Galore #####"
echo "####################################################################"
cd ../results
mkdir 2_trimmed_output
cd ../fastq

ls *.fastq.gz | while read id;
do
trim_galore \
--cores $threads \
--quality 20 \
--fastqc \
--length $length \
--output_dir ../results/2_trimmed_output/ \
$id
done

echo "####################################################################"
echo "########## Step 3: Aligning to Genome with Hisat2-aligner ##########"
echo "####################################################################"
cd ../results
mkdir 3_aligned_sequences
cd 2_trimmed_output

conda_path=$(which conda)
source ${conda_path%bin/conda}etc/profile.d/conda.sh
conda activate python2.7

ls *trimmed.fq.gz | while read id;
do
hisat2 -t -p $threads \
-x $hisat2_index \
-U $id \
-S ../3_aligned_sequences/${id%.fq.gz}.sam \
--summary-file ../3_aligned_sequences/${id%.fq.gz}.hisat2.report
done
wait

conda deactivate

cd ../3_aligned_sequences

ls *sam | while read id; do 
samtools view  -@ 12 -S $id -b | \
samtools sort -@ 12 \
> ${id%.*}.sorted.bam
done
wait

ls *sorted.bam | while read id; do 
samtools index ${id} 
done
wait

rm *.sam

echo "####################################################################"
echo "########## Step 4:  Generate BW files with bamCoverage #############"
echo "####################################################################"
cd ../
mkdir 4_bw_files
cd 3_aligned_sequences

ls *sorted.bam | while read id;
do
bamCoverage -b ${id} \
-o ${id%.sorted.bam}.RPKM.bw \
-of bigwig \
-p $threads \
--normalizeUsing RPKM
done
wait

echo "####################################################################"
echo "################## Step 5: RNAseq QC with RSeQC ####################"
echo "####################################################################"
cd ../
mkdir 5_RSeQC_report
cd 3_aligned_sequences

ls *sorted.bam | while read id;
do
bam_stat.py \
-i ${id} \
> ../5_RSeQC_report/${id%.sorted.bam}.RSeQC.bam_stat.txt

read_distribution.py \
-i ${id} -r $genes_bed \
> ../5_RSeQC_report/${id%.sorted.bam}.RSeQC.reads_distribution.txt
done
wait

echo "####################################################################"
echo "######## Step 6: Summarizing Gene Counts with featureCounts ########"
echo "####################################################################"
cd ../
mkdir 6_final_counts
cd 3_aligned_sequences

dirlist=$(ls -t ./*.bam | tr '\n' ' ')

featureCounts \
-a $genes_gtf \
-F GTF \
-o ../6_final_counts/final_counts.txt \
-g 'gene_name' \
-T 6 \
-M \
$dirlist


echo "####################################################################"
echo "######## Step 7: Generating analysis report with multiQC  ##########"
echo "####################################################################"
cd ../
mkdir 7_multiQC
cd ../

multiqc results \
--outdir results/7_multiQC


echo "####################################################################"
echo "Congratulations! Whole processs finished! "
echo "RNA-Seq Pipline Finished At:"
time=$(date "+%Y-%m-%d %H:%M:%S")
echo "${time}"
echo "####################################################################"
