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

################################################################
## It assumes fastqc,trim_galore,STAR,bedtools,samtools,RSeQC,
## featurecounts and multiQC are all available.
## We have tested this script on an Centos system.
################################################################

threads=24
length=25

# Augument Parsing
print_usage_and_exit(){
    echo "Usage: $0 [-t <number of threads> =24 default] [-m <Galore minimum length> =$cutadapt_m] -i <star index> =path/to/STAR/index -b <genes bed file> =path/to/genes.bed -g <genes gtf file> =path/to/genes.gtf "
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
            star_index="$OPTARG"
            echo "-i <STAR index> = $star_index"
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
if [ -z "$star_index" -o -z "$genes_gtf" -o -z "$genes_bed" ]
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

if [ -d ${star_index%\/*} ]
then
    echo "STAR index dir found! "
else
    echo "STAR index dir not exist! Please input correct STAR index dir! "
    print_usage_and_exit
fi




echo "####################################################################"
echo "RNA-Seq Pipline (STAR + FeatureCounts) for SE data Starts At:"
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

fastqc *.gz \
-t $threads \
-o ../results/1_initial_qc/ \
--noextract \
${id} 
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
echo "########### Step 3: Aligning to Genome with STAR-aligner ###########"
echo "####################################################################"
cd ../results
mkdir 3_aligned_STAR
cd 3_aligned_STAR
star_path=$(pwd)
cd ../2_trimmed_output

ls *trimmed.fq.gz | while read id;
do
STAR \
--readFilesCommand zcat \
--genomeDir $star_index \
--readFilesIn $id  \
--runThreadN 16 \
--outSAMtype BAM SortedByCoordinate \
--quantMode GeneCounts \
--outFileNamePrefix $star_path/${id%trimmed.fq.gz}
done
wait

cd ../3_aligned_STAR

ls *.bam | while read id; do 
samtools index ${id} 
done
wait


echo "####################################################################"
echo "########## Step 4:  Generate BW files with bamCoverage #############"
echo "####################################################################"
cd ../
mkdir 4_bw_files
cd 3_aligned_STAR

ls *.bam | while read id;
do
bamCoverage -b ${id} \
-o ../4_bw_files/${id%.bam}.RPKM.bw \
-of bigwig \
-p $threads \
--normalizeUsing RPKM
done
wait

num_bam=`ls -l *.bam | wc -l`

if [ $num_bam -ge 3 ]
then
    multiBamSummary bins \
    --bamfiles *.bam \
    --minMappingQuality 30 \
    -out readCounts.npz --outRawCounts readCounts.tab

    plotCorrelation \
        -in readCounts.npz \
        --corMethod pearson --skipZeros \
        --plotTitle "Pearson Correlation of Read Counts" \
        --whatToPlot scatterplot \
        -o scatterplot_PearsonCorr_BAMScores.png   \
        --outFileCorMatrix PearsonCorr_BAMScores.tab

    plotCorrelation \
        -in readCounts.npz \
        --corMethod spearman --skipZeros \
        --plotTitle "Spearman Correlation of Read Counts" \
        --whatToPlot scatterplot \
        -o scatterplot_SpearmanCorr_BAMScores.png   \
        --outFileCorMatrix SpearmanCorr_BAMScores.tab

    plotCorrelation \
        -in readCounts.npz \
        --corMethod spearman --skipZeros \
        --plotTitle "Spearman Correlation of Read Counts" \
        --whatToPlot heatmap --colorMap RdYlBu --plotNumbers \
        -o heatmap_SpearmanCorr_readCounts.png   \
        --outFileCorMatrix SpearmanCorr_readCounts.tab

    plotCorrelation \
        -in readCounts.npz \
        --corMethod pearson --skipZeros \
        --plotTitle "Pearson Correlation of Read Counts" \
        --whatToPlot heatmap --colorMap RdYlBu --plotNumbers \
        -o heatmap_PearsonCorr_readCounts.png   \
        --outFileCorMatrix PearsonCorr_readCounts.tab
fi &

echo "####################################################################"
echo "################## Step 5: RNAseq QC with RSeQC ####################"
echo "####################################################################"
cd ../
mkdir 5_RSeQC_report
cd 3_aligned_STAR

ls *.bam | while read id;
do
bam_stat.py \
-i ${id} \
> ../5_RSeQC_report/${id%.bam}.RSeQC.bam_stat.txt &
done

ls *.bam | while read id;
do
read_distribution.py \
-i ${id} -r $genes_bed \
> ../5_RSeQC_report/${id%.bam}.RSeQC.reads_distribution.txt &
done

ls *.bam | while read id;
do
echo ${id} >> bam.files.txt
done

geneBody_coverage.py -r $genes_bed \
-i bam.files.txt -o ../5_RSeQC_report/GeneBodyCoverage &


task=bam_stat.py
while ps -ef | grep $task | grep -v 'grep'; do
sleep 10
done
echo "bam_stat.py finished! "


task=read_distributi
while ps -ef | grep $task | grep -v 'grep'; do
sleep 10
done
echo "read_distributiion.py finished! "


task=geneBody_cove
while ps -ef | grep $task | grep -v 'grep'; do
sleep 10
done
echo "read_distributiion.py finished! "

rm bam.files.txt

task=multiBamSummary
while ps -ef | grep $task | grep -v 'grep'; do
sleep 10
done
echo "multiBamSummary finished! "

task=plotCorrelation
while ps -ef | grep $task | grep -v 'grep'; do
sleep 10
done
echo "plotCorrelation finished! "

echo "####################################################################"
echo "######## Step 6: Summarizing Gene Counts with featureCounts ########"
echo "####################################################################"
cd ../
mkdir 6_final_counts
cd 3_aligned_STAR

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
trap : 0
echo >&2 '####################################################################'
