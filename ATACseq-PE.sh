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

###############################################
## It assumes cutadapt, bowtie2 bedtools and samtools is 
## available.
## We have tested this script with 
## cutadapt 2.5, bowtie2 2.2.5 and samtools v1.9 
## on an Centos system.
###############################################
# change default parameters here!
bowtie2_index=/home/xilab/reference/UCSC/dm6/Sequence/Bowtie2Index/genome
blacklist=/data0/reference/encode-atac-seq-pipeline/dm6/dm6-blacklist.v2.bed
cutadapt_m=10
threads=12
output_folder=results
fastq=fastq

# Augument Parsing
print_usage_and_exit(){
	echo "Usage: $0 [-f <path of fastq files> ] [-j <number of threads> =12 default] [-m <cutadapt minimum length> =$cutadapt_m] -i <bowtie2 index> =path/to/bowtie2/index -b <blacklist bed file>=path/to/blacklist.bed -o <output directory>"
	exit 1
}

echo "****************"
echo "***for SE data**"
echo "****************"
while getopts ":f:o:j:m:i:b:h:" opt; do
    case $opt in
    	f)
            fastq="$OPTARG"
            echo "-d <fastq dir> = $fastq"
            ;;
        j)
            threads="$OPTARG"
            echo "-j <threads used> = $threads"
            ;;
        m)
            cutadapt_m="$OPTARG"
            echo "-m <cutadapt minimum length> = $cutadapt_m"
            ;;
        i)
            bowtie2_index="$OPTARG"
            echo "-i <bowtie2 index> = $bowtie2_index"
            ;;
        o)
            output_folder="$OPTARG"
            echo "-o <Output files Path> = $output_folder"
            ;;
        b)
            blacklist="$OPTARG"
            echo "-b <blacklist bed file> = $blacklist"
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
if [ -z "$output_folder" -o -z "$bowtie2_index" -o -z "$blacklist" ]
then
    echo "Error: missing required argument(s)"
    print_usage_and_exit
fi

mkdir -p $output_folder
cd $output_folder
# cutadapt remove nextera transposase sequence
echo "#******************************************************************#"
echo "ATAC-Seq Pipline (cutadapt + bowtie2) for PE data Starts At:"
starttime=`date +'%Y-%m-%d %H:%M:%S'`
echo "程序开始于： "${starttime}""
echo ""
echo "                  Step 1. Sequence QC and Trimming                  "
echo "#******************************************************************#"
mkdir -p fastqc
mkdir -p cutadapt
cd fastqc
qc_path=`pwd`
cd ../cutadapt/
trim_path=`pwd`
cd ../../${fastq}

num_fastq=`ls -l *.gz | wc -l`
fastqc -t $num_fastq *.gz -o $qc_path


ls *.gz | grep R1 | while read id;
do
fastq_name=${id%%.*}.cutadapt.fastq
cutadapt \
-a CTGTCTCTTATACACATCT \
-A CTGTCTCTTATACACATCT \
-m $cutadapt_m \
-j 12 \
-o $trim_path/${fastq_name} \
-p $trim_path/${fastq_name/R1/R2} \
>& $trim_path/${id%gz}cutadapt.log \
${id} ${id/R1/R2}
done

# bowtie2 mapping reads to reference
cd ../$output_folder
mkdir -p bowtie2
cd cutadapt/

num_fastq=`ls -l *cutadapt.fastq | wc -l`
fastqc -t $num_fastq *cutadapt.fastq -o $qc_path

echo "***********************************************************************************"
echo "Step 2: using bowtie2 to map reads to reference and using samtools to sort bamfiles"
echo "***********************************************************************************"
ls *.cutadapt.fastq | grep R1 | while read id;
do
bowtie2 \
-p $threads \
--very-sensitive \
-x $bowtie2_index \
-1 $id \
-2 ${id/R1/R2} \
2> ../bowtie2/${id%.cutadapt.fastq}.bowtie2.log \
| samtools view -u -@ 4 \
| samtools sort -@ 4 > ../bowtie2/${id%.cutadapt.fastq}.sorted.bam
done
wait

# count reads mapping to mitochondrial
cd ../bowtie2
echo "************************************************"
echo "Step 3: counting reads mapping to mitochondrial!"
echo "************************************************"
ls *.sorted.bam | while read id;
do
samtools idxstats $id >& ${id%.sorted.bam}.idxstat.log
done
wait

# remove Mitochondrial reads;sort reads
echo "**************************************************"
echo "Step 4: remove  Mitochondrial reads and sort files"
echo "**************************************************"
ls *sorted.bam | while read id;do
samtools view -@ 4 -h ${id} \
| grep -v chrM \
| samtools sort -@ 4 -O bam -o ${id%.bam}.rmChrM.bam
done
wait

rm *sorted.bam

echo "******************************"
echo "Step 5: remove PCR duplicates!"
echo "******************************"
#ls *rmChrM.bam | while read id;
#do
#samtools sort -n -o ${id%bam}namesort.bam $id &
#done
#wait
#
#ls *namesort.bam | while read id;
#do
#samtools fixmate -m $id ${id%bam}fixmate.bam &
#done
#wait
#
#ls *fixmate.bam | while read id;
#do
#samtools sort -o ${id%bam}positionsort.bam $id &
#done
#wait
#
#ls *positionsort.bam | while read id;
#do
#samtools markdup -r -@ 8 ${id} ${id%.bam}.rmDup.bam &
#done
#wait
#
#rm *namesort.bam *fixmate.bam *positionsort.bam


ls *rmChrM.bam | while read id;
do
picard MarkDuplicates I=${id} O=${id%.bam}.rmDup.bam \
REMOVE_DUPLICATES=true VALIDATION_STRINGENCY=LENIENT \
METRICS_FILE=${id%.*}.rmDup.txt
done


# Remove multi-mapped reads
echo "**********************************"
echo "Step 6: remove multi-mapped reads!"
echo "**********************************"
ls *rmDup.bam | while read id;
do
samtools view -h -@ 4 -q 30 ${id} > ${id%.bam}.rmMulti.bam
done
wait

# sort bam file and remove temporary files
echo "*************************************************"
echo "Step 7: sort bam file and remove temporary files!"
echo "*************************************************"
ls *rmMulti.bam | while read id;
do
samtools sort -@ 4 ${id}  -O bam -o ${id%.*}.sorted.bam
done
wait

rm *rmMulti.bam
rm *rmDup.bam
rm *rmChrM.bam


# index bam files and remove blacklist regions
echo "****************************************************"
echo "Step 8: index bam files and remove blacklist regions"
echo "****************************************************"
ls *rmMulti.sorted.bam | while read id; do 
samtools index ${id} -@ 4
done
wait

ls *rmMulti.sorted.bam | while read id;
do
bedtools intersect -nonamecheck -v -abam ${id} -b $blacklist > ${id%.bam}.rmBlacklist.bam
done
wait

# index bam files and generate bw files
echo "*************************************************"
echo "Sept 9: index the bam files and generate bw files"
echo "*************************************************"
ls *rmBlacklist.bam | while read id; do 
samtools index ${id} -@ 4
done
wait

rm *rmMulti.sorted.bam.bai
rm *rmMulti.sorted.bam

cd ../
mkdir -p bw
cd bowtie2
ls *rmBlacklist.bam | while read id;
do
bamCoverage -b ${id} -o ../bw/${id}.bw -of bigwig -p $threads --normalizeUsing BPM -bs 1
done
wait

# plot correlations of bam files
echo "***************************************"
echo "Step 10: plot correlations of bam files"
echo "***************************************"
num_bam=`ls -l *rmBlacklist.bam | wc -l` 
if [ $num_bam -ge 3 ]
then
	# using all BAM files in the folder
	multiBamSummary bins \
		--bamfiles *rmBlacklist.bam \
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
else
	echo "plot correlations of bam files escaped for bam files less than 3!"
fi

echo "####################################################################"
echo "######## Step 11: Generating analysis report with multiQC  ##########"
echo "####################################################################"
cd ../
mkdir -p multiQC
cd ../

multiqc $output_folder \
--outdir $output_folder/multiQC

echo "#******************************************************************#"
echo "Congratulations! Whole processs finished! "
echo "ATAC-Seq Pipline For PE Data Finished At:"
endtime=`date +'%Y-%m-%d %H:%M:%S'`
echo "${endtime}"
start_seconds=$(date --date="$starttime" +%s);
end_seconds=$(date --date="$endtime" +%s);
secondsdiff=$((end_seconds-start_seconds))
echo "程序总运行时间： "$((secondsdiff/3600))" h "$((secondsdiff%3600/60))" m"
trap : 0
echo >&2 '#******************************************************************#'
