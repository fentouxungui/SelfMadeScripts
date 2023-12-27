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

cutadapt_m=10
threads=12

# Augument Parsing
print_usage_and_exit(){
	echo "Usage: $0 [-j <number of threads> =12 default] [-m <cutadapt minimum length> =$cutadapt_m] -i <bowtie2 index> =path/to/bowtie2/index -b <blacklist bed file>=path/to/blacklist.bed -o <output directory>"
	exit 1
}

echo "****************"
echo "***for SE data**"
echo "****************"
while getopts ":j:m:i:b:o:h:" opt; do
    case $opt in
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
echo "************************************************************************"
echo "Step 1: using cutadapt to remove reads with length less than 10(default)"
echo "************************************************************************"
mkdir -p cutadapt
cd ../fastq
num_fastq=`ls -l | wc -l`
ls *.fastq | while read id;
do
cutadapt -a CTGTCTCTTATACACATCT -m 10 -j 12 \
-o ../$output_folder/cutadapt/${id%.*}.cutadapt.fastq ${id}
done

# bowtie2 mapping reads to reference
cd ../$output_folder
mkdir -p bowtie2
cd cutadapt/
echo "***********************************************************************************"
echo "Step 2: using bowtie2 to map reads to reference and using samtools to sort bamfiles"
echo "***********************************************************************************"
ls *.cutadapt.fastq | while read id;
do
bowtie2 -p $threads --very-sensitive  -x $bowtie2_index -U ${id}  \
| samtools view -u -@ $threads \
| samtools sort -@ $threads > ../bowtie2/${id}.sorted.bam
done
wait

# count reads mapping to mitochondrial
cd ../bowtie2
echo "************************************************"
echo "Step 3: counting reads mapping to mitochondrial!"
echo "************************************************"
ls *.sorted.bam | while read id;
do
touch ${id%.sorted.bam}.mito.log
samtools view $id | cut -f 3 | sort | uniq -c > ${id%.sorted.bam}.mito.log
done
wait

# remove Mitochondrial reads;sort reads
echo "**************************************************"
echo "Step 4: remove  Mitochondrial reads and sort files"
echo "**************************************************"
ls *sorted.bam | while read id;do
samtools view -@ $threads -h ${id} \
| grep -v chrM \
| samtools sort -@ $threads -O bam -o ${id%.bam}.rmChrM.bam
done
wait

# remove PCR duplicates
echo "******************************"
echo "Step 5: remove PCR duplicates!"
echo "******************************"
ls *rmChrM.bam | while read id;
do
samtools markdup -r -@ $threads ${id} ${id%.bam}.rmdup.bam
done
wait

# Remove multi-mapped reads
echo "**********************************"
echo "Step 6: remove multi-mapped reads!"
echo "**********************************"
ls *rmdup.bam | while read id;
do
samtools view -h -@ $threads -q 30 ${id}>  ${id%.bam}.rmMulti.bam
done
wait

# sort bam file and remove temporary files
echo "*************************************************"
echo "Step 7: sort bam file and remove temporary files!"
echo "*************************************************"
ls *rmMulti.bam | while read id;
do
samtools sort -@ 12 ${id}  -O bam -o ${id%.*}.sorted.bam
done
wait

rm *rmMulti.bam
rm *rmdup.bam
rm *rmChrM.bam
rm *fastq.sorted.bam

# index bam files and remove blacklist regions
echo "****************************************************"
echo "Step 8: index bam files and remove blacklist regions"
echo "****************************************************"
ls *rmMulti.sorted.bam | while read id; do 
samtools index ${id} -@ $threads
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
samtools index ${id} -@ $threads
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
# using all BAM files in the folder
multiBamSummary bins \
	--bamfiles *rmBlacklist.bam \
	--minMappingQuality 30 \
	-out readCounts.npz --outRawCounts readCounts.tab

if [ $num_fastq -ge 3 ]
then
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
fi

trap : 0
echo >&2 '
**********************************************
*** DONE Processing ...
*** You can use files ${output_folder}/bowtie2
*** ${output_folder}/cutadapt and ${output_folder}/bw
**********************************************
'
