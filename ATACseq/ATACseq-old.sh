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
# Analysis Roadmap and Softwares/Functions needed:
# 1. FASTQ QC: fastqc;
# 2. Trim Adapter: Cutadapt;
# 3. Align Reads: Bowtie2
# 4. Remove reads mapped to Mitochondria: Samtools and grep
# 5. Remove PCR duplicates By picard MarkDuplicates(PE) or Samtools markdup(SE)
# 6. Remove multi-mapped reads By samtools view
# 7.1 Remove blacklist regions By bedtools intersect and plot correlation of bam files By multiBamSummary + plotCorrelation 
# 7.2 plot insert size by bamPEFragmentSize(PE) and picard CollectInsertSizeMetrics
# 7.3 Bam QC by bam_stat.py and read_distribution.py
# 8. Shift bam file by alignmentSieve and Generage BW files By bamCoverage
# 9. Analysis Report: multiqc

# 案例
# 果蝇
# sh ATACseq.sh \
# -t true \
# -p true
# 小鼠
# sh ATACseq.sh \
# -i  /data0/reference/bowtie2_index/mm10/mm10 \
# -g /data0/reference/UCSC/mm10/mm10.knownGene.bed \
# -b /data0/reference/encode-atac-seq-pipeline/mm10/mm10.blacklist.bed \
# -p true
# 人
# sh ATACseq.sh \
# -i /data0/reference/encode-atac-seq-pipeline/hg38/bowtie2_index/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta \
# -g /home/xilab/reference/encode-atac-seq-pipeline/hg38/hg38.knownGene.bed \
# -b /data0/reference/encode-atac-seq-pipeline/hg38/hg38.blacklist.bed \
# -p true \
# -t true

# 备注:
# alignmentSieve shift reads 感觉有些问题，被丢掉很多reads，没有核实是否是improper paired reads。
# bamCoverage 有些成对reads无法被正常延伸。

# install softwares
# conda activate py37
# conda install -c bioconda fastqc
# conda install -c bioconda cutadapt
# conda install -c bioconda bowtie2
# conda install -c bioconda samtools
# conda install -c bioconda deeptools
# conda install -c bioconda rseqc
# conda install -c bioconda multiqc
###############################################
# change default parameters here!
bowtie2_index=/home/xilab/reference/UCSC/dm6/Sequence/Bowtie2Index/genome
blacklist=/data0/reference/encode-atac-seq-pipeline/dm6/dm6-blacklist.v2.bed
genes_bed=/data0/reference/STAR_index/Drosophlia/BDGP6_ensembl/genes.bed
reads_min_length=25
threads=24
output=results
fastq=fastq
paired=false
# test模式，也可以帮助检查是否有接头序列、引物序列啥的
test=false

# Augument Parsing
print_usage_and_exit(){
	echo "Usage: $0 
        [-p <paired data> = $paired default]
        [-f <path of fastq files> = $fastq]
        [-j <number of threads> = $threads default]
        [-m <cutadapt minimum length> =$reads_min_length]
        [-i <bowtie2 index> = $bowtie2_index]
        [-g <genes bed file> = $genes_bed]    
        [-b <blacklist bed file> = $blacklist]
        [-o <output directory> = $output default]
        [-t <test mode> = Default false - using all reads]"
	exit 1
}


while getopts ":t:g:p:f:o:j:m:i:b:h:" opt; do
    case $opt in
        t)
            test=$"$OPTARG"
            echo "-t <Test Mode> = $test"
            ;;
        p)
            paired="$OPTARG"
            echo "-p <Paired Data> = $paired"
            ;;
    	f)
            fastq="$OPTARG"
            echo "-f <fastq dir> = $fastq"
            ;;
        j)
            threads="$OPTARG"
            echo "-j <threads used> = $threads"
            ;;
        m)
            reads_min_length="$OPTARG"
            echo "-m <cutadapt minimum length> = $reads_min_length"
            ;;
        i)
            bowtie2_index="$OPTARG"
            echo "-i <bowtie2 index> = $bowtie2_index"
            ;;
        o)
            output="$OPTARG"
            echo "-o <Output files Path> = $output"
            ;;
        g)
            genes_bed="$OPTARG"
            echo "-g <genes bed file> = $genes_bed"
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
if [ -z "$output" -o -z "$bowtie2_index" -o -z "$blacklist" -o -z "$genes_bed" ]
then
    echo "Error: missing required argument(s)"
    print_usage_and_exit
fi

if [ -d ${fastq%\/*} ]
then
    echo "fastq file directory found! "
else
    echo "fastq file directory not exist! Please input correct fastq file directory!"
    print_usage_and_exit
fi


# functions
function plot_bam_function(){
    plotCorrelation \
        -in readCounts.npz \
        --corMethod pearson --skipZeros \
        --plotTitle "Pearson Correlation of Read Counts" \
        --whatToPlot scatterplot \
        -o scatterplot_PearsonCorr_BAMScores.png   \
        --removeOutliers \
        --outFileCorMatrix PearsonCorr_BAMScores.tab &
    plotCorrelation \
        -in readCounts.npz \
        --corMethod spearman --skipZeros \
        --plotTitle "Spearman Correlation of Read Counts" \
        --whatToPlot scatterplot \
        --removeOutliers \
        -o scatterplot_SpearmanCorr_BAMScores.png   \
        --outFileCorMatrix SpearmanCorr_BAMScores.tab &
    plotCorrelation \
        -in readCounts.npz \
        --corMethod spearman --skipZeros \
        --plotTitle "Spearman Correlation of Read Counts" \
        --whatToPlot heatmap --colorMap RdYlBu --plotNumbers \
        -o heatmap_SpearmanCorr_readCounts.png   \
        --removeOutliers \
        --outFileCorMatrix SpearmanCorr_readCounts.tab &
    plotCorrelation \
        -in readCounts.npz \
        --corMethod pearson --skipZeros \
        --plotTitle "Pearson Correlation of Read Counts" \
        --whatToPlot heatmap --colorMap RdYlBu --plotNumbers \
        -o heatmap_PearsonCorr_readCounts.png   \
        --removeOutliers \
        --outFileCorMatrix PearsonCorr_readCounts.tab
}

# Define the dirs
wroking_folder=`pwd`
## remove old results
if [ -d ./${output} ]
then
    echo "Removing the old analysis results..."
    rm -rf ./${output}
fi

# check and rename fastqs: SampleName_R1_001.fastq.gz/SampleName_R1_001.fq.gz or 
# SampleName_1.fq.gz/_R1.fq.gz/_1.fastq.gz 
# to SampleName_R1.fastq.gz
if [ -h $fastq ]
then
    # fastq is a link file
    samples_fastq=($(ls -l $(readlink ${fastq}) | grep "^[dl]" | awk '{print $9}'))
else
    # fastq in not a link file, and the subdir are links
    # fastq
    # |---- AD2_1_young -> /home/xilab/ATACseq/GSE164317/AD2_1_young
    # |---- AD2_2_young -> /home/xilab/ATACseq/GSE164317/AD2_2_young
    samples_fastq=($(ls -l ${fastq} | grep "^[dl]" | awk '{print $9}'))
fi


for ASample in ${samples_fastq[*]};
do
    ls ${wroking_folder}/${fastq}/${ASample}/*.gz | while read id
    do
        new_name_a=${id/.fq.gz/.fastq.gz}
        new_name_b=${new_name_a/_001.fastq.gz/.fastq.gz}
        new_name_c=${new_name_b/_1.fastq.gz/_R1.fastq.gz}
        new_name_d=${new_name_c/_2.fastq.gz/_R2.fastq.gz}
        if [[ "$id" != "$new_name_d" ]]; then
            mv ${id} ${new_name_d}
            echo ">>> Renamed the fastq file: from $(basename ${id}) to $(basename ${new_name_d})..."
        fi
    done
done

## test mode
if [[ "$test" = true ]]
then
    cd ${wroking_folder} 
    if [ -d ./test ]; then
        echo ">>> Removing the old test dir..."
        rm -rf test
    fi
    mkdir -p test && cd test
    source_fq=${wroking_folder}/${fastq}
    target_fq=`pwd`
    test_lines=100000
    if [ -h $source_fq ]
    then
        samples_fq=($(ls -l $(readlink ${source_fq}) | grep "^[dl]" | awk '{print $9}'))
    else
        samples_fq=($(ls -l ${source_fq} | grep "^[dl]" | awk '{print $9}'))
    fi
    for Asample in ${samples_fq[*]};
    do
        echo $Asample
        mkdir -p ${target_fq}/$Asample
        ls ${source_fq}/${Asample}/*q.gz | while read id
        do
            echo ">>> Subseting fastq: $(basename ${id})..."
            zcat $id | head -n ${test_lines} > ${target_fq}/${Asample}/$(basename ${id%.gz})
            gzip ${target_fq}/${Asample}/$(basename ${id%.gz})
        done
    done
    fastq=test
fi

cd ${wroking_folder}
mkdir -p ./${output} && cd ./${output} && output_path=`pwd` # results dir
cd ${wroking_folder} && cd ${fastq} && samples_folder=`pwd` && samples=($(ls -l  | grep "^[dl]" | awk '{print $9}')) # samples dir and the sample names



## generage new meta.data file
### tittle
touch ${output_path}/sample_fastq_meta.txt
if [[ "$paired" = true ]]
then
    echo -e 'SampleName\tR1\tR2' >>  ${output_path}/sample_fastq_meta.txt
else
    echo -e 'SampleName\tfastq' >>  ${output_path}/sample_fastq_meta.txt
fi
### content
for Asample in ${samples[*]};
do
    # echo $Asample
    cd ${samples_folder}/$Asample && Asample_fastqs=($(ls *.gz))
    if [[ "$paired" = true ]]
    then
        # for paired data
        ## should have 2 fastq files
        if [[ ${#Asample_fastqs[@]} -ne 2 ]]; then
            echo "An error occurred, Please check the if fastqs are paired. Exiting..." >&2
            exit 1
        fi
        ## file names should be sampleName_R1.fastq.gz sampleName_R2.fastq.gz format.
        if [[ ${Asample_fastqs[0]/_R1/_R2} != ${Asample_fastqs[1]} ]]
        then
            echo "An error occurred, Fastqs should be sampleName_R1.fastq.gz and sampleName_R2.fastq.gz. Exiting..." >&2
            exit 1
        fi
        ## content
        echo -e $Asample'\t'${samples_folder}/$Asample/${Asample_fastqs[0]}'\t'${samples_folder}/$Asample/${Asample_fastqs[1]} >> ${output_path}/sample_fastq_meta.txt
    else
        # for single data
        ## only one fastq file allowed
        if [[ ${#Asample_fastqs[@]} -ne 1 ]]; then
            echo "An error occurred, Only one fastq is allowed for each sample. Exiting..." >&2
            exit 1
        fi
        ## content
        echo -e $Asample'\t'${samples_folder}/$Asample/${Asample_fastqs} >> ${output_path}/sample_fastq_meta.txt
    fi    
done


# some code to check duplicated sample names! not necessary for now!

echo "####################################################################"
echo "ATAC-Seq Pipline Starts At:"
starttime=`date +'%Y-%m-%d %H:%M:%S'`
echo "程序开始于： "${starttime}""
echo "####################################################################"

echo ">>>>>>>>>>> Step 1. Analysing Sequence Quality with FastQC <<<<<<<<<<<"
cd ${output_path} && mkdir -p ./1_initial_qc && cd ./1_initial_qc && qc_path=`pwd` # fastq qc results dir
cd ${output_path} && touch fastq.path
if [[ "$paired" = true ]]
then
    sed '1d' sample_fastq_meta.txt |  awk '{print $1"_R1\t"$2}' >> fastq.path
    sed '1d' sample_fastq_meta.txt |  awk '{print $1"_R2\t"$3}' >> fastq.path
else
    sed '1d' sample_fastq_meta.txt |  awk '{print $1"\t"$2}' >> fastq.path
fi
cut -f 2 fastq.path | xargs fastqc -o $qc_path -t $threads

# Rename the fastqc report files 
# For SE, Change the Filename to sample_fastqc.zip/sample/fastqc_data.txt
# For PE, Change the Filename to sample_R1_fastqc.zip  sample_R2_fastqc.zip
cat fastq.path | while read sample_new_name fastq_path
do
    old_name=`basename $fastq_path`
    if [[ "${old_name%.fastq.gz}" != "$sample_new_name" ]]; then
        unzip ${qc_path}/${old_name/.fastq.gz/_fastqc.zip} -d ${qc_path} && rm ${qc_path}/${old_name/.fastq.gz/_fastqc.zip}
        mv ${qc_path}/${old_name/.fastq.gz/_fastqc} ${qc_path}/${sample_new_name}_fastqc
        cd ${qc_path}/${sample_new_name}_fastqc && sed -i "s/${old_name}/${sample_new_name}.fastq.gz/g" fastqc_data.txt && cd ../
        # echo ${old_name} ${sample_new_name}.fastq.gz
        zip -r ${sample_new_name}_fastqc.zip ${sample_new_name}_fastqc && rm -rf ${sample_new_name}_fastqc
        # mv ${qc_path}/${old_name/.fastq.gz/_fastqc.zip} ${qc_path}/${sample_new_name}_fastqc.zip
        mv ${qc_path}/${old_name/.fastq.gz/_fastqc.html} ${qc_path}/${sample_new_name}_fastqc.html
        cd ${output_path}
    fi
done



echo ">>>>>>>>>> Step 2:  Removing Adapter and Low Quality Sequences with Cutadapt <<<<<<<<<<"
cd ${output_path} && mkdir -p ./2_trimmed_output && cd ./2_trimmed_output && trim_path=`pwd` # trimed reads results dir
cd ${output_path}

export global_trim_path=$trim_path && export global_reads_min_length=$reads_min_length
if [[ "$paired" = true ]]
then
    cat sample_fastq_meta.txt | sed '1d' | xargs -P $(($threads/6)) -n 3 sh -c \
    'cutadapt -a CTGTCTCTTATACACATCT -A CTGTCTCTTATACACATCT -m $global_reads_min_length -j 6 -o $global_trim_path/${0}.cutadapt_R1.fastq -p $global_trim_path/${0}.cutadapt_R2.fastq >& $global_trim_path/${0}.cutadapt.log $1 $2'
else
    cat sample_fastq_meta.txt | sed '1d' | xargs -P $(($threads/6)) -n 2 sh -c \
    'cutadapt -a CTGTCTCTTATACACATCT -m $global_reads_min_length -j 6 -o $global_trim_path/${0}.cutadapt.fastq >& $global_trim_path/${0}.cutadapt.log $1'
fi
export -n global_trim_path && export -n global_reads_min_length

cd $trim_path
ls *.fastq | xargs fastqc -o ./ -t $threads




# in case cutadapt  not output equal reads in paired fastq.
# ls *.gz | grep '_R1'| while read id;
# do
# trim_galore \
# --cores $threads \
# --fastqc \
# --nextera \
# --length $length \
# --output_dir $trim_path/$abpath \
# --paired $id ${id/_R1/_R2}
# done

echo ">>>>>>>>>> Step 3: Aligning to Genome with Bowtie2-aligner <<<<<<<<<<"
cd ${output_path} && mkdir -p ./3_aligned_Bowtie2 && cd ./3_aligned_Bowtie2 && bowtie2_mapping_path=`pwd`
cd ${trim_path}

# some problem in threads assign. bowtie2 and samtools
if [[ "$paired" = true ]]
then
    ls *.cutadapt_R1.fastq | while read id;
    do
        bowtie2 \
            -p $threads \
            --very-sensitive \
            -x $bowtie2_index \
            -1 $id \
            -2 ${id/_R1/_R2} \
            2> ${bowtie2_mapping_path}/${id%.cutadapt_R1.fastq}.bowtie2.log \
        | samtools view -u -@ 4 \
        | samtools sort -@ 4 > ${bowtie2_mapping_path}/${id%.cutadapt_R1.fastq}.sorted.bam
    done
else
    ls *cutadapt.fastq | while read id;
    do
        bowtie2 \
            -p $threads \
            --very-sensitive \
            -x $bowtie2_index \
            -U $id \
            2> ${bowtie2_mapping_path}/${id%.cutadapt.fastq}.bowtie2.log \
        | samtools view -u -@ 4 \
        | samtools sort -@ 4 > ${bowtie2_mapping_path}/${id%.cutadapt.fastq}.sorted.bam
    done
fi


# ls *R1_val_1.fq.gz | grep _1 | while read id;
# do
# samplep=${id/_1/_2}
# bowtie2 \
# -p $threads \
# --very-sensitive \
# -x $bowtie2_index \
# -1 $id \
# -2 ${samplep/_R1/_R2} \
# 2> ../bowtie2/${id%_R1_val_1.fq.gz}.bowtie2.log \
# | samtools view -u -@ 4 \
# | samtools sort -@ 4 > ../bowtie2/${id%_R1_val_1.fq.gz}.sorted.bam
# done
# wait

echo ">>>>>>>>>> Step 4: Counting and Remove reads mapped to Mitochondrial <<<<<<<<<<"
cd ${bowtie2_mapping_path}
ls *bam | xargs -P $(($threads/4)) -n 1 sh -c 'samtools index $0 -@ 4'
ls *.sorted.bam | xargs -P $threads -n 1 sh -c 'samtools idxstats $0 >& ./${0%.sorted.bam}.idxstat.log'

ls *.sorted.bam | xargs -P $(($threads/9)) -n 1 sh -c 'samtools view -@ 4 -h $0 | grep -v chrM | samtools sort -@ 4 -O bam -o ./${0%.sorted.bam}.rmChrM.bam'

mkdir raw_bam && mv *sorted.bam *sorted.bam.bai ./raw_bam
ls *.bam | xargs -P $(($threads/4)) -n 1 sh -c 'samtools index $0 -@ 4'

echo ">>>>>>>>>> Step 5: Remove PCR duplicates By picard MarkDuplicates or Samtools markdup <<<<<<<<<<"
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

# 貌似这一步的picard MarkDuplicates多线程有些问题，貌似picard MarkDuplicates本身就会调用多线程？
if [[ "$paired" = true ]]
then
    ls *rmChrM.bam | xargs -P $threads -n 1 sh -c 'picard MarkDuplicates I=$0 O=${0%.bam}.rmDup.bam REMOVE_DUPLICATES=true VALIDATION_STRINGENCY=LENIENT METRICS_FILE=${0%.rmChrM.bam}.rmDup.log'
else
    ls *rmChrM.bam | xargs -P $(($threads/4)) -n 1 sh -c 'samtools markdup -r -@ 4 -s $0 ${0%bam}rmDup.bam >& ${0%rmChrM.bam}rmDup.log'
fi

echo ">>>>>>>>>> Step 6: Remove multi-mapped reads and Remove temporary files <<<<<<<<<<"
ls *rmDup.bam | xargs -P $(($threads/4)) -n 1 sh -c 'samtools view -h -@ 4 -q 30 $0 > ${0%.bam}.rmMulti.bam'
ls *rmMulti.bam | xargs -P $(($threads/4)) -n 1 sh -c 'samtools sort -@ 4 $0 -O bam -o ${0%.bam}.sorted.bam'


rm *rmMulti.bam
rm *rmDup.bam
rm *rmChrM.bam

echo ">>>>>>>>>> Step 7: Remove blacklist regions and plot correlation of bam files <<<<<<<<<<"
ls *rmMulti.sorted.bam | xargs -P $(($threads/4)) -n 1 sh -c 'samtools index ${0} -@ 4'
export global_blacklist=$blacklist
ls *rmMulti.sorted.bam | xargs -P $threads -n 1 sh -c 'bedtools intersect -nonamecheck -v -abam $0 -b ${global_blacklist} > ${0%.bam}.rmBlacklist.bam'
export -n global_blacklist


ls *rmBlacklist.bam | xargs -P $(($threads/4)) -n 1 sh -c 'samtools index $0 -@ 4'

rm *rmMulti.sorted.bam.bai
rm *rmMulti.sorted.bam


bowtie2_bam_files=$(ls *rmChrM.rmDup.rmMulti.sorted.rmBlacklist.bam)
bowtie2_sample_Labels=$(ls *rmChrM.rmDup.rmMulti.sorted.rmBlacklist.bam | xargs -n1 sh -c 'echo ${0//.rmChrM.rmDup.rmMulti.sorted.rmBlacklist.bam}')


if [[ "$paired" = true ]]
then
    multiBamSummary bins \
        --bamfiles $bowtie2_bam_files \
        --minMappingQuality 30 \
        --labels $bowtie2_sample_Labels \
        --numberOfProcessor $threads \
        --outFileName readCounts.npz \
        --outRawCounts readCounts.tab \
        --extendRead

    bamPEFragmentSize \
        -hist fragmentSize.png \
        -T "Fragment size of PE data" \
        --maxFragmentLength 1000 \
        -b $bowtie2_bam_files \
        --numberOfProcessors $threads \
        --samplesLabel $bowtie2_sample_Labels \
        --table fragmentSize.table.txt \
        --outRawFragmentLengths fragmentSize.raw.txt


else
    multiBamSummary bins \
        --bamfiles $bowtie2_bam_files \
        --minMappingQuality 30 \
        --labels $bowtie2_sample_Labels \
        --numberOfProcessor $threads \
        --outFileName readCounts.npz \
        --outRawCounts readCounts.tab
fi

bowtie2_bam_files_number=`ls *.bam | wc -l`
if [ $bowtie2_bam_files_number -ge 3 ]
then
    plot_bam_function
fi

plotPCA  \
    -in readCounts.npz \
    --plotTitle "PCA plot of Read Counts" \
    -o PCA_readCounts.pdf


export global_genes_bed=${genes_bed}
ls *.rmChrM.rmDup.rmMulti.sorted.rmBlacklist.bam | xargs -P $(($threads/2)) -n 1 \
sh -c 'bam_stat.py -i $0 > ./${0%rmChrM.rmDup.rmMulti.sorted.rmBlacklist.bam}RSeQC.bam_stat.txt & \
read_distribution.py -i ${0} -r $global_genes_bed >  ./${0%rmChrM.rmDup.rmMulti.sorted.rmBlacklist.bam}RSeQC.reads_distribution.txt'
export -n global_genes_bed

if [[ "$paired" = true ]]
then
ls *.rmChrM.rmDup.rmMulti.sorted.rmBlacklist.bam | xargs -P $(($threads/2)) -n 1 \
sh -c 'picard CollectInsertSizeMetrics I=$0 O=${0%rmChrM.rmDup.rmMulti.sorted.rmBlacklist.bam}insert_size_metrics.txt H=${0%rmChrM.rmDup.rmMulti.sorted.rmBlacklist.bam}insert_size_histogram.pdf M=0.5'
fi

echo ">>>>>>>>>> Step 8: Generage BW files and Shift bam files <<<<<<<<<<"
# some potential problems
# inproper pairs will be removed
# 但是我去看IGV的时候，有些看起来propered reads 也被去掉了，不知道咋回事，对alignmentSieve的结果保留疑问！
ls *rmBlacklist.bam | while read id;
do
alignmentSieve -b $id -o ${id%bam}Shifted.bam --ATACshift --numberOfProcessor $threads
# the bed file can be directly used for macs2
alignmentSieve -b $id -o ${id%bam}Shifted.bed --ATACshift --numberOfProcessor $threads --BED
done



ls *Shifted.bam | xargs -P $(($threads/8)) -n 1 sh -c  'samtools view -u -@ 4 $0 | samtools sort -@ 4 > ${0%bam}Sorted.bam'
rm *Shifted.bam
ls *Sorted.bam | xargs -P $(($threads/4)) -n 1 sh -c 'samtools index $0 -@ 4'

cd ${output_path} && mkdir -p ./4_bw_Bowtie2 && cd ./4_bw_Bowtie2 && bw_path=`pwd`
cd ${bowtie2_mapping_path}

if [[ "$paired" = true ]]
then
    ls *rmBlacklist.bam | while read id;
    do
        bamCoverage -b ${id} -o ${bw_path}/${id%bam}RPKM.bw -of bigwig -p $threads --extendReads --normalizeUsing RPKM -bs 1 >& ${id%rmChrM.rmDup.rmMulti.sorted.rmBlacklist.bam}bamCoverage.RPKM.log
        bamCoverage -b ${id} -o ${bw_path}/${id%bam}BPM.bw -of bigwig -p $threads --extendReads --normalizeUsing BPM -bs 1 >& ${id%rmChrM.rmDup.rmMulti.sorted.rmBlacklist.bam}bamCoverage.BPM.log
    done
else
    ls *rmBlacklist.bam | while read id;
    do
        bamCoverage -b ${id} -o ${bw_path}/${id%bam}RPKM.bw -of bigwig -p $threads --normalizeUsing RPKM -bs 1 >& ${id%rmChrM.rmDup.rmMulti.sorted.rmBlacklist.bam}bamCoverage.log
        bamCoverage -b ${id} -o ${bw_path}/${id%bam}BPM.bw -of bigwig -p $threads --normalizeUsing BPM -bs 1 >& ${id%rmChrM.rmDup.rmMulti.sorted.rmBlacklist.bam}bamCoverage.BPM.log
    done
fi    

echo ">>>>>>>>>> Step 9: Generating analysis report with multiQC <<<<<<<<<<<"
cd ${output_path} && mkdir -p ./5_multiQC && cd ./5_multiQC && multiqc_path=`pwd`
multiqc -f ${output_path} --outdir ${multiqc_path} -c /home/xilab/reference/Scripts/config_example.yaml

echo ">>>>>>>>>> Step-8: Delete unused files <<<<<<<<<<<"
if [[ "$test" = true ]]
then
    # rm -rf $fastq
    rm -rf ${wroking_folder}/test
fi

cd $trim_path && rm *.fastq
cd ${wroking_folder}

echo "#******************************************************************#"
echo "Congratulations! Whole processs finished! "
echo "ATAC-Seq Pipline Finished At:"
endtime=`date +'%Y-%m-%d %H:%M:%S'`
echo "${endtime}"
start_seconds=$(date --date="$starttime" +%s);
end_seconds=$(date --date="$endtime" +%s);
secondsdiff=$((end_seconds-start_seconds))
echo "程序总运行时间： "$((secondsdiff/3600))" h "$((secondsdiff%3600/60))" m"
trap : 0
echo >&2 '#******************************************************************#'
