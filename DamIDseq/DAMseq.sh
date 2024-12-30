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

####################################################################
## It assumes cutadapt, bowtie2 bedtools and samtools is available.
## We have tested this script with cutadapt 2.5, bowtie2 2.2.5 and 
## samtools v1.9 on an Centos 7 system.
####################################################################

##  输出结果
# 1. bam 文件及mapping logs文档目录： 2_mapping_bowtie2
# 2. bam correlation plots/data目录： 2_mapping_bowtie2
# 3. duplicated reads 统计信息 （当-r <remove duplicated reads> 为 true）目录： 2_mapping_bowtie2/rmDup_logs目录
# 4. RSeQC reports目录： 2_mapping_bowtie2/RSeQC_report
# 5. bw 文件目录：3_bw_bamCoverage
# 6. 综合报告目录：4_QC_multiQC/multiqc_report.html


## 经验总结
# 1. reads duplication 水平很高98%左右，peak很高，几乎无背景，可能是由于PCR循环数目太高导致的，本质是富集的DNA太少。
# 2. T7 外切酶和内切酶（内切酶Alw1）得到的reads起始位置不一样，内切酶不含GATC序列，并且在固定位置（GATC+5base）开始。
# 3. [DamIDseq] peak两端reads很多，中间少，应该是建库前没做超声打断，或者打断不充分（貌似DAMIDseq抓到的DNA目标片段就是很短）。强烈建议要做打断。
# 4. [DamIDseq] 正常情况，不需要不去除duplicated reads，也不需要延伸单端reads（得到的peaks更sharp，不会两边高中间低），如果未打断DNA，建议延伸单端reads。
# 5. PCR 引物序列为： GGTCGCGGCCGAGGATC， 可检查引物是否切除干净。
# 参考： 
# zcat SRR8955177_S1_L001_R1_001.fastq.gz | head -n 40000 | grep "GGTCGCGGCCGAGGATC" | wc -l # 正向
# zcat SRR8955177_S1_L001_R1_001.fastq.gz | head -n 40000 | grep "GATCCTCGGCCGCGACC" | wc -l # 反向


# 环境配置
# conda activate /home/xilab/software/miniconda-envs/bioinfo
# conda install bioconda::fastqc
# conda install bioconda::cutadapt
# conda install bioconda::bowtie2
# conda install bioconda::samtools
# conda install bioconda::deeptools
# conda install bioconda::rseqc
# conda install bioconda::picard
# conda install bioconda::multiqc
# conda install bioconda::homer # for downstream analysis
# perl /home/xilab/software/miniconda-envs/bioinfo/share/homer/.//configureHomer.pl -list
# perl /home/xilab/software/miniconda-envs/bioinfo/share/homer/.//configureHomer.pl -install dm6


# change default parameters here!
bowtie2_index=/home/xilab/reference/bowtie2_index/drosophlia/ucsc-dm6/dm6
genes_bed=/home/xilab/reference/Genome/Drosophlia/UCSC/dm6/dm6.refGene.bed
threads=24
output=results
fastq=fastq
paired=false
# test模式，也可以帮助检查是否有接头序列、引物序列啥的，本流程不包含去接头的代码
test=false
rm_dup=false
# 注意，默认SE reads会被延伸至300bp，这可能导致peaks跨越GATC位点
extend_SEReads=false

# Augument Parsing
print_usage_and_exit(){
	echo "Usage: $0 
        [-p <paired data> = $paired default]
        [-f <path of fastq files> = $fastq default]
        [-j <number of threads> = $threads default]
        [-i <bowtie2 index> = $bowtie2_index default]
        [-g <genes bed file> = $genes_bed default]
        [-o <output directory> = $output default]
        [-t <test mode> = $test default - using all reads]
        [-r <remove duplicated reads> = $rm_dup default]
        [-e <extend the SE reads to 300bp> = $extend_SEReads default]"
	exit 1
}


while getopts ":t:p:f:j:i:o:g:r:e:" opt; do
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
            echo "-b <genes bed file> = $genes_bed"
            ;;
        r)
            rm_dup="$OPTARG"
            echo "-r <remove duplicated reads> = $rm_dup"
            ;;
        e)
            extend_SEReads="$OPTARG"
            echo "-e <extend the SE reads to 300bp> = $extend_SEReads"
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
if [ -z "$output" -o -z "$bowtie2_index" -o -z "$genes_bed" ]
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


# 检查fastq文件名
## 格式为: SampleName_S1_L001_R1_001.fastq.gz
## 如不符合以上命名规则，可尝试使用测序数据管理目录下的prepare-fastqs.R脚本进行重命名
if [ -h $fastq ]
then
    # 如果fastq目录是一个软连接
    samples_fastq=($(ls -l $(readlink ${fastq}) | grep "^[dl]" | awk '{print $9}'))
else
    # fastq目录不是软连接，但该目录下的子目录均为软连接，如：
    # fastq
    # |---- AD2_1_young -> /home/xilab/ATACseq/GSE164317/AD2_1_young
    # |---- AD2_2_young -> /home/xilab/ATACseq/GSE164317/AD2_2_young
    samples_fastq=($(ls -l ${fastq} | grep "^[dl]" | awk '{print $9}'))
fi
## 列出所有fastq.gz文件，并 
echo "FASTQ检查"
for ASample in ${samples_fastq[*]};
do
    echo "找到样本: ${ASample}"
    ls ${wroking_folder}/${fastq}/${ASample}/*_001.fastq.gz | while read id
    do
        echo "fastq文件: ${id}"
    done
done
## 手动检查样本及其对应的fastq文件！
while true; do
    read -p "FASTQ检查： 请检查样本及其对应的fastq文件是否完整！文件名需为SampleName_S1_L001_R1_001.fastq.gz格式！" yn
    case $yn in
        [Yy]* ) echo "人工检查通过，继续..."; break;;
        [Nn]* ) echo "请检查你的样本及其fastq文件，退出中..."; exit 1;;
        * ) echo "yes - 检查通过； no - 有问题";;
    esac
done


## 测试模式
if [[ "$test" = true ]]
then
    cd ${wroking_folder} 
    if [ -d ./test ]; then
        echo "测试模式：删除旧的test目录..."
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
        echo "测试模式： 处理样本 - $Asample"
        mkdir -p ${target_fq}/$Asample
        ls ${source_fq}/${Asample}/*_001.fastq.gz | while read id
        do
            echo "测试模式：Subset 文件 $(basename ${id})..."
            zcat $id | head -n ${test_lines} > ${target_fq}/${Asample}/$(basename ${id%.gz})
            gzip ${target_fq}/${Asample}/$(basename ${id%.gz})
        done
    done
    fastq=test
fi


cd ${wroking_folder} && mkdir -p ./${output} && cd ./${output} && output_path=`pwd` # results dir
cd ${wroking_folder} && cd ${fastq} && samples_folder=`pwd` && samples=($(ls -l  | grep "^[dl]" | awk '{print $9}')) # samples dir and the sample names

## 生成新的样本和fastq文件的meta文件
### 表头
echo "Meta文件：生成表头"
touch ${output_path}/sample_fastq_meta.txt
if [[ "$paired" = true ]]
then
    echo -e 'SampleName\tR1\tR2' >>  ${output_path}/sample_fastq_meta.txt
else
    echo -e 'SampleName\tfastq' >>  ${output_path}/sample_fastq_meta.txt
fi
### 内容
echo "Meta文件：生成样本和FASTQ信息"
for Asample in ${samples[*]};
do
    cd ${samples_folder}/$Asample && Asample_fastqs=($(ls *.gz))
    if [[ "$paired" = true ]]
    then
        # 双端数据
        ## 要有两个fastq文件
        if [[ ${#Asample_fastqs[@]} -ne 2 ]]; then
            echo "Meta文件： 请检查您的fastq文件名是否为配对的，退出..." >&2
            exit 1
        fi
        ## 文件名要以_R1_001.fastq.gz结尾
        if [[ ${Asample_fastqs[0]/_R1_001.fastq.gz/_R2_001.fastq.gz} != ${Asample_fastqs[1]} ]]
        then
            echo "Meta文件： 请检查您的fastq文件名是否以_R1_001.fastq.gz结尾，退出......" >&2
            exit 1
        fi
        ## 保存
        echo -e $Asample'\t'${samples_folder}/$Asample/${Asample_fastqs[0]}'\t'${samples_folder}/$Asample/${Asample_fastqs[1]} >> ${output_path}/sample_fastq_meta.txt
    else
        # 单端数据
        ## 仅允许1个fastq文件
        if [[ ${#Asample_fastqs[@]} -ne 1 ]]; then
            echo "Meta文件： 单端模式，每个样本仅允许1个fastq文件，退出......" >&2
            exit 1
        fi
        ## 保存
        echo -e $Asample'\t'${samples_folder}/$Asample/${Asample_fastqs} >> ${output_path}/sample_fastq_meta.txt
    fi    
done

# some code to check duplicated sample names! not necessary for now!


echo "####################################################################"
echo "Dam-Seq Pipline (fstqc + bowtie2 + bamCoverage) for PE data Starts At:"
starttime=`date +'%Y-%m-%d %H:%M:%S'`
echo "程序开始于： "${starttime}""
echo "####################################################################"


echo ">>>>>>>>>>> Step-1: Analysing Sequence Quality with FastQC <<<<<<<<<<<"
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

# 重命名 fastqc report files 
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


rm fastq.path

echo ">>>>>>>>>> Step-2: Aligning to Genome with Bowtie2-aligner <<<<<<<<<<"
cd ${output_path} && mkdir -p ./2_mapping_bowtie2 && cd ./2_mapping_bowtie2 && bowtie2_mapping_path=`pwd`
cd ${bowtie2_mapping_path} && mkdir -p ./bowtie2_logs && bowtie2_log_path=`pwd`
cd ${output_path}

if [[ "$paired" = true ]]
then
    cat sample_fastq_meta.txt | sed '1d' | while read sample_name read1 read2;
    do
        bowtie2 \
            -p $threads \
            -x $bowtie2_index \
            -1 ${read1} \
            -2 ${read2} \
            2> ${bowtie2_log_path}/${sample_name}.log \
        | samtools view -u -@ 4 \
        | samtools sort -@ 4 > ${bowtie2_mapping_path}/${sample_name}.sorted.bam
    done
else
    cat sample_fastq_meta.txt | sed '1d' | while read sample_name read;
    do
        bowtie2 \
            -p $threads \
            -x $bowtie2_index \
            -U ${read} \
            2> ${bowtie2_log_path}/${sample_name}.log \
        | samtools view -u -@ 4 \
        | samtools sort -@ 4 > ${bowtie2_mapping_path}/${sample_name}.sorted.bam
    done
fi



echo ">>>>>>>>>> Step-3: Remove low mapping quality reads <<<<<<<<<<"
cd ${bowtie2_mapping_path}
ls *.sorted.bam | while read id;
do
samtools index $id
samtools idxstats $id >& ./${id%.sorted.bam}.idxstat.log
samtools view -S $id -q 20 -@ $threads -b > ./${id%.sorted.bam}.bam
done

rm *.sorted.bam
rm *.sorted.bam.bai

ls *.bam | while read id;
do
samtools sort ${id} -@ $threads -o ${id%.bam}.sorted.bam && samtools index ${id%.bam}.sorted.bam
rm $id
done

## 去除PCR重复的reads
if [[ "$rm_dup" = false ]]
then
    echo ">>>>>>>>>> Step-4: Remove Duplicated reads - Escaped <<<<<<<<<<"
else
    echo ">>>>>>>>>> Step-4: Remove Duplicated reads <<<<<<<<<<"
    mkdir ./rmDup_logs
    if [[ "$paired" = true ]]
    then
        ls *sorted.bam | xargs -P $threads -n 1 sh -c 'picard MarkDuplicates I=$0 O=${0%sorted.bam}bam REMOVE_DUPLICATES=true VALIDATION_STRINGENCY=LENIENT METRICS_FILE=./rmDup_logs/${0%sorted.bam}log'
    else
        ls *sorted.bam | xargs -P $(($threads/4)) -n 1 sh -c 'samtools markdup -r -@ 4 -s $0 ${0%sorted.bam}bam >& ./rmDup_logs/${0%sorted.bam}.log'
    fi
    rm *sorted.bam
    ls *.bam | while read id;
    do
        samtools sort ${id} -@ $threads -o ${id%bam}sorted.bam && samtools index ${id%bam}sorted.bam
        rm $id
    done
fi

echo ">>>>>>>>>>> Step-5:  mapping statics <<<<<<<<<<<"
bowtie2_bam_files=$(ls *sorted.bam)
bowtie2_sample_Labels=$(ls *sorted.bam | xargs -n1 sh -c 'echo ${0//.sorted.bam}')


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
    -o PCA_readCounts.png

mkdir -p ./RSeQC_report
export global_genes_bed=${genes_bed}
ls *sorted.bam | xargs -P $(($threads/2)) -n 1 \
sh -c 'bam_stat.py -i $0 > ./RSeQC_report/${0%sorted.bam}bam_stat.txt & \
read_distribution.py -i ${0} -r $global_genes_bed >  ./RSeQC_report/${0%sorted.bam}reads_distribution.txt'
export -n global_genes_bed && export -n global_rseqc_log_path=${rseqc_log_path}

if [[ "$paired" = true ]]
then
    ls *sorted.bam | xargs -P $(($threads/2)) -n 1 \
    sh -c 'picard CollectInsertSizeMetrics I=$0 O=${0%sorted.bam}insert_size_metrics.txt H=${0%sorted.bam}insert_size_histogram.pdf M=0.5'
fi

echo ">>>>>>>>>>> Step-6:  Generate BW files <<<<<<<<<<<"
cd ${output_path} && mkdir -p ./3_bw_bamCoverage && cd ./3_bw_bamCoverage && bw_path=`pwd`
cd ${bowtie2_mapping_path}

export global_bw_path=$bw_path

if [[ "$paired" = true ]]
then
    # for PE data, Reads with mates are always extended to match the fragment size defined by the two read mates.
    ls *sorted.bam | xargs -P $(($threads/6)) -n 1 \
    sh -c 'bamCoverage -b $0 -o ${global_bw_path}/${0%sorted.bam}Bowtie2.RPKM.bw -of bigwig -p 6 -bs 1 --extendReads --normalizeUsing RPKM'

    ls *sorted.bam | xargs -P $(($threads/6)) -n 1 \
    sh -c 'bamCoverage -b $0 -o ${global_bw_path}/${0%sorted.bam}Bowtie2.BPM.bw -of bigwig -p 6 -bs 1 --extendReads --normalizeUsing BPM'
else
    # for SE data, Extended reads to 300bp will improve the Peak shape.
    # 300bp: Default in damidseq_pipeline, suppose fragments size is 300bp
    # ref: https://owenjm.github.io/damidseq_pipeline/
    if [[ "$extend_SEReads" = true ]]
    then
        ls *sorted.bam | xargs -P $(($threads/6)) -n 1 \
        sh -c 'bamCoverage -b $0 -o ${global_bw_path}/${0%sorted.bam}Bowtie2.RPKM.bw -of bigwig -p 6 -bs 1 --extendReads 300 --normalizeUsing RPKM'
        ls *sorted.bam | xargs -P $(($threads/6)) -n 1 \
        sh -c 'bamCoverage -b $0 -o ${global_bw_path}/${0%sorted.bam}Bowtie2.BPM.bw -of bigwig -p 6 -bs 1 --extendReads 300 --normalizeUsing BPM'
    else
        ls *sorted.bam | xargs -P $(($threads/6)) -n 1 \
        sh -c 'bamCoverage -b $0 -o ${global_bw_path}/${0%sorted.bam}Bowtie2.RPKM.bw -of bigwig -p 6 -bs 1 --normalizeUsing RPKM'
        ls *sorted.bam | xargs -P $(($threads/6)) -n 1 \
        sh -c 'bamCoverage -b $0 -o ${global_bw_path}/${0%sorted.bam}Bowtie2.BPM.bw -of bigwig -p 6 -bs 1 --normalizeUsing BPM'
    fi
fi
export -n global_bw_path


echo ">>>>>>>>>> Step-7: Generating analysis report with multiQC <<<<<<<<<<<"
cd ${output_path} && mkdir -p ./4_QC_multiQC && cd ./4_QC_multiQC && multiqc_path=`pwd`
cd ${output_path}

multiqc -f ${output_path} --outdir ${multiqc_path} -c /home/xilab/reference/Scripts/multiqc_config.yaml

echo ">>>>>>>>>> Step-8: Delete unused files <<<<<<<<<<<"
if [[ "$test" = true ]]
then
    rm -rf ${wroking_folder}/test
fi


echo "#******************************************************************#"
echo "祝贺，任务染成啦"
endtime=`date +'%Y-%m-%d %H:%M:%S'`
echo "程序结束于： ${endtime}"
start_seconds=$(date --date="$starttime" +%s);
end_seconds=$(date --date="$endtime" +%s);
secondsdiff=$((end_seconds-start_seconds))
echo "程序总运行时间： "$((secondsdiff/3600))" h "$((secondsdiff%3600/60))" m"
trap : 0
echo >&2 '#******************************************************************#'
