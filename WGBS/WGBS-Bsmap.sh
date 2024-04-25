#!/bin/bash
# Using getopt

# 该流程尚未搭建完成！！！

source ~/.bashrc
conda activate py37

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

######################################################################################################################
# install softwares
# conda activate py37
# conda install -c bioconda fastqc
# conda install -c bioconda trim-galore
# conda install -c bioconda star
# conda install -c bioconda samtools
# conda install -c bioconda deeptools
# conda install -c bioconda rseqc
# conda install -c bioconda subread
# conda install -c bioconda gffcompare
# conda install -c bioconda stringtie
# conda install -c bioconda multiqc
######################################################################################################################
# fastq dir should be a link!
fastq=fastq
output=results
threads=24
reads_min_length=20
paired=false
hisat2_index=/home/xilab/reference/hisat2_index/dm6-UCSC/genome_tran/genome_tran
star_index=/home/xilab/reference/STAR_index/Drosophlia/dm6-UCSC/STARindex_100bp
genes_bed=/home/xilab/reference/Genome/Drosophlia/UCSC/dm6/dm6.refGene.bed
genes_gtf=/home/xilab/reference/Genome/Drosophlia/UCSC/dm6/dm6.refGene.gtf
rsem_index=/home/xilab/reference/RSEM_index/dm6-UCSC/dm6-UCSC
# test模式，也可以帮助检查是否有接头序列、引物序列啥的
test=false

############################### 案例
################## UCSC - dm6
# 1.1 hisat2_index 可以从hisat2官网下载 http://daehwankimlab.github.io/hisat2/download/
# 构建 Hisat2 index: dm6 genome_tran(hisat2 官网未提供)
# cp /home/xilab/reference/Genome/Drosophlia/UCSC/dm6/dm6.fa ./genome.fa
# cp /home/xilab/reference/Genome/Drosophlia/UCSC/dm6/dm6.refGene.gtf  ./genome.gtf
# hisat2_extract_splice_sites.py genome.gtf > genome.ss
# hisat2_extract_exons.py genome.gtf > genome.exon
# hisat2-build -p 32 --exon genome.exon --ss genome.ss genome.fa genome_tran
#
# 1.2 构建STAR index
# STAR \
# --runThreadN 24 \
# --runMode genomeGenerate \
# --genomeDir STARindex_100bp \
# --genomeFastaFiles /home/xilab/reference/Genome/Drosophlia/UCSC/dm6/dm6.fa \
# --sjdbGTFfile /home/xilab/reference/Genome/Drosophlia/UCSC/dm6/dm6.refGene.gtf \
# --sjdbOverhang 100 \
# --genomeSAindexNbases 12
# genomeSAindexNbases 12 为STAR推荐的针对于果蝇基因组大小的参数
# sjdbOverhang: 测序的reads长度-1， 默认是100bp，貌似默认的这个效果也不错！
#
# 1.3 构建RSEM index
# rsem-prepare-reference \
# --gtf ./dm6.refGene.gtf \
# ./dm6.fa \
# dm6-UCSC -p 24
# Augument Parsing
print_usage_and_exit(){
    echo -e "Usage: $0 
                [-p <paired data> = $paired default]
                [-j <number of threads> = $threads default] 
                [-m <Trimgalore minimum reads length> = $length default] 
                [-i <star index> = path/to/STAR-index-dir]
                [-a <hisat2 index> = path/to/Hsiat2-index-dir] 
                [-r <rsem index> = path/to/rsem-index-dir ]
                [-g <genes gtf file> = path/to/genes.gtf] 
                [-b <genes bed file> = path/to/genes.bed]             
                [-f <fastq file directory> = $fastq default]
                [-o <results directory> = $output default]
                [-t <test mode> = Default false - using all reads]

Examples:
    Drosophlia dm6:
        sh RNAseq-NEW.sh \\
        -p true
    Mouse Gencode-GRCm39-vM31:
        sh RNAseq-NEW.sh \\
        -i /home/xilab/reference/STAR_index/Mouse/Gencode-GRCm39-vM31/STARindex_100bp \\
        -a /home/xilab/reference/hisat2_index/GRCm39-Gencode-vM31/genome_tran/genome_tran \\
        -r /home/xilab/reference/RSEM_index/mm10-Gencode-GRCm39-vM31/Genecode-GRCm39-vM31 \\
        -g /home/xilab/reference/Genome/Mouse/GENCODE/GRCm39-vM31/gencode.vM31.primary_assembly.annotation.gtf \\
        -b /home/xilab/reference/Genome/Mouse/GENCODE/GRCm39-vM31/gencode.vM31.primary_assembly.annotation.bed \\
        -p true     
"
    exit 1
}

while getopts ":t:p:j:m:i:a:r:g:b:f:o:h:" opt; do
    case $opt in
        t)
            test=$"$OPTARG"
            echo "-t <Test Mode> = $test"
            ;;
        p)
            paired="$OPTARG"
            echo "-p <Paired Data> = $paired"
            ;;
        j)
            threads="$OPTARG"
            echo "-j <threads used> = $threads"
            ;;
        m)
            reads_min_length="$OPTARG"
            echo "-m <Trimgalore minimum reads length> = $reads_min_length"
            ;;
        i)
            star_index="$OPTARG"
            echo "-i <STAR index> = $star_index"
            ;;
        a)
            hisat2_index="$OPTARG"
            echo "-a <Hisat2 index> = $hisat2_index"
            ;;
        r)
            rsem_index="$OPTARG"
            echo "-a <Rsem index> = $rsem_index"
            ;;
        g)
            genes_gtf="$OPTARG"
            echo "-g <genes gtf file> = $genes_gtf"
            ;;
        b)
            genes_bed="$OPTARG"
            echo "-b <genes bed file> = $genes_bed"
            ;;
        f)
            fastq="$OPTARG"
            echo "-f <fastq dir> = $fastq"
            ;;
        o)
            output="$OPTARG"
            echo "-o <output dir> = $output"
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
    echo "Genes bed file found!"
else
    echo "Genes bed file not exist! Please input correct Genes bed file!"
    print_usage_and_exit
fi

if [ -f $genes_gtf ]
then
    echo "Genes GTF file found!"
else
    echo "Genes GTF file not exist! Please input correct Genes GTF file!"
    print_usage_and_exit
fi

if [ -d ${star_index%\/*} ]
then
    echo "STAR index dir found! "
else
    echo "STAR index dir not exist! Please input correct STAR index dir!"
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


echo ">>>>>>>>>>> Step-0: Prepare the Analysis <<<<<<<<<<<"
# Define the dirs
wroking_folder=`pwd`
## remove old results
if [ -d ./${output} ]
then
    echo ">>> Removing the old analysis results..."
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

## test mode
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

cd ${wroking_folder}
mkdir -p ./${output} && cd ./${output} && output_path=`pwd` # results dir
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
echo "RNA-Seq Pipline Starts At:"
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



echo ">>>>>>>>>> Step-2:  Removing Low Quality Sequences with Trim_Galore <<<<<<<<<<"
cd ${output_path} && mkdir -p ./2_trimmed_output && cd ./2_trimmed_output && trim_path=`pwd` # trimed reads results dir
cd ${output_path}

# problem: trim_galore report naming is not adjustable.
# $(($threads/6)) # echo $((14/6)) return 2
# Rename the cutadapt report files
export global_trim_path=$trim_path && export global_reads_min_length=$reads_min_length
if [[ "$paired" = true ]]
then
    cat sample_fastq_meta.txt | sed '1d' | xargs -P $(($threads/6)) -n 3 sh -c \
    'trim_galore --cores 6 --fastqc --quality 20 --length $global_reads_min_length --paired $1 $2 --basename $0 --output_dir $global_trim_path &&
    mv ${global_trim_path}/`basename ${1}`_trimming_report.txt ${global_trim_path}/${0}_R1_trimming_report.txt &&
    mv ${global_trim_path}/`basename ${2}`_trimming_report.txt ${global_trim_path}/${0}_R2_trimming_report.txt'
else
    cat sample_fastq_meta.txt | sed '1d' | xargs -P $(($threads/6)) -n 2 sh -c \
    'trim_galore --cores 6 --fastqc --quality 20 --length $global_reads_min_length --basename $0 --output_dir $global_trim_path $1 &&
    mv ${global_trim_path}/`basename ${1}`_trimming_report.txt ${global_trim_path}/${0}_trimming_report.txt'
fi
export -n global_trim_path && export -n global_reads_min_length


echo ">>>>>>>>>> Step-3:  mapping reads with Bsmap <<<<<<<<<<"
# bsmap 最新更新为2013年。
cd ${output_path} && mkdir -p ./3_aligned_BSMAP && cd ./3_aligned_BSMAP && bsmap_mapping_path=`pwd`
cd ${trim_path}

if [[ "$paired" = true ]]
then
    ls *_val_1.fq.gz | while read id;
    do
        bsmap \
        -a ${id} \
        -b ${id/_val_1.fq.gz/_val_2.fq.gz} \
        -d /home/xilab/zhangyc/DataAnalysis/chicken/chicken-reference/galGal6.fa \
        -q 20 \
        -f 5 \
        -p 24 \
        -w 20 \
        -r 0 \
        -v 0.05 \
        -s 16 \
        -S 1 \
        -n 0 \
        2> ${bsmap_mapping_path}/${id/_val_1.fq.gz/.BSMAP_report.txt} | \
        samtools view -b \
        -o ${bsmap_mapping_path}/${id/_val_1.fq.gz/.aligned.bam}
    done
fi
# 20min for 80m reads

cd ${bsmap_mapping_path}
ls *.aligned.bam | while read id
do
    sambamba sort \
    -m 32GB \
    --tmpdir tmp \
    -t 24 \
    -o ${id/.aligned.bam/_sorted.bam} $id
done
# -m: 所有线程的存储上限
# --tmpdir: 临时文件夹
# -t: 线程数
# -o: 输出文件的路径名

echo ">>>>>>>>>> Remove dulicated reads <<<<<<<<<<"
ls *_sorted.bam | while read id
do
    sambamba markdup \
    --overflow-list-size 3000000 \
    --tmpdir tmp \
    -t 12 \
    ${id} \
    ${id/_sorted.bam/_DupRemoved.bam} \
    2> ${id/_sorted.bam/.MarkDup_report.txt}
done

echo ">>>>>>>>>> Step-4:  methylation calling with MethylDackel <<<<<<<<<<"
# 提取甲基化信息 per-base methylation metrics
# https://github.com/dpryan79/MethylDackel
# MethylDackel (formerly named PileOMeth, which was a temporary name derived due to it using a PILEup to extract METHylation metrics) 
# will process a coordinate-sorted and indexed BAM or CRAM file containing some form of BS-seq alignments and extract per-base methylation 
# metrics from them. MethylDackel requires an indexed fasta file containing the reference genome as well.
cd ${output_path} && mkdir -p ./4_methylation_MethylDackel && cd ./4_methylation_MethylDackel && MethylDackel_path=`pwd`
cd ${bsmap_mapping_path}
ls *_DupRemoved.bam | while read id
do
    MethylDackel extract \
    /home/xilab/zhangyc/DataAnalysis/chicken/chicken-reference/galGal6.fa \
    ${id} --opref ${MethylDackel_path}/${id/_DupRemoved.bam}
done
# 输出的bedgraph需要转为bw！
# 此外还可计算CHH and CHG metrics
# MethylDackel mbias 可用于绘图
# --opref: 要保存文件的前缀
# 如果只是这样后面不再加其他参数，则默认提取CpG位点，输出文件为'sample_CpG.bedGraph'
# 如果后面再加上'--CHG'，则为提取CHG位点，输出文件为'sample_CHG.bedGraph'
# 同理：'--CHH' 对应 'sample_CHH.bedGraph'
############################ bsmap
# -a: read1
# -b: read2
# -d: 基因组文件
# -q: 质量阈值，低于该值的舍弃
# -f: 去除低质量的reads，将包含Ns > n (n指设的值)的reads舍弃
# -p: 参与处理的处理器个数
# -r: 如何报告重复的Hits，不知道啥意思
# -v: 所容许的错配率
# -s: seed大小，WGBS模式为16，RRBS模式为12，最小8最大16
# -S: 用于随机数生成的种子在选择多个命中时使用其他种子值根据读取的索引编号生成伪随机数，以允许再现映射结果。默认值为0。（从系统时钟获取种子，映射结果不可重现。）翻译的，同不知道啥意思
# -n: set mapping strand information
# -n  [0,1]   set mapping strand information:
#             -n 0: only map to 2 forward strands, i.e. BSW(++) and BSC(-+)    (i.e. the "Lister protocol")
#             for PE sequencing, map read#1 to ++ and -+, read#2 to +- and --. 
#             -n 1: map SE or PE reads to all 4 strands, i.e. ++, +-, -+, --    (i.e. the "Cokus protocol")
# 2> 标准错误重定向，在这里是保存BSMAP的报告内容
# samtools view -b : 将SAM文件输出为BAM文件
# -o: 输出的BAM文件名
# If use Hiseq2000, -z should be set 64;

echo ">>>>>>>>>> Step-8: Generating analysis report with multiQC <<<<<<<<<<<"
cd ${output_path} && mkdir -p ./7_multiQC && cd ./7_multiQC && multiqc_path=`pwd`
multiqc -f ${output_path} --outdir ${multiqc_path} -c /home/xilab/reference/Scripts/config_example.yaml

# cat /home/xilab/reference/Scripts/config_example.yaml
# extra_fn_clean_exts:
#   - .gz
#   - .fastq
#   - .fq
#   - .bam
#   - .sam
#   - .sra
#   - _tophat
#   - _trimming_report.txt
#   - .cutadapt_R1.fastq
#   - .cutadapt_R2.fastq
#   - _RSeQC.reads_distribution
#   - _RSeQC.junction_saturation
#   - _RSeQC.junction_annotation
#   - _star_aligned
#   - _fastqc
#   - _R1
#   - _R2
#   - type: remove
#     pattern: ".sorted"
#   - type: regex
#     pattern: '^Sample_\d+'
#   - type: regex_keep
#     pattern: "[A-Z]{3}[1-9]{4}[A,B][1-9]"

# use_filename_as_sample_name:
#   - cutadapt
# table_columns_visible:
#   QualiMap: False

echo ">>>>>>>>>> Step-9: Delete unused files <<<<<<<<<<<"
cd ${output_path}
cd $trim_path
echo "Delete the trimmed fastqs..."
rm *.fq.gz

echo "Delete the Hisat2 generaged sam files..."
cd ${hisat2_mapping_path}
rm *.sam

if [[ "$test" = true ]]
then
    # rm -rf $fastq
    rm -rf ${wroking_folder}/test
fi


cd ${wroking_folder}
echo "#******************************************************************#"
echo "Congratulations! Whole processs finished! "
echo "RNA-Seq Pipline Finished At:"
endtime=`date +'%Y-%m-%d %H:%M:%S'`
echo "${endtime}"
start_seconds=$(date --date="$starttime" +%s);
end_seconds=$(date --date="$endtime" +%s);
secondsdiff=$((end_seconds-start_seconds))
echo "程序总运行时间： "$((secondsdiff/3600))" h "$((secondsdiff%3600/60))" m"
trap : 0
echo >&2 '#******************************************************************#'
