#!/bin/bash
# Using getopt

source ~/.bashrc
conda activate /home/xilab/software/miniconda-envs/bioinfo

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

# 注意
# 1. 无论双端单端，reads长度最好大于50bp，因为rsem-calculate-expression --seed-length 默认25bp，短于25bp的会警告。
# 2. 使用流程时，注意fastq文件大小不要大于2.8G！详见Debug信息第一条。
# 3. Hisat2 比对很慢很慢，可能是有大量reads比对到多个位置，估计是有大量核糖体mRNA。

######################################################################################################################
# Analysis Roadmap
# 1. STAR + FeatureCounts (counts matrix - Gene level)     > DESeq2/edgeR/limma（counts）
# 2. STAR + RSEM (counts, TPM, FPKM matrix - Gene level)   > DESeq2/edgeR/limma（counts）
#                                                                                |-  Known             > DESeq2/edgeR/limma（counts）
# 3. Hisat2 + StringTie (counts - Gene/Transcript）, TPM, FPKM matrix - Gene)   -| 
#                                                                                |-  Known and de novo > DESeq2/edgeR/limma （counts）or balgown
# 注： 1. ./results/6_Counts_StringTie/Make-Reference/gffcmp.annotated.gtf 为从头方式得到的全部转录本和外显子
# 2. Known and de novo 的结果可考虑导入到 IsoformSwitchAnalyzeR R包进行分析。未尝试！ 不同处理下，isform的切换。
# 3. Known and de novo 的结果包含以MSTRG为前缀的新基因/转录本，Known的没有，这一点分析时需要注意，对于Known and de novo 的结果可以用IsoformSwitchAnalyzeR导入，
#   参考（https://www.biostars.org/p/282817/），说是可以解决50%的MSTRG基因的注释，然后再导出count matrix，用于DESeq2/edgeR/limma分析。
#
# 之后考虑搭建balgown分析流程！
#
# Softwares/Functions needed:
# 1. FASTQ QC: fastqc;
# 2. Trim Adapter: Trimgalore;
# 3. Align reads: STAR; Hisat2; Samtools;
# 4. BW files: deepTools: bamCoverage, plotCorrelation, multiBamSummary; 
# 5. RNAseq QC: RSeQC: bam_stat.py, read_distribution.py, geneBody_coverage.py(Not Used);
# 6. Summarizing Gene/Transcrpt/Exon Counts: Subread: featureCounts; GffCompare + StringTie;
# 7. Analysis Report: multiqc


# install softwares
# conda activate py37
# conda install bioconda::fastqc
# conda install bioconda::trim-galore
# conda install bioconda::star
# conda install bioconda::hisat2
# conda install bioconda::samtools
# conda install bioconda::deeptools
# conda install bioconda::rseqc
# conda install bioconda::qualimap
# conda install bioconda::subread
# conda install bioconda::gffcompare
# conda install bioconda::stringtie==2.2.1
# stringtie==2.2.1 because of https://github.com/gpertea/stringtie/issues/238 
# conda install bioconda::rsem
# conda install bioconda::multiqc

# check package avaliable and list the version
fastqc --version
trim_galore --version
STAR --version
hisat2 --version
samtools --version
deeptools --version
# rseqc version
bam_stat.py --version
qualimap --version
# subread version
subread-align -v
gffcompare --version
stringtie --version
# rsem version
rsem-calculate-expression --version
multiqc --version


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
merge_stringtie_results_code=/data0/reference/Scripts/RNAseq/merge-stringtie-results.R
merge_RSEM_results_code=/data0/reference/Scripts/RNAseq/merge-RSEM-results.R

# test模式，也可以帮助检查是否有接头序列、引物序列啥的
test=false

############################### 果蝇案例
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


############################### 小鼠案例
############# GENCODE - GRCm39 vM31
# 1.1 hisat2_index 可以从hisat2官网下载，Hisat2官网未提供GRCm39(mm39)的index
# 构建Hisat2 index
# 参考： http://daehwankimlab.github.io/hisat2/howto/#building-indexes
# cp /home/xilab/reference/Genome/Mouse/GENCODE/GRCm39-vM31/GRCm39.genome.fa ./genome.fa
# cp /home/xilab/reference/Genome/Mouse/GENCODE/GRCm39-vM31/gencode.vM31.primary_assembly.annotation.gtf  ./genome.gtf
# hisat2_extract_splice_sites.py genome.gtf > genome.ss
# hisat2_extract_exons.py genome.gtf > genome.exon
# hisat2-build -p 32 --exon genome.exon --ss genome.ss genome.fa genome_tran
#
# 1.2 构建STAR index
# ln -s /home/xilab/reference/Genome/Mouse/GENCODE/GRCm39-vM31/gencode.vM31.primary_assembly.annotation.gtf 
# ln -s /home/xilab/reference/Genome/Mouse/GENCODE/GRCm39-vM31/GRCm39.primary_assembly.genome.fa ./
# STAR \
# --runThreadN 24 \
# --runMode genomeGenerate \
# --genomeDir STARindex_100bp \
# --genomeFastaFiles ./GRCm39.primary_assembly.genome.fa \
# --sjdbGTFfile ./gencode.vM31.primary_assembly.annotation.gtf \
# --sjdbOverhang 100
# sjdbOverhang: 测序的reads长度-1， 默认是100bp，貌似默认的这个效果也不错！
#
# 1.3 构建RSEM index
# rsem-prepare-reference \
# --gtf /home/xilab/reference/Genome/Mouse/GENCODE/GRCm39-vM31/gencode.vM31.primary_assembly.annotation.gtf  \
# /home/xilab/reference/Genome/Mouse/GENCODE/GRCm39-vM31/GRCm39.primary_assembly.genome.fa \
#  Genecode-GRCm39-vM31 -p 24
#
############################### 人案例
############### GENCODE - GRCh38.p13.release42 Primary
# 1.1 hisat2_index 可以从hisat2官网下载
# http://daehwankimlab.github.io/hisat2/download/
# 1.2 RSEM indx
# rsem-prepare-reference \
# --gtf /data0/reference/Genome/Human/GENCODE/GRCh38.p13.release42/gencode.v42.primary_assembly.annotation.gtf \
# /data0/reference/Genome/Human/GENCODE/GRCh38.p13.release42/GRCh38.primary_assembly.genome.fa \
#  GENCODE-GRCh38-p13-release42-Primary -p 24
# 1.3 STAR index
# path：/data0/reference/Genome/Human/GENCODE/GRCh38.p13.release42
# cd /data0/reference/STAR_index/Human/GENCODE-GRCh38-p13-release42-primary
# ln -s /data0/reference/Genome/Human/GENCODE/GRCh38.p13.release42/gencode.v42.primary_assembly.annotation.gtf ./
# ln -s /data0/reference/Genome/Human/GENCODE/GRCh38.p13.release42/GRCh38.primary_assembly.genome.fa ./
# STAR \
# --runThreadN 24 \
# --runMode genomeGenerate \
# --genomeDir STARindex_100bp \
# --genomeFastaFiles ./GRCh38.primary_assembly.genome.fa \
# --sjdbGTFfile ./gencode.v42.primary_assembly.annotation.gtf \
# --sjdbOverhang 100
# sjdbOverhang: 测序的reads长度-1， 默认是100bp，貌似默认的这个效果也不错！
################# UCSC - T2T-CHM13-v2.0-hs1
# 2.1 Hisat2 index
# cp /home/xilab/reference/Genome/Human/UCSC/T2T-CHM13-v2.0-hs1/hs1.fa ./genome.fa
# cp /home/xilab/reference/Genome/Human/UCSC/T2T-CHM13-v2.0-hs1/catLiftOffGenesV1.gtf  ./genome.gtf
# hisat2_extract_splice_sites.py genome.gtf > genome.ss
# hisat2_extract_exons.py genome.gtf > genome.exon
# hisat2-build -p 32 --exon genome.exon --ss genome.ss genome.fa genome_tran
#
# 2.2 RSEM indx
# rsem-prepare-reference \
# --gtf /home/xilab/reference/Genome/Human/UCSC/T2T-CHM13-v2.0-hs1/catLiftOffGenesV1.gtf \
# /home/xilab/reference/Genome/Human/UCSC/T2T-CHM13-v2.0-hs1/hs1.fa \
#  Human-UCSC-T2T-CHM13-v2.0-hs1 -p 24
# 2.3 STAR index
# path：/data0/reference/Genome/Human/UCSC/T2T-CHM13-v2.0-hs1
# cd /data0/reference/STAR_index/Human/UCSC-T2T-CHM13-v2.0-hs1
# ln -s /home/xilab/reference/Genome/Human/UCSC/T2T-CHM13-v2.0-hs1/catLiftOffGenesV1.gtf ./
# ln -s /home/xilab/reference/Genome/Human/UCSC/T2T-CHM13-v2.0-hs1/hs1.fa ./
# STAR \
# --runThreadN 24 \
# --runMode genomeGenerate \
# --genomeDir STARindex_100bp \
# --genomeFastaFiles ./hs1.fa \
# --sjdbGTFfile ./catLiftOffGenesV1.gtf \
# --sjdbOverhang 100

# 补充分析
# 1. reads 比对率低，判定是污染了什么物种，注意others，一般是的是接头序列，引物序列等。
# cd results
# mkdir S_Kraken2
# cat sample_fastq_meta.txt | sed '1d' | while read sample fastq
# do
# kraken2 --db /data0/reference/MetaGenomics/kraken2/kraken2_db/k2_pluspf_fly_mouse_20220501 \
# --threads 24 --use-names $fastq --report ./S_Kraken2/$(basename ${fastq}).kraken2.report \
# > ./S_Kraken2/$(basename ${fastq}).kraken2.output
# done
# 再写个代码，提取top10的结果，或者图形展示。

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
    Human GENCODE-GRCh38-p13-release42:
        sh RNAseq-NEW.sh \\
        -i /data0/reference/STAR_index/Human/GENCODE-GRCh38-p13-release42-primary/STARindex_100bp \\
        -a /home/xilab/reference/hisat2_index/grch38/genome_tran/genome_tran \\
        -r /home/xilab/reference/RSEM_index/Human-GENCODE-GRCh38-p13-release42-primary/GENCODE-GRCh38-p13-release42-Primary \\
        -g /data0/reference/Genome/Human/GENCODE/GRCh38.p13.release42/gencode.v42.primary_assembly.annotation.gtf \\
        -b /data0/reference/Genome/Human/GENCODE/GRCh38.p13.release42/gencode.v42.primary_assembly.annotation.bed \\
        -p true
    Human UCSC-T2T-CHM13-v2.0-hs1:
        sh RNAseq-NEW.sh \\
        -i /data0/reference/STAR_index/Human/UCSC-T2T-CHM13-v2.0-hs1/STARindex_100bp \\
        -a /home/xilab/reference/hisat2_index/T2T-CHM13-v2.0-hs1-UCSC/genome_tran \\
        -r /home/xilab/reference/RSEM_index/Human-UCSC-T2T-CHM13-v2.0-hs1/Human-UCSC-T2T-CHM13-v2.0-hs1 \\
        -g /home/xilab/reference/Genome/Human/UCSC/T2T-CHM13-v2.0-hs1/catLiftOffGenesV1.gtf \\
        -b /home/xilab/reference/Genome/Human/UCSC/T2T-CHM13-v2.0-hs1/catLiftOffGenesV1.bed \\
        -p true
    Supplementary Analysis:
    1. rRNA contamination rates
        cd ./resulsts && mkdir -p S_bbmap_bbduk
        # PE data
        cat sample_fastq_meta.txt | sed '1d' | while read samplename r1 r2
        do
            ~/software/BBTools/bbmap/bbduk.sh \\
            in1=\${r1} \\
            in2=\${r2} \\
            stats=./S_bbmap_bbduk/\${samplename}.stats \\
            ref=/home/xilab/reference/Biology-Database/rRNA_Database/BDGP6.32_rRNA_Unspliced_cDNA.fa
        done
            #out=./8_rDNA/myolA-HNF4-IR1_val_R1_rRNAremoved.fq.gz \\ # 输出去除污染后的reads
            #out2=./8_rDNA/myolA-HNF4-IR1_val_R2_rRNAremoved.fq.gz \\
        # SE data
        cat sample_fastq_meta.txt | sed '1d' | while read samplename fastq
        do
            ~/software/BBTools/bbmap/bbduk.sh \\
            in=\${fastq} \\
            stats=./S_bbmap_bbduk/\${samplename}.stats \\
            ref=/home/xilab/reference/Biology-Database/rRNA_Database/BDGP6.32_rRNA_Unspliced_cDNA.fa
        done
        # 或者基于counts矩阵判定rRNA的污染比例，但是计算出来的值很低！可能是这种方法不适合？multimapped reads不被计数？使得比例偏低？

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



echo ">>>>>>>>>> Step-3: Aligning to Genome <<<<<<<<<<"
echo ">>>>>>>>>> Step-3.1: Aligning to Genome with STAR-aligner <<<<<<<<<<"
cd ${output_path} && mkdir -p ./3_aligned_STAR && cd ./3_aligned_STAR && star_mapping_path=`pwd`
cd ${trim_path}

# -limitBAMsortRAM 30000000000 \ # in case of large bam files
# --quantMode TranscriptomeSAM will output alignments translated into  tran-script coordinates,为了使用RSEM 进行定量分析做准备；
# If you pass one file, STAR will consider these as single-end reads: --readFilesIn single_reads.fastq.
# If you pass two files, STAR will consider these as paired reads: --readFilesIn pair_1.fastq pair_2.fastq.

if [[ "$paired" = true ]]
then
    ls *_val_1.fq.gz | while read id;
    do
        STAR \
            --readFilesCommand zcat \
            --genomeDir ${star_index} \
            --readFilesIn ${id} ${id/_val_1.fq.gz/_val_2.fq.gz} \
            --runThreadN $threads \
            --outSAMtype BAM SortedByCoordinate \
            --quantMode TranscriptomeSAM GeneCounts \
            --limitBAMsortRAM 30000000000 \
            --outFileNamePrefix ${star_mapping_path}/${id/_val_1.fq.gz/.}
        done
else
    ls *trimmed.fq.gz | while read id;
    do
        STAR \
            --readFilesCommand zcat \
            --genomeDir $star_index \
            --readFilesIn $id  \
            --runThreadN $threads \
            --outSAMtype BAM SortedByCoordinate \
            --quantMode TranscriptomeSAM GeneCounts \
            --limitBAMsortRAM 30000000000 \
            --outFileNamePrefix ${star_mapping_path}/${id%trimmed.fq.gz}
    done
fi


cd ${star_mapping_path} && ls *Aligned.sortedByCoord.out.bam | xargs -P $threads -n 1 samtools index
# Aligned.toTranscriptome.out.bam can not be indexed!

echo ">>>>>>>>>>> Step-3.2: Aligning to Genome with Hisat2-aligner <<<<<<<<<<<"
cd ${output_path} && mkdir -p ./3_aligned_Hisat2 && cd ./3_aligned_Hisat2 && hisat2_mapping_path=`pwd`
cd ${trim_path}

if [[ "$paired" = true ]]
then
    ls *_val_1.fq.gz | while read id;
    do
        hisat2 \
            -t --dta \
            -p $threads \
            -x ${hisat2_index} \
            -1 ${id} \
            -2 ${id/_val_1.fq.gz/_val_2.fq.gz} \
            -S ${hisat2_mapping_path}/${id%_val_1.fq.gz}.sam \
            --summary-file ${hisat2_mapping_path}/${id%_val_1.fq.gz}.hisat2.report
    done
else
    ls *trimmed.fq.gz | while read id
    do
        echo "Mapping reads to reference:" ${id%_trimmed.fq.gz}
        hisat2 \
            -t --dta \
            -p $threads \
            -x ${hisat2_index} \
            -U $id \
            -S ${hisat2_mapping_path}/${id%_trimmed.fq.gz}.sam \
            --summary-file ${hisat2_mapping_path}/${id%_trimmed.fq.gz}.hisat2.report
    done
fi


cd ${hisat2_mapping_path}
ls *.sam | xargs -P $threads -n 1 \
sh -c 'samtools view $0 -F 4 -Su | samtools sort -T ${0%.sam}.accepted_hits -o ${0%.sam}.accepted_hits.bam && samtools index ${0%.sam}.accepted_hits.bam'


echo ">>>>>>>>>>> Step-4:  Bamfiles correlation Plots <<<<<<<<<<<"
# STAR bam files
cd ${star_mapping_path}
star_bam_files=$(ls *Aligned.sortedByCoord.out.bam)
star_sample_Labels=$(ls *Aligned.sortedByCoord.out.bam | xargs -n1 sh -c 'echo ${0//.Aligned.sortedByCoord.out.bam}')
star_bam_files_number=`ls *.bam | wc -l`

multiBamSummary bins \
    --bamfiles $star_bam_files \
    --minMappingQuality 30 \
    --labels $star_sample_Labels \
    --numberOfProcessor $threads \
    --outFileName readCounts.npz \
    --outRawCounts readCounts.tab


if [ $star_bam_files_number -ge 3 ]
then
    plot_bam_function
    plotPCA  \
    -in readCounts.npz \
    --plotTitle "PCA plot of Read Counts" \
    -o PCA_readCounts.png
fi



# Hisat2 Bam files
cd ${hisat2_mapping_path}
hisat2_bam_files=$(ls *.bam)
hisat2_sample_Labels=$(ls *.bam | xargs -n1 sh -c 'echo ${0//.accepted_hits.bam}')
hisat2_bam_files_number=`ls *.bam | wc -l`

multiBamSummary bins \
    --bamfiles $hisat2_bam_files \
    --minMappingQuality 30 \
    --labels $hisat2_sample_Labels \
    --numberOfProcessor $threads \
    --outFileName readCounts.npz \
    --outRawCounts readCounts.tab



if [ $hisat2_bam_files_number -ge 3 ]
then
    plot_bam_function
    plotPCA  \
    -in readCounts.npz \
    --plotTitle "PCA plot of Read Counts" \
    -o PCA_readCounts.png
fi


echo ">>>>>>>>>>> Step-5:  Generate BW files <<<<<<<<<<<"
# STAR Bam files
cd ${output_path} && mkdir -p ./4_bw_STAR && cd ./4_bw_STAR && STAR_bw_path=`pwd`
cd ${star_mapping_path}

# bamCoverage did not seem to run fast, though I specified the -p $threads, so try xargs
export global_STAR_bw_path=$STAR_bw_path
ls *sortedByCoord.out.bam | xargs -P $(($threads/6)) -n 1 \
sh -c 'bamCoverage -b $0 -o ${global_STAR_bw_path}/${0%Aligned.sortedByCoord.out.bam}STAR.RPKM.bw -of bigwig -p 6 -bs 10 --normalizeUsing RPKM'
export -n global_STAR_bw_path

# Hisat2 Bam files
cd ${output_path} && mkdir -p ./4_bw_Hisat2 && cd ./4_bw_Hisat2 && hisat2_bw_path=`pwd`
cd ${hisat2_mapping_path}

export global_Hisat2_bw_path=$hisat2_bw_path
ls *.accepted_hits.bam | xargs -P $(($threads/6)) -n 1 \
sh -c 'bamCoverage -b $0 -o ${global_Hisat2_bw_path}/${0%accepted_hits.bam}Hisat2.RPKM.bw -of bigwig -p 6 -bs 10 --normalizeUsing RPKM'
export -n global_Hisat2_bw_path


echo ">>>>>>>>>> Step-6A: RNAseq QC with RSeQC <<<<<<<<<<"
echo ">>>>>>>>>> Step-6A.1: STAR align resulsts <<<<<<<<<<"
cd ${output_path} && mkdir -p ./5_RSeQC_report && cd ./5_RSeQC_report && rseqc_path=`pwd`
cd ${star_mapping_path}

# not tested code bellow!!!
# run in background and only one cpu is needed.
# running time = samples numbers * 20 mins
# nohup geneBody_coverage.py -r $genes_bed -i $hisat2_mapping_path -o $rseqc_path/GeneBodyCoverage &


export global_rseqc_path=${rseqc_path} && export global_genes_bed=${genes_bed}
ls *sortedByCoord.out.bam | xargs -P $(($threads/4)) -n 1 \
sh -c 'bam_stat.py -i $0 > ${global_rseqc_path}/${0%Aligned.sortedByCoord.out.bam}RSeQC.bam_stat.txt & \
read_distribution.py -i ${0} -r $global_genes_bed >  ${global_rseqc_path}/${0%Aligned.sortedByCoord.out.bam}STAR_RSeQC.reads_distribution.txt & \
junction_saturation.py -i $0 -r $global_genes_bed -o ${global_rseqc_path}/${0%Aligned.sortedByCoord.out.bam}STAR_RSeQC.junction_saturation & \
junction_annotation.py -i $0 -r $global_genes_bed -o ${global_rseqc_path}/${0%Aligned.sortedByCoord.out.bam}STAR_RSeQC.junction_annotation'


echo ">>>>>>>>>> Step-6A.2: Hisat2 align resulsts <<<<<<<<<<"
cd ${hisat2_mapping_path}
# read_distribution.py: no more than 5 mins for each sample 
# junction_annotation.py: no more than 3 mins for each sample
ls *.accepted_hits.bam | xargs -P $(($threads/4)) -n 1 \
sh -c  'bam_stat.py -i $0 > ${global_rseqc_path}/${0%.accepted_hits.bam}_Hisat2_RSeQC.bam_stat.txt & \
read_distribution.py -i ${0} -r $global_genes_bed >  ${global_rseqc_path}/${0%.accepted_hits.bam}_Hisat2_RSeQC.reads_distribution.txt & \
junction_saturation.py -i $0 -r $global_genes_bed -o ${global_rseqc_path}/${0%.accepted_hits.bam}_Hisat2_RSeQC.junction_saturation & \
junction_annotation.py -i $0 -r $global_genes_bed -o ${global_rseqc_path}/${0%.accepted_hits.bam}_Hisat2_RSeQC.junction_annotation'



echo ">>>>>>>>>> Step-6B: RNAseq QC with Qualimap <<<<<<<<<<"
echo ">>>>>>>>>> Step-6B.1: STAR align resulsts <<<<<<<<<<"
cd ${output_path} && mkdir -p ./5_Qualimap_report && cd ./5_Qualimap_report && qualimap_path=`pwd`
cd ${star_mapping_path}
export global_qualimap_path=${qualimap_path} && export global_genes_gtf=${genes_gtf}
ls *sortedByCoord.out.bam | xargs -P $(($threads/6)) -n 1 \
sh -c  'qualimap bamqc -bam $0 -outdir ${global_qualimap_path}/${0%.Aligned.sortedByCoord.out.bam}_STAR -nt 6 && \
qualimap rnaseq -bam $0 -gtf $global_genes_gtf -outdir ${global_qualimap_path}/${0%.Aligned.sortedByCoord.out.bam}_STAR'

echo ">>>>>>>>>> Step-6B.2: Hisat2 align resulsts <<<<<<<<<<"
cd ${hisat2_mapping_path}
ls *.accepted_hits.bam | xargs -P $(($threads/6)) -n 1 \
sh -c  'qualimap bamqc -bam $0 -outdir ${global_qualimap_path}/${0%.accepted_hits.bam}_Hisat2 -nt 6 && \
qualimap rnaseq -bam $0 -gtf $global_genes_gtf -outdir ${global_qualimap_path}/${0%.accepted_hits.bam}_Hisat2'

export -n global_rseqc_path & export -n global_genes_bed & export -n global_qualimap_path

echo ">>>>>>>>>> Step-7: Summarizing Gene/Isform Expression Values <<<<<<<<<<<"
echo ">>>>>>>>>> Step-7.1: Summarizing Gene Counts with featureCounts/SUBREAD <<<<<<<<<<<"
cd ${output_path} && mkdir -p ./6_Counts_featureCounts && cd ./6_Counts_featureCounts && featurecounts_path=`pwd`
cd ${star_mapping_path}

dirlist=$(ls -t ./*sortedByCoord.out.bam | tr '\n' ' ')

if [[ "$paired" = true ]]
then
    featureCounts \
        -a $genes_gtf \
        -F GTF \
        -o ${featurecounts_path}/final_counts.txt \
        -g 'gene_name' \
        -T 6 \
        -M \
        -p \
        $dirlist
else
    featureCounts \
        -a $genes_gtf \
        -F GTF \
        -o ${featurecounts_path}/final_counts.txt \
        -g 'gene_name' \
        -T 6 \
        -M \
        $dirlist
fi
# Never recommend the use of featureCounts to produce transcript-level counts! so use RSEM here

echo ">>>>>>>>>> Step-7.2: Summarizing Gene/Transcrpt/Exon Counts with Stringtie <<<<<<<<<<<"
cd ${output_path} && mkdir -p ./6_Counts_StringTie && cd ./6_Counts_StringTie && stringtie_path=`pwd`
mkdir -p ./Make-Reference && cd ./Make-Reference && stringtie_ref_path=`pwd`

echo "StringTie: Prepare the transcripts reference..."
cd ${hisat2_mapping_path}

# 组装转录本
export global_genes_gtf=${genes_gtf} && export global_stringtie_ref_path=${stringtie_ref_path}
ls *.accepted_hits.bam | xargs -P $(($threads/2)) -n 1 \
sh -c 'stringtie $0 -G ${global_genes_gtf} -o ${global_stringtie_ref_path}/${0%.accepted_hits.bam}.transcripts.gtf'

# 合并样本的gtf文件，得到非冗余转录本集
cd ${stringtie_ref_path}
stringtie --merge \
    -G ${genes_gtf} \
    -F 0.1 \
    -T 0.1 \
    -i \
    -o StringTie_merged.gtf *.gtf

# 将各样品的注释文件与参考基因组+参考注释比较，得到新转录本，其实基本没有啥变化！
gffcompare -r ${genes_gtf} StringTie_merged.gtf -o gffcmp

echo "StringTie: Quantify transcripts..."
cd $hisat2_mapping_path
export global_stringtie_path=$stringtie_path && export global_stringtie_ref_path=$stringtie_ref_path
ls *.accepted_hits.bam | xargs -P $(($threads/2)) -n 1 \
sh -c 'mkdir -p ${global_stringtie_path}/results/${0%.accepted_hits.bam} && \
    stringtie $0 -eB -G ${global_stringtie_ref_path}/gffcmp.annotated.gtf -o  ${global_stringtie_path}/results/${0%.accepted_hits.bam}/stringTie_asm.gtf -p 2 \
    -A  ${global_stringtie_path}/results/${0%.accepted_hits.bam}/gene_abundence.tab -C  ${global_stringtie_path}/results/${0%.accepted_hits.bam}/known.cov_refs.gtf'
export -n global_stringtie_path=$stringtie_path && export -n global_stringtie_ref_path=$stringtie_ref_path

# 输出counts矩阵
cd $stringtie_path/results
ls | while read id; do echo -e ${id}'\t'./${id}/stringTie_asm.gtf >> sample.list; done
prepDE.py -i sample.list

echo "StringTie: Quantify transcripts based on the known transcripts..."
cd ${output_path} && mkdir -p ./6_Counts_StringTie_OnlyKnownTranscripts && cd ./6_Counts_StringTie_OnlyKnownTranscripts && StringTie_Known_path=`pwd`
cd $hisat2_mapping_path

export global_StringTie_Known_path=$StringTie_Known_path
ls *.accepted_hits.bam | xargs -P $(($threads/2)) -n 1 \
sh -c 'mkdir -p ${global_StringTie_Known_path}/${0%.accepted_hits.bam} && \
    stringtie $0 -eB -G ${global_genes_gtf} -o  ${global_StringTie_Known_path}/${0%.accepted_hits.bam}/stringTie_asm.gtf -p 2 \
    -A  ${global_StringTie_Known_path}/${0%.accepted_hits.bam}/gene_abundence.tab -C  ${global_StringTie_Known_path}/${0%.accepted_hits.bam}/known.cov_refs.gtf'
export -n global_StringTie_Known_path

# 输出counts矩阵
cd $StringTie_Known_path
ls | while read id; do echo -e ${id}'\t'./${id}/stringTie_asm.gtf >> sample.list; done
prepDE.py -i sample.list

# 输出两种方式得到基因水平的TPM TPKM矩阵
cd ${output_path}
if [[ -e ${merge_stringtie_results_code} ]]; then
    Rscript ${merge_stringtie_results_code} -p ${stringtie_path}/results
    Rscript ${merge_stringtie_results_code} -p ${StringTie_Known_path}
fi


echo ">>>>>>>>>> Step-7.3: Summarizing Gene/Transcrpt Counts with RSEM <<<<<<<<<<<"
# sampleName.genes.results: gene_id transcript_id(s)        length  effective_length        expected_count  TPM     FPKM
# sampleName.isoforms.results: transcript_id   gene_id length  effective_length        expected_count  TPM     FPKM    IsoPct
# sampleName.transcript.bam
cd ${output_path} && mkdir -p ./6_Counts_RSEM && cd ./6_Counts_RSEM && RSEM_Merge=`pwd`
mkdir -p ./Single && cd ./Single && RSEM_path=`pwd`
cd $star_mapping_path

# length < seed length (= 25) will throw an error from rsem-run-em
# rsem: 18 CPUs at most
export global_RSEM_path=$RSEM_path && export global_rsem_index=$rsem_index
if [[ "$paired" = true ]]
then
    ls *Aligned.toTranscriptome.out.bam | xargs -P $(($threads/8)) -n 1 \
    sh -c 'rsem-calculate-expression --time -p 8 --bam --paired-end $0 $global_rsem_index ${global_RSEM_path}/${0%_Aligned.toTranscriptome.out.bam}'
else
    ls *Aligned.toTranscriptome.out.bam | xargs -P $(($threads/8)) -n 1 \
    sh -c 'rsem-calculate-expression --time -p 8 --bam $0 $global_rsem_index ${global_RSEM_path}/${0%_Aligned.toTranscriptome.out.bam}'
fi

export -n global_RSEM_path && export -n global_rsem_index

# 合并多样本的RSEM输出结果
cd ${output_path}
if [[ -e ${merge_RSEM_results_code} ]]; then
    Rscript ${merge_RSEM_results_code} -g ${genes_gtf} -p ${RSEM_path} -o ${RSEM_Merge}
fi

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
