#!/bin/bash
# Using getopt

source ~/.bashrc
conda activate /home/xilab/software/miniconda-envs/msPIPE

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
# Analysis Roadmap
# bismark + ?
#
# install softwares
# conda activate py37
# conda install -c bioconda fastqc
# conda install -c bioconda trim-galore
# conda install -c bioconda bismark
# conda install -c bioconda samtools
# conda install -c bioconda multiqc
######################################################################################################################
# fastq dir should be a link!
fastq=fastq
output=results
threads=24
test_reads=400000
bismark_parallel=8 # n parallel jobs for bismark,use about  8 * 5 =  about 40 threads？ and 40G RAM?
bismark_methylation_parallel=12 # n parallel jobs for bismark_methylation_extractor,use about 12*3 threads!
reads_min_length=20
paired=true
bismark_index=/home/xilab/software/msPIPE/msPIPE/reference/galGal6
# test模式，也可以帮助检查是否有接头序列、引物序列啥的
test=false

############################### 果蝇案例
################## UCSC - galGal6
# bismark
# 说明书： https://rawgit.com/FelixKrueger/Bismark/master/Docs/Bismark_User_Guide.html
# 该软件最近一直有更新！并且会有一些漂亮的统计图表（软件官网）！ 推荐！
# https://www.bioinformatics.babraham.ac.uk/projects/bismark/
# https://www.jianshu.com/p/418e4e29dd79
# 使用bismark_genome_preparation来对基因组进行处理，输入是给定的目录，目录里面是.fa或者.fasta文件。
# 下载的fasta文件位于/home/xilab/zhangyc/DataAnalysis/chicken/chicken-reference中。
# bismark_genome_preparation /home/xilab/zhangyc/DataAnalysis/chicken/chicken-reference
# bismark会生成两个目录，一个是C->T转换的基因组，一个是G->A转换的基因组共两套基因组，之后可以用bowtie2-build来建立索引。
# Bisulfite_Genome
# ├── CT_conversion
# │   ├── BS_CT.1.bt2
# │   ├── BS_CT.2.bt2
# │   ├── BS_CT.3.bt2
# │   ├── BS_CT.4.bt2
# │   ├── BS_CT.rev.1.bt2
# │   ├── BS_CT.rev.2.bt2
# │   └── genome_mfa.CT_conversion.fa
# └── GA_conversion
#     ├── BS_GA.1.bt2
#     ├── BS_GA.2.bt2
#     ├── BS_GA.3.bt2
#     ├── BS_GA.4.bt2
#     ├── BS_GA.rev.1.bt2
#     ├── BS_GA.rev.2.bt2
#     └── genome_mfa.GA_conversion.fa
# 
#
# Augument Parsing
print_usage_and_exit(){
    echo -e "Usage: $0 
                [-p <paired data> = $paired default]
                [-j <number of threads> = $threads default]
                [-i <directory of bismark index> = $bismark_index default]
                [-n <number of parallels for bismark> = $bismark_parallel default]
                [-m <Trimgalore minimum reads length> = $length default] 
                [-f <fastq file directory> = $fastq default]
                [-o <results directory> = $output default]
                [-t <test mode> = Default false - using all reads]

Examples:
    Attention:
        Rename the fastq dirs with their sample names before running this script.
        Remove the short and not used sequence (not show in gtf file)in fasta file will improve the run speed.
        10x galGal6 数据6个，跑了72h。
    Activate run env:
        conda activate /home/xilab/software/miniconda-envs/msPIPE
    Chicken galGal6:
        sh WGBS-Bismark.sh \\
        -p true
    msPIPE pipeline:
        ~/software/msPIPE/msPIPE/msPIPE.py -p params_galGal6.conf -c 12 --skip_calling --calling_data ./results/5_msPIPE/methylCALL -o ./results/5_msPIPE
        ~/software/msPIPE/msPIPE/msPIPE.py -p params_galGal6.conf -c 12 --skip_calling --calling_data ./results/5_msPIPE/methylCALL -o ./results/5_msPIPE --bsmooth --skip_GMA
"
    exit 1
}

while getopts ":t:p:j:m:i:f:o:h:" opt; do
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
            bismark_index="$OPTARG"
            echo "-i <Bismark index> = $bismark_index"
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
if [ -z "$bismark_index" ]
then
    echo "Error: missing required argument(s)"
    print_usage_and_exit
fi

if [ -d ${bismark_index%\/*} ]
then
    echo "Bismark index dir found! "
else
    echo "Bismark index dir not exist! Please input correct Bismark index dir!"
    print_usage_and_exit
fi

if [ -d ${fastq%\/*} ]
then
    echo "fastq file directory found! "
else
    echo "fastq file directory not exist! Please input correct fastq file directory!"
    print_usage_and_exit
fi

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
    test_lines=${test_reads}
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


echo ">>>>>>>>>> Step-3:  mapping and call methylation with bismark <<<<<<<<<<"
echo ">>>>>>>>>> Step-3.1:  mapping reads with bismark <<<<<<<<<<"
# bismark [options] --genome <genome_folder> {-1 <mates1> -2 <mates2> | <singles>}
# for single data: bismark --genome /data/genomes/homo_sapiens/GRCh37/ test_dataset.fastq
# 左边的C->T，是原始链的转化方式，然后跟两种参考基因组进行比对，得到OT和OB的结果；
# 右边的G->A，是扩增链的转化方式，然后跟两种参考基因组进行比对，得到CTOT和CTOB的结果；
# (一种是先片段化，然后末修，加上特殊处理的接头（甲基化的）后，进行CT转化，之后进行扩增，如下图)上面示意图提到建库方式是有方向性的,链特异性文库，因此bismark只需要进行两次比对即可(使用bismark默认参数比对)，而不需要进行4次的比对，如果文库是非链特异性文库，
# 需要加上--non_directional的参数，bismark会进行4次比对，选择最优的结果进行输出。
# 针对先转化后加接头的文库（PBAT），常见如单细胞甲基化文库，需要添加--pbat参数进行比对。
# 针对链特异性文库，若进行单端比对read1使用bismark单端参数即可，read2需要添加--pbat参数；而非链特异性文库--pbat参数的使用则相反。
cd ${output_path} && mkdir -p ./3_aligned_BISMARK && cd ./3_aligned_BISMARK && bismark_mapping_path=`pwd`
cd ${trim_path}
# 一定要问清楚文库类型！
if [[ "$paired" = true ]]
then
    ls *_val_1.fq.gz | while read id;
    do
        bismark \
            --bowtie2 \
            --parallel ${bismark_parallel} \
            --score_min L,0,-0.6 \
            --genome ${bismark_index} \
            -N 0 \
            -L 20 \
            -1 ${id} \
            -2 ${id/_val_1.fq.gz/_val_2.fq.gz} \
            -o ${bismark_mapping_path} \
            2> ${bismark_mapping_path}/${id/_val_1.fq.gz/.Bismark.mapping.log}
    done
else
    ls *trimmed.fq.gz | while read id
    do
        bismark \
            --bowtie2 \
            --score_min L,0,-0.6 \
            --parallel  ${bismark_parallel} \
            --genome ${bismark_index} \
            -N 0 \
            -L 20 \
            -o ${bismark_mapping_path} \
            ${id} \
            2> ${bismark_mapping_path}/${id/_val_1.fq.gz/.Bismark.mapping.log}
    done
fi
# --score_min L,0,-0.2(Default). This is quite a stringent setting, but you shoule not expect a big difference!
# 时间很长，至少一小时每个样本？
# 说是不建议使用-p大参数！建议用--parallel，分解成几份！参考https://github.com/FelixKrueger/Bismark/issues/106，
# --parallel为4，会同时运行16个bowtie2-align-s和8个samtools和8个perl，所以不要设置的太大。 40M reads大概1h15min
# 该步骤会生成文件：test_dataset_bismark_bt2.bam (包含比对结果和甲基化水平) bam文件：可以利用SeqMonk进行可视化；也可以用于继续的去除dup等
# test_dataset_bismark_SE_report.txt (包含比对和甲基化水平的统计结果)
# 以及未比对到参考基因组上的reads文件（--un参数）
# 若测序数据为双端数据文件中会有PE或pe标签，若为单端数据则会有SE或se字符加以区别。
# 报告展示了唯一比对上的read对数目等信息（Final Alignment report部分）以及有reads支持发生甲基化的C的数目统计（Final Cytosine Methylation Report部分），注意这里并非真正意义上的甲基化C，仍需要下游的分析。
# 运行 deduplicate_bismark
# 对于单选测序序列，去dup的时候，会根据染色体、起始坐标和链的方向来确定是否是duplication；
# 对于双端测序序列，Read1的染色体、起始坐标，Read2的起始坐标会用于去除duplication。需要注意的是，去dup的时候需要注意read1和read2是相邻的，而不是按照坐标排序的。



cd ${bismark_mapping_path} 
# 尝试在这一步修改bam和report的文件名！  未测试！！
if [[ "$paired" = true ]]
then
    ls *_pe.bam | xargs -P $threads -n 1 sh -c 'mv ${0} ${0/_val_1_bismark_bt2_pe.bam/_bismark_bt2_pe.bam}'
else
    ls *_se.bam | xargs -P $threads -n 1 sh -c 'mv ${0} ${0/_val_1_bismark_bt2_se.bam/_bismark_bt2_se.bam}'
fi

if [[ "$paired" = true ]]
then
    ls *bismark_bt2_PE_report.txt | xargs -P $threads -n 1 sh -c 'mv ${0} ${0/_val_1_bismark_bt2_PE_report.txt/_bismark_bt2_PE_report.txt}'
else
    ls *bismark_bt2_SE_report.txt | xargs -P $threads -n 1 sh -c 'mv ${0} ${0/_val_1_bismark_bt2_SE_report.txt/_bismark_bt2_SE_report.txt}'
fi



echo ">>>>>>>>>> Step-3.2:  remove duplicated reads by  'deduplicate_bismark' <<<<<<<<<<"

# if [[ "$paired" = true ]]
# then
#     ls *_pe.bam | while read id
#     do 
#         deduplicate_bismark ${id}  \
#             --paired \
#             --bam   \
#             --output_dir ./ \
#             2> ./${id/_val_1_bismark_bt2_pe.bam/.Bismark.deduplicate.log}
#     done
# else
#     ls *_se.bam | while read id
#     do 
#         deduplicate_bismark ${id}  \
#             --single \
#             --bam   \
#             --output_dir ./ \
#             2> ./${id/_trimmed_bismark_bt2_pe.bam/.Bismark.deduplicate.log}
#     done
# fi
if [[ "$paired" = true ]]
then
    ls *_pe.bam | xargs -P $threads -n 1 sh -c 'deduplicate_bismark ${0}  \
            --paired \
            --bam   \
            --output_dir ./ \
            2> ./${0/_bismark_bt2_pe.bam/.Bismark.deduplicate.log}'
else
    ls *_se.bam | xargs -P $threads -n 1 sh -c 'deduplicate_bismark ${0}  \
            --single \
            --bam   \
            --output_dir ./ \
            2> ./${0/_bismark_bt2_pe.bam/.Bismark.deduplicate.log}'
fi
# 这一步有些慢，大概30min每个样本
# 对 Bismark 的比对结果去重，去除以相同方向比对到相同位置的 reads。建议应用于 WGBS 数据，不适用与 RRBS ，amplicon，以及 target 富集文库的数据。

echo ">>>>>>>>>> Step-3.3:  methylation calling by bismark_methylation_extractor <<<<<<<<<<"
# 提取甲基化位点, 运行 bismark_methylation_extractor
# 使用bismark_methylation_extractor来提取每个C位点的甲基化的情况，结果文件可以导入到SeqMonk中。
# 可以使用--bedGraph，来产生bedgraph和coverage的文件；--count可以产生count文件，输出每个CpG位点的深度文件。
mkdir -p ./methylation && cd ./methylation && methylation=`pwd`
cd ${bismark_mapping_path} 
if [[ "$paired" = true ]]
then    
    ls *deduplicated.bam | while read id
    do 
        bismark_methylation_extractor \
            --paired-end  \
            --comprehensive \
            --no_overlap \
            --bedGraph \
            --counts \
            --gzip \
            --ample_memory \
            --CX_context \
            --report \
            --parallel ${bismark_methylation_parallel} \
            --cytosine_report  \
            --buffer_size 50% \
            --genome_folder ${bismark_index} \
            -o ${methylation} \
            ${id}
    done
else
    ls *deduplicated.bam | while read id
    do 
        bismark_methylation_extractor \
            --single-end  \
            --comprehensive \
            --bedGraph \
            --counts \
            --gzip \
            --ample_memory \
            --report \
            --CX_context \
            --parallel ${bismark_methylation_parallel} \
            --cytosine_report  \
            --buffer_size 50% \
            --genome_folder ${bismark_index} \
            -o ${methylation} \
            ${id}
    done
fi

# 很慢 1h每个样本
# --pair-end         # 指定双端数据
# --comprehensive    # 输出CHG CHH CpG的甲基化信息
# --no-overlap       # 去除reads重叠区域的bias
# --bedGraph         # 输出bedGraph文件
# --counts           # 每个C上甲基化reads和非甲基化reads的数目
# --report           # 一个甲基化summay
# --cytosine_report  # 输出全基因组所有CpG
# --genome_folder    # <path_to_reference_genome>
# input.bam          # 输入文件
# A typical command to extract context-dependent (CpG/CHG/CHH) methylation could look like this: 
# 生成三个甲基化文件：
# CpG_context_test_dataset_bismark_bt2.txt.gz
# CHG_context_test_dataset_bismark_bt2.txt.gz
# CHH_context_test_dataset_bismark_bt2.txt.gz
# 以及一个 bedGraph 和一个 Bismark coverag文件。
# bismark对不同状态的胞嘧啶进行了约定并记录到了BAM文件和中间文件中：
# X 代表CHG中甲基化的C
# x  代笔CHG中非甲基化的C
# H 代表CHH中甲基化的C
# h  代表CHH中非甲基化的C
# Z  代表CpG中甲基化的C
# z  代表CpG中非甲基化的C
# U 代表其他情况的甲基化C(CN或者CHN)
# u  代表其他情况的非甲基化C (CN或者CHN)
# bismark_methylation_extractor的一个非常重要的输出文件为全基因组胞嘧啶报告文件（*CX_report.txt），该文件是一个非常核心的记录胞嘧啶深度信息的文件，
# 可进行甲基化水平的计算、区域甲基化信号文件的转化、差异甲基化水平的计算等。文件示例如下：
# chr1    14716   +       1       6       CG      CGT
# chr1    14717   -       10      1       CG      CGA
# chr1    14741   +       2       1       CG      CGG
# chr1    14742   -       0       1       CG      CGC
# 第7列 : C背景序列
# 基于该文件可进行单个胞嘧啶位点的甲基化水平（C/C+T）计算，比如第一行：甲基化是水平为1/(1+6)；针对指定区域的甲基化水平计算则取区域内的胞嘧啶的总C/(总C+总T)即可。
# 四、bismark的优缺点
# 4.1、优点
# （1）比对和mC Calling比较准确，受众广泛，使用率很高
# （2）软件维护性好，更新快，bug修复及时
# （3）软件成体系，集成了比对、去重、mC Calling分析
# 4.2、缺点
# （1）比对和mC Calling速度慢，真的很慢，可以通过拆分输入数据为小份作为输出和增大计算资源解决。
# （2）bismak相匹配的直接进行下游甲基化水平、图形可视化等分析的套件较少，需要自行计算或进行格式转化作为第三方软件的输入。
# （3）说明文档不够系统和详细，在”外行“看来理解从建库到mC Calling有一定的挑战。
# bismark的结果只是研究甲基化信息的第一步也是关键的一步。未来，期望bismark在速度上，和在对甲基化水平层面分析的支持上多一下改进
# --ample_memory           Using this option will not sort chromosomal positions using the UNIX 'sort' command, but will
#                          instead use two arrays to sort methylated and unmethylated calls. This may result in a faster
#                          sorting process of very large files, but this comes at the cost of a larger memory footprint
#                          (two arrays of the length of the largest human chromosome 1 (~250M bp) consume around 16GB
#                          of RAM). Due to overheads in creating and looping through these arrays it seems that it will
#                          actually be *slower* for small files (few million alignments), and we are currently testing at
#                          which point it is advisable to use this option. Note that --ample_memory is not compatible
#                          with options '--scaffolds/--gazillion' (as it requires pre-sorted files to begin with).


# Bismark Nucleotide Coverage report 
cd ${bismark_mapping_path} 
ls *deduplicated.bam | while read id
do 
    bam2nuc --genome_folder ${bismark_index} ${id} 
done


echo ">>>>>>>>>> Step-3.4:  generate bismark analysis report <<<<<<<<<<"
# 运行 bismark2report
cd ${bismark_mapping_path} 
bismark2report


# 总结 report 结果，可以在运行完 bam2nuc 之后
bams=`ls */*_bismark_bt2_pe.bam|tr '\n' ' '`
bismark2summary $bams
# This command scans the current working directory for different Bismark alignment, 
# deduplication and methylation extraction (splitting) reports to produce a graphical summary HTML report, 
# as well as a data table, for all files in a directory. Here is a sample Bismark Summary Report.

# msPIPE
# https://github.com/jkimlab/msPIPE
# 貌似可以做bismark的下游分析！

 # R 包 methylKit  R v4，本地电脑上做吧！
 # 适用于Bismark aligner结果
 # https://www.bioconductor.org/packages/release/bioc/vignettes/methylKit/inst/doc/methylKit.html

echo ">>>>>>>>>> Step-4: Generating analysis report with multiQC <<<<<<<<<<<"
cd ${output_path} && mkdir -p ./4_multiQC && cd ./4_multiQC && multiqc_path=`pwd`
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

echo ">>>>>>>>>> Step-5: Delete unused files <<<<<<<<<<<"
cd ${output_path}
cd $trim_path
echo "Delete the trimmed fastqs..."
rm *.fq.gz

if [[ "$test" = true ]]
then
    # rm -rf $fastq
    rm -rf ${wroking_folder}/test
fi


echo ">>>>>>>>>> Step-6: Rearrange and prepare files for msPIPE pipeline <<<<<<<<<<<"
cd ${output_path} && mkdir -p ./5_msPIPE/methylCALL && cd ./5_msPIPE/methylCALL && mspipe_path=`pwd`
mv ${methylation}/* ${mspipe_path}/

cat ${output_path}/sample_fastq_meta.txt | sed '1d' | while read sample r1 r2
do 
    mkdir -p ./${sample}/methylcontext && mkdir -p ./${sample}/data
    ls -l | grep "^-" | awk '{print $9}' | grep ${sample} | while read id; do mv ${id} ./${sample}/methylcontext; done
done

ls | while read id
do
    bismark2bedGraph -o CpG.cov --dir ${id}/methylcontext ${id}/methylcontext/CpG_context_*  --buffer_size 50% --ample_memory
    bismark2bedGraph -o CHG.cov --dir ${id}/methylcontext --CX ${id}/methylcontext/CHG_context_*  --buffer_size 50% --ample_memory
    bismark2bedGraph -o CHH.cov --dir ${id}/methylcontext --CX ${id}/methylcontext/CHH_context_*  --buffer_size 50% --ample_memory
    ## split files
    gzip -dc ${id}/methylcontext/CpG.cov.gz.bismark.cov.gz 1> ${id}/methylcontext/${id}_CpG.cov.txt
    gzip -dc ${id}/methylcontext/CHG.cov.gz.bismark.cov.gz 1> ${id}/methylcontext/${id}_CHG.cov.txt
    gzip -dc ${id}/methylcontext/CHH.cov.gz.bismark.cov.gz 1> ${id}/methylcontext/${id}_CHH.cov.txt

    mkdir -p ${id}/methylcontext/CpG_chr && mkdir -p ${id}/methylcontext/CHG_chr && mkdir -p ${id}/methylcontext/CHH_chr

    ~/software/msPIPE/msPIPE/bin/script/splitF_bychr.py ${id}/methylcontext/${id}_CpG.cov.txt ${id}/methylcontext/CpG_chr/${id}
    ~/software/msPIPE/msPIPE/bin/script/splitF_bychr.py ${id}/methylcontext/${id}_CHG.cov.txt ${id}/methylcontext/CHG_chr/${id}
    ~/software/msPIPE/msPIPE/bin/script/splitF_bychr.py ${id}/methylcontext/${id}_CHH.cov.txt ${id}/methylcontext/CHH_chr/${id}
done


cd ${wroking_folder}
echo "#******************************************************************#"
echo "Congratulations! Whole processs finished! "
echo "WGBS Bismark Pipline Finished At:"
endtime=`date +'%Y-%m-%d %H:%M:%S'`
echo "${endtime}"
start_seconds=$(date --date="$starttime" +%s);
end_seconds=$(date --date="$endtime" +%s);
secondsdiff=$((end_seconds-start_seconds))
echo "程序总运行时间： "$((secondsdiff/3600))" h "$((secondsdiff%3600/60))" m"
trap : 0
echo >&2 '#******************************************************************#'
