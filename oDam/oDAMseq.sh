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
## fastq naming: SRR9305626_R2.fastq.gz and SRR9305626_R1.fastq.gz
###############################################
# change default parameters here!
bowtie2_index=/home/xilab/reference/UCSC/dm6/Sequence/Bowtie2Index/genome
genes_bed=/data0/reference/UCSC/dm6/Annotation/Archives/archive-2015-07-24-09-25-49/Genes/genes.bed
threads=24
output=results
fastq=fastq
paired=false
# test模式，也可以帮助检查是否有接头序列、引物序列啥的，本流程不包含去接头的代码
test=false

# Augument Parsing
print_usage_and_exit(){
	echo "Usage: $0 
        [-p <paired data> = $paired default]
        [-f <path of fastq files> = $fastq default]
        [-j <number of threads> = $threads default]
        [-i <bowtie2 index> = path/to/bowtie2/index]
        [-g <genes bed file> = $genes_bed]
        [-o <output directory> = $output default]
        [-t <test mode> = Default false - using all reads]"
	exit 1
}


while getopts ":t:p:f:o:g:j:i:h:" opt; do
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
        rm -rf test
    fi
    mkdir -p test && cd test
    source_fq=${wroking_folder}/${fastq}
    target_fq=`pwd`
    test_lines=100000
    if [[ -h $source_fq ]]
    then
        samples_fq=($(ls -l $(readlink ${source_fq}) | grep "^[dl]" | awk '{print $9}'))
    else
        samples_fq=($(ls -l ${source_fq} | grep "^[dl]" | awk '{print $9}'))
    fi
    for Asample in ${samples_fq[*]};
    do
        mkdir -p ${target_fq}/$Asample
        ls ${source_fq}/${Asample}/*fastq.gz | while read id
            do
            echo Subset fastq: ${id}...
            zcat $id | head -n ${test_lines} > ${target_fq}/${Asample}/$(basename ${id%.gz})
            gzip ${target_fq}/${Asample}/$(basename ${id%.gz})
            done
    done
    fastq=test
fi


cd ${wroking_folder} && mkdir -p ./${output} && cd ./${output} && output_path=`pwd` # results dir
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

# Rename the fastqc outputs and Change the Filename in sample_fastqc.zip/sample/fastqc_data.txt
cat fastq.path | while read sample_new_name fastq_path
do
    old_name=`basename $fastq_path`
    if [[ "${old_name%.fastq.gz}" != "$sample_new_name" ]]; then
        unzip ${qc_path}/${old_name/.fastq.gz/_fastqc.zip} -d ${qc_path} && rm ${qc_path}/${old_name/.fastq.gz/_fastqc.zip}
        mv ${qc_path}/${old_name/.fastq.gz/_fastqc} ${qc_path}/${sample_new_name}_fastqc
        cd ${qc_path}/${sample_new_name}_fastqc && sed -i "s/${old_name}/${sample_new_name}.fastq.gz/g" fastqc_data.txt && cd ../
        echo ${old_name} ${sample_new_name}.fastq.gz
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
    ls *sorted.bam | xargs -P $(($threads/6)) -n 1 \
    sh -c 'bamCoverage -b $0 -o ${global_bw_path}/${0%sorted.bam}Bowtie2.RPKM.bw -of bigwig -p 6 -bs 1 --extendReads 300 --normalizeUsing RPKM'

    ls *sorted.bam | xargs -P $(($threads/6)) -n 1 \
    sh -c 'bamCoverage -b $0 -o ${global_bw_path}/${0%sorted.bam}Bowtie2.BPM.bw -of bigwig -p 6 -bs 1 --extendReads 300 --normalizeUsing BPM'
fi
export -n global_bw_path


echo ">>>>>>>>>> Step-7: Generating analysis report with multiQC <<<<<<<<<<<"
cd ${output_path} && mkdir -p ./4_QC_multiQC && cd ./4_QC_multiQC && multiqc_path=`pwd`
cd ${output_path}

multiqc -f ${output_path} --outdir ${multiqc_path} -c /home/xilab/reference/Scripts/config_example.yaml

echo ">>>>>>>>>> Step-8: Delete unused files <<<<<<<<<<<"
if [[ "$test" = true ]]
then
    # rm -rf $fastq
    rm -rf ${wroking_folder}/test
fi


echo "#******************************************************************#"
echo "Congratulations! Whole processs finished! "
echo "Dam-Seq Pipline For PE Data Finished At:"
endtime=`date +'%Y-%m-%d %H:%M:%S'`
echo "${endtime}"
start_seconds=$(date --date="$starttime" +%s);
end_seconds=$(date --date="$endtime" +%s);
secondsdiff=$((end_seconds-start_seconds))
echo "程序总运行时间： "$((secondsdiff/3600))" h "$((secondsdiff%3600/60))" m"
trap : 0
echo >&2 '#******************************************************************#'
