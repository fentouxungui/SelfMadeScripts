#!/bin/bash
# Using getopt

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
# Analysis Roadmap
# 1. STAR + FeatureCounts > Gene level expression matrix - Counts > DESeq2
# 2. STAR + RSEM > DESeq2
# 3. Hisat2 + StringTie > Gene/Transcript(novel and known)/Exon level expression matrix - Counts and TPM > Ballgown 

# Softwares needed:
# 1. FASTQ QC: fastqc;
# 2. Trim Adapter: Trimgalore;
# 3. Align reads: STAR; Hisat2; Samtools;
# 4. BW files: deepTools: bamCoverage, plotCorrelation, multiBamSummary; 
# 5. RNAseq QC: RSeQC: bam_stat.py, read_distribution.py, geneBody_coverage.py(Not Used);
# 6. Summarizing Gene/Transcrpt/Exon Counts: Subread: featureCounts; GffCompare + StringTie;
# 7. Analysis Report: multiqc

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
hisat2_index=/home/xilab/reference/hisat2_index/dm6-UCSC/dm6/genome
star_index=/data0/reference/STAR_index/Drosophlia/BDGP6_ensembl/STARindex_ReadsLengthDefault
genes_bed=/home/xilab/reference/Genome/Drosophlia/UCSC/dm6/dm6.refGene.bed
genes_gtf=/home/xilab/reference/Genome/Drosophlia/UCSC/dm6/dm6.refGene.gtf
rsem_index=/data0/reference/RSEM_index/BDGP6/rsem-BDGP6
# testÄ£Ê½£¬Ò²¿ÉÒÔ°ïÖú¼ì²éÊÇ·ñÓÐ½ÓÍ·ÐòÁÐ¡¢ÒýÎïÐòÁÐÉ¶µÄ
test=false

# rsem-prepare-reference \
# --gtf /data0/reference/STAR_index/Mouse/m23.GRCm38.p6.gencode/gencode.vM23.primary_assembly.annotation.gtf \
# /data0/reference/STAR_index/Mouse/m23.GRCm38.p6.gencode/GRCm38.primary_assembly.genome.fa \
#  rsem-Mouse-genecode -p 8

# hisat2_index ¿ÉÒÔ´Óhisat2g¹ÙÍøÏÂÔØ
# http://daehwankimlab.github.io/hisat2/download/

# STAR \
# --runThreadN 18 \
# --runMode genomeGenerate \
# --genomeDir STARindex_149bp \
# --genomeFastaFiles ./GRCm38.primary_assembly.genome.fa \
# --sjdbGTFfile ./gencode.vM23.primary_assembly.annotation.gtf \
# --sjdbOverhang 149
# genomeSAindexNbases 12 ÎªSTARÍÆ¼öµÄÕë¶ÔÓÚ¹ûÓ¬»ùÒò×é´óÐ¡µÄ²ÎÊý
# sjdbOverhang: ²âÐòµÄreads³¤¶È-1£¬ Ä¬ÈÏÊÇ100bp£¬Ã²ËÆÄ¬ÈÏµÄÕâ¸öÐ§¹ûÒ²²»´í£¡


# Example use for fly data
# sh RNAseq-STAR-FeatureCounts+Hisat2-Stringtie.sh

# Example use for mouse data
# sh RNAseq-STAR-FeatureCounts+Hisat2-Stringtie.sh \
# -i /data0/reference/STAR_index/Mouse/m23.GRCm38.p6.gencode/MouseGenome149bp \
# -a /data0/reference/hisat2_index/mm10-UCSC \
# -r /data0/reference/RSEM_index/mm10/rsem-Mouse-genecode \
# -g /data0/reference/STAR_index/Mouse/m23.GRCm38.p6.gencode/gencode.vM23.primary_assembly.annotation.gtf \
# -b /data0/reference/STAR_index/Mouse/m23.GRCm38.p6.gencode/gencode.vM23.primary_assembly.annotation.bed \
# -t true \
# -p true \
# -f 2022.03.28_RNAseq_PE_Mouse-Crypt-CD24-low-Runx1-CKO-Male-Female


# Augument Parsing
print_usage_and_exit(){
    echo "Usage: $0 
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
            [-t <test mode> = Default false - using all reads]"
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
        mkdir -p ${target_fq}/$Asample
        ls ${source_fq}/${Asample}/*fastq.gz | while read id
        do
            echo ">>> Subseting fastq: $(basename ${id})..."
            zcat $id | head -n ${test_lines} > ${target_fq}/${Asample}/$(basename ${id%.gz})
            gzip ${target_fq}/${Asample}/$(basename ${id%.gz})
        done
    done
    fastq=test
fi


cd ${wroking_folder} && mkdir -p ./${output} && cd ./${output} && output_path=`pwd` # results dir
cd ${wroking_folder} && cd ${fastq} && samples_folder=`pwd` && samples=($(ls -l  | grep "^[dl]" | awk '{print $9}')) # samples dir and the sample names


## generage new meta.data file
echo ">>> Preparing the sample meta data file..."
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
echo "RNA-Seq Pipline Starts At:"
starttime=`date +'%Y-%m-%d %H:%M:%S'`
echo "³ÌÐò¿ªÊ¼ÓÚ£º "${starttime}""
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
# Hisat2 Bam files
cd ${hisat2_mapping_path}
hisat2_bam_files=$(ls *.bam)
hisat2_sample_Labels=$(ls *.bam | xargs -n1 sh -c 'echo ${0//.accepted_hits.bam}')
multiBamSummary bins \
    --bamfiles $hisat2_bam_files \
    --minMappingQuality 30 \
    --labels $hisat2_sample_Labels \
    --numberOfProcessor $threads \
    --outFileName readCounts.npz \
    --outRawCounts readCounts.tab

hisat2_bam_files_number=`ls *.bam | wc -l`
if [ $hisat2_bam_files_number -ge 3 ]
then
    plot_bam_function
fi
plotPCA  \
    -in readCounts.npz \
    --plotTitle "PCA plot of Read Counts" \
    -o PCA_readCounts.png

echo ">>>>>>>>>>> Step-5:  Generate BW files <<<<<<<<<<<"
echo ">>>>>>>>>> Step-5.2: Hisat2 align resulsts <<<<<<<<<<"
# Hisat2 Bam files
cd ${output_path} && mkdir -p ./4_bw_Hisat2 && cd ./4_bw_Hisat2 && hisat2_bw_path=`pwd`
cd ${hisat2_mapping_path}

export global_Hisat2_bw_path=$hisat2_bw_path
ls *.accepted_hits.bam | xargs -P $(($threads/6)) -n 1 \
sh -c 'bamCoverage -b $0 -o ${global_Hisat2_bw_path}/${0%accepted_hits.bam}Hisat2.RPKM.bw -of bigwig -p 6 -bs 10 --normalizeUsing RPKM'
export -n global_Hisat2_bw_path




echo ">>>>>>>>>> Step-6A: RNAseq QC with RSeQC <<<<<<<<<<"
cd ${output_path} && mkdir -p ./5_RSeQC_report && cd ./5_RSeQC_report && rseqc_path=`pwd`
export global_rseqc_path=${rseqc_path} && export global_genes_bed=${genes_bed}

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
cd ${output_path} && mkdir -p ./5_Qualimap_report && cd ./5_Qualimap_report && qualimap_path=`pwd`
export global_qualimap_path=${qualimap_path} && export global_genes_gtf=${genes_gtf}
echo ">>>>>>>>>> Step-6B.2: Hisat2 align resulsts <<<<<<<<<<"
cd ${hisat2_mapping_path}
ls *.accepted_hits.bam | xargs -P $(($threads/6)) -n 1 \
sh -c  'qualimap bamqc -bam $0 -outdir ${global_qualimap_path}/${0%.accepted_hits.bam}_Hisat2 -nt 6 && \
qualimap rnaseq -bam $0 -gtf $global_genes_gtf -outdir ${global_qualimap_path}/${0%.accepted_hits.bam}_Hisat2'

export -n global_rseqc_path & export -n global_genes_bed & export -n global_qualimap_path

echo ">>>>>>>>>> Step-7: Summarizing Gene/Isform Expression Values <<<<<<<<<<<"
echo ">>>>>>>>>> Step-7.1: Summarizing Gene Counts with featureCounts/SUBREAD <<<<<<<<<<<"
cd ${output_path} && mkdir -p ./6_Counts_featureCounts && cd ./6_Counts_featureCounts && featurecounts_path=`pwd`
cd ${hisat2_mapping_path}

dirlist=$(ls -t ./*accepted_hits.bam | tr '\n' ' ')

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
echo "³ÌÐò×ÜÔËÐÐÊ±¼ä£º "$((secondsdiff/3600))" h "$((secondsdiff%3600/60))" m"
trap : 0
echo >&2 '#******************************************************************#'
