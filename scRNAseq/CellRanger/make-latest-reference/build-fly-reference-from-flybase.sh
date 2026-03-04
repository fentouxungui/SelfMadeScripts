#!/bin/bash
# refer to: https://www.10xgenomics.com/cn/support/software/cell-ranger/tutorials/cr-tutorial-mr
project="flybase-dmel-r6.55" # 目录名
ref_dir="reference" # reference存储目录
fa_url="https://ftp.flybase.net/releases/current/dmel_r6.55/fasta/dmel-all-chromosome-r6.55.fasta.gz"
gtf_url="https://ftp.flybase.net/releases/current/dmel_r6.55/gtf/dmel-all-r6.55.gtf.gz"
cellranger="/home/xilab/software/cellranger7/cellranger-7.2.0/bin/cellranger"
genome="dmel-r6_55" #cellranger 基因名

working_dir=`pwd`
if [  -d  ${project} ]; then
	echo "项目目录已经存在！"
else
	mkdir -p ${project}
fi

cd ${working_dir}/${project}&& mkdir -p ${ref_dir} && cd ${ref_dir}

# download fasta file
wget ${fa_url}
# download gtf file
wget ${gtf_url}
gunzip *.gz

fa_file=$(basename ${fa_url})
gtf_file=$(basename ${gtf_url})

# some errors in gtf file from flybase
# $ cat dmel-all-r6.55.gtf | cut -f 7 | sort | uniq -c
#  279249 -
#      9 .
#  269873 +

awk 'BEGIN{FS=OFS="\t"} {if ($7 == ".") $7="-";}1' ${gtf_file%.gz} > ${gtf_file%gtf.gz}modified.gtf



# filter the gtf, remove entries for non-polyA transcripts, for These entries can cause reads to be flagged as mapped to multiple genes (multi-mapped) 
# because of the overlapping annotations. In the case where reads are flagged as multi-mapped, they are not counted
# 注意，没有去除非蛋白编码基因！
# ${cellranger} mkgtf  ${gtf_name} ${gtf_name%gtf}filtered.gtf --attribute=gene_biotype:protein_coding

cd ${working_dir}/${project}
${cellranger} mkref --genome=${genome} --fasta=${ref_dir}/${fa_file%.gz} --genes=${ref_dir}/${gtf_file%gtf.gz}modified.gtf