#!/bin/bash
# refer to: https://www.10xgenomics.com/cn/support/software/cell-ranger/tutorials/cr-tutorial-mr
project="ensembl-release-110"
ref_dir="reference"
fa_name="sm.primary_assembly.fa"
gtf_name="BDGP6.46.110.gtf"
cellranger="/home/xilab/software/cellranger7/cellranger-7.2.0/bin/cellranger"
genome="BDGP6-46_110"

working_dir=`pwd`
if [  -d  ${project} ]; then
	echo "项目目录已经存在！"
else
	mkdir -p ${project}
fi

cd ${working_dir}/${project}&& mkdir -p ${ref_dir} && cd ${ref_dir}

# download ensembl sm.primaryl_assembl fasta files, code refer to: https://useast.ensembl.org/info/data/ftp/rsync.html
rsync -av rsync://ftp.ensembl.org/ensembl/pub/release-110/fasta/drosophila_melanogaster/dna/Drosophila_melanogaster.BDGP6.46.dna_sm.primary_assembly* .
rsync -av rsync://ftp.ensembl.org/ensembl/pub/release-110/fasta/drosophila_melanogaster/dna/Drosophila_melanogaster.BDGP6.46.dna_sm.nonchromosomal.fa.gz .
zcat *.gz > ${fa_name} && rm *.gz
# download gtf file
rsync -av rsync://ftp.ensembl.org/ensembl/pub/release-110/gtf/drosophila_melanogaster/Drosophila_melanogaster.BDGP6.46.110.gtf.gz .
gunzip Drosophila_melanogaster.BDGP6.46.110.gtf.gz && mv Drosophila_melanogaster.BDGP6.46.110.gtf ${gtf_name}

# filter the gtf, remove entries for non-polyA transcripts, for These entries can cause reads to be flagged as mapped to multiple genes (multi-mapped) 
# because of the overlapping annotations. In the case where reads are flagged as multi-mapped, they are not counted
${cellranger} mkgtf  ${gtf_name} ${gtf_name%gtf}filtered.gtf --attribute=gene_biotype:protein_coding


cd ${working_dir}/${project}
${cellranger} mkref --genome=${genome} --fasta=${ref_dir}/${fa_name} --genes=${ref_dir}/${gtf_name%gtf}filtered.gtf 