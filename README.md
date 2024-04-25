# Scripts

自己学习和使用的一些分析脚本。

## 1. 测序数据管理

### 按样本归类fastq文件

同一样本的所有fastq文件放到同一目录下，并且以样本名命名该目录，后续分析会继续延用此样本名。

### 自动化修改fastq文件名

``Rename-fastqs.R``可自动识别一些类型的文件名，并统一格式为：``SampleName_S1_L001_R1_001.fastq.gz`` and/or ``SampleName_S1_L001_R2_001.fastq.gz``。

``` sh
./Rename-fastqs.r -h
# 测试
./Rename-fastqs.r
# 执行
./Rename-fastqs.r -r
```

### 生成md5sum值

``` sh
# check samples and md5 files: 
sh generate-md5sum-txt.sh -c true
# first time use or just want to regenerate all
sh generate-md5sum-txt.sh -f true
# just run on newlly added files: 
sh generate-md5sum-txt.sh
```

注意： 目录不能是软连接。后期再改bug!

## 2. ATACseq/Cut&Tag

###  ATACseq.sh

**using**:

``` sh
# help
sh ATACseq.sh -h

# 果蝇
sh ATACseq.sh -p true

# 小鼠
sh ATACseq.sh \
-i  /data0/reference/bowtie2_index/mm10/mm10 \
-g /data0/reference/UCSC/mm10/mm10.knownGene.bed \
-b /data0/reference/encode-atac-seq-pipeline/mm10/mm10.blacklist.bed \
-p true

# 人
sh ATACseq.sh \
-i /data0/reference/encode-atac-seq-pipeline/hg38/bowtie2_index/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta \
-g /home/xilab/reference/encode-atac-seq-pipeline/hg38/hg38.knownGene.bed \
-b /data0/reference/encode-atac-seq-pipeline/hg38/hg38.blacklist.bed \
-p true \
-t true
```

**Parameters**:

- -p: paired data, True or False
- -f: path of fastq files
- -j: number of threads
- -m: cutadapt minimum length
- -i: bowtie2 index
- -g: genes bed file
- -b: blacklist bed file
- -o: output directory
- -t: test mode, if True, will only use the top 100000 lines in the fastq files

### Chipseeker.r

For Drosophlia!

input: MACS2 opuput (narrowpeaks and summits.bed)

output: plots(feature distribution, DistributionRelativeToTSS, GO-BP, KEGG), Annotated peaks, extracted peaks for Homer(downstream Motif analysis)

Before using, pleas change the Rscript PATH in shebang: `` which Rscript``.

``` sh
./chipseeker.r -p xxx.narrowPeak \
               -s xxx.summits.bed \
               -n prefix
```      

### ChIPseeker.Rmd

same with ``Chipseeker.r``, should runs inside of rstudio IDE.


## 3. RNAseq

more info please refer to [twbattaglia/RNAseq-workflow](https://github.com/twbattaglia/RNAseq-workflow)

### RNAseq.sh

Files located in ``RNAseq/upstream-analysis``, including two other files:``merge-RSEM-results.R`` and ``merge-stringtie-results.R``，这两个文件的路径需要在``RNAseq.sh``中指定。

**Prepare**:

``` sh
# Help
############################### 果蝇案例
################## UCSC - dm6
# 1.1 hisat2_index 可以从hisat2官网下载 http://daehwankimlab.github.io/hisat2/download/
# 构建 Hisat2 index: dm6 genome_tran(hisat2 官网未提供)
cp /home/xilab/reference/Genome/Drosophlia/UCSC/dm6/dm6.fa ./genome.fa
cp /home/xilab/reference/Genome/Drosophlia/UCSC/dm6/dm6.refGene.gtf  ./genome.gtf
hisat2_extract_splice_sites.py genome.gtf > genome.ss
hisat2_extract_exons.py genome.gtf > genome.exon
hisat2-build -p 32 --exon genome.exon --ss genome.ss genome.fa genome_tran

# 1.2 构建STAR index
STAR \
--runThreadN 24 \
--runMode genomeGenerate \
--genomeDir STARindex_100bp \
--genomeFastaFiles /home/xilab/reference/Genome/Drosophlia/UCSC/dm6/dm6.fa \
--sjdbGTFfile /home/xilab/reference/Genome/Drosophlia/UCSC/dm6/dm6.refGene.gtf \
--sjdbOverhang 100 \
--genomeSAindexNbases 12
genomeSAindexNbases 12 为STAR推荐的针对于果蝇基因组大小的参数
sjdbOverhang: 测序的reads长度-1， 默认是100bp，貌似默认的这个效果也不错！

# 1.3 构建RSEM index
rsem-prepare-reference \
--gtf ./dm6.refGene.gtf \
./dm6.fa \
dm6-UCSC -p 24

############################### 小鼠案例
############# GENCODE - GRCm39 vM31
# 1.1 hisat2_index 可以从hisat2官网下载，Hisat2官网未提供GRCm39(mm39)的index
# 构建Hisat2 index
# 参考： http://daehwankimlab.github.io/hisat2/howto/#building-indexes
cp /home/xilab/reference/Genome/Mouse/GENCODE/GRCm39-vM31/GRCm39.genome.fa ./genome.fa
cp /home/xilab/reference/Genome/Mouse/GENCODE/GRCm39-vM31/gencode.vM31.primary_assembly.annotation.gtf  ./genome.gtf
hisat2_extract_splice_sites.py genome.gtf > genome.ss
hisat2_extract_exons.py genome.gtf > genome.exon
hisat2-build -p 32 --exon genome.exon --ss genome.ss genome.fa genome_tran

# 1.2 构建STAR index
ln -s /home/xilab/reference/Genome/Mouse/GENCODE/GRCm39-vM31/gencode.vM31.primary_assembly.annotation.gtf 
ln -s /home/xilab/reference/Genome/Mouse/GENCODE/GRCm39-vM31/GRCm39.primary_assembly.genome.fa ./
STAR \
--runThreadN 24 \
--runMode genomeGenerate \
--genomeDir STARindex_100bp \
--genomeFastaFiles ./GRCm39.primary_assembly.genome.fa \
--sjdbGTFfile ./gencode.vM31.primary_assembly.annotation.gtf \
--sjdbOverhang 100
sjdbOverhang: 测序的reads长度-1， 默认是100bp，貌似默认的这个效果也不错！

# 1.3 构建RSEM index
rsem-prepare-reference \
--gtf /home/xilab/reference/Genome/Mouse/GENCODE/GRCm39-vM31/gencode.vM31.primary_assembly.annotation.gtf  \
/home/xilab/reference/Genome/Mouse/GENCODE/GRCm39-vM31/GRCm39.primary_assembly.genome.fa \
Genecode-GRCm39-vM31 -p 24

############################### 人案例
############### GENCODE - GRCh38.p13.release42 Primary
# 1.1 hisat2_index 可以从hisat2官网下载
# http://daehwankimlab.github.io/hisat2/download/
# 1.2 RSEM indx
rsem-prepare-reference \
--gtf /data0/reference/Genome/Human/GENCODE/GRCh38.p13.release42/gencode.v42.primary_assembly.annotation.gtf \
/data0/reference/Genome/Human/GENCODE/GRCh38.p13.release42/GRCh38.primary_assembly.genome.fa \
GENCODE-GRCh38-p13-release42-Primary -p 24
# 1.3 STAR index
# path：/data0/reference/Genome/Human/GENCODE/GRCh38.p13.release42
cd /data0/reference/STAR_index/Human/GENCODE-GRCh38-p13-release42-primary
ln -s /data0/reference/Genome/Human/GENCODE/GRCh38.p13.release42/gencode.v42.primary_assembly.annotation.gtf ./
ln -s /data0/reference/Genome/Human/GENCODE/GRCh38.p13.release42/GRCh38.primary_assembly.genome.fa ./
STAR \
--runThreadN 24 \
--runMode genomeGenerate \
--genomeDir STARindex_100bp \
--genomeFastaFiles ./GRCh38.primary_assembly.genome.fa \
--sjdbGTFfile ./gencode.v42.primary_assembly.annotation.gtf \
--sjdbOverhang 100
# sjdbOverhang: 测序的reads长度-1， 默认是100bp，貌似默认的这个效果也不错！
################# UCSC - T2T-CHM13-v2.0-hs1
# 2.1 Hisat2 index
cp /home/xilab/reference/Genome/Human/UCSC/T2T-CHM13-v2.0-hs1/hs1.fa ./genome.fa
cp /home/xilab/reference/Genome/Human/UCSC/T2T-CHM13-v2.0-hs1/catLiftOffGenesV1.gtf  ./genome.gtf
hisat2_extract_splice_sites.py genome.gtf > genome.ss
hisat2_extract_exons.py genome.gtf > genome.exon
hisat2-build -p 32 --exon genome.exon --ss genome.ss genome.fa genome_tran

# 2.2 RSEM indx
rsem-prepare-reference \
--gtf /home/xilab/reference/Genome/Human/UCSC/T2T-CHM13-v2.0-hs1/catLiftOffGenesV1.gtf \
/home/xilab/reference/Genome/Human/UCSC/T2T-CHM13-v2.0-hs1/hs1.fa \
Human-UCSC-T2T-CHM13-v2.0-hs1 -p 24
# 2.3 STAR index
# path：/data0/reference/Genome/Human/UCSC/T2T-CHM13-v2.0-hs1
cd /data0/reference/STAR_index/Human/UCSC-T2T-CHM13-v2.0-hs1
ln -s /home/xilab/reference/Genome/Human/UCSC/T2T-CHM13-v2.0-hs1/catLiftOffGenesV1.gtf ./
ln -s /home/xilab/reference/Genome/Human/UCSC/T2T-CHM13-v2.0-hs1/hs1.fa ./
STAR \
--runThreadN 24 \
--runMode genomeGenerate \
--genomeDir STARindex_100bp \
--genomeFastaFiles ./hs1.fa \
--sjdbGTFfile ./catLiftOffGenesV1.gtf \
--sjdbOverhang 100
```

**Run**:

``` sh
# Drosophlia dm6:
sh RNAseq-NEW.sh \
-p true
# Mouse Gencode-GRCm39-vM31:
sh RNAseq-NEW.sh \
-i /home/xilab/reference/STAR_index/Mouse/Gencode-GRCm39-vM31/STARindex_100bp \
-a /home/xilab/reference/hisat2_index/GRCm39-Gencode-vM31/genome_tran/genome_tran \
-r /home/xilab/reference/RSEM_index/mm10-Gencode-GRCm39-vM31/Genecode-GRCm39-vM31 \
-g /home/xilab/reference/Genome/Mouse/GENCODE/GRCm39-vM31/gencode.vM31.primary_assembly.annotation.gtf \
-b /home/xilab/reference/Genome/Mouse/GENCODE/GRCm39-vM31/gencode.vM31.primary_assembly.annotation.bed \
-p true 
# Human GENCODE-GRCh38-p13-release42:
sh RNAseq-NEW.sh \
-i /data0/reference/STAR_index/Human/GENCODE-GRCh38-p13-release42-primary/STARindex_100bp \
-a /home/xilab/reference/hisat2_index/grch38/genome_tran/genome_tran \
-r /home/xilab/reference/RSEM_index/Human-GENCODE-GRCh38-p13-release42-primary/GENCODE-GRCh38-p13-release42-Primary \
-g /data0/reference/Genome/Human/GENCODE/GRCh38.p13.release42/gencode.v42.primary_assembly.annotation.gtf \
-b /data0/reference/Genome/Human/GENCODE/GRCh38.p13.release42/gencode.v42.primary_assembly.annotation.bed \
-p true
# Human UCSC-T2T-CHM13-v2.0-hs1:
sh RNAseq-NEW.sh \
-i /data0/reference/STAR_index/Human/UCSC-T2T-CHM13-v2.0-hs1/STARindex_100bp \
-a /home/xilab/reference/hisat2_index/T2T-CHM13-v2.0-hs1-UCSC/genome_tran \
-r /home/xilab/reference/RSEM_index/Human-UCSC-T2T-CHM13-v2.0-hs1/Human-UCSC-T2T-CHM13-v2.0-hs1 \
-g /home/xilab/reference/Genome/Human/UCSC/T2T-CHM13-v2.0-hs1/catLiftOffGenesV1.gtf \
-b /home/xilab/reference/Genome/Human/UCSC/T2T-CHM13-v2.0-hs1/catLiftOffGenesV1.bed \
-p true
```

上游分析相关的其他分析请见：``sh RNAseq.sh -h``; 分析流程Debug信息请见文件：``RNAseq-Pepeline-Debug.md``。

### 上游分析： 无参转录组

见``无参转录组分析-202304.md``和``无参转录组分析-参考维科盟.md``。

### 下游分析： DESeq2 + ClusterProfiler

``DESeq2-ClusterProfiler.Rmd``

### 下游分析： EBSeq

支持无生物学重复。

``EBSeq-NoReplicates.R``


## 4. scRNAseq/spatialRNAseq

目前仅支持10×Genomics数据。使用CellRanger跑上游分析。

### Prepare the latest CellRanger reference

``scRNAseq/CellRanger/cellranger-make-latest-reference``

### 下游分析：Seurat

见目录： ``SelfMadeScripts/scRNAseq/SeuratV5``

### 结果展示： ShinyApps

见目录： ``spatialRNAseq/paceRanger+Seurat/Shiny-APP``和***ZhangShinyApps***里的scRNAseq Shiny Apps。A scRNAseq Demo: https://xilab.shinyapps.io/database/。 Thanks to https://www.shinyapps.io/ for holding this app.


## 5. CrisprScreen

### Mageck Based Pipeline

``Mageck-codes.sh``

### Negative screen

``Crispr-Negative-Screen.Rmd``

### my Rpackage: ``fentouxungui/CrisprNS``

``` r
if(!require("CrisprNS")){
  library("devtools")
  devtools::install_github("fentouxungui/CrisprNS")
}
library(CrisprNS)
```

## 6. MetaGenomics

### GrimmLab

support **Amplicon** and **Metagenomics** analysis.

## 7. DamIDseq

### 上游分析 - Bowtie2 Based

``DAMseq.sh``

### 下游分析 - iDEAR

``iDEAR-Drosophlia-for-SE.Rmd`` And ``iDEAR-Drosophlia-for-PE.Rmd``

## 8. RIPseq

上游分析流程参考是转录组分析。下游用的``RIPSeeker``R包。

## 9. Others

- ``extract_gene_summary_from_gbff_file.py``: 从gbff文件中提取Gene summary。
- ``extract_gene_summary_from_XML_file.py``: 从XML文件中提取Gene summary。
- ``AlphaFold-Download-glb-files.py``: 从AlphaFold官网下载蛋白的3D结构文件glb。
- ``batch_download_Entrez_gene_summary.py``：从Entrez自动下载Gene summary。
- ``PeptideAtlas-spider.py``： 从PeptideAtlas下载蛋白。
- 获取gRNA序列对应的基因信息
- 使用fasta和gtf文件提取基因的promoter序列，并统计特定序列出现的次数

## 10. WGBS

### 上游分析: WGBS-Bismark + msPIPE

``` sh
sh WGBS-Bismark.sh -p true
# msPIPE pipeline:
~/software/msPIPE/msPIPE/msPIPE.py -p params_galGal6.conf -c 12 --skip_calling --calling_data ./results/5_msPIPE/methylCALL \
-o ./results/5_msPIPE
~/software/msPIPE/msPIPE/msPIPE.py -p params_galGal6.conf -c 12 --skip_calling --calling_data ./results/5_msPIPE/methylCALL \
-o ./results/5_msPIPE --bsmooth --skip_GMA
```

### 上游分析： Bsmap 未完成

``WGBS-Bsmap.sh``

### 下游分析: methykit

``WGBS/methylKit/methylKit.Rmd``

