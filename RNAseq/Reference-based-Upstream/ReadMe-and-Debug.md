# ReadMe

> 注意:
>
> 无论双端单端，reads长度最好大于50bp，因为rsem-calculate-expression --seed-length 默认25bp，短于25bp的会警告。
>
> 使用流程时，注意fastq文件大小不要大于2.8G！详见Debug信息第一条。
>
> Hisat2 比对很慢，可能是有大量reads比对到多个位置，估计是有大量核糖体mRNA。
>
> 计划：
>
> 添加重启功能？

## Analysis Roadmap

1. STAR + FeatureCounts (counts matrix - Gene level)     > DESeq2/edgeR/limma（counts）
2. STAR + RSEM (counts, TPM, FPKM matrix - Gene level)   > DESeq2/edgeR/limma（counts）
                                                                               |-  Known             > DESeq2/edgeR/limma（counts）
3. Hisat2 + StringTie (counts - Gene/Transcript）, TPM, FPKM matrix - Gene)   -| 
                                                                                |-  Known and de novo > DESeq2/edgeR/limma （counts）or balgown
                                                                            注： 1. ./results/6_Counts_StringTie/Make-Reference/gffcmp.annotated.gtf 为从头方式得到的全部转录本和外显子
2. Known and de novo 的结果可考虑导入到 IsoformSwitchAnalyzeR R包进行分析。未尝试！ 不同处理下，isform的切换。
3. Known and de novo 的结果包含以MSTRG为前缀的新基因/转录本，Known的没有，这一点分析时需要注意，对于Known and de novo 的结果可以用IsoformSwitchAnalyzeR导入，
参考（https://www.biostars.org/p/282817/），说是可以解决50%的MSTRG基因的注释，然后再导出count matrix，用于DESeq2/edgeR/limma分析。

之后考虑搭建balgown分析流程！

Softwares/Functions needed:
1. FASTQ QC: fastqc;
2. Trim Adapter: Trimgalore;
3. Align reads: STAR; Hisat2; Samtools;
4. BW files: deepTools: bamCoverage, plotCorrelation, multiBamSummary; 
5. RNAseq QC: RSeQC: bam_stat.py, read_distribution.py, geneBody_coverage.py(Not Used);
6. Summarizing Gene/Transcrpt/Exon Counts: Subread: featureCounts; GffCompare + StringTie;
7. Analysis Report: multiqc



```shell
# 果蝇案例
## UCSC - dm6
### 1.1 hisat2_index 可以从hisat2官网下载 http://daehwankimlab.github.io/hisat2/download/
#### 构建 Hisat2 index: dm6 genome_tran(hisat2 官网未提供)
cp /home/xilab/reference/Genome/Drosophlia/UCSC/dm6/dm6.fa ./genome.fa
cp /home/xilab/reference/Genome/Drosophlia/UCSC/dm6/dm6.refGene.gtf  ./genome.gtf
hisat2_extract_splice_sites.py genome.gtf > genome.ss
hisat2_extract_exons.py genome.gtf > genome.exon
hisat2-build -p 32 --exon genome.exon --ss genome.ss genome.fa genome_tran

### 1.2 构建STAR index
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

### 1.3 构建RSEM index
rsem-prepare-reference \
--gtf ./dm6.refGene.gtf \
./dm6.fa \
dm6-UCSC -p 24


# 小鼠案例
## GENCODE - GRCm39 vM31
1.1 hisat2_index 可以从hisat2官网下载，Hisat2官网未提供GRCm39(mm39)的index
### 构建Hisat2 index
### 参考： http://daehwankimlab.github.io/hisat2/howto/#building-indexes
cp /home/xilab/reference/Genome/Mouse/GENCODE/GRCm39-vM31/GRCm39.genome.fa ./genome.fa
cp /home/xilab/reference/Genome/Mouse/GENCODE/GRCm39-vM31/gencode.vM31.primary_assembly.annotation.gtf  ./genome.gtf
hisat2_extract_splice_sites.py genome.gtf > genome.ss
hisat2_extract_exons.py genome.gtf > genome.exon
hisat2-build -p 32 --exon genome.exon --ss genome.ss genome.fa genome_tran

### 1.2 构建STAR index
ln -s /home/xilab/reference/Genome/Mouse/GENCODE/GRCm39-vM31/gencode.vM31.primary_assembly.annotation.gtf 
ln -s /home/xilab/reference/Genome/Mouse/GENCODE/GRCm39-vM31/GRCm39.primary_assembly.genome.fa ./
STAR \
--runThreadN 24 \
--runMode genomeGenerate \
--genomeDir STARindex_100bp \
--genomeFastaFiles ./GRCm39.primary_assembly.genome.fa \
--sjdbGTFfile ./gencode.vM31.primary_assembly.annotation.gtf \
--sjdbOverhang 100
# sjdbOverhang: 测序的reads长度-1， 默认是100bp，貌似默认的这个效果也不错！

### 1.3 构建RSEM index
rsem-prepare-reference \
--gtf /home/xilab/reference/Genome/Mouse/GENCODE/GRCm39-vM31/gencode.vM31.primary_assembly.annotation.gtf  \
/home/xilab/reference/Genome/Mouse/GENCODE/GRCm39-vM31/GRCm39.primary_assembly.genome.fa \
Genecode-GRCm39-vM31 -p 24

# 人案例
## GENCODE - GRCh38.p13.release42 Primary
### 1.1 hisat2_index 可以从hisat2官网下载
### http://daehwankimlab.github.io/hisat2/download/
### 1.2 RSEM indx
rsem-prepare-reference \
--gtf /data0/reference/Genome/Human/GENCODE/GRCh38.p13.release42/gencode.v42.primary_assembly.annotation.gtf \
/data0/reference/Genome/Human/GENCODE/GRCh38.p13.release42/GRCh38.primary_assembly.genome.fa \
GENCODE-GRCh38-p13-release42-Primary -p 24

### 1.3 STAR index
### path：/data0/reference/Genome/Human/GENCODE/GRCh38.p13.release42
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

## UCSC - T2T-CHM13-v2.0-hs1
### 2.1 Hisat2 index
cp /home/xilab/reference/Genome/Human/UCSC/T2T-CHM13-v2.0-hs1/hs1.fa ./genome.fa
cp /home/xilab/reference/Genome/Human/UCSC/T2T-CHM13-v2.0-hs1/catLiftOffGenesV1.gtf  ./genome.gtf
hisat2_extract_splice_sites.py genome.gtf > genome.ss
hisat2_extract_exons.py genome.gtf > genome.exon
hisat2-build -p 32 --exon genome.exon --ss genome.ss genome.fa genome_tran

### 2.2 RSEM indx
rsem-prepare-reference \
--gtf /home/xilab/reference/Genome/Human/UCSC/T2T-CHM13-v2.0-hs1/catLiftOffGenesV1.gtf \
/home/xilab/reference/Genome/Human/UCSC/T2T-CHM13-v2.0-hs1/hs1.fa \
Human-UCSC-T2T-CHM13-v2.0-hs1 -p 24

### 2.3 STAR index
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

# 补充分析
## 1. reads 比对率低，判定是污染了什么物种，注意others，一般是的是接头序列，引物序列等。
cd results
mkdir S_Kraken2
cat sample_fastq_meta.txt | sed '1d' | while read sample fastq
do
kraken2 --db /data0/reference/MetaGenomics/kraken2/kraken2_db/k2_pluspf_fly_mouse_20220501 \
--threads 24 --use-names $fastq --report ./S_Kraken2/$(basename ${fastq}).kraken2.report \
> ./S_Kraken2/$(basename ${fastq}).kraken2.output
done
# 再写个代码，提取top10的结果，或者图形展示。
```



# Debug

## 滞留在Cutadapt步骤

停止在cutadapt步骤，等候一天无任何反应。
逐行运行代码发现，原因是停留在fastqc一步，尤其是对于较大的样本（目前发现的是双端150bp，文件大小各自为2.9G），报错信息：

>Now running FastQC on the validated data BratOE-Rep2_val_1.fq.gz<<<
>
>Started analysis of BratOE-Rep2_val_1.fq.gz
>Exception in thread "Thread-1" java.lang.OutOfMemoryError: Java heap space
>	at uk.ac.babraham.FastQC.Utilities.QualityCount.<init>(QualityCount.java:33)
>	at uk.ac.babraham.FastQC.Modules.PerTileQualityScores.processSequence(PerTileQualityScores.java:281)
>	at uk.ac.babraham.FastQC.Analysis.AnalysisRunner.run(AnalysisRunner.java:89)
>	at java.lang.Thread.run(Thread.java:745)

经查发现：https://github.com/s-andrews/FastQC/issues/86
使用代码 'fastqc BratOE-Rep2_val_1.fq.gz -t 24'错误会消失！
总结： 使用流程时，**注意双端测序数据的每个文件大小不要大于2.8G！**

