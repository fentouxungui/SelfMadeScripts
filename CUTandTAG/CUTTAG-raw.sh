mkdir -p 1_QC_fastqc
cd 1_QC_fastqc
fastqc -t 2 ../fastq/*.gz -o ./

mkdir -p 2_Trim-TrimGalore
cd 2_Trim-TrimGalore

trim_galore --cores 6 --fastqc --quality 20 --length 35 \
--paired ../fastq/86_1.fq.gz ../fastq/86_2.fq.gz \
--basename 86 --output_dir ./

bowtie2 --local --very-sensitive --no-mixed --no-discordant --phred33 -I 10 -X 700 \
-p 24 -x /data0/reference/bowtie2_index/hg38/GRCh38_noalt_as/GRCh38_noalt_as \
-1 ../2_Trim-TrimGalore/86_val_1.fq.gz \
-2 ../2_Trim-TrimGalore/86_val_2.fq.gz \
-S 86_bowtie2.sam \
&> 86_bowtie2.log


bowtie2 --local --very-sensitive --no-mixed --no-discordant --phred33 -I 10 -X 700 \
-p 24 -x /data0/reference/bowtie2_index/hg38/GRCh38_noalt_as/GRCh38_noalt_as \
-1 ../2_Trim-TrimGalore/C_val_1.fq.gz \
-2 ../2_Trim-TrimGalore/C_val_2.fq.gz \
-S C_bowtie2.sam \
&> C_bowtie2.log


# bowtie2-build E.coli-gennome-sequence.fasta E.coli-genome


# bowtie2 --local --very-sensitive --no-mixed --no-discordant --phred33 -I 10 -X 700 \
# -p 24 -x /home/xilab/reference/bowtie2_index/E.coli/NCBI/strainK12-substrainMG1655/E.coli-genome \
# -1 ../2_Trim-TrimGalore/86_val_1.fq.gz \
# -2 ../2_Trim-TrimGalore/86_val_2.fq.gz \
# -S 86_Ecoli.sam \
# &> 86_Ecoli.txt

# seqdepth=$(head -n1 86_Ecoli.txt | cut -d ' ' -f1)

# bowtie2 --local --very-sensitive --no-mixed --no-discordant --phred33 -I 10 -X 700 \
# -p 24 -x /data0/reference/bowtie2_index/hg38/GRCh38_noalt_decoy_as/GRCh38_noalt_decoy_as \
# -1 ../2_Trim-TrimGalore/C_val_1.fq.gz \
# -2 ../2_Trim-TrimGalore/C_val_2.fq.gz \
# -S C_Ecoli.sam \
# &> C_Ecoli.txt


# ls *.sam | while read id
# do
#     samtools view -u -@ 4 $id | samtools sort -@ 4 > ${id%sam}sorted.bam
# done

ls *.sam | while read id
do
    picard SortSam -I $id -O ${id%sam}sorted.bam -SO coordinate -CREATE_INDEX TRUE
    samtools idxstats ${id%sam}sorted.bam >& ${id%sam}idxstat
done

grep "chrM" *idxstat


# Error in MarkDuplicates
# > https://gatk.broadinstitute.org/hc/en-us/community/posts/24452455341467-Error-in-MarkDuplicates
ls *sorted.bam | while read id
do
    samtools addreplacerg -@ 8 -r "@RG\tID:RG1\tSM:SampleName\tPL:Illumina\tLB:Library.fa" -o ${id%bam}modified.bam $id
    picard MarkDuplicates -QUIET true -INPUT ${id%bam}modified.bam -OUTPUT ${id%sorted.bam}marked.bam -METRICS_FILE \
    ${id%sorted.bam}markDup.metrics -REMOVE_DUPLICATES false -CREATE_INDEX true -VALIDATION_STRINGENCY LENIENT -TMP_DIR .
done    

ls *markDup.metrics | while read id
do
    head -n 8 $id | cut -f 7,9 | grep -v ^# | tail -n 2
done

ls *marked.bam | while read id
do
    qualimap bamqc -bam $id -gd hg38 -outdir . -outfile ${id%sorted.modified.bam}qualimap.report -outformat html
done


ls *marked.bam | while read id
do
    samtools index $id -@ 4
    samtools idxstats $id >& ${id%marked.bam }idxstat
    samtools view -@ 4 -h $id | grep -v chrM | samtools sort -@ 4 -O bam -o ${id%bam}rmChrM.bam
    samtools view -@ 4 -h -b -F 1024 ${id%bam}rmChrM.bam > ${id%bam}rmChrM.rmDup.bam
    samtools view -@ 4 -q 30 -b ${id%bam}rmChrM.rmDup.bam > ${id%bam}rmChrM.rmDup.filtered.bam
done


ls *filtered.bam | while read id
do
    picard SortSam -I $id -O ${id%bam}sorted.bam -SO queryname -CREATE_INDEX TRUE
    bedtools bamtobed -i ${id%bam}sorted.bam -bedpe > ${id%bam}bed
done

# # 这里有问题，.打头的条目是怎么来的？
# ls *rmChrM.rmDup.filtered.bed | while read id
# do
#     cat $id | grep -v "^\\." > ${id%bed}modified.bed
# done

ls *filtered.bed | while read id
do
    bedtools sort -i $id > ${id%bed}sorted.bed
    cut -f 1,2,6  ${id%bed}sorted.bed | sort -k1,1 -k2,2n -k3,3n > ${id%bed}sorted-converted.bed
done

ls *sorted-converted.bed | while read id
do
    bedtools genomecov -bg -i $id -g /home/xilab/reference/bowtie2_index/hg38/GRCh38_noalt_as/chrom.sizes > ${id%bed}bedGraph
done

conda activate /home/xilab/software/miniconda-envs/bioinfo2

ls *bedGraph | while read id
do
    bedGraphToBigWig $id /home/xilab/reference/bowtie2_index/hg38/GRCh38_noalt_as/chrom.sizes ${id%bedGraph}bigWig
done

conda deactivate

# warnings: some chromosome names did not match 
computeMatrix scale-regions -S <sample1>.bigWig <sample2>.bigWig -R hg38-gene-coordinates.bed \
--beforeRegionStartLength 3000 --regionBodyLength 5000 --afterRegionStartLength 3000 --missingDataAsZero --skipZeros -o matrix.mat.gz


plotHeatmap -m matrix.mat.gz -out ExampleHeatmap1.png 


#>https://github.com/FredHutch/SEACR/blob/master/SEACR_1.3.sh
SEACR_1.3.sh ../3_Mapping-Bowtie2/86_bowtie2.rmChrM.rmDup.filtered.sorted-converted.bedGraph \
../3_Mapping-Bowtie2/C_bowtie2.rmChrM.rmDup.filtered.sorted-converted.bedGraph non stringent 86-vs-C_seacr.peaks

SEACR_1.3.sh ../3_Mapping-Bowtie2/86_bowtie2.rmChrM.rmDup.filtered.sorted-converted.bedGraph \
../3_Mapping-Bowtie2/C_bowtie2.rmChrM.rmDup.filtered.sorted-converted.bedGraph norm stringent 86-vs-C_seacr-norm.peaks



conda activate /home/xilab/software/miniconda-envs/bioinfo2

macs2 callpeak -t ../3_Mapping-Bowtie2/86_bowtie2.rmChrM.rmDup.filtered.sorted.bam \
-c ../3_Mapping-Bowtie2/C_bowtie2.rmChrM.rmDup.filtered.sorted.bam \
-g hs -f BAMPE -n 86-vs-C_macs2-peak-q0.1 --outdir ./ -q 0.1 2>./86-vs-C_macs2-peak-q0.1_summary.txt

conda deactivate




