# use ensembl id
# use R to compare gtf file from ucsc and ensembl, and prepare the promoter bed
# 1kb from the tss
# 添加DMR、HyperDMC和HypoDMC对应的基因列表

bed=/home/xilab/zhangyc/WGBS/results/5_msPIPE/Analysis-Ensembl/annotations/promoter.bed
output=/home/xilab/zhangyc/WGBS/results/5_msPIPE/Analysis-Ensembl/DMR/

cd /home/xilab/zhangyc/WGBS/results/5_msPIPE/Analysis-Ensembl/DMR
ls -d */ | while read id
do
### DMC
# bedtools intersect -wa -wb -a $bed -b ${output}${id}/reform.DMC_q0.5.bed > ${output}${id}/intersection.DMC2Promoter.Ensembl.txt
# cut -f 4 ${output}${id}/intersection.DMC2Promoter.Ensembl.txt | cut -d ';' -f1| cut -d ':' -f2 | sort -u > ${output}${id}/DMC_genelist.txt

### DMR
sed '1d' ${output}${id}/DMR_0.5.bed > ${output}${id}/reform.DMR_0.5.bed
bedtools intersect -wa -wb -a $bed -b ${output}${id}/reform.DMR_0.5.bed > ${output}${id}/intersection.DMR2Promoter.txt
cut -f 4 ${output}${id}/intersection.DMR2Promoter.txt | cut -d ';' -f1| cut -d ':' -f2 | sort -u > ${output}${id}/DMR_genelist.txt


################  Hyper and Hypo
cut -f 4 ${output}${id}/hyperDMC_detailed_count_methyl.txt | cut -d ':' -f2 | sort -u > ${output}${id}/hyperDMC_genelist.txt
cut -f 4 ${output}${id}/hypoDMC_detailed_count_methyl.txt | cut -d ':' -f2 | sort -u > ${output}${id}/hypoDMC_genelist.txt

done

