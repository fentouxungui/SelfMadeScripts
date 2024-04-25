# use ensembl id
# use R to compare gtf file from ucsc and ensembl, and prepare the promoter bed
# 1kb from the tss
# 添加DMR、HyperDMC和HypoDMC对应的基因列表

# bed=/home/xilab/zhangyc/WGBS/results/5_msPIPE/R-Ensembl-prepared.promoter.bed
bed=/home/xilab/zhangyc/WGBS/results/5_msPIPE/R-Ensembl-Release100.2k.promoter.bed
output=/data3/Xilab-Data-Analysis/zhangyc/WGBS/results/5_msPIPE/Analysis/DMR/
cd /home/xilab/zhangyc/WGBS/results/5_msPIPE/Analysis/DMR
ls -d */ | while read id
do
### DMC
# q0.5
bedtools intersect -wa -wb -a $bed -b ${output}${id}q0.5/reform.DMC_q0.5.bed > ${output}${id}q0.5/intersection.DMC2Promoter.Ensembl.2k.txt

cut -f 4 ${output}${id}q0.5/intersection.DMC2Promoter.Ensembl.2k.txt | cut -d ';' -f1| cut -d ':' -f2 | sort -u > ${output}${id}q0.5/DMC_genelist_Ensembl.2k.txt

# q0.01
bedtools intersect -wa -wb -a $bed -b ${output}${id}q0.01/reform.DMC_q0.01.bed > ${output}${id}/q0.01/intersection.DMC2Promoter.Ensembl.2k.txt

cut -f 4 ${output}${id}q0.01/intersection.DMC2Promoter.Ensembl.2k.txt | cut -d ';' -f1| cut -d ':' -f2 |sort -u > ${output}${id}q0.01/DMC_genelist_Ensembl.2k.txt

### DMR
##### 0.5

bedtools intersect -wa -wb -a $bed -b ${output}${id}q0.5/reform.DMR_0.5.bed > ${output}${id}q0.5/intersection.DMR2Promoter.Ensembl.2k.txt

cut -f 4 ${output}${id}q0.5/intersection.DMR2Promoter.Ensembl.2k.txt | cut -d ';' -f1| cut -d ':' -f2 | sort -u > ${output}${id}q0.5/DMR_genelist.Ensembl.2k.txt

###### 0.01

bedtools intersect -wa -wb -a $bed -b ${output}${id}q0.01/reform.DMR_0.01.bed > ${output}${id}q0.01/intersection.DMR2Promoter.Ensembl.2k.txt

cut -f 4 ${output}${id}q0.01/intersection.DMR2Promoter.Ensembl.2k.txt | cut -d ';' -f1| cut -d ':' -f2 | sort -u > ${output}${id}q0.01/DMR_genelist.Ensembl.2k.txt



################  Hyper and Hypo


/home/xilab/software/msPIPE/msPIPE/bin/script/getStat.pl ${output}${id}q0.5/intersection.DMC2Promoter.Ensembl.2k.txt ${output}${id}test
ls ${output}${id}test | while read sample
do
mv ${output}${id}test/${sample} ${output}${id}q0.5/${sample/.txt/_Ensembl.2k.txt}
done


/home/xilab/software/msPIPE/msPIPE/bin/script/getStat.pl ${output}${id}q0.01/intersection.DMC2Promoter.Ensembl.2k.txt ${output}${id}test
ls ${output}${id}test | while read sample
do
mv ${output}${id}test/${sample} ${output}${id}q0.01/${sample/.txt/_Ensembl.2k.txt}
done

cut -f 4 ${output}${id}q0.5/hyperDMC_detailed_count_methyl_Ensembl.2k.txt | cut -d ':' -f2 | sort -u > ${output}${id}q0.5/hyperDMC_genelist_Ensembl.2k.txt
cut -f 4 ${output}${id}q0.5/hypoDMC_detailed_count_methyl_Ensembl.2k.txt | cut -d ':' -f2 | sort -u > ${output}${id}q0.5/hypoDMC_genelist_Ensembl.2k.txt
cut -f 4 ${output}${id}q0.01/hyperDMC_detailed_count_methyl_Ensembl.2k.txt | cut -d ':' -f2 | sort -u > ${output}${id}q0.01/hyperDMC_genelist_Ensembl.2k.txt
cut -f 4 ${output}${id}q0.01/hypoDMC_detailed_count_methyl_Ensembl.2k.txt | cut -d ':' -f2 | sort -u > ${output}${id}q0.01/hypoDMC_genelist_Ensembl.2k.txt

done

