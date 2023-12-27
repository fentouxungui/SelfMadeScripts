head all_methylCalls.txt 
# chr location total reads readsWithMethylation
chr1	1442	14	8
chr1	8830	10	2
chr1	10459	10	4
chr1	13272	11	11
chr1	13310	16	10
chr1	13369	14	6
chr1	13388	12	3
chr1	13451	13	3
chr1	13475	22	7
chr1	13504	23	13

find ../ -name *_CpG.cov.txt | while read id
do
new=$(basename $id)
awk '{print $1"\t"$2"\t"$5+$6"\t"$5}' $id >  ${new/CpG.cov.txt/all_methylCalls.txt}
done



ls *_all_methylCalls.txt | while read id
do
~/software/msPIPE/msPIPE/bin/GMA/GMA.Union_txt2bed.pl ${id} > ${id/all_methylCalls.txt/CpG_methylCalls.bed}
done

ls *CpG_methylCalls.bed | while read id
do
bedtools intersect -wo -a /home/xilab/zhangyc/WGBS/results/5_msPIPE/Analysis/annotations/Genomic_context.bed -b ${id} > ${id/CpG_methylCalls.bed/methylation.Genomic_Context.CpG.txt}
done



ls *methylation.Genomic_Context.CpG.txt | while read id
do
old=$(basename $id)
new=${old%_methylation.Genomic_Context.CpG.txt}
Rscript /home/xilab/software/msPIPE/msPIPE/bin/vis_script/genomic_context_levels.R $id ${new} ${new}.Genomic_Context_CpG.pdf 1> ${new}.Avg_Genomic_Context_CpG.txt
done