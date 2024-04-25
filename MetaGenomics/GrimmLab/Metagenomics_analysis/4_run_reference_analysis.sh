#!/usr/bin/evn bash

### Refereence based Analysis


ref_analysis_main(){
   # check_and_install
   run_kraken
   run_metaphlan
   run_diamond
   run_humann3
}


run_kraken(){
   echo "Running kraken2"
   mkdir -p ${ANALYSIS_FOLDER}/reference_classification/kraken
   # cd ${ANALYSIS_FOLDER}/reference_classification/kraken

   # running kraken
   finallist=$(ls -d  ${ANALYSIS_FOLDER}/QC/final_QC_output/*.1.unmerged.final.clean.fq | awk '{print $NF}')
   for i in ${finallist}
   do
      # 非常快，几秒钟
      kraken2 \
      --db /data0/reference/MetaGenomics/kraken2/kraken2_db/k2_pluspf_20210517 \
      --threads $THREADS \
      --use-names \
      --paired \
      ${i} \
      ${i%1*}2.unmerged.final.clean.fq \
      --report ${ANALYSIS_FOLDER}/reference_classification/kraken/$(basename ${i%1*}kreport) > \
      ${ANALYSIS_FOLDER}/reference_classification/kraken/$(basename ${i%1*}kraken)
   done
   echo "DONE running kraken2"
}


run_metaphlan(){
   echo "Running metaphlan"

   mkdir -p ${ANALYSIS_FOLDER}/reference_classification/metaphlan3
   mkdir -p ${ANALYSIS_FOLDER}/reference_classification/metaphlan3/profiled_samples
   mkdir -p ${ANALYSIS_FOLDER}/reference_classification/metaphlan3/output_images
   # cd ${ANALYSIS_FOLDER}/reference_classification/metaphlan3

   finallist=$(ls -d ${ANALYSIS_FOLDER}/QC/final_QC_output/*.1.unmerged.final.clean.fq | awk '{print $NF}')
   echo "Running metaphlan"
   for i in ${finallist}
   do
      metaphlan \
      --input_type fastq \
      ${i},${i%1*}2.unmerged.final.clean.fq,${i%1*}u.final.clean.fq \
      --bowtie2db /data0/reference/MetaGenomics/MetaPhlAn-3/mpa_v31 \
      --bt2_ps very-sensitive \
      --nproc $THREADS \
      --bowtie2out ${ANALYSIS_FOLDER}/reference_classification/metaphlan3/$(basename ${i%1*}bt2out) > \
      ${ANALYSIS_FOLDER}/reference_classification/metaphlan3/profiled_samples/$(basename ${i%1*}txt)
   done


   echo "Creating sample profiles using metaphlan"
   profilelist=$(ls -d ${ANALYSIS_FOLDER}/reference_classification/metaphlan3/profiled_samples/* | awk '{print $NF}')

   #Merge the output
   merge_metaphlan_tables.py \
   ${ANALYSIS_FOLDER}/reference_classification/metaphlan3/profiled_samples/*.txt > \
   ${ANALYSIS_FOLDER}/reference_classification/metaphlan3/merged_abundance_table.txt

   # create a species abundance table
   grep -E "(s__)|(^clade)" ${ANALYSIS_FOLDER}/reference_classification/metaphlan3/merged_abundance_table.txt \
   | grep -v "t__" | sed 's/^.*s__//g' > \
   ${ANALYSIS_FOLDER}/reference_classification/metaphlan3/merged_abundance_table_species.txt

   echo "Creating heatmap"
   profilelist=$(ls -d ${ANALYSIS_FOLDER}/reference_classification/metaphlan3/profiled_samples/* | awk '{print $NF}')
   # 由于样本数目少于2，未能测试以下代码，另外，是否需要去除第二列？
   # https://github.com/biobakery/biobakery/wiki/metaphlan3#overview： 原教程去除了第二列
   hclust2.py \
   -i ${ANALYSIS_FOLDER}/reference_classification/metaphlan3/merged_abundance_table_species.txt \
   -o ${ANALYSIS_FOLDER}/reference_classification/metaphlan3/output_images/abundance_heatmap_species.png \
   --ftop 25 \
   --f_dist_f braycurtis \
   --s_dist_f braycurtis \
   --cell_aspect_ratio 0.5 \
   -l \
   --flabel_size 6 \
   --slabel_size 6 \
   --max_flabel_len 100 \
   --max_slabel_len 100 \
   --minv 0.1 \
   --dpi 300


   # creating files for graphlan
   echo "Running graphlan"
   # conda install export2graphlan
   # conda install graphlan
   conda activate py27
   mkdir -p ${ANALYSIS_FOLDER}/reference_classification/metaphlan3/graphlan_output
   mkdir -p ${ANALYSIS_FOLDER}/reference_classification/metaphlan3/graphlan_output/output_images
   export2graphlan.py \
   --skip_rows 1,2 \
   -i ${ANALYSIS_FOLDER}/reference_classification/metaphlan3/merged_abundance_table.txt \
   --tree ${ANALYSIS_FOLDER}/reference_classification/metaphlan3/graphlan_output/merged_abundance.tree.txt \
   --annotation ${ANALYSIS_FOLDER}/reference_classification/metaphlan3/graphlan_output/merged_abundance.annot.txt \
   --most_abundant 100 \
   --abundance_threshold 1 \
   --least_biomarkers 10 \
   --annotations 5,6 \
   --external_annotations 7 \
   --min_clade_size 1

   graphlan_annotate.py \
   --annot ${ANALYSIS_FOLDER}/reference_classification/metaphlan3/graphlan_output/merged_abundance.annot.txt \
   ${ANALYSIS_FOLDER}/reference_classification/metaphlan3/graphlan_output/merged_abundance.tree.txt \
   ${ANALYSIS_FOLDER}/reference_classification/metaphlan3/graphlan_output/merged_abundance.xml

   graphlan.py \
   --dpi 300 \
   --size 15 \
   --format pdf \
   ${ANALYSIS_FOLDER}/reference_classification/metaphlan3/graphlan_output/merged_abundance.xml \
   ${ANALYSIS_FOLDER}/reference_classification/metaphlan3/graphlan_output/output_images/merged_abundance.pdf
   #  --external_legends
   conda deactivate
   echo "DONE running metaphlan!"
}

run_diamond(){
   echo "Running diamond"

   mkdir -p ${ANALYSIS_FOLDER}/reference_classification/diamond_output

   finallist=$(ls -d ${ANALYSIS_FOLDER}/QC/final_QC_output/*.12.interleaved.final.clean.fa | awk '{print $NF}')
   for i in ${finallist}
   do
      diamond blastx \
      --threads $THREADS \
      --query ${i} \
      --db /data0/reference/MetaGenomics/diamond_nr/nr.dmnd \
      --daa ${ANALYSIS_FOLDER}/reference_classification/diamond_output/$(basename ${i%.12*})

      # 问题： shell中，daa2rma命令无法运行，也无任何报错
      daa2rma \
      --in ${ANALYSIS_FOLDER}/reference_classification/diamond_output/$(basename ${i%.12*}.daa) \
      --mapDB /data0/reference/MetaGenomics/megan_ref/megan-map-Feb2022.db \
      -fwa true \
      --threads $THREADS \
      --out ${ANALYSIS_FOLDER}/reference_classification/diamond_output/$(basename ${i%.12*}.rma)
      # 旧版软件，不赞成使用
      # daa2rma \
      # --in ${ANALYSIS_FOLDER}/reference_classification/diamond_output/$(basename ${i%12*}.daa) \
      # --acc2taxa ${REFERENCE_FOLDER}/reference_database/megan_ref/nucl_acc2tax-Nov2018.abin \
      # --acc2interpro2go ${REFERENCE_FOLDER}/reference_database/megan_ref/acc2interpro-June2018X.abin \
      # --acc2seed  ${REFERENCE_FOLDER}/reference_database/megan_ref/acc2seed-May2015XX.abin \
      # --acc2eggnog ${REFERENCE_FOLDER}/reference_database/megan_ref/acc2eggnog-Oct2016X.abin \
      # -fwa true \
      # --out ${ANALYSIS_FOLDER}/reference_classification/diamond_output/$(basename ${i%.12*}.rma)
   done

   echo "DONE running diamond!"
}


run_humann3(){
   echo "Running humann3"

   mkdir -p ${ANALYSIS_FOLDER}/reference_classification/humann3
   finallist=$(ls -d ${ANALYSIS_FOLDER}/QC/final_QC_output/*.interleaved.final.clean.fa | awk '{print $NF}')
   for i in ${finallist}
   do
      humann3 \
      --input ${i} \
      # --metaphlan ${TOOLS_FOLDER}/metaphlan3/ \
      --output ${ANALYSIS_FOLDER}/reference_classification/humann3 \
      --nucleotide-database /data0/reference/MetaGenomics/HUMAnN3/chocophlan \
      --protein-database /data0/reference/MetaGenomics/HUMAnN3/uniref \
      --threads $THREADS
   done

   genefamilieslist=$(ls -d ${ANALYSIS_FOLDER}/reference_classification/humann3/*_genefamilies.tsv | awk '{print $NF}')
   for i in ${genefamilieslist}
   do
      humann_renorm_table \
      --input $i \
      --output ${i%_genefamilies*}_genefamilies_relab.tsv \
      --units relab

      humann_join_tables \
      --input ${ANALYSIS_FOLDER}/reference_classification/humann3 \
      --output ${i%_genefamilies*}_genefamilies.tsv \
      --file_name genefamilies_relab

      humann_join_tables \
      --input ${ANALYSIS_FOLDER}/reference_classification/humann3 \
      --output ${i%_genefamilies*}_pathcoverage.tsv \
      --file_name pathcoverage

      humann_join_tables \
      --input ${ANALYSIS_FOLDER}/reference_classification/humann3 \
      --output ${i%genefamilies*}_pathabundance.tsv \
      --file_name pathabundance_relab
   done

   echo "DONE running humann3!"
}

ref_analysis_main
