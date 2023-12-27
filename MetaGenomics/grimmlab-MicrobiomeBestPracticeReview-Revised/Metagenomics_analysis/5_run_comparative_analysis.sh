#!/usr/bin/evn bash

### Comparative Analysis


comparative_analysis_main(){
   # check_and_install
   contigs_classification_with_kraken
   comparative_functional_annotation_prokka
   gene_classification_with_diamond_megan
   predict_metabolic_pathway_minpath
}

contigs_classification_with_kraken(){

   echo "Running contigs classification with kraken"
   # running kraken
   mkdir -p ${ANALYSIS_FOLDER}/contigs_taxonomic_classification/kranken_output
   finallist=$(ls -d ${ANALYSIS_FOLDER}/QC/final_QC_output/*.1.unmerged.final.clean.fq | awk '{print $NF}')
   for i in ${finallist}
   do
      kraken2 \
      --db /data0/reference/MetaGenomics/kraken2/kraken2_db/k2_pluspf_20210517 \
      --threads $THREADS \
      --use-names \
      ${ANALYSIS_FOLDER}/assembly/spades/$(basename ${i%.1*})/contigs_1000_filtered.fasta \
      --report ${ANALYSIS_FOLDER}/contigs_taxonomic_classification/kranken_output/$(basename ${i%.1*}.kreport) > \
      ${ANALYSIS_FOLDER}/contigs_taxonomic_classification/kranken_output/$(basename ${i%.1*}.kraken)
   done

   echo "DONE running contigs classification with kraken"
}


comparative_functional_annotation_prokka(){
   echo "Running gene functional annotation using prokka"

   #Run prokka
   mkdir -p ${ANALYSIS_FOLDER}/Gene_based_analysis_onContigs/functional_classification/prokka_output

   finallist=$(ls -d ${ANALYSIS_FOLDER}/QC/final_QC_output/*.1.unmerged.final.clean.fq | awk '{print $NF}')
   for i in ${finallist}
   do
      #ln -s ~/mg-workshop/results/assembly/$SAMPLE/${SAMPLE}_${kmer}/contigs.fa .
      prokka ${ANALYSIS_FOLDER}/assembly/spades/$(basename ${i%.1*})/contigs_1000_filtered.fasta \
      --outdir ${ANALYSIS_FOLDER}/Gene_based_analysis_onContigs/functional_classification/prokka_output \
      --locustag 'PROKKA' \
      --kingdom 'Bacteria' \
      --prefix $(basename ${i%.1*}) \
      --norrna \
      --notrna \
      --metagenome \
      --cpus $THREADS \
      --addgenes \
      --centre X \
      --force

      grep \
         "eC_number=" ${ANALYSIS_FOLDER}/Gene_based_analysis_onContigs/functional_classification/prokka_output/$(basename ${i%.1*}).gff | cut -f9 | cut -f1,3 -d ';' | sed "s/ID=//g" | sed "s/;eC_number=/\t/g" > ${ANALYSIS_FOLDER}/Gene_based_analysis_onContigs/functional_classification/prokka_output/$(basename ${i%.1*}).ec
   done

   echo "DONE running gene functional annotation using prokka!"

}

gene_classification_with_diamond_megan(){
   echo "Running contigs classification with diamond"

   mkdir -p ${ANALYSIS_FOLDER}/contigs_taxonomic_classification/diamond_output
   finallist=$(ls -d ${ANALYSIS_FOLDER}/QC/final_QC_output/*.1.unmerged.final.clean.fq | awk '{print $NF}')
   for i in ${finallist}
   do
      diamond blastp \
      --threads $THREADS \
      --query ${ANALYSIS_FOLDER}/Gene_based_analysis_onContigs/functional_classification/prokka_output/$(basename ${i%.1*}).faa \
      --db /data0/reference/MetaGenomics/diamond_nr/nr.dmnd \
      --daa ${ANALYSIS_FOLDER}/contigs_taxonomic_classification/diamond_output/$(basename ${i%.1*}.daa)


      daa2rma \
      --in ${ANALYSIS_FOLDER}/contigs_taxonomic_classification/diamond_output/$(basename ${i%.1*}.daa) \
      --mapDB /data0/reference/MetaGenomics/megan_ref/megan-map-Feb2022.db \
      -fwa true \
      --threads $THREADS \
      --out ${ANALYSIS_FOLDER}/contigs_taxonomic_classification/diamond_output/$(basename ${i%.1*}.rma)
      # daa2rma \
      # --in ${ANALYSIS_FOLDER}/contigs_taxonomic_classification/diamond_output/$(basename ${i%.1*}.daa) \
      # --acc2taxa ${REFERENCE_FOLDER}/reference_database/megan_ref/prot_acc2tax-Nov2018X1.abin \
      # --acc2interpro2go ${REFERENCE_FOLDER}/reference_database/megan_ref/acc2interpro-June2018X.abin \
      # --acc2seed  ${REFERENCE_FOLDER}/reference_database/megan_ref/acc2seed-May2015XX.abin \
      # --acc2eggnog ${REFERENCE_FOLDER}/reference_database/megan_ref/acc2eggnog-Oct2016X.abin \
      # -fwa true \
      # --out ${ANALYSIS_FOLDER}/contigs_taxonomic_classification/diamond_output/$(basename ${i%.1*}.rma)
   done
   
   echo "DONE running contigs classification with diamond"
}


predict_metabolic_pathway_minpath(){
   echo "Running predict metabolic pathway"

   #awk -F "\t" '{print $1 "\t" $3}' 121832.assembled.EC > for_MinPath.ec
   #sed -i "s/EC://g" for_MinPath.ec
   #MinPath1.2.py -any for_MinPath.ec -map /home/rprops/MinPath/data/ec2path -report report.metacyc.minpath

   mkdir -p ${ANALYSIS_FOLDER}/Gene_based_analysis_onContigs/functional_classification/minpath_output

   finallist=$(ls -d ${ANALYSIS_FOLDER}/QC/final_QC_output/*.1.unmerged.final.clean.fq | awk '{print $NF}')
   for i in ${finallist}
   do
      python `which MinPath.py` \
         -any  ${ANALYSIS_FOLDER}/Gene_based_analysis_onContigs/functional_classification/prokka_output/$(basename ${i%.1*}).ec\
         -map /home/xilab/software/MinPath/MinPath/data/ec2path \
         -report ${ANALYSIS_FOLDER}/Gene_based_analysis_onContigs/functional_classification/minpath_output/$(basename ${i%.1*}).report.metacyc.minpath
   done

   echo "DONE running predict metabolic pathway!"
}

comparative_analysis_main
