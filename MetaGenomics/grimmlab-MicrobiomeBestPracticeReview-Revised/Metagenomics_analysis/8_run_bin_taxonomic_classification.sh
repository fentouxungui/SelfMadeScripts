#!/usr/bin/evn bash

### coverage and bining 

bin_taxonomic_classification_main(){
   # check_and_install
   run_binclassification
}




run_binclassification(){
   echo "Running Bin Classification"
   
   mkdir -p ${ANALYSIS_FOLDER}/bin_taxonomic_classification/kraken_out/metabat
   mkdir -p ${ANALYSIS_FOLDER}/bin_taxonomic_classification/kraken_out/maxbin
   # Run gene prediction using prodigal on
   #finallist=$(ll -d ${ANALYSIS_FOLDER}/QC/final_QC_output/*.1.unmerged.final.clean.fq | awk '{print $NF}')

   echo "Running Bin Classification on Metabat"
   metabatbinlist=$(ls -d ${ANALYSIS_FOLDER}/binning/metabat/*.fa | awk '{print $NF}')
   for bin in ${metabatbinlist}
   do
      bin_name=$(basename ${bin%.*})
      kraken2 \
      --db /data0/reference/MetaGenomics/kraken2/kraken2_db/k2_pluspf_20210517 \
      --threads $THREADS \
      --use-names \
      ${bin} \
      --report ${ANALYSIS_FOLDER}/bin_taxonomic_classification/kraken_out/metabat/${bin_name}.kreport > \
      ${ANALYSIS_FOLDER}/bin_taxonomic_classification/kraken_out/metabat/${bin_name}.kraken
   done

   echo "Running Bin Classification on Maxbin"
   maxbinlist=$(ls -d ${ANALYSIS_FOLDER}/binning/maxbin/*.fasta | awk '{print $NF}')
   for bin in ${maxbinlist}
   do
      bin_name=$(basename ${bin%.*})
      kraken2 \
      --db /data0/reference/MetaGenomics/kraken2/kraken2_db/k2_pluspf_20210517 \
      --threads $THREADS \
      --use-names \
      ${bin} \
      --report ${ANALYSIS_FOLDER}/bin_taxonomic_classification/kraken_out/maxbin/${bin_name}.kreport > \
      ${ANALYSIS_FOLDER}/bin_taxonomic_classification/kraken_out/maxbin/${bin_name}.kraken
   done

   echo "DONE running bin classification"
}

bin_taxonomic_classification_main
