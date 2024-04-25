#!/usr/bin/evn bash

### coverage and bining 

bin_functional_classification_main(){
   # check_and_install
   run_binclassification
}



run_binclassification(){
   echo "Running Bin Classification"
   
   mkdir -p ${ANALYSIS_FOLDER}/bin_functional_annotation/prokka_out/metabat
   mkdir -p ${ANALYSIS_FOLDER}/bin_functional_annotation/prokka_out/maxbin
   # Run gene prediction using prodigal on
   #finallist=$(ll -d ${ANALYSIS_FOLDER}/QC/final_QC_output/*.1.unmerged.final.clean.fq | awk '{print $NF}')

   echo "Running prokka on maxbin"
   maxbinlist=$(ls -d ${ANALYSIS_FOLDER}/binning/maxbin/*.fasta | awk '{print $NF}')
   for i in $(ls ${maxbinlist})
   do
      bin_name=$(basename ${i%.*})

      prokka ${i} \
      --quiet \
      --cpus $THREADS \
      --outdir ${ANALYSIS_FOLDER}/bin_functional_annotation/prokka_out/maxbin/$bin_name \
      --prefix $bin_name \
      --metagenome \
      --kingdom 'Bacteria' \
      --locustag 'PROKKA' \
      --force
      
   done  

   echo "Running prokka on metabat"
   metabatbinlist=$(ls -d ${ANALYSIS_FOLDER}/binning/metabat/*.fa | awk '{print $NF}')
   for bin in ${metabatbinlist}
   do
      bin_name=$(basename ${i%.*})

      prokka ${i} \
      --quiet \
      --cpus $THREADS \
      --outdir ${ANALYSIS_FOLDER}/bin_functional_annotation/prokka_out/metabat/$bin_name \
      --prefix $bin_name \
      --metagenome \
      --kingdom 'Bacteria' \
      --locustag 'PROKKA' \
      --force
   done  

   echo "DONE running bin classification"
}

bin_functional_classification_main
