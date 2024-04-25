#!/usr/bin/evn bash

### binrefinement

binrefinement_main(){
   # check_and_install
   run_binrefiner
}


run_binrefiner(){
   echo "Running Binrefiner"

   mkdir -p ${ANALYSIS_FOLDER}/bin_refinement/binrefiner_output/
   mkdir -p ${ANALYSIS_FOLDER}/bin_refinement/input/maxbin
   mkdir -p ${ANALYSIS_FOLDER}/bin_refinement/input/metabat

   cp ${ANALYSIS_FOLDER}/binning/maxbin/*.fasta ${ANALYSIS_FOLDER}/bin_refinement/input/maxbin
   cp ${ANALYSIS_FOLDER}/binning/metabat/*.fa ${ANALYSIS_FOLDER}/bin_refinement/input/metabat

   # cd  ${ANALYSIS_FOLDER}/bin_refinement/
   Binning_refiner \
   -i ${ANALYSIS_FOLDER}/bin_refinement/input/ \
   -p BR

   #mv BR* ${ANALYSIS_FOLDER}/bin_refinement/binrefiner_output/

   #running CheckM
   mkdir -p ${ANALYSIS_FOLDER}/bin_refinement/checkM/maxbin/
   #running for maxbin
   checkm lineage_wf \
   --threads $THREADS \
   -x fasta ${ANALYSIS_FOLDER}/binning/maxbin/ \
   ${ANALYSIS_FOLDER}//bin_refinement/checkM/maxbin/

   checkm taxonomy_wf \
   domain \
   Bacteria \
   ${ANALYSIS_FOLDER}/binning/maxbin/ \
   ${ANALYSIS_FOLDER}/bin_refinement/checkM/maxbin/ \
   -x fasta \
   --threads $THREADS

   # Python 2 to 3: Removed Functionality
   # bin_qa_plot: non-critical, rarely used plot which does not scale to the large numbers 
   # of MAGs now being recovered
   # mkdir -p ${ANALYSIS_FOLDER}/bin_refinement/checkM/maxbin/plots/
   # checkm bin_qa_plot \
   # -x fasta \
   # ${ANALYSIS_FOLDER}/bin_refinement/checkM/ \
   # ${ANALYSIS_FOLDER}/binning/maxbin/ \
   # ${ANALYSIS_FOLDER}/bin_refinement/checkM/maxbin/plots


   mkdir -p ${ANALYSIS_FOLDER}/bin_refinement/checkM/metabat
   #runrning for metabat
   checkm lineage_wf \
   -x fa \
   --threads $THREADS \
   ${ANALYSIS_FOLDER}/binning/metabat/ \
   ${ANALYSIS_FOLDER}/bin_refinement/checkM/metabat

   mkdir -p ${ANALYSIS_FOLDER}/bin_refinement/checkM/metabat/plots/
   # 注意，与maxbin不同，后缀为fa
   checkm taxonomy_wf \
   domain \
   Bacteria \
   ${ANALYSIS_FOLDER}/binning/metabat/ \
   ${ANALYSIS_FOLDER}/bin_refinement/checkM/metabat/ \
   -x fa \
   --threads $THREADS

   # checkm bin_qa_plot \
   # -x fasta \
   # ${ANALYSIS_FOLDER}/bin_refinement/checkM/ \
   # ${ANALYSIS_FOLDER}/binning/metabat/ \
   # ${ANALYSIS_FOLDER}/bin_refinement/checkM/metabat/plots

   echo "DONE running Binrefiner"

}

binrefinement_main
