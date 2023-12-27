#!/usr/bin/evn bash

THREADS=32

### Assembly Analysis
coassembly_main(){
   check_and_install
   run_coassembly
   #assembly_stats
   filter_on_assembly_length
}


# Run assembly 
run_coassembly(){
   echo "Running all coassemblies"
   echo "Running coassembly"
   mkdir -p ${ANALYSIS_FOLDER}/coassembly
   run_coassem_megahit
   run_coassem_spades
   # run_coassem_ray
   #run_coassem_idba
   echo "DONE running all coassemblies!"
}

run_coassem_idba(){
   #Assembly with IDBA
   #first covert paired into single
   #1. Merged the paired end data and covert to fasta format.
   echo "Running idba/fq2fa"
   mkdir -p ${ANALYSIS_FOLDER}/coassembly/idba_ud
   finallist=$(ls -d ${ANALYSIS_FOLDER}/QC/final_QC_output/*.1.unmerged.final.clean.fq | awk '{print $NF}')
   for i in ${finallist}
   do
      fq2fa \
      --merge ${i} ${i%1*}2.unmerged.final.clean.fq \
      ${i%1*}12.final.clean.fa

      idba_ud \
      -l ${i%1*}12.final.clean.fa \
      --mink 20 \
      --maxk 124 \
      --num_threads 16 \
      -o ${ANALYSIS_FOLDER}/coassembly/idba_ud/$(basename ${i%.1*})
   done
}

run_coassem_ray(){
   mkdir -p ${ANALYSIS_FOLDER}/coassembly/ray
   #Assembly with Ray
   echo "Running ray"
   finallist=$(ls -d ${ANALYSIS_FOLDER}/QC/final_QC_output/*.1.unmerged.final.clean.fq | awk '{print $NF}')
   for i in ${finallist}
   do
      mpiexec -n 16 \
      ray \
      - k 31 \
      -p ${i} ${i%1*}2.unmerged.final.clean.fq \
      -s ${i%1*}u.final.clean.fq  \
      -o ${ANALYSIS_FOLDER}/coassembly/ray/ray_31_assembly/$(basename ${i%.1*})
   done
}

run_coassem_spades(){
   mkdir -p ${ANALYSIS_FOLDER}/coassembly/spades
   # Assembly with spades
   finallist=$(ls -d ${ANALYSIS_FOLDER}/QC/final_QC_output/*.1.unmerged.final.clean.fq | awk '{print $NF}')
   echo "Running spades"
   #   spades.py --pe1-1 lib1_forward_1.fastq --pe1-2 lib1_reverse_1.fastq \
   #   --pe1-1 lib1_forward_2.fastq --pe1-2 lib1_reverse_2.fastq \
   #   -o spades_output
   cospd_1="spades.py --meta "
   for i in ${finallist}
   do
      cospd="$cospd --pe1-1 ${i} --pe1-2 ${i%1*}2.unmerged.final.clean.fq --pe-s ${i%1*}u.final.clean.fq "
   done
   cospd_2=" -k auto -o ${ANALYSIS_FOLDER}/coassembly/spades  -t 16"
   $cospd_1 $cospd $cospd_2
}

run_coassem_megahit(){
   mkdir -p ${ANALYSIS_FOLDER}/coassembly/megahit
   # Assembly Megahit
   echo "Running megahit"
   R1s=`ls -d ${ANALYSIS_FOLDER}/QC/final_QC_output/*.1.unmerged.final.clean.fq | python -c 'import sys; print ",".join([x.strip() for x in sys.stdin.readlines()])'`
   R2s=`ls -d ${ANALYSIS_FOLDER}/QC/final_QC_output/*.2.unmerged.final.clean.fq | python -c 'import sys; print ",".join([x.strip() for x in sys.stdin.readlines()])'`
   Rs=`ls -d ${ANALYSIS_FOLDER}/QC/final_QC_output/*.u.final.clean.fq | python -c 'import sys; print ",".join([x.strip() for x in sys.stdin.readlines()])'`

   megahit \
   --pe-1 ${R1s} \
   --pe-2 ${R2s} \
   --pe-s ${Rs} \
   -t $THREADS \
   --presets meta-large \
   --mem-flag 2 \
   -o ${ANALYSIS_FOLDER}/coassembly/megahit/
}

assembly_stats(){
   echo "Running stats"
   # creating stats of all the assemblers
   #Stats with quast
   mkdir -p ${ANALYSIS_FOLDER}/coassembly/coassembly_stats

   metaquast.py \
   --threads $THREADS \
   --gene-finding \
   -o $ANALYSIS_FOLDER/coassembly/coassembly_stats \
   -l megahit,SPAdes,IDBA-UD \
   $ANALYSIS_FOLDER/coassembly/megahit/megahit_final.contigs.fa \
   $ANALYSIS_FOLDER/coassembly/spades/contigs.fasta \
   ${ANALYSIS_FOLDER}/coassembly/idba_ud/contig.fa
   #${ANALYSIS_FOLDER}/coassembly/ray/$(basename ${i%.1*})/ray_31_assembly/  \

   echo "DONE running stats!"
}


filter_on_assembly_length(){
   echo "Running assembly length filter"
   #running megahit
   #spades
   reformat.sh \
   in=${ANALYSIS_FOLDER}/coassembly/spades/contigs.fasta \
   out=${ANALYSIS_FOLDER}/coassembly/spades/contigs_1000_filtered.fasta \
   minlength=1000

   #megahit
   reformat.sh \
   in=${ANALYSIS_FOLDER}/coassembly/megahit/megahit_final.contigs.fa \
   out=${ANALYSIS_FOLDER}/coassembly/megahit/megahit_final_1000_filtered.fasta \
   minlength=1000

   #idba_ud
   #${TOOLS_FOLDER}/bbmap/reformat.sh \
   #in=${ANALYSIS_FOLDER}/coassembly/idba_ud/contig.fa \
   #out=${ANALYSIS_FOLDER}/coassembly/idba_ud/contig_1000_filtered.fa \
   #minlength=1000

   echo "DONE running assembly length filter!"

}


coassembly_main
