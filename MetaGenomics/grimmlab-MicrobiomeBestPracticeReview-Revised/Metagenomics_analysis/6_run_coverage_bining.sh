#!/usr/bin/evn bash

### coverage and bining 


coverage_and_bining_main(){
   # check_and_install
   create_coveragefile
   binning_maxbin
   binning_metabat
}


create_coveragefile(){
   echo "Running Create coverage file"

   mkdir -p ${ANALYSIS_FOLDER}/coverage
   mkdir -p ${ANALYSIS_FOLDER}/coverage/bam
   finallist=$(ls -d ${ANALYSIS_FOLDER}/QC/final_QC_output/*.1.unmerged.final.clean.fq | awk '{print $NF}')
   conda activate py37
   for i in ${finallist}
   do
      echo "Running bbwrap in coverage file"
      bbwrap.sh \
      ref=${ANALYSIS_FOLDER}/assembly/spades/$(basename ${i%.1*})/contigs_1000_filtered.fasta \
      nodisk \
      in=${i},${i%1*}u.final.clean.fq \
      in2=${i%1*}2.unmerged.final.clean.fq,NULL \
      t=24 \
      kfilter=22 \
      subfilter=15 \
      maxindel=80 \
      out=${ANALYSIS_FOLDER}/coverage/bam/$(basename ${i%1*}sam) \
      covstats=${ANALYSIS_FOLDER}/coverage/bam/$(basename ${i%1*}coverage)
      
      #converting sam file to bam
      samtools view \
      -S \
      -b \
      -u \
      ${ANALYSIS_FOLDER}/coverage/bam/$(basename ${i%1*}sam) \
      > ${ANALYSIS_FOLDER}/coverage/bam/$(basename ${i%1*}bam)

      #sorting the bam file
      samtools sort \
      -m 16G \
      -@ 3 \
      ${ANALYSIS_FOLDER}/coverage/bam/$(basename ${i%1*}bam) \
      -o ${ANALYSIS_FOLDER}/coverage/bam/$(basename ${i%1*}sorted.bam)

      #indexing the sorted bam file
      samtools index \
      ${ANALYSIS_FOLDER}/coverage/bam/$(basename ${i%1*}sorted.bam)
   done
   conda deactivate

   jgi_summarize_bam_contig_depths \
   --outputDepth ${ANALYSIS_FOLDER}/coverage/depth.txt \
   --pairedContigs ${ANALYSIS_FOLDER}/coverage/paired.txt \
   --minContigLength 1000 \
   --minContigDepth 2 \
   ${ANALYSIS_FOLDER}/coverage/bam/*sorted.bam

   tail -n+2 ${ANALYSIS_FOLDER}/coverage/depth.txt | cut -f 1,3 > ${ANALYSIS_FOLDER}/coverage/maxbin.cov

   echo "DONE running create coverage file!"

}

binning_maxbin(){
   echo "Running Binning MaxBin"

   #chose the contigs
   #samplename=${i%1*}
   mkdir -p ${ANALYSIS_FOLDER}/binning/
   mkdir -p ${ANALYSIS_FOLDER}/binning/maxbin

   finallist=$(ls -d ${ANALYSIS_FOLDER}/QC/final_QC_output/*.1.unmerged.final.clean.fq | awk '{print $NF}')
   for i in ${finallist}
   do
      run_MaxBin.pl \
      -thread $THREADS \
      -contig ${ANALYSIS_FOLDER}/assembly/spades/$(basename ${i%.1*})/contigs_1000_filtered.fasta \
      -abund ${ANALYSIS_FOLDER}/coverage/maxbin.cov \
      -out ${ANALYSIS_FOLDER}/binning/maxbin/$(basename ${i%.1*})_maxbin
   done
   echo "DONE running Binning MaxBin"
}

binning_metabat(){
   echo "Running Binning Metabat"

   mkdir -p ${ANALYSIS_FOLDER}/binning/metabat
   finallist=$(ls -d ${ANALYSIS_FOLDER}/QC/final_QC_output/*.1.unmerged.final.clean.fq | awk '{print $NF}')
   for i in ${finallist}
   do
      metabat \
      -i ${ANALYSIS_FOLDER}/assembly/spades/$(basename ${i%.1*})/contigs.fasta \
      -a ${ANALYSIS_FOLDER}/coverage/depth.txt \
      -o ${ANALYSIS_FOLDER}/binning/metabat/$(basename ${i%.1*})_metabat \
      -t $THREADS \
      -m 1500 \
      -v \
      --unbinned
   done
   #In order to check which contigs were grouped together into separate bins by MetaBat example
   #grep ">" E01452_L001_to_L004_metabat.1.fa | sed 's/>//g' > E01452_L001_to_L004_metabat.1.contigNames

   echo "DONE running Binning Metabat!"
}

coverage_and_bining_main
