#!/usr/bin/evn bash


# conda install -c bioconda idba
THREADS=24

### Assembly Analysis
assembly_main(){
   # check_and_install
   conda activate metagenomic
   run_assembly
   assembly_stats
   conda activate py37
   filter_on_assembly_length
}

# Run assembly 
run_assembly(){
   echo "Running all assemblies"
   echo "Running assembly"
   mkdir -p ${ANALYSIS_FOLDER}/assembly
   run_assem_megahit
   run_assem_spades
   # run_assem_ray
   run_assem_idba
   echo "DONE running all assemblies!"
}

run_assem_idba(){
   # Assembly with IDBA
   # first covert paired into single
   # idba_ud -b y /usr/bin/idba_ud -r reads12.fas --num_threads 14 -o ${ANALYSIS_FOLDER}/assembly/idba_ud_out/${i%.1*}
   # 1. Merged the paired end data and covert to fasta format.
   # 注意事项： https://github.com/loneknightpy/idba/issues/26
   # 注意reads的长度，-r适用于<128bp的reads，-l适用于大于128bp的reads
   # seqkit stats sample.12.final.clean.fa # check your reads length
   # 当reads长度大于128bp时，用-l选项。我还没有尝试用-l选项来处理长度大于128bp的reads！
   # 注意--maxk不要大于reads的长度
   # 有人说：Even more strangely, I get a segfault with 12 cores, no segfault with 16 cores, segfault again at 24 cores, 
   # no segfault at 3 cores. Not sure what's going on here. I can also confirm this on real datasets and that it does not segfault at 40 cores
   # from <https://groups.google.com/g/hku-idba/c/RzTkrVTod8o/m/UWAPxhKQKSEJ>
   # 运行时间很长~，建议过夜跑
   # 如果既有大于128bp的reads，也有小于的，是不是分成两个文件？
   echo "Running idba/fq2fa"
   mkdir -p ${ANALYSIS_FOLDER}/assembly/idba_ud
   finallist=$(ls -d ${ANALYSIS_FOLDER}/QC/bbmap/*.1.unmerged.final.clean.fq | awk '{print $NF}')
   for i in ${finallist}
   do
      fq2fa --merge ${i} ${i%1*}2.unmerged.final.clean.fq ${i%1*}12.final.clean.fa

      idba_ud \
      # -l ${i%1*}12.final.clean.fa \
      -r ${i%1*}12.final.clean.fa \
      --mink 20 \
      --maxk 120 \
      --num_threads 16 \
      -o ${ANALYSIS_FOLDER}/assembly/idba_ud/$(basename ${i%.1*})
   done
}

# run_assem_ray(){
#    mkdir -p ${ANALYSIS_FOLDER}/assembly/ray
#    #Assembly with Ray
#    echo "Running ray"
#    finallist=$(ls -d ${ANALYSIS_FOLDER}/QC/final_QC_output/*.1.unmerged.final.clean.fq | awk '{print $NF}')
#    for i in ${finallist}
#    do
#       mpiexec -n 16 \
#       ray \
#       -k 31 \
#       -p ${i} ${i%1*}2.unmerged.final.clean.fq \
#       -s ${i%1*}u.final.clean.fq  \
#       -o ${ANALYSIS_FOLDER}/assembly/ray/ray_31_assembly/$(basename ${i%.1*})
#    done
# }

run_assem_spades(){
   mkdir -p ${ANALYSIS_FOLDER}/assembly/spades
   # Assembly with spades
   finallist=$(ls -d ${ANALYSIS_FOLDER}/QC/final_QC_output/*.1.unmerged.final.clean.fq | awk '{print $NF}')
   echo "Running spades"
   for i in ${finallist}
   do
      spades.py \
      --meta \
      -1 ${i} \
      -2 ${i%1*}2.unmerged.final.clean.fq \
      -s ${i%1*}u.final.clean.fq \
      -k auto \
      -o ${ANALYSIS_FOLDER}/assembly/spades/$(basename ${i%.1*}) \
      -t $THREADS
   done
}

run_assem_megahit(){
   mkdir -p ${ANALYSIS_FOLDER}/assembly/megahit
   # Assembly Megahit
   finallist=$(ls -d ${ANALYSIS_FOLDER}/QC/final_QC_output/*.1.unmerged.final.clean.fq | awk '{print $NF}')
   #running megahit
   echo "Running megahit"
   for i in ${finallist}
   do
     #To run a co-assembly using multiple samples, you group together the R1 and R2 reads as shown for three libraries (six read files):
     #megahit [options] {-1 pe1_R1.fq,pe2_R1.fq,pe3_R1.fq -2 pe1_R2.fq,pe2_R2.fq,pe3_R2.fq} -o outdir
      megahit \
      -1 ${i} \
      -2 ${i%1*}2.unmerged.final.clean.fq \
      -r ${i%1*}u.final.clean.fq \
      -t $THREADS \
      --presets meta-large \
      --mem-flag 2 \
      -o ${ANALYSIS_FOLDER}/assembly/megahit/$(basename ${i%.1*})
   done

   cat ${ANALYSIS_FOLDER}/assembly/megahit/$(basename ${i%.1*})/final.contigs.fa | cut -d ' ' -f 1 > ${ANALYSIS_FOLDER}/assembly/megahit/$(basename ${i%.1*})/megahit_final.contigs.fa
}

assembly_stats(){
   echo "Running stats"
   finallist=$(ls -d ${ANALYSIS_FOLDER}/QC/final_QC_output/*.1.unmerged.final.clean.fq | awk '{print $NF}')
   # creating stats of all the assemblers
   # Stats with quast
   # http://quast.sourceforge.net/docs/manual.html
   # If you run MetaQUAST without providing reference genomes, the tool will try to identify genome content of the metagenome. 
   # MetaQUAST uses BLASTN for aligning contigs to SILVA 16S rRNA database, i.e. FASTA file containing small subunit ribosomal RNA sequences. 
   # For each assembly, 50 reference genomes with top scores are chosen. Maximum number of references to download can be specified with --max-ref-number.
   mkdir -p ${ANALYSIS_FOLDER}/assembly/assembly_stats
   for i in ${finallist}
   do
      metaquast.py \
      --threads $THREADS \
      --gene-finding \
      -o $ANALYSIS_FOLDER/assembly/assembly_stats/$(basename ${i%.1*}) \
      -l megahit,SPAdes,IDBA-UD \
      $ANALYSIS_FOLDER/assembly/megahit/$(basename ${i%.1*})/megahit_final.contigs.fa \
      $ANALYSIS_FOLDER/assembly/spades/$(basename ${i%.1*})/contigs.fasta \
      ${ANALYSIS_FOLDER}/assembly/idba_ud/$(basename ${i%.1*})/contig.fa
      #${ANALYSIS_FOLDER}/assembly/ray/$(basename ${i%.1*})/ray_31_assembly/  \
   done
   echo "DONE running stats"
}



filter_on_assembly_length(){
   conda activate py37
   echo "Running assembly length filter"
   finallist=$(ls -d ${ANALYSIS_FOLDER}/QC/final_QC_output/*.1.unmerged.final.clean.fq | awk '{print $NF}')
   #running megahit
   for i in ${finallist}
   do
      #spades
      reformat.sh \
      in=${ANALYSIS_FOLDER}/assembly/spades/$(basename ${i%.1*})/contigs.fasta \
      out=${ANALYSIS_FOLDER}/assembly/spades/$(basename ${i%.1*})/contigs_1000_filtered.fasta \
      minlength=1000

      #megahit
      reformat.sh \
      in=${ANALYSIS_FOLDER}/assembly/megahit/$(basename ${i%.1*})/megahit_final.contigs.fa \
      out=${ANALYSIS_FOLDER}/assembly/megahit/$(basename ${i%.1*})/megahit_final_1000_filtered.fasta \
      minlength=1000

      #idba_ud
      reformat.sh \
      in=${ANALYSIS_FOLDER}/assembly/idba_ud/$(basename ${i%.1*})/contig.fa \
      out=${ANALYSIS_FOLDER}/assembly/idba_ud/$(basename ${i%.1*})/contig_1000_filtered.fa \
      minlength=1000
   done
   echo "DONE running assembly length filter"
}

assembly_main
