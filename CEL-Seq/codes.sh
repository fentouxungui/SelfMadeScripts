STAR \
--runMode alignReads \
--genomeDir /home/xilab/Data_Backup/Data_From_Paper/2021.03.09_CellReports_PMID33691112_Human-gastrointestinal-epithelia-of-the-esophagus-stomach-and-duodenum-resolved-at-single-cell-resolution/data/reference_genome/star_index \
--readFilesIn /home/xilab/Data_Backup/Data_From_Paper/2021.03.09_CellReports_PMID33691112_Human-gastrointestinal-epithelia-of-the-esophagus-stomach-and-duodenum-resolved-at-single-cell-resolution/data/fastq/SRR12615656_2.fastq.gz /home/xilab/Data_Backup/Data_From_Paper/2021.03.09_CellReports_PMID33691112_Human-gastrointestinal-epithelia-of-the-esophagus-stomach-and-duodenum-resolved-at-single-cell-resolution/data/fastq/SRR12615656_1.fastq.gz \
--runThreadN 16 \
--soloType CB_UMI_Simple \
--soloCBstart 7 \
--soloCBlen 8 \
--soloUMIstart 1 \
--soloUMIlen 6 \
--soloBarcodeReadLength 0 \
--readFilesCommand zcat \
--soloCBwhitelist GSE157694_CELseq2_barcodes.txt \
--outFileNamePrefix output/SRR12615656

# use for seurat
# tree filtered/
# filtered/
# ├── barcodes.tsv
# ├── features.tsv
# └── matrix.mtx

