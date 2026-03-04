library(DWLS)
library(Seurat)

file_path <- "/data3/Xilab-Data-Analysis/paper-data-analysis/2017.11.08_Nature_PMID29144463_Mouse-Adult-Small-Intestine-EpCAM-Epithelium-Atalas-scRNAseq/scRNAseq/velocyto-notebooks-Demo/loom-to-Seurat/haber.epithelium.from.loom.rds"
cds <- readRDS(file = file_path)
cds <- UpdateSeuratObject(cds)
cds

T <- read.delim("~/guoxt/RNAseq/2024.12.30_RNAseq_PE_Mouse-intestine-epithelial-Irx3-cko-homozygous-vs-heterozygous-mutation/results/6_Counts_featureCounts/final_counts.txt",stringsAsFactors = FALSE, skip = 1)
rownames(T) <- T$Geneid
T <- T[,-c(1:6)]
head(T)
colnames(T) <- gsub(".Aligned.sortedByCoord.out.bam", "", colnames(T), fixed = TRUE)
colnames(T) <- gsub("..", "", colnames(T), fixed = TRUE)

counts_mat <- as.matrix(cds@assays$RNA@counts)
counts_mat <- counts_mat[rownames(counts_mat) %in% rownames(T),]

labels <- as.character(cds@meta.data$Maintype)

Signature <- buildSignatureMatrixMAST(counts_mat,labels, path= "results",diff.cutoff = 0.5,pval.cutoff = 0.01)
Signature = as.matrix(Signature)



RESULTS <- apply(T, 2, function(x){
  b = setNames(x, rownames(T))
  tr <- DWLS::trimData(Signature, b)
  RES <- t(DWLS::solveDampenedWLS(tr$sig, tr$bulk))
})

rownames(RESULTS) <- as.character(unique(colnames(Signature)))
RESULTS = apply(RESULTS,2,function(x) ifelse(x < 0, 0, x)) #explicit non-negativity constraint
RESULTS = apply(RESULTS,2,function(x) x/sum(x)) #explicit STO constraint
print(head(RESULTS))

write.csv(RESULTS,file = "cell-type-percentage.csv")



