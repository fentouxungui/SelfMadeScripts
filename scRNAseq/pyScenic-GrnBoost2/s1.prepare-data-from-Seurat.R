library(Seurat)
cds <- readRDS('../subset-EE.rds')
cds
DimPlot(cds)
# expr_df <- as.data.frame(cds@assays$RNA@layers$data)
expr_df <- as.data.frame(cds@assays$RNA@data)
rownames(expr_df) <- rownames(cds)
colnames(expr_df) <- colnames(cds)
expr_df[1:5,1:5]
write.csv(t(expr_df),file = 'expression_data.csv', row.names = TRUE)

# 确保TF基因list存在于单细胞的基因名中
tf <- read.table('./mm_mgi_tfs.txt', stringsAsFactors = FALSE, col.names = 'gene')
head(tf)
table(tf$gene %in% rownames(expr_df))
library(BiologyDBLight)
head(AnimalTFDB_Mouse_TF)
table(AnimalTFDB_Mouse_TF$Symbol %in% rownames(expr_df))
table(rownames(expr_df) %in% c(tf$gene, AnimalTFDB_Mouse_TF$Symbol))
TFs <- rownames(expr_df)[rownames(expr_df) %in% c(tf$gene, AnimalTFDB_Mouse_TF$Symbol)]
TFs
write.table(TFs,file = 'mm_mgi_tfs_all.txt', quote = FALSE, col.names = FALSE, row.names = FALSE)
