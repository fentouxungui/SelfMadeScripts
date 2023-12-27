library(Seurat)

cds <- readRDS("../../Subset-Endothelial-seuratHarmony.rds")
cds <- subset(cds,subset = type == "H")
table(cds@meta.data$type)

all(rownames(cds@meta.data) == Cells(cds))

cells <- paste(cds@meta.data$orig.ident,gsub("_.*","x",Cells(cds)),sep = ":")
write.csv(cells, file = "cellID_obs.csv")

embedding.umap <- Embeddings(cds, reduction = "umap")
all(rownames(embedding.umap) == Cells(cds))
rownames(embedding.umap) <- cells
write.csv(embedding.umap, file = "cell_embeddings.csv")

cluster <- cds@meta.data
cluster$cells <- cells
write.csv(cluster, file = "clusters.csv")
