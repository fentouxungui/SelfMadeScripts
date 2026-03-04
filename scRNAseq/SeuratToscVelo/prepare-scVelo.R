library(Seurat)
cds <- readRDS("./jiang-pdgfra-from-adult-SI-Crypt.rds")

Results.dir <- "scVelo/"
if ( !dir.exists(Results.dir)) {
  dir.create(Results.dir)
}
cells <- Cells(cds)
cells <- gsub("-\\d+$","",cells)
write.csv(cells, file = paste(Results.dir,"cellID_obs.csv",sep = ""))

embedding.umap <- Embeddings(cds, reduction = "umap")
all(rownames(embedding.umap) == Cells(cds))
rownames(embedding.umap) <- cells
write.csv(embedding.umap, file = paste(Results.dir,"cell_embeddings_umap.csv",sep = ""))

embedding.umap <- Embeddings(cds, reduction = "tsne")
all(rownames(embedding.umap) == Cells(cds))
rownames(embedding.umap) <- cells
write.csv(embedding.umap, file = paste(Results.dir,"cell_embeddings_tsne.csv",sep = ""))

cluster <- cds@meta.data
cluster$cells <- cells
write.csv(cluster, file = paste(Results.dir,"clusters.csv",sep = ""))
