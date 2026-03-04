library(Seurat)

cds <- readRDS("./Merge-Total-SeuratHarmony.rds")
DimPlot(cds,label = TRUE)
Embeddings(cds)


pbmc <- readRDS("scvi-tools/All_merged_scvi_tools.rds")
DimPlot(pbmc,label = TRUE)
Embeddings(pbmc,reduction = "umap")
Reductions(pbmc)

Reductions(cds)


## add new reductions
cds[["scvi_umap"]] <- pbmc@reductions$umap
DimPlot(cds, reduction = "scvi_umap", pt.size = 0.5)
cds[["scvi_tsne"]] <- pbmc@reductions$tsne
DimPlot(cds, reduction = "scvi_tsne", pt.size = 0.5)
cds[["scvi"]] <- pbmc@reductions$scvi
DimPlot(cds, reduction = "scvi", pt.size = 0.5)

## add new cluster infos
head(cds@meta.data)
head(pbmc@meta.data)

cds@meta.data[,paste("scvi",grep("^RNA",colnames(pbmc@meta.data),value = TRUE),sep = "_")] <- 
  pbmc@meta.data[,grep("^RNA",colnames(pbmc@meta.data),value = TRUE)]

DimPlot(cds, reduction = "scvi_umap", pt.size = 0.5,group.by = "scvi_RNA_snn_res.0.6")
DimPlot(pbmc,reduction = "umap",pt.size = 0.5,group.by = "RNA_snn_res.0.6")
DimPlot(cds, reduction = "scvi_tsne", pt.size = 0.5,group.by = "scvi_RNA_snn_res.0.6")
DimPlot(pbmc,reduction = "tsne",pt.size = 0.5,group.by = "RNA_snn_res.0.6")
Reductions(cds)

DimPlot(cds,label = TRUE)
cds <- ProjectDim(cds, reduction = "scvi",dims.print = 1:10)
DimHeatmap(cds, reduction = "scvi", dims = 1, cells = 500, projected = TRUE, balanced = TRUE)

saveRDS(cds,file = "Merge-Total-SeuratHarmony.rds")
