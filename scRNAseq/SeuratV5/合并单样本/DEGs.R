.libPaths(.libPaths()[2])
library(Seurat)
cds <- readRDS("./merged.rds")
cds <- JoinLayers(cds)
Idents(cds) <- "cca_clusters_res_0.6"
markers <- FindAllMarkers(cds,only.pos = TRUE)
write.csv(markers,file = "markers-of-cca_clusters_res_0.6.csv")
library(dplyr)
top_markers <- markers %>%
  group_by(cluster) %>%
  top_n(n = -20, p_val_adj)
pdf(file = "plot-markers.pdf", width = 9,height = 4)
for (i in unique(top_markers$gene)) {
  print(FeaturePlot(cds, features = i, split.by = "orig.ident", reduction = "umap.cca"))
}
dev.off()

cds$celltype <- Idents(cds)
Group.by <- "orig.ident"
group.levels <- c("ST72","ST_IL972")
cds@meta.data[,Group.by][is.na(cds@meta.data[,Group.by])] <- "ST_IL972"
Results.dir <-"DEGs"
Idents(cds) <- paste(cds@meta.data[,Group.by],Idents(cds), sep = "-")

cds@meta.data[,Group.by] <- factor(as.character(cds@meta.data[,Group.by]), levels = group.levels)

top_diff <- c()
if ( !dir.exists(Results.dir)) {
  dir.create(Results.dir, recursive = TRUE)
}

for (i in (1:length(levels(cds$celltype)))) {
  samplesSelected <- paste(levels(cds@meta.data[,Group.by]),i-1,sep = "-")
  message(paste("Using",samplesSelected[2],"as Case!",sep = " "))
  if (sum(as.character(Idents(cds)) %in% samplesSelected[1]) < 3 | sum(as.character(Idents(cds)) %in% samplesSelected[2]) < 3) {
    message(paste("Escaped cluster",i-1,"for very few cells!",sep = " "))
    next()
  }
  b.interferon.response <- suppressWarnings(suppressMessages(FindMarkers(cds, ident.1 = samplesSelected[2], ident.2 = samplesSelected[1], verbose = FALSE)) )
  write.csv(b.interferon.response, file = paste(Results.dir,"/DEGs-of-cluster-",i-1,"-",samplesSelected[2],"-vs-",samplesSelected[1],".csv",sep = ""))
  top_diff <- append(top_diff,rownames(b.interferon.response)[1:10])
}

Idents(cds) <- "celltype"

pdf(file = "plot-DEGs.pdf", width = 9,height = 4)
for ( gene in unique(top_diff)) {
  print(FeaturePlot(cds, features = gene, split.by = Group.by, pt.size = 0.5, cols = c("grey", "red"), reduction = "umap.cca"))
}
dev.off()
