library(Seurat)

ee <- Read10X(data.dir = "../ee_filtered_gene_bc_matrices/BDGP6/")
ee <- CreateSeuratObject(raw.data = ee, project = "ee_data", min.cells = 5)
ee@meta.data$status <- "ee"
ee <- FilterCells(ee, subset.names = "nGene", low.thresholds = 500, high.thresholds = Inf)
ee <- NormalizeData(ee)
ee <- ScaleData(ee, display.progress = F)
# Set up stimulated object
esg <- Read10X(data.dir = "../esg_filtered_gene_bc_matrices/BDGP6/")
esg <- CreateSeuratObject(raw.data = esg, project = "esg_data", min.cells = 5)
esg@meta.data$status <- "esg"
esg <- FilterCells(esg, subset.names = "nGene", low.thresholds = 500, high.thresholds = Inf)
esg <- NormalizeData(esg)
esg <- ScaleData(esg, display.progress = F)

# Gene selection for input to CCA
ee <- FindVariableGenes(ee, do.plot = F)
esg <- FindVariableGenes(esg, do.plot = F)
g.1 <- head(rownames(ee@hvg.info), 1000)
g.2 <- head(rownames(esg@hvg.info), 1000)
genes.use <- unique(c(g.1, g.2))
genes.use <- intersect(genes.use, rownames(ee@scale.data))
genes.use <- intersect(genes.use, rownames(esg@scale.data))
head(genes.use)
length(genes.use)

# Perform a canonical correlation analysis (CCA)
combined <- RunCCA(ee, esg, genes.use = genes.use, num.cc = 30,add.cell.id1 = "ee",add.cell.id2 = "esg")
# visualize results of CCA plot CC1 versus CC2 and look at a violin plot
p1 <- DimPlot(object = combined, reduction.use = "cca", group.by = "status", 
              pt.size = 0.5, do.return = TRUE)
p2 <- VlnPlot(object = combined, features.plot = "CC1", group.by = "status", 
              do.return = TRUE)
plot_grid(p1, p2)


PrintDim(object = combined, reduction.type = "cca", dims.print = 1:2, 
         genes.print = 10)

p3 <- MetageneBicorPlot(combined, grouping.var = "status", dims.eval = 1:30, 
                        display.progress = FALSE)

DimHeatmap(object = combined, reduction.type = "cca", cells.use = 500, 
           dim.use = 1:9, do.balanced = TRUE)

# Align the CCA subspaces
immune.combined <- AlignSubspace(combined, reduction.type = "cca", grouping.var = "status", 
                                 dims.align = 1:20)

p1 <- VlnPlot(object = immune.combined, features.plot = "ACC1", group.by = "status", 
              do.return = TRUE)
p2 <- VlnPlot(object = immune.combined, features.plot = "ACC2", group.by = "status", 
              do.return = TRUE)
plot_grid(p1, p2)

# Perform an integrated analysis
# t-SNE and Clustering
immune.combined <- RunTSNE(immune.combined, reduction.use = "cca.aligned", dims.use = 1:20, 
                           do.fast = T)
immune.combined <- FindClusters(immune.combined, reduction.type = "cca.aligned", 
                                resolution = 0.6, dims.use = 1:20)

# Visualization
p1 <- TSNEPlot(immune.combined, do.return = T, pt.size = 0.5, group.by = "status")
p2 <- TSNEPlot(immune.combined, do.label = T, do.return = T, pt.size = 0.5)
plot_grid(p1, p2)








