require("Matrix")
# import ee and esg data, construct SingleCellExperiment Object
cellbarcodes <- read.table("../ee_filtered_gene_bc_matrices/BDGP6/barcodes.tsv")
genenames <- read.table("../ee_filtered_gene_bc_matrices/BDGP6/genes.tsv")
molecules <- Matrix::readMM("../ee_filtered_gene_bc_matrices/BDGP6/matrix.mtx")

head(cellbarcodes)
head(genenames)
colnames(molecules) <- paste("ee", cellbarcodes[,1], sep="_")
head(molecules)

ee <- molecules

cellbarcodes <- read.table("../esg_filtered_gene_bc_matrices/BDGP6/barcodes.tsv")
genenames <- read.table("../esg_filtered_gene_bc_matrices/BDGP6/genes.tsv")
molecules <- Matrix::readMM("../esg_filtered_gene_bc_matrices/BDGP6/matrix.mtx")

head(cellbarcodes)
head(genenames)
colnames(molecules) <- paste("esg", cellbarcodes[,1], sep="_")
head(molecules)

esg <- molecules

library(monocle)
HSMM_expr_matrix <- cbind(ee,esg)
head(HSMM_expr_matrix)
length(colnames(ee))
length(colnames(esg))
HSMM_sample_sheet <- data.frame(type = c(rep("ee",4769),rep("esg",3121)))



rownames(HSMM_sample_sheet) <- colnames(HSMM_expr_matrix)
head(HSMM_sample_sheet)
HSMM_gene_annotation <- data.frame(gene_short_name = genenames$V2)
rownames(HSMM_gene_annotation) <- genenames$V1


pd <- new("AnnotatedDataFrame", data = HSMM_sample_sheet)
fd <- new("AnnotatedDataFrame", data = HSMM_gene_annotation)
HSMM <- newCellDataSet(as.matrix(HSMM_expr_matrix),
                       phenoData = pd, featureData = fd,lowerDetectionLimit = 0.5)
############ expressionFamily not checked!!!
HSMM <- estimateSizeFactors(HSMM)
HSMM <- estimateDispersions(HSMM)
HSMM <- detectGenes(HSMM, min_expr = 0.1)
head(pData(HSMM))
expressed_genes <- row.names(subset(fData(HSMM),
                                    num_cells_expressed >= 10))
pData(HSMM)$Total_mRNAs <- Matrix::colSums(exprs(HSMM))
upper_bound <- 10^(mean(log10(pData(HSMM)$Total_mRNAs)) +
                     2*sd(log10(pData(HSMM)$Total_mRNAs)))
lower_bound <- 10^(mean(log10(pData(HSMM)$Total_mRNAs)) -
                     2*sd(log10(pData(HSMM)$Total_mRNAs)))
qplot(Total_mRNAs, data = pData(HSMM), color = type, geom =
        "density") +
  geom_vline(xintercept = lower_bound) +
  geom_vline(xintercept = upper_bound)
head(pData(HSMM))


HSMM <- HSMM[,pData(HSMM)$Total_mRNAs > lower_bound &
               pData(HSMM)$Total_mRNAs < upper_bound]
HSMM <- detectGenes(HSMM, min_expr = 0.1)

# Log-transform each value in the expression matrix.
L <- log(exprs(HSMM[expressed_genes,]))

library(reshape2)
# Standardize each gene, so that they are all on the same scale,
# Then melt the data with plyr so we can plot it easily
melted_dens_df <- melt(Matrix::t(scale(Matrix::t(L))))

# Plot the distribution of the standardized gene expression values.
qplot(value, geom = "density", data = melted_dens_df) +
  stat_function(fun = dnorm, size = 0.5, color = 'red') +
  xlab("Standardized log(FPKM)") +
  ylab("Density")
## failed! Warning messages:
#1: Removed 70796970 rows containing non-finite values (stat_density). 
#2: Computation failed in `stat_function()`:
#  'from' must be a finite numbe

# Clustering cells without marker genes 
disp_table <- dispersionTable(HSMM)
unsup_clustering_genes <- subset(disp_table, mean_expression >= 0.1)
HSMM <- setOrderingFilter(HSMM, unsup_clustering_genes$gene_id)
plot_ordering_genes(HSMM)

# HSMM@auxClusteringData[["tSNE"]]$variance_explained <- NULL
plot_pc_variance_explained(HSMM, return_all = F) # norm_method='log'

HSMM <- reduceDimension(HSMM, max_components = 2, num_dim = 6,
                        reduction_method = 'tSNE', verbose = T)
HSMM <- clusterCells(HSMM, num_clusters = 2)
plot_cell_clusters(HSMM, 1, 2, color = "type",
                   markers = c("NPF", "N"))

plot_cell_clusters(HSMM, 1, 2, color = "type")


HSMM <- reduceDimension(HSMM, max_components = 2, num_dim = 2,
                        reduction_method = 'tSNE',
                        residualModelFormulaStr = "~type + num_genes_expressed",
                        verbose = T)
HSMM <- clusterCells(HSMM, num_clusters = 2)
plot_cell_clusters(HSMM, 1, 2, color = "type")
plot_cell_clusters(HSMM, 1, 2, color = "Cluster") +
  facet_wrap(~type)

# Constructing Single Cell Trajectories
# Trajectory step 1: choose genes that define a cell's progress
diff_test_res <- differentialGeneTest(HSMM[expressed_genes,],
                                      fullModelFormulaStr = "~type")
ordering_genes <- row.names (subset(diff_test_res, qval < 0.01))

HSMM_myo <- setOrderingFilter(HSMM, ordering_genes)
plot_ordering_genes(HSMM)

HSMM_myo <- reduceDimension(HSMM, max_components = 2,
                            method = 'DDRTree')
HSMM_myo <- orderCells(HSMM_myo)

plot_cell_trajectory(HSMM_myo, color_by = "State")
plot_cell_trajectory(HSMM_myo, color_by = "type")






























