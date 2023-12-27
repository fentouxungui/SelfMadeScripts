library(monocle)
HSMM <- importCDS(pbmc, import_all = TRUE)

HSMM
HSMM <- detectGenes(HSMM, min_expr = 0.1)
HSMM <- estimateSizeFactors(HSMM)
HSMM <- estimateDispersions(HSMM)
expressed_genes <- row.names(subset(fData(HSMM),
                                    num_cells_expressed >= 10))
# Trajectory step 1: choose genes that define a cell's progress
pData(HSMM)$res.0.2 <- as.factor(pData(HSMM)$res.0.2)
diff_test_res <- differentialGeneTest(HSMM[expressed_genes,],
                                      fullModelFormulaStr = "~res.0.2")
ordering_genes <- row.names (subset(diff_test_res, qval < 0.01))
HSMM <- setOrderingFilter(HSMM, ordering_genes)
plot_ordering_genes(HSMM)
# Trajectory step 2: reduce data dimensionality
HSMM <- reduceDimension(HSMM, max_components = 2,
                            method = 'DDRTree')
# Trajectory step 3: order cells along the trajectory
HSMM_myo <- orderCells(HSMM)
plot_cell_trajectory(HSMM, color_by = "res.0.2")
plot_cell_trajectory(HSMM, color_by = "res.0.4")
plot_cell_trajectory(HSMM, color_by = "res.0.3")