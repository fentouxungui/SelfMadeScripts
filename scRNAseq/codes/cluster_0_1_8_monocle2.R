HSMM <- HSMM[,pData(HSMM)$res.0.2 %in% c("0","1","8")]
HSMM
head(pData(HSMM)$res.0.2)
table(pData(HSMM)$res.0.2)
#set up the color for each experiment
cols <- c("0" = "#edf8fb", "1" = "#ccece6", "8" = "#99d8c9")
HSMM <- estimateSizeFactors(HSMM)
HSMM <- estimateDispersions(HSMM)
#Running dpFeature for selecting ordering gene
#1. set ordering genes
ordering_genes <- c("N","Dl","AstC","TK","NPF","NPFR")
HSMM <- setOrderingFilter(HSMM, ordering_genes = ordering_genes)
URMM_pc_variance <- plot_pc_variance_explained(HSMM, return_all = T, norm_method = 'log')
#2. run reduceDimension with tSNE as the reduction_method
set.seed(2017)
URMM_all_fig1b <- reduceDimension(HSMM, max_components=2, norm_method = 'log', reduction_method = 'tSNE', num_dim = 12,  verbose = F, check_duplicates = FALSE)

#3. initial run of clusterCells
URMM_all_fig1b <- clusterCells(URMM_all_fig1b, verbose = F)

#4. check the clusters
options(repr.plot.width=4, repr.plot.height=3)
plot_cell_clusters(URMM_all_fig1b, color_by = 'as.factor(Cluster)') + theme (legend.position="left", legend.title=element_blank())# show_density = F,
plot_cell_clusters(URMM_all_fig1b, color_by = 'res.0.4') + theme (legend.position="left", legend.title=element_blank())
plot_cell_clusters(URMM_all_fig1b, color_by = 'res.0.2') + theme (legend.position="left", legend.title=element_blank())
plot_rho_delta(URMM_all_fig1b) 

URMM_all_fig1b@expressionFamily <- negbinomial.size()
pData(URMM_all_fig1b)$Cluster <- factor(pData(URMM_all_fig1b)$Cluster)
URMM_clustering_DEG_genes <- differentialGeneTest(URMM_all_fig1b, fullModelFormulaStr = '~Cluster', cores = detectCores() - 2)

#use all DEG gene from the clusters
URMM_ordering_genes <- row.names(URMM_clustering_DEG_genes)[order(URMM_clustering_DEG_genes$qval)][1:1000]

# Reconstruct the developmental trajectory for the wild-type data
URMM_all_fig1b <- setOrderingFilter(URMM_all_fig1b, ordering_genes = c(URMM_ordering_genes))
URMM_all_fig1b <- reduceDimension(URMM_all_fig1b, verbose = F, scaling = T, max_components = 4, maxIter = 100, norm_method = 'log',  lambda = 20 * ncol(URMM_all_fig1b)) 
URMM_all_fig1b <- orderCells(URMM_all_fig1b)
options(repr.plot.width=3, repr.plot.height=3)
plot_cell_trajectory(URMM_all_fig1b, color_by = 'res.0.2')
options(repr.plot.width=8, repr.plot.height=8)
plot_cell_trajectory(URMM_all_fig1b, color_by = 'Cluster') + facet_wrap(~Cluster)
options(repr.plot.width=8, repr.plot.height=8)
plot_cell_trajectory(URMM_all_fig1b, color_by = 'Cluster', x = 1, y = 3) + facet_wrap(~Cluster)

options(repr.plot.width=8, repr.plot.height=8)
plot_cell_trajectory(URMM_all_fig1b, color_by = 'res.0.2', x = 1, y = 3) + facet_wrap(~Cluster)

