# seurat 2.3.4版本
library(Seurat)
library(dplyr)
# Load the PBMC dataset
pbmc.data <- Read10X(data.dir = "hg19/")
pbmc <- CreateSeuratObject(raw.data = pbmc.data,
                           min.cells = 3, min.genes = 200, 
                           project = "10X_PBMC")

library(monocle)
ee_wild <- importCDS(pbmc.data)
# Estimate size factor and dispersion
DelayedArray:::set_verbose_block_processing(TRUE)
options(DelayedArray.block.size=1000e6)
ee_wild <- estimateSizeFactors(ee_wild)
ee_wild <- estimateDispersions(ee_wild)
# Reduce the dimensionality of ee_wild dataset
library(dplyr)

disp_table = dispersionTable(ee_wild)
disp_table = disp_table %>% mutate(excess_disp =
                                     (dispersion_empirical - dispersion_fit) / dispersion_fit) %>%
  arrange(plyr::desc(excess_disp))
top_subset_genes = as.character(head(disp_table, 2500)$gene_id)

ee_wild = setOrderingFilter(ee_wild, top_subset_genes)
ee_wild <- preprocessCDS(ee_wild,  method = 'PCA',
                         norm_method = 'log',
                         num_dim = 50,
                         verbose = T)
ee_wild <- reduceDimension(ee_wild, max_components = 2,
                           reduction_method = 'UMAP',
                           metric="correlation",
                           min_dist = 0.75,
                           n_neighbors = 50,
                           verbose = T)

# Group cells into different clusters
ee_wild <- clusterCells(ee_wild,
                        method = 'louvain',
                        res = 1e-6,
                        louvain_iter = 1,
                        verbose = T)
    # -Number of clusters: 3 

# Visualize the clustering results
col_vector_origin <- c("#db83da",
                       "#53c35d",
                       "#a546bb",
                       "#83b837",
                       "#a469e6",
                       "#babb3d",
                       "#4f66dc",
                       "#e68821",
                       "#718fe8",
                       "#d6ac3e",
                       "#7957b4",
                       "#468e36",
                       "#d347ae",
                       "#5dbf8c",
                       "#e53e76",
                       "#42c9b8",
                       "#dd454a",
                       "#3bbac6",
                       "#d5542c",
                       "#59aadc",
                       "#cf8b36",
                       "#4a61b0",
                       "#8b8927",
                       "#a24e99",
                       "#9cb36a",
                       "#ca3e87",
                       "#36815b",
                       "#b23c4e",
                       "#5c702c",
                       "#b79add",
                       "#a55620",
                       "#5076af",
                       "#e38f67",
                       "#85609c",
                       "#caa569",
                       "#9b466c",
                       "#88692c",
                       "#dd81a9",
                       "#a35545",
                       "#e08083",
                       "#17becf",
                       "#9edae5")


#col_vector <-
#  col_vector_origin[1:length(unique(as.character(pData(ee_wild)$orig.ident)))]
#names(col_vector) <- unique(as.character(pData(ee_wild)$orig.ident))
#options(repr.plot.width = 11)
#options(repr.plot.height = 8)
#plot_cell_clusters(ee_wild,
#                   color_by = 'orig.ident',
#                   cell_size = 0.1,
#                   show_group_id = T)  +
#  theme(legend.text=element_text(size=6)) + #set the size of the text
#  theme(legend.position="right") #put the color legend on the right


options(repr.plot.width = 11)
options(repr.plot.height = 8)
cluster_col_vector <- col_vector_origin[1:length(unique(as.character(pData(ee_wild)$Cluster)))]
names(cluster_col_vector) <- unique(as.character(pData(ee_wild)$Cluster))
plot_cell_clusters(ee_wild,
                   color_by = 'Cluster',
                   cell_size = 0.5,
                   show_group_id = T) +
  scale_color_manual(values = c("red","green","blue","yellow")) +
  theme(legend.text=element_text(size=15)) + #set the size of the text
  theme(legend.position="top") #put the color legend on the right


# Identify genes that are differentially expressed between clusters
start <- Sys.time()
spatial_res <- principalGraphTest(ee_wild, relative_expr = TRUE, k = 3, cores = detectCores() - 2, verbose = FALSE)
end <- Sys.time()
end - start

# Finding cluster-specific marker genes
cluster_marker_res <-
  find_cluster_markers(ee_wild,
                       spatial_res,
                       group_by = 'Cluster',
                       morans_I_threshold = 0.25)

genes <- (cluster_marker_res %>%
            dplyr::filter(mean > 0.5, percentage > 0.1) %>%
            dplyr::group_by(Group) %>% dplyr::slice(which.max(specificity)))
options(repr.plot.width=22, repr.plot.height=12)
plot_markers_by_group(ee_wild, genes$gene_short_name, group_by = 'Cluster', ordering_type = 'maximal_on_diag')

write.csv(genes,"cluster-specific_marker_genes.csv")
# Visualizing marker expression across cell clusters
genes <- (cluster_marker_res %>%
     dplyr::filter(mean > 0.5, percentage > 0.1) %>%
     dplyr::group_by(Group) %>%
     dplyr::top_n(50, wt = specificity))
#pdf("heatmap_top50_markers_by_specificity.pdf",height = 25,width = 8)
plot_markers_cluster(ee_wild, as.character(genes$gene_short_name),
                     minimal_cluster_fraction = 0.01)

#dev.off()
  #按照specificity > 0.7
options(repr.plot.width=22, repr.plot.height=12)
genes <- cluster_marker_res %>%
  dplyr::filter(mean > 0.5, percentage > 0.1, specificity > 0.7) %>%
  dplyr::group_by(Group) %>%
  dplyr::arrange(Group, dplyr::desc(specificity))
pdf("heatmap_top50_markers_filtered_specificity_.pdf",height = 25,width = 8)
plot_markers_cluster(ee_wild,
                     as.character(genes$gene_short_name),
                     minimal_cluster_fraction = 0.01,
                     show_rownames = T)
dev.off()
# Visualizing marker expression in UMAP plots
top_4_genes <- genes %>%
  dplyr::filter(Group == 1) %>%
  dplyr::top_n(n = 4, wt  = specificity)
top_4_genes

options(repr.plot.width = 8)
options(repr.plot.height = 9)
marker_genes <- top_4_genes$Gene
pdf("top4_of_Group1_UMAP_plots.pdf",height = 10,width = 8)
plot_cell_clusters(ee_wild,
                   markers = as.character(marker_genes),
                   show_group_id = T, cell_size = 0.5)
dev.off()


# 描述各个基因在UMAP上的分布
cluster_marker_res %>%
  dplyr::filter(gene_short_name %in% c('CG42826', 'CG13229'), Group == 7)
options(repr.plot.width = 8)
options(repr.plot.height = 5.5)
plot_cell_clusters(ee_wild,
                   markers = c('CG42826','CG13229'),
                   show_group_id = T, cell_size = 0.5)


# 构建trajectories
# 个人感觉，应该多分出几个群，然后再做此分析，3个群不好解释，
# 或者把最大的那个群拿出来，再做分析
# 尝试一：3个群直接做
#Step 1: Noramlize and pre-process the data
DelayedArray:::set_verbose_block_processing(TRUE)

# Passing a higher value will make some computations faster but use more memory. Adjust with caution!
options(DelayedArray.block.size=1000e6)
#ee_wild <- estimateSizeFactors(ee_wild)
#ee_wild <- estimateDispersions(ee_wild)
ee_wild <- preprocessCDS(ee_wild, num_dim = 20)
# Reduce the dimensionality of the MCA dataset
ee_wild <- reduceDimension(ee_wild, reduction_method = 'UMAP')
# Step 3: Partition the cells into supergroups
ee_wild <- partitionCells(ee_wild)
# Step 4: Learn the principal graph
ee_wild <- learnGraph(ee_wild,  RGE_method = 'SimplePPT')
# Step 5: Visualize the trajectory
plot_cell_trajectory(ee_wild,
                     color_by = "Cluster") +
  scale_color_manual(values = c("red","green","blue"))
######## bad results


# Try split to more clusters

ee_wild2 <- clusterCells(ee_wild,
                         method = 'louvain',
                         res = 1e-8,
                         louvain_iter = 1,
                         verbose = T)
        ### -Number of clusters: 12
options(repr.plot.width = 11)
options(repr.plot.height = 8)
cluster_col_vector <- col_vector_origin[1:length(unique(as.character(pData(ee_wild2)$Cluster)))]
names(cluster_col_vector) <- unique(as.character(pData(ee_wild2)$Cluster))
plot_cell_clusters(ee_wild2,
                   color_by = 'Cluster',
                   cell_size = 0.5,
                   show_group_id = T) +
  scale_color_manual(values = cluster_col_vector) +
  theme(legend.text=element_text(size=10)) + #set the size of the text
  theme(legend.position="top") #put the color legend on the right
  

# Identify genes that are differentially expressed between 12 clusters
start <- Sys.time()
spatial_res <- principalGraphTest(ee_wild2, relative_expr = TRUE, k = 12, cores = detectCores() - 2, verbose = FALSE)
end <- Sys.time()
end - start

# Finding cluster-specific marker genes
cluster_marker_res <-
  find_cluster_markers(ee_wild2,
                       spatial_res,
                       group_by = 'Cluster',
                       morans_I_threshold = 0.25)

genes <- (cluster_marker_res %>%
            dplyr::filter(mean > 0.5, percentage > 0.1) %>%
            dplyr::group_by(Group) %>% dplyr::slice(which.max(specificity)))
options(repr.plot.width=22, repr.plot.height=12)
plot_markers_by_group(ee_wild2, genes$gene_short_name, group_by = 'Cluster', ordering_type = 'maximal_on_diag')

write.csv(genes,"cluster-specific_marker_genes.csv")
# Visualizing marker expression across cell clusters
genes <- (cluster_marker_res %>%
            dplyr::filter(mean > 0.5, percentage > 0.1) %>%
            dplyr::group_by(Group) %>%
            dplyr::top_n(50, wt = specificity))
#pdf("heatmap_top50_markers_by_specificity.pdf",height = 25,width = 8)
plot_markers_cluster(ee_wild2, as.character(genes$gene_short_name),
                     minimal_cluster_fraction = 0.01)

#dev.off()
#按照specificity > 0.7
options(repr.plot.width=22, repr.plot.height=12)
genes <- cluster_marker_res %>%
  dplyr::filter(mean > 0.5, percentage > 0.1, specificity > 0.7) %>%
  dplyr::group_by(Group) %>%
  dplyr::arrange(Group, dplyr::desc(specificity))
pdf("heatmap_top50_markers_filtered_specificity_.pdf",height = 25,width = 8)
plot_markers_cluster(ee_wild2,
                     as.character(genes$gene_short_name),
                     minimal_cluster_fraction = 0.01,
                     show_rownames = T)
dev.off()
# Visualizing marker expression in UMAP plots
top_4_genes <- genes %>%
  dplyr::filter(Group == 8) %>%
  dplyr::top_n(n = 4, wt  = specificity)
top_4_genes

options(repr.plot.width = 8)
options(repr.plot.height = 9)
marker_genes <- top_4_genes$Gene
pdf("top4_of_Group8_UMAP_plots.pdf",height = 6,width = 8)
plot_cell_clusters(ee_wild2,
                   markers = as.character(marker_genes),
                   show_group_id = T, cell_size = 0.5)
dev.off()


# 描述各个基因在UMAP上的分布
cluster_marker_res %>%
  dplyr::filter(gene_short_name %in% c('CG42826', 'CG13229'), Group == 7)
options(repr.plot.width = 8)
options(repr.plot.height = 5.5)
plot_cell_clusters(ee_wild2,
                   markers = c('CG42826','CG13229'),
                   show_group_id = T, cell_size = 0.5)



# 构建trajectories
#Step 1: Noramlize and pre-process the data
DelayedArray:::set_verbose_block_processing(TRUE)
# Passing a higher value will make some computations faster but use more memory. Adjust with caution!
options(DelayedArray.block.size=1000e6)
#ee_wild <- estimateSizeFactors(ee_wild)
#ee_wild <- estimateDispersions(ee_wild)
ee_wild2 <- preprocessCDS(ee_wild2, num_dim = 20)
# Reduce the dimensionality of the MCA dataset
ee_wild2 <- reduceDimension(ee_wild2, reduction_method = 'UMAP')
# Step 3: Partition the cells into supergroups
ee_wild2 <- partitionCells(ee_wild2)
# Step 4: Learn the principal graph
ee_wild2 <- learnGraph(ee_wild2,  RGE_method = 'SimplePPT')
# Step 5: Visualize the trajectory
plot_cell_trajectory(ee_wild2,
                     color_by = "Cluster") +
  scale_color_manual(values = c("red","green","blue"))





















