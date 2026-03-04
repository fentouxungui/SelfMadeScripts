library(Seurat)
library(qs2)
library(ggplot2)

cds <- qs_read("../featal-SI_default.qs2")
table(cds@meta.data$RCTD_2_spot_class)
Idents(cds) <- "RCTD_2_spot_class"
cds <- subset(cds, idents = "reject", invert = TRUE)
options(future.globals.maxSize = 30 * 1024^3)

# day15.5
Idents(cds) <- "group"
cds <- subset(cds, idents = "E12.5_SI")
cds

ImageDimPlot(cds, cols = "polychrome", axes = TRUE)



df <- as.data.frame(cds@images$zoom$centroids@coords)
plot(df$x, df$y)

# subset_region
region_cells <- colnames(cds)[df$y > 8570 & df$y < 9170]
pbmc <- subset(cds, cells = region_cells)




df_sub <- as.data.frame(pbmc@images$zoom$centroids@coords)
meta_data <- pbmc@meta.data
meta_data <- meta_data[order(df_sub$y),]
df_sub <- df_sub[order(df_sub$y),]
# bin y values, for each bin, calculate the mean of x and y
df_sub$group <- floor(c(1:nrow(df_sub))/200)
df_sub_list <- split(df_sub,f = df_sub$group)
df_sub_points <- as.data.frame(Reduce(rbind, lapply(df_sub_list, function(x)apply(x, 2,mean))))
model <- lm(df_sub_points$y ~ df_sub_points$x)
model

plot(df_sub$x, df_sub$y)
points(x = df_sub_points$x, y=df_sub_points$y, col = 'red', pch = 15)
abline(model, col = 'blue', cex = 12)
# attention: not support vertical line, if center points are vertical, the lm not work!

cfs <- coef(model)
# > https://blog.csdn.net/qq_32867925/article/details/113799382
cal_distance <- function(x, y, slope, intercept){
  return((abs(slope * x - y + intercept))/sqrt(slope^2 + (-1)^2))
}
cal_distance(x = df_sub_points$x, y = df_sub_points$y, slope=unname(cfs[2]), intercept=unname(cfs[1]))
df_sub$celltype <- as.character(meta_data$RCTD_2_first_type)

plot_order <- rev(c('Enterocytes','EECs', 'Goblet cells', 'Paneth cells','IPC-1', 'IPC-2', 'IPC-3', 'ISCs',
                    'MFPs', 'Telocytes','Intermediate','Trophocytes','Endothelial cells','Pericytes','Ebf1+Postn+ FLCs', 
                    'Smooth muscle progenitors','Smooth muscle cells-1','Smooth muscle cells-2', 'Smooth muscle cells-3',
                    'Neural cells','Mesothelial cells','Dkk2+ FLCs', 'Frzb+ FLCs', 'Ccl19+ FLCs',  'Macrophages', 'T cells', 
                    'Dendritic cells'))

# get color from heatmap
group.colors <- rev(c("#023ea5", "#7d87b9","#bec1d4","#d6bcc0","#bb7783","#8e063a","#4a70e3","#8594e1","#b5bbe3",
                      "#e6afb9", "#e07b91", "#d33f69", "#11c638", "#8dd593", "#c6dec7", "#ead3c6", "#f0ba8d", "#ef9608",
                      "#0fcfbf", "#9cded6", "#d5eae7", "#f3e1eb", "#f6c4e1", "#f79cd4", "#9f9f9f", "#d8d8d8","#1ce5ff"))

names(group.colors) <- plot_order
df_sub$celltype <- factor(as.character(df_sub$celltype), levels = rev(plot_order))

ggplot(df_sub, aes(x = x, y = y, color = celltype)) +
  geom_point(size = 0.3) +
  geom_point(data = df_sub_points, aes(x = x, y = y), size = 1, alpha = 1, color = 'red') +
  #xlim(14600, 15500) +
  #ylim(3500, 9500) +
  theme_bw() +
  geom_abline(intercept = unname(cfs[1]), slope = unname(cfs[2]), color="blue", size=0.5) +
  # theme(legend.position="none") +
  scale_color_manual(values=group.colors) +
  coord_fixed(ratio = 1, xlim = NULL, ylim = NULL, expand = TRUE, clip = "on") +
  theme(axis.text.x = element_text(angle = -45, hjust=0))
  
# library(ggplot2)
# ImageDimPlot(pbmc, cols = "polychrome", axes = TRUE,group.by = "RCTD_2_first_type") +
#   coord_flip()




df_celltype <- split(df_sub, f = as.character(meta_data$RCTD_2_first_type))
df_celltype_list <- lapply(df_celltype, function(x)data.frame(dis = cal_distance(x=x$x,y=x$y,slope=unname(cfs[2]), intercept=unname(cfs[1])),
                                                              celltype = x$celltype))
df_celltype_res <- Reduce(rbind, df_celltype_list)
# only keep cells with distance < 350
df_celltype_res <- df_celltype_res[df_celltype_res$dis < 350, ]

# change cell type to celltype + cellcounts
table_res <- c(table(df_celltype_res$celltype))
name_mapping <- paste0(names(table_res), " [", unname(table_res), "]")
names(name_mapping) <- names(table_res)
df_celltype_res$celltype <- name_mapping[df_celltype_res$celltype]
df_celltype_res$celltype <- factor(df_celltype_res$celltype, levels = name_mapping[plot_order])
names(group.colors) <- unname(name_mapping[plot_order])

# plot all cell types
ggplot(df_celltype_res, aes(x = celltype, y = dis, fill = celltype)) +
  geom_boxplot(outlier.colour="gray", outlier.size = 0.2) + 
  coord_flip() +
  # geom_dotplot(binaxis='y', stackdir='center', dotsize = 0.01, alpha = 0.01)  + 
  # geom_jitter(shape=16, position=position_jitter(0.2)) +
  ylim(0,160) + 
  # theme(legend.position="none") +
  labs(title = "Distance to center of E12.5") +
  ylab("Distance to Center (μm)") +
  xlab("Cell Type") +
  scale_fill_manual(values=group.colors) +
  scale_x_discrete(drop = FALSE)

# plot only epithelial cells
epi_celltypes <- rev(c('Enterocytes','EECs', 'Goblet cells', 'Paneth cells','IPC-1', 'IPC-2', 'IPC-3', 'ISCs'))
df_celltype_res <- Reduce(rbind, df_celltype_list)
df_celltype_res_sub <- df_celltype_res[df_celltype_res$celltype %in% epi_celltypes, ]
df_celltype_res_sub$celltype <- factor(df_celltype_res_sub$celltype, levels = epi_celltypes)
table_res <- c(table(df_celltype_res_sub$celltype))
name_mapping <- paste0(names(table_res), " [", unname(table_res), "]")
names(name_mapping) <- names(table_res)
df_celltype_res_sub$celltype <- name_mapping[df_celltype_res_sub$celltype]
df_celltype_res_sub$celltype <- factor(df_celltype_res_sub$celltype, levels = name_mapping[epi_celltypes])
group.colors <- rev(c("#1f77b4", "#ff7e0e", "#2ca02c", "#d62727", "#9467bd", "#8c564b", "#e377c3", "#9f9f9f"))
names(group.colors) <- unname(name_mapping)

# plot all cell types
ggplot(df_celltype_res_sub, aes(x = celltype, y = dis, fill = celltype)) +
  geom_boxplot(outlier.colour="gray", outlier.size = 0.2) + 
  coord_flip() +
  # geom_dotplot(binaxis='y', stackdir='center', dotsize = 0.01, alpha = 0.01)  + 
  # geom_jitter(shape=16, position=position_jitter(0.2)) +
  ylim(0,130) + 
  # theme(legend.position="none") +
  labs(title = "Distance to center of E14.5") +
  ylab("Distance to Center (μm)") +
  xlab("Cell Type") +
  scale_fill_manual(values=group.colors) +
  scale_x_discrete(drop = FALSE)







