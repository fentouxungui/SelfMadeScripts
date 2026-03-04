library(Seurat)
library(qs2)
cds <- qs_read("../../default/featal-SI_default.qs2")
df <- cds@meta.data
df <- df[df$RCTD_2_spot_class != "reject",]
dim(df)
df <- df[df$RCTD_2_first_type %in% c('Enterocytes','EECs', 'Goblet cells', 'Paneth cells','IPC-1', 'IPC-2', 'IPC-3', 'ISCs'),]
write.csv(df, file = "epi-metadata-fro-scanpy.csv")
getwd()
