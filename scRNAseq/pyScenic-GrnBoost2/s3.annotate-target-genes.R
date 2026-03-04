library(rtracklayer)
gtf_data <- import("~/reference/cellranger9_reference/customref-gex-mm10-2024-A/genes.gtf", format = "gtf")
gtf_data <- as.data.frame(gtf_data)
head(gtf_data)
table(gtf_data$type)
gtf_data <- gtf_data[gtf_data$type == 'gene',]
head(gtf_data)
gtf_data <- gtf_data[,c('gene_name', 'gene_type', 'gene_id','mgi_id', 'seqnames', 'start', 'end')]
head(gtf_data)

res <- read.csv('./pyscenic-GRNBoost2.csv', stringsAsFactors = FALSE)
res$X <- NULL
head(res)
res <- merge(res, gtf_data, by.x = 'target', by.y = 'gene_name', all.x = TRUE)
res <- res[order(res$importance, decreasing = TRUE),]
head(res)
write.csv(res, file = './pyscenic-GRNBoost2-annotated.csv',row.names = FALSE)
