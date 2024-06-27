promoter.len <- 2000
transcripts <- read.table("./transcripts.gtf",stringsAsFactors = FALSE)
chrom.sizes <- read.table("./chrom.sizes",stringsAsFactors = FALSE)
genes <- read.table("./genes.gtf",stringsAsFactors = FALSE)
#################### filter transcripts based on DEGs
degs <- read.delim("../results_gene_annotated_significant.txt",stringsAsFactors = FALSE)
table(degs$X %in% transcripts$V13)
table(degs$symbol %in% transcripts$V13)
table(degs$ensembl %in% transcripts$V10)

Genes <- read.table("geneInfo.tab",stringsAsFactors = FALSE)
table(degs$symbol %in% Genes$V2)
degs$FBgn <- Genes$V1[match(degs$X, Genes$V2)]
table(degs$FBgn %in% transcripts$V10)
write.csv(degs,file = "for-flybase-update-FBgn-IDs.csv")
updated <- read.delim("FlyBase_Fields_download.txt",stringsAsFactors = FALSE)
table(updated$FBID_KEY != "")
transcripts <- transcripts[transcripts$V10 %in% updated$FBID_KEY,]
genes <- genes[genes$V10 %in% updated$FBID_KEY,]
########### for promoter
# bed start from 0
transcripts$V4 <- transcripts$V4 -1
transcripts$V5 <- transcripts$V5 -1
table(transcripts$V7)
transcripts[!transcripts$V7 %in% c("-","+"),]
transcripts[!transcripts$V7 %in% c("-","+"),"V7"] <- "-"
transcripts[transcripts$V7== "-","V4"] <- transcripts[transcripts$V7== "-","V5"]
transcripts[transcripts$V7== "-","V5"] <- transcripts[transcripts$V7== "-","V5"] + promoter.len
transcripts[transcripts$V7== "+","V5"] <- transcripts[transcripts$V7== "+","V4"]
transcripts[transcripts$V7== "+","V4"] <- transcripts[transcripts$V7== "+","V4"] - promoter.len
transcripts$name <- paste(transcripts$V10, transcripts$V13, transcripts$V16, sep = ";")
transcripts <- transcripts[,c("V1","V4","V5","name")]
transcripts$V4[transcripts$V4 < 0] <- 0
size.mapping <- chrom.sizes$V2
names(size.mapping) <- chrom.sizes$V1
transcripts$V5 <- pmin(transcripts$V5,size.mapping[transcripts$V1]) # parallel minimal value of two or more given vectors
transcripts$V1[transcripts$V1 == "mitochondrion_genome"] <- "M" # 修改染色体命名
transcripts$V1 <- paste0("chr", transcripts$V1)
write.table(transcripts, file = "promoter.2kb.bed",row.names = FALSE,col.names = FALSE,quote = FALSE,sep = "\t")
############## for gene body
genes$name <- paste(genes$V10, genes$V13, sep = ";")
genes <- genes[,c("V1","V4","V5","name")]
write.table(genes, file = "genebody.bed",row.names = FALSE,col.names = FALSE,quote = FALSE,sep = "\t")

