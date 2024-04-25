library(rtracklayer)
gtf <- as.data.frame(import("~/software/msPIPE/msPIPE/reference/galGal6/galGal6.ncbiRefSeq.gtf"))
head(gtf)
gtf <- gtf[gtf$type == "transcript",]
head(gtf)
gtf$Newend <- ifelse(gtf$strand == "-",gtf$end + 1000,gtf$start)
gtf$Newstart <- ifelse(gtf$strand == "-",gtf$end,gtf$start -1000)
head(gtf)
gtf$name <- paste0("upstream1K:",gtf$gene_id,";",gtf$gene_name)
gtf$point <- "."
write.table(gtf[,c("seqnames","Newstart","Newend","name","point","strand")],row.names = FALSE,col.names = FALSE,quote = FALSE,sep = "\t",file = "R-prepared.promoter.bed")
dim(gtf)

# compare fasta 
library(seqinr)
ucsc <- read.fasta("~/software/msPIPE/msPIPE/reference/galGal6/galGal6.fa")
ensembl <- read.fasta("~/zhangyc/WGBS/reference/fasta/Gallus_gallus.GRCg6a.dna_sm.toplevel.fa")



ucsc_len <- sort(unlist(lapply(ucsc, function(x)length(x))))
ensembl_len <- sort(unlist(lapply(ensembl, function(x)length(x))))
table(ucsc_len == ensembl_len)
ucsc_len[duplicated(ucsc_len)]

chr_mapping <- names(ucsc_len)
names(chr_mapping) <- names(ensembl_len)
# top_100 <- lapply(ensembl, function(x)x[1:1000])
# for (i in names(ucsc)) {
#   print(names(top_100)[which(unlist(lapply(top_100,function(x)all(x==ucsc[[i]][1:1000]))))])
# }
# 
chr_mapping

names(ucsc_len[ucsc_len %in% (ucsc_len[duplicated(ucsc_len)])])
names(ensembl_len[ensembl_len %in% (ensembl_len[duplicated(ensembl_len)])])

table(as.character(ucsc$chrUn_NW_020109980v1) == as.character(ensembl$`AADN05001305.1`))
chr_mapping["AADN05001305.1"] <- "chrUn_NW_020109980v1"
table(as.character(ucsc$chrUn_NW_020109977v1) == as.character(ensembl$`AADN05001173.1`))
chr_mapping["AADN05001173.1"] <- "chrUn_NW_020109977v1"
table(as.character(ucsc$chrUn_NW_020110031v1) == as.character(ensembl$`AADN05001446.1`))
chr_mapping["AADN05001446.1"]
table(as.character(ucsc$chrW_NW_020109797v1_random) == as.character(ensembl$`AADN05001453.1`))
chr_mapping["AADN05001453.1"]

chr_mapping

write.table(data.frame(ensembl = names(chr_mapping),
                       ucsc = unname(chr_mapping)),row.names = FALSE,col.names = TRUE,file = "chrosome-mapping-of-ucsc-and-ensembl.txt")

length(as.character(ucsc$chr1))
length(as.character(ensembl$`1`))
table(as.character(ucsc$chr1) == as.character(ensembl$`1`))
table(as.character(ucsc$chr2) == as.character(ensembl$`2`))
table(as.character(ucsc$chr3) == as.character(ensembl$`3`))
table(as.character(ucsc$chr4) == as.character(ensembl$`4`))
table(as.character(ucsc$chr5) == as.character(ensembl$`5`))
table(as.character(ucsc$chr6) == as.character(ensembl$`6`))
table(as.character(ucsc$chr7) == as.character(ensembl$`7`))
table(as.character(ucsc$chr8) == as.character(ensembl$`8`))
table(as.character(ucsc$chr9) == as.character(ensembl$`9`))
table(as.character(ucsc$chr10) == as.character(ensembl$`10`))
table(as.character(ucsc$chr11) == as.character(ensembl$`11`))
table(as.character(ucsc$chr12) == as.character(ensembl$`12`))
table(as.character(ucsc$chr13) == as.character(ensembl$`13`))
table(as.character(ucsc$chr14) == as.character(ensembl$`14`))
table(as.character(ucsc$chr15) == as.character(ensembl$`15`))
table(as.character(ucsc$chr16) == as.character(ensembl$`16`))
table(as.character(ucsc$chr17) == as.character(ensembl$`17`))
table(as.character(ucsc$chr18) == as.character(ensembl$`18`))
table(as.character(ucsc$chr19) == as.character(ensembl$`19`))
table(as.character(ucsc$chr20) == as.character(ensembl$`20`))
table(as.character(ucsc$chr21) == as.character(ensembl$`21`))
table(as.character(ucsc$chr22) == as.character(ensembl$`22`))
table(as.character(ucsc$chr23) == as.character(ensembl$`23`))
table(as.character(ucsc$chr24) == as.character(ensembl$`24`))
table(as.character(ucsc$chr25) == as.character(ensembl$`25`))
table(as.character(ucsc$chr26) == as.character(ensembl$`26`))
table(as.character(ucsc$chr27) == as.character(ensembl$`27`))
table(as.character(ucsc$chr28) == as.character(ensembl$`28`))
table(as.character(ucsc$chr30) == as.character(ensembl$`30`))


################## for Ensembl gtf file
library(rtracklayer)
gtf <- as.data.frame(import("~/zhangyc/WGBS/reference/Gallus_gallus.GRCg6a.100.gtf"))
head(gtf)
gtf <- gtf[gtf$type == "transcript",]
gtf$seqnames <- as.character(gtf$seqnames)
gtf$seqnames <- chr_mapping[gtf$seqnames]
# gtf <- gtf[!is.na(gtf$gene_name),] # too much NA in gene name
# gtf <- gtf[as.character(gtf$seqnames) %in% levels(gtf$seqnames)[1:35],]
# gtf$seqnames <- paste0("chr",as.character(gtf$seqnames))
# gtf$seqnames[gtf$seqnames == "chrMT"] <- "chrM"
head(gtf)
gtf$Newend <- ifelse(gtf$strand == "-",gtf$end + 1000,gtf$start)
gtf$Newstart <- ifelse(gtf$strand == "-",gtf$end,gtf$start -1000)
head(gtf)
gtf$name <- paste0("upstream1K:",gtf$gene_id,";",gtf$gene_name)
gtf$point <- "."
gtf$Newstart[gtf$Newstart < 0] <- 0 
write.table(gtf[,c("seqnames","Newstart","Newend","name","point","strand")],row.names = FALSE,col.names = FALSE,quote = FALSE,sep = "\t",file = "R-Ensembl-prepared.promoter.bed")
dim(gtf)

dmc <- read.table("./Analysis/DMR/D1.D63/q0.5/reform.DMC_q0.5.bed",stringsAsFactors = FALSE)
head(dmc)
unique(gtf$seqnames)[!unique(gtf$seqnames) %in% unique(dmc$V1)]

# compare DMC gene list
ncbi <- read.table("./Analysis/DMR/D1.D63/q0.5/DMC_genelist.txt",stringsAsFactors = FALSE)
ensembl <- read.table("./Analysis/DMR/D1.D63/q0.5/DMC_genelist_Ensembl_GeneSymbol.txt",stringsAsFactors = FALSE)
table(ncbi$V1 %in% ensembl$V1)
table(ensembl$V1 %in% ncbi$V1 )

