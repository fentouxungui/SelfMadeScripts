m1 <- read.delim("./match-AGGGCGG.txt",stringsAsFactors = FALSE)
m2 <- read.delim("./match-AGGGTGG.txt",stringsAsFactors = FALSE)
head(m1)
m1$seqID <- gsub("\\(.*","",m1$seqID)
m2$seqID <- gsub("\\(.*","",m2$seqID)
head(m2)

meta <- read.table("./genes.txt",stringsAsFactors = FALSE)
head(meta)
gene.meta <- meta[,c(2,5,11)]
head(gene.meta)
colnames(gene.meta) <- c("FBgn","Symbol","Type")

temp <- as.data.frame(as.matrix(table(m1$seqID)))
mapping <- temp$V1
names(mapping) <- rownames(temp)

gene.meta$AGGGCGG <- mapping[gene.meta$FBgn]

temp <- as.data.frame(as.matrix(table(m2$seqID)))
mapping <- temp$V1
names(mapping) <- rownames(temp)

gene.meta$AGGGTGG <- mapping[gene.meta$FBgn]
head(gene.meta)
gene.meta[is.na(gene.meta)] <- 0
sum(gene.meta$AGGGCGG)
sum(gene.meta$AGGGTGG)

counts.df <- gene.meta
head(counts.df)
degs.df <- read.csv("./Supply-1.csv",stringsAsFactors = FALSE)
head(degs.df)
degs.df <- degs.df[apply(degs.df[,3:5], 1, mean) < apply(degs.df[,9:11], 1, mean),]
dim(degs.df)
table(degs.df$FBgn.latest %in% counts.df$FBgn)
table(degs.df$gene.symbol %in% counts.df$Symbol)

counts.df$motif <- apply(counts.df[,4:5],1,sum) > 0

degs.df <- degs.df[degs.df$FBgn.latest %in% counts.df$FBgn,]
dim(degs.df)
table(counts.df$motif)
counts.df.sub <- counts.df[counts.df$FBgn %in% degs.df$FBgn.latest,]
table(counts.df.sub$motif)

fisher.test(x = as.vector(table(counts.df$motif)), y =as.vector(table(counts.df.sub$motif)))

###############

m1 <- read.delim("./match-AGGACAA.txt",stringsAsFactors = FALSE)
m2 <- read.delim("./match-AGGATAA.txt",stringsAsFactors = FALSE)
head(m1)
m1$seqID <- gsub("\\(.*","",m1$seqID)
m2$seqID <- gsub("\\(.*","",m2$seqID)
head(m2)

meta <- read.table("./genes.txt",stringsAsFactors = FALSE)
head(meta)
gene.meta <- meta[,c(2,5,11)]
head(gene.meta)
colnames(gene.meta) <- c("FBgn","Symbol","Type")

temp <- as.data.frame(as.matrix(table(m1$seqID)))
mapping <- temp$V1
names(mapping) <- rownames(temp)

gene.meta$AGGACAA <- mapping[gene.meta$FBgn]

temp <- as.data.frame(as.matrix(table(m2$seqID)))
mapping <- temp$V1
names(mapping) <- rownames(temp)

gene.meta$AGGATAA <- mapping[gene.meta$FBgn]
head(gene.meta)
gene.meta[is.na(gene.meta)] <- 0
sum(gene.meta$AGGACAA)
sum(gene.meta$AGGATAA)

counts.df <- gene.meta
head(counts.df)
degs.df <- read.csv("./Supply-1.csv",stringsAsFactors = FALSE)
head(degs.df)
degs.df <- degs.df[apply(degs.df[,3:5], 1, mean) < apply(degs.df[,9:11], 1, mean),]

table(degs.df$FBgn.latest %in% counts.df$FBgn)
table(degs.df$gene.symbol %in% counts.df$Symbol)

counts.df$motif <- apply(counts.df[,4:5],1,sum) > 0

degs.df <- degs.df[degs.df$FBgn.latest %in% counts.df$FBgn,]
dim(degs.df)
table(counts.df$motif)
counts.df.sub <- counts.df[counts.df$FBgn %in% degs.df$FBgn.latest,]
table(counts.df.sub$motif)

fisher.test(x = as.vector(table(counts.df$motif)), y =as.vector(table(counts.df.sub$motif)))

###############################

###############

m1 <- read.delim("./match-AAGCT.txt",stringsAsFactors = FALSE)
m1$seqID <- gsub("\\(.*","",m1$seqID)
head(m1)

meta <- read.table("./genes.txt",stringsAsFactors = FALSE)
head(meta)
gene.meta <- meta[,c(2,5,11)]
head(gene.meta)
colnames(gene.meta) <- c("FBgn","Symbol","Type")

temp <- as.data.frame(as.matrix(table(m1$seqID)))
mapping <- temp$V1
names(mapping) <- rownames(temp)

gene.meta$AAGCT<- mapping[gene.meta$FBgn]
head(gene.meta)
gene.meta[is.na(gene.meta)] <- 0
sum(gene.meta$AAGCT)
counts.df <- gene.meta
head(counts.df)
degs.df <- read.csv("./Supply-1.csv",stringsAsFactors = FALSE)
head(degs.df)
degs.df <- degs.df[apply(degs.df[,3:5], 1, mean) < apply(degs.df[,9:11], 1, mean),]

table(degs.df$FBgn.latest %in% counts.df$FBgn)
table(degs.df$gene.symbol %in% counts.df$Symbol)

counts.df$motif <- counts.df[,4] > 0

degs.df <- degs.df[degs.df$FBgn.latest %in% counts.df$FBgn,]
dim(degs.df)
table(counts.df$motif)
counts.df.sub <- counts.df[counts.df$FBgn %in% degs.df$FBgn.latest,]
table(counts.df.sub$motif)

fisher.test(x = as.vector(table(counts.df$motif)), y =as.vector(table(counts.df.sub$motif)))


