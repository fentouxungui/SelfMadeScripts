library(scran)
require("Matrix")

# import ee and esg data, construct SingleCellExperiment Object
cellbarcodes <- read.table("../ee_filtered_gene_bc_matrices/BDGP6/barcodes.tsv")
genenames <- read.table("../ee_filtered_gene_bc_matrices/BDGP6/genes.tsv")
molecules <- Matrix::readMM("../ee_filtered_gene_bc_matrices/BDGP6/matrix.mtx")

head(cellbarcodes)
head(genenames)
rownames(molecules) <- genenames[,2]
colnames(molecules) <- paste("ee", cellbarcodes[,1], sep="_")
head(molecules)

ee <- SingleCellExperiment(list(counts=molecules))

cellbarcodes <- read.table("../esg_filtered_gene_bc_matrices/BDGP6/barcodes.tsv")
genenames <- read.table("../esg_filtered_gene_bc_matrices/BDGP6/genes.tsv")
molecules <- Matrix::readMM("../esg_filtered_gene_bc_matrices/BDGP6/matrix.mtx")

head(cellbarcodes)
head(genenames)
rownames(molecules) <- genenames[,2]
colnames(molecules) <- paste("esg", cellbarcodes[,1], sep="_")
head(molecules)

esg <- SingleCellExperiment(list(counts=molecules))

# check if there are spike-in transcripts - ERCC
is.spike <- grepl("^ERCC", genenames$V2)
sum(is.spike) #0 No ERCC

library(scater)
library(scran)
# Quality control and normalization by scater and scran
# ee
ee <- calculateQCMetrics(ee, compact=TRUE)
QC <- ee$scater_qc
low.lib <- isOutlier(QC$all$log10_total_counts, type="lower", nmad=3)
low.genes <- isOutlier(QC$all$log10_total_features_by_counts, type="lower", nmad=3)
high.spike <- isOutlier(QC$feature_control_ERCC$pct_counts, type="higher", nmad=3)
data.frame(LowLib=sum(low.lib), LowNgenes=sum(low.genes), 
           HighSpike=sum(high.spike, na.rm=TRUE))
# high.spike is 0, remove from the discard
discard <- low.lib | low.genes
ee <- ee[,!discard]
summary(discard)

set.seed(1000)    
clusters <- quickCluster(ee, method="igraph", min.mean=0.1)
table(clusters)

ee <- computeSumFactors(ee, min.mean=0.1, clusters=clusters)
summary(sizeFactors(ee))

ee <- normalize(ee)

# Identifying highly variable genes
fit <- trendVar(ee, parametric=TRUE,use.spikes = FALSE)
dec <- decomposeVar(ee, fit)

plot(dec$mean, dec$total, xlab="Mean log-expression", 
     ylab="Variance of log-expression", pch=16)
curve(fit$trend(x), col="dodgerblue", add=TRUE)

dec.ee <- dec
dec.ee$Symbol <- rowData(ee)$Symbol
dec.ee <- dec.ee[order(dec.ee$bio, decreasing=TRUE),]
head(dec.ee) # bio factors for choosing the feature
write.csv(dec.ee,file = "ee_bio_factors_decresing.csv")

# esg
esg <- calculateQCMetrics(esg, compact=TRUE)
QC <- esg$scater_qc
low.lib <- isOutlier(QC$all$log10_total_counts, type="lower", nmad=3)
low.genes <- isOutlier(QC$all$log10_total_features_by_counts, type="lower", nmad=3)
high.spike <- isOutlier(QC$feature_control_ERCC$pct_counts, type="higher", nmad=3)
data.frame(LowLib=sum(low.lib), LowNgenes=sum(low.genes), 
           HighSpike=sum(high.spike, na.rm=TRUE))
# high.spike is 0, remove from the discard
discard <- low.lib | low.genes
esg <- esg[,!discard]
summary(discard)

set.seed(1000)    
clusters <- quickCluster(esg, method="igraph", min.mean=0.1)
table(clusters)

esg <- computeSumFactors(esg, min.mean=0.1, clusters=clusters)
summary(sizeFactors(esg))

esg <- normalize(esg)

# Identifying highly variable genes
fit <- trendVar(esg, parametric=TRUE,use.spikes = FALSE)
dec <- decomposeVar(esg, fit)

plot(dec$mean, dec$total, xlab="Mean log-expression", 
     ylab="Variance of log-expression", pch=16)
curve(fit$trend(x), col="dodgerblue", add=TRUE)

dec.esg <- dec
dec.esg$Symbol <- rowData(esg)$Symbol
dec.esg <- dec.esg[order(dec.esg$bio, decreasing=TRUE),]
head(dec.esg) # bio factors for choosing the feature
write.csv(dec.esg,file = "esg_bio_factors_decresing.csv")

#  it is often possible to obtain good results when applying MNN correction 
# to batches of data generated with different technologies.

# Feature selection across batches
universe <- intersect(rownames(dec.ee), rownames(dec.esg))
mean.bio <- (dec.ee[universe,"bio"] + dec.esg[universe,"bio"])/2
chosen <- universe[mean.bio > 0]  # should be >0,but the gene numbers is too short, only 4776
length(chosen)

rescaled <- multiBatchNorm(ee[universe,], esg[universe,])
ee <- rescaled[[1]]
esg <- rescaled[[2]]


# Performing MNN-based correction
set.seed(100) 
original <- list(
  ee=logcounts(ee)[chosen,],
  esg=logcounts(esg)[chosen,]
)
# Slightly convoluted call to avoid re-writing code later.
# Equivalent to fastMNN(GSE81076, GSE85241, k=20, d=50, approximate=TRUE)
mnn.out <- do.call(fastMNN, c(original, list(k=30,d=50, approximate=TRUE,auto.order=TRUE)))
dim(mnn.out$corrected) # Row:cells  Col:top 50 Principle compoments
mnn.out$batch
mnn.out$pairs # the number of MNN pairs

# Examining the effect of correction
omat <- do.call(cbind, original)
sce <- SingleCellExperiment(list(logcounts=omat))
reducedDim(sce, "MNN") <- mnn.out$corrected
sce$Batch <- as.character(mnn.out$batch)
sce

set.seed(100)
# Using irlba to set up the t-SNE, for speed.
osce <- runPCA(sce, ntop=Inf, method="irlba")
osce <- runTSNE(osce, use_dimred="PCA")
ot <- plotTSNE(osce, colour_by="Batch") + ggtitle("Original")

set.seed(100)
csce <- runTSNE(sce, use_dimred="MNN")
ct <- plotTSNE(csce, colour_by="Batch") + ggtitle("Corrected")

multiplot(ot, ct, cols=2)

ct.gcg <- plotTSNE(csce, colour_by="Tk") + 
  ggtitle("Tk")
ct.ins <- plotTSNE(csce, colour_by="AstC") + 
  ggtitle("AstC")
ct.sst <- plotTSNE(csce, colour_by="NPF") + 
  ggtitle("NPF")
ct.ppy <- plotTSNE(csce, colour_by="CG14989") + 
  ggtitle("CG14989")
multiplot(ct.gcg, ct.ins, ct.sst, ct.ppy, cols=2)


## Using the corrected values in downstream analyses
snn.gr <- buildSNNGraph(sce, use.dimred="MNN")
clusters <- igraph::cluster_walktrap(snn.gr)
table(clusters$membership, sce$Batch)

csce$Cluster <- factor(clusters$membership)
plotTSNE(csce, colour_by="Cluster")


m.out <- findMarkers(sce, clusters$membership, block=sce$Batch,
                     direction="up")        
demo <- m.out[["4"]] # looking at cluster 4 (probably alpha cells).
demo <- demo[demo$Top <= 5,]
head(demo)

# Setting up the design matrix (we remove intercept for full rank
# in the final design matrix with the cluster-specific terms).
design <- model.matrix(~sce$Batch)
design <- design[,-1,drop=FALSE]

m.alt <- findMarkers(sce, clusters$membership, design=design,
                     direction="up")
demo <- m.alt[["4"]]
demo <- demo[demo$Top <= 5,]
head(demo)
















