library(EBSeq)
counts <- read.table("../6_Counts_featureCounts/final_counts.txt",stringsAsFactors = FALSE,row.names = 1, header = TRUE)
head(counts)

GeneMat.norep <- as.matrix(counts[,c(6,7)])
colnames(GeneMat.norep) <- c("control","cko")
Sizes.norep=MedianNorm(GeneMat.norep)
EBOut.norep=EBTest(Data=GeneMat.norep,
                   Conditions=as.factor(rep(c("C1","C2"))),
                   sizeFactors=Sizes.norep, maxround=5)
EBDERes.norep=GetDEResults(EBOut.norep)
GeneFC.norep=PostFC(EBOut.norep)

str(EBOut.norep)
all(names(EBOut.norep$C1Mean[[1]]) == names(GeneFC.norep$PostFC))
all(names(EBOut.norep$C2Mean[[1]]) == names(GeneFC.norep$PostFC))
all(rownames(EBOut.norep$PPMat) == names(GeneFC.norep$PostFC))

str(EBDERes.norep)
all(rownames(EBDERes.norep$PPMat) == rownames(counts))

str(GeneFC.norep)

res <- as.data.frame(EBOut.norep$PPMat)
res$C1Mean <- unname(EBOut.norep$C1Mean[[1]])
res$C2Mean <- unname(EBOut.norep$C2Mean[[1]])
res$PostFC <- GeneFC.norep$PostFC
res$RealFC <- GeneFC.norep$RealFC
res <- cbind(GeneMat.norep[rownames(res),],res)
write.csv(res,file = "EBSeq_without_any_replicates.csv",quote = FALSE)


