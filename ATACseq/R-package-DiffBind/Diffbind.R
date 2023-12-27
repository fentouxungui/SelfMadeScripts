# ls *.bam | while read id
# do
# macs2 callpeak \
# -t $id \
# -f BAMPE \
# -g dm \
# -n ${id%.bam} \
# --outdir ../2_macs2_6samples --shift -100 --extsize 200
# done



library(DiffBind, quietly = TRUE)
library(dplyr, quietly = TRUE)
library(parallel, quietly = TRUE)

# https://www.jieandze1314.com/post/cnposts/214/ 
samples <- read.csv("samplesheet.csv", header = TRUE)
samples
sampleDba <- dba(sampleSheet = samples)
sampleCount <- dba.count(sampleDba)

sampleNor <- dba.normalize(sampleCount)
sampleCon <- dba.contrast(sampleNor, reorderMeta = list(Condition="WT"))
sampleCon
sampleDiff <- dba.analyze(sampleCon, bBlacklist = FALSE, bGreylist = FALSE)

dba.show(sampleDiff, bContrasts = TRUE)
diffPeak <- dba.report(sampleDiff,bCalled = TRUE)
diffPeak
table(diffPeak$Fold > 0)

allPeak <- dba.report(sampleDiff,th=1,bCalled = TRUE)
allPeak

# save as bed file
library(rtracklayer)
export.bed(object = allPeak,con = "AllPeaks.bed")
export.bed(object = diffPeak[diffPeak$Fold > 0 ],con = "26-vs-esg-Peaks.Up.bed")
export.bed(object = diffPeak[diffPeak$Fold < 0 ],con = "26-vs-esg-Peaks.Down.bed")



dba.plotVolcano(sampleDiff)


dba.plotPCA(sampleDiff,DBA_CONDITION,label=DBA_CONDITION)
dba.plotVenn(sampleDiff, contrast=1, bDB=TRUE,bGain=TRUE, bLoss=TRUE, bAll=FALSE)

# USE ALL SITES
corvals <- dba.plotHeatmap(sampleDiff)
# RPKM fold (RPKM of the ChIP reads divided by RPKM of the control reads)
dba.plotHeatmap(sampleDiff, score=DBA_SCORE_RPKM_FOLD)

# ONLY USE DE PEAKS
hmap <- colorRampPalette(c("red", "black", "green"))(n = 13)
readscores <- dba.plotHeatmap(sampleDiff, contrast=1, correlations=FALSE,scale="row", colScheme = hmap)

# BOX PLOT
pvals <- dba.plotBox(sampleDiff)

library(ChIPseeker)
# can also inlude a list of grange objects
# peak_list = list(All = allPeak, UpSigRegulate =  diffPeak[diffPeak$Fold > 0], DownSigRegulate =  diffPeak[diffPeak$Fold < 0])


covplot(allPeak,weightCol="Fold")

library(TxDb.Dmelanogaster.UCSC.dm6.ensGene)
txdb <- TxDb.Dmelanogaster.UCSC.dm6.ensGene
promoter <- getPromoters(TxDb=txdb, upstream=3000, downstream=3000)
tagMatrix <- lapply(list(allPeak), getTagMatrix, windows=promoter)
tagHeatmap(tagMatrix, xlim=c(-3000, 3000), color="red")
plotAvgProf(tagMatrix, xlim=c(-3000, 3000),xlab="Genomic Region (5'->3')", ylab = "Read Count Frequency")

allPeak.annotated <- annotatePeak(allPeak, tssRegion=c(-3000, 3000),TxDb=txdb, annoDb="org.Dm.eg.db", verbose=FALSE)
write.table(as.data.frame(allPeak.annotated), file = "./AllPeaks.Annotated.txt",quote=F, row.names=F, sep="\t")


