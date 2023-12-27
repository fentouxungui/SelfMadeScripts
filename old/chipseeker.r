#!/PATH/TO/bin/Rscript
# Fentouxungui
# NIBS
# 2019-11-12
# Using Chipseeker to analysis MACS2 outs and prepare data for homer analysis

# usage: ./Chipseeker.r -p xxx.narrowPeak -s xxx.summits.bed -n ProjectName
#        ./Chipseeker.r --peak xxx.narrowPeak --summit xxx.summits.bed --name ProjectName


suppressPackageStartupMessages(require(optparse)) # don't say "Loading required package: optparse"

option_list = list(
  make_option(c("-p", "--peak"), action="store", default=NA, type='character',
              help="Please specify the narrowpeak file from MACS2 output!"),
  make_option(c("-s", "--summit"), action="store", default=NA, type='character',
              help="Please specify the summit file from MACS2 output!"),
  make_option(c("-n", "--projectname"), action="store", default="ChipSeeker",type='character',
              help="Please specify the project name! [default is ChipSeeker]")
)
opt = parse_args(OptionParser(option_list=option_list))

if (!is.na(opt$peak) & !is.na(opt$summit)) {
  library(ChIPseeker)
  library(TxDb.Dmelanogaster.UCSC.dm6.ensGene)
  library(org.Dm.eg.db)
  library(clusterProfiler)
  
  peak.file <- opt$peak
  summit.file <- opt$summit
  if(is.na(opt$projectname)){project.name <- "ChipSeeker"}
  else {project.name <- opt$projectname}
  
  txdb <- TxDb.Dmelanogaster.UCSC.dm6.ensGene 
  nanog <- readPeakFile(peak.file,as = "GRanges",header = FALSE)
  
  file.dir <- getwd()
  new.dir <- paste("ChipSeeker-",project.name,sep = "")
  dir.create(new.dir)
  setwd(paste(file.dir,"/",new.dir,sep = ""))
  
  peaks <- list(Nanog=nanog)
  # promotor区间范围可以自己设定
  promoter <- getPromoters(TxDb=txdb, upstream=3000, downstream=3000)
  tagMatrixList <- lapply(peaks, getTagMatrix, windows=promoter)
  #annotatePeak传入annoDb参数,可进行基因ID转换（Entrez，ENSEMBL，SYMBOL，GENENAME）
  peakAnnoList <- lapply(peaks, annotatePeak, TxDb=txdb,tssRegion=c(-3000, 3000),
                         verbose=FALSE,addFlankGeneInfo=TRUE, flankDistance=5000,annoDb="org.Dm.eg.db")
  
  pdf(file = paste(project.name,"FeatureDistribution.pdf",sep = "-"),width = 9,height = 3)
  print(plotAnnoBar(peakAnnoList))
  dev.off()
  
  pdf(file = paste(project.name,"DistributionRelativeToTSS.pdf",sep = "-"),width = 9,height = 3)
  print(plotDistToTSS(peakAnnoList, title="Distribution of transcription factor-binding loci \n relative to TSS"))
  dev.off()
  
  
  # Create a list with genes from each sample
  gene = lapply(peakAnnoList, function(i) as.data.frame(i)$geneId)
  gene <- unlist(gene)
  ENTREZID <- bitr(gene, fromType = "FLYBASE",
                   toType = c("ENTREZID"),
                   OrgDb = org.Dm.eg.db)
  # Run GO enrichment analysis 
  ego <- enrichGO(gene = ENTREZID$ENTREZID, 
                  #keytype = "ENTREZID", 
                  OrgDb = org.Dm.eg.db, 
                  ont = "BP", 
                  pAdjustMethod = "BH", 
                  qvalueCutoff = 0.05, 
                  readable = TRUE)
  # Dotplot visualization
  pdf(file = paste(project.name,"GO-BP.pdf",sep = "-"),width = 16,height = 12)
  print(dotplot(ego, showCategory=30))
  dev.off()
  
  
  ENTREZID <- bitr(gene, fromType = "FLYBASE",
                   toType = c("FLYBASECG"),
                   OrgDb = org.Dm.eg.db)
  
  ENTREZID$FLYBASECG <- paste0("Dmel_",ENTREZID$FLYBASECG,sep="")
  kk <- enrichKEGG(gene = ENTREZID$FLYBASECG,
                   organism = 'dme',
                   #keyType = 'FLYBASE',
                   pvalueCutoff = 0.05)
  
  pdf(file = paste(project.name,"KEGG.pdf",sep = "-"),width = 12,height = 10)
  print(dotplot(kk, showCategory=15))
  dev.off()
  
  # Output peakAnnolist file
  save(peakAnnoList,file=paste(project.name,"peakAnnolist.rda",sep = "-"))
  results <- as.data.frame(peakAnnoList$Nanog)
  colnames(results)[6:12] <- c("peak name","int(-10*log10qvalue)","strand","fold change","peaks -log10pvalue","peaks -log10qvalue","relative-summit-position-to-peak-start")
  write.csv(results,file=paste(project.name,"Peaks.Annotated.csv",sep = "-"),row.names = FALSE)
  # Output results from GO analysis to a table
  cluster_summary <- data.frame(ego)
  write.csv(cluster_summary, paste(project.name,"clusterProfiler.GO.csv",sep = "-"))
  # Output results from KEGG analysis to a table
  kegg_summary <- data.frame(kk)
  write.csv(kegg_summary, paste(project.name,"KEGG.csv",sep = "-"))
  
  
  
  ####### export for homer:  findMotifsGenome.pl
  ## 1. exact size
  library(dplyr)
  peaks <- read.delim(file=paste("../",peak.file,sep = ""),header = FALSE)
  colnames(peaks)[c(1:4,6)] <- c("chromosome","starting position","ending position","Peak ID","Strand")
  peaks <- select(peaks,"Peak ID","chromosome","starting position","ending position","Strand")
  peaks$Strand <- rep("+",length(peaks$Strand))
  write.table(peaks,file = paste(project.name,"tmp.file.peaks_for_homer.txt",sep = "-"),quote=FALSE,sep="\t",row.names=FALSE)
  
  ## 2. summits
  peaks <- read.delim(file=paste("../",summit.file,sep = ""),header = FALSE)
  head(peaks)
  colnames(peaks) <- c("chromosome","starting position","ending position","Peak ID","Strand")
  peaks <- select(peaks,"Peak ID","chromosome","starting position","ending position","Strand")
  peaks$Strand <- rep("+",length(peaks$Strand))
  write.table(peaks,file = paste(project.name,"tmp.file_summits_for_homer.txt",sep = "-"),quote=FALSE,sep="\t",row.names=FALSE)
} else {
  cat("you didn't specify both variables p (narrowpeak file) and s (summits file)\n!", file=stderr()) # print error messages to stderr
}
