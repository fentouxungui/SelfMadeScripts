#参考： https://www.jianshu.com/p/7a04647f6e48?utm_campaign=hugo&utm_content=note&utm_medium=seo_notes&utm_source=recommendation

library(GenomicRanges)
# library(TxDb.Ggallus.UCSC.galGal6.refGene)
# download CpG islands
# https://genome.ucsc.edu/cgi-bin/hgTables?hgsid=578954849_wF1QP81SIHdfr8b0kmZUOcsZcHYr&clade=vertebrate&org=Chicken&db=galGal6&hgta_group=regulation&hgta_track=gold&hgta_table=0&hgta_regionType=genome&position=chr4%3A45%2C667%2C317-45%2C670%2C831&hgta_outputType=primaryTable&hgta_outFileName=  
# cgGR <- GRanges() read the cpg islands

vFoler <- list.files("../4_methylation_metrics_MethylDackel",pattern = "bedGraph$",full.names = TRUE)
vFoler
for(i in vFoler){   
  
  tag = gsub("_CpG.bedGraph","",basename(vFoler), fixed = TRUE)
  cat("Processing file", basename(i),"\n")    
  
  Tx = read.table(file = i, row.names =NULL, header =FALSE, sep ="\t", skip =1)    
  chr = as.character(Tx[,1])    
  ST = as.numeric(Tx[,2])+1    
  ED = as.numeric(Tx[,3])    
  methyCount = as.numeric(Tx[,5])    
  readCount = methyCount + as.numeric(Tx[,6])    
  gr = GRanges(seqnames = Rle(chr),  ranges = IRanges(ST,ED),  methyCount = methyCount, readCount = readCount)    
  
  gr = gr[countOverlaps(gr, cgGR, type ="any", ignore.strand =TRUE)==1]    
  x = findOverlaps(gr, cgGR, type ="any", ignore.strand =TRUE)    
  cM = data.frame(mcols(gr))    
  agT = aggregate(cM, by = list(subjectHits(x)),FUN="sum")    
  mgr = cgGR[as.integer(agT[,1])]    
  ols(mgr)= agT[,2:3]    
  mgr$SID = tag    
  grList[[tag]]= mgr
}
CpG_Batch = unlist(GRangesList(grList))
save(CpG_Batch, file = paste0("RData/", args$tag,"_CpG_Batch.RData"))