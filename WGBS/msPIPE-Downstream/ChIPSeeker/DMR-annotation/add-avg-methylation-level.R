dmr <- read.delim("./D7.D63_DMR_0.5.anno.txt",stringsAsFactors = FALSE)
head(dmr)

dmc <- read.delim("../../DMC-annotation/q0.5/D7.D63_DMC_q0.5.anno.txt",stringsAsFactors = FALSE)
head(dmc)

# dmr$avg.g1.methyl <- 0
# dmr$avg.g2.methyl <- 0
# 
# library(dplyr)
# system.time(for (i in 1:1000) {
#   temp <- filter(dmc,seqnames == dmr$seqnames[i], start >= dmr$start[i],end <= dmr$end[i] )
#   dmr[i,"avg.g1.methyl"] <- mean(temp$g1.methyl.per)
#   dmr[i,"avg.g2.methyl"] <- mean(temp$g2.methyl.per)
# })
library(parallel)
library(dplyr)

library(parallel)
system.time(g1 <- unlist(mclapply(1:nrow(dmr), function(x){
  temp <- filter(dmc,seqnames == dmr$seqnames[x], start >= dmr$start[x],end <= dmr$end[x] )
  return(mean(temp$g1.methyl.per))
},mc.cores = 24)))

system.time(g2 <- unlist(mclapply(1:nrow(dmr), function(x){
  temp <- filter(dmc,seqnames == dmr$seqnames[x], start >= dmr$start[x],end <= dmr$end[x] )
  return(mean(temp$g2.methyl.per))
},mc.cores = 24)))

dmr$avg.g1.methyl <- g1
dmr$avg.g2.methyl <- g2
table(sign(dmr$avg.g2.methyl - dmr$avg.g1.methyl),dmr$case.control)
write.table(dmr,file = "./D7.D63_DMR_0.5.anno.withMeanMethylationLevel.txt",quote = FALSE,sep = "\t",row.names = FALSE)
