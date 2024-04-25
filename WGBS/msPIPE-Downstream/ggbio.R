library(ggbio)
library(genomation)
peaks.files <- list.files("./",pattern = "txt$")
i <- peaks.files[1]
i
temp <- read.delim(i,stringsAsFactors = FALSE)
dim(temp)
#temp <- temp[temp$annotation != "Distal Intergenic",]
temp <- temp[!(grepl("^chrUn",temp$seqnames)),]
temp <- temp[!(grepl("random",temp$seqnames)),]
temp$seqnames <- gsub("chr","",temp$seqnames)
temp$start <- temp$start - 1

# temp <- temp[sample(1:nrow(temp),40000),]
temp[temp$q_val < 1e-30,"q_val"] <-  1e-30

library(GenomicRanges)
gr <- GRanges(
  seqnames = factor(temp$seqnames,levels = c(1:28,30:33,"W","Z")),
  ranges = IRanges(temp$start, temp$end, names = temp$name),
  strand = temp$strand,
  score = temp$score,
  qvalue = temp$q_val,
  methdiff = temp$meth_diff)
gr



plotGrandLinear(gr, aes(y = -log10(qvalue)), 
                # color = c("#3399ff", "#707070"),
                alpha = 1)  + 
  theme(axis.text.x = element_text(vjust = 0.5, hjust = 0.5, angle = 90)) + 
  labs(y= "-log10(q value)", x = "Chrosome") +
  scale_y_continuous(breaks = c(0,2,5,10,20,50,100))
  # scale_y_continuous(trans = scales::log10_trans(),breaks = c(0,2,5,10,20,50,100))



# head(temp[order(temp$q_val),],10)











