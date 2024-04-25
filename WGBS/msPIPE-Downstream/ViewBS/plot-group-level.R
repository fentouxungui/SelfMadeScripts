meth <- "./MergedSingle_MethOverRegion_CG.txt"
tab <- read.table(meth,head=T,stringsAsFactors = FALSE)
tab <- tab[order(tab[,3]),]

head(tab)
tab$group <- factor(gsub("-.*$","",tab$sample_name),levels = c("D1","D7","D63"))
library(dplyr)

## optional
tab <- group_by(tab, group, bin_num) %>%
  summarise(sample_name = group[1],
            Methylation_level = mean(Methylation_level),
            bin_num = bin_num[1])




library(ggplot2)
p=ggplot(tab, aes(x=bin_num, y=Methylation_level, group=sample_name, col=sample_name)) +
  geom_line() + xlab("Gene")
p
min <- min(tab$bin_num)  ## by default the lowest value is -19
max <- abs(max(tab$bin_num))
#p = p + scale_x_continuous(breaks=c(min(tab$Pos), 0.5,max(tab$Pos)+min(tab$Pos), max(tab$Pos)), labels=c(abs(min(tab$Pos))/10, "TSS", "TTS", abs(min(tab$Pos))/10))

#p = p + scale_x_continuous(breaks=c(min/2, (max + min)/2, max + min/2), labels=c("Upstream", xlab, "Downstream"))
#p = p + scale_x_continuous(breaks=c(min/2, (max + min)/2, max + min/2), labels=c("Upstream", xlab, "Downstream"))
#p = p + theme(axis.text.x = element_blank())

# flank = paste( -(min -1)/adjustXaxis, "kb", sep = " ")
flank = "2kb"

p = p + scale_x_continuous(breaks=c(min, max), labels=c(flank, flank));

p <- p + theme(legend.title=element_blank()) ## no legend title
#p = p + scale_fill_continuous(guide = guide_legend(title = legend_title)) # title text

## 1 means the first bin in the gene body, max + min + 1 means the last bean in the gene body.

p = p + geom_vline(xintercept = c(1, max + min - 1), linetype = "dashed")
p = p + expand_limits(y=0)
p = p + ylab("Methylation level")
p
p + ylim(c(0.4,0.7))
p + facet_wrap(~group) + ylim(c(0.4,0.7))

