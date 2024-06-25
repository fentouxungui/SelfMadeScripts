# 火山图

```
library(MAGeCKFlute)
library(ggplot2)

# gdata = ReadRRA("./Enriched-vs-CTR.gene_summary.txt")
source("ReadRRA.R") 
gdata = ReadRRA_2("./Enriched-vs-CTR-Rep2.gene_summary.txt") # 使用pvalue， not fdr,
# RearRRA中，神奇的地方是，neglfc和poslfc仅使用绝对值大的！

VolcanoView(gdata, x = "Score", y = "Pvalue", Label = "id", top = 20, x_cut = 2, 
            max.overlaps = 1e6) +
  scale_x_continuous(breaks=seq(-8, 8, 2)) + 
  ylab("-log10(Pvalue)") + xlab("Score") + geom_vline(xintercept = 0, linetype="dotted", 
                                                      color = "black")
```

# ![image-20240507154938148](D:\BaiduNetdiskWorkspace\Learning-bioinformatics-生物信息学笔记\NGS数据分析\Figures\image-20240507154938148.png)

```
# ReadRRA.R
ReadRRA_2 <- function (gene_summary, score = c("lfc", "rra")[1]) 
{
  if (is.null(dim(gene_summary))) {
    gene_summary = read.table(file = gene_summary, sep = "\t", 
                              header = TRUE, quote = "", comment.char = "", check.names = FALSE, 
                              stringsAsFactors = FALSE)
  }
  if (all(c("id", "Score", "FDR") %in% colnames(gene_summary))) {
    dd = as.data.frame(gene_summary[, c("id", "Score", "FDR")], 
                       stringsAsFactors = FALSE)
    dd$id = as.character(dd$id)
    return(dd)
  }
  # gene_summary = gene_summary[, c(1, 3, 9, 8, 14, 5, 11)]
  gene_summary = gene_summary[, c(1, 3, 9, 8, 14, 4, 10)]
  colnames(gene_summary) = c("id", "negscore", "poscore", 
                             "neglfc", "poslfc", "negfdr", "posfdr")
  dd = gene_summary
  if ("lfc" %in% tolower(score)) {
    dd$LFC = dd$poslfc
    dd$FDR = dd$posfdr
    dd$LFC[abs(dd$neglfc) > dd$poslfc] = dd$neglfc[abs(dd$neglfc) > 
                                                     dd$poslfc]
    dd$FDR[abs(dd$neglfc) > dd$poslfc] = dd$negfdr[abs(dd$neglfc) > 
                                                     dd$poslfc]
    dd = dd[, c("id", "LFC", "FDR")]
  }  else if ("rra" %in% tolower(score)) {
    idx_neg = dd$negscore < dd$poscore
    dd$LFC = apply(-log10(dd[, 2:3]), 1, max)
    dd$LFC[idx_neg] = -dd$LFC[idx_neg]
    dd$FDR = dd$posfdr
    dd$FDR[idx_neg] = dd$negfdr[idx_neg]
    dd = dd[, c("id", "LFC", "FDR")]
  }
  colnames(dd) = c("id", "Score", "Pvalue")
  dd$id = as.character(dd$id)
  return(dd)
}

```

