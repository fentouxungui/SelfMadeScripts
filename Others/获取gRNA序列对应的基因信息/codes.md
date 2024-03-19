# 注释基因信息

>  Crispr 双gRNA基因激活的果蝇库

## 1. 比对到基因body上（不靠谱）

后来发现，激活位点是位于promoter区域的！

```R
# 数据清洗
df <- read.csv("./fly-stock.csv",stringsAsFactors = FALSE)
head(df)
anyDuplicated(df$CV.number)
df$primer_sequence_forward
df$Sequence_forward <- gsub("TATATAGGAAAGATATCCGGGTGAACTTCG","", df$primer_sequence_forward)
df$Sequence_forward <- gsub("GTTTTAGAGCTAGAAATAGCAAG", "", df$Sequence_forward)
df$Sequence_forward <- unname(sapply(df$Sequence_forward, function(x){
  if (nchar(x) > 20) {
    gsub("^ttcg", "", x)
  }else{
    x
  }
}))
df$Sequence_forward <- unname(sapply(df$Sequence_forward, function(x){
  if (nchar(x) > 20) {
    gsub("^GTCG", "", x)
  }else{
    x
  }
}))
table(nchar(df$Sequence_forward))

df$Sequence_reverse <- gsub("ATTTTAACTTGCTATTTCTAGCTCTAAAAC","", df$primer_sequence_reverse)
df$Sequence_reverse <- gsub("CGACGTTAAATTGAAAATAGGTC", "", df$Sequence_reverse)
df$Sequence_reverse <- unname(sapply(df$Sequence_reverse, function(x){
  if (nchar(x) > 20) {
    gsub("^aaac", "", x)
  }else{
    x
  }
}))
df$Sequence_reverse <- unname(sapply(df$Sequence_reverse, function(x){
  if (nchar(x) > 20) {
    gsub("^AAAC", "", x)
  }else{
    x
  }
}))
table(nchar(df$Sequence_reverse))
write.csv(df,file = "fly-stock-clean-data.csv",row.names = FALSE)

file.create("input.forward.fa")
for (i in 1:nrow(df)) {
  cat("> ", df$CV.number[i], "\n", file = "input.forward.fa",append = TRUE)
  cat(df$Sequence_forward[i],"\n", file = "input.forward.fa",append = TRUE)
}

file.create("input.reverse.fa")
for (i in 1:nrow(df)) {
  cat("> ", df$CV.number[i], "\n", file = "input.reverse.fa",append = TRUE)
  cat(df$Sequence_reverse[i],"\n", file = "input.reverse.fa",append = TRUE)
}
```

```{shell}
cd /home/xilab/reference/Genome/Drosophlia/Flybase/r6.51
# 仅仅保留type为基因的条目
cat dmel-all-r6.51.gtf |awk '$3 == "gene"' > genes.gtf
```

```{r}
# 得到gene bed文件，问题是：不是从0开始的
df <- read.table("./genes.gtf")
df <- df[,c("V1","V4","V5","V10")]
write.table(df,file = "genes.clean.df",row.names = FALSE,col.names = FALSE,quote = FALSE,sep = "\t")
```

```{shell}
# 获取基因的序列
bedtools getfasta -fi dmel-all-chromosome-r6.51.fasta -bed genes.clean.gtf -name > genes.fa
# 获取基因序列文件中的基因长度
seqkit fx2tab --length genes.fa --name > genes.length.txt
```

```{r}
# 依据长度，来检查得到的基因序列是不是正确
l.1 <- read.table("./genes.clean.df",stringsAsFactors = FALSE)
l.2 <- read.table("./genes.length.txt",stringsAsFactors = FALSE)
l.1$length <- abs(l.1$V2 - l.1$V3)
head(l.1)
head(l.2)
```

```{shell}
# 构建基因序列的blast数据库，比对gRNA序列，从而获取对应的基因
makeblastdb -in genes.fa -dbtype nucl -out db/genedb
blastn -query input.forward.fa -db db/genedb -outfmt 6 -out result.forward.txt -max_target_seqs 1 -word_size 11 -dust no
blastn -query input.reverse.fa -db db/genedb -outfmt 6 -out result.reverse.txt -max_target_seqs 1 -word_size 11 -dust no
```

```R
# 将比对结果整合到结果中
res.f <- read.table("./result.forward.txt",stringsAsFactors = FALSE)
res.r <- read.table("./result.reverse.txt",stringsAsFactors = FALSE)

colNames <- c("query_id",	"refer_id",	"identity",	"alignment_length",	"mismatches",	"gap_openings",	
              "q.start",	"q.end",	"s.start",	"s.end",	"e-value",	"bit_score")
colnames(res.f) <- colNames
colnames(res.r) <- colNames

res.f$refer_id <- gsub("::.*","",res.f$refer_id)
res.r$refer_id <- gsub("::.*","",res.r$refer_id)

res.f.final <- res.f[res.f$q.start == 1 & res.f$identity == 100 & res.f$mismatches == 0 & res.f$`e-value` < 0.1,]
res.r.final <- res.r[res.r$q.start == 1 & res.r$identity == 100 & res.r$mismatches == 0 & res.r$`e-value` < 0.1,]

table(duplicated(res.f.final$refer_id))
table(duplicated(res.f.final$query_id))

table(duplicated(res.r.final$refer_id))
table(duplicated(res.r.final$query_id))

res.f.final[res.f.final$query_id %in% res.f.final$query_id[duplicated(res.f.final$query_id)],]
res.r.final[res.r.final$query_id %in% res.r.final$query_id[duplicated(res.r.final$query_id)],]

res.f.final <- res.f.final[!duplicated(res.f.final$query_id),]
res.r.final <- res.r.final[!duplicated(res.r.final$query_id),]

df <- read.csv("./fly-stock-clean-data.csv",stringsAsFactors = FALSE)
mapping.f <- res.f.final$refer_id
names(mapping.f) <- res.f.final$query_id

mapping.r <- res.r.final$refer_id
names(mapping.r) <- res.r.final$query_id

df$forward_mapped <- mapping.f[df$CV.number]
df$reverse_mapped <- mapping.r[df$CV.number]
write.csv(df, file = "fly-stock-annotated.csv")
```

## 2. 比对到基因组上（只能注释得到最近的基因）

```shell
# 去除空格
# cat input.forward.fa | sed 's/ //g' > input.forward2.fa
# cat input.reverse.fa | sed 's/ //g' > input.reverse2.fa
# 获取gRNA序列在基因组上的位置
seqkit locate -f input.forward.fa dmel-all-chromosome-r6.51.fasta > output.forward.location.txt
seqkit locate -f input.reverse.fa dmel-all-chromosome-r6.51.fasta > ouput.reverse.location.txt
```

```R
# 合并双向gRNA的基因组比对结果

forward <- read.delim("./output.forward.location.txt",stringsAsFactors = FALSE)
reverse <- read.delim("./output.reverse.location.txt",stringsAsFactors = FALSE)
forward$matchBy <- "forward"
reverse$matchBy <- "reverse"
table(duplicated(forward$patternName))

# res: a list,用于存储每对gRNA的比对结果
df <- read.csv("../fly-stock-clean-data.csv",stringsAsFactors = FALSE)
res <- list()
for (i in df$CV.number) {
  forward.sub <- forward[forward$patternName == i,]
  reverse.sub <- reverse[reverse$patternName == i,]
  res[[i]] <- rbind(forward.sub, reverse.sub)
}
table(unlist(lapply(res, function(x){nrow(x)})))

# 后面是对于每一对结果进行处理，理论上gRNA的比对位置不能离得太远！
# 第一步，得到的行数为0，或2，对于非0，2行的，进行人工矫正
res[which(!unlist(lapply(res, function(x){nrow(x)})) %in% c(0,2))]
## 手动去除一些错误的比对位置
res$CV221 <- res$CV221[-2, ]
res$CV236 <- res$CV236[-2, ]
res$CV347 <- res$CV347[-1, ]
res$CV753 <- res$CV753[c(1,4), ]
res$CV844 <- res$CV844[-1, ]
res$CV884 <- res$CV884[c(2,4), ]
res$CV1298 <- res$CV1298[c(1,3), ]
## 输出有多种比对可能的条目
unsolved <- c("CV1001","CV1015","CV1028","CV1051","CV1068","CV1197")
library(dplyr)
for (i in unsolved) {
  write.csv(arrange(res[[i]], seqID,start), file = paste0("unsolved-",i,".csv"))
}
res[unsolved] <- NULL # 去除这几个

# 第二步，行数为2的，检查是否forward和reverse各为1
res.matched <- res[unlist(lapply(res, function(x){nrow(x)})) == 2]
res.matched[!unlist(lapply(res.matched, function(x)all(c("forward", "reverse") %in% x$matchBy)))]
## 检查染色体是不是一样的
res.matched[!unlist(lapply(res.matched, function(x){x$seqID[1] == x$seqID[2]}))]

# 没有问题
# 第三步，检查长度
hist(unlist(lapply(res.matched, function(x){max(c(x$start,x$end)) - min(x$start,x$end)})))

# 把res.matched数据添加到df中
res.matched.arranged <- lapply(res.matched, function(x){
  data.frame(cvID = x$patternName[1],
             seqID = x$seqID[1],
             forward.pattern = x$pattern[x$matchBy == "forward"],
             forward.strand = x$strand[x$matchBy == "forward"],
             forward.start = x$start[x$matchBy == "forward"],
             forward.end = x$end[x$matchBy == "forward"],
             reverse.pattern = x$pattern[x$matchBy == "reverse"],
             reverse.strand = x$strand[x$matchBy == "reverse"],
             reverse.start = x$start[x$matchBy == "reverse"],
             reverse.end = x$end[x$matchBy == "reverse"],
             range.min = min(c(x$start, x$end)),
             range.max = max(c(x$start, x$end)),
             length = max(c(x$start, x$end)) - min(c(x$start, x$end)))
})
results <- Reduce(rbind, res.matched.arranged)
write.csv(results,file = "matched.results.csv")

df.location <- merge(df, results, by.x = "CV.number", by.y = "cvID",all.x = TRUE)
df.location <- df.location[match(df$CV.number, df.location$CV.number),]
write.csv(df.location,file = "../fly-stock-clean-data-with-location.csv")

# 输出bed文件，用于可视化gRNA对的位置和注释
temp <- df.location[!is.na(df.location$length),c("seqID","range.min","range.max","CV.number")]
temp$anno <- "fly-stock-cripsr"
temp$Strand <- "+"
temp$seqID <- paste0("chr",temp$seqID)
write.table(temp,file = "crispr.bed",row.names = FALSE,col.names = FALSE,quote = FALSE,sep = "\t")
```

```shell
# 使用homer注释bed
annotatePeaks.pl crispr.bed dm6 > crispr.annotated.txt
```

```R
# 或者使用Chipseeker注释
library(ChIPseeker)
library(TxDb.Dmelanogaster.UCSC.dm6.ensGene)
txdb <- TxDb.Dmelanogaster.UCSC.dm6.ensGene
library(clusterProfiler)
library(org.Dm.eg.db)
peak <- readPeakFile("../crispr.bed")
peakAnno <- annotatePeak(peak, tssRegion=c(-3000, 3000),level = "gene",
                         TxDb=txdb, annoDb="org.Dm.eg.db")
as.data.frame(peakAnno)
```

```R
# 后续使用HOMER注释的结果进行分析
anno <- read.delim("./crispr.annotated.txt",stringsAsFactors = FALSE)
anno <- anno[,c(1,8,10,11,12,13,14,16,18,19)]
colnames(anno)[1] <- "cvID"

df <- read.csv("../fly-stock-clean-data-with-location.csv",stringsAsFactors = FALSE)
df$X <- NULL
res <- merge(df, anno, by.x = "CV.number", by.y = "cvID", all.x = TRUE)
res <- res[order(as.numeric(gsub("CV","",res$CV.number))),]
write.csv(res,file = "temp.csv")

# homer得到的ncbi ref ID, nm_xxxx类型的，上传到flybase基因进行ID的转换
ref <- as.vector(na.omit(unique(res$Nearest.Refseq)))
ref <- ref[ref != ""]
write.csv(data.frame(ref = ref),file = "Nearest-Ref-for-flybase-batch-download.csv")
# 发现尽管nearest promoter列与nearest refseq列不一样的，但是最终对应的flybase fbgn id是全部一样的！
ref <- as.vector(na.omit(unique(res$Nearest.PromoterID)))
ref <- ref[ref != ""]
write.csv(data.frame(ref = ref),file = "Nearest-Promoter-for-flybase-batch-download.csv")

## NM2FBgn.txt 是nearest refseq 到 flybase转换后的结果
res.anno <- read.delim("./NM2FBgn.txt",stringsAsFactors = FALSE)
res.anno$X.SUBMITTED.ID[duplicated(res.anno$X.SUBMITTED.ID)]
res$Nearest.Refseq_FBgn <- res.anno$FBID_KEY[match(res$Nearest.Refseq, res.anno$X.SUBMITTED.ID)]
res$Nearest.Refseq_Symbol <- res.anno$SYMBOL[match(res$Nearest.Refseq, res.anno$X.SUBMITTED.ID)]

## Promoter-NM2-FBgn.txt 是 nearest promoter 到 flybase转换后的结果
res.anno <- read.delim("./Promoter-NM2-FBgn.txt",stringsAsFactors = FALSE)
res.anno$X.SUBMITTED.ID[duplicated(res.anno$X.SUBMITTED.ID)]
res$Nearest.PromoterID_FBgn <- res.anno$FBID_KEY[match(res$Nearest.PromoterID, res.anno$X.SUBMITTED.ID)]
res$Nearest.PromoterID_Symbol <- res.anno$SYMBOL[match(res$Nearest.PromoterID, res.anno$X.SUBMITTED.ID)]
head(res[,c("Nearest.PromoterID", "Nearest.Refseq","Nearest.Refseq_FBgn", "Nearest.Refseq_Symbol", "Nearest.PromoterID_FBgn", "Nearest.PromoterID_Symbol")], 100)

## 发现尽管nearest promoter列与nearest refseq列不一样的，但是最终对应的flybase fbgn id是全部一样的！
table(res$Nearest.PromoterID == res$Nearest.Refseq)
table(res$Nearest.PromoterID_FBgn == res$Nearest.Refseq_FBgn)

## 原数据提供了部分基因的FBgn和CG信息
genes <- res$FBgn
genes <- trimws(genes)
genes[genes == ""] <- res$Target.gene[genes == ""] # 有两个条目有两个FBgn id
# [1117] "FBgn0026063; FBgn0026064"
# [985] "FBgn0027363; FBgn0024227"
genes[985] <- "FBgn0027363"
genes[1117] <- "FBgn0026064"
write.csv(data.frame(ref = unique(genes[genes != ""])), file = "for-flybase-batch-download.csv")

# ref.txt： 原数据提供的基因提交到flybase转换后的结果
res.anno <- read.delim("./ref.txt",stringsAsFactors = FALSE)
res.anno <- res.anno[!duplicated(res.anno),]
duplicates <- res.anno$X.SUBMITTED.ID[duplicated(res.anno$X.SUBMITTED.ID)]
duplicates.rows <- res.anno[res.anno$X.SUBMITTED.ID %in% duplicates, ]
duplicates.rows[order(duplicates.rows$X.SUBMITTED.ID),] # 有对应多条FBgn id的基因
# X.SUBMITTED.ID    FBID_KEY   SYMBOL
# 34         CG10924 FBgn0003067   Pepck1
# 35         CG10924 FBgn0034356   Pepck2
# 129            PCB FBgn0027580      Pcb
# 130            PCB FBgn0286227 Hsap\\PC
res.anno[c("34","130"),] # 手动排查去除这两行
res.anno <- res.anno[!rownames(res.anno) %in% c("34","130"),]

# 将最新的注释信息加到结果中
res$Old <- res$FBgn
res$Old <- trimws(res$Old)
res$Old[res$Old == ""] <- res$Target.gene[res$Old == ""]
res$Old <- trimws(res$Old)
transfer_name <- function(input, ref.old, ref.new){
  results <- c()
  for (i in input) {
    old.names <- trimws(unlist(strsplit(i,split = ";")))
    new.names <- ref.new[match(old.names, ref.old)]
   results <- append(results, paste(new.names,collapse = ";"))
  }
  return(results)
}

res$FBgn.latest <- transfer_name(input = res$Old, ref.old = res.anno$X.SUBMITTED.ID, ref.new = res.anno$FBID_KEY)
res$Symbol.latest <- transfer_name(input = res$Old, ref.old = res.anno$X.SUBMITTED.ID, ref.new = res.anno$SYMBOL)

# 保存最终结果
write.csv(res,file = "final.results.csv")
```



问题： 发现很多基因不对应【注释得到的与原有的注释不一样】。经检查发现，此类绝大多数是因为gRNA对的匹配区域位于两个或多个基因的上游！



## 3. 优先比对到promoter上，然后再比对到基因body上（推荐）

### 提取promoter和genebody区域为bed文件

```shell
# extract all transcript from gtf file for promoter analysis
cat dmel-all-r6.51.gtf | awk '$3=="mRNA"' > transcripts.gtf
# extract all genes from gtf file for gene body analysis
cat dmel-all-r6.51.gtf |awk '$3 == "gene"' > genes.gtf
# get chrosome size
samtools faidx dmel-all-chromosome-r6.51.fasta 
cat dmel-all-chromosome-r6.51.fasta.fai | cut -f1,2 > chrom.sizes
```

```R
# gtf2bed.R
promoter.len <- 1000
transcripts <- read.table("./transcripts.gtf",stringsAsFactors = FALSE)
chrom.sizes <- read.table("./chrom.sizes",stringsAsFactors = FALSE)
genes <- read.table("./genes.gtf",stringsAsFactors = FALSE)
########### for promoter
# bed start from 0
transcripts$V4 <- transcripts$V4 -1
transcripts$V5 <- transcripts$V5 -1
table(transcripts$V7)
transcripts[!transcripts$V7 %in% c("-","+"),] # 从IGV人工检查一下基因，发现mod(mdg4)应该是位于负链
transcripts[!transcripts$V7 %in% c("-","+"),"V7"] <- "-"
transcripts[transcripts$V7== "-","V4"] <- transcripts[transcripts$V7== "-","V5"]
transcripts[transcripts$V7== "-","V5"] <- transcripts[transcripts$V7== "-","V5"] + promoter.len
transcripts[transcripts$V7== "+","V5"] <- transcripts[transcripts$V7== "+","V4"]
transcripts[transcripts$V7== "+","V4"] <- transcripts[transcripts$V7== "+","V4"] - promoter.len
transcripts$name <- paste(transcripts$V10, transcripts$V13, transcripts$V16, sep = ";")
transcripts <- transcripts[,c("V1","V4","V5","name")]
transcripts$V4[transcripts$V4 < 0] <- 0
size.mapping <- chrom.sizes$V2
names(size.mapping) <- chrom.sizes$V1
transcripts$V5 <- pmin(transcripts$V5,size.mapping[transcripts$V1]) # parallel minimal value of two or more given vectors
write.table(transcripts, file = "promoter.1kb.bed",row.names = FALSE,col.names = FALSE,quote = FALSE,sep = "\t")
############## for gene body
genes$name <- paste(genes$V10, genes$V13, sep = ";")
genes <- genes[,c("V1","V4","V5","name")]
write.table(genes, file = "genebody.bed",row.names = FALSE,col.names = FALSE,quote = FALSE,sep = "\t")
```

### 清华果蝇库数据清洗，并输出正反向gRNA序列

```R
df <- read.csv("./fly-stock.csv",stringsAsFactors = FALSE)
anyDuplicated(df$CV.number)
df$primer_sequence_forward
df$Sequence_forward <- gsub("TATATAGGAAAGATATCCGGGTGAACTTCG","", df$primer_sequence_forward)
df$Sequence_forward <- gsub("GTTTTAGAGCTAGAAATAGCAAG", "", df$Sequence_forward)
df$Sequence_forward <- unname(sapply(df$Sequence_forward, function(x){
  if (nchar(x) > 20) {
    gsub("^ttcg", "", x)
  }else{
    x
  }
}))
df$Sequence_forward <- unname(sapply(df$Sequence_forward, function(x){
  if (nchar(x) > 20) {
    gsub("^GTCG", "", x)
  }else{
    x
  }
}))
table(nchar(df$Sequence_forward))

df$Sequence_reverse <- gsub("ATTTTAACTTGCTATTTCTAGCTCTAAAAC","", df$primer_sequence_reverse)
df$Sequence_reverse <- gsub("CGACGTTAAATTGAAAATAGGTC", "", df$Sequence_reverse)
df$Sequence_reverse <- unname(sapply(df$Sequence_reverse, function(x){
  if (nchar(x) > 20) {
    gsub("^aaac", "", x)
  }else{
    x
  }
}))
df$Sequence_reverse <- unname(sapply(df$Sequence_reverse, function(x){
  if (nchar(x) > 20) {
    gsub("^AAAC", "", x)
  }else{
    x
  }
}))
table(nchar(df$Sequence_reverse))
write.csv(df,file = "fly-stock-clean-data.csv",row.names = FALSE)

file.create("input.forward.fa")
for (i in 1:nrow(df)) {
  cat("> ", df$CV.number[i], "\n", file = "input.forward.fa",append = TRUE)
  cat(df$Sequence_forward[i],"\n", file = "input.forward.fa",append = TRUE)
}

file.create("input.reverse.fa")
for (i in 1:nrow(df)) {
  cat("> ", df$CV.number[i], "\n", file = "input.reverse.fa",append = TRUE)
  cat(df$Sequence_reverse[i],"\n", file = "input.reverse.fa",append = TRUE)
}
```

### 基因组水平搜索gRNA序列的位置，并检查成对gRNA的匹配位置

```shell
# 去除空格
# cat input.forward.fa | sed 's/ //g' > input.forward2.fa
# cat input.reverse.fa | sed 's/ //g' > input.reverse2.fa
# 获取gRNA序列在基因组上的位置
seqkit locate -f input.forward.fa dmel-all-chromosome-r6.51.fasta > output.forward.location.txt
seqkit locate -f input.reverse.fa dmel-all-chromosome-r6.51.fasta > ouput.reverse.location.txt
```

```R
# check-gRNA.R
# 合并双向gRNA的基因组比对结果
forward <- read.delim("./output.forward.location.txt",stringsAsFactors = FALSE)
reverse <- read.delim("./output.reverse.location.txt",stringsAsFactors = FALSE)
forward$matchBy <- "forward"
reverse$matchBy <- "reverse"
table(duplicated(forward$patternName))

# res: a list,用于存储每对gRNA的比对结果
df <- read.csv("../fly-stock-clean-data.csv",stringsAsFactors = FALSE)
res <- list()
for (i in df$CV.number) {
  forward.sub <- forward[forward$patternName == i,]
  reverse.sub <- reverse[reverse$patternName == i,]
  res[[i]] <- rbind(forward.sub, reverse.sub)
}
table(unlist(lapply(res, function(x){nrow(x)})))

# 后面是对于每一对结果进行处理，理论上gRNA的比对位置不能离得太远！
# 第一步，得到的行数为0，或2，对于非0，2行的，进行人工矫正
res[which(!unlist(lapply(res, function(x){nrow(x)})) %in% c(0,2))]
## 手动去除一些错误的比对位置
res$CV221 <- res$CV221[-2, ]
res$CV236 <- res$CV236[-2, ]
res$CV347 <- res$CV347[-1, ]
res$CV753 <- res$CV753[c(1,4), ]
res$CV844 <- res$CV844[-1, ]
res$CV884 <- res$CV884[c(2,4), ]
res$CV1298 <- res$CV1298[c(1,3), ]
## 输出有多种比对可能的条目
unsolved <- c("CV1001","CV1015","CV1028","CV1051","CV1068","CV1197")
library(dplyr)
for (i in unsolved) {
  write.csv(arrange(res[[i]], seqID,start), file = paste0("unsolved-",i,".csv"))
}
res[unsolved] <- NULL # 去除这几个

# 第二步，行数为2的，检查是否forward和reverse各为1
res.matched <- res[unlist(lapply(res, function(x){nrow(x)})) == 2]
res.matched[!unlist(lapply(res.matched, function(x)all(c("forward", "reverse") %in% x$matchBy)))]
## 检查染色体是不是一样的
res.matched[!unlist(lapply(res.matched, function(x){x$seqID[1] == x$seqID[2]}))]

# 没有问题
# 第三步，检查长度
hist(unlist(lapply(res.matched, function(x){max(c(x$start,x$end)) - min(x$start,x$end)})))

# 把res.matched数据添加到df中
res.matched.arranged <- lapply(res.matched, function(x){
  data.frame(cvID = x$patternName[1],
             seqID = x$seqID[1],
             forward.pattern = x$pattern[x$matchBy == "forward"],
             forward.strand = x$strand[x$matchBy == "forward"],
             forward.start = x$start[x$matchBy == "forward"],
             forward.end = x$end[x$matchBy == "forward"],
             reverse.pattern = x$pattern[x$matchBy == "reverse"],
             reverse.strand = x$strand[x$matchBy == "reverse"],
             reverse.start = x$start[x$matchBy == "reverse"],
             reverse.end = x$end[x$matchBy == "reverse"],
             range.min = min(c(x$start, x$end)),
             range.max = max(c(x$start, x$end)),
             length = max(c(x$start, x$end)) - min(c(x$start, x$end)))
})
results <- Reduce(rbind, res.matched.arranged)
write.csv(results,file = "matched.results.csv")

df.location <- merge(df, results, by.x = "CV.number", by.y = "cvID",all.x = TRUE)
df.location <- df.location[match(df$CV.number, df.location$CV.number),]
write.csv(df.location,file = "../fly-stock-clean-data-with-location.csv")

# 输出bed文件，用于可视化gRNA对的位置和注释
temp <- df.location[!is.na(df.location$length),c("seqID","range.min","range.max","CV.number")]
temp$anno <- "fly-stock-cripsr"
temp$Strand <- "+"
temp$range.min <- temp$range.min - 1
temp$range.max <- temp$range.max - 1
write.table(temp, file = "targets.bed",row.names = FALSE,col.names = FALSE,quote = FALSE,sep = "\t")
```

### 分别使用promoter和genebody区域，注释成对gRNA组成的区域

```shell
conda activate py37
# https://bedtools.readthedocs.io/en/latest/content/tools/intersect.html
bedtools intersect -wao -a targets.bed -b promoter.1kb.bed  > promoter.intersect.txt
bedtools intersect -wao -a targets.bed -b genebody.bed  > genebody.intersect.txt
```
### 将注释结果添加到excel里

```R
## by promoter
mapping.promoter <- read.table("./promoter.intersect.txt",stringsAsFactors = FALSE)
table(mapping.promoter$V10 == ".") # 没有注释上的条目数
mapping.promoter <- mapping.promoter[mapping.promoter$V10 != ".",] # 去除没有比对上的条目，有可能是比对到了基因区域，而非promoter区域
mapping.promoter$overlaped_base_percentage <- mapping.promoter$V11/(mapping.promoter$V3 - mapping.promoter$V2)
head(mapping.promoter)
write.csv(mapping.promoter,file = "promoter比对详情.csv")

mapping.promoter.GeneLevel <- lapply(split(mapping.promoter, f = mapping.promoter$V4), function(x)paste0(unique(unlist(lapply(strsplit(x$V10,split = ";"),"[",1))),collapse = ";"))
promoter2FBgn <- unlist(mapping.promoter.GeneLevel)
mapping.promoter.GeneLevel <- lapply(split(mapping.promoter, f = mapping.promoter$V4), function(x)paste0(unique(unlist(lapply(strsplit(x$V10,split = ";"),"[",2))),collapse = ";"))
promoter2Symbol <- unlist(mapping.promoter.GeneLevel)

######## by gene body
mapping.genebody <- read.table("./genebody.intersect.txt",stringsAsFactors = FALSE)
table(mapping.genebody$V10 == ".") # 没有注释上的条目数
mapping.genebody <- mapping.genebody[mapping.genebody$V10 != ".",] # 去除没有比对上的条目，有可能是比对到了基因区域，而非promoter区域
mapping.genebody$overlaped_base_percentage <- mapping.genebody$V11/(mapping.genebody$V3 - mapping.genebody$V2)
head(mapping.genebody)
write.csv(mapping.genebody,file = "genebody比对详情.csv")

mapping.genebody.GeneLevel <- lapply(split(mapping.genebody, f = mapping.genebody$V4), function(x)paste0(unique(unlist(lapply(strsplit(x$V10,split = ";"),"[",1))),collapse = ";"))
genebody2FBgn <- unlist(mapping.genebody.GeneLevel)
mapping.genebody.GeneLevel <- lapply(split(mapping.genebody, f = mapping.genebody$V4), function(x)paste0(unique(unlist(lapply(strsplit(x$V10,split = ";"),"[",2))),collapse = ";"))
genebody2Symbol <- unlist(mapping.genebody.GeneLevel)

df <- read.csv("../fly-stock-clean-data-with-location.csv",stringsAsFactors = FALSE)
df <- df[,c(1:11,23)]
df$Target.Gene.FBgn.By.Promoter <- promoter2FBgn[df$CV.number]
df$Target.Gene.Symbol.By.Promoter <- promoter2Symbol[df$CV.number]
df$Target.Gene.FBgn.By.Genebody <- genebody2FBgn[df$CV.number]
df$Target.Gene.Symbol.By.Genebody <- genebody2Symbol[df$CV.number]
write.csv(df,file = "清华CV(CRISPR-activation)-fly-stock-注释结果.csv", row.names = FALSE)
```







