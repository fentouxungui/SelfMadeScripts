#!/home/xilab/miniconda3/bin/Rscript
# https://github.com/fentouxungui
# zhangyongchao@nibs.ac.cn
# NIBS
# 2022-07-25
# 在基因水平和isform水平，合并RSEM counts/TPM/FPKM outputs，并依据GTF文件里的gene_id和gene_name列添加基因注释。
# usage: 
# Rscript merge-RSEM-results.R \
# -g /data0/reference/Genome/Human/GENCODE/GRCh38.p13.release42/gencode.v42.primary_assembly.annotation.gtf \
# -p ./results/6_Counts_RSEM


suppressPackageStartupMessages(require(optparse)) # don't say "Loading required package: optparse"

option_list = list(
  make_option(c("-p", "--path"), action="store", default='./results/6_Counts_RSEM/Single', type='character',
              help="请输入RSEM结果输出目录！[默认为 ./results/6_Counts_RSEM/Single]"),
  make_option(c("-g", "--gtf"), action="store", default='', type='character',
              help="请输入RNAseq分析对应的gtf文件！"),
  make_option(c("-o", "--out"), action="store", default='', type='character',
              help="请输入保存位置，默认与path一致！")
)
opt = parse_args(OptionParser(option_list=option_list))


# 检查工作路径
if (opt$path == "./") {
  cat("检测消息：工作路径为默认路径，即当前目录!\n", file=stderr()) # print error messages to stderr
}else{
  cat("注意：工作路径为", opt$path, "!\n", file=stderr()) 
  if (!dir.exists(opt$path)) {
    cat("检测消息：工作路径不存在!\n", file=stderr()) 
  }
}
# 检查gtf文件
if (opt$gtf != "" & file.exists(opt$gtf)) {
  cat("检测消息：找到了gtf文件！\n", file=stderr()) # print error messages to stderr
}else{
  cat("注意：gtf文件未找到，或未设置！", opt$path, "!\n", file=stderr()) 
  stop("请检查gtf文件的设置！")
}
# 检查结果存储目录
if (opt$out == "") {
  opt$out <- opt$path
  cat("检测消息：结果存储目录未设定，默认为：", opt$out,"！\n", file=stderr()) # print error messages to stderr
}else{
  if(!dir.exists(opt$out)){dir.create(opt$out)}
  cat("检测消息：结果存储目录设定为：", opt$out,"！\n", file=stderr()) # print error messages to stderr
}

# Prepare Annotation
message(">>> Prepare gtf file for annotation!")
message(paste0("Using gtf file: ", opt$gtf))
suppressMessages(library(rtracklayer)) 
gtf <- rtracklayer::import(opt$gtf)



message(">>> merge results at Gene level...")
prefix <- "RSEM_Gene_"
gene.res <- grep("genes.results",list.files(opt$path,full.names = TRUE),value = TRUE)
# check if all file have same column - gene id
message("> Check Gene Names...")
for (gene.file in gene.res) {
  sampleName <- gsub(".genes.results","",basename(gene.file),fixed = TRUE)
  temp <- read.table(gene.file,header = TRUE,stringsAsFactors = FALSE)
  if (match(gene.file,gene.res) == 1) {
    gene.id <<- temp$gene_id
  }else{
    if(all(temp$gene_id == gene.id)){
      message(paste0("基因id一致性检查通过：",sampleName))
    }else{
      stop("Different Gene Order!")
    }
  }
}

gtf.gene <- as.data.frame(gtf)[,c("gene_id","gene_name")]
gtf.gene <- gtf.gene[!duplicated(gtf.gene$gene_id),]
mapping <- gtf.gene$gene_name
names(mapping) <- gtf.gene$gene_id

# merge
message(">>> Start to merge at Gene level...")
for (key in c("expected_count", "TPM", "FPKM")) {
  for (gene.file in gene.res) {
    sampleName <- gsub(".genes.results","",basename(gene.file),fixed = TRUE)
    if ( gene.file == gene.res[1] ) {
      temp <- read.table(gene.file,header = TRUE,stringsAsFactors = FALSE)
      temp <- temp[c("gene_id","length", "effective_length", key)]
      colnames(temp)[4] <- sampleName
    } else {
      temp <- read.table(gene.file,header = TRUE,stringsAsFactors = FALSE)
      temp <- temp[c("gene_id",key)]
      colnames(temp)[2] <- sampleName
    }
    
    if ( gene.file == gene.res[1] ) { 
      FPKM.res <- temp 
    } else {
      FPKM.res <- merge(FPKM.res,temp,by = "gene_id",all = TRUE)
    }
  }
  
  # head(FPKM.res)
  res <- cbind(FPKM.res[,1,drop = FALSE],
               data.frame(GeneName = mapping[FPKM.res$gene_id]), 
               FPKM.res[,2:length(colnames(FPKM.res))])
  
  write.csv(res,file = paste0(opt$out,"/", prefix, "Annotated_", key, ".csv"),quote = FALSE,row.names = FALSE)
  message(paste0("Gene Level - Merge ",key," done!"))
}

message(">>> merge results at Isform level...")
prefix <- "RSEM_Isform_"
gene.res <- grep("isoforms.results",list.files(opt$path,full.names = TRUE),value = TRUE)
# check if all file have same column - gene id
message("> Check Isform Names...")
# check if all file have same column - gene id
for (gene.file in gene.res) {
  sampleName <- gsub(".isoforms.results","",basename(gene.file),fixed = TRUE)
  temp <- read.table(gene.file,header = TRUE,stringsAsFactors = FALSE)
  if (match(gene.file,gene.res) == 1) {
    transcript.id <<- temp$transcript_id
  }else{
    if(all(temp$transcript.id ==transcript.id)){
      message(paste0("转录本id一致性检查通过：",sampleName))
    }else{
      stop("Different Gene Order!")
    }
  }
}


gtf.isform <- as.data.frame(gtf)[,c("transcript_id","gene_name")]
gtf.isform <- gtf.isform[!duplicated(gtf.isform$transcript_id),]
mapping <- gtf.isform$gene_name
names(mapping) <- gtf.isform$transcript_id


# merge
message(">>> Start to merge at Isform level...")
for (key in c("expected_count", "TPM", "FPKM")) {
  for (gene.file in gene.res) {
    sampleName <- gsub(".isoforms.results","",basename(gene.file),fixed = TRUE)
    if ( gene.file == gene.res[1] ) {
      temp <- read.table(gene.file,header = TRUE,stringsAsFactors = FALSE)
      temp <- temp[c("transcript_id","gene_id","length", "effective_length", key)]
      colnames(temp)[5] <- sampleName
    } else {
      temp <- read.table(gene.file,header = TRUE,stringsAsFactors = FALSE)
      temp <- temp[c("transcript_id",key)]
      colnames(temp)[2] <- sampleName
    }
    
    if ( gene.file == gene.res[1] ) { 
      FPKM.res <- temp 
    } else {
      FPKM.res <- merge(FPKM.res,temp,by = "transcript_id",all = TRUE)
    }
  }
  
  res <- cbind(FPKM.res[,1,drop = FALSE],
               data.frame(GeneName = mapping[FPKM.res$transcript_id]), 
               FPKM.res[,2:length(colnames(FPKM.res))])
  write.csv(res,file = paste0(opt$out,"/", prefix, "Annotated_", key, ".csv"),quote = FALSE,row.names = FALSE)
  message(paste0("Isform Level - Merge ",key," done!"))
}

message("All Done!")

