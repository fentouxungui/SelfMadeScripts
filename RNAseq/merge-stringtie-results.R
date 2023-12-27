#!/home/xilab/miniconda3/bin/Rscript
# https://github.com/fentouxungui
# zhangyongchao@nibs.ac.cn
# NIBS
# 2022-07-25
# 在基因水平和isform水平，合并stringtie TPM/FPKM outputs from gene_abundence.tab
# usage: 
# Rscript merge-stringtie-results.R \
# -p ./results/6_Counts_StringTie_OnlyKnownTranscripts


suppressPackageStartupMessages(require(optparse)) # don't say "Loading required package: optparse"

option_list = list(
  make_option(c("-p", "--path"), action="store", default='./results/6_Counts_StringTie_OnlyKnownTranscripts', type='character',
              help="请输入RSEM结果输出目录！[默认为 ./results/6_Counts_RSEM/Single]")
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

message(">>> merge results at Gene level...")
prefix <- "StringTie_Gene_"
sample.dirs <- list.dirs(opt$path,full.names = TRUE)[-1]
# check if all file have same column - gene id
message("> Check Gene Names...")
for (sample.dir in sample.dirs) {
  sampleName <- paste0(sample.dir,"/","gene_abundence.tab")
  temp <- read.delim(sampleName,header = TRUE,stringsAsFactors = FALSE)
  if (any(duplicated(temp$Gene.ID))) {
    warnings("Duplicated gene IDs found in gene_abundence.tab!")
    print(temp[temp$Gene.ID %in% temp$Gene.ID[duplicated(temp$Gene.ID)],])
    message("Please check the terms carefully! and last duplicates will be removed!")
    temp <- temp[!duplicated(temp$Gene.ID),]
  }
  
  if (match(sample.dir,sample.dirs) == 1) {
    gene.id <<- temp$Gene.ID
  }else{
    if(all(temp$Gene.ID == gene.id)){
      message(paste0("基因id一致性检查通过：",sampleName))
    }else if(all(temp$Gene.ID %in% gene.id) & all(gene.id %in% temp$Gene.ID )){
      message(paste0("基因id check passed with different order：",sampleName))
    }else{
      stop("check the gene name!")
    }
  }
}

# merge
message(">>> Start to merge at Gene level...")
for (key in c("TPM", "FPKM")) {
  for (sample.dir in sample.dirs) {
    sampleName <- paste0(sample.dir,"/","gene_abundence.tab")
    temp <- read.delim(sampleName,header = TRUE,stringsAsFactors = FALSE)
    temp <- temp[!duplicated(temp$Gene.ID),]
    if ( sample.dir == sample.dirs[1] ) {
      temp <- temp[c("Gene.ID", "Gene.Name", "Reference", "Strand", "Start",   "End" , key)]
      colnames(temp)[7] <- basename(sample.dir)
    } else {
      temp <- temp[c("Gene.ID",key)]
      colnames(temp)[2] <- basename(sample.dir)
    }
    
    if ( sample.dir == sample.dirs[1] ) { 
      FPKM.res <- temp 
    } else {
      FPKM.res <- merge(FPKM.res,temp,by = "Gene.ID",all = TRUE)
    }
  }
  write.csv(FPKM.res,file = paste0(opt$path,"/", prefix, key, ".csv"),quote = FALSE,row.names = FALSE)
  message(paste0("Gene Level - Merge ",key," done!"))
}

message("All Done!")