#!/home/xilab/miniconda3/bin/Rscript
# https://github.com/fentouxungui
# zhangyongchao@nibs.ac.cn
# NIBS
# 2022-07-25
# arrange the fastq database: 
# 1. renaming the fastq files to a unified format;
# 2. generate the md5sum file for each directory containning the fastq files.

# 1. Rename fastq files
# 1.1 标准名称
# 参考illumina测序仪下机FASTQ命名: https://support.illumina.com/help/BaseSpace_OLH_009008/Content/Source/Informatics/BS/NamingConvention_FASTQ-files-swBS.htm
# 例： SampleName_S1_L001_R1_001.fastq.gz
# 特别支持： SampleName_S1_L006_I1_002.fastq.gz
# SampleName  The sample name provided in the sample sheet. If a sample name is not provided, the file name includes the sample ID, which is a required field in the sample sheet and must be unique.
# S1  The sample number based on the order that samples are listed in the sample sheet starting with 1. In this example, S1 indicates that this sample is the first sample listed in the sample sheet.
# NOTE： Reads that cannot be assigned to any sample are written to a FASTQ file for sample number 0, and excluded from downstream analysis.
#	L001  The lane number.
# R1  The read. In this example, R1 means Read 1. For a paired-end run, there is at least one file with R2 in the file name for Read 2. When generated, index reads are I1 or I2.
# 001  The last segment is always 001.
# 1.2 支持的转换
# 按以下顺序逐步替换
# SampleName_S1_L001_R1_001.fq.gz # fq to fastq
# SampleName.1.fastq.gz and SampleName_1.fastq.gz and SampleName_1_fastq.gz # [._]1[._]fastq.gz to _R1.fastq.gz
# SampleName_S1_L001_R1.fastq.gz and  SampleName_S1_L001.R1.fastq.gz # [_.]R1[_.]fastq.gz to _R1_001.fastq.gz
# SampleName_R1_001.fastq.gz # _R1_001.fastq.gz to _S1_L001_R1_001.fastq.gz


# usage: ./prepare-fastqs.r
#        ./prepare-fastqs.r -p './' -v
#        ./prepare-fastqs.r --path './' -r


suppressPackageStartupMessages(require(optparse)) # don't say "Loading required package: optparse"

option_list = list(
  make_option(c("-p", "--path"), action="store", default='./', type='character',
              help="请输入您要整理的数据目录！"),
  make_option(c("-e", "--escaped-dirs"), action="store", default='scRNAseq,SpatialRNAseq', type='character',
              help="需要跳过的子数据目录，此路径下的fastq文件名将不会被修改！ [Default scRNAseq,SpatialRNAseq]。注意不适用于非标准目录[比如Data_From_Paper目录]!\n
              标准目录格式为：/data1/data_backup/CheMinhua/RNAseq/2023.08.07_RNAseq_PE_Fly-Gut-esg-Brat-IR-And-OE/BratIR-Rep1/Brat-IR1_S1_L001_R1_001.fastq.gz"),
  make_option(c("-r", "--run"), action="store_true", default=FALSE,type='logical',
              help="执行重命名"),
  make_option(c("-v", "--verbose"), action="store_true", default=FALSE,type='logical',
              help="输出更多信息 [默认不输出]")
)
opt = parse_args(OptionParser(option_list=option_list))

# 检查模式
# 1. 列出所有可被重命名的文件
# 2. 列出不可被重命名的文件
# 3. 列出跳过的目录

# 检查工作路径
if (opt$path == "./") {
  cat("检测消息：工作路径为默认路径，即当前目录!\n", file=stderr()) # print error messages to stderr
}else{
  cat("注意：工作路径为", opt$path, "!\n", file=stderr()) 
  if (!dir.exists(opt$path)) {
    cat("检测消息：工作路径不存在!\n", file=stderr()) 
  }
}
# 列出所有fastq文件
fastq.files <-  list.files(path = opt$path, all.files = FALSE, recursive = TRUE, pattern = "(*fq.gz$)|(*fastq.gz$)")
if (length(fastq.files) == 0) {
  stop("检测消息：错误！未发现以fq.gz或fastq.gz结尾的文件！")
}
cat(paste0("检测消息：共发现",length(fastq.files),"个样本以fq.gz或fastq.gz结尾！\n"), file=stderr()) # print error messages to stderr 
# 解析fastq文件名
format.standard <- ".*_S\\d+_L00[1-8]_(I1|R1|R2)_00[1-8]\\.fastq\\.gz" # 标准格式
res <- data.frame(fastq = basename(fastq.files), path = dirname(fastq.files),IsStandard = "",NewName = "",stringsAsFactors = FALSE)
# 去除指定目录下的文件
dirs.escaped <- trimws(unlist(strsplit(opt$'escaped-dirs',split = ",")))
category.dirs <- dirname(dirname(res$path))
category.dirs.escaped <- category.dirs[category.dirs %in% dirs.escaped]
if (length(category.dirs.escaped) != 0) {
  for (i in unique(category.dirs.escaped)) {
    cat(paste0("检测消息：被跳过的目录",paste0(opt$path,"/",i),"\n"), file=stderr()) # print error messages to stderr
  }
  cat(paste0("检测消息：总共",length(category.dirs.escaped), "个fastq文件被跳过！\n"), file=stderr()) # print error messages to stderr
  res <- res[!category.dirs %in% dirs.escaped,]
  cat(paste0("检测消息：最终",nrow(res),"个样本用于本次分析！\n"), file=stderr()) # print error messages to stderr 
}

## 判定是否为标准格式
res$IsStandard <- grepl(pattern = format.standard, res$fastq)
cat(paste0("检测消息：共发现",sum(!res$IsStandard),"个非标准格式的fastq样本名！\n"), file=stderr()) # print error messages to stderr
if (opt$verbose) { # 输出所有非标准格式的样本
  for (i in paste0(res$path, "/", res$fastq)[!res$IsStandard]) {
    cat(paste0(">>> ", i, "\n"), file=stderr()) # print error messages to stderr
  }
}
## 尝试解析非标准的格式
trans_name <- function(astring, format = format.standard){
  tmp <- gsub("fq\\.gz$","fastq.gz", astring)  # fq to fastq
  tmp <- gsub("\\.clean\\.fastq\\.gz$",".fastq.gz", tmp)  # .clean.fastq.gz to .fastq.gz
  if (!grepl(pattern = format, tmp)) {
    tmp <- gsub("[\\._]([12])[\\._]fastq.gz$","_R\\1.fastq.gz",tmp) # [._]1[._]fastq.gz to _R1.fastq.gz And [._]2[._]fastq.gz to _R2.fastq.gz
    if (!grepl(pattern = format, tmp)) {
      tmp <- gsub("[\\._](I1|R1|R2)[\\._]fastq.gz$","_\\1_001.fastq.gz",tmp) # [_.]R1[_.]fastq.gz to _R1_001.fastq.gz And [_.]R2[_.]fastq.gz to _R2_001.fastq.gz And [_.]I1[_.]fastq.gz to _I1_001.fastq.gz
      if (!grepl(pattern = format, tmp)) {
        tmp <- gsub("_(L00[1-8])_(I1|R1|R2)_(00[1-8]).fastq.gz$","_S1_\\1_\\2_\\3.fastq.gz",tmp) # _L001_R1_001.fastq.gz to _S1_L001_R1_001.fastq.gz And _L001_R1_002.fastq.gz to _S1_L001_R1_002.fastq.gz etc.
        if (!grepl(pattern = format, tmp)) {
          tmp <- gsub("_(S[1-96])_([I1|R1|R2)_(00[1-8]).fastq.gz$","_\\1_L001_\\2_\\3.fastq.gz",tmp) # _S1_R1_001.fastq.gz to _S1_L001_R1_001.fastq.gz etc.
          if (!grepl(pattern = format, tmp)) {
            tmp <- gsub("_(I1|R1|R2)_(00[1-8]).fastq.gz$","_S1_L001_\\1_\\2.fastq.gz",tmp) # _R1_001.fastq.gz to _S1_L001_R1_001.fastq.gz etc.
            if (!grepl(pattern = format, tmp)) {
              tmp <- gsub("\\.fastq.gz$","_S1_L001_R1_001.fastq.gz",tmp) # .fastq.gz to _S1_L001_R1_001.fastq.gz
            }
          }
        }
      }
    }
  }
  if (grepl(pattern = format, tmp)) {
    return(tmp)
  }else{
    return("")
  }
}
trans_name_vector <- function(vector, format = format.standard){ # 向量化trans_name函数
  res <- character(0)
  for (i in vector) {
    res <- append(res, trans_name(i, format = format))
  }
  return(res)
}
res$NewName[!res$IsStandard] <- trans_name_vector(res$fastq[!res$IsStandard])
if (sum(!res$IsStandard) == 0) {
  stop("检测消息：所有文件均为标准格式！")
}
res.solved <- res[!res$IsStandard,c("fastq","NewName")]
cat("检测消息：可以被解析的文件名有", nrow(res.solved),"个！\n", file=stderr()) 
cat("检测消息：请手动检查以下重命名是否正确！\n", file=stderr()) 
if (nrow(res.solved) != 0) {
  for (i in 1:nrow(res.solved)) {
    cat(paste0(res.solved$fastq[i], " >>>\n",res.solved$NewName[i],"\n\n"), file=stderr()) 
  }
}
res.sub <- res[!(res$IsStandard) & (res$NewName == ""),"fastq"] # 代码未测试
cat("检测消息：不可以被解析的文件名有", length(res.sub),"个！\n", file=stderr())
if (length(res.sub) != 0) {
  print(res.sub, file=stderr())
}
if (opt$run) {
  res.sub <- res[!(res$IsStandard) & (res$NewName != ""),]
  old.name <- paste0(res.sub$path,"/",res.sub$fastq)
  new.name <- paste0(res.sub$path,"/",res.sub$NewName)
  if (length(old.name) !=0) {
    file.rename(from = old.name, to = new.name)
  }
  message(paste0("执行消息：共计",length(old.name),"个文件名被修改！"))
}
message("代码运行完毕！")
