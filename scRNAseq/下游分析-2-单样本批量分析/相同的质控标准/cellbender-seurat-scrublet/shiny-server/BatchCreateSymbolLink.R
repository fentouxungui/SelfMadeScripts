Pipeline.dir <- paste("../Analysis-Pipeline",basename(getwd()),sep = "-")

if(!dir.exists(Pipeline.dir)){
  dir.create(Pipeline.dir)
}else{
  system(paste("rm -rf",Pipeline.dir,sep = " "))
}


# 1. CellRanger and Cellbender Report
save.dir <- "1.CellRanger-And-Cellbender-Report"

cellranger.dir <- "/home/xilab/jinz/Mouse-Large-Intestine-scRNAseq/CellRanger-Outputs"
cellranger.samples <- basename(list.dirs(cellranger.dir,recursive = FALSE))

cellbender.filename <- "cellbender.output.pdf"
cellranger.filename <- "web_summary.html"

ln.cellbender.from <- paste(cellranger.dir,cellranger.samples,"outs",cellbender.filename,sep = "/")
ln.cellranger.from <- paste(cellranger.dir,cellranger.samples,"outs",cellranger.filename,sep = "/")
ln.to <- paste(Pipeline.dir,save.dir,cellranger.samples,sep = "/")

for (i in 1:length(ln.cellbender.from)) {
  system(paste("mkdir -p",ln.to[i],sep = " "))
  ln.command.cellbender <- paste("ln -s",ln.cellbender.from[i],ln.to[i],sep = " ")
  ln.command.cellranger <- paste("ln -s",ln.cellranger.from[i],ln.to[i],sep = " ")
  # print(ln.command.cellranger)
  # print(ln.command.cellranger)
  system(ln.command.cellbender)
  system(ln.command.cellranger)
}


# 2. Downstream Analysis report
save.dir <- "2.Downstream-analysis-report"
system(paste("mkdir -p ",Pipeline.dir,"/",save.dir,sep = ""))

analysis.dir <- "/home/xilab/jinz/Mouse-Large-Intestine-scRNAseq/Analysis/SingleSample"

html.files <- grep("html$",list.files(analysis.dir,recursive = FALSE),value = TRUE) # attention! only top levels files will be found!

for ( htmlfile in html.files) {
  
  ln.command.from <- paste(analysis.dir, htmlfile,sep = "/")
  ln.command.to <- paste(Pipeline.dir,save.dir,htmlfile,sep = "/")
  
  # print(paste("ln -s",ln.command.from,ln.command.to))
  system(paste("ln -s",ln.command.from,ln.command.to))
}

