sample.names <- basename(list.dirs("../../../Cellranger-Outputs",recursive = FALSE))
for (i in sample.names) {
  if (!dir.exists(i)) {
    dir.create(paste("./",i,sep = ""))
  }
}

