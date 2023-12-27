single.dirs <- list.dirs("../../SingleSample/Cellbender-Seurat-Scrublet",recursive = FALSE)
day.dirs.table <- table(gsub("_B\\d+.*$","",basename(single.dirs)))
day.dirs <- names(day.dirs.table)[unname(day.dirs.table) > 1]
for (i in day.dirs) {
  if( !dir.exists(i) ){
    dir.create(i)
  }
}

# create ln -s for samples without replicates
day.dirs.single <- names(day.dirs.table)[unname(day.dirs.table) == 1]
ln.commnds <- paste("ln -s",single.dirs[gsub("_B\\d+.*$","",basename(single.dirs)) %in% day.dirs.single],"./")
for (i in ln.commnds) {
  system(i)
}
