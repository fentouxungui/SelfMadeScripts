samples.dirs <- list.dirs(".",recursive = FALSE)
res.list <- list()
for (i in samples.dirs) {
  qc <- readRDS(paste0(i,"/QC.rds"))
  res <- as.data.frame(t(qc$parameters)[1,,drop = FALSE])
  res$cells.before <- unname(qc$cells.info[1])
  res$cells.after <- unname(qc$cells.info[2])
  rownames(res) <- basename(i)
  res.list[[basename(i)]] <- res
}
res <- Reduce(rbind, res.list)
write.csv(res,file = "QC.summary.csv")