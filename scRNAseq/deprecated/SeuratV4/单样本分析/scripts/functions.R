qc_filter <- function(obj,
                      nFeature_RNA.min = 300,
                      nFeature_RNA.max = +Inf,
                      nCount_RNA.min = 3000,
                      nCount_RNA.max = +Inf,
                      percent.mt.max = 20,
                      log10GenesPerUMI.min = 0.8){
  require(Seurat)
  require(SeuratObject)
  # create a data.frame to save parameters
  filter.df <- data.frame(row.names = c("nFeature_RNA.min","nFeature_RNA.max","nCount_RNA.min","nCount_RNA.max","percent.mt.max","log10GenesPerUMI.min"),
                          cutoff = c(nFeature_RNA.min, nFeature_RNA.max, nCount_RNA.min, nCount_RNA.max, percent.mt.max, log10GenesPerUMI.min),
                          good.cells = 0,
                          bad.cells = 0)
  for (i in rownames(filter.df)) {
    if (grepl("min$",i)) {
      filter.df[i, "good.cells"] <- sum(obj@meta.data[,gsub("\\.min$","",i)] > filter.df[i,"cutoff"])
      filter.df[i, "bad.cells"] <- sum(obj@meta.data[,gsub("\\.min$","",i)] <= filter.df[i,"cutoff"])
    }else{
      filter.df[i, "good.cells"] <- sum(obj@meta.data[,gsub("\\.max$","",i)] < filter.df[i,"cutoff"])
      filter.df[i, "bad.cells"] <- sum(obj@meta.data[,gsub("\\.max$","",i)] >= filter.df[i,"cutoff"])
    }
  }
  cells.before <- ncol(obj)
  obj <- subset(obj, nFeature_RNA > nFeature_RNA.min & nFeature_RNA < nFeature_RNA.max & nCount_RNA > nCount_RNA.min & nCount_RNA < nCount_RNA.max &
                          percent.mt < percent.mt.max & log10GenesPerUMI > log10GenesPerUMI.min)
  cells.after <- ncol(obj)
  cells.info <- c(cells.before = cells.before, cells.after = cells.after)
  saveRDS(list(parameters = filter.df, cells.info = cells.info), file = "QC.rds")
  return(obj)
}