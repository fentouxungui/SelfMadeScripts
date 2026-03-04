sample.paths <- list.dirs(c("/data2/shared_data_backup/jinz/Mouse-Small-Intestine-scRNAseq/Analysis/SingleSample/Cellbender-Seurat-Scrublet",
                            "/data2/shared_data_backup/jinz/Mouse-Large-Intestine-scRNAseq/Analysis/SingleSample/Cellbender-Seurat-Scrublet",
                            "/data3/Xilab-Data-Analysis/paper-data-analysis/2021.08.03_CellReports_PMID34348153_Mouse-Embryo-Colon-Development-Damage-E145-E155-E185-DSS-scRNAseq/scRNAseq/from-fastq/Analysis/SingleSample/Cellbender-Seurat-Scrublet"),
                          recursive = FALSE)
samples.used <- list.files("../../CellRanger-Outputs/", recursive = FALSE)
sample.used.paths <- sample.paths[basename(sample.paths) %in% samples.used]
length(samples.used)
length(sample.paths)
length(sample.used.paths)

for (i in sample.used.paths) {
  if (!dir.exists(basename(i))) {
    ln.command <- paste("ln -s", i, basename(i))
    print(ln.command)
    system(ln.command)
  }
}

