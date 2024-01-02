# 1. 自动化为每个样本拆创建分析目录
Cellranger_output_dir <- "/home/xilab/Data_Backup/Data_From_Paper/2023.03.15_DevelopmentalCell_PMID36924771_Smooth-muscle-contributes-to-the-development-and-function-of-a-layered-intestinal-stem-cell-niche/GSE184158/Supplementary-file/cellranger-outputs"
samples_dirs <- list.dirs(Cellranger_output_dir,recursive = FALSE)
# file_postfix <- "/outs/filtered_feature_bc_matrix"
file_postfix <- "/filtered_feature_bc_matrix"
for (i in samples_dirs ) {
  dir.create(basename(i))
  cellranger_mat_path <- paste0(i,file_postfix)
  print(cellranger_mat_path)
  saveRDS(cellranger_mat_path, file = paste0(basename(i),"/CellRanger_Matrix_path.rds"))
}


# 2. 自动化为每个目录拷贝分析代码
## 注意，所有代码都放在同一层中，不可以建立子目录
scripts.dir <- "/home/xilab/reference/Scripts/SelfMadeScripts/scRNAseq/CellRanger-SeuratV4/单样本水平批量分析/scripts"
scripts.pattern <- "(Rmd$)|(ipynb$)|(\\.R$)"
scripts <- list.files(path = scripts.dir, pattern = scripts.pattern,full.names = TRUE)
for (i in basename(samples_dirs)) {
  file.copy(from = scripts, to = i)
}
