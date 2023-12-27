total.looms <- grep(".loom$",list.files("../../CellRanger-Outputs/",recursive = TRUE),value = TRUE)
looms <- grep("E_SI",total.looms,value = TRUE)

total.rds <- grep("rds$",list.files("../SingleSample/",recursive = TRUE),value = TRUE)
grep("^E_SI",total.rds,value = TRUE)


suppressMessages(library(reticulate))
suppressMessages(library(Seurat))
suppressMessages(library(dplyr))
sc <- reticulate::import("scanpy")
# remove bad samples
bad.samples <- c("E_LI_D135","E_SI_D135","E_SI_D135_B4")
all.files <- c()
# import all files
for (loom in looms) {
  fileName <- unlist(lapply(strsplit(loom,split = "/"),"[",1))[1]
  if (!fileName %in% bad.samples) {
    all.files <- append(all.files,fileName)
    filePath <- paste("../../CellRanger-Outputs/",loom,sep = "")
    adata = sc$read_loom(filePath)
    print(adata)
    
    ## Setup the Seurat Object
    # Get the expression matrix
    exprs <- t(as.matrix(adata$X))
    colnames(exprs) <- adata$obs_names$to_list()
    rownames(exprs) <- adata$var_names$to_list()
    # Create the Seurat object
    pbmc <- CreateSeuratObject(exprs)
    print(pbmc)
    assign(fileName,pbmc)
    rm(adata,pbmc,fileName,exprs,filePath)
  }
}

# QC and merge all
for (seurat.file in all.files[1]) {
  seurat <- get(seurat.file)
  print(seurat)
  seurat[["percent.mt"]] <- PercentageFeatureSet(seurat, pattern = "^mt-")
  # Visualize QC metrics as a violin plot
  print(VlnPlot(seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size = 0.001))
}















