source("./ui-preparation-functions/prepare-data.R")
source("./ui-preparation-functions/prepare-html-reports.R")

Encrypted.App <- FALSE # is it a encrypted app? Public Data or Private Data
Universal.code <- TRUE # can app be opened with universal account
Datalist.mode <- "customized" # customized or automated

# Analysis.Main.Dir is needed for html pdf csv etc reports, can be set as a vector
Analysis.Main.Dir <- "/data2/shared_data_backup/qiannn/SpatialRNAseq/SpatialRNAseq_Mouse/Downstream-Analysis/Single-Sample"
# Analysis.Main.Dir <- c("/home/xilab/qiannn/scRNA_Mouse/DownstreamAnalysis/Merge-Total-Colon-Samples", "/data2/shared_data_backup/qiannn/scRNA_Mouse/DownstreamAnalysis/Merge-Total-Colon-Samples/Merge-All-Colon-Samples/Subset")

if (Datalist.mode == "automated") {
  # prepare data
  data.list <- prepare_data(Analysis.Main.Dir = Analysis.Main.Dir,
                            Rds.dir.keyword = c("(3.3-Clean-Cells-Rds)|(4.1-Replicates-Merged-Cells-Rds)"),
                            Markers.dir.keyword = list(c("Downstream-Analysis","(Cluster-Markers)|(ClusterMarkers)")), # the first keyword is for the supreme dir, the second is for the direct upper dir of the marker file.
                            Markers.file.keyword = c("res."),
                            Check.marker.file.format = FALSE)
  # Demo
  # data.list <- prepare_data(Analysis.Main.Dir = Analysis.Main.Dir,
  #                           Rds.dir.keyword = c("(3.3-Clean-Cells-Rds)|(4.1-Replicates-Merged-Cells-Rds)|(4.1-Merged-Rds)", "5.1-Subset-Rds"),
  #                           Markers.dir.keyword = list(c("Downstream-Analysis","(Cluster-Markers)|(ClusterMarkers)"), c("Downstream-Analysis","(Cluster-Markers)|(ClusterMarkers)")), # the first keyword is for the supreme dir, the second is for the direct upper dir of the marker file.
  #                           Markers.file.keyword = rep("res.", 2),
  #                           Check.marker.file.format = rep(FALSE, 2))
}else if (Datalist.mode == "customized") {
  rds.files <- list.files(Analysis.Main.Dir,pattern = "rds$",recursive = TRUE)
  data.list <- data.frame(Name =  unlist(lapply(strsplit(rds.files,split = "/"),"[",1)),
                          AnalysisDir =  paste0(Analysis.Main.Dir,"/",dirname(rds.files)),
                          Rds = paste0(Analysis.Main.Dir,"/",rds.files),
                          Project.Name =  unlist(lapply(strsplit(rds.files,split = "/"),"[",1)),
                          Markers.File = gsub("rds$","cluster.markers.csv",paste0(Analysis.Main.Dir,"/",rds.files)),
                          Cluster.Anno = rep(NA,length(rds.files)),stringsAsFactors = FALSE)
  # data.list <- data.frame(Name = c("Gut-10Genomics"),
  #                         AnalysisDir = Analysis.Main.Dir,
  #                         Rds = c(paste0(Analysis.Main.Dir,"/","s_fca_biohub_gut_10x.raw.seurat.rds")),
  #                         Project.Name = c("Gut-10Genomics"),
  #                         Markers.File = NA,
  #                         Cluster.Anno = rep(NA,1),stringsAsFactors = FALSE)
  # Name                                                                                                            AnalysisDir
  # E_LI_D115_B1 /home/xilab/jinz/scRNAseq/Mouse-Large-Intestine-scRNAseq/Analysis/SingleSample/Cellbender-Seurat-Scrublet/E_LI_D115_B1
  # Rds Project.Name
  # /home/xilab/jinz/scRNAseq/Mouse-Large-Intestine-scRNAseq/Analysis/SingleSample/Cellbender-Seurat-Scrublet/E_LI_D115_B1/3.3-Clean-Cells-Rds/E_LI_D115_B1.rds E_LI_D115_B1
  # Markers.File Cluster.Anno
  # /home/xilab/jinz/scRNAseq/Mouse-Large-Intestine-scRNAseq/Analysis/SingleSample/Cellbender-Seurat-Scrublet/E_LI_D115_B1/3.4-Downstream-Analysis/ClusterMarkers/E_LI_D115_B1.Cluster-Markers.CellBender_snn_res.0.6.csv           NA
}else{
  stop("Check the Datalist.mode, only allow customized or automated!")
}

# prepare or update html reports
# Demo: Mode:Merged suitable for a merged analysis and a simple analysis
prepare_reports(Analysis.Main.Dir = Analysis.Main.Dir, 
                Mode = "Merged") 
# Demo: For Single Sample
# prepare_reports(Analysis.Main.Dir = Analysis.Main.Dir, 
#                 Mode = "Single",
#                 Cellranger.Cellbender.Outputs.Dir = "/home/xilab/jinz/Mouse-Large-Intestine-scRNAseq/CellRanger-Outputs/")
Genes.qc <- c("nFeature_RNA", "percent.mt", "nCount_RNA", "Phase") # show in Note:
Markers.column.order <- c("cluster", "gene", "avg_logFC", "p_val_adj", "p_val", "pct.1", "pct.2")
Default.reduction.choice <- "umap" # in UPPER CASE
# Cluster.Annotation.grepl <- "(^RNA_snn)|(^Annotation)|(^CellBender)|(^integrated_)|(^annotation)|(^R_annotation)|(^leiden_res)"
Cluster.Annotation.grepl <- NA # using all factor columns as cluster resolution choices
Cluster.Annotation.Remove.grepl <- "(snn_)|(Annotation.)" # words not show in annotation choices
Default.Annotation.choice <- "Annotation.SCSA.Cluster.level"
# SplitBy.Choice.total <- c("By Sample" = "orig.ident",  # **Attention**: split choice should be factor! please check this before running the app.
#                           "By Type" = "type", 
#                           "By Subject" = "subject",
#                           "By Cell Cycle" = "Phase",
#                           "By Group" = "group") # Letters are not case-sensitive: type == Type; columns with two levels can used to generage splited Vlnplot
SplitBy.Choice.total <- NA
n_top_markers <- 25
SplitBy.levels.max <- 15 # guess factors to be splited according to the levels numbers 
# Attention: FeaturePlot in Seurat can not all customized group label! usually not necessary to add the label...

