# **************** Define functions **************** #
## filter strings by keyword
filter_file_path_by_name_keyword <- function(PathVector,keywords){
  for (i in keywords) {
    PathVector <- PathVector[grepl(i,basename(PathVector))]
  }
  return(PathVector)
}

## finding the rds file
find_Rds <- function(MainDir = "each.sample.analysis.dir",
                     Rds.dir.keyword = "Rds.dir.keyword"){
  message(paste("Finding Rds file from",MainDir, "...",sep = " "))
  Analysis.dirs <- list.dirs(MainDir,recursive = FALSE,full.names = TRUE)
  ## finding the Rds dir
  Rds.dir <- filter_file_path_by_name_keyword(Analysis.dirs, Rds.dir.keyword)
  if (length(Rds.dir) == 1) {
    message("Found the Rds dir...")
    Analysis.Rds.path <- Rds.dir
    if (length(list.files(Analysis.Rds.path,pattern = "rds$")) == 1) {
      Analysis.Rds.file <- list.files(Analysis.Rds.path,pattern = "rds$",full.names = TRUE)
    }else{
      stop("Check the if the Rds file exist, and only one Rds file is allowed!")
    }
  }else{
    stop("Check the Rds dir and the keyword, and only one Rds dir is allowed!")
  }
  return(Analysis.Rds.file)
}

## finding the cluster markers file
find_marker <- function(MainDir = "each.sample.analysis.dir",
                        Markers.dir.keyword = c("upper.dir","lower.dir"),
                        Markers.file.keyword = "marker.file.keyword",
                        check.marker.file.format = FALSE # if FALSE, compatiable with old analysis naming format(only one file are allowed), 
                        # if TRUE, the markers files should be: ProjectName+Resolution+Method.csv(Can be multiple files).
                        ){
  message(paste("Finding marker file from",MainDir, "...",sep = " "))
  Analysis.dirs <- list.dirs(MainDir,recursive = FALSE,full.names = TRUE)
  # filter the upper dir
  upper.dirs <- filter_file_path_by_name_keyword(Analysis.dirs, Markers.dir.keyword[1])
  #Analysis.dirs.names <- basename(Analysis.dirs)
  if (length(upper.dirs) == 1) {
    message("Found the Supreme dir for markers file...")
    Markers.supreme.path <- upper.dirs
    Markers.supreme.path.dirs <- list.dirs(Markers.supreme.path,recursive = TRUE)
    # filter the lower dir
    lower.dirs <- filter_file_path_by_name_keyword(Markers.supreme.path.dirs, Markers.dir.keyword[2])
    if (length(lower.dirs) == 1) {
      Markers.dir <- lower.dirs
      csv.files <- list.files(Markers.dir,pattern = "csv$")
      # filter csv file by keyword
      csv.files <- filter_file_path_by_name_keyword(csv.files, Markers.file.keyword)
      if (check.marker.file.format) {
        # Filter files by format: two underlines
        csv.files <- csv.files[grepl(".*\\+.*\\+.*csv$",basename(csv.files))]
        if (length(csv.files) == 0) {
          warning("Check the if the marker file format!")
          warning("the marker file will be set to NA!")
          Markers.file <- NA
        }else{
          Markers.file <- paste(paste(Markers.dir, csv.files, sep = "/"),collapse = ";")
        }
      }else{
        if (length(csv.files) == 1) {
          Markers.file <- paste(Markers.dir, csv.files,sep = "/")
        }else{
          warning("Check the if the marker file exist by using the defined keywords and only only one marker file is allowed!")
          warning("the marker file will be set to NA!")
          Markers.file <- NA
        }
      }
    }else{
      warning("Check the marker file dir keywords!")
      warning("the marker file will be set to NA!")
      Markers.file <- NA
    }
    
    
  }else{
    warning("Check the marker dir and the keyword, and only one marker dir is allowed!")
    warning("the marker file is set to NA!")
    Markers.file <- NA
  }
  return(Markers.file)
}

# prepare the data list - from results of cellbender seurat scrubelt pipeline
prepare_data <- function(Analysis.Main.Dir = c("/home/xilab/jinz/Mouse-Large-Intestine-scRNAseq/Analysis/SingleSample/Cellbender-Seurat-Scrublet"),
                         Rds.dir.keyword = c("3.3-Clean-Cells-Rds"),
                         Markers.dir.keyword = list(c("Downstream-Analysis","ClusterMarkers")), # the first keyword is for the supreme dir, the second is for the direct upper dir of the marker file.
                         Markers.file.keyword = c("res.0.6"),
                         Check.marker.file.format = c(FALSE)){
  # Check parameters
  if(!all(c(is.vector(Analysis.Main.Dir), is.vector(Rds.dir.keyword), is.list(Markers.dir.keyword), is.vector(Markers.file.keyword), is.vector(Check.marker.file.format)))){ stop("Check the parameter types!")}
  if(!all(sapply(Markers.dir.keyword,function(x){if (is.vector(x) & length(x) == 2){return(TRUE)}else{return(FALSE)}}))){stop("Check the Markers.dir.keyword, each must be length = 2!")}
  if (any(c(length(Analysis.Main.Dir) != length(Rds.dir.keyword),
            length(Analysis.Main.Dir) != length(Markers.dir.keyword),
            length(Analysis.Main.Dir) != length(Markers.file.keyword),
            length(Analysis.Main.Dir) != length(Check.marker.file.format)))) {
    stop("Check the parameters length!")
  }
  ## Search sample Dirs in each main directory
  final.list <- list()
  for (i in 1:length(Analysis.Main.Dir)) {
    Analysis.Samples.Dir <- list.dirs(Analysis.Main.Dir[i],recursive = FALSE)
    print(Analysis.Samples.Dir)
    res <- list(as.data.frame(t(sapply(Analysis.Samples.Dir, function(x){
      project.name <- basename(x)
      analysis.dir <- x
      rds.file <- find_Rds(MainDir =  x, Rds.dir.keyword = Rds.dir.keyword[i])
      markers.file <- find_marker(MainDir = x, Markers.dir.keyword =  Markers.dir.keyword[[i]], Markers.file.keyword = Markers.file.keyword[i], check.marker.file.format = Check.marker.file.format[i])
      return(c("project.name" = project.name, "analysis.dir" = analysis.dir, "rds.file" = rds.file, "markers.file" = markers.file))
    },USE.NAMES = FALSE)),stringsAsFactors = FALSE))
    final.list <- append(final.list, res)
  }
  final.list <- Reduce(rbind,final.list)

  # data infos
  data.list <- data.frame("Name" = final.list$project.name,
                          "AnalysisDir" = final.list$analysis.dir,
                          "Rds" = final.list$rds.file,
                          "Project.Name" = final.list$project.name,
                          "Markers.File" = final.list$markers.file, # try list all marker files~. marker file nameing principle: project name + assay used + resolution + method used for calculate markers + .csv; files are sepatated by commas.
                          # such as: ProjectName_Resolution_Method.csv Method could be:"Default","MAST"...
                          "Cluster.Anno" = rep(NA,length(final.list$project.name)), # if NA, will use a default empty plot - "www/default.Cluster.Annotation.jpg"; else, you need offer a path to annotated cluster plot. Attention, when split is not null, will output a vlnplot.
                          stringsAsFactors = FALSE) 
  return(data.list)
}
