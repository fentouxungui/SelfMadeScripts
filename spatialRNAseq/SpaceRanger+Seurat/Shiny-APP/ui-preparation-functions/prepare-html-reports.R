prepare_reports <- function(Analysis.Main.Dir = c("/home/xilab/qiannn/scRNA_Mouse/DownstreamAnalysis/Merge-Total-Colon-Samples"), 
                            # such as  "/home/xilab/qiannn/scRNA_Mouse/DownstreamAnalysis/Single-Sample/Cellbender-Seurat-Scrublet" for single samle analysis
                            # "/home/xilab/qiannn/scRNA_Mouse/DownstreamAnalysis/Merge-Total-Colon-Samples" for merged data.
                            # or multipe single analysis or multiple merged analysis, such as: 
                            # c("/home/xilab/qiannn/scRNA_Mouse/DownstreamAnalysis/Merge-Total-Colon-Samples", "/data2/shared_data_backup/qiannn/scRNA_Mouse/DownstreamAnalysis/Merge-Total-Colon-Samples/Merge-All-Colon-Samples/Subset")
                            Mode = "Merged", # Single or Merged, Attention, Merged and Single Analysis can not be put together.
                            Cellranger.Cellbender.Outputs.Dir = "/home/xilab/qiannn/scRNA_Mouse/Cellranger-Outputs/" # Only used for single mode
                            ){
  # remove old html reports
  Pipeline.dir <- paste("../Analysis-Pipeline",basename(getwd()),sep = "-")
  if(!dir.exists(Pipeline.dir)){
    dir.create(Pipeline.dir)
  }else{
    unlink(Pipeline.dir, recursive = TRUE)
    dir.create(Pipeline.dir)
  }
  # create bash file
  bash.command.file <- "prepare-reports.sh"
  file.create(bash.command.file)
  cat("#!/bin/bash\n",file = bash.command.file) # shebang
  # Define create links function
  create_link_command <- function(source.dir = "", target.dir = "", file.types = "(\\.html$)|(\\.tiff$)|(\\.csv$)|(\\.pdf$)|(\\.jpg$)|(\\.jpeg$)|(\\.png$)|(\\.bmp$)|(\\.svg$)", bash.command.file = "prepare-reports.sh"){
    # Downstream data analysis directorys (could include subset analysis)
    htmls <- list.files(source.dir, recursive = TRUE, pattern = file.types)
    # 3.1 create dirs
    dirs.new <- paste(target.dir, unique(dirname(htmls)[dirname(htmls) != "."]), sep = "/")  
    if (length(dirs.new) != 0) {
      for (i in dirs.new) { # create main sample directory
        cat(paste("mkdir -p", i, "\n"), file =  bash.command.file, append = TRUE)
      }
    }
    # 3.2 Create html lns
    ln.html.from <- paste(source.dir, htmls, sep = "/")
    ln.to <- paste(target.dir, htmls, sep = "/")
    for (i in 1:length(ln.html.from)) { # running too much time
      if (!dirname(ln.to[i]) %in% dirs.new) {
        cat(paste("mkdir -p", dirname(ln.to[i]), "\n"), file =  bash.command.file, append = TRUE) # create sub directory
      }
      ln.command.html <- paste("ln -s",ln.html.from[i], ln.to[i])
      cat(paste(ln.command.html, "\n"), file =  bash.command.file, append = TRUE)
    }
  }

  if (Mode == "Single") { # For Single Sample Analysis
    # 1. CellRanger and Cellbender Report
    save.dir <- "1.CellRanger-And-Cellbender-Report"
    cellranger.samples <- basename(list.dirs(Cellranger.Cellbender.Outputs.Dir,recursive = FALSE))
    cellbender.filename <- "cellbender.output.pdf"
    cellranger.filename <- "web_summary.html"
    ln.cellbender.from <- paste(Cellranger.Cellbender.Outputs.Dir,cellranger.samples,"outs",cellbender.filename,sep = "/")
    ln.cellranger.from <- paste(Cellranger.Cellbender.Outputs.Dir,cellranger.samples,"outs",cellranger.filename,sep = "/")
    ln.to <- paste(Pipeline.dir,save.dir,cellranger.samples,sep = "/")
    for (i in 1:length(ln.cellbender.from)) {
      cat(paste("mkdir -p",ln.to[i],"\n"), file = bash.command.file, append = TRUE)
      ln.command.cellbender <- paste("ln -s",ln.cellbender.from[i],ln.to[i])
      ln.command.cellranger <- paste("ln -s",ln.cellranger.from[i],ln.to[i])
      cat(paste(ln.command.cellbender,"\n"), file = bash.command.file, append = TRUE)
      cat(paste(ln.command.cellranger,"\n"), file = bash.command.file, append = TRUE)
    }
    # 2. Downstream Analysis report
    save.dir <- "2.Downstream-analysis-report"
    Downstream.Analysis.path <- paste(Pipeline.dir, save.dir, sep = "/")
    cat(paste("mkdir -p", Downstream.Analysis.path, "\n"), file = bash.command.file, append = TRUE)
    for (i in Analysis.Main.Dir) {
      create_link_command(source.dir = i, target.dir = Downstream.Analysis.path, file.types = "(\\.html$)|(\\.tiff$)|(\\.csv$)|(\\.pdf$)|(\\.jpg$)|(\\.jpeg$)|(\\.png$)|(\\.bmp$)|(\\.svg$)", bash.command.file = bash.command.file)
    }
  } else if(Mode == "Merged"){ # For Merged Analysis
    # Downstream Analysis report
    for (i in Analysis.Main.Dir) {
      create_link_command(source.dir = i, target.dir = Pipeline.dir, file.types = "(\\.html$)|(\\.tiff$)|(\\.csv$)|(\\.pdf$)|(\\.jpg$)|(\\.jpeg$)|(\\.png$)|(\\.bmp$)|(\\.svg$)", bash.command.file = bash.command.file)
    }
  }else{
    stop("Mode should be set as 'Single' or 'Merged'!")
  }
  # in case some files not exist, such as cellbender.output.pdf and other mistakes in creating links, please check the information carefully.
  cat(paste0('for file in $(find ', Pipeline.dir,' -type l); do if [ ! -e $file ]; then echo "rm $file"; rm -f $file; fi; done'), file = bash.command.file, append = TRUE)
  system(paste("sh", bash.command.file),wait = FALSE)
}

