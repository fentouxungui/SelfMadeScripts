library(shiny)
library(Seurat)
library(ggplot2)
library(RColorBrewer)
library(dplyr)
library(DT)
library(ggpubr)
library(stringr)
library(patchwork)

# load customized functions
source("./server-functions/functions.R")

# define some variables
cache.rds.list <- list() # save the imported rds file
cache.markers.list <- list() # save the imported markers file path

server <- function(input, output, session) {
  
  if (Encrypted.App) {
    # credentials infos
    credentials <- extract_credentials()
    # call the server part
    # check_credentials returns a function to authenticate users
    # timeout:  Timeout session (minutes) before logout if sleeping, Defaut to 15. 0 to disable.
    res_auth <- secure_server( check_credentials = check_credentials(credentials), timeout = 180)
    output$auth_output <- renderPrint({ reactiveValuesToList(res_auth) })
  }

  # Data list show in UI: Analysis Report
  output$UIDataList <- DT::renderDataTable(DT::datatable(data.list,
                                                         class = 'cell-border stripe',
                                                         caption = htmltools::tags$caption(
                                                           style = 'caption-side: bottom; text-align: center;',
                                                           'Table 1: ', htmltools::em('Files used for running this App.')
                                                         ),
                                                         options = list(searching = FALSE, scrollX = T)))
  
  
  # define function: import the Rds file from path or list
  dataset.input <- reactive({
    message("Data: Loading scRNAseq Data - ", input$Dataset, " ...")
    if (is.null(names(cache.rds.list)) | ! input$Dataset %in% names(cache.rds.list)) {
      Rds.raw <- readRDS(data.list$Rds[match(input$Dataset, data.list$Project.Name)])
      Rds.raw@project.name <- input$Dataset
      # change some columns in meta.data to factors
      columns.to.factor <- colnames(Rds.raw@meta.data)[check_df_level(Rds.raw@meta.data) & !check_df_factor(Rds.raw@meta.data)]
      Rds.raw@meta.data[columns.to.factor] <- lapply(Rds.raw@meta.data[columns.to.factor],as.factor)
      cache.rds.list[[input$Dataset]] <<- Rds.raw
    }else{
      Rds.raw <- cache.rds.list[[input$Dataset]]
    }
    return(Rds.raw)
  })
 
  
  # define function: import the markers file from path or list
  markers.path.input <- reactive({
    message("Data: Loading the Markers file of dataset - ", input$Dataset, " ...")
    if (is.null(names(cache.markers.list)) | ! input$Dataset %in% names(cache.markers.list)) {
      marker.info <- data.list$Markers.File[match(input$Dataset, data.list$Project.Name)]
      if (is.na(marker.info)) {
        markers.path <- NA
      }else{
        markers.path <- unlist(strsplit(marker.info, split = ";"))
      }
      cache.markers.list[[input$Dataset]] <<- markers.path
    }else{
      markers.path <- cache.markers.list[[input$Dataset]]
    }
    return(markers.path)
  })
  
  # define actionbutton - OpenLinkedDimplot UI
  output$LinkedDimplot.UI <- renderUI({
    system(paste0("mkdir -p ../LinkedDimplot-SingleSample/", input$Dataset))
    system(paste0("cp ./LinkedDimplotApp/app.R ../LinkedDimplot-SingleSample/", input$Dataset))
    app.path <- paste0("../LinkedDimplot-SingleSample/", input$Dataset,"/app.R")
    x <- readLines(app.path)
    x[8] <- paste0('object <- readRDS("', data.list$Rds[match(input$Dataset, data.list$Project.Name)], '")')
    cat(x, file = app.path, sep = "\n")
    actionButton(inputId='OpenLinkedDimplot', label="Check Linked Dimplot", icon = icon("th"), onclick = paste("window.open('",
                                                                                                               # "window.open('http://192.168.13.45",
                                                                                                               gsub("/data3/Xilab-Data-Analysis/xilab/","",dirname(getwd())),
                                                                                                               paste0("LinkedDimplot-SingleSample/", input$Dataset ,"')"),
                                                                                                               sep = "/"))
  })
  
  
  
  
  # define markers choices UI
  output$MarkerList.UI <- renderUI({
    message("UI: Loading Markers List UI...")
    if (all(is.na(markers.path.input()))) {
      # no markers path
      markers.choices <- "Not Found"
    } else {
      markers.choices <- markers.path.input()
      if (length(markers.choices) > 1 ) {
        names(markers.choices) <- unlist(lapply(lapply(strsplit(gsub(".csv", "", basename(markers.choices)), split = "\\+"), "[", c(2,3)), function(x){paste(x, collapse = "_")}))
      } else if(grepl("\\+", markers.choices)) {
        names(markers.choices) <- unlist(lapply(lapply(strsplit(gsub(".csv", "", basename(markers.choices)), split = "\\+"), "[", c(2,3)), function(x){paste(x, collapse = "_")}))
      }else{
        names(markers.choices) <- gsub(".csv", "", basename(markers.choices))
      }
    }
    selectInput("MarkerList","Cluster Marker Choices:", choices = markers.choices)
  })
  
  # define reductions choices UI
  output$Reductions.UI <- renderUI({
    message("UI: Loading Reductions Choice UI...")
    reduction.choice <- rev(Reductions(dataset.input()))
    names(reduction.choice) <- toupper(reduction.choice)
    # Remove PCA from reduction choices
    reduction.choice <- reduction.choice[names(reduction.choice) != "PCA"]
    if (Default.reduction.choice %in% unname(reduction.choice)) { # in case the customized reduction choice not exist!
      selectInput("DimensionReduction", "Dimension Reduction:", choices = reduction.choice, selected = Default.reduction.choice) # set default reduction
    }else{
      selectInput("DimensionReduction", "Dimension Reduction:", choices = reduction.choice) # set default reduction
    }
  })
  
  # define Cluster Annotation choice
  output$ClusterResolution.UI <- renderUI({
    message("UI: Loading Resolution Choice UI...")
    # without setting cluster resolution: using all factor columns
    if (is.na(Cluster.Annotation.grepl)) {
      full.res <- colnames(dataset.input()@meta.data)[check_df_factor(dataset.input()@meta.data)]
      names(full.res) <- full.res
      selectInput("ClusterResolution","Cluster Resolution:", choices = full.res)
    }else{ # using customized resolutions
      full.res <- rev(grep(Cluster.Annotation.grepl, colnames(dataset.input()@meta.data), value = TRUE))
      names(full.res) <- gsub(Cluster.Annotation.Remove.grepl, "", full.res)
      if (Default.Annotation.choice %in% unname(full.res)) { # in case the customized Annotation choice not exist!
        selectInput("ClusterResolution","Cluster Resolution:", choices = full.res, selected = Default.Annotation.choice)
      }else{
        selectInput("ClusterResolution","Cluster Resolution:", choices = full.res)
      }
    }
  })
  
  # return the markers of each cluster
  ClusterMarkers <- reactive({
    message("Data: Loading the selected Marker file...")
    if (length(is.na(input$MarkerList)) == 0) { # when app initiation, length(is.na(input$MarkerList)) = 0
      message("Data: Marker file Loaded!")
      setNames(data.frame(matrix(ncol = length(Markers.column.order), nrow = 0)), Markers.column.order) # Return an empty data.frame
    }else if (is.na(input$MarkerList)) { # for app initiation
      message("Data: Marker file Loaded!")
      setNames(data.frame(matrix(ncol = length(Markers.column.order), nrow = 0)), Markers.column.order) # Return an empty data.frame
    }else if(input$MarkerList != "Not Found"){
      message("Data: Marker file Loaded!")
      read.csv(input$MarkerList, stringsAsFactors = FALSE, row.names = 1)
    }else{
      message("Data: Marker file Loaded!")
      setNames(data.frame(matrix(ncol = length(Markers.column.order), nrow = 0)), Markers.column.order) # Return an empty data.frame
    }
  })
  

  #************** Cluster Annotation Plot and Vlnplot (when two levels) ********************#
  # define all available features
  GeneLibrary <- reactive({
    message("Action-1: Generating Features Library...")
    c(rownames(dataset.input()), Genes.qc)
  })
  # Check the input gene
  Gene.Revised <- reactive({
    message("Action-2: Revising the input Features...")
    CheckGene(InputGene = ifelse(is.na(input$GeneSymbol), NA, input$GeneSymbol), GeneLibrary = GeneLibrary())
  })
  # # slow down render speed
  # Gene.Revised <- Gene.Revised.init %>% debounce(2000)
  
  # Generate cluster annotation plot / Vlnplot (when two levels)
  PlotClusterAnnotation <- reactive({
    message("Action-4: Generating the Cluster Annotation Plot / Vlnpot...")
    if (!is.na(Gene.Revised())) { 
        print(paste("{ Dataset:",input$Dataset, "}; { Action: Vlnplot }; { Feature:", Gene.Revised(), "}; { Resolution: ", input$ClusterResolution, "}",sep = " "))
        return(VlnPlot(dataset.input(), features = Gene.Revised(), pt.size = 0, group.by = input$ClusterResolution))
    }else{
      Annotation.figure.path <- data.list$Cluster.Anno[match(input$Dataset, data.list$Project.Name)]
      if (is.na(Annotation.figure.path)) {
        list(src = "./www/default.Empty.Cluster.Annotation.jpg",
             alt = "Cluster Annotation")
      }else{
        list(src = Annotation.figure.path,
             alt = "Cluster Annotation")
      }
    }
  })

  output$ClusterAnnotationPlot <- renderPlot({PlotClusterAnnotation()}, height = function() {session$clientData$output_ClusterAnnotationPlot_width * input$ClusterAnnotationPlotHWRatio}) # set plot height = width

  
  #************** Scatter Plot ********************#
  # Scatter Plot
  PlotFeature <- reactive({
    message("Action-4: Generating the Scatter Plot...")
    if(is.na(Gene.Revised())) { # Dimplot, Attention: wrong input gene will also output the dimplot!
        print(paste("{ Dataset:",input$Dataset, "}; { Action: DimPlot }; { Reduction:", input$DimensionReduction, "}; { Resolution: ", input$ClusterResolution, "}",sep = " "))
        p1 <- DimPlot(dataset.input(), reduction = input$DimensionReduction, label = input$ShowLabel, pt.size = input$PointSize, label.size = input$LabelSize, group.by = input$ClusterResolution)
        p2 <- SpatialDimPlot(dataset.input(), label = input$ShowLabel, label.size = input$LabelSize, group.by = input$ClusterResolution, pt.size.factor = input$SpatialPointSize)
        p1 + p2
    }else { # Feature plot
      temp.rds <- dataset.input()
      Idents(temp.rds) <- input$ClusterResolution
      # Gene.Revised()[1]: could set only use the first gene (input DL, get Dl and dl)
      print(paste("{ Dataset:",input$Dataset, "}; { Action: FeaturePlot }; { Reduction:", input$DimensionReduction, "}; { Feature: ", Gene.Revised(), "}",sep = " "))
      p1 <- FeaturePlot(temp.rds, features = Gene.Revised(), pt.size = input$PointSize, label = input$ShowLabel, label.size = input$LabelSize, reduction = input$DimensionReduction) + NoLegend()
      p2 <- SpatialFeaturePlot(temp.rds, features = Gene.Revised(), pt.size.factor = input$SpatialPointSize, alpha = c(input$SpatialLowerTran, input$SpatialHighTran)) + theme(legend.position = "right")
      g2 <- ggplot_build(p2)
      g1 <- ggplot_build(p1)
      g1$data[[1]]["colour"] <- g2$data[[1]]["fill"]
      gtable <- ggplot_gtable(g1)
      wrap_plots(list(gtable, p2))
    }
  })
  
  output$FeaturePlot <- renderPlot({PlotFeature()}, height = function(){session$clientData$output_FeaturePlot_width * input$FeaturePlotHWRatio}) # box plot: height = width

  
  output$DownloadPlot <- downloadHandler( # Unsolved Problems: The Label size in Donwloaded pdf is larger than what we see in web.
    filename = function() { ifelse(is.na(Gene.Revised()), 
                                   paste(dataset.input()@project.name, "_Dimplot_", input$DimensionReduction, '.pdf', sep=''), 
                                   paste(dataset.input()@project.name, "_FeaturePlot_", input$DimensionReduction, "_", Gene.Revised(), '.pdf', sep=''))},
    content = function(filename) {ggsave(filename, plot = PlotFeature(), 
                                         # plot size are same with what client see in website, 1 pixel (X)	== 0.0104166667 inch, session$clientData$output_FeaturePlot_width is in pixel.
                                         width = session$clientData$output_FeaturePlot_width * 0.0104166667, height= ( session$clientData$output_FeaturePlot_width * 0.0104166667 )* input$FeaturePlotHWRatio, 
                                         units = "in", 
                                         device = "pdf", dpi = 300)
    })
  
  Markers <- eventReactive(input$CalculateDEG, {
    message("Action-4: Calculating DEGs...")
    if(Check.Group()$Check & !Check.Group()$Default) {
      showModal(modalDialog(title = "Calculating DEGs...", "Please wait for a few minutes!", footer= NULL, size = "l"))
      # ifelse(NA, NA, 1) will cause an error
      if(is.na(Check.Group()$G2)){ Ident.2 <- NULL } else { Ident.2 <- as.numeric(Check.Group()$G2)}
      Markers.DEG <- FindMarkers(dataset.input(), ident.1 = as.numeric(Check.Group()$G1), ident.2 = Ident.2, min.pct = 0.25, group.by = input$ClusterResolution)
      removeModal()
      return(Markers.DEG )
    } else{
      return(ClusterMarkers())
    }
  })
  
  observeEvent(input$CalculateDEG, {
    if(!is.na(input$CalculateDEG)) {
      removeUI("#DownloadDEG")}
  })
  
  
  Check.Group <- eventReactive(input$CalculateDEG, {
    message("Action-3: Checking DEGs Calculating Parameters...")
    # check parameters
    if (trimws(input$Group1) == "" & trimws(input$Group2) == "") {
      showModal(modalDialog(title = "Attention", "Both group are set to blank, Cluster Markers will be downloaded by default!", size = "l"))
      return(list(Check = TRUE, Default = TRUE, G1 = NA, G2 = NA))
    } else {
      g1 <- trimws(unlist(strsplit(as.character(input$Group1), ",")))
      g2 <- trimws(unlist(strsplit(as.character(input$Group2), ",")))
      g2 <- ifelse(length(g2) == 0, NA, g2) # why g2 == character(0), when g2 reset to "", will cause an error
      if (length(g1) == 0) { # g1 can not be empty
        showModal(modalDialog(title = "Error", "Group 1 can not be empty!", size = "l"))
        return(list(Check = FALSE, Default = TRUE, G1 = NA, G2 = NA))
      } else {
        groups.all <- levels(dataset.input()@meta.data[, input$ClusterResolution]) # check if the groups are right(in the selected resolution levels)
        if (!all(g1 %in% groups.all)) {
          showModal(modalDialog(title = "Error", "Please Check if Group 1 is in the selected Cluster Annotation!", size = "l"))
          return(list(Check = FALSE, Default = TRUE, G1 = NA, G2 = NA))
        }else if(!is.na(g2) & !all(g2 %in% groups.all)){
          showModal(modalDialog(title = "Error", "Please Check if Group 2 is in the selected Cluster Annotation!", size = "l"))
          return(list(Check = FALSE, Default = TRUE, G1 = NA, G2 = NA))
        }else{
          return(list(Check = TRUE, Default = FALSE, G1 = g1, G2 = g2))
        }
      }
    }
  })
  
  observeEvent(input$CalculateDEG, {
    if (Check.Group()$Check) {
      insertUI(selector = "#CalculateDEG", where = "afterEnd", ui = downloadButton("DownloadDEG",label = "Download DEGs!"))
    }
  })
  
  output$DownloadDEG <- downloadHandler(
    filename = function(){paste(dataset.input()@project.name, "_DEGs-of-cluster-", gsub(",", "-", input$Group1), "-Vs-", gsub(",", "-", input$Group2),"_Resolution-", input$ClusterResolution, ".csv", sep = "")},
    content = function(filename) {write.csv(Markers(), filename)}
  )
  

  output$MarkersTable <- DT::renderDataTable(
    ClusterMarkers() %>%  
      dplyr::select(Markers.column.order) %>% 
      dplyr::group_by(cluster) %>% 
      dplyr::top_n(n = n_top_markers, wt = avg_logFC) %>% 
      DT::datatable(options = list(scrollX = TRUE)) %>% 
      DT::formatRound(c('avg_logFC')) %>% 
      DT::formatSignif(c('p_val','p_val_adj'))
  )
  #>>>>>>>>>>>>>>>>>>> Notes Part <<<<<<<<<<<<<<<<<<<<<<<<<<<<<#
  notes.file <- "notes.rds"
  notesExist <- reactiveVal(file.exists(notes.file))
  
  Check.Notes.Para <- eventReactive(input$NotesSubmit, {
      message("Action-4: Checking Notes Parameters...")
    # check parameters
    if (trimws(input$NotesName) == "" | trimws(input$NotesContent) == "" | trimws(input$NotesContent) == "Notes, such as Data quality, Analysis Evaluation, Biology Interpretation.") {
      showModal(modalDialog(title = "Error", "Please Input Your Name and Note Content!", size = "l"))
      return(FALSE)
    } else {
      return(TRUE)
    }
  })
  
  format_content <- function(df){
    paste(sapply(nrow(df):1,function(x){paste(paste(colnames(df), df[x,],sep = ": "),collapse = "\n")}),collapse = "\n\n")
  }
  
  observeEvent(input$NotesSubmit,{
    if(Check.Notes.Para()) {
      if(notesExist()) {
        notes.all <- readRDS(notes.file)
      }else{
       notes.all <- data.frame(Date = character(0),
                              Data = character(0),
                              User = character(0),
                              Content = character(0),
                              stringsAsFactors = FALSE)
      }
      notes.new <- data.frame(Date = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
                              Data = input$NotesDataset,
                              User = input$NotesName,
                              Content = input$NotesContent,
                              stringsAsFactors = FALSE)
      notes.final <- rbind(notes.all, notes.new)
      saveRDS(notes.final, file = notes.file)
      showModal(modalDialog(title = "Cons", "New notes added!", size = "l"))
      output$OldNotes <- renderText({  format_content(notes.final) })
  }})
  
  output$OldNotes <- renderText({  
    if(notesExist()) {
      notes.all <- readRDS(notes.file)
      format_content(notes.all)
    }else{
      "No Notes Avaliable!"
    }
  })
  
  
}