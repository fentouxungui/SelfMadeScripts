# define some credentials
lastElement <- function(Avector){
  Avector[length(Avector)]
}
AppID <- unlist(lapply(strsplit(getwd(),"/"),FUN = lastElement))

AppOwerCredentials <- readRDS("/srv/shiny-server/scRNAseq/initiated.account.and.password.rds")[[AppID]][,c(1:6)]
AdminCredentials <- readRDS("/srv/shiny-server/scRNAseq/initiated.account.and.password.rds")[["Admin"]][,c(1:6)]
credentials <- rbind(AppOwerCredentials,AdminCredentials)

# load library
library(shinymanager)
library(Seurat)
library(stringr)
library(shiny)
library(ggplot2)
library(jpeg)
library(RColorBrewer)
library(dplyr)
library(shinydashboard)
library(DT)
library(ggpubr)
library(shinycssloaders)
#********************************************************************************************************************#
# @Varibles
# importdata: ClusterMarkers,FileName
# ui: FeaturePlot, DimensionReduction, GeneSymbol, UpdatePlot, DownloadPlot, Group1, Group2, ClusterInfo,CalculateDEG
# server: PlotFeature, Reaction, DownloadDEG, Markers, MarkersTable

# @Attention
# 1. no ".","-" in variable names! and alse the dashboard name!
# 2. RunApp with external website will give a correct filename!
#********************************************************************************************************************#


batch.analysis.dir <- "/home/xilab/qiannn/scRNA_Mouse/DownstreamAnalysis/Single-Sample/cellbender-scrublet-seurat/"

# rds.files <- list.files(batch.analysis.dir,pattern = "(removed.rds$)|(scrublet.rds$)")
rds.files <- list.files(batch.analysis.dir,pattern = "removed.rds$")
markers.files <- gsub("rds$","cluster.markers.csv",rds.files)
project.names <- gsub("\\.rds$","",rds.files)

# data settings
data.list <- data.frame("Name" = project.names,
                        "Rds" = paste(batch.analysis.dir,rds.files,sep = ""),
                        "Project.Name" = project.names ,
                        "Markers.File" = paste(batch.analysis.dir,markers.files,sep = ""),
                        "Cluster.Anno" = rep("VlnPlot",length(rds.files)),
                        stringsAsFactors = FALSE) 
data.list


cache.list <- list()

selectInput.Dataset <- data.list$Project.Name
color.vector <- colorRampPalette(c("#638ED0", "#ffff33","#ff3300"))(15)

# some function

FeaturePlotSingle<- function(obj, feature, metadata_column, ...){
  all_cells<- colnames(obj)
  groups<- levels(obj@meta.data[, metadata_column])
  print(groups)
  # the minimal and maximal of the value to make the legend scale the same. 
  minimal<- min(obj[['RNA']]@data[feature, ])
  maximal<- max(obj[['RNA']]@data[feature, ])
  ps<- list()
  for (group in groups) {
    subset_indx<- obj@meta.data[, metadata_column] == group
    subset_cells<- all_cells[subset_indx]
    print(length(subset_cells))
    if (length(subset_cells) < 2) { # if cells less than 2 will cause an error
      next()
    }else {
      if (length(unique(obj@assays$RNA@data[feature,subset_cells])) == 1) {
        minimal.color <- color.vector[14*unique(obj@assays$RNA@data[feature,subset_cells])/maximal+1] # color decided by expr value
        p<- FeaturePlot(obj, features = feature, cols = c(minimal.color, color.vector[length(color.vector)]), cells= subset_cells, ...) +
          ggtitle(group) +
          theme(plot.title = element_text(size = 10, face = "bold"))
      }else{
        p<- FeaturePlot(obj, features = feature, cells= subset_cells, ...) +
          scale_color_viridis_c(limits=c(minimal, maximal), direction = 1) +
          scale_colour_gradientn(colours = color.vector) +
          ggtitle(group) +
          theme(plot.title = element_text(size = 10, face = "bold"))
      }
      ps[[group]]<- p
    }
  }
  return(ps)
}

DimPlotSingle<- function(obj, metadata_column, ...){
  all_cells<- colnames(obj)
  groups<- levels(obj@meta.data[, metadata_column])
  ps<- list()
  for (group in groups) {
    print(group)
    subset_indx<- obj@meta.data[, metadata_column] == group
    print(length(subset_indx))
    subset_cells<- all_cells[subset_indx]
    if (length(subset_cells) < 2) { # if cells less than 2 will cause an error
      next()
    }else {
      p <- DimPlot(obj, cells= subset_cells, ...) +
        ggtitle(group) +
        theme(plot.title = element_text(size = 10, face = "bold"))
        ps[[group]]<- p
      }

  }
  return(ps)
}




## app.R ##
ui <- secure_app(
    dashboardPage(
    dashboardHeader(title = "Lab of Xi Rongwen"),
    dashboardSidebar(
      sidebarMenu(
        menuItem("About This Site", tabName = "widgets", icon = icon("th")),
        menuItem("Mouse scRNAseq", tabName = "dashboard", icon = icon("dashboard"))
      )
    ),
    dashboardBody(
      tabItems(
        # First tab content
        tabItem(tabName = "dashboard",
                fluidRow(
                  box(title = "Scatter Plot",
                      shinycssloaders::withSpinner(plotOutput("FeaturePlot",height = "550px")),
                      width = 9,status = "primary",collapsible = TRUE,solidHeader = TRUE),
                  box(title = "Settings",solidHeader = TRUE,width = 3,
                      selectInput("Dataset", "Choose a dataset:", choices = selectInput.Dataset),
                      shinycssloaders::withSpinner(uiOutput("Reductions.UI"),proxy.height = "10px"),
                      # selectInput("DimensionReduction","Dimension Reduction:",choices = c("UMAP" = "umap","tSNE" = "tsne")),
                      # selectInput("ClusterResolution","Cluster Resolution:",choices = ClusterResolution.Choice ),
                      shinycssloaders::withSpinner(uiOutput("ClusterResolution.UI"),proxy.height = "10px"),
                      # selectInput("Split","Split by:", choices = c("No" = "No",SplitBy.Choice)),
                      shinycssloaders::withSpinner(uiOutput("Split.UI"),proxy.height = "10px"),
                      textInput("GeneSymbol","Gene Symbol:",value = ""),
                      helpText("Note: also support for nFeature_RNA, percent.mt, nCount_RNA, Phase."),
                      checkboxInput("ShowLable",label = "Show cluster label",TRUE),
                      sliderInput("LabelSize",label = "Label Size:",min = 0,max = 10,value = 7),
                      sliderInput("PointSize",label = "Point Size",min = 0.001,max = 2,value = 0.05),
                      downloadButton('DownloadPlot', label = 'Download Plot'),
                      textInput("Group1","Group 1:",value = "1"),
                      textInput("Group2","Group 2:",value = "2"),
                      actionButton("CalculateDEG", "Calculate DEGs!"),
                      textInput("Cluster","Cluster:",value = "0"),
                      actionButton("CalculateDEGWithin", "Calculate DEGs Within A Cluster"),
                      status = "primary",collapsible = TRUE)
                ),
                fluidRow(
                  box(height = "850px",solidHeader = TRUE,
                      title = "Cluster Info",
                      shinycssloaders::withSpinner(plotOutput("ClusterInfo")),
                      collapsible = TRUE,status = "primary"
                  ),
                  box(solidHeader = TRUE,
                      title = "Top 10 Markers of each Cluster By log(FC) Using Resolution 0.6",
                      shinycssloaders::withSpinner(DT::dataTableOutput("MarkersTable")),
                      collapsible = TRUE,status = "primary"
                  )
                )   
        ) ,
        # About this site
        tabItem(tabName = "widgets",
                h1("Welcome to Xi Lab!",align = "center"),
                h2(""),
                h5("A Database Only for Internal Use.",align = "center"),
                h2(""),
                h3(strong("1. CellRanger and Cellbender Report")),
                p(tags$a(href="http://192.168.13.45//shiny-server/PrivateData/xilab-jin/Mouse/Analysis-Pipeline-Mouse-SingleSample-Large-Intestine-Jinzhen/1.CellRanger-And-Cellbender-Report/","Here")),
                br(),
                h3(strong("1. Downstream Analysis Report")),
                p(tags$a(href="http://192.168.13.45//shiny-server/PrivateData/xilab-jin/Mouse/Analysis-Pipeline-Mouse-SingleSample-Large-Intestine-Jinzhen/2.Downstream-analysis-report/","Here"))
          )
        )
    )
  )
)


server <- function(input, output,session) {
  
  # call the server part
  # check_credentials returns a function to authenticate users
  res_auth <- secure_server(
    check_credentials = check_credentials(credentials)
  )
  
  output$auth_output <- renderPrint({
    reactiveValuesToList(res_auth)
  })
  
  dataset.input <- reactive({
    if (is.null(names(cache.list)) | ! input$Dataset %in% names(cache.list)) {
      Rds.raw <- readRDS(data.list$Rds[match(input$Dataset,selectInput.Dataset)])
      Rds.raw@project.name <- input$Dataset
      cache.list[[input$Dataset]] <<- Rds.raw
    }else{
      Rds.raw <- cache.list[[input$Dataset]]
    }
    return(Rds.raw)}
  )
    
  output$Reductions.UI <- renderUI({
    reduction.choice <- rev(Reductions(dataset.input()))
    names(reduction.choice) <- rev(toupper(Reductions(dataset.input())))
    selectInput("DimensionReduction","Dimension Reduction:",choices = reduction.choice)
  })
  
  
  output$ClusterResolution.UI <- renderUI({
    full.res <- rev(grep("RNA_snn",colnames(dataset.input()@meta.data),value = TRUE))
    part.res <- gsub("RNA_snn_","",full.res)
    names(full.res) <- part.res
    ClusterResolution.Choice <- full.res
    selectInput("ClusterResolution","Cluster Resolution:",choices = ClusterResolution.Choice )
  })
  
  
  output$Split.UI <- renderUI({
    SplitBy.Choice.total <- c("orig.ident", "type", "subject","Phase")
    names(SplitBy.Choice.total) <- c("By Sample", "By Type","By Subject","By Cell Cycle")
    SplitBy.Choice <- names(dataset.input()@meta.data)[toupper(names(dataset.input()@meta.data)) %in% toupper(SplitBy.Choice.total)]
    selectInput("Split","Split by:", choices = c("No" = "No",SplitBy.Choice))
  })
  
  

  # return the markers of each cluster
  ClusterMarkers <- reactive({
    read.csv(data.list$Markers.File[match(input$Dataset,selectInput.Dataset)])
  })
  
  # return cluster info plot
  output$ClusterInfo <- renderPlot({
    dataset.input.isolate <- isolate(dataset.input())
    if (!is.null(input$ClusterResolution)) {
      Idents(dataset.input.isolate) <- input$ClusterResolution
    }
    
    fig.dir <- data.list$Cluster.Anno[match(input$Dataset,selectInput.Dataset)]
    if (is.na(fig.dir)) {
      clusterinfo <- readJPEG("www/figure_4C.jpg")
      if(exists("rasterImage")){
        plot(1:2, type='n',yaxt="n",xaxt="n",xlab="",ylab="")
        rasterImage(clusterinfo,1,1,2,2)}
    } else if(fig.dir == "VlnPlot"){
      if (input$GeneSymbol != "") {
        gene <- tolower(input$GeneSymbol)
        gene_library <- tolower(rownames(dataset.input.isolate@assays$RNA@data))
        if (input$GeneSymbol %in% c(rownames(dataset.input.isolate@assays$RNA@data),"nFeature_RNA", "percent.mt", "nCount_RNA","Phase")) {
          gene.input <- input$GeneSymbol
        }else if( gene %in% gene_library){
          gene.input <-  rownames(dataset.input.isolate@assays$RNA@data)[gene_library == gene][1]
        }else{
          gene.input <- NULL
        }
        
        if (!is.null(gene.input)) {
          # whether to split feature plots
          if( is.null(input$Split)){  # (is.null(input$Split) | input$Split == "No") will cause an error!
            split <- NULL
          }else if( input$Split == "No" ){
            split <- NULL
          }else{
            split <- input$Split
            if (!is.factor(dataset.input.isolate@meta.data[,split])) {
              dataset.input.isolate@meta.data[,split] <- factor(dataset.input.isolate@meta.data[,split])
            }
          }
          
          
          if (split  == "type") {
            VlnPlot(dataset.input.isolate, features = gene.input,split.by = "type", pt.size = 0)
          }else{
            VlnPlot(dataset.input.isolate, features = gene.input,pt.size = 0)
          }
        }else{
          clusterinfo <- readJPEG("www/figure_4C.jpg")
          if(exists("rasterImage")){
            plot(1:2, type='n',yaxt="n",xaxt="n",xlab="",ylab="")
            rasterImage(clusterinfo,1,1,2,2)}
        }
      }else{
        clusterinfo <- readJPEG("www/figure_4C.jpg")
        if(exists("rasterImage")){
          plot(1:2, type='n',yaxt="n",xaxt="n",xlab="",ylab="")
          rasterImage(clusterinfo,1,1,2,2)}
      }
      
    } else {
      clusterinfo <- readJPEG(fig.dir)
      if(exists("rasterImage")){
        plot(1:2, type='n',yaxt="n",xaxt="n",xlab="",ylab="")
        rasterImage(clusterinfo,1,1,2,2)}
    }


  },
  height = function() {
    session$clientData$output_ClusterInfo_width}
  )
  
  
  # Scatter Plot
  PlotFeature <- reactive({
    dataset.input.isolate <- isolate(dataset.input())
    print(input$ClusterResolution)
    # set id
    if (!is.null(input$ClusterResolution)) {
      Idents(dataset.input.isolate) <- input$ClusterResolution
    }
    
    
    # whether to split feature plots
    if( is.null(input$Split)){  # (is.null(input$Split) | input$Split == "No") will cause an error!
      split <- NULL
    }else if( input$Split == "No" ){
      split <- NULL
    }else{
      split <- input$Split
      if (!is.factor(dataset.input.isolate@meta.data[,split])) {
        dataset.input.isolate@meta.data[,split] <- factor(dataset.input.isolate@meta.data[,split])
      }
    }
  print(split)
  

  
  
  
    if (input$GeneSymbol == "") { # Dimplot
      print(paste("scRNAseq-",dataset.input.isolate@project.name,"-DimPlot-",input$DimensionReduction,sep=""))
      if (is.null(split)) {
        DimPlot(dataset.input.isolate,reduction = input$DimensionReduction,label = input$ShowLable,split.by = split, pt.size = input$PointSize, label.size = input$LabelSize)
      }else{
        p_list<- DimPlotSingle(dataset.input.isolate, metadata_column = split, reduction = input$DimensionReduction,label = input$ShowLable,pt.size = input$PointSize, label.size = input$LabelSize, order =TRUE)
        ggarrange(plotlist = p_list,common.legend = TRUE, legend="right")
      }
      
    }else {
      gene <- tolower(input$GeneSymbol)
      gene_library <- tolower(rownames(dataset.input.isolate@assays$RNA@data))
      if (input$GeneSymbol %in% c(rownames(dataset.input.isolate@assays$RNA@data),"nFeature_RNA", "percent.mt", "nCount_RNA","Phase")) {
        print(paste("scRNAseq-",dataset.input.isolate@project.name,"-FeaturePlot-",input$GeneSymbol,sep=""))
        if (is.null(split)) {
          FeaturePlot(dataset.input.isolate, features = input$GeneSymbol, split.by = NULL,pt.size = input$PointSize, label = input$ShowLable,label.size = input$LabelSize,reduction = input$DimensionReduction) +
            scale_colour_gradientn(colours = color.vector)
        } else{
          p_list<- FeaturePlotSingle(dataset.input.isolate, feature = input$GeneSymbol, metadata_column = split, pt.size = input$PointSize, label = input$ShowLable,label.size = input$LabelSize,reduction = input$DimensionReduction, order =TRUE)
          ggarrange(plotlist = p_list,common.legend = TRUE, legend="right")
        }
       
      }else if(gene %in% gene_library) {
        gene_name <- rownames(dataset.input.isolate@assays$RNA@data)[gene_library == gene][1]
        print(paste("scRNAseq-",dataset.input.isolate@project.name,"-FeaturePlot-",gene_name,sep=""))
        if (is.null(split)) {
          FeaturePlot(dataset.input.isolate, features = gene_name, split.by = NULL,pt.size = input$PointSize, label = input$ShowLable,label.size = input$LabelSize,reduction = input$DimensionReduction) +
            scale_colour_gradientn(colours = color.vector)
        }else{
          p_list<- FeaturePlotSingle(dataset.input.isolate, feature = gene_name, metadata_column = split, pt.size = input$PointSize, label = input$ShowLable,label.size = input$LabelSize,reduction = input$DimensionReduction, order =TRUE)
          ggarrange(plotlist = p_list,common.legend = TRUE, legend="right")
        }

      }else{
        error <- readJPEG("www/error.jpg")
        if(exists("rasterImage")){
          plot(1:2, type='n',yaxt="n",xaxt="n",xlab="",ylab="")
          rasterImage(error,1,1,2,2)
        }
      }
    }
  })
  
  output$FeaturePlot <- renderPlot({PlotFeature()})
  
  output$DownloadPlot <- downloadHandler(
    print( paste("scRNAseq-",dataset.input()@project.name, "_",input$GeneSymbol, "_",input$DimensionReduction,'.png', sep='') ),
    filename = function() { paste("scRNAseq-",dataset.input()@project.name, "_",input$GeneSymbol, "_",input$DimensionReduction,'.png', sep='') },
    content = function(filename) {ggsave(filename,plot = PlotFeature(), 
                                         # plot size are same with what client see in website
                                         width = 10, height= 10,units="in",device = "png")
    })
  
  Markers <- eventReactive(input$CalculateDEG, {
    if (input$Group1 == "" & input$Group2 == "") {
      print(paste("Download markers of each cluser from data - ",dataset.input()@project.name,sep = ""))
      return(ClusterMarkers())}
    else{
      g1 <- unlist(strsplit(as.character(input$Group1),","))
      g2 <- unlist(strsplit(as.character(input$Group2),","))
      print(paste("Calculate DEGs between cluster",g1,"-and-",g2," using data-",dataset.input()@project.name,sep=""))
      FindMarkers(dataset.input(), ident.1 = as.numeric(g1), ident.2 = as.numeric(g2), min.pct = 0.25)
    }
  })
  
  observeEvent(input$CalculateDEG,{
    if(!is.na(input$CalculateDEG)) {
      removeUI("#DownloadDEG")}
  })
  
  observeEvent(input$CalculateDEG,{
    insertUI(selector = "#CalculateDEG", where = "afterEnd", ui = downloadButton("DownloadDEG",label = "Download DEGs!"))
  })
  
  output$DownloadDEG <- downloadHandler(
    filename = function(){paste("scRNAseq-",dataset.input()@project.name, "_cluster",gsub(",","-",input$Group1),"_Vs_",gsub(",","-",input$Group2),".csv", sep = "")},
    content = function(filename) {write.csv(Markers(), filename)}
  )
  
  DEGs <- eventReactive(input$CalculateDEGWithin, {
    g3 <- unlist(strsplit(as.character(input$Cluster),","))
    FindMarkers(dataset.input(), ident.1 = levels(dataset.input()@meta.data$orig.ident)[1],
                ident.2 = levels(dataset.input()@meta.data$orig.ident)[2],group.by = "orig.ident",subset.ident = g3)
  })
  
  observeEvent(input$CalculateDEGWithin,{
    if(!is.na(input$CalculateDEGWithin)) {removeUI("#DownloadDEGWithin")}
  })
  
  observeEvent(input$CalculateDEGWithin,{
    insertUI("#CalculateDEGWithin","afterEnd",ui = downloadButton("DownloadDEGWithin",label = "Download DEGs Within A Cluster"))
  })
  
  output$DownloadDEGWithin <- downloadHandler(
    filename = function(){paste("scRNAseq-",dataset.input()@project.name,"-", levels(dataset.input()@meta.data$orig.ident)[1],"-Vs-",
                                levels(dataset.input()@meta.data$orig.ident)[2],"_cluster_",unlist(strsplit(as.character(input$Cluster),",")),".csv", sep = "")},
    content = function(filename) {write.csv(DEGs(), filename)}
  )
  
  output$MarkersTable <- DT::renderDataTable(
    ClusterMarkers() %>%  group_by(cluster) %>% top_n(n = 10, wt = avg_logFC) %>% DT::datatable(options = list(scrollX = TRUE)) %>% 
      DT::formatRound(c('avg_logFC')) %>% DT::formatSignif(c('p_val','p_val_adj'))
  )
}

shinyApp(ui, server)