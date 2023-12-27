library(Seurat)
library(stringr)
library(shiny)
library(ggplot2)
library(jpeg)
library(dplyr)
library(shinydashboard)
library(DT)
#********************************************************************************************************************#
# @Varibles
# importdata: ClusterMarkers,FileName
# ui: FeaturePlot, DimensionReduction, GeneSymbol, UpdatePlot, DownloadPlot, Group1, Group2, ClusterInfo,CalculateDEG
# server: PlotFeature, Reaction, DownloadDEG, Markers, MarkersTable

# @Attention
# 1. no ".","-" in variable names! and alse the dashboard name!
# 2. RunApp with external website will give a correct filename!
#********************************************************************************************************************#

# data settings
EC <- readRDS("EC-dims15-res0.5.rds")
EC@project.name <- "EC"
ECMarkers <- read.csv("markers-EC-dims15.csv")

APC <- readRDS("scRNA_APC.rds")
APC@project.name <- "APC"
APCMarkers <- read.csv("markers_APC_cluster.csv")

Ttk5d <- readRDS("Ttk-5d-dims15-res0.5.rds")
Ttk5d@project.name <- "Ttk-5d"
Ttk5dMarkers <- read.csv("markers-Ttk-5d-cluster.csv")

Ttk14drep1 <- readRDS("Ttk-14d-rep1-dims11-res0.5.rds")
Ttk14drep1@project.name <- "Ttk-14d-rep1"
Ttk14drep1dMarkers <- read.csv("markers-Ttk-14d-rep1-cluster.csv")

Ttk14drep2 <- readRDS("Ttk-14d-rep2-dims15-res0.5.rds")
Ttk14drep2@project.name <- "Ttk-14d-rep2"
Ttk14drep2Markers <- read.csv("markers-Ttk-14d-rep2-cluster.csv")

## app.R ##
ui <- dashboardPage(
  dashboardHeader(title = "Lab of Xi Rongwen"),
  dashboardSidebar(
    sidebarMenu(
      menuItem("About This Site", tabName = "widgets", icon = icon("th")),
      menuItem("Fly gut scRNAseq", tabName = "dashboard", icon = icon("dashboard"))
    )
  ),
  dashboardBody(
    tabItems(
      # First tab content
      tabItem(tabName = "dashboard",
              fluidRow(
                box(title = "Scatter Plot",plotOutput("FeaturePlot",height = "550px"),status = "primary",collapsible = TRUE,solidHeader = TRUE),
                box(title = "Settings",solidHeader = TRUE,
                    selectInput("dataset", "Choose a dataset:", 
                                choices = c(
                                  "EC", 
                                  "APC","Ttk-5d",
                                  "Ttk-14d-rep1",
                                  "Ttk-14d-rep2"
                                  )),
                    helpText(""),
                    selectInput("DimensionReduction","Dimension Reduction:",choices = c("UMAP" = "umap","tSNE" = "tsne")),
                    textInput("GeneSymbol","Gene Symbol:",value = ""),
                    helpText("Note: please input a gene symbol(not case-sensitive,except for genes like Dl,dl,H,hald,Ald,and so on), 
                             if no gene symbols, using CG number instead. Also support for nFeature_RNA, percent.mt, nCount_RNA."),
                    helpText(""),
                    helpText(""),
                    actionButton("UpdatePlot", "Update View"),
                    helpText(""),
                    downloadButton('DownloadPlot', label = 'Download Plot'),
                    helpText(""),
                    textInput("Group1","Group 1:",value = "1"),
                    textInput("Group2","Group 2:",value = "2"),
                    actionButton("CalculateDEG", "Calculate DEGs!"),
                    helpText(""),
                    status = "primary",collapsible = TRUE)
                ),
              fluidRow(
                box(height = "850px",solidHeader = TRUE,
                    title = "Cluster Info",
                    plotOutput("ClusterInfo"),collapsible = TRUE,status = "primary"
                ),
                box(solidHeader = TRUE,
                  title = "Top 20 Markers of each Cluster By log(FC)",
                  DT::dataTableOutput("MarkersTable"),collapsible = TRUE,status = "primary"
                )
              )   
      ) ,
      # About this site
      tabItem(tabName = "widgets",
              h1("Welcome to Xi Lab!",align = "center"),
              h2(""),
              h5("A Database Only for Internal Use.",align = "center"),
              h2("")
      )   
    )
  )
)



server <- function(input, output,session) {
  # Return the requested dataset
  datasetInput <- reactive({
    switch(input$dataset,
           "EC" = EC,
           "APC" = APC,
           "Ttk-5d" = Ttk5d,
           "Ttk-14d-rep1" = Ttk14drep1,
           "Ttk-14d-rep2" = Ttk14drep2
           )
  })
  # return the markers of each cluster
  ClusterMarkers <- reactive({
    if ( datasetInput()@project.name == "EC" ) { return(ECMarkers) }
    if ( datasetInput()@project.name == "APC" ) { return(APCMarkers) }
    if ( datasetInput()@project.name == "Ttk-5d" ) { return(Ttk5dMarkers) }
    if ( datasetInput()@project.name == "Ttk-14d-rep1" ) { return(Ttk14drep1dMarkers) }
    if ( datasetInput()@project.name == "Ttk-14d-rep2" ) { return(Ttk14drep2Markers) }
  })
  
  # return cluster info plot
  output$ClusterInfo <- renderPlot({
    if ( datasetInput()@project.name == "EC" ) { clusterinfo <- readJPEG("EC-anno-small.jpg") }
    if ( datasetInput()@project.name == "APC" ) { clusterinfo <- readJPEG("APC-anno-small.jpg") }
    if ( datasetInput()@project.name == "Ttk-5d" ) { clusterinfo <- readJPEG("Ttk-5d-anno.jpg") }
    if ( datasetInput()@project.name == "Ttk-14d-rep1" ) { clusterinfo <- readJPEG("Ttk-14d-rep1-anno.jpg") }
    if ( datasetInput()@project.name == "Ttk-14d-rep2" ) { clusterinfo <- readJPEG("Ttk-14d-rep2-anno.jpg") }
    
    if(exists("rasterImage")){
      plot(1:2, type='n',yaxt="n",xaxt="n",xlab="",ylab="")
      rasterImage(clusterinfo,1,1,2,2)}
  },
  height = function() {
    session$clientData$output_ClusterInfo_width}
  )
  
  # Scatter Plot
  PlotFeature <- reactive({
    if (input$GeneSymbol == "") {
      print(paste("scRNAseq-",datasetInput()@project.name,"-DimPlot-",input$DimensionReduction,sep=""))
      DimPlot(datasetInput(),reduction = input$DimensionReduction,label = TRUE,coord.fixed = TRUE)
    }else {
      gene <- tolower(input$GeneSymbol)
      gene_library <- tolower(rownames(datasetInput()@assays$RNA@data))
      if (input$GeneSymbol %in% c(rownames(datasetInput()@assays$RNA@data),"nFeature_RNA", "percent.mt", "nCount_RNA")) {
        print(paste("scRNAseq-",datasetInput()@project.name,"-FeaturePlot-",input$GeneSymbol,sep=""))
        FeaturePlot(datasetInput(),features = input$GeneSymbol,cols = c("gray", "red"),
                    label = TRUE,label.size = 7,coord.fixed = TRUE,reduction = input$DimensionReduction)
      }else if(gene %in% gene_library) {
        gene_name <- rownames(datasetInput()@assays$RNA@data)[gene_library == gene]
        print(paste("scRNAseq-",datasetInput()@project.name,"-FeaturePlot-",gene_name,sep=""))
        FeaturePlot(datasetInput(),features = gene_name,cols = c("gray", "red"),
                    label = TRUE,label.size = 7,coord.fixed = TRUE,reduction = input$DimensionReduction)
      }else{
        error <- readJPEG("www/error.jpg")
        if(exists("rasterImage")){
          plot(1:2, type='n',yaxt="n",xaxt="n",xlab="",ylab="")
          rasterImage(error,1,1,2,2)
        }
      }
    }
  })
  
  Reaction <- eventReactive(input$UpdatePlot,PlotFeature())
  
  output$FeaturePlot <- renderPlot({
    if(input$UpdatePlot == 0){return(PlotFeature())}else{return(Reaction())}
  })
  
  output$DownloadPlot <- downloadHandler(
    print( paste("scRNAseq-",datasetInput()@project.name, "_",input$GeneSymbol, "_",input$DimensionReduction,'.png', sep='') ),
    filename = function() { paste("scRNAseq-",datasetInput()@project.name, "_",input$GeneSymbol, "_",input$DimensionReduction,'.png', sep='') },
    content = function(filename) {ggsave(filename,plot = PlotFeature(), 
                                         # plot size are same with what client see in website
                                         width = 10, height= 10,units="in",device = "png")
    })
  
  Markers <- eventReactive(input$CalculateDEG, {
    if (input$Group1 == "" & input$Group2 == "") {
      print(paste("Download markers of each cluser from data - ",datasetInput()@project.name,sep = ""))
      return(ClusterMarkers())}
    else{
      g1 <- unlist(strsplit(as.character(input$Group1),","))
      g2 <- unlist(strsplit(as.character(input$Group2),","))
      print(paste("Calculate DEGs between cluster",g1,"-and-",g2," using data-",datasetInput()@project.name,sep=""))
      FindMarkers(datasetInput(), ident.1 = as.numeric(g1), ident.2 = as.numeric(g2), min.pct = 0.25)
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
    filename = function(){paste("scRNAseq-",datasetInput()@project.name, "_cluster",gsub(",","-",input$Group1),"_Vs_",gsub(",","-",input$Group2),".csv", sep = "")},
    content = function(filename) {write.csv(Markers(), filename)}
  )
  
  output$MarkersTable <- DT::renderDataTable(
    ClusterMarkers() %>%  group_by(cluster) %>% top_n(n = 20, wt = avg_logFC) %>% DT::datatable(options = list(scrollX = TRUE)) %>% 
      DT::formatRound(c('avg_logFC')) %>% DT::formatSignif(c('p_val','p_val_adj'))
  )
}

shinyApp(ui, server)