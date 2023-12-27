# library
library(Seurat)
library(stringr)
library(shiny)
library(ggplot2)
library(jpeg)
library(shinydashboard)

# import rds file of the scRNAseq data
pbmc <- readRDS("Fly-Gut-EEs-scRNAseq.rds")

# define UI
ui <- dashboardPage(
  dashboardHeader(title = "Flygut-EEs"),
  dashboardSidebar(
    sidebarMenu(
      menuItem("About the Database", tabName = "widgets", icon = icon("th")),
      menuItem("Flygut-EEs", tabName = "dashboard", icon = icon("dashboard")),
      #menuItem("National Institute of Biological Sciences, Beijing", icon = icon("home"),href = "http://www.nibs.ac.cn/en/index.php"),
      menuItem("Contact Us", tabName = "contact",icon = icon("address-book"))
    )
  ),
  dashboardBody(
    tabItems(
      # First tab content - Fly gut EEs scRNAseq
      tabItem(tabName = "dashboard",
              fluidRow(
                box(title = "Feature Plot",solidHeader = TRUE,
                    plotOutput("FeaturePlot",height = "680px"),
                    status = "primary",
                    collapsible = TRUE),
                box(title = "Settings",solidHeader = TRUE,
                    selectInput("resolution","Choose Resolution(Seurat parameter):",choices = c("Res0.4 Default in paper" = "0.4","Res0.6" = "0.6","Res1.0" = "1.0")),
                    selectInput("dimension_reduction","Dimension Reduction Method:",choices = c("tSNE" = "tsne","UMAP" = "umap")),
                    textInput("GeneSymbol","Gene Symbol:",value = ""),
                    helpText("Note: Please input a gene symbol ( not case-sensitive, except for genes below ), if no gene symbols exist, use CG number instead. If the input gene is not found in our data, this will return the default cluster map."),
                    helpText(""),
                    helpText("Exception List: ald, Ald, caps, Caps, cbs, Cbs, crc, Crc, Dip1, DIP1, Dip2, DIP2, dl, Dl, dlp, DLP, dor, DOR, dys, Dys, eap, Eap, faf, Faf, 
                              grp, Grp, h, H, hk, Hk, hop, Hop, Lkr, LKR, mio, Mio, Nat1, NAT1, pcl, Pcl, pen, Pen, pfk, Pfk, r, R, RpL37a, RpL37A, sip1, Sip1, smg, SmG, 
                              sug, Sug, trn, Trn."),
                    #actionButton("updateplot", "Update View"),
                    helpText(""),
                    downloadButton('downloadPlot', label = 'Download Plot'),
                    helpText(""),
                    downloadButton('downloadMarker', label = 'Download Markers of Each Cluster'),
                    helpText(""),
                    helpText("Note: This will download the markers of each cluster defined by the choosen Resolution."),
                    #textInput("group1","Group 1:",value = "1"),
                    #textInput("group2","Group 2:",value = "2"),
                    #helpText(""),
                    #helpText("To calculate differentially expressed genes (DEGs) between two groups, Please input the cluter number separted by comma, such as, Group1: 1,3,4 
                    #         and Group2: 2,5,6. leave Both Group blank to calculate markers of each cluter. Download DEGs may take a few miniutes."),
                    #helpText(""),
                    #actionButton("CalculateDEG", "Calculate DEGs"),
                    status = "info",
                    collapsible = TRUE)
             ),
             
              fluidRow(
                box(title = "Location of each cell cluster (Resolution 0.4)",width = 8,solidHeader = TRUE,
                    plotOutput("ClusterLocation",width = "960px",height = "600px"),
                    collapsible = TRUE,
                    status = "success")
             )
          ),
      
      # Second tab content - About the Database
      tabItem(tabName = "widgets",
              h2(strong("Welcome to Flygut-EEs"),align = "center"),
              p("This website is built based on R package - shiny and shinydashboard",align = "center"),
              h2(""),
              h3(strong("Flygut EEs scRNA-seq")),
              h4("Single cell transcriptome analysis of Drosophila enteroendocrine cells. EEs were labeled and sorted out by FACS, and single cell analysis was 
                 done by 10xGenomics and Seurat."),
              h2(""),
              h3("Reference:"),
              h4("Guo, Xingting and Yin, Chang and Yang, Fu and Zhang, Yongchao and Huang, Huanwei and Wang, Jiawen and Deng, Bowen and Cai, Tao and Rao, Yi and Xi, Rongwen, 
                 The Cellular Diversity and Transcription Factor Code of Drosophila Enteroendocrine Cells (June 25, 2019). CELL-REPORTS-D-19-02252. Available at SSRN: ",
                 tags$a(href = "https://ssrn.com/abstract=3409458","https://ssrn.com/abstract=3409458") ,"or",
                 tags$a(href = "http://dx.doi.org/10.2139/ssrn.3409458","http://dx.doi.org/10.2139/ssrn.3409458")),
              h3(""),
              h3("")
              ),
      # Third tab content - Contact US
      tabItem(tabName = "contact",
              h2(""),
              h2("Scientific Contacts"),
              h2(""),
              h4("Rongwen Xi"),
              h4(tags$a(href="http://www.nibs.ac.cn/en/yjsjyimgshow.php?cid=5&sid=6&id=774","National Institute of Biological Sciences, Beijing")),
              h4(tags$a(href = "mailto:xirongwen@nibs.ac.cn","xirongwen@nibs.ac.cn")),
              h2(""),
              h4("Xingting Guo"),
              h4(tags$a(href = "mailto:guoxingting@nibs.ac.cn","guoxingting@nibs.ac.cn")),
              h2(""),
              h2(""),
              h2(""),
              h2("Website Contact"),
              h2(""),
              h4("Yongchao Zhang"),
              h4(tags$a(href="mailto:zhangyongchao@nibs.ac.cn","zhangyongchao@nibs.ac.cn")),
              h2(""),
              h3("")
      )
            )
          )
        )

# define server
server <- function(input, output) {
 featureplot <- reactive({
   # change idents
    if (input$resolution == "0.4") {Idents(pbmc) <- pbmc@meta.data$res.0.4}
    if (input$resolution == "0.6") {Idents(pbmc) <- pbmc@meta.data$RNA_snn_res.0.6}
    if (input$resolution == "1.0") {Idents(pbmc) <- pbmc@meta.data$RNA_snn_res.1}
  # draw the plot
    if (input$GeneSymbol == "") {
        DimPlot(pbmc,reduction = input$dimension_reduction,label = TRUE,coord.fixed = TRUE)}
    else {
        gene <- tolower(input$GeneSymbol)
        gene_library <- tolower(rownames(pbmc@assays$RNA@data))
        if (input$GeneSymbol %in% c(rownames(pbmc@assays$RNA@data),"nFeature_RNA", "percent.mt", "nCount_RNA")) {
            FeaturePlot(pbmc,features = input$GeneSymbol,cols = c("gray", "red"),
                        label = TRUE,label.size = 7,coord.fixed = TRUE,reduction = input$dimension_reduction)}
        else if(gene %in% gene_library) {
            gene_name <- rownames(pbmc@assays$RNA@data)[gene_library == gene]
            FeaturePlot(pbmc,features = gene_name,cols = c("gray", "red"),
                        label = TRUE,label.size = 7,coord.fixed = TRUE,reduction = input$dimension_reduction)}
        else {DimPlot(pbmc,reduction = input$dimension_reduction,label = TRUE,coord.fixed = TRUE)}
      }
    })

  output$FeaturePlot <- renderPlot(featureplot())
  
  output$downloadPlot <- downloadHandler(
    filename = function() { paste("Fly-Gut-EEs-scRNAseq_Gene-",input$GeneSymbol, '.png', sep='') },
    content = function(filename) {ggsave(filename,plot = featureplot(), device = "png")}
    )
      

  markers <- reactive({
      if (input$resolution == "0.4"){return(read.csv(file = "Fly-Gut-EEs-scRNAseq-Resolution0.4-Cluster-Markers.csv"))}
      if (input$resolution == "0.6"){return(read.csv(file = "Fly-Gut-EEs-scRNAseq-Resolution0.6-Cluster-Markers.csv"))}
      if (input$resolution == "1.0"){return(read.csv(file = "Fly-Gut-EEs-scRNAseq-Resolution1.0-Cluster-Markers.csv"))}
      })
  
  output$downloadMarker <- downloadHandler(
    filename = function() {paste("Fly-Gut-EEs-scRNAseq_Resolution",input$resolution, "-Cluster-Markers.csv", sep = "")},
    content = function(filename) {write.csv(markers(), filename)}
  )
  
  output$ClusterLocation <- renderPlot({
    clusterlocation <- readJPEG("www/figure_4C.jpg")
    if(exists("rasterImage")){
      plot(1:2, type='n',yaxt="n",xaxt="n",xlab="",ylab="")
      rasterImage(clusterlocation,1,1,2,2)}
  })
}

shinyApp(ui, server)
