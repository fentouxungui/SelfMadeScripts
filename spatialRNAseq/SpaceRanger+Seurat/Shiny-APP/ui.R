library(shiny)
library(shinymanager)
library(shinydashboard)
library(shinycssloaders)
library(DT)

# Initialize/Update the function: Attention: The new code will take effect next time.
## update codes while keep the old Parameters.R
update_code <- function(source.code.path = "/data0/reference/Scripts/Shiny-spatialRNAseq-Browse-App-Code/", 
                        file.backup = "./Parameters.R"){
  # backup parameters.R file
  file.copy(file.backup, paste0(file.backup,".backup"), overwrite = TRUE)
  system(paste("cp -r",  paste0(source.code.path,"*"), "./"))
  file.copy(paste0(file.backup,".backup"),file.backup, overwrite = TRUE)
}
# Commented out lines bellow before debuging!!!  ****** Attention ********
update_code()

# load parameters -  data.list and init the html reports
suppressMessages(source("./Parameters.R", local = FALSE) )
# Rearrange data list by rds file size
file.size(data.list$Rds)
data.list <- data.list[order(file.size(data.list$Rds)),]
  
# If you use the default value of local = FALSE, then the file will be sourced in the global environment.
# functions and objects will be used in ui.r and server.r
# do not put this code in app.R, for it will cause an error!
app.ui <-   dashboardPage(
    dashboardHeader(title = "Lab of Xi Rongwen"),
    dashboardSidebar(
      sidebarMenu(
        menuItem("Analysis Report", tabName = "widgets", icon = icon("th")),
        menuItem("Browse Data", tabName = "dashboard", icon = icon("dashboard")),
        menuItem("Notes", tabName = "notes", icon = icon("dashboard"))
      )
    ),
    dashboardBody(
      tabItems(
        # About this site - Analysis Report
        tabItem(tabName = "widgets",
                h1("Welcome to Xi Lab!",align = "center"),
                p(),
                h5("A Database Only for Internal Use.",align = "center"),
                br(),
                h3(strong("Analysis Reports")),
                p(tags$a(href=paste("http://192.168.13.45//",gsub("/data3/Xilab-Data-Analysis/xilab/","",dirname(getwd())),paste("Analysis-Pipeline-",basename(getwd()),"/",sep = ""),sep = "/"),"Here")),
                br(),
                h3(strong("Files info")),
                br(),
                fluidRow(
                  box(title = "Included Data", width = 12, status = "primary",
                      DT::dataTableOutput("UIDataList"))
                )
        ),
        # Browse Data
        tabItem(tabName = "dashboard",
                fluidRow(
                  box(title = "Scatter Plot",
                      shinycssloaders::withSpinner(plotOutput("FeaturePlot",height = "auto")), # Add a spinner that shows when an output is recalculating
                      width = 9, status = "primary", collapsible = TRUE, solidHeader = TRUE),
                  box(title = "Settings", solidHeader = TRUE, width = 3,
                      selectInput("Dataset", "Choose a dataset:", choices = data.list$Project.Name),
                      shinycssloaders::withSpinner(uiOutput("LinkedDimplot.UI"), proxy.height = "10px"), # proxy.height: spinner height
                      shinycssloaders::withSpinner(uiOutput("Reductions.UI"), proxy.height = "10px"), # proxy.height: spinner height
                      shinycssloaders::withSpinner(uiOutput("ClusterResolution.UI"), proxy.height = "10px"),
                      # shinycssloaders::withSpinner(uiOutput("Split.UI"), proxy.height = "10px"),
                      textInput("GeneSymbol", "Gene Symbol:", value = ""),
                      # helpText(strong(paste("Support: ", paste(Genes.qc, collapse = ", "), ".",sep = "")),style = "font-size:5px;"),
                      sliderInput("FeaturePlotHWRatio", label = "Scater Plot Height/Width Ratio", min = 0.1, max = 2, value = 0.5), # adjust the Ratio of width and height of scater plot.
                      checkboxInput("ShowLabel",label = "Show cluster label", TRUE),
                      sliderInput("LabelSize", label = "Label Size:", min = 0, max = 10, value = 7),
                      sliderInput("PointSize", label = "Left Point Size", min = 0.001, max = 2, value = 0.05),
                      sliderInput("SpatialPointSize", label = "Right Point Size", min = 0.5, max = 3, value = 1.6),
                      sliderInput("SpatialLowerTran", label = "Tranparency of points with lower expression", min = 0, max = 1, value = 1),
                      sliderInput("SpatialHighTran", label = "Tranparency of points with higher expression", min = 0, max = 1, value = 1),
                      downloadButton('DownloadPlot', label = 'Download Scater Plot - Not Recommmended'),
                      textInput("Group1","Group 1:", value = "1"),
                      textInput("Group2","Group 2:", value = "2"),
                      actionButton("CalculateDEG", "Calculate DEGs!"),
                      # textInput("Cluster", "Cluster:", value = "0"),
                      # actionButton("CalculateDEGWithin", "Calculate DEGs Within A Cluster"),
                      # shinycssloaders::withSpinner(uiOutput("Vlnplot.Cluster.Selection.UI"), proxy.height = "10px"),
                      sliderInput("ClusterAnnotationPlotHWRatio", label = "Cluster Annotation Plot Height/Width Ratio", min = 0.1, max = 2, value = 1), # adjust the Ratio of width and height of scater plot.
                      shinycssloaders::withSpinner(uiOutput("MarkerList.UI"), proxy.height = "10px"),
                      status = "primary", collapsible = TRUE)
                ),
                fluidRow(
                  box(solidHeader = TRUE,
                      title = "Cluster Annotation / Vlnplot",
                      shinycssloaders::withSpinner(plotOutput("ClusterAnnotationPlot", height = "auto")),
                      collapsible = TRUE, status = "primary"
                  ),
                  box(solidHeader = TRUE,
                      title = "Top Cluster Markers",
                      shinycssloaders::withSpinner(DT::dataTableOutput("MarkersTable")),
                      collapsible = TRUE, status = "primary"
                  )
                )   
        ),
        # Notes
        tabItem(tabName = "notes",
                h3(strong("Take New Notes")),
                br(),
                selectInput("NotesDataset", "Select dataset:", choices = data.list$Project.Name),
                textInput("NotesName",label = "Your Name"),
                textAreaInput(inputId = "NotesContent", "New Note", width = "100%",
                              placeholder = "Notes, such as Data quality, Analysis Evaluation, Biology Interpretation."),
                actionButton("NotesSubmit", "Submit"),
                br(),
                h3(strong("Previous Notes")),
                br(),
                verbatimTextOutput(outputId = "OldNotes", placeholder = TRUE)
        )
      )
    )
  )

if (Encrypted.App) {
  ui <- secure_app(tags_bottom = tags$div(
                       tags$p(
                         "For any question, please  contact ",
                         tags$a(
                           href = "mailto:zhangyongchao@nibs.ac.cn?Subject=Shiny%20aManager",
                           target="_top", "ZhangYongchao"
                         )
                       )
                     ),
                     background  = "linear-gradient(rgba(0, 0, 255, 0.5), rgba(255, 255, 0, 0.5))", 
                   app.ui)
}else{
  ui <- app.ui
}