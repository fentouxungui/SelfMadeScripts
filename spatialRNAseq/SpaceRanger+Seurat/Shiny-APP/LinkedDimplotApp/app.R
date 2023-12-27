library(Seurat)
library(ggplot2)
library(shiny)
library(miniUI)
library(magrittr)
library(DT)
library(patchwork)
object <- readRDS("")
# SpatialDimPlot(object, interactive = TRUE)
# SpatialFeaturePlot(object, features = "Epcam", interactive = TRUE)
# LinkedFeaturePlot(object = object, feature = "Epcam")

dims = 1:2
alpha = c(0.1, 1)
combine = TRUE

# image <- Seurat::DefaultImage(object = object)
image <- names(object@images)[1]
cells.use <- Seurat::Cells(x = object[[image]])
reduction <- "umap"
dims <- dims[1:2]
dims <- paste0(Key(object = object[[reduction]]), dims)
group.by <- "seurat_clusters"
group.data <- Seurat::FetchData(object = object, vars = group.by, 
                                cells = cells.use)
coords <- Seurat::GetTissueCoordinates(object = object[[image]])
embeddings <- Seurat::Embeddings(object = object[[reduction]])[cells.use, 
                                                               dims]
plot.data <- cbind(coords, group.data, embeddings)
plot.data$selected_ <- FALSE
Idents(object = object) <- group.by


ui <- miniPage(gadgetTitleBar(title = "LinkedDimPlot", 
                              left = miniTitleBarButton(inputId = "reset", label = "Reset")), 
               miniContentPanel(fillRow(plotOutput(outputId = "spatialplot",height = "100%", click = clickOpts(id = "spclick", clip = TRUE), 
                                                   hover = hoverOpts(id = "sphover", delay = 10, nullOutside = TRUE)), 
                                        plotOutput(outputId = "dimplot", height = "100%", brush = brushOpts(id = "brush",delay = 10, clip = TRUE, resetOnNew = FALSE), 
                                                   click = clickOpts(id = "dimclick",  clip = TRUE), hover = hoverOpts(id = "dimhover", delay = 10, nullOutside = TRUE)), 
                                        height = "97%"), 
                                verbatimTextOutput(outputId = "info"))
)

server <- function(input, output, session) {
  library(Seurat)
  click <- reactiveValues(pt = NULL, invert = FALSE)
  plot.env <- reactiveValues(data = plot.data, alpha.by = NULL)
  observeEvent(eventExpr = input$done, handlerExpr = {
    plots <- list(plot.env$spatialplot, plot.env$dimplot)
    if (combine) {
      plots <- wrap_plots(plots, ncol = 2)
    }
    stopApp(returnValue = plots)
  })
  observeEvent(eventExpr = input$reset, handlerExpr = {
    click$pt <- NULL
    click$invert <- FALSE
    session$resetBrush(brushId = "brush")
  })
  observeEvent(eventExpr = input$brush, handlerExpr = click$pt <- NULL)
  observeEvent(eventExpr = input$spclick, handlerExpr = {
    click$pt <- input$spclick
    click$invert <- TRUE
  })
  observeEvent(eventExpr = input$dimclick, handlerExpr = {
    click$pt <- input$dimclick
    click$invert <- FALSE
  })
  observeEvent(eventExpr = c(input$brush, input$spclick, 
                             input$dimclick), handlerExpr = {
                               plot.env$data <- if (is.null(x = input$brush)) {
                                 clicked <- nearPoints(df = plot.data, coordinfo = if (click$invert) {
                                   Seurat:::InvertCoordinate(x = click$pt)
                                 }
                                 else {
                                   click$pt
                                 }, threshold = 10, maxpoints = 1)
                                 if (nrow(x = clicked) == 1) {
                                   cell.clicked <- rownames(x = clicked)
                                   group.clicked <- plot.data[cell.clicked, group.by, 
                                                              drop = TRUE]
                                   idx.group <- which(x = plot.data[[group.by]] == 
                                                        group.clicked)
                                   plot.data[idx.group, "selected_"] <- TRUE
                                   plot.data
                                 }
                                 else {
                                   plot.data
                                 }
                               }
                               else if (input$brush$outputId == "dimplot") {
                                 brushedPoints(df = plot.data, brush = input$brush, 
                                               allRows = TRUE)
                               }
                               else if (input$brush$outputId == "spatialplot") {
                                 brushedPoints(df = plot.data, brush = Seurat:::InvertCoordinate(x = input$brush), 
                                               allRows = TRUE)
                               }
                               plot.env$alpha.by <- if (any(plot.env$data$selected_)) {
                                 "selected_"
                               }
                               else {
                                 NULL
                               }
                             })
  output$spatialplot <- renderPlot(expr = {
    plot.env$spatialplot <- Seurat::SingleSpatialPlot(data = plot.env$data, 
                                                      image = object[[image]], col.by = group.by, pt.size.factor = 1.6, 
                                                      crop = TRUE, alpha.by = plot.env$alpha.by) + 
      scale_alpha_ordinal(range = alpha) + NoLegend()
    plot.env$spatialplot
  })
  output$dimplot <- renderPlot(expr = {
    plot.env$dimplot <- Seurat::SingleDimPlot(data = plot.env$data, 
                                              dims = dims, col.by = group.by, alpha.by = plot.env$alpha.by, label = TRUE) + 
      scale_alpha_ordinal(range = alpha) + guides(alpha = FALSE)
    plot.env$dimplot
  })
  output$info <- renderPrint(expr = {
    cell.hover <- rownames(x = nearPoints(df = plot.data, 
                                          coordinfo = if (is.null(x = input[["sphover"]])) {
                                            input$dimhover
                                          }
                                          else {
                                            Seurat:::InvertCoordinate(x = input$sphover)
                                          }, threshold = 10, maxpoints = 1))
    if (length(x = cell.hover) == 1) {
      paste(cell.hover, paste("Group:", plot.data[cell.hover, 
                                                  group.by, drop = TRUE]), collapse = "<br />")
    }
    else {
      NULL
    }
  })
}
shinyApp(ui = ui, server = server)
