#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library (ggplot2)
library (ggraph)
library (igraph)
library (ggrepel)
library (data.table)


# ---- startup data ----

# load apms data
# baitPrey.dt <- fread ("./data/D_Final_v3.txt")
# 
# #standardize
# # keep it simple and gene based bait, prey mostly. 
# # possible weights will have weight prefixed to their name.
# # other types of fields should have other prefix, maybe attribute.cluster for example...
# baitPrey.dt <- baitPrey.dt[, .(bait = Bait, prey = PreyGene, prey.uniprot = Prey, 
#                                weight.AvgSpec = AvgSpec, weight.rank_WD = rank_WD, weight.inverse_BFDR = `1-BFDR`)]

baitPrey.dt <- fread ("./data/baitPrey.dt.csv")

allBaits <- sort(unique(baitPrey.dt$bait))
allPreys <- unique(baitPrey.dt$prey)


## ---- other edge types ----
otherEdgeDir <- "./data/otherEdges/"
edge.files <- list.files(path = otherEdgeDir ) # , pattern = "(tsv|csv|txt)$"
#print (edge.files)
if(length(edge.files) > 0){
  names(edge.files) <- tstrsplit(basename(edge.files), "\\.")[[1]]
  otherData.list <- sapply(edge.files, function(f)fread (file.path(otherEdgeDir, f) ))
}else{
  otherData.list <- list()
}

otherData.sliderColumn <- lapply(otherData.list, function(dt)grep("^weight\\.", colnames(dt), value = TRUE)[1])

makeSliderOtherData <- function (name){
  print (name)
  if(!is.na(otherData.sliderColumn[[name]])){
    quants <- quantile( otherData.list[[name]][[otherData.sliderColumn[[name]]]], c(0,0.9, 1.0))
    print(quants)
    return (shiny::sliderInput(sprintf("slider.%s", name), 
                               min = quants[1], max = quants[3], value = quants[2], 
                               label = sprintf ("%s : %s", 
                                                name, 
                                                gsub("^weight.",
                                                     "",
                                                     otherData.sliderColumn[[name]]) )))
  }else{
    return(NULL)
  }
}




colors <- c("purple", "#ff7777", "forestgreen", "orange", "blue", "magenta")
otherData.colors <- colors[1:length(otherData.list)]
names(otherData.colors) <- names(otherData.list)
edgeColors <- c(c(apms = "gray30", apms2ndDegree = "gray30"),
                otherData.colors)

# corum.dt <- fread ("./data/corum.dt.csv")
# corum.sets.dt <- corum.dt[, .(gene = unique(c(gene.x, gene.y))), by= name][, complexSize := length(unique(gene)), by = name][]


#alphaFold <- fread ("")




# ---- functions ----
# sub network functions

makeEdgeTable <- function(baitsOI, allBaits, baitPrey.dt, otherData.list){
  
  subEdges.dt <- baitPrey.dt[bait %in% baitsOI]
  
  preySet <- setdiff(unique(subEdges.dt$prey), baitsOI)
  
  allNodes <- c(preySet, baitsOI)
  
  setnames(subEdges.dt, new = c("gene.x", "gene.y"), old = c("bait", "prey"))
  
  # secondary edges bait1->bait2->preyOfBait1
  # so basically any prey-prey edge that is already in the baitPrey.dt table
  secondaryAPMS <- baitPrey.dt[bait %in% subEdges.dt$gene.y & prey %in% subEdges.dt$gene.y]
  setnames(secondaryAPMS, new = c("gene.x", "gene.y"), old = c("bait", "prey"))

  # subset based on nodes, reorder edge direction, take unique only  
  otherData.subSet <-lapply (otherData.list, function(dt)unique(dt[gene.x %in% allNodes & gene.y %in% allNodes,
                                                                   .(gene.x = ifelse(gene.x < gene.y, gene.x, gene.y),
                                                                     gene.y = ifelse(gene.x < gene.y, gene.y, gene.x))]) )
  
  
  
  combinedEdges.dt <- rbindlist(c(list(apms = subEdges.dt,
                                       apms2ndDegree = secondaryAPMS),
                                  otherData.subSet),
                                use.names = TRUE, fill = TRUE,
                                idcol = "edge.source")
  
  
  combinedEdges.dt[, edgeColor := edgeColors[edge.source]]
  combinedEdges.dt[, edgeLineType := ifelse (edge.source == "apms", "solid", "dashed")] #:= c(apms = "solid", corum = "dashed", apms2ndDegree = "dashed")[edge.source]]
  
  setcolorder(combinedEdges.dt, c("gene.x", "gene.y"))
  return (combinedEdges.dt[])  
}


makeNodeTable <- function(subEdges.dt, baitsOI, allBaits, allPreys){
  nodes.dt <- data.table(gene = unique(c(subEdges.dt$gene.x, subEdges.dt$gene.y)))
  
  nodes.dt[gene %in% baitsOI, type := "bait"] # main bait
  nodes.dt[is.na(type) & gene %in% allBaits, type := "other.bait"]
  nodes.dt[is.na(type) & gene %in% allPreys, type := "prey"]
  
  # format nodes by type
  nodes.dt[, shape := c(bait = "diamond", other.bait = "diamond open", prey = "circle small")[type]]
  nodes.dt[, nodeSize := c(bait = 6, other.bait = 4, prey = 3)[type]]
  nodes.dt[, color := c(bait = "black", other.bait = "black", prey = "red")[type]]
  
  return(nodes.dt[])
}


makeSubGraphPlot <- function(g.layout, labelSize){

  
  #g.layout <- graph.reactive()
  
  # define limits for coord_fixed bug...
  limits <- range(c(g.layout$x, g.layout$y))
  
  
  
  p <- ggraph::ggraph(g.layout) + 
    #ggraph::geom_edge_link(aes(color = edgeColor, lty = edgeLineType), alpha = 0.5) +
    ggraph::geom_edge_fan(aes(color = edgeColor, lty = edgeLineType), edge_width = 0.5) + # don't do fans until I deal with 2nd degree edges on multi-baits
    ggraph::geom_node_point(aes(shape = shape, color = color, size = nodeSize), alpha = 0.5) +
    scale_shape_identity() +
    scale_size_identity() +
    scale_color_identity() +
    ggraph::scale_edge_color_identity() +
    ggraph::scale_edge_linetype_identity() + 
    ggraph::geom_node_text(aes(label = name, color = color), repel = TRUE, size = labelSize, point.padding = max(g.layout$nodeSize)) +
    theme_void() +
    coord_fixed(xlim = limits, ylim = limits)
  #coord_fixed() # causes an error when range(y) is vyer msal, like -2.424677e-17  6.234884e-17
  p
  
}

# ---- ui ----
# Define UI for application that draws a histogram
ui <- fluidPage(
  
  # Application title
  titlePanel("APMS Subnetwork Viewer"),
  
  sidebarLayout(
    sidebarPanel(
      # tags$head(tags$script('$(document).on("shiny:connected", function(e) {
      #                     Shiny.onInputChange("innerHeight", window.innerHeight);
      #                     });
      #                     $(window).resize(function(e) {
      #                     Shiny.onInputChange("innerHeight", window.innerHeight);
      #                     });
      #                     ')),
      
      selectInput(inputId = "baitsOI",
                  label = "Enter baits to show",
                  choices= allBaits,
                  selected = sample(allBaits, 3),
                  multiple=TRUE, 
                  selectize=TRUE),
      selectInput(inputId = "layout",
                  label = "Network layout algorithm",
                  choices = c("kk", "fr", "stress", "dh", "lgl", "nicely"),
                  multiple = FALSE),
      checkboxInput(inputId =  "hideAPMS",
                    label = "hide main APMS edges",
                    value = FALSE),
      shiny::h3("Other Edge Types:"),
      uiOutput("dataSelectors"),
      shiny::hr(),
      shiny::h2("Other settings:"),
      sliderInput(inputId =  "labelSize",label = "Label Size", min = 1, max = 10, value = 3.5, step = 0.1),
      
      downloadLink("downloadData", "Download Network Image")
      
    ),
    
    
    # Show a plot of the generated distribution
    shiny::mainPanel(
      #plotOutput("networkPlot", width = "100%")
      fillPage(
        tags$style(type = "text/css", "#networkPlot {height: calc(100vh - 80px) !important;}"),
        plotOutput("networkPlot", width = "100%",height = "100%")
      )
    )
  )
)
# ---- server ----

# Define server logic required to draw a histogram
server <- function(input, output) {
  
  # set up the data-type selectors
  # doesn't yet need to live here in the server.
  output$dataSelectors <- renderUI({
    ui.list <- list()
    for (n in names(otherData.list)){
      ui.list[[length(ui.list) + 1]] <-  checkboxInput(paste0("select.", n), n)
      ui.list[[length(ui.list) + 1]] <-  colourpicker::colourInput(paste0("color.", n),label = NULL, value = otherData.colors[n], width = 15)
      slider <- makeSliderOtherData(n)
      if (!is.null(slider))
        ui.list[[length(ui.list) + 1]] <-  slider
      ui.list[[length(ui.list) + 1]] <-  shiny::hr()
    }
    ui.list
    
    # lapply(names(otherData.list),function(n){
      # list(
      #   checkboxInput(paste0("select.", n), n),
      #   makeSliderOtherData(n))
         #makeSliderOtherData(n)
      # list(
      # ,
      # shiny::textOutput("dummy")
      # )
      
    # })
  })
  

  combinedEdges.reactive <- reactive({
    makeEdgeTable(input$baitsOI, allBaits, baitPrey.dt, otherData.reactive())
  })
  
  otherData.reactive <- reactive({
    od <- list()
    for (n in names(otherData.list)){
      nn <- paste0("select.", n)
      if (input[[nn]] == TRUE){
        
        # check for slider
        sliderValue <- input[[sprintf("slider.%s", n)]]
        if(!is.null(sliderValue)){
          subSet <- otherData.list[[n]][[otherData.sliderColumn[[n]]]] > sliderValue
          od[[n]] <- otherData.list[[n]][subSet == TRUE,]
          
        }else{
          od[[n]] <- otherData.list[[n]]
        }
      }
    }
    od
  })
  
  nodes.reactive <- reactive({
    makeNodeTable(combinedEdges.reactive(), input$baitsOI, allBaits, allPreys)
  })
  
  plot.reactive <- reactive({
    makeSubGraphPlot(graph.reactive(),
                     labelSize = input$labelSize
    )
    
  })
  
  edgeColors.reactive <- reactive({
    otherDataNames <- names(otherData.reactive())
    #get colors from color pickers
    colors <- sapply(otherDataNames, function(odn)input[[sprintf ("color.%s", odn)]])
    names(colors) <- otherDataNames
    return (colors)
    
  })
  
  
  graph.reactive <- reactive({
    edges <- combinedEdges.reactive()
    nodes <- nodes.reactive()
    edgeColorMap <- edgeColors.reactive()
    for(odn in names(edgeColorMap))
      edges[edge.source == odn, edgeColor := edgeColorMap[odn]]
    g <- igraph::graph_from_data_frame(edges,
                                       vertices = nodes,
                                       directed = FALSE) # not fully undirected, but for simplicity sake 
    
    if (input$hideAPMS){
      g <- igraph::delete.edges(g, which(igraph::get.edge.attribute(g, "edge.source") == "apms"))       
    }
    
    
    g.layout <- ggraph::create_layout(g, layout = input$layout)
    
    # center
    g.layout$x <- g.layout$x - mean(range(g.layout$x))
    g.layout$y <- g.layout$y - mean(range(g.layout$y))
    
    g.layout
    
  })
  
    output$networkPlot <- renderPlot({
      plot.reactive()
      }
    )
    
    
    output$downloadData <- downloadHandler(
      filename = function() {
        paste("PPI_", paste0(input$baitsOI, collapse = "_"), ".pdf", sep="")
      },
      content = function(file) {
        pdf(file=file, width = 10, height = 10) #res = 300, height = 6, width = 6, units = "in")
        print(plot.reactive())
        dev.off()
      }
    )
}

# Run the application 
shinyApp(ui = ui, server = server)
