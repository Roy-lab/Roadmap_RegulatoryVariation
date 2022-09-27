library(shiny)
library(tools)
library(tidyverse)
library(networkD3)
library(Matrix)
library(scales)

## Load in aux functions and data
load('phenotypeFiles_binaryInput_fixID/Coronaryarterydisease.Rdata')
source('aux_functions2.R')
MyClickScript <- 'Shiny.setInputValue("node_name", d.name)'
cell_line_file=read_tsv("cell_lines_nocom.txt", col_names=FALSE)
phenotype_file=read_tsv("phenos.txt", col_names=FALSE)
#netList=list(Net, TGS, Module, SNPs,transitioning_sets)
################################# ui ###########################################
ui <- fluidPage(
    title = "CD56 TEST EXAMPLE",
    forceNetworkOutput("network"),

    fluidRow(
      column(4,
             h4("Search Network"),
             selectInput(inputId = "method", "Search Method",
                  c("Transitioning Gene Set" = "geneset",Modules = "module")
             ),
             # Only show this panel if the plot type is a histogram
             conditionalPanel(
              condition = "input.method == 'module'",
                selectInput(
                inputId = "module_id", "Module ID",
                sort(module_ids)
                ),
              ),
	     conditionalPanel(
              condition = "input.method == 'geneset'",
                selectInput(
                inputId = "transitioning_set", "Transitioning Gene Set",
                sort(transitioning_sets)
                )
              )
	),
    column(4,
              h4("SNPs linked to Gene"), 
              htmlOutput(outputId = "node_info")
      ),
    column(4,
             h4("Change Cell Type"),
             selectInput(inputId = "celltype", "Cell type",
                  cell_line_file$X1
          ),
             h4("Change Phenotype"),
                selectInput(inputId = "phenotype", "Phenotype",
                  phenotype_file$X1
          )
        )
)
)

server <- function(input, output, session) {
   updateSelectizeInput(session, 'gene', choices = genes, server = TRUE)
   #### Variables ################
   node_name_info <- reactiveVal(value = NA)
  module_id_info <- reactiveVal(value = NA)
   gene_name <-reactiveVal()
   gene_list <- reactiveVal()
 #  ph='Coronaryarterydisease'
 # load("phenotypeFiles_binaryInput_fixID/Coronaryarterydisease.Rdata")
#     react <- reactiveValues(Net = Net, TGS = TGS,Module = Module,SNPs = SNPs,transitioning_sets = transitioning_sets) 
    react <- reactiveValues('foo' = list(Net, TGS,Module,SNPs,transitioning_sets))

    observeEvent(input$phenotype,{
      rm(Net,TGS,SNPs,Module) 
      react$foo = NULL
      print(input$phenotype)
      req(input$phenotype)
      filename <- paste0('phenotypeFiles_binaryInput_fixID/', input$phenotype, ".Rdata")
      load(filename)
      updateSelectInput(session, "transitioning_set",
                           label = "Transitioning Gene Set",
                           choices = transitioning_sets)
      react$foo[[1]] <- Net
      react$foo[[2]] <- TGS
      react$foo[[3]] <- Module
      react$foo[[4]] <- SNPs
      react$foo[[5]] <- transitioning_sets
      },ignoreInit = TRUE,ignoreNULL=FALSE)	
 
   sub_net <- reactive({
   rm(sub_net)
   if(input$method=="module"){
      module <- input$module_id
      print(module)
      module_id_info(input$module_id)
      node_name_info(NA)
      #moduleSubgraph(netList[[1]], netList[[3]], module,input$celltype)
      moduleSubgraph(react$foo[[1]], react$foo[[3]], module,input$celltype)
    }else if (input$method=="geneset"){
      gene_set <- input$transitioning_set
      print(gene_set)
      #print(react$Net)
      module_id_info(input$transitioning_set)
      node_name_info(NA)
      #genesetSubgraph(netList[[1]], netList[[2]], gene_set,input$celltype)
      genesetSubgraph(react$foo[[1]], react$foo[[2]], gene_set,input$celltype)
   }
   #HTML(printSNPInfo(netList[[4]],input$node_name))
 })

   ############## Main Render ######################
  output$network <- renderForceNetwork({
  print(input$celltype)
  #doesn't print on the third time 
  #issue where sub_net is called
  #print(react$Net)
  S<-sub_net()
  #no issue
  S_tables <- graph2NodeEdgeTables(S,input$celltype)
  S_nodes <- S_tables[[1]]
  S_edges <- S_tables[[2]]
  #print(S_nodes)
  #print(S_edges)
  #ColInds = which(S_nodes$Input > 0, arr = TRUE)[, 1]
  #NodeCols = ifelse(1:nrow(MisLinks) %in% ColInds, "#bf3eff", "#666")
  ### Add superficial edge to overcome empty array 
  S_edges <- S_edges %>% add_row(from = 0, to = 0)
  ColourScale <- 'd3.scaleOrdinal()
            .domain([0, 1])
           .range(["#FF6900", "#694489"]);'
  S_nodes$Diffused = S_nodes$Diffused * 5;
  forceNetwork(Links = S_edges, Nodes = S_nodes,
                 Source = "from", Target = "to",
                 NodeID = "feature",Nodesize='Diffused',
                 Group = "Input", opacity = 1, zoom = TRUE, fontSize=40,
		colourScale = ColourScale,
                 charge = -10, clickAction = MyClickScript)
  })

  ############### observe Events ##########################################
  output$node_info <- renderUI({
        #netList <- load_new()
	#HTML(printSNPInfo(netList[[4]],input$node_name,input$celltype))
  	HTML(printSNPInfo(isolate(react$foo[[4]]),input$node_name,input$celltype))
  })
}

shinyApp(ui, server)
