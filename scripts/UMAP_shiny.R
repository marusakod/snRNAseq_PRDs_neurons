# SHINY APP FOR EXPLORING EXPRESSION OFF DIFFERENT GENES IN snRNAseq DATASET

library(shiny)
library(esquisse)
library(Seurat)
library(SingleCellExperiment) 
library(shinyWidgets)
library(tidyverse)
library(shinyjs)

all_SCE_objects <- readRDS("all_SCE_objects.rds")
NBH_noHashed_Seurat_combined <- readRDS("NBH_noHashed_Seurat_combined.rds")
RML_noHashed_Seurat_combined <- readRDS("RML_noHashed_Seurat_combined.rds")  
NBH_noHashed_Seurat_comb_filt <- readRDS("NBH_noHashed_Seurat_comb_filt.rds")
RML6_noHashed_Seurat_comb_filt <- readRDS("RML6_noHashed_Seurat_comb_filt.rds")

# append Replicate column to metadata of nonhashed samples
addReplicate_to_Seurat <- function(seurat){
  # add metadata column assigning cells to replicate 1 or replicate 2
  cells <- rownames(seurat@meta.data)
  replicate_vec <- vector("character", length = length(cells))
  replicate_vec[grepl("Replicate1.*", cells)] <- "Replicate 1"
  replicate_vec[grepl("Replicate2.*", cells)] <- "Replicate 2"
  
  rep_df <- data.frame(Replicate = replicate_vec)
  rownames(rep_df) <- cells
  
  s1 <- AddMetaData(seurat, metadata = rep_df)
  s1
}

all_noHashed_Seurat_comb_filt <- sapply(list(NBH_noHashed = NBH_noHashed_Seurat_comb_filt, RML_noHashed = RML6_noHashed_Seurat_comb_filt),
                                        FUN = addReplicate_to_Seurat, 
                                        simplify = FALSE)


all_Hashed_demultiplexed <- readRDS("all_Hashed_demultiplexed.rds")
all_Hashed_demultiplexed_filt2x <- readRDS("all_Hashed_demultiplexed_filt2x.rds")

# remove doublets and negatives from nonfiltered Hashed Seurat objects

all_HTO_names <- rownames(all_Hashed_demultiplexed$RML6_Hashed[["HTO"]])

# remove all negatives and doublets and assign cells to replicates 
filter_dumux_seurat <- function(demux_seurat){
  s1 <- subset(x = demux_seurat, subset = hash.ID %in% all_HTO_names)
  
  # add metadata column assigning cells to replicate 1 or replicate 2
  singlets <- rownames(s1@meta.data)
  singlets_hashIDs <- s1@meta.data$hash.ID
  replicate_vec <- vector("character", length = length(singlets_hashIDs))
  replicate_vec[grepl(".*REP1", singlets_hashIDs)] <- "Replicate 1"
  replicate_vec[grepl(".*REP2", singlets_hashIDs)] <- "Replicate 2"
  
  rep_df <- data.frame(Replicate = replicate_vec)
  rownames(rep_df) <- singlets
  
  s1 <- AddMetaData(s1, metadata = rep_df)
  s1
}

all_Hashed_demultiplexed <- sapply(all_Hashed_demultiplexed, FUN = filter_dumux_seurat, simplify = FALSE)


# gene Symbols for selectInputs
all_genes <- rowData(all_SCE_objects$NBH_1_NoHashed)$Symbol
all_genes_Hashed <- rowData(all_SCE_objects$NBH_Hashed)$Symbol

# USER INTERFACE ==========================================================================
ui <- fluidPage(sidebarLayout(
  sidebarPanel = sidebarPanel(width = 2, 
                              radioButtons(
                                inputId = "pickSample",
                                label = "Sample type",
                                choices = c("Nonhashed", "Hashed"),
                                selected = "Nonhashed"
                               
                              ),
                              conditionalPanel("input.pickSample == 'Nonhashed'", 
                              selectizeInput("geneSymbol",
                                             label = "Select Gene Symbol", 
                                             multiple = FALSE,
                                             size = 20,
                                             choices = NULL,
                                             width = "200px",
                                             options = list(
                                               placeholder = 'eg. Agap2',
                                               onInitialize = I('function() { this.setValue("");}')
                                             ))),
                              conditionalPanel("input.pickSample == 'Hashed'", 
                                               selectizeInput("geneSymbolHashed",
                                                              label = "Select Gene Symbol", 
                                                              multiple = FALSE,
                                                              size = 20,
                                                              choices = NULL,
                                                              width = "200px",
                                                              options = list(
                                                                placeholder = 'eg. Agap2',
                                                                onInitialize = I('function() { this.setValue("");}')
                                                              ))),
                
                              radioButtons("colorWhat", 
                                                "Colour cells by:", 
                                                choices = c("Normalized counts", "UMI counts (raw)")),
                              radioButtons("replicatesDisp",
                                           "Display replicates:",
                                           choices = c("merged", "split")),
                              esquisse::colorPicker("colors",
                                                    choices = list(grey1 = "#EAECEE",
                                                                   grey2 = "#BFC9CA", 
                                                                   yellow = "#F7DC6F",
                                                                   orange = "#F39C12",
                                                                   brown = "#A04000", 
                                                                   lightGreen = "#56EB4F",
                                                                   darkGreen = "#196F3D", 
                                                                   lightBlue = "#4FCFEF", 
                                                                   darkBlue = "#2874A6", 
                                                                   lightPurple = "#D2B4DE",
                                                                   darkPurple = "#8E44AD",
                                                                   pink = "#F596E6",
                                                                   red1 = "#F00B1C",
                                                                   red2 = "#8D161F"),
                                                    label = "Pick color palette",
                                                    width = "200px", 
                                                    selected = list(grey1 = "#EAECEE", yellow = "#F7DC6F", orange = "#F39C12",red1 = "#F00B1C"),
                                                    multiple = TRUE,
                                                    textColor = "#000",
                                                    plainColor = FALSE),
                              numericInput("pointSize", 
                                           label = "Select point size", 
                                           value = 1,
                                           min = 0.1, 
                                           max = 10, 
                                           step = 0.1, 
                                           width = "200px"),
                              numericInput("alpha",
                                           label = "Select point transparency", 
                                           value = 0.5, 
                                           min = 0,
                                           max = 1,
                                           step = 0.1, 
                                           width = "200px"),
                              downloadButton(outputId = "downloadNBH",
                                           label = "Download NBH plot",
                                           icon = icon("download")),
                             downloadButton(outputId = "downloadRML",
                                                        label = "Download RML plot",
                                                        icon = icon("download"))),
  mainPanel = mainPanel(width = 10, 
                        uiOutput("get_plot_NBH"), 
                        uiOutput("get_plot_RML"))
  
))


# SERVER ======================================================================================================================

server <- function(input, output, session){
  
  updateSelectizeInput(session  = session,
                       inputId  = "geneSymbol",
                       server   = TRUE,
                       choices  = all_genes,
                       selected = character(0),
                       options = list(
                         placeholder = 'eg. Agap2',
                         onInitialize = I('function() { this.setValue("");}')
                       ))
  
  updateSelectizeInput(session  = session,
                       inputId  = "geneSymbolHashed",
                       server   = TRUE,
                       choices  = all_genes_Hashed,
                       selected = character(0),
                       options = list(
                         placeholder = 'eg. Agap2',
                         onInitialize = I('function() { this.setValue("");}')
                       ))
  
  observeEvent(input$pickSample, {
    shinyjs::reset("geneSymbol")
    shinyjs::reset("geneSymbolHashed")
  })
  
  plot_width <- reactive({
    if(input$replicatesDisp == "merged"){
      w <- "700px"
    }else{
      w <- "1400px"
    }
    w
  })
  

  Entry_variable <- reactiveVal()
  
  observeEvent(input$geneSymbol, {
    Entry_variable(input$geneSymbol)
  }) # manual single gene selection (symbol) #
  
  observeEvent(input$geneSymbolHashed, {
    Entry_variable(input$geneSymbolHashed)
  }) # manual single gene selection (Hashed symbol) 
  
  # get ensembl id for selected gene symbol
  ensembl <- reactive({
    if(input$pickSample == "Nonhashed"){
      data <- rowData(all_SCE_objects$NBH_1_NoHashed)
    }else{
      data <- rowData(all_SCE_objects$NBH_Hashed)
    }
    
    ensembl <- as.data.frame(data) %>% filter(Symbol == Entry_variable()) %>% dplyr::select(ID) %>% flatten_chr()
    ensembl
  })
  
  
  # plot expression of a particular gene to a UMAP
  make_UMAP_df <- function(ensembl, unfilt_seurat, filt_seurat){
    # get normalized counts and raw counts from unfiltered seurat object
    raw_counts <-  GetAssayData(unfilt_seurat, slot = "counts")[ensembl, ]
    norm_counts <- GetAssayData(unfilt_seurat, slot = "data")[ensembl, ]
    
    # append counts to UMAP dataframe
    UMAP_df <- as.data.frame(filt_seurat[["umap"]]@cell.embeddings)
    UMAP_df$replicate <-  filt_seurat@meta.data$Replicate
    UMAP_df$raw_counts <- raw_counts
    UMAP_df$norm_counts <- norm_counts
    
    UMAP_df
  }

  
NBH_filt_Seurat <- reactive({
  if(input$pickSample == "Nonhashed"){
    x <- all_noHashed_Seurat_comb_filt$NBH_noHashed
  }else{
    x <- all_Hashed_demultiplexed_filt2x$NBH_Hashed
  }
  x
})

RML_filt_Seurat <- reactive({
  if(input$pickSample == "Nonhashed"){
    x <- all_noHashed_Seurat_comb_filt$RML_noHashed
  }else{
    x <- all_Hashed_demultiplexed_filt2x$RML6_Hashed
  }
  x
})

NBH_unfilt_seurat <- reactive({
  if(input$pickSample == "Nonhashed"){
    x <- NBH_noHashed_Seurat_combined
  }else{
    x <- all_Hashed_demultiplexed$NBH_Hashed
  }
  x
  
})

RML_unfilt_seurat <- reactive({
  if(input$pickSample == "Nonhashed"){
    x <- RML_noHashed_Seurat_combined
  }else{
    x <- all_Hashed_demultiplexed$RML6_Hashed
  }
  x
})

NBH_UMAP_df <- reactive({
  df <- make_UMAP_df(ensembl = ensembl(), 
                     unfilt_seurat = NBH_unfilt_seurat(), 
                     filt_seurat = NBH_filt_Seurat())
  
  # if Hashed samples are to be displayed add hash ID to UMAP_df
  if(input$pickSample == "Hashed"){
    df$hashID <- NBH_filt_Seurat()@meta.data$hash.ID
    return(df)
  }else{
    return(df)
  }
})

RML_UMAP_df <- reactive({
  df <- make_UMAP_df(ensembl = ensembl(), 
                     unfilt_seurat = RML_unfilt_seurat(), 
                     filt_seurat = RML_filt_Seurat())
  
  if(input$pickSample == "Hashed"){
    df$hashID <- RML_filt_Seurat()@meta.data$hash.ID
    return(df)
  }else{
    return(df)
  }
})

# variable determining what should be colored 
inColor <- reactive({
  if(input$colorWhat == "Normalized counts"){
    x <- "norm_counts"
  }else{
    x <- "raw_counts"
  }
  x
})



make_UMAP <- function(UMAP_df, p_size, trans, which_color, Title, cols, Subtit, legendName){
    
    p <- ggplot(UMAP_df, aes(x = UMAP_1, y = UMAP_2, color = get(which_color))) +
      geom_point(alpha = trans, size = p_size) +
      theme_light() +
      labs(x = "UMAP 1", y = "UMAP 2", title = Title, subtitle = Subtit, color = legendName) +
      theme(panel.grid.major = element_blank(),
            strip.text = element_text(size = 16, colour = "black"),
            axis.text = element_text(size = 12), 
            axis.title = element_text(size = 14), 
            legend.text = element_text(size = 12), 
            legend.title = element_text(size = 14), 
            plot.title = element_text(size = 18),
            plot.subtitle = element_text(size = 16), 
            panel.grid.minor = element_blank()) +
      scale_color_gradientn(colours = cols)
    
    p

}

empty_umap <- ggplot(data.frame(x = c(1,2,3), y = c(1,2,3)), aes(x = x, y = y)) +
              geom_blank() +
              labs(x = "UMAP 1", y = "UMAP 2") +
              theme_light() +
  theme(panel.grid.major = element_blank(),
        axis.text = element_text(size = 12), 
        axis.title = element_text(size = 14),
        panel.grid.minor = element_blank())

#vals <- reactiveValues()

NBH_UMAP <- reactive({
  if(!ensembl() %in% rownames(GetAssayData(NBH_filt_Seurat(), slot = "counts"))){
    note <- "This gene has been FILTERED OUT!!!"
  }else{
    note <- ""
  }
  
  if(input$pickSample == "Nonhashed"){
    sub <- input$geneSymbol
  }else{
    sub <- input$geneSymbolHashed
  }
  
  p <-  make_UMAP(UMAP_df = NBH_UMAP_df(),
                  p_size = input$pointSize,
                  trans = input$alpha,
                  which_color = inColor(),
                  Title = "NBH",
                  legendName = input$colorWhat, 
                  cols = input$colors, 
                  Subtit = paste(sub, note, sep = " "))
  
  if(input$replicatesDisp == "split"){
    p <- p + facet_wrap(~replicate)
    #vals$NBH <- p
    return(p)
  }else{
    #vals$NBH <- p
    return(p)
  }
})


RML_UMAP <- reactive({
  if(!ensembl() %in% rownames(GetAssayData(RML_filt_Seurat(), slot = "counts"))){
  note <- "This gene has been FILTERED OUT!!!"
}else{
  note <- ""
}
  
  if(input$pickSample == "Nonhashed"){
    sub <- input$geneSymbol
  }else{
    sub <- input$geneSymbolHashed
  }

p <-  make_UMAP(UMAP_df = RML_UMAP_df(),
                p_size = input$pointSize,
                trans = input$alpha,
                which_color = inColor(),
                Title = "RML",
                legendName = input$colorWhat, 
                cols = input$colors,
                Subtit = paste(sub, note, sep = " "))

if(input$replicatesDisp == "split"){
  p <- p + facet_wrap(~replicate)
  #vals$RML <- p
  return(p)
}else{
  #vals$RML <- p 
  return(p)
}

})




output$NBH <- renderPlot({
  
  if(input$pickSample == "Nonhashed"){
    if(length(ensembl()) == 0){
      return(empty_umap)
    }else{
      return(NBH_UMAP())
    }
  }else{
  if(length(ensembl()) == 0){
   return(empty_umap)
  }else{
    return(NBH_UMAP())
  }
  }
})



output$RML <- renderPlot({
  
  if(input$pickSample == "Nonhashed"){
    if(length(ensembl()) == 0){
      return(empty_umap)
    }else{
      return(RML_UMAP())
    }
  }else{
    if(length(ensembl()) == 0){
      return(empty_umap)
    }else{
      return(RML_UMAP())
    } 
  }
  })

 # output$downloadNBH <- downloadHandler(
  #  filename = function(){paste(input$geneSymbol, "_NBH_UMAP", '.png', sep = "")},
    
  #  content = function(file){
  #    png(file, width = as.numeric(gsub("px", "", plot_width())), height = 700, units = "px")
    #  print(vals$NBH)
   #   dev.off()
   # }
 # )
  
  
  #output$downloadRML <- downloadHandler(
   # filename = function(){paste(input$geneSymbol, "_RML_UMAP", '.png', sep = "")},
    
   # content = function(file){
    #  png(file, width = as.numeric(gsub("px", "", plot_width())), height = 700, units = "px")
  #    print(vals$RML)
  #    dev.off()
  #  }
#  )
  


output$get_plot_NBH <- renderUI({
    
    plotOutput("NBH", width = plot_width(), height = "700px")
    
    
  })
  
output$get_plot_RML <- renderUI({
    
    plotOutput("RML", width = plot_width(), height = "700px")
    
  })
  
}

shinyApp(ui = ui, server = server)

