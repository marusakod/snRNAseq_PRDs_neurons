# SHINY APP FOR EXPLORING EXPRESSION OFF DIFFERENT GENES IN snRNAseq DATASET

library(shiny)
library(esquisse)
library(Seurat)
library(SingleCellExperiment) # install with BiocManager::install("SingleCellExperiment")
library(shinyWidgets)
library(tidyverse)

all_SCE_objects <- readRDS("all_SCE_objects.rds")
NBH_noHashed_Seurat_combined <- readRDS("NBH_noHashed_Seurat_combined.rds")
RML_noHashed_Seurat_combined <- readRDS("RML_noHashed_Seurat_combined.rds")  
NBH_noHashed_Seurat_comb_filt <- readRDS("NBH_noHashed_Seurat_comb_filt.rds")
RML6_noHashed_Seurat_comb_filt <- readRDS("RML6_noHashed_Seurat_comb_filt.rds")


all_genes <- rowData(all_SCE_objects$NBH_1_NoHashed)$Symbol

# USER INTERFACE ==========================================================================
ui <- fluidPage(sidebarLayout(
  sidebarPanel = sidebarPanel(width = 2, 
                              selectizeInput("geneSymbol",
                                             label = "Select Gene Symbol", 
                                             multiple = FALSE,
                                             size = 20,
                                             choices = NULL,
                                             width = "200px",
                                             options = list(
                                               placeholder = 'eg. Agap2',
                                               onInitialize = I('function() { this.setValue("");}')
                                             )),
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
  
  plot_width <- reactive({
    if(input$replicatesDisp == "merged"){
      w <- "700px"
    }else{
      w <- "1400px"
    }
    w
  })
  
  # get ensembl ID for selected gene symbol 
  ensembl <- reactive({
    as.data.frame(rowData(all_SCE_objects$NBH_1_NoHashed)) %>% filter(Symbol == input$geneSymbol) %>% dplyr::select(ID) %>% flatten_chr()
  })
  
  
  # plot expression of a particular gene to a UMAP
  make_UMAP_df <- function(ensembl, unfilt_seurat, filt_seurat){
    # get normalized counts and raw counts from unfiltered seurat object
    raw_counts <-  GetAssayData(unfilt_seurat, slot = "counts")[ensembl, ]
    norm_counts <- GetAssayData(unfilt_seurat, slot = "data")[ensembl, ]
    
    # append counts to UMAP dataframe
    UMAP_df <- as.data.frame(filt_seurat[["umap"]]@cell.embeddings)
    UMAP_df$replicate <-  gsub("_cell_.*", "", rownames(UMAP_df))
    UMAP_df$raw_counts <- raw_counts
    UMAP_df$norm_counts <- norm_counts
    
    UMAP_df
  }


NBH_UMAP_df <- reactive({
  df <- make_UMAP_df(ensembl = ensembl(), 
                     unfilt_seurat = NBH_noHashed_Seurat_combined, 
                     filt_seurat = NBH_noHashed_Seurat_comb_filt)
  df
})

RML_UMAP_df <- reactive({
  df <- make_UMAP_df(ensembl = ensembl(), 
                     unfilt_seurat = RML_noHashed_Seurat_combined, 
                     filt_seurat = RML6_noHashed_Seurat_comb_filt)
  df
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

make_UMAP <- function(UMAP_df, p_size, trans, which_color, Title, cols, Subtit){
    
    p <- ggplot(UMAP_df, aes(x = UMAP_1, y = UMAP_2, color = get(which_color))) +
      geom_point(alpha = trans, size = p_size) +
      theme_light() +
      labs(x = "UMAP 1", y = "UMAP 2", title = Title, subtitle = Subtit) +
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

vals <- reactiveValues()

output$NBH <- renderPlot({
  if(identical(input$geneSymbol, character(0))){
    vals$NBH <- empty_umap
    return(empty_umap)
  }else{
    
  # if ensemblID is is not in the filtered seurat object add notification to plot title
  if(!ensembl() %in% rownames(GetAssayData(NBH_noHashed_Seurat_comb_filt, slot = "counts"))){
    note <- "This gene has been FILTERED OUT!!!"
  }else{
    note <- ""
  }
    
  p <-  make_UMAP(UMAP_df = NBH_UMAP_df(),
                  p_size = input$pointSize,
                  trans = input$alpha,
                  which_color = inColor(),
                  Title = "NBH",
                  cols = input$colors, 
                  Subtit = paste(input$geneSymbol, note, sep = " "))
  
  if(input$replicatesDisp == "split"){
    p <- p + facet_wrap(~replicate)
    vals$NBH <- p
    return(p)
  }else{
    vals$NBH <- p
    return(p)
  }
  
  } 
})
    
output$RML <- renderPlot({
  if(identical(input$geneSymbol, character(0))){
    vals$RML <- empty_umap
    return(empty_umap)
  
  }else{
    
    # if ensemblID is is not in the filtered seurat object add notification to plot title
    if(!ensembl() %in% rownames(GetAssayData(RML6_noHashed_Seurat_comb_filt, slot = "counts"))){
      note <- "This gene has been FILTERED OUT!!!"
    }else{
      note <- ""
    }
    
    p <-  make_UMAP(UMAP_df = RML_UMAP_df(),
                    p_size = input$pointSize,
                    trans = input$alpha,
                    which_color = inColor(),
                    Title = "RML",
                    cols = input$colors,
                    Subtit = paste(input$geneSymbol, note, sep = " "))
    
    if(input$replicatesDisp == "split"){
      p <- p + facet_wrap(~replicate)
      vals$RML <- p
      return(p)
    }else{
      vals$RML <- p 
      return(p)
    }
    
  } 
  
  output$downloadNBH <- downloadHandler(
    filename = function(){paste(input$geneSymbol, "_NBH_UMAP", '.png', sep = "")},
    
    content = function(file){
      png(file, width = as.numeric(gsub("px", "", plot_width())), height = 700, units = "px")
      print(vals$NBH)
      dev.off()
    }
  )
  
  
  output$downloadRML <- downloadHandler(
    filename = function(){paste(input$geneSymbol, "_RML_UMAP", '.png', sep = "")},
    
    content = function(file){
      png(file, width = as.numeric(gsub("px", "", plot_width())), height = 700, units = "px")
      print(vals$RML)
      dev.off()
    }
  )
  
})


output$get_plot_NBH <- renderUI({
    
    plotOutput("NBH", width = plot_width(), height = "700px")
    
    
  })
  
output$get_plot_RML <- renderUI({
    
    plotOutput("RML", width = plot_width(), height = "700px")
    
  })
  
}

shinyApp(ui = ui, server = server)

