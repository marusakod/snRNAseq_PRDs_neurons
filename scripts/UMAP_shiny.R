# SHINY APP FOR EXPLORING EXPRESSION OFF DIFFERENT GENES IN snRNAseq DATASET

library(shiny)
library(esquisse)
library(Seurat)
library(SingleCellExperiment) 
library(shinyWidgets)
library(tidyverse)
library(shinyjs)

all_SCE_objects <- readRDS("all_SCE_objects.rds") # just for ensembl IDs

ensembl_symbol_table <- as.data.frame(rowData(all_SCE_objects$NBH_1_NoHashed)) %>% dplyr::select(ID, Symbol)

# all samples merged seurat objects
all_srat_Merged_filtered <- readRDS("all_srat_Merged_filtered2.rds")

# all samples merged ambient corrected seurat objects
all_srat_Merged_ambCorr_filtered <- readRDS("all_srat_Merged_ambCorr_filtered2.rds")

# make a UMAP df 

make_UMAP_df  <- function(srat){
  umap_df <- as.data.frame(srat[["umap"]]@cell.embeddings)
  umap_df2 <- cbind(umap_df, srat@meta.data)
  
  # remove columns that start with pANN
  to_remove <- colnames(umap_df2)[grepl("pANN_.*", colnames(umap_df2))]
  umap_df2 <- umap_df2 %>% dplyr::select(-all_of(to_remove)) %>%
    replace_na(list(nCount_HTO = 0, 
                    nFeature_HTO = 0, 
                    HTO_maxID = "noId", 
                    HTO_secondID = "noId", 
                    HTO_margin = 0, 
                    HTO_classification = "noClass", 
                    HTO_classification.global = "noClass", 
                    hash.ID = "noId"))
 
  umap_df3 <- umap_df2 %>% dplyr::filter(type == "Hashed")
  
  singlets_hashIDs <- umap_df3$hash.ID
  replicate_vec <- vector("character", length = length(singlets_hashIDs))
  replicate_vec[grepl(".*REP1", singlets_hashIDs)] <- "Replicate 1"
  replicate_vec[grepl(".*REP2", singlets_hashIDs)] <- "Replicate 2"
  replicate_vec[grepl("Negative", singlets_hashIDs)] <- "Negative"
  
  Nonhashed_replicates <- umap_df2%>% dplyr::filter(type == "Nonhashed") %>% dplyr::select(Replicate) %>% flatten_chr()
  all_replicate_vec <- c(Nonhashed_replicates, replicate_vec)
  
  
  umap_df2$Replicate2 <- all_replicate_vec
  umap_df2 <- umap_df2 %>% mutate(cond_rep = paste(Replicate2, condition, sep = " "))
  umap_df2$cond_rep <- factor(umap_df2$cond_rep, levels = c("Replicate 1 NBH", "Replicate 2 NBH", "Replicate 1 RML6", "Replicate 2 RML6",
                                                            "Negative NBH", "Negative RML6"))
  
  umap_df2
}

UMAP_df <- make_UMAP_df(all_srat_Merged_filtered)
UMAP_df_ambCorr <- make_UMAP_df(all_srat_Merged_ambCorr_filtered)


# gene Symbols for selectInputs
all_genes <- rownames(GetAssayData(all_srat_Merged_filtered, slot = "data", assay = "SCT"))
all_symbols <- ensembl_symbol_table %>% dplyr::filter(ID %in% all_genes) %>% dplyr::select(Symbol) %>% flatten_chr()
all_genes_ambCorr <- rownames(GetAssayData(all_srat_Merged_ambCorr_filtered, slot = "data", assay = "SCT"))
all_symbols_ambCorr <- ensembl_symbol_table %>% dplyr::filter(ID %in% all_genes_ambCorr) %>% dplyr::select(Symbol) %>% flatten_chr()

# USER INTERFACE ==========================================================================
ui <- fluidPage(sidebarLayout(
  sidebarPanel = sidebarPanel(width = 2, 
                              radioButtons(
                                inputId = "pickSample",
                                label = "Data",
                                choices = c("Original", "Ambient corrected"),
                                selected = "Original"
                                
                              ),
                              conditionalPanel("input.pickSample == 'Original'", 
                                               selectizeInput("geneSymbol",
                                                              label = "Select Gene Symbol", 
                                                              multiple = FALSE,
                                                              size = 20,
                                                              
                                                              choices = NULL,
                                                              width = "200px",
                                                              # options = list(
                                                              # placeholder = 'eg. Agap2',
                                                              #   onInitialize = I('function() { this.setValue("");}')
                                               )),
                              conditionalPanel("input.pickSample == 'Ambient corrected'", 
                                               selectizeInput("geneSymbol_ambCorr",
                                                              label = "Select Gene Symbol", 
                                                              multiple = FALSE,
                                                              size = 20,
                                                              
                                                              choices = NULL,
                                                              width = "200px",
                                                              #options = list(
                                                              #  placeholder = 'eg. Agap2',
                                                              # onInitialize = I('function() { this.setValue("");}')
                                               )),
                              radioButtons("splitBy",
                                           "Split by:",
                                           choices = c("condition", "sample", "replicate"), 
                                           selected = "condition"),
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
                                           width = "200px")),
  mainPanel = mainPanel(width = 10, 
                        uiOutput("get_plot")
  )
  
))



# SERVER ======================================================================================================================

server <- function(input, output, session){
  
  updateSelectizeInput(session  = session,
                       inputId  = "geneSymbol",
                       server   = TRUE,
                       choices  = all_symbols,
                       selected = "Gpnmb",
                       #options = list(
                        # placeholder = 'eg. Agap2',
                         #onInitialize = I('function() { this.setValue("");}')
                       )
  
  updateSelectizeInput(session  = session,
                       inputId  = "geneSymbol_ambCorr",
                       server   = TRUE,
                       choices  = all_symbols_ambCorr,
                       selected = "Gpnmb",
                      # options = list(
                         #placeholder = 'eg. Agap2',
                        # onInitialize = I('function() { this.setValue("");}')
                       )
  
  observeEvent(input$pickSample, {
    shinyjs::reset("geneSymbol")
    shinyjs::reset("geneSymbol_ambCorr")
  })
  
  plot_height <- reactive({
    if(input$splitBy == "condition"){
      h <- "1200px"
    }else if(input$splitBy == "type"){
      h <- "1200"
    }else{
    h <- "2200px"
    }
    h
  })
  
  
  selected_gene <- reactive({
    if(input$pickSample == "Original"){
      gene <- input$geneSymbol
    }else{
      gene <- input$geneSymbol_ambCorr
    }
    gene
  })
  
  ensembl <- reactive({
    ensembl_symbol_table %>% dplyr::filter(Symbol == selected_gene()) %>% dplyr::select(ID) %>% flatten_chr()
  })
  
  seurat <- reactive({
    if(input$pickSample == "Original"){
      seu <- all_srat_Merged_filtered
    }else{
      seu <- all_srat_Merged_ambCorr_filtered
    }
    seu
  })
  
  
  original_umap_df <- reactive({
    if(input$pickSample == "Original"){
      u <- UMAP_df
    }else{
      u <- UMAP_df_ambCorr
    }
    u
  })
  
  umap_df_for_plot <- reactive({
    norm_counts <- GetAssayData(seurat(), slot = "data", assay = "SCT")[ensembl(), ]
    umap1 <- original_umap_df()
    umap1$normCounts <- norm_counts
    umap1
  })
  
  

  output$umap_plot <- renderPlot({
    table <- umap_df_for_plot()
    if(input$splitBy == "replicate"){
      table <- table %>% dplyr::filter(Replicate2 != "Negative")
    }
    
    p <- ggplot(table, aes(x = UMAP_1, y = UMAP_2, color = normCounts)) +
      geom_point(alpha = input$alpha, size = input$pointSize) +
      theme_light() +
      labs(x = "UMAP 1", y = "UMAP 2", title = paste("Selected gene:", selected_gene(), sep = " "), color = "Normalized \ncounts")
    
  if(input$splitBy == "condition"){
    p <- p + facet_grid(type~condition, scales = "free")
  }else if(input$splitBy == "sample"){
    p <- p + facet_grid(Replicate~condition, scales = "free")
  }else{
    p <- p + facet_grid(cond_rep~type, scales ="free")
  }
    
   p <-  p + theme(panel.grid.major = element_blank(), 
              panel.grid.minor = element_blank(), 
              strip.text = element_text(size = 16, color = "black"), 
              axis.text = element_text(size = 12), 
              axis.title = element_text(size = 14), 
              legend.title = element_text(size = 14), 
              legend.text = element_text(size = 14), 
              plot.title = element_text(size = 18)
              ) +
      
      scale_color_gradientn(colours = input$colors)
    
    p
    
  })



output$get_plot <- renderUI({
  
  plotOutput("umap_plot", height = plot_height())
  
})

}


shinyApp(ui = ui, server = server)

