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
all_srat_Merged_ambCorr_filtered <- readRDS("/Users/marusa/Documents/snRNAseq_PRDs_neurons/all_srat_Merged_ambCorr_filtered2.rds")

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
ui <- fluidPage(tabsetPanel(
  tabPanel("Gene search", 
  sidebarLayout(
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
                                           choices = c("condition", "condition and type",  "sample", "replicate"), 
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
  ))), 
  tabPanel("Explore Metadata", 
           sidebarLayout(
             sidebarPanel(width = 2, 
                          radioButtons(
                            inputId = "pickSampleMeta",
                            label = "Data",
                            choices = c("Original", "Ambient corrected"),
                            selected = "Original"
                            
                          ), 
                          selectInput("metadata", 
                                      label = "Select metadata column", 
                                      choices = c(
                                                  "Clusters (resolution 0.1)", 
                                                  "Clusters (resolution 0.2)", 
                                                  "Clusters (resolution 0.3)", 
                                                  "Clusters (resolution 0.4)", 
                                                  "Clusters (resolution 0.5)", 
                                                  "Ribosomal protein genes content", 
                                                  "Mouse RNAseq reference classification", 
                                                  "Chen et al. classification", 
                                                  "Campbell et al. classification"
                                                  ) , 
                                      selected = "Clusters (resolution 0.1)",
                                      multiple = FALSE, 
                                      width = "200px"
                                
                                      ),
                          uiOutput("get_meta_cat"), 
                       
                          radioButtons("splitBymeta",
                                      "Split by:",
                                      choices = c("none", "condition", "condition and type", "sample", "replicate"),
                                      selected = "condition"), 
                          numericInput("pointSizeMeta", 
                                       label = "Select point size", 
                                       value = 0.1,
                                       min = 0.1, 
                                       max = 10, 
                                       step = 0.1, 
                                       width = "200px"),
                          numericInput("alphaMeta",
                                       label = "Select point transparency", 
                                       value = 1, 
                                       min = 0,
                                       max = 1,
                                       step = 0.1, 
                                       width = "200px")
                          ), 
             mainPanel(uiOutput("get_meta_plot"))
           ))
  
))



# SERVER ======================================================================================================================

server <- function(input, output, session){
  
########################## GENE SEARCH TAB ##############################################
  
  updateSelectizeInput(session  = session,
                       inputId  = "geneSymbol",
                       server   = TRUE,
                       choices  = all_symbols,
                       selected = character(0),
                       options = list(
                        placeholder = 'eg. Agap2',
                        onInitialize = I('function() { this.setValue("");}')
                       ))
  
  updateSelectizeInput(session  = session,
                       inputId  = "geneSymbol_ambCorr",
                       server   = TRUE,
                       choices  = all_symbols_ambCorr,
                       selected = character(0),
                      options = list(
                        placeholder = 'eg. Agap2',
                      onInitialize = I('function() { this.setValue("");}')
                       ))
  
  observeEvent(input$pickSample, {
    shinyjs::reset("geneSymbol")
    shinyjs::reset("geneSymbol_ambCorr")
  })
  
  plot_height <- reactive({
    
    if(selected_gene() == ""){
      h <- "800px"
    }else{
    if(input$splitBy == "condition"){
      h <- "800px"
    }else if(input$splitBy == "condition and type"){
      h <- "1200px"
    }else if(input$splitBy == "sample"){
      h <- "1700px"
    }else{
    h <- "2400px"
    }
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
    
    if(selected_gene() == ""){
     umap1 <- data.frame(matrix(ncol = 0, nrow = 0))
     
    }else{

    norm_counts <- GetAssayData(seurat(), slot = "data", assay = "SCT")[ensembl(), ]
    umap1 <- original_umap_df()
    umap1$normCounts <- norm_counts
    }
    
    umap1
  })
  
  

  output$umap_plot <- renderPlot({
    table <- umap_df_for_plot()
    
    if(selected_gene() == ""){ # if no genes are selected print an empty umap
      p <- ggplot(data.frame(x = c(1,2,3), y = c(1,2,3)), aes(x = x, y = y)) +
        theme_light() +
        labs(x = "UMAP 1", y = "UMAP 2") +
            theme(axis.title = element_text(size = 14),
              panel.grid.minor = element_blank(), 
              panel.grid.major = element_blank())
    
      
    }else{
      
    if(input$splitBy == "replicate"){
      table <- table %>% dplyr::filter(Replicate2 != "Negative")
    }
    
    p <- ggplot(table, aes(x = UMAP_1, y = UMAP_2, color = normCounts)) +
      geom_point(alpha = input$alpha, size = input$pointSize) +
      theme_light() +
      labs(x = "UMAP 1", y = "UMAP 2", title = paste("Selected gene:", selected_gene(), sep = " "), color = "Normalized \ncounts")
    
  if(input$splitBy == "condition"){
    p <- p + facet_wrap(~condition, scales = "free")
  }else if(input$splitBy == "condition and type"){
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
    
    
    
    }
    
    return(p)
    
  })



output$get_plot <- renderUI({
  
  plotOutput("umap_plot", height = plot_height())
  
})


####################### EXPLORE METADATA TAB ##################################################

# pick the correct metadata column based on metadata select input
meta_column <- reactive({
  if(input$metadata == "Clusters (resolution 0.1)"){
    col <- "SCT_snn_res.0.1"
  }else if(input$metadata == "Clusters (resolution 0.2)"){
    col <-  "SCT_snn_res.0.2"
  }else if(input$metadata ==  "Clusters (resolution 0.3)"){
    col <- "SCT_snn_res.0.3"
  }else if(input$metadata == "Clusters (resolution 0.4)"){
    col <- "SCT_snn_res.0.4" 
  }else if(input$metadata ==  "Clusters (resolution 0.5)"){
    col <-  "SCT_snn_res.0.5"
  }else if(input$metadata == "Ribosomal protein genes content"){
    col <- "high_subsets_Ribo_percent"
  }else if(input$metadata == "Mouse RNAseq reference classification"){
    col <- "mouseRNAseq_labels" 
  }else if(input$metadata == "Chen et al. classification"){
    col <- "Chen_labels"
  }else if(input$metadata == "Campbell et al. classification"){
    col <- "Campbell_labels"
  }

  col
})

original_umap_df_meta <- reactive({
  if(input$pickSampleMeta == "Original"){
    u <- UMAP_df
  }else{
    u <- UMAP_df_ambCorr
  }
  u
})

# select the categories of metadata columns to be displayed in selectCat input
meta_categories <- reactive({
  umap_df <- original_umap_df_meta()
  meta_cats  <- umap_df %>% dplyr::select(meta_column()) 
  meta_cats <- as.character(sort(unique(meta_cats[, 1])))
  meta_cats <- c("all", meta_cats)
  meta_cats
})

output$get_meta_cat <- renderUI({
  
  selectInput("metaCat",
              label = "Colored metadata category", 
              multiple = FALSE,
              #size = 20,
              choices = meta_categories(),
              width = "200px"
             
  )
  
})

# add cat_color column to umap_df based on selected metadata categories 

umap_for_metaTab <- reactive({
  
  # get selected metadata column as vector
  umap_df <- original_umap_df_meta()
  meta_cats  <- umap_df %>% dplyr::select(meta_column()) 
  meta_cats <- as.character(meta_cats[, 1])
  
  # if metaCat is "all" keep all metata categories the same, if not keep only selected category as it is and assign others to the same category
  
  if(input$metaCat == "all"){
    meta_color <- meta_cats
  }else{
    cat_to_kepp <- input$metaCat
    meta_cats[which(meta_cats != cat_to_kepp)] <- "Others"
    meta_color <- meta_cats
    
  }
  
  umap_df$metaColor <- meta_color
  umap_df
  
})


plot_height_meta <- reactive({
  
  if(input$splitBymeta == "condition"|input$splitBymeta == "none"){
    h <- "750px"
  }else if(input$splitBymeta == "condition and type"){
    h <- "1200px"
  }else if(input$splitBymeta == "sample"){
    h <- "1650px"
  }else{
    h <- "2400px"
  }
  
  h
})



output$umap_metadata_plot <- renderPlot({
  
  table <- umap_for_metaTab()
  
  if(input$splitBymeta == "replicate"){
    table <- table %>% dplyr::filter(Replicate2 != "Negative")
  }
  
  p <- ggplot(table, aes(x = UMAP_1, y = UMAP_2, color = metaColor)) +
    geom_point(alpha = input$alphaMeta, size = input$pointSizeMeta) +
    theme_light() +
    labs(x = "UMAP 1", y = "UMAP 2", title = paste("Selected metadata category:", input$metaCat, sep = " "), color = input$metadata)
  
  if(input$splitBymeta == "condition"){
    p <- p + facet_wrap(~condition, scales = "free")
  }else if(input$splitBymeta == "condition and type"){
    p <- p + facet_grid(type~condition, scales = "free")
  }else if(input$splitBymeta == "sample"){
    p <- p + facet_grid(Replicate~condition, scales = "free")
  }else if(input$splitBymeta == "replicate"){
    p <- p + facet_grid(cond_rep~type, scales ="free")
  }else{
    p <- p
  }
  
  
  
  p <-  p + theme(panel.grid.major = element_blank(), 
                  panel.grid.minor = element_blank(), 
                  strip.text = element_text(size = 16, color = "black"), 
                  axis.text = element_text(size = 12), 
                  axis.title = element_text(size = 14), 
                  legend.title = element_text(size = 14), 
                  legend.text = element_text(size = 14), 
                  plot.title = element_text(size = 18))+
                    
       guides(colour = guide_legend(override.aes = list(size= 5)))
  
  if(input$metaCat != "all"){
    p <- p + scale_color_manual(values = c("#E5E7E9", "#E74C3C"), breaks = c("Others", input$metaCat))
  }
  
  
 p 
  
})


output$get_meta_plot <- renderUI({
  
  if(input$splitBymeta == "none"){
    w <- "950px"
  }else{
    w <- "1350px"
  }
  
  plotOutput("umap_metadata_plot", height = plot_height_meta(), width = w)
  
})


}


shinyApp(ui = ui, server = server)

