# --- Package Management ---
list_of_packages <- c("shinyFiles", "fs", "shiny", "bslib", "Seurat", "ggplot2", "dplyr", "bsicons", "plotly", "DT")

new_packages <- list_of_packages[!(list_of_packages %in% installed.packages()[,"Package"])]
if (length(new_packages)) install.packages(new_packages)

if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")
if (!requireNamespace("presto", quietly = TRUE)) devtools::install_github("immunogenomics/presto")

invisible(lapply(c(list_of_packages, "presto"), library, character.only = TRUE))

# --- Global Configurations ---
options(shiny.maxRequestSize = 500 * 1024^2) # 500MB max file size

# Define the User Interface ####
ui <- fluidPage(
  titlePanel("Gene Expression Visualization"),
  theme = bs_theme(version = 5, bootswatch = "lumen"),
  
  # Define the sidebar
  sidebarLayout(
    sidebarPanel(
      tags$img(
        src = "https://static.wixstatic.com/media/cf7f0d_baf87b2bcb6d47a687a25cad51b2ccd5~mv2.png",
        style = "width: 100%;", 
        alt = "Developmental Genetics"
      ),
      br(), br(),
      
      # --- FILE SELECTION (shinyFiles) ---
      h4("Load Dataset"),
      shinyFilesButton(
        id = "file_picker", 
        label = "Select .rds File", 
        title = "Choose Seurat Object File", 
        multiple = FALSE,
        class = "btn-primary"
      ),
      br(), br(),
      verbatimTextOutput("loaded_file_info"),
      hr(),
      
      downloadButton("downloadDimplot", "Download Dimplot"),
      br(), br(), br(),
      
      h4("Differential Expression"),
      uiOutput("identity_test_ui"),
      uiOutput("genotype1_ui"),      
      uiOutput("genotype2_ui"),
      actionButton("run_de", "Run Differential Expression"),
      br(), br(),
      downloadButton("downloadDE", "Download DE Results"),
      downloadButton("downloadDEVolcano", "Download Volcano Plot"),
      br(), br(), br(), br(),
      
      h4("FeaturePlot & VlnPlot Inputs"),
      textInput("gene", "Enter Gene Name:"),
      actionButton("updatePlot", "Update Plot"), 
      br(), br(),
      downloadButton("downloadFeatureplot", "Download FeaturePlot"),
      downloadButton("downloadVlnplot", "Download VlnPlot"),
      br(), br(), br(),
      
      h4("Scatter Plot Inputs"),
      textInput("gene_x", "X-axis Gene:"),
      textInput("gene_y", "Y-axis Gene:"),
      textInput("gene_color", "Color by Gene:"),
      h4("Cells expressing 3 genes"),
      br(),
      tableOutput("triple_cell_count_table"),
      br(),
      downloadButton("downloadScatterplot", "Download Scatterplot"),
      br(), br(), br(),
      
      h4("Total Cells per Genotype"),
      tableOutput("total_cell_count_table"),
      h4("Cell Count by Genotype"),
      textInput("gene_count", "Enter Gene Name:"),
      tableOutput("cell_count_table"),
      width = 3
    ),
    
    # Define the cards
    mainPanel(
      fluidRow(
        column(6,
               card(
                 card_header(div(style = "text-align: center;", strong("Dimplot"))),
                 plotOutput("dimPlot"),
                 actionButton("showHeatmap", "Show Heatmap"),
                 downloadButton("downloadMarkers", "Download Markers")
               )
        ),
        column(6,
               card(
                 card_header(div(style = "text-align: center;", strong("Volcano plot"))),
                 plotOutput("de_volcano"),
                 hr(),
                 fluidRow(
                   column(4, numericInput("xlim_min", "X-axis Min:", value = -1)),
                   column(4, numericInput("xlim_max", "X-axis Max:", value = 1)),
                   column(4, numericInput("ylim_max", "Y-axis Max:", value = 5))
                 )
               )
        )
      ),
      fluidRow(
        column(6,
               card(
                 card_header(div(style = "text-align: center;", strong("FeaturePlot"))),
                 plotOutput("featurePlot"),
                 actionButton("showSplitFeaturePlot", "Split view by genotype")
               )
        ),
        column(6,
               card(
                 card_header(div(style = "text-align: center;", strong("VlnPlot"))),
                 plotOutput("vlnPlot"),
                 actionButton("showSplitVlnPlot", "Split view by genotype")
               )
        )
      ),
      fluidRow(
        column(12,
               card(
                 card_header(div(style = "text-align: center;", strong("Scatterplot"))),
                 plotOutput("scatterPlot")
               )
        )
      ),
      fluidRow(
        column(12,
               card(
                 card_header(div(style = "text-align: center;", strong("Blended FeaturePlot"))),
                 plotOutput("blendedFeaturePlot"),
                 hr(), 
                 fluidRow(
                   column(3, textInput("gene1", "Gene 1:")),
                   column(3, textInput("gene2", "Gene 2:")),
                   column(3, textInput("color1", "Color Gene 1:", value = "red")),
                   column(3, textInput("color2", "Color Gene 2:", value = "blue"))
                 ),
                 br(),
                 actionButton("updateBlendedPlot", "Update Blended Plot"),
                 actionButton("showBlendedPlot", "Split view by genotype")
               )
        )
      ),
      fluidRow(
        column(12,
               card(
                 card_header(div(style = "text-align: center;", strong("Differential Expression Results"))),
                 conditionalPanel(
                   condition = "input.run_de > 0",
                   DT::dataTableOutput("de_results")
                 )
               )
        )
      ),
      width = 9
    )
  )
)

# Define the server ####
server <- function(input, output, session) {
  
  # --- 1. FILE PICKER SETUP (shinyFiles) ---
  roots <- c(Home = fs::path_home(), getVolumes()())
  
  shinyFileChoose(
    input, 
    'file_picker', 
    roots = roots, 
    session = session,
    filetypes = c('rds', 'RData', 'rda')
  )
  
  file_path <- reactive({
    req(input$file_picker)
    parseFilePaths(roots, input$file_picker)$datapath
  })
  
  output$loaded_file_info <- renderText({
    path <- file_path()
    if (length(path) == 0 || path == "") {
      "No file selected"
    } else {
      paste("Loaded:", basename(path))
    }
  })
  
  # --- 2. REACTIVE SEURAT OBJECT LOADING & PREPROCESSING ---
  sc_obj <- reactive({
    path <- file_path()
    req(length(path) > 0 && path != "")
    
    ext <- tools::file_ext(path)
    
    withProgress(message = "Loading Seurat Object...", value = 0.3, {
      if (tolower(ext) == "rds") {
        obj <- readRDS(path)
      } else {
        env <- new.env()
        load(path, envir = env)
        obj_name <- ls(env)[1]
        obj <- env[[obj_name]]
      }
      
      incProgress(0.4, detail = "Configuring cluster levels...")
      
      obj <- SetIdent(
        obj, 
        value = factor(
          Idents(obj),
          levels = c("C1", "C2", "C3", "M1", "M2", "M3", "M4", "P1", "P2", "P3", "P4")
        )
      )
      
      incProgress(0.3, detail = "Done!")
      return(obj)
    })
  })

  genes_to_test <- reactive({
    req(sc_obj())
    obj <- sc_obj()
    rna_genes <- rownames(obj@assays$RNA)
    
    Rik.genes <- grep(pattern = "\\d+.*Rik$", x = rna_genes, value = TRUE)
    Rps.genes <- grep(pattern = "^Rps", x = rna_genes, value = TRUE)
    Rpl.genes <- grep(pattern = "^Rpl", x = rna_genes, value = TRUE)
    mito.genes <- grep(pattern = "^mt-", x = rna_genes, value = TRUE)
    pseudo.genes <- grep(pattern = "^Gm", x = rna_genes, value = TRUE)
    blood.genes <- c('Hbb-y', 'Hba-a1', 'Hbb-x', 'Hbb-bh1', 'Hba-a2', 'Hbb-bh1', "Hba-x", "Hbb-bs", "Hbb-bt")
    
    ychrom <- if (exists("removal") && "V2" %in% colnames(removal)) removal$V2[1:1555] else c()
    xchrom <- if (exists("removal") && "V3" %in% colnames(removal)) removal$V3 else c()
    
    removing <- c(Rps.genes, Rpl.genes, mito.genes, pseudo.genes, blood.genes, ychrom, xchrom, Rik.genes)
    setdiff(rna_genes, removing)
  })

  # --- 3. DYNAMIC UI RENDERERS ---
  output$identity_test_ui <- renderUI({
    req(sc_obj())
    identities <- unique(sc_obj()$identity)
    selected_val <- if ("C1" %in% identities) "C1" else identities[1]
    selectInput("identity_test", "Select Identity", choices = identities, selected = selected_val)
  })
  
  output$genotype1_ui <- renderUI({
    req(sc_obj())
    genotypes <- unique(sc_obj()$genotype)
    selectInput("genotype1", "Genotype 1", choices = genotypes, selected = genotypes[1])
  })
  
  output$genotype2_ui <- renderUI({
    req(sc_obj())
    genotypes <- unique(sc_obj()$genotype)
    selected_val <- if (length(genotypes) > 1) genotypes[2] else genotypes[1]
    selectInput("genotype2", "Genotype 2", choices = genotypes, selected = selected_val)
  })

  top10_markers <- reactiveVal(NULL)
  
  # --- 4. FEATURE PLOT AND VLN PLOT UPDATES ---
  observeEvent(input$updatePlot, {
    req(sc_obj(), input$gene)
    
    validate(
      need(input$gene %in% rownames(sc_obj()),
           paste("Gene", input$gene, "not found in Seurat object."))
    )
    
    output$featurePlot <- renderPlot({
      req(sc_obj())
      FeaturePlot(sc_obj(), features = input$gene, order = TRUE)
    })
    
    output$featurePlotSplit <- renderPlot({
      req(sc_obj())
      FeaturePlot(sc_obj(), features = input$gene, split.by = "genotype", order = TRUE)
    })
    
    output$vlnPlot <- renderPlot({
      req(sc_obj())
      VlnPlot(sc_obj(), features = input$gene, pt.size = 0.1)
    })
    
    output$vlnPlotSplit <- renderPlot({
      req(sc_obj())
      VlnPlot(sc_obj(), features = input$gene, pt.size = 0.1, split.by = 'genotype')
    })
  })
  
  # --- SHOW MODAL WITH SPLIT FEATUREPLOT ---
  observeEvent(input$showSplitFeaturePlot, {
    req(sc_obj(), input$gene)
    
    validate(
      need(input$gene %in% rownames(sc_obj()),
           paste("Gene", input$gene, "not found in Seurat Object."))
    )
    
    if (input$gene %in% rownames(sc_obj())){
      showModal(
        modalDialog(
          title = "FeaturePlot Split by Genotype",
          plotOutput("featurePlotSplitModal", height = "300px", width = "900px"),
          easyClose = TRUE,
          footer = tagList(
            modalButton("Close"),
            downloadButton("downloadSplitFeature", "Download PDF")
          )
        )
      )
    }
  })
  
  output$featurePlotSplitModal <- renderPlot({
    req(sc_obj(), input$gene)
    FeaturePlot(sc_obj(), features = input$gene, split.by = "genotype", order = TRUE)
  })
  
  output$downloadSplitFeature <- downloadHandler(
    filename = function() {
      paste("Featureplot_Split_", input$gene, "_", Sys.Date(), ".pdf", sep = "")
    },
    content = function(file) {
      req(sc_obj(), input$gene)
      plot_to_save <- FeaturePlot(sc_obj(), features = input$gene, split.by = "genotype", order = TRUE)
      pdf(file, width = 10, height = 5)
      print(plot_to_save)
      dev.off()
    }
  )
  
  # --- SHOW MODAL WITH SPLIT VLNPLOT ---
  observeEvent(input$showSplitVlnPlot, {
    req(sc_obj(), input$gene)
    
    validate(
      need(input$gene %in% rownames(sc_obj()), "Gene not found")
    )
    
    if (input$gene %in% rownames(sc_obj())){
      showModal(modalDialog(
        title = "Violin Plot - Split by Genotype",
        plotOutput("vlnPlotSplitModal", height = "300px", width = "900px"),
        easyClose = TRUE,
        footer = tagList(
          modalButton("Close"),
          downloadButton("downloadSplitVln", "Download PDF")
        )
      ))
    }
  })
  
  output$vlnPlotSplitModal <- renderPlot({
    req(sc_obj(), input$gene)
    VlnPlot(sc_obj(), features = input$gene, pt.size = 0.1, split.by = 'genotype')
  })
  
  output$downloadSplitVln <- downloadHandler(
    filename = function() {
      paste("Vlnplot_Split_", input$gene, "_", Sys.Date(), ".pdf", sep = "")
    },
    content = function(file) {
      req(sc_obj(), input$gene)
      plot_to_save <- VlnPlot(sc_obj(), features = input$gene, pt.size = 0.1, split.by = 'genotype')
      pdf(file, width = 10, height = 5)
      print(plot_to_save)
      dev.off()
    }
  )
  
# --- 5. DIFFERENTIAL EXPRESSION ---
  observeEvent(input$run_de, {
    req(sc_obj(), input$identity_test, input$genotype1, input$genotype2)
    
    withProgress(message = "Calculating Differential Expression...", value = 0, {
      # Subset cells for the selected cluster
      subset_obj <- subset(sc_obj(), idents = c(input$identity_test))
      
      incProgress(0.2, detail = "Preparing data...")
      
      tryCatch({
        # Set genotype as active ident for the subset
        Idents(subset_obj) <- subset_obj$genotype
        
        incProgress(0.4, detail = "Running Wilcoxon test...")
        
        # Native Seurat v5 DE test
        de_res <- FindMarkers(
          subset_obj,
          ident.1 = input$genotype2,
          ident.2 = input$genotype1,
          test.use = "wilcox",
          logfc.threshold = 0,
          min.pct = 0.1
        )
        
        # Format columns for volcano plot and table
        de_wilcox <- de_res %>%
          tibble::rownames_to_column(var = "feature") %>%
          rename(logFC = avg_log2FC, padj = p_val_adj) %>%
          mutate(
            padj = ifelse(padj == 0, .Machine$double.xmin, padj),
            DE = abs(logFC) > log2(1.1) & padj < 0.01,
            DEG = ifelse(DE, feature, NA)
          )
        
        incProgress(0.8, detail = "Rendering results...")
        
        initial_xlim <- c(min(de_wilcox$logFC, na.rm = TRUE), max(de_wilcox$logFC, na.rm = TRUE))
        initial_ylim <- c(0, max(-log10(de_wilcox$padj), na.rm = TRUE))
        
        updateNumericInput(session, "xlim_min", value = round(initial_xlim[1], 2))
        updateNumericInput(session, "xlim_max", value = round(initial_xlim[2], 2))
        updateNumericInput(session, "ylim_max", value = round(initial_ylim[2], 2))
        
        output$de_results <- DT::renderDataTable({
          de_wilcox %>% arrange(padj)
        })
        
        output$de_volcano <- renderPlot({
          ggplot(de_wilcox, aes(x = logFC, y = -log10(padj), col = DE, label = DEG)) +
            geom_point() +
            ggrepel::geom_text_repel() +
            geom_vline(xintercept = c(-log2(1.1), log2(1.1), 0), linetype = "dotted") +
            geom_hline(yintercept = -log10(0.01), linetype = "dotted") +
            scale_color_manual(values = c("#909090", "red")) +
            theme_minimal() +
            xlim(input$xlim_min, input$xlim_max) +
            ylim(0, input$ylim_max)
        })
        incProgress(1.0, detail = "Done!")
        
      }, error = function(e) {
        output$de_volcano <- renderPlot({
          ggplot() +
            geom_text(aes(x = 0.5, y = 0.5, label = paste("Error in DE:", e$message)), size = 5) +
            theme_void()
        })
        output$de_results <- DT::renderDataTable(NULL)
      })
    })
  })
  
  # --- DE DOWNLOAD HANDLERS ---
  output$downloadDEVolcano <- downloadHandler(
    filename = function() {
      paste("VolcanoPlot_", input$identity_test, "_", input$genotype1, "_vs_", input$genotype2, "_", Sys.Date(), ".pdf", sep = "")
    },
    content = function(file) {
      req(sc_obj(), input$identity_test, input$genotype1, input$genotype2)
      subset_obj <- subset(sc_obj(), idents = c(input$identity_test))
      Idents(subset_obj) <- subset_obj$genotype
      
      de_res <- FindMarkers(
        subset_obj,
        ident.1 = input$genotype2,
        ident.2 = input$genotype1,
        test.use = "wilcox",
        logfc.threshold = 0,
        min.pct = 0.1
      )
      
      de_wilcox <- de_res %>%
        tibble::rownames_to_column(var = "feature") %>%
        rename(logFC = avg_log2FC, padj = p_val_adj) %>%
        mutate(
          padj = ifelse(padj == 0, .Machine$double.xmin, padj),
          DE = abs(logFC) > log2(1.1) & padj < 0.01,
          DEG = ifelse(DE, feature, NA)
        )
        
      pdf(file, width = 8, height = 6)
      print(ggplot(de_wilcox, aes(x = logFC, y = -log10(padj), col = DE, label = DEG)) +
              geom_point() +
              ggrepel::geom_text_repel() +
              geom_vline(xintercept = c(-log2(1.1), log2(1.1), 0), col = "#303030", linetype = "dotted") +
              geom_hline(yintercept = -log10(0.01), col = "#303030", linetype = "dotted") +
              scale_color_manual(values = c("#909090", "red")) +
              theme_minimal() +
              xlim(input$xlim_min, input$xlim_max) +
              ylim(0, input$ylim_max)
      )
      dev.off()
    }
  )
  
  output$downloadDE <- downloadHandler(
    filename = function() {
      paste("DE_results_", input$identity_test, "_", input$genotype1, "_vs_", input$genotype2, "_", Sys.Date(), ".csv", sep = "")
    },
    content = function(file) {
      req(sc_obj(), input$identity_test, input$genotype1, input$genotype2)
      subset_obj <- subset(sc_obj(), idents = c(input$identity_test))
      Idents(subset_obj) <- subset_obj$genotype
      
      de_res <- FindMarkers(
        subset_obj,
        ident.1 = input$genotype2,
        ident.2 = input$genotype1,
        test.use = "wilcox",
        logfc.threshold = 0,
        min.pct = 0.1
      )
      
      de_wilcox <- de_res %>%
        tibble::rownames_to_column(var = "feature") %>%
        rename(logFC = avg_log2FC, padj = p_val_adj) %>%
        mutate(
          padj = ifelse(padj == 0, .Machine$double.xmin, padj),
          DE = abs(logFC) > log2(1.1) & padj < 0.01,
          DEG = ifelse(DE, feature, NA)
        )
      write.csv(de_wilcox, file, row.names = FALSE)
    }
  )
  # --- 6. SCATTERPLOT & DIMPLOT DOWNLOADS ---
  output$downloadScatterplot <- downloadHandler(
    filename = function() {
      paste("Scatterplot_", Sys.Date(), "_", input$gene_x, "_", input$gene_y, "_", input$gene_color, ".pdf", sep = "")
    },
    content = function(file) {
      req(sc_obj(), input$gene_x, input$gene_y, input$gene_color)
      validate(
        need(all(c(input$gene_x, input$gene_y, input$gene_color) %in% rownames(sc_obj())),
             "One or more selected genes not found in Seurat object.")
      )
      pdf(file, width = 7.5, height = 4.5)
      expr_data <- FetchData(sc_obj(), vars = c(input$gene_x, input$gene_y, input$gene_color, "genotype"))
      print(ggplot(expr_data, aes_string(x = input$gene_x, y = input$gene_y)) +
              geom_point(aes_string(color = input$gene_color), size = 2, alpha = 0.8) +
              scale_color_gradient(low = "grey", high = "red") +
              labs(color = paste(input$gene_color, "Expression")) +
              geom_vline(xintercept = 0, linetype = "dashed") +
              geom_hline(yintercept = 0, linetype = "dashed") +
              theme_minimal() +
              facet_wrap(~ genotype))
      dev.off()
    }
  )
  
  output$downloadDimplot <- downloadHandler(
    filename = function() {
      paste("Dimplot_", Sys.Date(), ".pdf", sep = "")
    },
    content = function(file) {
      req(sc_obj())
      pdf(file, width = 7.5, height = 4.5)
      print(DimPlot(sc_obj(), group.by = 'identity', reduction = "umap", label = TRUE, label.size = 6))
      dev.off()
    }
  )
  
  output$downloadFeatureplot <- downloadHandler(
    filename = function() {
      paste("Featureplot_", input$gene, "_", Sys.Date(), ".pdf", sep = "")
    },
    content = function(file) {
      req(sc_obj(), input$gene)
      if (input$gene %in% rownames(sc_obj())) {
        pdf(file, width = 7.5, height = 4.5)
        print(FeaturePlot(sc_obj(), features = input$gene, order = TRUE))
        dev.off()
      }
    }
  )
  
  output$downloadVlnplot <- downloadHandler(
    filename = function() {
      paste("Vlnplot_", input$gene, "_", Sys.Date(), ".pdf", sep = "")
    },
    content = function(file) {
      req(sc_obj(), input$gene)
      if (input$gene %in% rownames(sc_obj())) {
        pdf(file, width = 7.5, height = 4.5)
        print(VlnPlot(sc_obj(), features = input$gene, pt.size = 0.1))
        dev.off()
      }
    }
  )
  
  # --- 7. TABLES ---
  output$cell_count_table <- renderTable({
    req(sc_obj(), input$gene_count)
    if (input$gene_count %in% rownames(sc_obj())) {
      gene_expr <- FetchData(sc_obj(), vars = c(input$gene_count, "genotype"))
      gene_sym <- input$gene_count
      cell_count <- table(subset(gene_expr, gene_expr[[gene_sym]] > 0)$genotype)
      t(as.matrix(cell_count))
    }
  }, rownames = TRUE, colnames = TRUE)
  
  output$total_cell_count_table <- renderTable({
    req(sc_obj())
    total_cell_count <- table(sc_obj()$genotype)
    t(as.matrix(total_cell_count))
  }, rownames = TRUE, colnames = TRUE)
  
  output$triple_cell_count_table <- renderTable({
    req(sc_obj(), input$gene_x, input$gene_y, input$gene_color)
    if (all(c(input$gene_x, input$gene_y, input$gene_color) %in% rownames(sc_obj()))) {
      gene_expr <- FetchData(sc_obj(), vars = c(input$gene_x, input$gene_y, input$gene_color, "genotype"))
      combined_expression <- gene_expr[[input$gene_x]] > 0 & 
        gene_expr[[input$gene_y]] > 0 & 
        gene_expr[[input$gene_color]] > 0
      
      cell_count <- table(subset(gene_expr, combined_expression)$genotype)
      t(as.matrix(cell_count))
    }
  }, rownames = TRUE, colnames = TRUE)
  
  # --- 8. INITIAL DIMPLOT & HEATMAP ---
  output$dimPlot <- renderPlot({
    req(sc_obj())
    DimPlot(sc_obj(), group.by = 'identity', reduction = "umap", label = TRUE, label.size = 6)
  })
  
  observeEvent(input$showHeatmap, {
    req(sc_obj())
    
    showModal(modalDialog(
      title = "Heatmap of Top 10 Markers per Cluster",
      "Calculating markers and rendering heatmap...",
      plotOutput("heatmap", height = "700px"),
      easyClose = TRUE,
      footer = tagList(
        modalButton("Close"),
        downloadButton("downloadHeatmap", "Download PDF")
      )
    ))
    
    withProgress(message = "Calculating Markers...", value = 0, {
      cluster_markers_identity <- FindAllMarkers(
        sc_obj(), 
        features = genes_to_test(),
        only.pos = TRUE, 
        min.pct = 0.25,
        logfc.threshold = log(1.2)
      )
      incProgress(0.6, detail = "Grouping Markers...")
      
      top10 <- cluster_markers_identity %>%
        group_by(cluster) %>%
        top_n(n = 10, wt = avg_log2FC)
        
      top10_markers(top10)
      incProgress(1.0, detail = "Done")
    })
    
    output$heatmap <- renderPlot({
      req(sc_obj(), top10_markers())
      DoHeatmap(sc_obj(), features = top10_markers()$gene) +
        theme(axis.text.y = element_text(size = 8))
    })
  })
  
  output$downloadHeatmap <- downloadHandler(
    filename = function() {
      paste("Heatmap_", Sys.Date(), ".pdf", sep = "")
    },
    content = function(file) {
      req(sc_obj(), top10_markers())
      pdf(file, width = 6, height = 8)
      print(DoHeatmap(sc_obj(), features = top10_markers()$gene) +
              theme(axis.text.y = element_text(size = 8)))
      dev.off()
    }
  )
  
  output$downloadMarkers <- downloadHandler(
    filename = function() {
      paste("Cluster_Markers_", Sys.Date(), ".csv", sep = "")
    },
    content = function(file) {
      req(top10_markers())
      write.csv(top10_markers(), file, row.names = FALSE)
    }
  )
  
  output$scatterPlot <- renderPlot({
    req(sc_obj(), input$gene_x, input$gene_y, input$gene_color)
    validate(
      need(all(c(input$gene_x, input$gene_y, input$gene_color) %in% rownames(sc_obj())),
           "One or more selected genes not found in Seurat object.")
    )
    if (all(c(input$gene_x, input$gene_y, input$gene_color) %in% rownames(sc_obj()))) {
      expr_data <- FetchData(sc_obj(), vars = c(input$gene_x, input$gene_y, input$gene_color, "genotype"))
      ggplot(expr_data, aes_string(x = input$gene_x, y = input$gene_y)) +
        geom_point(aes_string(color = input$gene_color), size = 2, alpha = 0.8) +
        scale_color_gradient(low = "grey", high = "red") +
        labs(color = paste(input$gene_color, "Expression")) +
        geom_vline(xintercept = 0, linetype = "dashed") +
        geom_hline(yintercept = 0, linetype = "dashed") +
        theme_minimal() +
        facet_wrap(~ genotype)
    }
  })

  # --- 9. BLENDED FEATUREPLOT ---
  observeEvent(input$updateBlendedPlot, {
    req(sc_obj(), input$gene1, input$gene2, input$color1, input$color2)
    
    validate(
      need(input$gene1 %in% rownames(sc_obj()) && input$gene2 %in% rownames(sc_obj()),
           "One or both of the specified genes were not found.")
    )
    
    output$blendedFeaturePlot <- renderPlot({
      default_color1 <- "red"
      default_color2 <- "blue"
      default_color3 <- "grey90"
      
      color1_val <- ifelse(is.null(input$color1) || input$color1 == "", default_color1, input$color1)
      color2_val <- ifelse(is.null(input$color2) || input$color2 == "", default_color2, input$color2)
      color3_val <- default_color3
      
      FeaturePlot(
        sc_obj(), 
        features = c(input$gene1, input$gene2), 
        blend = TRUE, 
        order = TRUE, 
        blend.threshold = 0.1, 
        cols = c(color3_val, color1_val, color2_val)
      ) 
    })
  })
  
  observeEvent(input$showBlendedPlot, {
    req(sc_obj(), input$gene1, input$gene2, input$color1, input$color2)
    
    validate(
      need(input$gene1 %in% rownames(sc_obj()) && input$gene2 %in% rownames(sc_obj()),
           "One or both of the specified genes were not found in Seurat Object.")
    )
    
    if (input$gene1 %in% rownames(sc_obj()) && input$gene2 %in% rownames(sc_obj())){
      default_color1 <- "red"
      default_color2 <- "blue"
      default_color3 <- "grey90"
      
      color1_val <- ifelse(is.null(input$color1) || input$color1 == "", default_color1, input$color1)
      color2_val <- ifelse(is.null(input$color2) || input$color2 == "", default_color2, input$color2)
      color3_val <- default_color3
      
      output$blendedPlotSplitModal <- renderPlot({
        FeaturePlot(
          sc_obj(), 
          features = c(input$gene1, input$gene2), 
          blend = TRUE,
          order = TRUE, 
          blend.threshold = 0.1,
          split.by = 'genotype',
          cols = c(color3_val, color1_val, color2_val)
        )
      })
      
      showModal(
        modalDialog(
          title = "BlendedPlot Split by Genotype",
          plotOutput("blendedPlotSplitModal", height = "600px", width = "900px"),
          easyClose = TRUE,
          footer = tagList(
            modalButton("Close"),
            downloadButton("downloadSplitPlot", "Download PDF")
          )
        )
      )
    }
  })
  
  output$downloadSplitPlot <- downloadHandler(
    filename = function() {
      req(input$gene1, input$gene2)
      if (!is.null(input$gene1) && !is.null(input$gene2)) {
        paste("BlendedPlot_Split_", input$gene1, "_", input$gene2, ".pdf", sep = "")
      } else {
        "BlendedPlot_Split.pdf"
      }
    },
    content = function(file) {
      req(sc_obj(), input$gene1, input$gene2, input$color1, input$color2)
      
      default_color1 <- "red"
      default_color2 <- "blue"
      default_color3 <- "grey90"
      
      color1_val <- ifelse(is.null(input$color1) || input$color1 == "", default_color1, input$color1)
      color2_val <- ifelse(is.null(input$color2) || input$color2 == "", default_color2, input$color2)
      color3_val <- default_color3
      
      plot_to_save <- FeaturePlot(
        sc_obj(), 
        features = c(input$gene1, input$gene2), 
        blend = TRUE,
        order = TRUE, 
        blend.threshold = 0.1,
        split.by = 'genotype',
        cols = c(color3_val, color1_val, color2_val)
      )
      
      pdf(file, width = 12, height = 8)
      print(plot_to_save)
      dev.off()
    }
  )
}

shinyApp(ui = ui, server = server)
