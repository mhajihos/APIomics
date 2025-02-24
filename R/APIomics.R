#Main-Shiny_App_Code

APIomics<-function()
{
#Packages

  suppressMessages(suppressWarnings(library(shiny)))
  suppressMessages(suppressWarnings(library(shinydashboard)))
  suppressMessages(suppressWarnings(library(shinyjs)))
  suppressMessages(suppressWarnings(library(limma)))
  suppressMessages(suppressWarnings(library(edgeR)))
  suppressMessages(suppressWarnings(library(clusterProfiler)))
  suppressMessages(suppressWarnings(library(org.Hs.eg.db)))
  suppressMessages(suppressWarnings(library(DT)))
  suppressMessages(suppressWarnings(library(shinyalert)))
  suppressMessages(suppressWarnings(library(gplots)))
  suppressMessages(suppressWarnings(library(RColorBrewer)))
  suppressMessages(suppressWarnings(library(plotly)))
  suppressMessages(suppressWarnings(library(ggplot2)))
  suppressMessages(suppressWarnings(library(ggheatmap)))
  suppressMessages(suppressWarnings(library(WGCNA)))
  suppressMessages(suppressWarnings(library(heatmaply)))
  suppressMessages(suppressWarnings(library(tidyverse)))
  suppressMessages(suppressWarnings(library(reshape2)))
  suppressMessages(suppressWarnings(library(ggrepel)))
  suppressMessages(suppressWarnings(library(ggraph)))
  suppressMessages(suppressWarnings(library(igraph)))
  suppressMessages(suppressWarnings(library(tidygraph)))
  suppressMessages(suppressWarnings(library(ComplexHeatmap)))
  suppressMessages(suppressWarnings(library(shinydashboardPlus)))
  suppressMessages(suppressWarnings(library(corto)))
  suppressMessages(suppressWarnings(library(dplyr)))
  suppressMessages(suppressWarnings(library(httr)))
  suppressMessages(suppressWarnings(library(jsonlite)))
  suppressMessages(suppressWarnings(library(rentrez)))
  suppressMessages(suppressWarnings(library(DOSE)))
  suppressMessages(suppressWarnings(library(visNetwork)))
  suppressMessages(suppressWarnings(library(xml2)))
allowWGCNAThreads() 


# Database search functions
search_chembl <- function(genes) {
  results <- data.frame(
    molecule_pref_name = character(),
    molecule_chembl_id = character(),
    target_pref_name = character(),
    target_chembl_id = character(),
    target_genes = character(),  # ✅ NEW: List of target genes for each compound
    standard_type = character(),
    standard_value = numeric(),
    standard_units = character(),
    max_phase = numeric(),
    indication_class = character(),
    pubchem_cid = character(),
    chembl_img_url = character(),
    stringsAsFactors = FALSE
  )
  
  total_genes <- length(genes)
  i <- 0
  
  for (gene in genes) {
    i <- i + 1
    incProgress((1 / total_genes), detail = paste("Processing", gene, "..."))
    
    tryCatch({
      url <- paste0("https://www.ebi.ac.uk/chembl/api/data/target/search.xml?q=", URLencode(gene))
      response <- GET(url)
      
      if (status_code(response) == 200) {
        xml_content <- read_xml(rawToChar(response$content))
        target_ids <- xml_find_all(xml_content, "//target_chembl_id") %>% xml_text()
        
        for (target_id in target_ids) {
          compounds_url <- paste0("https://www.ebi.ac.uk/chembl/api/data/activity.xml?target_chembl_id=", target_id)
          compounds_response <- GET(compounds_url)
          
          if (status_code(compounds_response) == 200) {
            compounds_xml <- read_xml(rawToChar(compounds_response$content))
            activities <- xml_find_all(compounds_xml, "//activity")
            
            for (activity in activities) {
              molecule_id <- xml_text(xml_find_first(activity, ".//molecule_chembl_id")) %||% NA
              molecule_name <- xml_text(xml_find_first(activity, ".//molecule_pref_name")) %||% NA
              
              # **1️⃣ Get Target Genes for This Molecule**
              target_gene_url <- paste0("https://www.ebi.ac.uk/chembl/api/data/target/", target_id, ".xml")
              target_gene_response <- GET(target_gene_url)
              
              target_genes <- "No target genes found"  # Default
              if (status_code(target_gene_response) == 200) {
                target_gene_xml <- read_xml(rawToChar(target_gene_response$content))
                gene_list <- xml_find_all(target_gene_xml, "//gene_symbol") %>% xml_text()
                
                if (length(gene_list) > 0) {
                  target_genes <- paste(gene_list, collapse = ", ")
                }
              }
              
              # **2️⃣ PubChem CID lookup**
              pubchem_url <- paste0("https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/", URLencode(molecule_name), "/cids/TXT")
              pubchem_response <- GET(pubchem_url)
              pubchem_cid <- if (status_code(pubchem_response) == 200) rawToChar(pubchem_response$content) else NA
              
              # **3️⃣ ChEMBL structure image URL**
              chembl_img_url <- paste0("https://www.ebi.ac.uk/chembl/api/data/image/", molecule_id, ".svg")
              
              temp_df <- data.frame(
                molecule_pref_name = molecule_name,
                molecule_chembl_id = molecule_id,
                target_pref_name = xml_text(xml_find_first(activity, ".//target_pref_name")) %||% NA,
                target_chembl_id = target_id,
                target_genes = target_genes,  # ✅ NEW: List of genes related to the target
                standard_type = xml_text(xml_find_first(activity, ".//standard_type")) %||% NA,
                standard_value = as.numeric(xml_text(xml_find_first(activity, ".//standard_value")) %||% NA),
                standard_units = xml_text(xml_find_first(activity, ".//standard_units")) %||% NA,
                max_phase = as.numeric(xml_text(xml_find_first(activity, ".//max_phase")) %||% NA),
                indication_class = xml_text(xml_find_first(activity, ".//indication_class")) %||% NA,
                pubchem_cid = trimws(pubchem_cid),
                chembl_img_url = chembl_img_url,
                stringsAsFactors = FALSE
              )
              results <- rbind(results, temp_df)
            }
          }
        }
      }
    }, error = function(e) {
      warning(paste("Error processing gene:", gene, "-", e$message))
    })
    
    Sys.sleep(0.2)  # Small delay to prevent API overload
  }
  return(results)
}



search_bindingdb <- function(genes) {
  results <- data.frame(
    compound_name = character(),
    target_name = character(),
    target_organism = character(),
    affinity_type = character(),
    affinity_value = numeric(),
    affinity_units = character(),
    assay_type = character(),
    bindingdb_id = character(),
    stringsAsFactors = FALSE
  )
  
  total_genes <- length(genes)
  i <- 0
  
  for (gene in genes) {
    i <- i + 1
    incProgress(1 / total_genes, detail = paste("Processing", gene, "..."))
    
    # **🔍 Debugging Step: Print API URL**
    url <- paste0("https://www.bindingdb.org/REST/ligand?target=", URLencode(gene))
    print(paste("Querying BindingDB:", url))
    
    response <- GET(url)
    
    if (status_code(response) == 200) {
      xml_content <- read_xml(rawToChar(response$content))
      
      # **🔍 Debugging Step: Check if XML has expected structure**
      print(xml_content)
      
      compounds <- xml_find_all(xml_content, "//ligand")
      
      if (length(compounds) == 0) {
        warning(paste("No results found for gene:", gene))
      }
      
      for (compound in compounds) {
        temp_df <- data.frame(
          compound_name = xml_text(xml_find_first(compound, ".//ligandName")) %||% NA,
          target_name = xml_text(xml_find_first(compound, ".//targetName")) %||% NA,
          target_organism = xml_text(xml_find_first(compound, ".//targetOrganism")) %||% NA,
          affinity_type = xml_text(xml_find_first(compound, ".//affinityType")) %||% NA,
          affinity_value = as.numeric(xml_text(xml_find_first(compound, ".//affinityValue")) %||% NA),
          affinity_units = xml_text(xml_find_first(compound, ".//affinityUnits")) %||% NA,
          assay_type = xml_text(xml_find_first(compound, ".//assayType")) %||% NA,
          bindingdb_id = xml_text(xml_find_first(compound, ".//bindingDBid")) %||% NA,
          stringsAsFactors = FALSE
        )
        results <- rbind(results, temp_df)
      }
    } else {
      warning(paste("BindingDB API request failed for:", gene, "Status:", status_code(response)))
    }
    
    Sys.sleep(0.2)  # Small delay to prevent API overload
  }
  
  if (nrow(results) == 0) {
    print("❌ No BindingDB results found for any gene.")
  } else {
    print("✅ Successfully retrieved BindingDB data.")
  }
  
  return(results)
}



# Helper function to safely handle NULL values
`%||%` <- function(x, y) if (is.null(x)) y else x


addResourcePath("static", system.file("www", package = "APIomics", mustWork = TRUE))

# 1000 MB (1 GB) file size limit
options(shiny.maxRequestSize = 1000 * 1024^2) 
ui <- dashboardPage(
  skin = "black",
  dashboardHeader(
    title = div(
      HTML('<svg width="40" height="40" viewBox="0 0 100 100" xmlns="http://www.w3.org/2000/svg">
                <circle cx="50" cy="50" r="45" fill="#007bff"/>
                <text x="50" y="60" text-anchor="middle" fill="white" font-size="30" font-weight="bold">API</text>
            </svg>'),
      span("APIomics v1.0", style = "margin-left: 10px; font-weight: 600; color: #333;"),
      style = "display: flex; align-items: center; font-size: 22px;"
    ),
    titleWidth = 250,
    tags$li(
      class = "dropdown",
      style = "padding: 10px;",
      actionButton("theme_toggle", label = NULL, icon = icon("adjust"), 
                   style = "background: none; border: none; color: #007bff;")
    )
  ),
  dashboardSidebar(
    width = 250,
    tags$head(
      tags$style(HTML("
                :root {
                    --primary-color: #007bff;
                    --secondary-color: #6c757d;
                    --background-light: #f8f9fa;
                    --background-dark: #343a40;
                }
                body {
                    font-family: 'Inter', 'Roboto', sans-serif;
                    transition: background-color 0.3s, color 0.3s;
                }
                .skin-black .sidebar {
                    background-color: var(--background-light);
                    border-right: 1px solid #e9ecef;
                }
                .sidebar-menu li a {
                    font-size: 15px;
                    font-weight: 500;
                    color: var(--secondary-color);
                    border-radius: 8px;
                    margin: 5px;
                    transition: all 0.3s ease;
                }
                .sidebar-menu li a:hover {
                    background-color: rgba(0,123,255,0.1);
                    color: var(--primary-color);
                }
                .sidebar-menu .active a {
                    background-color: var(--primary-color);
                    color: white !important;
                }
                .main-header {
                    box-shadow: 0 2px 4px rgba(0,0,0,0.1);
                }
                .box {
                    border-radius: 12px;
                    box-shadow: 0 4px 6px rgba(0,0,0,0.08);
                    border: none;
                }
                .box-header {
                    background-color: transparent;
                    border-bottom: 1px solid #e9ecef;
                }
            "))
    ),
    sidebarMenu(
      id = "sidebar",
      menuItem("Dashboard", tabName = "intro", icon = icon("gauge")),
      menuItem("Data Input", tabName = "data_input", icon = icon("upload")),
      menuItem("Preprocessing", tabName = "preprocessing", icon = icon("sliders")),
      menuItem("Differential Expression", tabName = "deg_analysis", icon = icon("chart-line")),
      menuItem("Gene Set Enrichment", tabName = "geneset_enrichment", icon = icon("microscope")),
      menuItem("Gene Regulators", tabName = "gene_regulators", icon = icon("dna")),
      menuItem("Master Regulators", tabName = "master_regulators", icon = icon("crown")),
      menuItem("Literature Search", tabName = "pubmed_search", icon = icon("search")),
      menuItem("Gene Disease Network", tabName = "gene_disease_network", icon = icon("network-wired")),
      menuItem("Drug Database Search", tabName = "database_search", icon = icon("pills"))
    )
  ),
  dashboardBody(
    useShinyjs(),
    useShinyalert(force = TRUE),
    tags$head(
      tags$link(href = "https://fonts.googleapis.com/css2?family=Inter:wght@300;400;600&display=swap", rel = "stylesheet")
    ),
    tabItems(
      # Data Input Tab
      tabItem(tabName = "data_input",
              fluidRow(
                box(title = "Upload Expression Data", status = "primary", solidHeader = TRUE, 
                    fileInput("expression_file", "Choose CSV File", 
                              accept = c(".csv", ".txt", ".tsv")),
                    checkboxInput("normalized_data", "Normalized Data", FALSE),
                    radioButtons("file_type", "File Type",
                                 choices = c("CSV" = "csv", 
                                             "Tab-separated" = "tsv", 
                                             "Space-separated" = "space")),
                    checkboxInput("header", "File has header", TRUE)
                ),
                box(title = "Data Preview", status = "info", solidHeader = TRUE,
                    DTOutput("data_preview")
                )
              )
      ),
      
      # Pre_processing Tab
      tabItem(tabName = "preprocessing",
              fluidRow(
                box(title = "Preprocessing Options", status = "primary", solidHeader = TRUE,
                    selectInput("normalization_method", "Normalization Method",
                                choices = c("TMM" = "tmm", 
                                            "RLE" = "rle", 
                                            "Upper Quartile" = "upperquartile")),
                    sliderInput("filter_low_counts", "Filter Low Counts",
                                min = 0, max = 10, value = 1),
                    actionButton("preprocess_data", "Run Preprocessing")
                ),
                box(title = "Preprocessing Results", status = "info", solidHeader = TRUE,
                    DTOutput("preprocessed_data"),
                    downloadButton("download_preprocessed", "Download Preprocessed Data")
                )
              )
      ),
      
      # Differential Expression Analysis Tab
      tabItem(tabName = "deg_analysis",
              fluidRow(
                useShinyjs(),
                box(title = "DEG Analysis Setup", status = "primary", solidHeader = TRUE,
                    selectInput("data_type", "Select Data Type", choices = c("Raw Count Data" = "raw", "Preprocessed Data" = "processed", "Normalized Data"="normalized")),
                    selectInput("deg_method", "DEG Method",
                                choices = c("edgeR-Count" = "edger",
                                            "limma" = "voom")),
                    selectInput("comparison_group", "Select Comparison Group",
                                choices = c("Choose after data input")),
                    
                    
                    
                    # Horizontal layout for the Run button and Progress bar
                    fluidRow(
                      column(6, 
                             actionButton("run_deg", "Run DEG Analysis"),
                             br(),
                             br(),
                             numericInput("logfc_threshold", "Log Fold Change Threshold for Plots", 
                                          value = 0.5, min = 0),
                             numericInput("qvalue_threshold", "Adjusted_P-value Threshold for Plots", 
                                          value = 0.05, min = 0, max = 1)
                      ),
                      column(6, 
                             div(id = "progress_bar_container", style = "display:none; width: 100%;",
                                 div(style = "height: 20px; background-color: lightgray; width: 100%;",
                                     div(id = "progress_bar", style = "height: 100%; width: 0%; background-color: blue;")
                                 )
                             )
                      )
                    )
                ),
                box(title = "Differential Expression Results", status = "info", solidHeader = TRUE,
                    DTOutput("deg_results"),
                    downloadButton("download_deg", "Download Results")
                )
              ),
              
              fluidRow(
                box(title = "Volcano Plot", status = "success", solidHeader = TRUE, collapsible = TRUE,
                    plotlyOutput("volcano_plot"),
                    sliderInput("volcano_width", "Width (inches)", min = 4, max = 20, value = 8,width='400px'),
                    sliderInput("volcano_height", "Height (inches)", min = 4, max = 20, value = 8,width='400px'),
                    numericInput("volcano_res", "Resolution (dpi)", value = 300, min = 100),
                    downloadButton("download_volcano_tiff", "Download Volcano Plot (TIFF)")
                ),
                box(title = "Heatmap (Top 20)", status = "success", solidHeader = TRUE, collapsible = TRUE,
                    plotlyOutput("heatmap_plot"),
                    sliderInput("heatmap_width", "Width (inches)", min = 4, max = 20, value = 8,width='400px'),
                    sliderInput("heatmap_height", "Height (inches)", min = 4, max = 20, value = 8,width='400px'),
                    numericInput("heatmap_res", "Resolution (dpi)", value = 300, min = 100),
                    downloadButton("download_heatmap_tiff", "Download Heatmap (TIFF)")
                )
              )
      ),
      
      # Gene Regulators Tab
      tabItem(tabName = "gene_regulators",
              fluidRow(
                box(title = "Regulatory Network Analysis", status = "primary", solidHeader = TRUE,
                    selectInput("regulator_method", "Regulatory Analysis Method",
                                choices = c("WGCNA" = "wgcna")),
                    selectInput("data_type_reg", "Select Data Type", choices = c("Preprocessed Data" = "processed", "Normalized Data"="normalized")),
                    selectInput("comparison_group2", "Select Comparison Group",
                                choices = c("Choose after data input")),
                    actionButton("find_regulators", "Find Regulators"),
                    br(),
                    br(),
                    radioButtons("module_selection", "Select Module", choices = c("To See Modules, Run Analysis First")),
                    numericInput("top_regulators", "Number of Top Genes for Plots", 
                                 value = 10, min = 1),
                    uiOutput("group_filter_radio")
                    
                ),
                box(title = "Regulatory Network Results", status = "info", solidHeader = TRUE,
                    DTOutput("regulators_table"),
                    downloadButton("download_module_genes", "Download Module Genes")
                )
              ),
              
              fluidRow(
                box(title = "Module-Specific Heatmap", status = "success", solidHeader = TRUE,collapsible = TRUE,
                    plotlyOutput("module_heatmap"),
                    sliderInput("heatmap_width", "Width (inches)", min = 4, max = 20, value = 8,width='400px'),
                    sliderInput("heatmap_height", "Height (inches)", min = 4, max = 20, value = 8,width='400px'),
                    numericInput("heatmap_res", "Resolution (dpi)", value = 300, min = 100),
                    downloadButton("download_module_heatmap", "Download Heatmap (TIFF)")
                ),
                box(title = "Module-Specific Network", status = "success", solidHeader = TRUE,collapsible = TRUE,
                    plotOutput("module_network"),
                    sliderInput("network_width", "Width (inches)", min = 4, max = 20, value = 8,width='400px'),
                    sliderInput("network_height", "Height (inches)", min = 4, max = 20, value = 8,width='400px'),
                    numericInput("network_res", "Resolution (dpi)", value = 300, min = 100),
                    downloadButton("download_module_network", "Download Network (TIFF)")
                )
              )
      ),
      
      # Gene Set Enrichment Tab
      tabItem(tabName = "geneset_enrichment",
              fluidRow(
                box(title = "Enrichment Analysis", status = "primary", solidHeader = TRUE,
                    selectInput("enrichment_type", "Enrichment Type",
                                choices = c("GO Biological Process" = "gobp",
                                            "GO Molecular Function" = "gomf",
                                            "GO Cellular Component" = "gocc",
                                            "KEGG Pathways" = "kegg")),
                    numericInput("GSEA_pvalue_threshold", "GSEA Pvalue Cutoff", 
                                 value = 0.1, min = 0, max = 1),
                    actionButton("run_enrichment", "Run Enrichment"),
                    br(),
                    br(),
                    numericInput("top_enriched", "Top Enriched Terms", 
                                 value = 10, min = 1)
                ),
                box(title = "Enrichment Results", status = "info", solidHeader = TRUE,
                    DTOutput("enrichment_table"),
                    downloadButton("download_enrichment_results", "Download GSEA Results")
                )),
              fluidRow(
                box(title = "Enrichment Plot", status = "success", solidHeader = TRUE,collapsible = TRUE,
                    plotOutput("enrichment_plot", height = 600),
                    sliderInput("enrichment_plot_width", "Width (inches)", min = 4, max = 20, value = 10,width='400px'),
                    sliderInput("enrichment_plot_height", "Height (inches)", min = 4, max = 20, value = 8,width='400px'),
                    numericInput("enrichment_plot_res", "Resolution (dpi)", value = 300, min = 100),
                    downloadButton("download_enrichment_plot", "Download Enrichment Plot (TIFF)")
                ),
                box(title = "Enrichment Network", status = "success", solidHeader = TRUE,collapsible = TRUE,
                    plotOutput("enrichment_network", height = 600),
                    sliderInput("network_width", "Width (inches)", min = 4, max = 20, value = 10,width='400px'),
                    sliderInput("network_height", "Height (inches)", min = 4, max = 20, value = 8,width='400px'),
                    numericInput("network_res", "Resolution (dpi)", value = 300, min = 100),
                    downloadButton("download_enrichment_network", "Download Network (TIFF)")
                )
                
              )
      ),
      
      # Gene Set MRA Tab
      tabItem(tabName = "master_regulators",
              fluidRow(
                box(title = "Master Regulators Analysis", status = "primary", solidHeader = TRUE, width = 12,
                    selectInput("MRA_type", "MRA Type",
                                choices = c("Corto" = "corto")),
                    selectInput("data_type_mra", "Select Data Type", choices = c("Raw Count Data" = "raw", "Preprocessed Data" = "processed", "Normalized Data"="normalized")),
                    selectInput("mra_group", "Select Group of Interest:",
                                choices = NULL),
                    actionButton("run_mra", "Run MRA")
                )
              ),
              fluidRow(
                column(width = 6,
                       box(title = "List of All Master Regulators", status = "success", solidHeader = TRUE, collapsible = TRUE, width = 12,
                           DTOutput("mra_table"),
                           downloadButton("download_mra_table", "Download MRA Results")
                       )
                ),
                column(width = 6,
                       box(title = "Top 5 Master Regulators", status = "success", solidHeader = TRUE, collapsible = TRUE, width = 12,
                           plotOutput("mra_plot1", height = "600px"),
                           sliderInput("mra_plot2_width", "Width (inches)", min = 4, max = 20, value = 8, width='400px'),
                           sliderInput("mra_plot2_height", "Height (inches)", min = 4, max = 20, value = 8, width='400px'),
                           numericInput("mra_plot2_res", "Resolution (dpi)", value = 300, min = 100),
                           downloadButton("download_mra_plot1", "Download MRA Plot (TIFF)")
                       )
                )
              )
      ),
      
      #Pubmed and CT
      tabItem(tabName = "pubmed_search",
              fluidRow(
                box(title = "Search Options", status = "primary", solidHeader = TRUE,
                    radioButtons("search_source", "Select Source of Terms:",
                                 choices = list("DEG Analysis (Top 20 Genes)" = "deg_analysis",
                                                "Master Regulators" = "master_regulators",
                                                "Gene Regulators (Top 20 Module-Genes)" = "gene_regulators")),
                    uiOutput("module_selection"),
                    textInput("disease_filter", "Specify Disease:", value = "Specify a disease"),
                    selectInput("search_field", "Search in:", 
                                choices = c("Title/Abstract" = "[Title/Abstract]", 
                                            "MeSH Terms" = "[MeSH Terms]"), 
                                selected = "[Title/Abstract]"),
                    radioButtons("module_selection2", "Select Module", choices = c("Only for Gene Regulators")),
                    checkboxInput("clinical_trial_only", "Limit PubMed Search to Clinical Trials", FALSE),
                    actionButton("search", "Search in Pubmed & Clinical Trials"),
                    downloadButton("download_results", "Download Pubmed Results as CSV"),
                    downloadButton("download_results1", "Download Clinical Trials Results as CSV"),
                    verbatimTextOutput("debug_info")
                )
              ),
              fluidRow(
                tabBox(id = "results_tabs", width = 12, 
                       tabPanel("PubMed Results", DTOutput("pubmed_table")),
                       tabPanel("Clinical Trials Results", DTOutput("clinical_trials_table"))
                )
              )
      ),
      
      
      #Gene Disease Network Tab
      tabItem(tabName = "gene_disease_network",
              fluidRow(
                box(title = "Gene Disease Network Analysis", status = "primary", solidHeader = TRUE, width = 12,
                    selectInput("gene_source", "Select Source of Genes:",
                                choices = c("DEG Analysis" = "deg",
                                            "Master Regulators" = "mra",
                                            "Gene Regulators (WGCNA)" = "wgcna")),
                    numericInput("qvalue_threshold_DEGEN", "q-value Threshold for the Results", 
                                 value = 0.05, min = 0, max = 1),
                    uiOutput("module_selection"),
                    radioButtons("module_selection3", "Select Module", choices = c("Only for Gene Regulators")),
                    sliderInput("centrality_filter", "Filter Nodes by Centrality:", min = 0, max = 1, value = 0.5, step = 0.05),
                    sliderInput("top_n_results", "Top N Diseases for Network:", min = 1, max = 50, value = 10, step = 1),
                    actionButton("run_disease_network", "Run Analysis")
                )
              ),
              fluidRow(
                column(width = 6,
                       box(title = "Disease Enrichment Results", status = "info", solidHeader = TRUE, width = 12,
                           DTOutput("disease_enrichment_table"),
                           downloadButton("download_disease_results", "Download Results")
                       )
                ),
                column(width = 6,
                       box(title = "Disease-Gene Network", status = "success", solidHeader = TRUE, width = 12,
                           visNetworkOutput("disease_network_plot", height = "600px"),
                           sliderInput("dn_plot_width", "Width (inches)", min = 4, max = 20, value = 8, width='400px'),
                           sliderInput("dn_plot_height", "Height (inches)", min = 4, max = 20, value = 8, width='400px'),
                           numericInput("dn_plot_res", "Resolution (dpi)", value = 300, min = 100),
                           downloadButton("download_disease_network", "Download Network (TIFF)")
                       )
                )
              )
      ),
      
      tabItem(tabName = "intro",
              fluidRow(
                box(title = "Welcome to APIomics", status = "primary", solidHeader = TRUE, width = 12, height = "100%",
                    p("APIomics is a bioinformatics analysis pipeline designed to process and analyze high-throughput expression data.",
                      style = "font-size: 16px; font-weight: bold;"),
                    p("This tool enables users to perform preprocessing, differential expression analysis, gene set enrichment, regulatory network analysis, pubmed and clinical trial search.",
                      style = "font-size: 14px;"),
                    p("Navigate through the tabs to explore various functionalities and start your analysis.",
                      style = "font-size: 14px;"),
                    p("If you need help, please refer to the user guide or contact support. \nEmail:morteza.hajihosseini@appliedpharma.ca",
                      style = "font-size: 14px;"),
                    p(br()),
                    p("Instructions:",
                      style = "font-size: 16px; font-weight: bold;"),
                    p("Input Dataset Format: Make sure your samples are in the rows and genes in the columns. It is important to have the condition (grouping) variable at the end of the file as the last column.",
                      style = "font-size: 14px"),
                    p("Preprocessing: The preprocessing step performs the normalization and low count filtering.",
                      style = "font-size: 14px; "),
                    p("Differential Expression: The DEG analysis can be done on the raw counts or normalized data. The combination of Log Fold Change Threshold (more than) and Adjusted_P-value Threshold (less than) will be used to summerize the data for the plots.",
                      style = "font-size: 14px; "),
                    p("Gene Set Enrichment (GSEA): The GSEA will use the gene list from DEG analysis to perform gene set enrichment analysis using four different human spicy databases.",
                      style = "font-size: 14px; "),
                    p("Gene Regulators (WGCNA): The weighted gene co-expression network analysis (WGCNA) is a widely used method for describing the correlation patterns of genes across a large set of samples.",
                      style = "font-size: 14px; "),
                    p("Master Regulators (Corto): The corto algorithm will run the gene network inference and master regulator analysis (MRA).",
                      style = "font-size: 14px; "),
                    p("Search in PubMed & Clinical Trials: In this section, you can call the result from DEG Analysis, Master Regulators, and WGCNA to search in the Pubmed and Clinical Trial databases.",
                      style = "font-size: 14px; "),
                    p("Gene Disease Network: In this section, you can call the result from DEG Analysis, Master Regulators, and WGCNA to find the connections for various diseases using the DISGENET database. \nDisGeNET is a discovery platform for the dynamical exploration of human diseases and their genes",
                      style = "font-size: 14px; "),
                    p("APIomics Pipeline",
                      style = "font-size: 16px; font-weight: bold;"),
                    p("The gray boxes are under development.",
                      style = "font-size: 16px; font-weight:"),
                    tags$img(src = "static/flowchart.jpg", height = "50%", width = "70%")
                    
                    
                )
              )
      ),
      
      tabItem(tabName = "database_search",
              fluidRow(
                box(title = "Search Options", status = "primary", solidHeader = TRUE,
                    selectInput("gene_source_db", "Select Gene List:", 
                                choices = list("DEG Analysis (Top 20 Genes)" = "deg_analysis",
                                               "Master Regulators" = "master_regulators",
                                               "Gene Regulators (Top 20 Module-Genes)" = "gene_regulators")),
                    radioButtons("module_selection_db", "Select Module", choices = c("Only for Gene Regulators")),
                    actionButton("search_database", "Search"),
                    actionButton("selected_molecule", "Showe Molecule"))
                 
              ),
              fluidRow(
                tabBox(id = "database_tabs", width = 12,
                       tabPanel("ChEMBL Results", DTOutput("chembl_results"),
                                downloadButton("download_chembl", "Download ChEMBL Results")),
                       tabPanel("BindingDB Results", DTOutput("bindingdb_results"),
                                downloadButton("download_bindingdb", "BindingDB Results"))
                )
              ),
              fluidRow(
                box(title = "Molecule 3D Viewer", width = 12, solidHeader = TRUE, status = "info",
                    uiOutput("molecule_3d_viewer")
                )
              )
      )
      
    )
  )
)





server <- function(input, output, session) {
  # Reactive values to store data across tabs



  rv <- reactiveValues(
    raw_data = NULL,
    preprocessed_data = NULL,
    Normalized_data=NULL,
    deg_data=NULL,
    deg_results = NULL,
    logcpm=NULL,
    reg_data=NULL,
    wgcna_modules=NULL,
    module_genes=NULL,
    mod_mat=NULL,
    group_data=NULL,
    moduleMembership=NULL,
    TOM=NULL,
    net_colors=NULL,
    predicted = NULL, 
    regulons = NULL,
    mra_data=NULL,
    mra_results=NULL,
    all_results=NULL
  )
  
  # Initialize shinyjs and disable tabs except "Data Input"
  #shinyjs::disable(selector = ".sidebar-menu a[data-value='preprocessing']")
  #shinyjs::disable(selector = ".sidebar-menu a[data-value='deg_analysis']")
  #shinyjs::disable(selector = ".sidebar-menu a[data-value='geneset_enrichment']")
  #shinyjs::disable(selector = ".sidebar-menu a[data-value='gene_regulators']")
  #shinyjs::disable(selector = ".sidebar-menu a[data-value='master_regulators']")
  
  # Data Input Tab
  observeEvent(input$expression_file, {
    req(input$expression_file)
    
    # Read file based on selected type
    if (input$file_type == "csv") {
      data <- read.csv(input$expression_file$datapath, 
                       header = input$header)
      rownames(data)<-data[,1]
      data<-data[,-1]
    } else if (input$file_type == "tsv") {
      data <- read.delim(input$expression_file$datapath, 
                         header = input$header)
    } else {
      data <- read.table(input$expression_file$datapath, 
                         header = input$header)
    }
    
    names(data)[dim(data)[2]]<-"Group"
    rv$raw_data <- data
    
    # Preview data
    output$data_preview <- renderDT({
      datatable(head(data, 50), 
                options = list(scrollX = TRUE))
    })
    
  })
  
  observe({
    if (input$normalized_data) {
      shinyjs::disable(selector = ".sidebar-menu a[data-value='preprocessing']")
      shinyjs::enable(selector = ".sidebar-menu a[data-value='deg_analysis']")
      shinyjs::enable(selector = ".sidebar-menu a[data-value='gene_regulators']")
      shinyjs::enable(selector = ".sidebar-menu a[data-value='pubmed_search']")
    } else {
      shinyjs::enable(selector = ".sidebar-menu a[data-value='preprocessing']")
    }
  })
  
  
  
  # Preprocessing Tab
  observeEvent(input$preprocess_data, {
    req(rv$raw_data)

    # preprocessing logic
    rv$raw_data2<-rv$raw_data[,-dim(rv$raw_data)[2]]
    rv$raw_data2<-rv$raw_data2[,colSums(rv$raw_data2)>0]
    rv$preprocessed_data <- DGEList(rv$raw_data2)
    rv$logcpm<-cpm(rv$preprocessed_data,log=T)
    keep<-rowSums(cpm(rv$preprocessed_data,log=T,keep.lib.sizes=TRUE)>input$filter_low_counts)>=3
    rv$preprocessed_data=rv$preprocessed_data[keep,]
    
    if(input$normalization_method=="tmm")
    {
      tryCatch(
        expr={rv$preprocessed_data <-calcNormFactors(rv$preprocessed_data, method="TMM")
        rv$preprocessed_data<-data.frame(cpm(rv$preprocessed_data$counts,log=TRUE))
        rv$Normalized_data<-data.frame(rv$preprocessed_data ,Group=rv$raw_data[,dim(rv$raw_data)[2]])
        rownames(rv$Normalized_data)=rownames(rv$raw_data)
        },
        error=function(e){
          shinyalert( title = "ERROR", text = "Try a different normalization method")
        }
      )
      
    }else if (input$normalization_method=="rle"){
      tryCatch(
        expr={rv$preprocessed_data <-calcNormFactors(rv$preprocessed_data, method="RLE")
        rv$preprocessed_data<-data.frame(cpm(rv$preprocessed_data$counts,log=TRUE))
        rv$Normalized_data<-data.frame(rv$preprocessed_data ,Group=rv$raw_data[,dim(rv$raw_data)[2]])
        rownames(rv$Normalized_data)=rownames(rv$raw_data)
        },
        error=function(e){
          shinyalert( title = "ERROR", text = "Try a different normalization method")
        }
      )
    }else if (input$normalization_method=="upperquartile"){
      tryCatch(
        expr={rv$preprocessed_data <-calcNormFactors(rv$preprocessed_data, method="upperquartile")
        rv$preprocessed_data<-data.frame(cpm(rv$preprocessed_data$counts,log=TRUE))
        rv$Normalized_data<-data.frame(rv$preprocessed_data ,Group=rv$raw_data[,dim(rv$raw_data)[2]])
        rownames(rv$Normalized_data)=rownames(rv$raw_data)
        },
        error=function(e){
          shinyalert( title = "ERROR", text = "Try a different normalization method")
        }
      )
    }
    
    output$preprocessed_data <- renderDT({
      datatable(head(rv$Normalized_data, 50), 
                options = list(scrollX = TRUE))
    })
  })
  
  output$download_preprocessed <- downloadHandler(
    filename = function() {
      paste("preprocessed_data", Sys.Date(), ".csv", sep = "")
    },
    content = function(file) {
      write.csv(rv$Normalized_data, file, row.names = TRUE)
    }
  )
  
  
  
  
  
  # DEG Analysis Tab
  deg_alert_shown <- reactiveVal(FALSE)
  
  observe({
    req(input$data_type)
    if (input$sidebar == "deg_analysis") {
      isolate({

          if (input$data_type == "raw" & !is.null(rv$raw_data) & !input$normalized_data) {
            non_numeric_cols <- names(rv$raw_data)[!sapply(rv$raw_data, is.numeric)]
            updateSelectInput(session, "comparison_group", choices = non_numeric_cols)
            rv$deg_data <- rv$raw_data
            shinyalert("Note", "Select edgeR-Count method for Raw Count Data.")
            deg_alert_shown(TRUE)
          } else if (input$data_type == "processed" & !is.null(rv$Normalized_data) & !input$normalized_data) {
            non_numeric_cols <- names(rv$Normalized_data)[!sapply(rv$Normalized_data, is.numeric)]
            updateSelectInput(session, "comparison_group", choices = non_numeric_cols)
            rv$deg_data <- rv$Normalized_data
            shinyalert("Note", "Select limma-voom method for Preprocessed Data.")
            deg_alert_shown(TRUE)
          } else if (input$data_type == "normalized" & input$normalized_data & !is.null(rv$raw_data)) {
            non_numeric_cols <- names(rv$raw_data)[!sapply(rv$raw_data, is.numeric)]
            updateSelectInput(session, "comparison_group", choices = non_numeric_cols)
            shinyalert("Note", "Normalized data provided in the Data Input tab. \nSelect limma method for Normalized Data.")
            rv$deg_data <- rv$raw_data
            deg_alert_shown(TRUE)
          }
          
      })
    }
  })



  
  
  observeEvent(input$run_deg, {
    req(rv$deg_data)
    
    shinyjs::enable(selector = ".sidebar-menu a[data-value='geneset_enrichment']")
    shinyjs::enable(selector = ".sidebar-menu a[data-value='master_regulators']")

    
    group_column <- input$comparison_group
    group_data <- rv$deg_data[,group_column]
    deg_data1<-rv$deg_data[,-dim(rv$deg_data)[2]]
    rv$group_data=data.frame(Samples=rownames(deg_data1),Group=group_data)
    
    # DEG analysis
    set.seed(123)
    
    if(input$deg_method=="voom"){
      group <- rv$deg_data[,input$comparison_group]
      data_matrix <- as.matrix(deg_data1[,!colnames(rv$deg_data) %in% "Group"])
      design <- model.matrix(~1+group)
      y <- voom(t(data_matrix), design)
      fit <- lmFit(y, design)
      fit2 <- eBayes(fit)
      top.table <- topTable(fit2, adjust.method = "BH", sort.by = "p", n = Inf)
      top.table<-data.frame(top.table)
      #top.table=round(top.table,4)
      
    }else if (input$deg_method=="edger"){
      raw_data2<-rv$deg_data[,-dim(rv$deg_data)[2]]
      raw_data2<-raw_data2[,colSums(raw_data2)>0]
      
      counts=data.frame(t(raw_data2))
      group=rv$deg_data[,dim(rv$deg_data)[2]]
      d <- DGEList(counts,group=c(group))
      d$samples$lib.size <- colSums(d$counts)
      d <- calcNormFactors(d)
      design.mat <- model.matrix(~ 0 + d$samples$group)
      colnames(design.mat) <- levels(d$samples$group)
      
      # Progress bar starts here
      withProgress(message = 'Performing DEG Analysis', value = 0, {
        incProgress(0.3, detail = "Estimating Common Dispersion")
        d2 <- estimateGLMCommonDisp(d, design.mat)
        
        incProgress(0.3, detail = "Estimating Tagwise Dispersion")
        d2 <- estimateGLMTagwiseDisp(d2, design.mat)
        
        incProgress(0.3, detail = "Fitting Model")
        fit <- glmFit(d2, design.mat)
        lrt12 <- glmLRT(fit, contrast = c(1, -1))
        top.table <- data.frame(topTags(lrt12, n = Inf))
      })
      
      names(top.table)[c(1,4,5)]=c("logFC","P.Value","adj.P.Val")
      
    }
    
    
    
    deg_results <- top.table
    
    rv$deg_results <- data.frame(deg_results)
    
    output$deg_results <- renderDT({
      datatable(data.frame(deg_results), 
                options = list(scrollX = TRUE,pageLength = 5))
    })
    
    # Volcano Plot
    output$volcano_plot <- renderPlotly({
      if(!is.null(deg_results)) {
        deg_results$Significance <- ifelse(
          deg_results$adj.P.Val < input$qvalue_threshold & 
            abs(deg_results$logFC) > input$logfc_threshold, 
          "Significant", "Not Significant"
        )
        plot_ly(
          data = deg_results,
          x = ~logFC, 
          y = ~-log10(adj.P.Val),
          text = ~paste(
            "Gene: ", rownames(deg_results), 
            "<br>Log2 Fold Change: ", round(logFC, 2),
            "<br>Adjusted P-value: ", round(adj.P.Val, 4)
          ),
          hoverinfo = 'text',
          mode = 'markers',
          color = ~Significance,
          colors = c('Significant' = 'red', 'Not Significant' = 'black'),
          type = 'scatter'
        ) %>%
          layout(
            title = "Interactive Volcano Plot",
            xaxis = list(title = "Log2 Fold Change"),
            yaxis = list(title = "-log10 Adjusted P-value"),
            shapes = list(
              list(
                type = "line",
                x0 = -input$logfc_threshold,
                x1 = -input$logfc_threshold,
                y0 = 0,
                y1 = max(-log10(deg_results$adj.P.Val)),
                line = list(color = "blue", dash = "dot")
              ),
              list(
                type = "line",
                x0 = input$logfc_threshold,
                x1 = input$logfc_threshold,
                y0 = 0,
                y1 = max(-log10(deg_results$adj.P.Val)),
                line = list(color = "blue", dash = "dot")
              )
            )
          )
      }
    })
    
    
    
    # Volcano Plot TIFF download handler
    output$download_volcano_tiff <- downloadHandler(
      filename = function() {
        paste("volcano_plot", Sys.Date(), ".tiff", sep = "")
      },
      content = function(file) {
        width <- input$volcano_width
        height <- input$volcano_height
        res <- input$volcano_res
        tiff(file, width = width, height = height, units = "in", res = res)
        deg_results$Significance <- ifelse(
          deg_results$adj.P.Val < input$qvalue_threshold & 
            abs(deg_results$logFC) > input$logfc_threshold, 
          "Significant", "Not Significant"
        )
        plot <- ggplot(deg_results, aes(x = logFC, y = -log10(adj.P.Val), color = Significance)) +
          geom_point() +
          theme_minimal() +
          labs(title = "Volcano Plot", x = "Log2 Fold Change", y = "-log10 Adjusted P-value") +
          scale_color_manual(values = c("Significant" = "red", "Not Significant" = "black"))
        print(plot)
        dev.off()
      }
    )
    
    
    output$heatmap_plot <- renderPlotly({
      req(rv$deg_results, rv$deg_data)
      
        # Validate data
      if(is.null(rv$deg_results) || is.null(rv$deg_data)) {
        return(NULL)
      }
      
      # Select top significant genes
      top_genes <- rv$deg_results %>%
        filter(adj.P.Val < input$qvalue_threshold & 
                 abs(logFC) > input$logfc_threshold) %>%
        arrange(adj.P.Val) %>%
        head(20)
      
      
      # Ensure we have genes
      if(nrow(top_genes) == 0) {
        shinyalert(
          title = "No Significant Genes Found",
          text = "No DEGs were found based on the selected thresholds. Consider adjusting the log fold change or p-value thresholds.",
          type = "warning",
          confirmButtonText = "OK"
        )
        return(NULL)
      }
 
      # Prepare data for plotting
      # Extract only numeric columns and the group column
      numeric_cols <- names(rv$deg_data)[sapply(rv$deg_data, is.numeric)]
      group_col <- input$comparison_group
      
      # Subset data to include only top genes
      heatmap_data <- rv$deg_data %>%
        dplyr::select(all_of(c(intersect(rownames(top_genes), numeric_cols), group_col)))
      heatmap_data=data.frame(heatmap_data)
      heatmap_data=heatmap_data[order(heatmap_data[,dim(heatmap_data)[2]]),]
      
      
      col <- colorRampPalette(rev(RColorBrewer::brewer.pal(n = 10, name = "RdYlBu")))(100)
      ByPal <- colorRampPalette(c('yellow','purple'))
      # Create plot
      heatmap_plot<-heatmaply(heatmap_data[,-dim(heatmap_data)[2]],column_text_angle=90,Rowv = F, Colv = F,colors =col ,scale = "column",
                              row_side_colors=heatmap_data[,dim(heatmap_data)[2]],row_side_palette= ByPal,
                              showticklabels = c(TRUE, FALSE))
    })
    
    
    
    # Heatmap TIFF download handler
    output$download_heatmap_tiff <- downloadHandler(
      filename = function() {
        paste("heatmap_plot", Sys.Date(), ".tiff", sep = "")
      },
      content = function(file) {
        width <- input$heatmap_width
        height <- input$heatmap_height
        res <- input$heatmap_res

        # Validate data
        if(is.null(rv$deg_results) || is.null(rv$deg_data)) {
          return(NULL)
        }
        
        # Select top significant genes
        top_genes <- rv$deg_results %>%
          filter(adj.P.Val < input$qvalue_threshold & 
                   abs(logFC) > input$logfc_threshold) %>%
          arrange(adj.P.Val) %>%
          head(20)
        
        # Ensure we have genes
        if(nrow(top_genes) == 0) {
          return(NULL)
        }
        
        # Prepare data for plotting
        # Extract only numeric columns and the group column
        numeric_cols <- names(rv$deg_data)[sapply(rv$deg_data, is.numeric)]
        group_col <- input$comparison_group
        
        # Subset data to include only top genes
        heatmap_data <- rv$deg_data %>%
          dplyr::select(all_of(c(intersect(rownames(top_genes), numeric_cols), group_col)))
        heatmap_data=data.frame(heatmap_data)
        heatmap_data=heatmap_data[order(heatmap_data[,dim(heatmap_data)[2]]),]
        
        
        col <- colorRampPalette(rev(RColorBrewer::brewer.pal(n = 10, name = "RdYlBu")))(100)
        
        
        # Generate unique colors for groups
        unique_groups <- unique(heatmap_data[,group_col])
        group_colors <- setNames(
          RColorBrewer::brewer.pal(min(length(unique_groups), 8), "Set2")[1:length(unique_groups)],
          unique_groups
        )
        row_side_colors <- group_colors[as.character(heatmap_data[,group_col])]
        
        tiff(file, width = width, height = height, units = "in", res = res)
        # Create heatmap using gplots
        heatmap.2(
          data_matrix, 
          scale = "column", 
          col = col, 
          Rowv = FALSE, 
          Colv = FALSE, 
          dendrogram = "none", 
          trace = "none", 
          margins = c(10, 15),
          main = "Gene Expression Heatmap",
          RowSideColors = row_side_colors,
          keysize = 1,  
          density.info = "none",  
          key = TRUE  
        )
        
        # Add a legend for group colors
        legend(
          "topright", 
          legend = names(group_colors), 
          col = group_colors, 
          lty = 1, 
          lwd = 10, 
          xpd = TRUE,
          title = group_col,
          cex = 0.8
        )
        
        dev.off()
      }
    )
    
    
  })
  
  
  
  # Download DEG Results
  output$download_deg <- downloadHandler(
    filename = function() {
      paste("deg_results_", Sys.Date(), ".csv", sep="")
    },
    content = function(file) {
      write.csv(rv$deg_results, file, row.names = TRUE)
    }
  )
  
  
  
 
  # Gene Regulators Tab
  deg_alert_shown2 <- reactiveVal(FALSE)
  
  observe({
    req(input$data_type_reg)
    if (input$sidebar == "gene_regulators") {
      isolate({
        if (!is.null(rv$raw_data) & !input$normalized_data & is.null(rv$Normalized_data)) {
          shinyalert("Note", "WGCNA needs preprocessed or normalized data.")
          deg_alert_shown(TRUE)
        } else if (input$data_type_reg == "processed" & !is.null(rv$raw_data) & !is.null(rv$Normalized_data) & !input$normalized_data) {
            non_numeric_cols <- names(rv$Normalized_data)[!sapply(rv$Normalized_data, is.numeric)]
            updateSelectInput(session, "comparison_group2", choices = non_numeric_cols)
            rv$reg_data <- rv$Normalized_data
          } else if (input$data_type_reg == "normalized" & input$normalized_data & !is.null(rv$raw_data)) {
            non_numeric_cols <- names(rv$raw_data)[!sapply(rv$raw_data, is.numeric)]
            updateSelectInput(session, "comparison_group2", choices = non_numeric_cols)
            rv$reg_data <- rv$raw_data
          }

      })
    }
  })
  
  
  
  
  # Gene Regulators Tab
  observeEvent(input$find_regulators, {
    
    req(rv$reg_data)
    group_column <- input$comparison_group2
    group_data <- rv$reg_data[, group_column]
    
    reg_data1<-rv$reg_data[,-dim(rv$reg_data)[2]]
    rv$group_data=data.frame(Samples=rownames(reg_data1),Group=group_data)
    
    
    
    if(input$regulator_method=="wgcna"){
      
      datExpr <- reg_data1[,sapply(reg_data1, is.numeric)]
      
      # Remove rows with NAs
      datExpr <- datExpr[complete.cases(datExpr), ]
      
      # Ensure numeric data
      datExpr=as.matrix(datExpr)
      
      
      withProgress(message = 'Performing WGCNA Network Analysis', value = 0, {
        
        incProgress(0.1, detail = "Selecting Soft Threshold")
        powers <- c(1:10, seq(from = 12, to = 20, by = 2))
        sft <- pickSoftThreshold(
          datExpr,
          dataIsExpr = TRUE,
          corFnc = cor,
          networkType = "signed"
        )
        
        
        # Select soft power
        softPower <- sft$powerEstimate
        if(is.na(softPower)) softPower <- 6
        
        incProgress(0.3, detail = "Preparing Network Construction")
        netwk <- blockwiseModules(datExpr,                # <= input here
                                  
                                  # == Adjacency Function ==
                                  power = softPower,                # <= power here
                                  networkType = "signed",
                                  
                                  # == Tree and Block Options ==
                                  deepSplit = 2,
                                  pamRespectsDendro = F,
                                  # detectCutHeight = 0.75,
                                  minModuleSize = 30,
                                  maxBlockSize = 4000,
                                  
                                  # == Module Adjustments ==
                                  reassignThreshold = 0,
                                  mergeCutHeight = 0.25,
                                  
                                  # == TOM == Archive the run results in TOM file (saves time)
                                  saveTOMs = F,
                                  saveTOMFileBase = "ER",
                                  
                                  # == Output Options
                                  numericLabels = T,
                                  randomSeed = 1234,
                                  verbose = 3)
        
        
        
        incProgress(0.3, detail = "Generating Module Eigengenes")
        
        module_eigengenes <- netwk$MEs
        module_colors <- labels2colors(netwk$colors)
        colnames(module_eigengenes) <- unique(module_colors)
        
        #For Network Plot
        rv$net_colors<-labels2colors(netwk$colors)
        
        rv$moduleMembership <- data.frame(cor(datExpr, module_eigengenes, use = "p"))
        
        # Prepare group data
        module_df <- data.frame(
          module_eigengenes, 
          rv$group_data
        ) 
        
        rv$module_colors <- unique(module_colors)
        
        
        module_df_melt=melt(module_df,id.vars=c("Samples","Group"))
       
      
      incProgress(0.3, detail = "Calculate the Topological Overlap Matrix")
      gene_module_key <- tibble::enframe(netwk$colors, name = "gene", value = "module") %>%
        # Let's add the `ME` part so its more clear what these numbers are and it matches elsewhere
        dplyr::mutate(module =labels2colors(netwk$colors))
      
      rv$wgcna_modules <-gene_module_key[order(gene_module_key$module),] 
      
      # Calculate the adjacency matrix
      adjacency <- adjacency(datExpr, power = softPower, type = "signed")
      
      # Calculate the Topological Overlap Matrix (TOM)
      rv$TOM <- TOMsimilarity(adjacency)
      
       incProgress(1, detail = "Complete")
      })
      
    }
    
    updateRadioButtons(
      session, 
      "module_selection", 
      choices = rv$module_colors
    )
    
    output$regulators_table <- renderDT({
      req(rv$wgcna_modules)
      # Filter genes for the selected module and create a data frame
      module_genes_df <- rv$wgcna_modules %>%
        filter(module == input$module_selection) %>%
        select(gene)
      
      datatable(module_genes_df, 
                colnames = c("Genes in Module"),
                options = list(
                  scrollX = TRUE, 
                  pageLength = 5  # Increase number of rows displayed
                ))
    })
    # Add download handler for module genes
    output$download_module_genes <- downloadHandler(
      filename = function() {
        paste0("module_genes_", input$module_selection, "_", Sys.Date(), ".csv")
      },
      content = function(file) {
        module_genes_df <- rv$wgcna_modules %>%
          filter(module == input$module_selection) %>%
          select(gene)
        write.csv(module_genes_df, file, row.names = FALSE)
      }
    )
    
  })
  
  output$group_filter_radio <- renderUI({
    req(rv$group_data)
    unique_groups <- unique(rv$group_data$Group)
    radioButtons("selected_group", "Filter by Group (Heatmap only):",
                 choices = c("All Groups", unique_groups),
                 selected = "All Groups")
  })
  
  # module heatmap rendering
  output$module_heatmap <- renderPlotly({
    req(rv$group_data)
    req(rv$wgcna_modules)
    req(rv$group_data)
    req(rv$reg_data)
    req(input$module_selection)
    req(rv$moduleMembership)
    req(input$top_regulators)
    
    reg_data1<-rv$reg_data[,-dim(rv$reg_data)[2]]
    
    module_genes_top <- rv$moduleMembership%>%
      select(input$module_selection)
    module_genes_top<-data.frame(Genes=rownames(module_genes_top),module_genes_top)
    
    module_genes_top2<-module_genes_top[order(module_genes_top[,2],decreasing =TRUE),]
    
    
    topGenes<-module_genes_top2[,1][1:input$top_regulators]
    
    # Filter by selected group if not "All Groups"
    if(input$selected_group != "All Groups") {
      group_indices <- which(rv$group_data$Group == input$selected_group)
      reg_data1 <- reg_data1[group_indices,]
    }
    
    
    # Prepare group information
    mod_mat <- reg_data1 %>%
      dplyr::select(all_of(topGenes)) %>%
      as.matrix()
    
    group_data <- rv$reg_data[rownames(mod_mat), ncol(rv$reg_data), drop = FALSE]
    
    # Set up the gene expression data frame
    
    
    mod_mat <- data.frame(mod_mat, group_data)
    mod_mat <- mod_mat[order(mod_mat[,dim(mod_mat)[2]]),]
    
    
    col <- colorRampPalette(rev(RColorBrewer::brewer.pal(n = 10, name = "RdYlBu")))(100)
    ByPal <- colorRampPalette(c('yellow','purple'))
    
    # Heatmap with Plotly
    heatmaply::heatmaply(
      mod_mat[,-dim(mod_mat)[2]],
      Rowv = FALSE,
      Colv = FALSE,
      colors =col,
      main = paste("Top", input$top_regulators ,"Genes, Module:", input$module_selection),column_text_angle=90,
      showticklabels = c(TRUE, FALSE)
    )%>%
      layout(
        legend = list(
          orientation = "h", # Horizontal legend placement
          x = 0.5, y = -0.5, # Adjust legend position
          xanchor = "center",
          yanchor = "botttom"
        )
      )
  })
  
  
  
  # Network with Plotly
  output$module_network <- renderPlot({
    req(rv$TOM)
    req(rv$wgcna_modules)
    req(input$module_selection)
    req(rv$net_colors)
    req(rv$reg_data)
    req(input$top_regulators)
    
    reg_data1<-rv$reg_data[,-dim(rv$reg_data)[2]]
    
    # Filter by selected group if not "All Groups"
    if(input$selected_group != "All Groups") {
      group_indices <- which(rv$group_data$Group == input$selected_group)
      reg_data1 <- reg_data1[group_indices,]
    }
    
    datExpr <- reg_data1[,sapply(reg_data1, is.numeric)]
    datExpr <- datExpr[complete.cases(datExpr), ]
    datExpr=as.matrix(datExpr)
    
    # Filter genes for the selected module and create a data frame
    module_genes <- which(rv$net_colors==input$module_selection)
    gene_names <- colnames(datExpr)[module_genes]
    TOM_module <- rv$TOM[module_genes, module_genes]
    network <- graph_from_adjacency_matrix(TOM_module, mode = "undirected", weighted = TRUE)
    V(network)$name <- gene_names
    # Calculate degree centrality for all nodes
    V(network)$degree <- degree(network, mode = "all")
    # Order nodes by degree centrality
    top_genes <- order(V(network)$degree, decreasing = TRUE)[1:input$top_regulators]
    filtered_network <- induced_subgraph(network, vids = top_genes)
    network<-filtered_network
    network_tidy <- as_tbl_graph(network)
    
    ggraph(network_tidy, layout = "fr") +
      geom_edge_link(aes(edge_alpha = weight), show.legend = FALSE) +
      geom_node_point(aes(color = input$module_selection), size = 3) +
      geom_node_text(aes(label = name), size = 4,repel=T, max.overlaps = 50) + 
      scale_color_manual(values =input$module_selection) +            
      theme_void() +theme(legend.position = "none")+
      ggtitle(paste("Network for Module:", input$module_selection,
                    "\nTop", input$top_regulators, "Genes",
                    if(input$selected_group != "All Groups") 
                      paste("\nGroup:", input$selected_group) else ""))
    
    
  })
  
  
  # Heatmap download handler
  output$download_module_heatmap <- downloadHandler(
    filename = function() {
      paste0("module_heatmap_", input$module_selection, "_", 
             if(input$selected_group != "All Groups") paste0(input$selected_group, "_"),
             Sys.Date(), ".tiff")
    },
    content = function(file) {
      # Get dimensions and resolution from UI inputs
      width <- input$heatmap_width
      height <- input$heatmap_height
      res <- input$heatmap_res
      
      # Create TIFF file
      tiff(file, width = width, height = height, units = "in", res = res)
      
      # Prepare data
      reg_data1 <- rv$reg_data[,-dim(rv$reg_data)[2]]
      
      module_genes_top <- rv$moduleMembership %>%
        select(input$module_selection)
      module_genes_top <- data.frame(Genes=rownames(module_genes_top), module_genes_top)
      
      module_genes_top2 <- module_genes_top[order(module_genes_top[,2], decreasing=TRUE),]
      
      topGenes <- module_genes_top2[,1][1:input$top_regulators]
      
      # Filter by selected group if not "All Groups"
      if(input$selected_group != "All Groups") {
        group_indices <- which(rv$group_data$Group == input$selected_group)
        reg_data1 <- reg_data1[group_indices,]
      }
      
      # Set up the gene expression data frame
      mod_mat <- reg_data1 %>%
        dplyr::select(all_of(topGenes)) %>%
        as.matrix()
      
      group_data <- rv$reg_data[rownames(mod_mat), ncol(rv$reg_data), drop = FALSE]
      
      # Prepare final matrix with ordered groups
      mod_mat <- data.frame(mod_mat, Group=group_data[,1])
      mod_mat <- mod_mat[order(mod_mat$Group),]
      
      # Extract the matrix without the group column
      plot_matrix <- as.matrix(mod_mat[, -ncol(mod_mat)])
      
      # Color schemes to match plotly display
      col <- colorRampPalette(rev(RColorBrewer::brewer.pal(n = 10, name = "RdYlBu")))(100)
      ByPal <- colorRampPalette(c('yellow', 'purple'))
      
      # Get unique groups for side colors
      groups <- mod_mat$Group
      unique_groups <- unique(groups)
      group_colors <- ByPal(length(unique_groups))
      names(group_colors) <- unique_groups
      row_colors <- group_colors[groups]
      
      # Create the heatmap
      par(oma=c(2,2,2,2))
      heatmap.2(plot_matrix,
                main = paste("Top", input$top_regulators, "Genes,",
                             "\nModule:", input$module_selection,
                             if(input$selected_group != "All Groups") 
                               paste("\nGroup:", input$selected_group) else ""),
                col = col,
                scale = "column",  # Scale columns like in plotly
                trace = "none",
                density.info = "none",
                key = TRUE,
                keysize = 1,
                margins = c(10, 10),  # Adjust margins for labels
                dendrogram = "none",  # No dendrograms to match plotly
                Rowv = FALSE,         # No row clustering
                Colv = FALSE,         # No column clustering
                srtCol = 90,          # Rotate column labels
                cexRow = 0.8,         # Adjust row label size
                cexCol = 0.8)         # Adjust column label size
      dev.off()
    }
  )
  
  
  
  # Network download handler
  output$download_module_network <- downloadHandler(
    filename = function() {
      paste0("module_network_", input$module_selection, "_",
             if(input$selected_group != "All Groups") paste0(input$selected_group, "_"),
             Sys.Date(), ".tiff")
    },
    content = function(file) {
      # Get dimensions and resolution from UI inputs
      width <- input$network_width
      height <- input$network_height
      res <- input$network_res
      
      # Create TIFF file
      tiff(file, width = width, height = height, units = "in", res = res)
      
      # Prepare data
      reg_data1 <- rv$reg_data[,-dim(rv$reg_data)[2]]
      
      # Filter by selected group if not "All Groups"
      if(input$selected_group != "All Groups") {
        group_indices <- which(rv$group_data$Group == input$selected_group)
        reg_data1 <- reg_data1[group_indices,]
      }
      
      datExpr <- reg_data1[,sapply(reg_data1, is.numeric)]
      datExpr <- datExpr[complete.cases(datExpr), ]
      datExpr <- as.matrix(datExpr)
      
      # Create network
      module_genes <- which(rv$net_colors==input$module_selection)
      gene_names <- colnames(datExpr)[module_genes]
      TOM_module <- rv$TOM[module_genes, module_genes]
      network <- graph_from_adjacency_matrix(TOM_module, mode = "undirected", weighted = TRUE)
      V(network)$name <- gene_names
      V(network)$degree <- degree(network, mode = "all")
      top_genes <- order(V(network)$degree, decreasing = TRUE)[1:input$top_regulators]
      filtered_network <- induced_subgraph(network, vids = top_genes)
      network <- filtered_network
      network_tidy <- as_tbl_graph(network)
      
      # Create plot
      print(ggraph(network_tidy, layout = "fr") +
              geom_edge_link(aes(edge_alpha = weight), show.legend = FALSE) +
              geom_node_point(aes(color = input$module_selection), size = 3) +
              geom_node_text(aes(label = name), size = 4, repel = T, max.overlaps = 50) + 
              scale_color_manual(values = input$module_selection) +            
              theme_void() + 
              theme(legend.position = "none") +
              ggtitle(paste("Network for Module:", input$module_selection,
                            "\nTop", input$top_regulators, "Genes",
                            if(input$selected_group != "All Groups") 
                              paste("\nGroup:", input$selected_group) else "")))
      
      dev.off()
    }
  )
  
  
  
  
  # Gene Set Enrichment Tab
  observeEvent(input$run_enrichment, {
    req(rv$deg_results)
    
    # Extract gene list from DEG results
    gene_list <- rv$deg_results$logFC
    names(gene_list) <- rownames(rv$deg_results)
    gene_list <- sort(gene_list, decreasing = TRUE)
    
    # Run enrichment analysis based on user selection
    withProgress(message = "Running Enrichment Analysis", value = 0, {
      enrichment_results <- NULL
      tryCatch({
        incProgress(0.3, detail = "Running GSEA")
        if (input$enrichment_type == "gobp") {
          enrichment_results <- gseGO(
            geneList = gene_list,
            OrgDb = org.Hs.eg.db,
            ont = "BP",
            keyType = "SYMBOL",
            pvalueCutoff = input$GSEA_pvalue_threshold
          )
        } else if (input$enrichment_type == "gomf") {
          enrichment_results <- gseGO(
            geneList = gene_list,
            OrgDb = org.Hs.eg.db,
            ont = "MF",
            keyType = "SYMBOL",
            pvalueCutoff = input$GSEA_pvalue_threshold
          )
        } else if (input$enrichment_type == "gocc") {
          enrichment_results <- gseGO(
            geneList = gene_list,
            OrgDb = org.Hs.eg.db,
            ont = "CC",
            keyType = "SYMBOL",
            pvalueCutoff = input$GSEA_pvalue_threshold
          )
        } else if (input$enrichment_type == "kegg") {
          enrichment_results <- gseKEGG(
            geneList = gene_list,
            organism = "hsa",
            pvalueCutoff = input$GSEA_pvalue_threshold
          )
        }
      }, error = function(e) {
        shinyalert("Error", "Enrichment analysis failed. Check your input or database selection.")
      })
      incProgress(1, detail = "Complete")
      
      
      # Check if results are empty
      if (is.null(enrichment_results) || nrow(as.data.frame(enrichment_results)) == 0) {
        shinyalert("No Results", "No significant enrichment terms were found. Try adjusting the p-value threshold or using a different gene set.", type = "warning")
        return(NULL)
      }
      
      # Save results and render outputs
      rv$enrichment_results <- enrichment_results
    })
    
    output$enrichment_table <- renderDT({
      req(rv$enrichment_results)
      datatable(as.data.frame(rv$enrichment_results), options = list(scrollX = TRUE,pageLength = 3))
    })
    
    output$enrichment_network <- renderPlot({
      req(rv$enrichment_results)
      goplot(rv$enrichment_results, showCategory  = input$top_enriched, title = "Enrichment Netwrok")
    })
    
    output$enrichment_plot <- renderPlot({
      req(rv$enrichment_results)
      dotplot(rv$enrichment_results, showCategory = input$top_enriched)+
        theme(axis.text.x=element_text(size=15),axis.text.y=element_text(size=15))
    })
    
  })
  
  # Download Handlers
  output$download_enrichment_results <- downloadHandler(
    filename = function() {
      paste("enrichment_results_", Sys.Date(), ".csv", sep = "")
    },
    content = function(file) {
      write.csv(as.data.frame(rv$enrichment_results), file, row.names = FALSE)
    }
  )
  
  
  output$download_enrichment_plot <- downloadHandler(
    filename = function() {
      paste("enrichment_plot_", Sys.Date(), ".tiff", sep = "")
    },
    content = function(file) {
      width <- input$enrichment_plot_width
      height <- input$enrichment_plot_height
      res <- input$enrichment_plot_res
      tiff(file, width = width, height = height, units = "in", res = res)
      P1<-dotplot(rv$enrichment_results, showCategory = input$top_enriched) +
        theme(axis.text.x = element_text(size = 15), axis.text.y = element_text(size = 15))
      print(P1)
      dev.off()
    }
  )  
  
  output$download_enrichment_network <- downloadHandler(
    filename = function() {
      paste("enrichment_network_", Sys.Date(), ".tiff", sep = "")
    },
    content = function(file) {
      width <- input$network_width
      height <- input$network_height
      res <- input$network_res
      tiff(file, width = width, height = height, units = "in", res = res)
      P2<-goplot(rv$enrichment_results, showCategory = input$top_enriched, title = "Enrichment Netwrok")
      print(P2)
      dev.off()
    }
  )  
  
  
  
  
  
  #MRA
  # MRA Tab
  observe({
    req(input$data_type_mra,req(rv$group_data))
    if (input$sidebar == "master_regulators") {
      isolate({
        if (input$data_type_mra == "raw" & !is.null(rv$raw_data) & !input$normalized_data) {
          updateSelectInput(session, "mra_group", choices = unique(rv$group_data$Group))
          rv$mra_data <- rv$raw_data
        } else if (input$data_type_mra == "processed" & !is.null(rv$Normalized_data) & !input$normalized_data) {
          updateSelectInput(session, "mra_group", choices = unique(rv$group_data$Group))
          rv$mra_data <- rv$Normalized_data

        } else if (input$data_type_mra == "normalized" & input$normalized_data & !is.null(rv$raw_data)) {
          updateSelectInput(session, "mra_group", choices = unique(rv$group_data$Group))
          rv$mra_data <- rv$raw_data

        }
        
      })
    }
  })
  
  
  # Update the observeEvent section for MRA
  observeEvent(input$run_mra, {
    req(rv$deg_results)
    req(rv$mra_data)
    req(input$mra_group)
    
    withProgress(message = 'Running Master Regulator Analysis', value = 0, {
      
      incProgress(0.3, detail = "Preparing data matrices")
      # Convert to character
      rv$mra_data$Group <- as.character(rv$mra_data$Group)
      
      # Prepare expression data matrices
      expr_data_mra <- rv$mra_data %>%
        dplyr::filter(Group == input$mra_group) %>%
        dplyr::select(-Group)  # Drop Group column
      
      other_data <- rv$mra_data %>%
        dplyr::filter(Group != input$mra_group) %>%
        dplyr::select(-Group)
      
      # Get significant DEGs
      deg_genes <- rownames(rv$deg_results)[rv$deg_results$adj.P.Val < input$qvalue_threshold & 
                                              abs(rv$deg_results$logFC) > input$logfc_threshold]
      
      if(length(deg_genes) == 0) {
        shinyalert("Warning", "No significant DEGs found. Using top 100 genes by p-value.")
        deg_genes <- rownames(rv$deg_results)[order(rv$deg_results$adj.P.Val)[1:100]]
      }
      
      # Initialize storage for results
      rv$regulons <- NULL
      rv$predicted <- NULL
      
      # Run corto with error handling
      tryCatch({
        
        incProgress(0.2, detail = "Running Corto")
        rv$regulons <- corto(as.matrix(t(expr_data_mra)), 
                             centroids = deg_genes,
                             nbootstraps = 100,
                             verbose = FALSE)
        
        if (!is.null(rv$regulons)) {
          incProgress(0.3, detail = "Preparing Results")
          rv$predicted <- mra(as.matrix(t(expr_data_mra)), 
                              as.matrix(t(other_data)),
                              regulon = rv$regulons)
        }
        
      }, error = function(e) {
        shinyalert("Error", paste("MRA Analysis failed:", e$message))
        return(NULL)
      })
      
      incProgress(1, detail = "Complete")
      
      # Store expression data results
      rv$mra_results <- as.matrix(names(rv$predicted$nes))
    })
    
    # Render table
    output$mra_table <- renderDT({
      req(rv$mra_results)
      datatable(rv$mra_results,
                options = list(scrollX = TRUE,
                               pageLength = 10))
    })
    
    # Create visualization
    output$mra_plot1 <- renderPlot({
      req(rv$predicted)  # Make sure predicted results exist
      mraplot(rv$predicted, 
              title = paste("Master Regulator Analysis using DEG Genes"))
    })
    
  })
  
  # Download handler for MRA results table
  output$download_mra_table <- downloadHandler(
    filename = function() {
      paste("mra_results_", Sys.Date(), ".csv", sep = "")
    },
    content = function(file) {
      write.csv(rv$mra_results, file, row.names = TRUE)
    }
  )
  
  # Download handler for MRA plot
  output$download_mra_plot1 <- downloadHandler(
    filename = function() {
      paste("mra_plot_", Sys.Date(), ".tiff", sep = "")
    },
    content = function(file) {
      width <- input$mra_plot2_width
      height <- input$mra_plot2_height
      res <- input$mra_plot2_res
      
      tiff(file, width = width, height = height, units = "in", res = res)
      mraplot(rv$predicted, 
              title = paste("Master Regulator Analysis using DEG Genes"))
      dev.off()
    })
  
  
  
  
  
  #Search in PubMed & Clinical Trials
  
  debug_log <- reactiveVal("")
  
  
  # Function to add to debug log
  add_debug <- function(msg) {
    current <- debug_log()
    debug_log(paste0(current, "\n", Sys.time(), ": ", msg))
  }
  
  # Safely extract data from PubMed summaries
  safe_extract <- function(summary, field) {
    tryCatch({
      if (field == "title") {
        return(summary$title %||% "Title not available")
      } else if (field == "year") {
        date_str <- summary$pubdate %||% "Year not available"
        if (nchar(date_str) == 4) {
          return(date_str)
        }
        return(sub("^(\\d{4}).*", "\\1", date_str))
      } else if (field == "abstract") {
        return(summary$summary %||% "Abstract not available")
      }
    }, error = function(e) {
      return("Error extracting field")
    })
  }
  
  # Replace `%||%` with a utility function for safer handling
  `%||%` <- function(a, b) if (!is.null(a) && length(a) > 0) a else b
  
  # Reactive value to store results
  search_results <- reactiveVal(data.frame())
  
  
  observe({
    if(input$search_source=="gene_regulators") 
    {
      req(rv$module_colors)
      updateRadioButtons(
        session, 
        "module_selection2", 
        choices = rv$module_colors
      )}
  })
  
  
  
  
  
  # First, add this at the beginning of your server function
  search_results <- reactiveVal()
  
  # Create a function to perform the search (place this before the observers)
  performSearch <- function() {
    # Validate required inputs
    req(input$search_source)
    req(input$disease_filter)
    req(input$search_field)
    
    if (is.null(input$disease_filter) || trimws(input$disease_filter) == "") {
      shinyalert(
        title = "Missing Disease!",
        text = "Please enter at least one disease before searching.",
        type = "warning"
      )
      return()
    }
    
    # Initialize genes variable
    genes <- NULL
    
    # Get genes based on selected source
    if(input$search_source == "deg_analysis") {
      req(rv$deg_results)
      genes <- rownames(rv$deg_results[rv$deg_results$adj.P.Val < input$qvalue_threshold & 
                                         abs(rv$deg_results$logFC) > input$logfc_threshold,]) %>%
        head(20)
    } else if(input$search_source == "master_regulators") {
      req(rv$mra_results)
      genes <- rv$mra_results[,1]  # Assuming first column contains gene names
    } else if(input$search_source == "gene_regulators") {
      req(rv$wgcna_modules)
      req(input$module_selection2)
      
      # Get genes from selected module
      module_genes <- rv$wgcna_modules %>%
        filter(module == input$module_selection2) %>%
        pull(gene) %>%
        head(20)
      genes <- module_genes
    }
    
    # Validate we have genes to search
    if(is.null(genes) || length(genes) == 0) {
      shinyalert(
        title = "Warning",
        text = "No genes available for search. Please check your selection criteria.",
        type = "warning"
      )
      return()
    }
    
    # Parse diseases
    diseases <- unlist(strsplit(input$disease_filter, "\\s*,\\s*"))
    if(length(diseases) == 0) {
      shinyalert(
        title = "Warning",
        text = "Please enter at least one disease term",
        type = "warning"
      )
      return()
    }
    
    # Initialize results dataframe
    results <- data.frame(
      Gene = character(),
      Disease = character(),
      PMID = character(),
      Title = character(),
      Year = character(),
      URL = character(),
      stringsAsFactors = FALSE
    )
    
    withProgress(message = 'Searching PubMed...', value = 0, {
      total_searches <- length(genes) * length(diseases)
      search_count <- 0
      
      for(gene in genes) {
        for(disease in diseases) {
          search_count <- search_count + 1
          incProgress(search_count/total_searches, 
                      detail = sprintf("Searching %s with %s", gene, disease))
          
          # Construct query
          query <- paste0(gene, "[Gene/Protein] AND ", disease, input$search_field)
          
          if (input$clinical_trial_only) {
            query <- paste0("Clinical Trial AND ", query)
          }
          
          tryCatch({
            search_result <- entrez_search(db = "pubmed", 
                                           term = query, 
                                           retmax = 10)
            
            if(length(search_result$ids) > 0) {
              summaries <- entrez_summary(db = "pubmed", 
                                          id = search_result$ids)
              
              for(id in search_result$ids) {
                summary <- summaries[[as.character(id)]]
                
                title <- tryCatch(
                  summary$title,
                  error = function(e) "Title not available"
                )
                
                year <- tryCatch(
                  substr(summary$pubdate, 1, 4),
                  error = function(e) "Year not available"
                )
                
                results <- rbind(results, data.frame(
                  Gene = gene,
                  Disease = disease,
                  PMID = id,
                  Title = title,
                  Year = year,
                  URL = paste0("https://pubmed.ncbi.nlm.nih.gov/", id, "/"),
                  stringsAsFactors = FALSE
                ))
              }
            }
          }, error = function(e) {
            warning(sprintf("Error searching for %s and %s: %s", 
                            gene, disease, e$message))
          })
          
          Sys.sleep(0.1)
        }
      }
    })
    
    # Store results in reactive value
    search_results(results)
    
    # Update table
    updateResultsTable()
  }
  
  # Function to update the results table
  updateResultsTable <- function() {
    results <- search_results()
    
    if(!is.null(results) && nrow(results) > 0) {
      
      results$Year <- as.numeric(results$Year)
      results <- results[order(-results$Year, na.last = TRUE), ]
      
      output$pubmed_table <- renderDT({
        datatable(
          results,
          options = list(
            pageLength = 10,
            scrollX = TRUE,
            dom = 'Bfrtip',
            buttons = c('copy', 'csv'),
            columnDefs = list(list(
              targets = which(colnames(results) == "URL") - 1,
              render = JS("
              function(data, type, row, meta) {
                if (type === 'display') {
                  return '<a href=\"' + data + '\" target=\"_blank\">PubMed</a>';
                }
                return data;
              }
            ")
            ))
          ),
          escape = FALSE,
          selection = 'none',
          rownames = FALSE
        ) %>%
          formatStyle('Gene', backgroundColor = 'lightblue')
      })
    } else {
      output$pubmed_table <- renderDT({
        datatable(
          data.frame(Message = "No results found for the specified search criteria"),
          options = list(dom = 't'),
          rownames = FALSE
        )
      })
    }
  }
  
  # Initialize reactive storage
  rv_ct <- reactiveValues(results = NULL)
  
  # Observer for search button
  observeEvent(input$search, {
    performSearch()
  
  
  # Observer for clinical trials checkbox
  observeEvent(input$clinical_trial_only, {
    # Only perform search if we already have some results
    if (!is.null(search_results())) {
      performSearch()
    }
  })
  
  # Download handler
  output$download_results <- downloadHandler(
    filename = function() {
      paste0("pubmed_results_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv")
    },
    content = function(file) {
      write.csv(search_results(), file, row.names = FALSE)
    }
  )
  
  
  
  # Clinical Trials

    req(input$disease_filter)
    
    if (is.null(input$disease_filter) || trimws(input$disease_filter) == "") {
      shinyalert(
        title = "Missing Disease!",
        text = "Please enter at least one disease before searching.",
        type = "warning"
      )
      return()
    }
    
    
    diseases <- unlist(strsplit(input$disease_filter, "\\s*,\\s*"))
    if (length(diseases) == 0) {
      shinyalert("Warning", "Please enter at least one disease term.", type = "warning")
      return()
    }
    
    # Create empty results dataframe
    results_ct <- data.frame(
      Disease = character(),
      NCT_ID = character(),
      Title = character(),
      Status = character(),
      Phase = character(),
      URL = character(),
      stringsAsFactors = FALSE
    )
    
    withProgress(message = 'Searching ClinicalTrials.gov...', value = 0, {
      for (i in seq_along(diseases)) {
        disease <- diseases[i]
        # Update progress
        incProgress(1/length(diseases), detail = paste("Processing", disease))
        
        query_url <- paste0(
          "https://clinicaltrials.gov/api/v2/studies?query.term=",
          URLencode(disease),
          "&fields=NCTId,BriefTitle,OverallStatus,Phase&format=json&pageSize=50"
        )
        
        response <- tryCatch({
          httr::GET(query_url)
        }, error = function(e) {
          print(paste("API request failed:", e$message))
          return(NULL)
        })
        
        if (!is.null(response) && httr::status_code(response) == 200) {
          api_content <- httr::content(response, as = "parsed")
          
          if (!is.null(api_content$studies) && length(api_content$studies) > 0) {
            # Process each study
            for (study in api_content$studies) {
              new_row <- data.frame(
                Disease = disease,
                NCT_ID = ifelse(!is.null(study$protocolSection$identificationModule$nctId), 
                                study$protocolSection$identificationModule$nctId, 
                                "N/A"),
                Title = ifelse(!is.null(study$protocolSection$identificationModule$briefTitle), 
                               study$protocolSection$identificationModule$briefTitle, 
                               "Title not available"),
                Status = ifelse(!is.null(study$protocolSection$statusModule$overallStatus), 
                                study$protocolSection$statusModule$overallStatus, 
                                "Unknown"),
                Phase = ifelse(!is.null(study$protocolSection$designModule$phases), 
                               paste(study$protocolSection$designModule$phases, collapse = ", "), 
                               "Not Specified"),
                URL = paste0("https://clinicaltrials.gov/study/", 
                             study$protocolSection$identificationModule$nctId),
                stringsAsFactors = FALSE
              )
              results_ct <- rbind(results_ct, new_row)
            }
          }
        }
        Sys.sleep(0.5)  # Prevent rate limits
      }
    })
    
    # Store results in reactive value
    rv_ct$results <- results_ct
    
    # Print debugging information
    print(paste("Number of results:", nrow(results_ct)))
    print("Results stored in rv_ct$results")
  })
  

  # Render Clinical Trials DataTable with Clickable Links
  output$clinical_trials_table <- renderDT({
    req(rv_ct$results)  # Ensure results exist
    
    if (!is.null(rv_ct$results) && nrow(rv_ct$results) > 0) {
      
      # Create a new data frame to avoid modifying reactive values directly
      results_ct_display <- rv_ct$results
      
      # Convert URLs into clickable hyperlinks
      results_ct_display$URL <- paste0(
        "<a href='", results_ct_display$URL, "' target='_blank'>View on ClinicalTrials.gov</a>"
      )
      
      # Render DataTable
      datatable(
        results_ct_display,
        options = list(
          pageLength = 10,
          scrollX = TRUE,
          dom = 'Bfrtip',
          buttons = c('copy', 'csv', 'excel')
        ),
        escape = FALSE,  # This allows the URL column to render as HTML
        selection = 'none'
      )
      
    } else {
      datatable(
        data.frame(Message = "No results found"),
        options = list(dom = 't')
      )
    }
  })
  
  
  # Download handler
  output$download_results1 <- downloadHandler(
    filename = function() {
      paste0("clinical_trials_results_", Sys.Date(), ".csv")
    },
    content = function(file) {
      req(rv_ct$results)
      write.csv(rv_ct$results, file, row.names = FALSE)
    }
  )
  
  
  
  
  
  observe({
    req(input$gene_source_db)
    if(input$gene_source_db=="wgcna") 
    {
      req(rv$module_colors)
      updateRadioButtons(
        session, 
        "module_selection3", 
        choices = rv$module_colors
      )}
  })
  
  
  #Gene Disease Network Analysis
  observeEvent(input$run_disease_network, {
    req(input$gene_source)
    
    
    withProgress(message = 'Running Gene Disease Network Analysis...', value = 0, {
    gene_list <- NULL
    if (input$gene_source == "deg") {
      req(rv$deg_results)
      gene_list <- rownames(rv$deg_results)
    } else if (input$gene_source == "mra") {
      req(rv$mra_results)
      gene_list <- rv$mra_results[,1]
    } else if (input$gene_source == "wgcna") {
      req(rv$wgcna_modules)
      req(input$module_selection3)
      gene_list <-  rv$wgcna_modules %>%
        filter(module == input$module_selection3) %>%
        pull(gene) 
    }
    

    incProgress(0.2, detail = "Mapping gene IDs")
    entrez_ids <- mapIds(org.Hs.eg.db, keys = gene_list, column = "ENTREZID", keytype = "SYMBOL", multiVals = "first")
    disease_enrichment <- enrichDGN(gene = entrez_ids, pvalueCutoff = 0.05, pAdjustMethod = "BH")
    
    incProgress(0.4, detail = "Performing disease enrichment")
    sig_results <- disease_enrichment@result[disease_enrichment@result$p.adjust < input$qvalue_threshold_DEGEN, ]
    sig_results <- sig_results[1:min(input$top_n_results, nrow(sig_results)), ]
    
    incProgress(0.6, detail = "Building network edges")
    edges <- data.frame()
    for (i in 1:nrow(sig_results)) {
      genes_in_term <- strsplit(sig_results$geneID[i], "/")[[1]]
      disease <- sig_results$Description[i]
      
      gene_symbols <- sapply(genes_in_term, function(id) {
        symbol <- names(entrez_ids)[which(entrez_ids == id)]
        if (length(symbol) == 0) return(id)
        return(symbol[1])
      })
      
      edge_data <- data.frame(
        source = gene_symbols,
        target = disease,
        weight = -log10(sig_results$p.adjust[i]),
        stringsAsFactors = FALSE
      )
      edges <- rbind(edges, edge_data)
    }
    
    edges <- na.omit(edges)  # Remove NA entries
    edges <- edges[!duplicated(edges), ]  # Ensure no duplicates
    edges$weight <- as.numeric(edges$weight)
    
    # Ensure edges length matches unique sources/targets
    if (nrow(edges) != length(unique(c(edges$source, edges$target)))) {
      edges <- edges[1:min(nrow(edges), length(unique(c(edges$source, edges$target)))), ]
    }
    
    incProgress(0.8, detail = "Creating graph")
    network <- graph_from_data_frame(edges[1:min(nrow(edges), nrow(edges)), ], directed = FALSE)
    centrality_metrics <- data.frame(
      node = V(network)$name,
      degree = degree(network),
      betweenness = betweenness(network),
      closeness = closeness(network)
    )
    network <- delete.vertices(network, V(network)$name[centrality_metrics$betweenness < input$centrality_filter])
    
    V(network)$type <- ifelse(V(network)$name %in% gene_list, "gene", "disease")
    V(network)$color <- ifelse(V(network)$type == "gene", "#A8D5E2", "#F9A7B0")
    V(network)$size <- ifelse(V(network)$type == "gene", 12, 18)
    V(network)$label.color <- "black"
   
    if (length(E(network)) == nrow(edges)) {
      E(network)$width <- edges$weight / max(edges$weight) * 2
    } else {
      warning("Mismatch between edges and network. Skipping weight assignment.")
    }
    
    
    incProgress(1, detail = "Rendering network")
    output$disease_network_plot <- renderVisNetwork({
      # Convert igraph to visNetwork format
      nodes <- data.frame(id = V(network)$name,
                          label = V(network)$name,
                          color = ifelse(V(network)$type == "gene", "#A8D5E2", "#F9A7B0"),
                          size = ifelse(V(network)$type == "gene", 12, 18))
      
      edges <- data.frame(from = edges$source,
                          to = edges$target,
                          width = edges$weight / max(edges$weight) * 5)  # Scale edge weights
      
      visNetwork(nodes, edges) %>%
        visOptions(highlightNearest = TRUE, nodesIdSelection = TRUE) %>%
        visInteraction(zoomView = TRUE, dragNodes = TRUE, dragView = TRUE) %>%
        visPhysics(enabled = TRUE)
    })
    
    
    output$disease_enrichment_table <- renderDT({
      datatable(disease_enrichment@result, options = list(scrollX = TRUE))
    })
    
    output$download_disease_results <- downloadHandler(
      filename = function() { paste("disease_enrichment_results", Sys.Date(), ".csv", sep = "") },
      content = function(file) { write.csv(disease_enrichment@result, file, row.names = FALSE) }
    )
    
    output$download_disease_network <- downloadHandler(
      filename = function() { paste("disease_network", Sys.Date(), ".tiff", sep = "") },
      content = function(file) {
        width <- input$dn_plot_width
        height <- input$dn_plot_height
        res <- input$dn_plot_res
        tiff(file, width = width, height = height, units = "in", res = res)
        plot(network, layout = layout_with_fr, vertex.label.cex = 0.9, vertex.frame.color = "darkgray", edge.curved = 0.3, edge.color = "gray70")
        dev.off()
      }
    )
  })
  })

  
  
  #Database Search
  observe({
    req(input$gene_source_db)
    if(input$gene_source_db=="gene_regulators") 
    {
      req(rv$module_colors)
      updateRadioButtons(
        session, 
        "module_selection_db", 
        choices = rv$module_colors
      )}
  })
  
   # 
  observeEvent(input$search_database, {
    # Get gene list based on selection
    req(input$gene_source_db)
    
    gene_list <- NULL
    if (input$gene_source_db == "deg_analysis") {
      req(rv$deg_results)
      gene_list <- rownames(rv$deg_results)[1:20]
    } else if (input$gene_source_db == "master_regulators") {
      req(rv$mra_results)
      gene_list <- rv$mra_results[,1]
    } else if (input$gene_source_db == "gene_regulators") {
      req(rv$wgcna_modules)
      req(input$module_selection_db)
      gene_list <- rv$wgcna_modules %>%
        filter(module == input$module_selection_db) %>%
        pull(gene) %>%
        head(20)
    }
    
    gene_query <- as.character(gene_list)
    
    withProgress(message = "Searching Drug Databases...", value = 0, {
      
      # Step 1: Search ChEMBL
      rv$chembl_data <- search_chembl(gene_query)
      
      rv$chembl_data$molecule_pref_name <- paste0(
        '<a href="#" onclick="Shiny.setInputValue(\'selected_molecule\', \'', 
        rv$chembl_data$molecule_pref_name, '\', {priority: \'event\'}); return false;">', 
        rv$chembl_data$molecule_pref_name, '</a>'
      )
      
      # Step 2: Search BindingDB 
      rv$bindingdb_data <- search_bindingdb(gene_query)
      
      if (nrow(rv$bindingdb_data) == 0) {
        showNotification("No BindingDB results found.", type = "warning", duration = 5)
      }
      
      output$chembl_results <- renderDT({
        datatable(rv$chembl_data, escape = FALSE, options = list(
          scrollX = TRUE, 
          pageLength = 5,
          autoWidth = TRUE,
          searchHighlight = TRUE
        )) %>%
          formatStyle("molecule_pref_name", cursor = "pointer", color = "blue", textDecoration = "underline")
      })
      
      output$bindingdb_results <- renderDT({
        datatable(rv$bindingdb_data, options = list(
          scrollX = TRUE, 
          pageLength = 5,
          autoWidth = TRUE,
          searchHighlight = TRUE
        ))
      })
    })      
  })

    # Download Handlers
    output$download_chembl <- downloadHandler(
      filename = function() { paste("ChEMBL_results_", Sys.Date(), ".csv", sep="") },
      content = function(file) { write.csv(rv$chembl_data, file, row.names = FALSE) }
    )
    
    output$download_bindingdb <- downloadHandler(
      filename = function() { paste("BindingDB_results_", Sys.Date(), ".csv", sep="") },
      content = function(file) { write.csv(rv$bindingdb_data, file, row.names = FALSE) }
    )
 
    observeEvent(input$selected_molecule, {
      req(rv$chembl_data)
      selected_row <- rv$chembl_data[rv$chembl_data$molecule_pref_name == input$selected_molecule, ]
      
      pubchem_cid <- selected_row$pubchem_cid
      chembl_img_url <- selected_row$chembl_img_url
      
      if (!is.na(pubchem_cid) && pubchem_cid != "") {
        viewer_url <- paste0("https://pubchem.ncbi.nlm.nih.gov/compound/", pubchem_cid, "#section=3D-Conformer")
        output$molecule_3d_viewer <- renderUI({
          tags$iframe(src = viewer_url, width = "100%", height = "600px", style = "border: none;")
        })
      } else {
        output$molecule_3d_viewer <- renderUI({
          tags$img(src = chembl_img_url, width = "400px", height = "400px")
        })
      }
    })
    
}

shinyApp(ui, server)

}

APIomics()


