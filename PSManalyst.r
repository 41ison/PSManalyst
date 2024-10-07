## PSM analyst dashboard for FragPipe search results
## The files psm.tsv and protein.tsv are the inputs for the viewwer
## It is possible to filter the PSMs by the hyperscore

# Load required libraries
library(shiny)            # from CRAN
library(shinydashboard)   # from CRAN
library(tidyverse)        # from CRAN
library(janitor)          # from CRAN
library(ggseqlogo)        # from CRAN
library(ggpointdensity)   # from GitHub
library(wordcloud2)       # from GitHub
library(ggtext)           # from CRAN
library(lsa)              # from CRAN
library(plotly)           # from CRAN
library(viridis)          # from CRAN
library(ggfortify)        # from CRAN

# Increase the maximum file size to 200 MB
options(shiny.maxRequestSize = 200 * 1024^2)

# set the general theme for the plots
theme_set(theme_bw())
theme_update(
    text = element_text(color = "black", size = 15),
    axis.text = element_text(color = "black"),
    axis.title = element_text(color = "black", face = "bold"),
    legend.title = element_text(face = "bold", hjust = 0.5),
    legend.title.position = "top"
  )

# to calculate the frequency of each amino acid in each column in percentage
aa_freq <- function(x) {
    table(x) / length(x) * 100
}

# to impute missing amino acids with zero and reorder if necessary
complete_and_reorder_amino_acids <- function(element) {
# List of 20 amino acids
twenty_amino_acids <- c('A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y')

# Complete missing amino acids with zero
  for (amino_acid in twenty_amino_acids) {
    if (!amino_acid %in% names(element)) {
      element[[amino_acid]] <- 0
    }
  }

# Reorder amino acids
  element <- element[match(twenty_amino_acids, names(element))]
  
  return(element)
}

# function to extract the matrix of amino acid frequencies
extract_matrix <- function(data) {
    fingerprint_protease <- c(data$fingerprint_Nterm,
            data$fingerprint_Cterm) %>%
        na.omit() %>%
        strsplit("")
 
# create a matrix with the list of peptide sequences
  mat_aa <- matrix(unlist(fingerprint_protease),
                        ncol = 8, byrow = TRUE)

# remove the rows containing "B" or any other unwanted amino acids in the matrix
 mat_aa <- mat_aa[!apply(mat_aa, 1,
                    function(x) any(x == "B|X|Y")), ]

# calculate the frequency of each amino acid in each column and plot a heatmap
  mat_aa_freq <- apply(mat_aa, 2, aa_freq)

  new_list <- lapply(mat_aa_freq, complete_and_reorder_amino_acids)

  final_matrix <- matrix(unlist(new_list,), ncol = 8, byrow = FALSE)

  colnames(final_matrix) <- c("P4", "P3", "P2", "P1", "P1'", "P2'", "P3'", "P4'")
  row.names(final_matrix) <- c("A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y")

  return(final_matrix)
}

# Define UI for application that reads a psm.tsv file and generates a PICS map report dashboard
ui <- dashboardPage(

  dashboardHeader(
      title = "PSM Analyst for FragPipe",
      titleWidth = "250",
      dropdownMenu(type = "messages",
                   messageItem(
                     from = "Support",
                     message = "felipealison@gmail.com",
                     icon = icon("envelope")
                   )
      )
  ),

  dashboardSidebar(
    sidebarMenu(
      menuItem("PSM viewer",
              tabName = "psm",
              icon = icon("barcode",
              lib = "glyphicon")),
      fileInput(inputId = "psm",
              label = "Choose the psm.tsv file",
               accept = ".tsv"),
      sliderInput("hyperscore",
              label = "PSM hyperscore filter",
              min = 0, max = 1000,
              value = 20, step = 10),
      menuItem("Protein viewer",
              tabName = "protein",
              icon = icon("equalizer",
              lib = "glyphicon")),
      fileInput(inputId = "protein",
              label = "Choose the protein.tsv file",
              accept = ".tsv"),
      fileInput(inputId = "combined_protein",
              label = "Choose the combined_protein.tsv file",
              accept = ".tsv"),
      selectInput("xcol", 
              "X Sample",
              choices = NULL),
      selectInput("ycol",
              "Y Sample",
              choices = NULL)
    )
  ),

  dashboardBody(
    tabItems(
      tabItem(tabName = "psm",
              fluidRow(
                  infoBoxOutput("info_box1", width = 12),
                  box(title = "Protease fingerprint", status = "primary", solidHeader = TRUE, plotOutput("plot1"), collapsible = TRUE),
                  box(title = "Peptide length", status = "primary", solidHeader = TRUE, plotOutput("plot2"), collapsible = TRUE),
                  box(title = "N-termini SeqLogo", status = "primary", solidHeader = TRUE, plotOutput("plot3"), collapsible = TRUE),
                  box(title = "C-termini SeqLogo", status = "primary", solidHeader = TRUE, plotOutput("plot4"), collapsible = TRUE),
                  box(title = "Charge state distribution", status = "primary", solidHeader = TRUE, plotOutput("plot5"), collapsible = TRUE),
                  box(title = "m/z over retention time", status = "primary", solidHeader = TRUE, plotOutput("plot6"), collapsible = TRUE),
                  box(title = "Number of missed cleavages", status = "primary", solidHeader = TRUE, plotOutput("plot7"), collapsible = TRUE),
                  box(title = "Number of enzymatic termini", status = "primary", solidHeader = TRUE, plotOutput("plot8"), collapsible = TRUE),
                  box(title = "Hyperscore distribution", status = "primary", solidHeader = TRUE, plotOutput("plot9"), collapsible = TRUE),
                  box(title = "Next Score distribution", status = "primary", solidHeader = TRUE, plotOutput("plot10"), collapsible = TRUE),
                  box(title = "PeptideProphet probability", status = "primary", solidHeader = TRUE, plotOutput("plot11"), collapsible = TRUE),
                  box(title = "Expectation (PeptideProphet)", status = "primary", solidHeader = TRUE, plotOutput("plot12"), collapsible = TRUE),
                  box(title = "Purity (Philosopher Freequant)", status = "primary", solidHeader = TRUE, plotOutput("plot13"), collapsible = TRUE),
                  box(title = "Uniqueness", status = "primary", solidHeader = TRUE, plotOutput("plot14"), collapsible = TRUE),
                  box(title = "Top 20 proteins with more PSMs", status = "primary", solidHeader = TRUE, plotOutput("plot15"), collapsible = TRUE),
                  box(title = "Word cloud of peptide sequences", status = "primary", solidHeader = TRUE, wordcloud2Output("plot16"), collapsible = TRUE)
      )
    ),

      tabItem(tabName = "protein",
            fluidRow(
                  infoBoxOutput("info_box2", width = 12),
                  box(title = "Protein coverage", status = "primary", solidHeader = TRUE, plotOutput("plot17"), collapsible = TRUE),
                  box(title = "Number of proteins by organim", status = "primary", solidHeader = TRUE, plotOutput("plot18"), collapsible = TRUE),
                  box(title = "Protein existence evidence", status = "primary", solidHeader = TRUE, plotOutput("plot19"), collapsible = TRUE),
                  box(title = "Protein probability (ProteinProphet)", status = "primary", solidHeader = TRUE, plotOutput("plot20"), collapsible = TRUE),
                  box(title = "Top Peptide Probability", status = "primary", solidHeader = TRUE, plotOutput("plot21"), collapsible = TRUE),
                  box(title = "Total peptides mapped to the proteins", status = "primary", solidHeader = TRUE, plotOutput("plot22"), collapsible = TRUE),
                  box(title = "Razor spectral count", status = "primary", solidHeader = TRUE, plotOutput("plot23"), collapsible = TRUE),
                  box(title = "Razor intensity", status = "primary", solidHeader = TRUE, plotOutput("plot24"), collapsible = TRUE),
                  box(title = "Top 20 proteins with higher razor intensity", status = "primary", solidHeader = TRUE, plotOutput("plot25"), collapsible = TRUE),
                  box(title = "MaxLFQ intensity distribution", status = "primary", solidHeader = TRUE, plotOutput("plot26"), collapsible = TRUE),
                  box(title = "Sample correlation - Non-normalized log2(Intensity)", status = "primary", height = 600, solidHeader = TRUE, plotlyOutput("plot27"), collapsible = FALSE),
                tabBox(
                  title = "Similarity metrics", side = "right", height = 600,
                  tabPanel("Cosine similarity", plotOutput("cosine_similarity")),
                  tabPanel("Euclidean distance", plotOutput("euclidean_distance")),
                  tabPanel("Jaccard similarity", plotOutput("jaccard_similarity"))
                )
            )
          )
      )
  )
)

# Define server logic required to read the psm.tsv file and generate the PICS map report
server <- function(input, output, session) {

# Information box to display the hyperscore filter
output$info_box1 <- renderInfoBox({
    infoBox("Filter the PSMs by hyperscore to select the best matches",
            paste("Showing PSMs with Hyperscore ≥ ", input$hyperscore),
            icon = icon("info"),
            color = "black"
    )
  })

  # Import and pre-process the uploaded psm.tsv file
  data <- reactive({
    req(input$psm)
    psm_file <- readr::read_tsv(input$psm$datapath) %>%
      janitor::clean_names() %>%
      dplyr::filter(.$hyperscore >= input$hyperscore) %>%
      dplyr::mutate(
        fingerprint_Nterm = case_when(
            str_detect(extended_peptide, "^\\.") ~ "NA",
            TRUE ~ substr(extended_peptide, 2, 16)
            ),
        fingerprint_Cterm = substr(extended_peptide, nchar(extended_peptide) - 15, nchar(extended_peptide) - 2),
        fingerprint_Nterm = str_extract(fingerprint_Nterm, ".{4}\\..{4}"),
        fingerprint_Nterm = str_remove_all(fingerprint_Nterm, "\\."),
        fingerprint_Cterm = str_extract(fingerprint_Cterm, ".{4}\\..{4}"),
        fingerprint_Cterm = str_remove_all(fingerprint_Cterm, "\\.")
        ) %>%
      dplyr::relocate(extended_peptide, .before = fingerprint_Nterm)
  })

# Extract the matrix of amino acid frequencies
frequency_matrix_of_aa <- reactive({
  req(data())
  extract_matrix(data())
})

  # Render plots for the PSM viewer
  output$plot1 <- renderPlot({
    frequency_matrix_of_aa() %>%
    as.data.frame() %>%
    rownames_to_column(var = "residue") %>%
    pivot_longer(cols = -residue, names_to = "position",
                    values_to = "frequency") %>%
    dplyr::mutate(
        position = factor(position, c("P4", "P3", "P2", "P1", "P1'", "P2'", "P3'", "P4'")),
        residue = factor(residue, c("A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y"))
        ) %>%
    ggplot(aes(x = position, y = residue, fill = frequency)) +
    geom_tile(color = "black") +
    scale_fill_gradient(low = "grey90", high = "dodgerblue4") +
    theme_void() +
    labs(
        title = "Cleavage Site Specificity",
        x = "Position",
        y = "Amino acid residue",
        fill = "Frequency (%)"
    ) +
    theme(text = element_text(size = 15, color = "black"),
        axis.text.x = element_text(hjust = 0.5),
        axis.text.y = element_text(hjust = 0.5),
        legend.title = element_text(hjust = 0.5),
        plot.title = element_text(hjust = 0.5),
        legend.position = "bottom",
        legend.key.width = unit(1.5, "cm"),
        legend.title.position = "top")
  })

  output$plot2 <- renderPlot({
    data() %>%
    as.data.frame() %>%
    ggplot() +
    geom_density(aes(x = peptide_length),
        fill = "dodgerblue4", alpha = 0.7) +
    labs(x = "Peptide Length",
        y = "Frequency (%)") +
    theme(text = element_text(size = 15, color = "black"),
        axis.text.x = element_text(hjust = 0.5),
        plot.title = element_text(hjust = 0.5),)
  })

  output$plot3 <- renderPlot({
    data() %>%
    as.data.frame() %>%
    dplyr::select(fingerprint_Nterm) %>%
    na.omit() %>%
    ggseqlogo::ggseqlogo(
  method = "bits",
  seq_type = "AA"
  ) +
  geom_hline(yintercept = 0, 
        color = "black", linetype = "dashed") +
    geom_vline(xintercept = 4.5, 
        color = "black", linetype = "dashed") +
  theme_bw() +
  theme(plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
    text = element_text(size = 15, color = "black"),
    legend.position = "bottom",
    legend.title.position = "top",
    legend.title = element_text(size = 12, hjust = 0.5)
  ) +
  labs(title = "SeqLogo of the N-termini fingerprint",
       x = "Amino acid position",
       y = "Bits")
  })

  output$plot4 <- renderPlot({
    data() %>%
    as.data.frame() %>%
    dplyr::select(fingerprint_Cterm) %>%
    na.omit() %>%
    ggseqlogo::ggseqlogo(
  method = "bits",
  seq_type = "AA"
  ) +
  geom_hline(yintercept = 0, 
        color = "black", linetype = "dashed") +
    geom_vline(xintercept = 4.5, 
        color = "black", linetype = "dashed") +
  theme_bw() +
  theme(plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
    text = element_text(size = 15, color = "black"),
    legend.position = "bottom",
    legend.title.position = "top",
    legend.title = element_text(size = 12, hjust = 0.5)
  ) +
  labs(title = "SeqLogo of the C-termini fingerprint",
       x = "Amino acid position",
       y = "Bits")
  })

  output$plot5 <- renderPlot({
    data() %>%
    as.data.frame() %>%
    ggplot() +
    geom_bar(aes(x = charge), 
      fill = "dodgerblue4", alpha = 0.7, color = "black") +
    labs(x = "Charge state",
        y = "Count")
  })

  output$plot6 <- renderPlot({
    data() %>%
    as.data.frame() %>%
    ggplot(aes(x = retention / 60, y = observed_m_z)) +
    ggpointdensity::geom_pointdensity(size = 0.25) +
    viridis::scale_color_viridis(option = "plasma") +
    labs(x = "Retention time (min)",
        y = "Scan range (m/z)") +
    theme(
    legend.position = "bottom",
        legend.key.width = unit(1.5, "cm"))
  })

  output$plot7 <- renderPlot({
    data() %>%
    as.data.frame() %>%
    ggplot() +
    geom_bar(aes(x = number_of_missed_cleavages), 
        fill = "dodgerblue4", alpha = 0.7, color = "black") +
    labs(x = "Number of Missed Cleavages",
        y = "Count",
        caption = "Number of potential enzymatic cleavage sites within the identified sequence")
  })

  output$plot8 <- renderPlot({
    data() %>%
    as.data.frame() %>%
    ggplot() +
    geom_bar(aes(x = number_of_enzymatic_termini), stat = "count",
        fill = "dodgerblue4", alpha = 0.7, color = "black") +
    labs(x = "Number of Enzymatic Termini",
        y = "Count",
        caption = "2 = fully-enzymatic, 1 = semi-enzymatic, 0 = non-enzymatic")
  })

  output$plot9 <- renderPlot({
    data() %>%
    as.data.frame() %>%
    ggplot() +
    geom_histogram(aes(x = hyperscore), 
        fill = "dodgerblue4", alpha = 0.7, color = "black") +
    labs(x = "Hyperscore",
        y = "Count",
        caption = "Similarity score between observed and theoretical spectra, higher values indicate greater similarity")
  })

  output$plot10 <- renderPlot({
    data() %>%
    as.data.frame() %>%
    ggplot() +
    geom_histogram(aes(x = nextscore), 
        fill = "dodgerblue4", alpha = 0.7, color = "black") +
    labs(x = "Nextscore",
        y = "Count",
        caption = "Similarity score (hyperscore) of the second-highest scoring match for the spectrum")
  })

  output$plot11 <- renderPlot({
    data() %>%
    as.data.frame() %>%
    ggplot() +
    geom_histogram(aes(x = probability), 
        fill = "dodgerblue4", alpha = 0.7, color = "black") +
    labs(x = "PeptideProphet Probability",
        y = "Count",
        caption = "Confidence score determined by PeptideProphet, higher values indicate greater confidence")
  })

  output$plot12 <- renderPlot({
    data() %>%
    as.data.frame() %>%
    ggplot() +
    geom_histogram(aes(x = expectation), 
        fill = "dodgerblue4", alpha = 0.7, color = "black") +
    labs(x = "Expectation value",
        y = "Count",
        caption = "Expectation value from statistical modeling with PeptideProphet, lower values indicate higher likelihood")
  })

  output$plot13 <- renderPlot({
    data() %>%
    as.data.frame() %>%
    ggplot() +
    geom_histogram(aes(x = purity), 
        fill = "dodgerblue4", alpha = 0.7, color = "black") +
    labs(x = "Purity",
        y = "Count",
        caption = "Proportion of total ion abundance in the inclusion window from the precursor (including precursor isotopic peaks)")
  })

  output$plot14 <- renderPlot({
    data() %>%
    dplyr::count(is_unique) %>%
    as.data.frame() %>%
    dplyr::mutate(uniqueness = case_when(
        is_unique == TRUE ~ "Unique",
        TRUE ~ "Shared"
    )) %>%
    ggplot(aes(x = uniqueness, y = n)) +
    geom_bar(stat = "identity", position = "dodge", show.legend = FALSE, 
        fill = "dodgerblue4", alpha = 0.7, color = "black") +
    geom_text(aes(label = n), vjust = -0.5, size = 5) +
    labs(x = "Unique peptides",
        y = "Count")
  })

  output$plot15 <- renderPlot({
    data() %>%
    as.data.frame() %>%
    dplyr::group_by(entry_name) %>%
    dplyr::summarize(n_psm = n()) %>%
    dplyr::arrange(desc(n_psm)) %>%
    dplyr::mutate(entry_name = factor(entry_name,
        levels = entry_name
    )) %>%
    head(20) %>%
    ggplot() +
    geom_bar(aes(x = n_psm, y = reorder(entry_name, n_psm)), 
        fill = "dodgerblue4", alpha = 0.7, color = "black", stat = "identity") +
    labs(x = "Number of PSMs",
        y = "Protein")
  })

  output$plot16 <- renderWordcloud2({
    data() %>%
    as.data.frame() %>%
    dplyr::count(peptide) %>%
    wordcloud2::wordcloud2(color = "random-dark",
        backgroundColor = "white",
        size = 1,
        shuffle = TRUE)
  })

# Import and pre-process the uploaded protein.tsv file
  protein_data <- reactive({
    req(input$protein)
    protein_file <- readr::read_tsv(input$protein$datapath) %>%
      janitor::clean_names()
  })

# Information box to display the hyperscore filter
  output$info_box2 <- renderInfoBox({
    infoBox("protein.tsv files contain FDR-filtered protein results, where each row is an identified protein group",
            icon = icon("info"),
            color = "black"
    )
  })
  
# Render plots for the protein viewer
  output$plot17 <- renderPlot({
    protein_data() %>%
    as.data.frame() %>%
    ggplot() +
    geom_histogram(aes(x = coverage), 
        fill = "dodgerblue4", alpha = 0.7, color = "black") +
    labs(x = "Protein coverage (%)",
        y = "Count")
  })

  output$plot18 <- renderPlot({
    protein_data() %>%
    as.data.frame() %>%
    dplyr::count(organism) %>%
    ggplot() +
    geom_bar(aes(x = n, y = reorder(organism, n)),
        fill = "dodgerblue4", alpha = 0.7, color = "black", stat = "identity") +
    geom_text(aes(x = n, y = reorder(organism, n), label = n), vjust = -0.5, size = 5) +
    labs(x = "Number of proteins",
        y = NULL) +
    theme(axis.text.y = element_text(face = "italic"))
  })

  output$plot19 <- renderPlot({
    protein_data() %>%
    as.data.frame() %>%
     dplyr::mutate(
      protein_existence = str_remove(protein_existence, ".*\\:"),
      protein_existence = factor(protein_existence,
        levels = c("Experimental evidence at protein level", "Experimental evidence at transcript level", "Protein inferred from homology", "Protein predicted"))
      ) %>%
    ggplot() +
    geom_bar(aes(y = protein_existence),
        fill = "dodgerblue4", alpha = 0.7, color = "black") +
    geom_text(aes(y = protein_existence, label = ..count..), stat = "count", vjust = -0.5, size = 5) +
    labs(y = NULL,
        x = "Count")
  })

  output$plot20 <- renderPlot({
    protein_data() %>%
    as.data.frame() %>%
    ggplot() +
    geom_histogram(aes(x = protein_probability),
        fill = "dodgerblue4", alpha = 0.7, color = "black") +
    labs(x = "Protein Probability",
        y = "Count")
  })

  output$plot21 <- renderPlot({
    protein_data() %>%
    as.data.frame() %>%
    ggplot() +
    geom_histogram(aes(x = top_peptide_probability),
        fill = "dodgerblue4", alpha = 0.7, color = "black") +
    labs(x = "Peptide Probability",
        y = "Count",
        caption = "Best peptide probability of supporting peptides")
  })

  output$plot22 <- renderPlot({
    protein_data() %>%
    as.data.frame() %>%
    ggplot() +
    geom_histogram(aes(x = total_peptides),
        fill = "dodgerblue4", alpha = 0.7, color = "black") +
    labs(x = "Total peptides mapped to proteins",
        y = "Count")
  })

  output$plot23 <- renderPlot({
    protein_data() %>%
    as.data.frame() %>%
    ggplot() +
    geom_histogram(aes(x = razor_spectral_count),
        fill = "dodgerblue4", alpha = 0.7, color = "black") +
    labs(x = "Razor Spectral Count",
        y = "Count",
        caption = "Number of PSMs corresponding to the razor peptides")
  })

  output$plot24 <- renderPlot({
    protein_data() %>%
    as.data.frame() %>%
    ggplot() +
    geom_histogram(aes(x = razor_intensity),
        fill = "dodgerblue4", alpha = 0.7, color = "black") +
    labs(x = "Razor Intensity",
        y = "Count",
        caption = "Protein intensity calculated using the unique peptides (from the top-N algorithm)")
  })

  output$plot25 <- renderPlot({
    protein_data() %>%
    as.data.frame() %>%
    dplyr::arrange(desc(razor_intensity)) %>%
    head(20) %>%
    ggplot() +
    geom_bar(aes(x = log2(razor_intensity), y = reorder(entry_name, razor_intensity)),
        fill = "dodgerblue4", alpha = 0.7, color = "black", stat = "identity") +
    labs(x = "log2 of Razor Intensity",
        y = "Protein")
  })

# Import and pre-process the uploaded combined_protein.tsv file
  combined_protein_data <- reactive({
    req(input$combined_protein)
    combined_protein_file <- readr::read_tsv(input$combined_protein$datapath) %>%
      janitor::clean_names() %>%
      dplyr::select(protein_id, ends_with("max_lfq_intensity")) %>%
      dplyr::mutate(sample = str_remove(sample, "_max_lfq_intensity")) %>%
      column_to_rownames("protein_id") %>%
      log2()
  })

  # Observe the uploaded file and update selectInput choices
  observe({
    req(combined_protein_data())
    colnames <- colnames(combined_protein_data())
    updateSelectInput(session, "xcol", choices = colnames)
    updateSelectInput(session, "ycol", choices = colnames)
  })

  output$plot26 <- renderPlot({
    combined_protein_data() %>%
    as.data.frame() %>%
    rownames_to_column(var = "protein_id") %>%
    tidyr::pivot_longer(
      cols = -protein_id,
      names_to = "sample",
      values_to = "maxlfq_intensity"
      ) %>%
    ggplot() +
    geom_violin(aes(x = sample, y = log2(maxlfq_intensity)),
        fill = "dodgerblue4", alpha = 0.7, color = "black") +
    geom_boxplot(aes(x = sample, y = log2(maxlfq_intensity)),
        fill = "white", alpha = 0.7, outliers = FALSE,
        color = "black", width = 0.1, show.legend = FALSE) +
    labs(x = NULL,
        y = "Log2(MaxLFQ intensity)") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
  })

output$plot27 <- renderPlotly({
    combined_protein_data() %>%
    as.data.frame() %>%
    ggplot(aes(x = !!sym(input$xcol), y = !!sym(input$ycol))) +
    geom_point(alpha = 0.7, show.legend = FALSE) +
    geom_smooth(method = "lm", se = FALSE,
        color = "darkblue") +
    labs(x = paste0("Log2(", input$xcol, ")"),
        y = paste0("Log2(", input$ycol, ")"))
  })

# calculate the cosine similarity in the matrix and plot the heatmap
output$cosine_similarity <- renderPlot({
    combined_protein_data() %>%
    as.matrix() %>%
    na.omit() %>%
    lsa::cosine() %>%
    as.data.frame() %>%
    rownames_to_column(var = "Sample") %>%
    pivot_longer(-Sample, names_to = "Match", values_to = "value") %>%
    dplyr::mutate(Similarity = "Cosine similarity") %>%
    ggplot() +
    geom_tile(aes(x = Sample, y = Match, fill = value)) +
    viridis::scale_fill_viridis(option = "E") +
    theme(text = element_text(size = 15),
        axis.text.x = element_text(angle = 90,
                        hjust = 1, vjust = 0.5),
        axis.text.y = element_text(angle = 0,
                        hjust = 1, vjust = 0.5),
        legend.position = "bottom",
        legend.key.width = unit(2.5, "cm")) +
    labs(x = NULL,
        y = NULL,
        fill = "Cosine similarity")
})

# calculate the euclidean distance in the matrix and plot the heatmap
output$euclidean_distance <- renderPlot({
    combined_protein_data() %>%
    t() %>%
    dist(method = "euclidean") %>%
    as.matrix() %>%
    as.data.frame() %>%
    rownames_to_column(var = "Sample") %>%
    pivot_longer(-Sample, names_to = "Match", values_to = "value") %>%
    dplyr::mutate(Similarity = "Euclidean distance") %>%
    ggplot() +
    geom_tile(aes(x = Sample, y = Match, fill = value)) +
    viridis::scale_fill_viridis(option = "E") +
    theme(text = element_text(size = 15),
    axis.text.x = element_text(angle = 90,
                        hjust = 1, vjust = 0.5),
        axis.text.y = element_text(angle = 0,
                        hjust = 1, vjust = 0.5),
        legend.position = "bottom",
        legend.key.width = unit(2.5, "cm")) +
    labs(x = NULL,
        y = NULL,
        fill = "Euclidean distance")
        })

# calculate the Jaccard similarity in the matrix and plot the heatmap
output$jaccard_similarity <- renderPlot({
    combined_protein_data() %>%
    t() %>%
    vegan::vegdist(method = "jaccard", na.rm = TRUE) %>%
    as.matrix() %>%
    as.data.frame(as.table(.)) %>% 
    dplyr::mutate(Sample = colnames(.)) %>%
    pivot_longer(-Sample, names_to = "Match", values_to = "value") %>%
    dplyr::mutate(Similarity = "Jaccard similarity") %>%
    ggplot() +
    geom_tile(aes(x = Sample, y = Match, fill = value)) +
    viridis::scale_fill_viridis(option = "E") +
    theme(text = element_text(size = 15),
    axis.text.x = element_text(angle = 90,
                        hjust = 1, vjust = 0.5),
        axis.text.y = element_text(angle = 0,
                        hjust = 1, vjust = 0.5),
        legend.position = "bottom",
        legend.key.width = unit(2.5, "cm")) +
    labs(x = NULL,
        y = NULL,
        fill = "Jaccard similarity")
})

}

# Run the application
shinyApp(ui = ui, server = server)
