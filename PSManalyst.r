## PSM analyst dashboard for FragPipe search results
## The files psm.tsv and protein.tsv are the inputs for the viewwer
## It is possible to filter the PSMs by the hyperscore

# Check if the required R libraries are installed and install them if necessary.
CRAN_packages <- c("shiny", "shinydashboard", "tidyverse", "janitor", "ggseqlogo", "ggtext", "lsa", "vegan", "plotly", "viridis", "ggfortify")
not_installed_CRAN <- CRAN_packages[!(CRAN_packages %in% installed.packages()[ , "Package"])]
if(length(not_installed_CRAN)) install.packages(not_installed_CRAN)

GitHub_packages <- c("ggpointdensity", "wordcloud2")
not_installed_GitHub <- GitHub_packages[!(GitHub_packages %in% installed.packages()[ , "Package"])]
if(length(not_installed_GitHub)) install.packages(not_installed_GitHub)

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
library(vegan)            # from CRAN

# Increase the maximum file size to 1000 MB
options(shiny.maxRequestSize = 1000 * 1024^2)

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

color_blue_seq <- c("#d4e6f1", "#a9cce3", "#7fb3d5", "#5499c7", "#2980b9", "#1f618d", "#154360")

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
              value = 15, step = 5),
      sliderInput("probability",
                    label = "PeptideProphet Probability",
                    min = 0, max = 1,
                    value = 0.95, step = 0.05),
      textInput("protein_pettern",
                label = "Remove an organism by entry name",
                value = "",
                placeholder = "HUMAN"),
      checkboxInput("case_sensitive",
                label = "Case sensitive",
                value = FALSE),
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
                  box(title = "Word cloud of peptide sequences", status = "primary", solidHeader = TRUE, wordcloud2Output("plot16"), collapsible = TRUE),
                  box(title = "N-termini SeqLogo", status = "primary", solidHeader = TRUE, plotOutput("plot3"), collapsible = TRUE),
                  box(title = "C-termini SeqLogo", status = "primary", solidHeader = TRUE, plotOutput("plot4"), collapsible = TRUE),
                  box(title = "m/z over retention time", status = "primary", solidHeader = TRUE, plotOutput("plot6"), collapsible = TRUE),
                  box(title = "Mass error (ppm)", status = "primary", solidHeader = TRUE, plotOutput("plot8"), collapsible = TRUE),
                  box(title = "Peptide length", status = "primary", solidHeader = TRUE, plotOutput("plot2"), collapsible = TRUE),
                  box(title = "Charge state distribution", status = "primary", solidHeader = TRUE, plotOutput("plot5"), collapsible = TRUE),
                  box(title = "Number of missed cleavages", status = "primary", solidHeader = TRUE, plotOutput("plot7"), collapsible = TRUE),
                  box(title = "Uniqueness", status = "primary", solidHeader = TRUE, plotOutput("plot14"), collapsible = TRUE),
                  box(title = "Hyperscore distribution", status = "primary", solidHeader = TRUE, plotOutput("plot9"), collapsible = TRUE),
                  box(title = "Next Score distribution", status = "primary", solidHeader = TRUE, plotOutput("plot10"), collapsible = TRUE),
                  box(title = "PeptideProphet probability", status = "primary", solidHeader = TRUE, plotOutput("plot11"), collapsible = TRUE),
                  box(title = "Expectation (PeptideProphet)", status = "primary", solidHeader = TRUE, plotOutput("plot12"), collapsible = TRUE),
                  box(title = "Assigned modifications", status = "primary", solidHeader = TRUE, plotOutput("plot13"), collapsible = TRUE),
                  box(title = "Top 20 proteins with more PSMs", status = "primary", solidHeader = TRUE, plotOutput("plot15"), collapsible = TRUE)
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
  filter_text <- paste("Showing PSMs with Hyperscore ≥", input$hyperscore,
                       "and PeptideProphet probability ≥", input$probability)
  
  if (!is.null(input$protein_pattern) && input$protein_pattern != "") {
    pattern_text <- if(inpt$case_sensitive) {
      paste("and organism entry name matching:", input$protein_pattern, "(case sensitive)")
    }
    else {
      paste("and organism entry name matching:", input$protein_pattern, "(case insensitive)")
    }
    filter_text <- paste(filter_text, pattern_text)
  }
    
    infoBox("Filter settings",
            filter_text,
            icon = icon("info"),
            color = "black"
    )
  })

  # Import and pre-process the uploaded psm.tsv file
  data <- reactive({
    req(input$psm)
  # Read the psm.tsv file and filter based on hyperscore and PeptideProphet probability
    psm_file <- readr::read_tsv(input$psm$datapath) %>%
      janitor::clean_names() %>%
      dplyr::filter(.$hyperscore >= input$hyperscore & .$probability >= input$probability) 
  # Filter by organism entry name if provided
    if (input$protein_pettern != "") {
      if (input$case_sensitive) {
        psm_file <- psm_file %>%
          dplyr::filter(str_detect(entry_name, input$protein_pettern, negate = TRUE))
      } else {
        psm_file <- psm_file %>%
          dplyr::filter(str_detect(tolower(entry_name), tolower(input$protein_pettern), negate = TRUE))
      }
    }
      
    psm_file <- psm_file %>%
      dplyr::mutate(
        fingerprint_Nterm = case_when(
            str_detect(extended_peptide, "^\\.") ~ "NA",
            TRUE ~ substr(extended_peptide, 2, 16)
            ),
        fingerprint_Cterm = substr(extended_peptide, nchar(extended_peptide) - 15, nchar(extended_peptide) - 2),
        fingerprint_Nterm = str_extract(fingerprint_Nterm, ".{4}\\..{4}"),
        fingerprint_Nterm = str_remove_all(fingerprint_Nterm, "\\."),
        fingerprint_Cterm = str_extract(fingerprint_Cterm, ".{4}\\..{4}"),
        fingerprint_Cterm = str_remove_all(fingerprint_Cterm, "\\."),
        delta_mass_ppm = (observed_m_z-calculated_m_z)/calculated_m_z*1e6
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
    scale_fill_gradient(low = "#d4e6f1", high = "#154360") +
    geom_vline(xintercept = 4.5, color = "black", linetype = "dashed") +
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
        legend.key.height = unit(0.25, "cm"),
        legend.title.position = "top")
  })

  output$plot2 <- renderPlot({
    data() %>%
    as.data.frame() %>%
    ggplot() +
    geom_density(aes(x = peptide_length),
        fill = "#5499c7") +
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
  scale_x_continuous(breaks = c(1, 2, 3, 4, 5, 6, 7, 8),
                     labels = c("1" = "P4", "2" = "P3", "3" = "P2", "4" = "P1", "5" = "P1'", "6" = "P2'", "7" =  "P3'", "8" = "P4'")) +
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
  scale_x_continuous(breaks = c(1, 2, 3, 4, 5, 6, 7, 8),
                     labels = c("1" = "P4", "2" = "P3", "3" = "P2", "4" = "P1", "5" = "P1'", "6" = "P2'", "7" =  "P3'", "8" = "P4'")) +
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
      fill = "#5499c7", color = "black") +
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
        y = "Scan range (m/z)",
        color = "Number of Neighborhoods") +
    theme(
    legend.position = "bottom",
        legend.key.width = unit(1.5, "cm"),
        legend.key.height = unit(0.25, "cm")
    )
  })

  output$plot7 <- renderPlot({
    data() %>%
    as.data.frame() %>%
    dplyr::count(number_of_missed_cleavages) %>%
    dplyr::mutate(number_of_missed_cleavages = factor(number_of_missed_cleavages)) %>%
    ggplot(aes(x = number_of_missed_cleavages, y = n)) +
    geom_bar(stat = "identity", position = "dodge", show.legend = FALSE, 
        fill = "#5499c7", color = "black") +
    geom_text(aes(label = n), vjust = -0.5, size = 5) +
    labs(x = "Number of Missed Cleavages",
        y = "Count")
  })

 output$plot8 <- renderPlot({
    data() %>%
    as.data.frame() %>%
      dplyr::filter(abs(delta_mass_ppm) < 100) %>% 
      ggplot(aes(x = retention/60,
                 y = delta_mass_ppm)
      ) +
      geom_point(alpha = 0.1, color = "black", size = 1) +
      geom_hline(yintercept = c(10, 0, -10), color = "red", linetype = "dashed", linewidth = 0.2) +
      labs(title = "Mass error in ppm",
           x = "Retention time (min)",
           y = "Mass error (ppm)",
          caption = "ppm error is calculated as ∆m/z over theoretical m/z * 1e6")
  })

  output$plot9 <- renderPlot({
    data() %>%
    as.data.frame() %>%
    ggplot() +
    geom_histogram(aes(x = hyperscore), 
        fill = "#5499c7", color = "black") +
    labs(x = "Hyperscore",
        y = "Count",
        caption = "Similarity score between observed and theoretical spectra, higher values indicate greater similarity")
  })

  output$plot10 <- renderPlot({
    data() %>%
    as.data.frame() %>%
    ggplot() +
    geom_histogram(aes(x = nextscore), 
        fill = "#5499c7", color = "black") +
    labs(x = "Nextscore",
        y = "Count",
        caption = "Similarity score (hyperscore) of the second-highest scoring match for the spectrum")
  })

  output$plot11 <- renderPlot({
    data() %>%
    as.data.frame() %>%
    ggplot() +
    geom_histogram(aes(x = probability), 
        fill = "#5499c7", color = "black") +
    labs(x = "PeptideProphet Probability",
        y = "Count",
        caption = "Confidence score determined by PeptideProphet, higher values indicate greater confidence")
  })

  output$plot12 <- renderPlot({
    data() %>%
    as.data.frame() %>%
    ggplot() +
    geom_histogram(aes(x = expectation), 
        fill = "#5499c7", color = "black") +
    labs(x = "Expectation value",
        y = "Count",
        caption = "Expectation value from statistical modeling with PeptideProphet, lower values indicate higher likelihood")
  })

  output$plot13 <- renderPlot({
    data() %>%
    as.data.frame() %>%
    tidyr::separate_rows(assigned_modifications, sep = ",") %>%
    dplyr::mutate(assigned_modifications = str_remove_all(assigned_modifications, ".*\\(|\\)"),
                  assigned_modifications = ifelse(is.na(assigned_modifications), 
                                                    "Unassigned modifications", 
                                                    assigned_modifications)) %>%
    dplyr::count(assigned_modifications) %>%
    ggplot(aes(y = assigned_modifications, x = n)) +
      geom_col(fill = "#5499c7", color = "black") +
      geom_text(aes(label = n), hjust = -0.1, size = 5) +
      labs(y = "Assigned Modifications",
           x = "Count",
           caption = "Number of modifications assigned to the peptide sequences")
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
        fill = "#5499c7", color = "black") +
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
        fill = "#5499c7", color = "black", stat = "identity") +
    labs(x = "Number of PSMs",
        y = "Protein")
  })

  output$plot16 <- renderWordcloud2({
    data() %>%
    as.data.frame() %>%
    dplyr::count(peptide) %>%
    dplyr::mutate(frequency = round(n / sum(n) * 100, 2)) %>%
    dplyr::select(-n) %>%
    wordcloud2::wordcloud2(color = rep_len(color_blue_seq, nrow(.)),
        backgroundColor = "white",
        size = 1,
        shuffle = TRUE,
        minRotation = -pi/6,
        maxRotation = pi/6,
        widgetsize = "100%")
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
        fill = "#5499c7", color = "black") +
    labs(x = "Protein coverage (%)",
        y = "Count")
  })

  output$plot18 <- renderPlot({
    protein_data() %>%
    as.data.frame() %>%
    dplyr::count(organism) %>%
    ggplot() +
    geom_bar(aes(x = n, y = reorder(organism, n)),
        fill = "#5499c7", color = "black", stat = "identity") +
    geom_text(aes(x = n, y = reorder(organism, n), label = n), hjust = -0.1, size = 5) +
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
        levels = c("Experimental evidence at protein level", 
                   "Experimental evidence at transcript level", 
                   "Protein inferred from homology", 
                   "Protein predicted"))
      ) %>%
    ggplot() +
    geom_bar(aes(y = protein_existence,
                 fill = protein_existence),
        color = "black", show.legend = FALSE) +
    scale_fill_manual(values = c("#5499c7", "#7fb3d5", "#a9cce3", "#d4e6f1")) +
    geom_text(aes(y = protein_existence, 
                  label = ..count..), show.legend = FALSE,
              stat = "count", vjust = -0.5, size = 7, fontface = "bold") +
    labs(y = NULL,
        x = "Count",
        fill = NULL)
  })

  output$plot20 <- renderPlot({
    protein_data() %>%
    as.data.frame() %>%
    ggplot() +
    geom_histogram(aes(x = protein_probability),
        fill = "#5499c7", color = "black") +
    labs(x = "Protein Probability",
        y = "Count")
  })

  output$plot21 <- renderPlot({
    protein_data() %>%
    as.data.frame() %>%
    ggplot() +
    geom_histogram(aes(x = top_peptide_probability),
        fill = "#5499c7", color = "black") +
    labs(x = "Peptide Probability",
        y = "Count",
        caption = "Best peptide probability of supporting peptides")
  })

  output$plot22 <- renderPlot({
    protein_data() %>%
    as.data.frame() %>%
    ggplot() +
    geom_histogram(aes(x = total_peptides),
        fill = "#5499c7", color = "black") +
    labs(x = "Total peptides mapped to proteins",
        y = "Count")
  })

  output$plot23 <- renderPlot({
    protein_data() %>%
    as.data.frame() %>%
    ggplot() +
    geom_histogram(aes(x = razor_spectral_count),
        fill = "#5499c7", color = "black") +
    labs(x = "Razor Spectral Count",
        y = "Count",
        caption = "Number of PSMs corresponding to the razor peptides")
  })

  output$plot24 <- renderPlot({
    protein_data() %>%
    as.data.frame() %>%
    ggplot() +
    geom_histogram(aes(x = razor_intensity),
        fill = "#5499c7", color = "black") +
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
        fill = "#5499c7", color = "black", stat = "identity") +
    labs(x = "log2 of Razor Intensity",
        y = "Protein")
  })

# Import and pre-process the uploaded combined_protein.tsv file
  combined_protein_data <- reactive({
    req(input$combined_protein)
    combined_protein_file <- readr::read_tsv(input$combined_protein$datapath) %>%
      janitor::clean_names() %>%
      dplyr::select(protein_id, ends_with("max_lfq_intensity")) %>%
      column_to_rownames("protein_id") %>%
      dplyr::rename_all(~str_remove(., "_max_lfq_intensity")) %>%
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
    geom_violin(aes(x = sample, y = maxlfq_intensity),
        fill = "#5499c7", alpha = 0.7, color = "black") +
    geom_boxplot(aes(x = sample, y = maxlfq_intensity),
        fill = "white", outliers = FALSE,
        color = "black", width = 0.1, show.legend = FALSE) +
    labs(x = NULL,
        y = "log2(MaxLFQ intensity)") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
  })

output$plot27 <- renderPlotly({
    combined_protein_data() %>%
    as.data.frame() %>%
    ggplot(aes(x = !!sym(input$xcol), y = !!sym(input$ycol))) +
    geom_point(alpha = 0.7, show.legend = FALSE) +
    geom_smooth(method = "lm", se = FALSE,
        color = "#5499c7") +
    labs(x = paste0("log2(", input$xcol, ")"),
        y = paste0("log2(", input$ycol, ")"))
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
        legend.key.width = unit(2.5, "cm"),
        legend.key.height = unit(0.25, "cm")
         ) +
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
        legend.key.width = unit(2.5, "cm"),
        legend.key.height = unit(0.25, "cm")) +
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
        legend.key.width = unit(2.5, "cm"),
        legend.key.height = unit(0.25, "cm")) +
    labs(x = NULL,
        y = NULL,
        fill = "Jaccard similarity")
})

}

# Run the application
shinyApp(ui, server)
