## PSM analyst dashboard for FragPipe search results
## The input are psm.tsv, protein.tsv and combined_protein.tsv files
## It is possible to filter the PSMs by the hyperscore and PeptideProphet probability, as well as for enzymatic specificity
## You can remove a contaminant organism as well
## You can customize the color of most of the plots and download all the PSM plots in high resolution

CRAN_packages <- c(
  "shiny",
  "shinydashboard",
  "tidyverse",
  "janitor",
  "ggseqlogo",
  "ggtext",
  "lsa",
  "vegan",
  "plotly",
  "viridis",
  "ggfortify",
  "colourpicker"
)
not_installed_CRAN <- CRAN_packages[
  !(CRAN_packages %in% installed.packages()[, "Package"])
]
if (length(not_installed_CRAN)) {
  install.packages(not_installed_CRAN)
}

GitHub_packages <- c("ggpointdensity", "wordcloud2")
not_installed_GitHub <- GitHub_packages[
  !(GitHub_packages %in% installed.packages()[, "Package"])
]
if (length(not_installed_GitHub)) {
  install.packages(not_installed_GitHub)
}

library(shiny)
library(shinydashboard)
library(tidyverse)
library(janitor)
library(ggseqlogo)
library(ggpointdensity)
library(wordcloud2)
library(ggtext)
library(lsa)
library(plotly)
library(viridis)
library(ggfortify)
library(vegan)
library(colourpicker)

# Increase the maximum file size to 1000 MB
options(shiny.maxRequestSize = 1000 * 1024^2)

theme_set(theme_bw())
theme_update(
  text = element_text(color = "black", size = 18),
  axis.text = element_text(color = "black"),
  axis.title = element_text(color = "black", face = "bold"),
  legend.title = element_text(face = "bold", hjust = 0.5),
  legend.title.position = "top"
)

aa_freq <- function(x) {
  table(x) / length(x) * 100
}

complete_and_reorder_amino_acids <- function(element) {
  twenty_amino_acids <- c(
    'A',
    'C',
    'D',
    'E',
    'F',
    'G',
    'H',
    'I',
    'K',
    'L',
    'M',
    'N',
    'P',
    'Q',
    'R',
    'S',
    'T',
    'V',
    'W',
    'Y'
  )

  for (amino_acid in twenty_amino_acids) {
    if (!amino_acid %in% names(element)) {
      element[[amino_acid]] <- 0
    }
  }
  element <- element[match(twenty_amino_acids, names(element))]

  return(element)
}

# Extract the matrix of amino acid frequencies
extract_matrix <- function(data) {
  fingerprint_protease <- c(data$fingerprint_Nterm, data$fingerprint_Cterm) %>%
    na.omit() %>%
    as.character() %>%
    strsplit("")
  mat_aa <- matrix(unlist(fingerprint_protease), ncol = 8, byrow = TRUE)
  mat_aa <- mat_aa[
    !apply(mat_aa, 1, function(x) any(x %in% c("B", "X", "Z", "U"))),
  ]
  mat_aa_freq <- apply(mat_aa, 2, aa_freq)
  new_list <- lapply(mat_aa_freq, complete_and_reorder_amino_acids)
  final_matrix <- matrix(unlist(new_list), ncol = 8, byrow = FALSE)

  colnames(final_matrix) <- c(
    "P4",
    "P3",
    "P2",
    "P1",
    "P1'",
    "P2'",
    "P3'",
    "P4'"
  )
  rownames(final_matrix) <- c(
    "A",
    "C",
    "D",
    "E",
    "F",
    "G",
    "H",
    "I",
    "K",
    "L",
    "M",
    "N",
    "P",
    "Q",
    "R",
    "S",
    "T",
    "V",
    "W",
    "Y"
  )

  return(final_matrix)
}

# GRAVY (Grand Average of Hydropathy)
# Ref: Kyte J, Doolittle RF. A simple method for displaying the hydropathic character of a protein. J Mol Biol. 1982 May 5;157(1):105-32. doi: 10.1016/0022-2836(82)90515-0
GRAVY <- function(sequence) {
  hydropathy_index <- c(
    A = 1.8,
    R = -4.5,
    N = -3.5,
    D = -3.5,
    C = 2.5,
    Q = -3.5,
    E = -3.5,
    G = -0.4,
    H = -3.2,
    I = 4.5,
    L = 3.8,
    K = -3.9,
    M = 1.9,
    F = 2.8,
    P = -1.6,
    S = -0.8,
    T = -0.7,
    W = -0.9,
    Y = -1.3,
    V = 4.2
  )
  scores <- sapply(strsplit(sequence, NULL)[[1]], function(aa) {
    hydropathy_index[aa]
  })
  return(mean(scores, na.rm = TRUE))
}

calculate_pI <- function(sequence) {
  pKa_values <- c(
    A = 2.34,
    R = 12.48,
    N = 10.76,
    D = 3.86,
    C = 8.33,
    Q = 10.76,
    E = 4.25,
    G = 2.34,
    H = 6.00,
    I = 6.04,
    L = 6.04,
    K = 9.74,
    M = 5.74,
    F = 5.48,
    P = 1.99,
    S = 2.21,
    T = 2.15,
    W = 9.39,
    Y = 10.07,
    V = 6.02
  )

  pI <- mean(
    sapply(strsplit(sequence, NULL)[[1]], function(aa) pKa_values[aa]),
    na.rm = TRUE
  )

  return(pI)
}

color_blue_seq <- c(
  "#d4e6f1",
  "#a9cce3",
  "#7fb3d5",
  "#5499c7",
  "#2980b9",
  "#1f618d",
  "#154360"
)

ui <- dashboardPage(
  dashboardHeader(
    title = "PSM Analyst for FragPipe",
    titleWidth = "250",
    dropdownMenu(
      type = "messages",
      messageItem(
        from = "Support",
        message = "felipealison@gmail.com",
        icon = icon("envelope")
      )
    )
  ),

  dashboardSidebar(
    sidebarMenu(
      menuItem(
        "PSM viewer",
        tabName = "psm",
        icon = icon("barcode", lib = "glyphicon")
      ),
      fileInput(
        inputId = "psm",
        label = "Choose the psm.tsv file",
        accept = ".tsv"
      ),
      sliderInput(
        "hyperscore",
        label = "PSM hyperscore filter",
        min = 0,
        max = 1000,
        value = 0,
        step = 5
      ),
      sliderInput(
        "probability",
        label = "PeptideProphet Probability",
        min = 0,
        max = 1,
        value = 0.95,
        step = 0.01
      ),
      selectInput(
        "specificity_filter",
        label = "Proteolysis fingerprinting specificity",
        choices = c(
          "All" = "all",
          "Fully specific" = "fully_specific",
          "Semi-specific at N-termini" = "semi_n_termini",
          "Semi-specific at C-termini" = "semi_c_termini",
          "Fully semi-specific" = "fully_semi_specific"
        ),
        selected = "all"
      ),
      textInput(
        "protein_pattern",
        label = "Remove an organism by entry name",
        value = "",
        placeholder = "HUMAN"
      ),
      checkboxInput("case_sensitive", label = "Case sensitive", value = FALSE),
      menuItem(
        "Protein viewer",
        tabName = "protein",
        icon = icon("equalizer", lib = "glyphicon")
      ),
      fileInput(
        inputId = "protein",
        label = "Choose the protein.tsv file",
        accept = ".tsv"
      ),
      fileInput(
        inputId = "combined_protein",
        label = "Choose the combined_protein.tsv file",
        accept = ".tsv"
      ),
      selectInput("xcol", "X Sample", choices = NULL),
      selectInput("ycol", "Y Sample", choices = NULL),
      colourpicker::colourInput(
        inputId = "plot_color",
        label = "Select plot color",
        value = "#5499c7"
      ),
      div(
        style = "text-align: center; margin-top: 10px;",
        downloadButton(
          outputId = "download_all_plots",
          label = "Download PSM plots",
          class = "butt"
        )
      ),
      tags$head(tags$style(".butt{background:grey;} .butt{color: #337ab7;}"))
    )
  ),

  dashboardBody(
    tabItems(
      tabItem(
        tabName = "psm",
        fluidRow(
          infoBoxOutput("info_box1", width = 12),
          box(
            title = "Protease fingerprint",
            status = "primary",
            solidHeader = TRUE,
            plotOutput("plot01"),
            collapsible = TRUE
          ),
          box(
            title = "Word cloud of peptide sequences",
            status = "primary",
            solidHeader = TRUE,
            wordcloud2Output("plot02"),
            collapsible = TRUE
          ),
          box(
            title = "N-termini SeqLogo",
            status = "primary",
            solidHeader = TRUE,
            plotOutput("plot03"),
            collapsible = TRUE
          ),
          box(
            title = "C-termini SeqLogo",
            status = "primary",
            solidHeader = TRUE,
            plotOutput("plot04"),
            collapsible = TRUE
          ),
          box(
            title = "m/z over retention time",
            status = "primary",
            solidHeader = TRUE,
            plotOutput("plot05"),
            collapsible = TRUE
          ),
          box(
            title = "Mass error (ppm)",
            status = "primary",
            solidHeader = TRUE,
            plotOutput("plot06"),
            collapsible = TRUE
          ),
          box(
            title = "Peptide length",
            status = "primary",
            solidHeader = TRUE,
            plotOutput("plot07"),
            collapsible = TRUE
          ),
          box(
            title = "GRAVY (Grand Average of Hydropathy)",
            status = "primary",
            solidHeader = TRUE,
            plotOutput("plot08"),
            collapsible = TRUE
          ),
          box(
            title = "Isoelectric Point (pI)",
            status = "primary",
            solidHeader = TRUE,
            plotOutput("plot09"),
            collapsible = TRUE
          ),
          box(
            title = "Charge state distribution",
            status = "primary",
            solidHeader = TRUE,
            plotOutput("plot10"),
            collapsible = TRUE
          ),
          box(
            title = "Number of missed cleavages",
            status = "primary",
            solidHeader = TRUE,
            plotOutput("plot11"),
            collapsible = TRUE
          ),
          box(
            title = "Uniqueness",
            status = "primary",
            solidHeader = TRUE,
            plotOutput("plot12"),
            collapsible = TRUE
          ),
          box(
            title = "Hyperscore distribution",
            status = "primary",
            solidHeader = TRUE,
            plotOutput("plot13"),
            collapsible = TRUE
          ),
          box(
            title = "Next Score distribution",
            status = "primary",
            solidHeader = TRUE,
            plotOutput("plot14"),
            collapsible = TRUE
          ),
          box(
            title = "PeptideProphet probability",
            status = "primary",
            solidHeader = TRUE,
            plotOutput("plot15"),
            collapsible = TRUE
          ),
          box(
            title = "Expectation (PeptideProphet)",
            status = "primary",
            solidHeader = TRUE,
            plotOutput("plot16"),
            collapsible = TRUE
          ),
          box(
            title = "Assigned modifications",
            status = "primary",
            solidHeader = TRUE,
            plotOutput("plot17"),
            collapsible = TRUE
          ),
          box(
            title = "Top 20 proteins with more PSMs",
            status = "primary",
            solidHeader = TRUE,
            plotOutput("plot18"),
            collapsible = TRUE
          )
        )
      ),

      tabItem(
        tabName = "protein",
        fluidRow(
          infoBoxOutput("info_box2", width = 12),
          box(
            title = "Protein coverage",
            status = "primary",
            solidHeader = TRUE,
            plotOutput("plot19"),
            collapsible = TRUE
          ),
          box(
            title = "Number of proteins by organim",
            status = "primary",
            solidHeader = TRUE,
            plotOutput("plot20"),
            collapsible = TRUE
          ),
          box(
            title = "Protein existence evidence",
            status = "primary",
            solidHeader = TRUE,
            plotOutput("plot21"),
            collapsible = TRUE
          ),
          box(
            title = "Protein probability (ProteinProphet)",
            status = "primary",
            solidHeader = TRUE,
            plotOutput("plot22"),
            collapsible = TRUE
          ),
          box(
            title = "Top Peptide Probability",
            status = "primary",
            solidHeader = TRUE,
            plotOutput("plot23"),
            collapsible = TRUE
          ),
          box(
            title = "Total peptides mapped to the proteins",
            status = "primary",
            solidHeader = TRUE,
            plotOutput("plot24"),
            collapsible = TRUE
          ),
          box(
            title = "Razor spectral count",
            status = "primary",
            solidHeader = TRUE,
            plotOutput("plot25"),
            collapsible = TRUE
          ),
          box(
            title = "Razor intensity",
            status = "primary",
            solidHeader = TRUE,
            plotOutput("plot26"),
            collapsible = TRUE
          ),
          box(
            title = "Top 20 proteins with higher razor intensity",
            status = "primary",
            solidHeader = TRUE,
            plotOutput("plot27"),
            collapsible = TRUE
          ),
          box(
            title = "MaxLFQ intensity distribution",
            status = "primary",
            solidHeader = TRUE,
            plotOutput("plot28"),
            collapsible = TRUE
          ),
          box(
            title = "Sample correlation - Non-normalized log2(Intensity)",
            status = "primary",
            height = 600,
            solidHeader = TRUE,
            plotlyOutput("plot29"),
            collapsible = FALSE
          ),
          tabBox(
            title = "Similarity metrics",
            side = "right",
            height = 600,
            tabPanel("Cosine similarity", plotOutput("cosine_similarity")),
            tabPanel("Euclidean distance", plotOutput("euclidean_distance")),
            tabPanel("Jaccard similarity", plotOutput("jaccard_similarity"))
          )
        )
      )
    )
  )
)

server <- function(input, output, session) {
  output$info_box1 <- renderInfoBox({
    filter_text <- paste(
      "Showing PSMs with Hyperscore ≥",
      input$hyperscore,
      "and PeptideProphet probability ≥",
      input$probability
    )

    if (input$specificity_filter != "all") {
      specificity_labels <- c(
        "fully_specific" = "Fully specific",
        "semi_n_termini" = "Semi-specific at N-termini",
        "semi_c_termini" = "Semi-specific at C-termini",
        "fully_semi_specific" = "Fully semi-specific"
      )
      filter_text <- paste(
        filter_text,
        "and",
        specificity_labels[input$specificity_filter],
        "peptides"
      )
    }

    infoBox(
      "Filter settings",
      filter_text,
      icon = icon("info"),
      color = "black"
    )
  })

  # Import and pre-process the uploaded psm.tsv file
  data <- reactive({
    req(input$psm)
    psm_file <- readr::read_tsv(input$psm$datapath) %>%
      janitor::clean_names() %>%
      dplyr::filter(
        .$hyperscore >= input$hyperscore & .$probability >= input$probability
      )

    if (input$protein_pattern != "") {
      if (input$case_sensitive) {
        psm_file <- psm_file %>%
          dplyr::filter(str_detect(
            entry_name,
            input$protein_pattern,
            negate = TRUE
          ))
      } else {
        psm_file <- psm_file %>%
          dplyr::filter(str_detect(
            tolower(entry_name),
            tolower(input$protein_pattern),
            negate = TRUE
          ))
      }
    }

    psm_file <- psm_file %>%
      dplyr::mutate(
        fingerprint_Nterm = case_when(
          str_detect(extended_peptide, "^\\.") ~ "NA",
          TRUE ~ substr(extended_peptide, 2, 16)
        ),
        fingerprint_Cterm = substr(
          extended_peptide,
          nchar(extended_peptide) - 15,
          nchar(extended_peptide) - 2
        ),
        fingerprint_Nterm = str_extract(fingerprint_Nterm, ".{4}\\..{4}"),
        fingerprint_Nterm = str_remove_all(fingerprint_Nterm, "\\."),
        fingerprint_Cterm = str_extract(fingerprint_Cterm, ".{4}\\..{4}"),
        fingerprint_Cterm = str_remove_all(fingerprint_Cterm, "\\."),
        delta_mass_ppm = (observed_m_z - calculated_m_z) / calculated_m_z * 1e6,
        gravy = sapply(peptide, GRAVY),
        isoelectric_point = sapply(peptide, calculate_pI),
        specificity = case_when(
          prev_aa %in% c("K", "R") & str_sub(peptide, -1) %in% c("K", "R") ~
            "Fully specific",
          prev_aa %in% c("K", "R") & !str_sub(peptide, -1) %in% c("K", "R") ~
            "Semi-specific at C-termini",
          !prev_aa %in% c("K", "R") & str_sub(peptide, -1) %in% c("K", "R") ~
            "Semi-specific at N-termini",
          !prev_aa %in% c("K", "R") & !str_sub(peptide, -1) %in% c("K", "R") ~
            "Fully semi-specific",
          TRUE ~ "Unknown"
        )
      ) %>%
      dplyr::relocate(extended_peptide, .before = fingerprint_Nterm) %>%
      dplyr::relocate(specificity, .after = peptide)

    if (input$specificity_filter != "all") {
      specificity_map <- c(
        "fully_specific" = "Fully specific",
        "semi_n_termini" = "Semi-specific at N-termini",
        "semi_c_termini" = "Semi-specific at C-termini",
        "fully_semi_specific" = "Fully semi-specific"
      )
      psm_file <- psm_file %>%
        dplyr::filter(
          specificity == specificity_map[input$specificity_filter]
        )
    }

    return(psm_file)
  })

  frequency_matrix_of_aa <- reactive({
    req(data())
    extract_matrix(data())
  })

  # Render plots for the PSM viewer
  output$plot01 <- renderPlot({
    frequency_matrix_of_aa() %>%
      as.data.frame() %>%
      rownames_to_column(var = "residue") %>%
      pivot_longer(
        cols = -residue,
        names_to = "position",
        values_to = "frequency"
      ) %>%
      dplyr::mutate(
        position = factor(
          position,
          c("P4", "P3", "P2", "P1", "P1'", "P2'", "P3'", "P4'")
        ),
        residue = factor(
          residue,
          c(
            "A",
            "C",
            "D",
            "E",
            "F",
            "G",
            "H",
            "I",
            "K",
            "L",
            "M",
            "N",
            "P",
            "Q",
            "R",
            "S",
            "T",
            "V",
            "W",
            "Y"
          )
        )
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
      theme(
        plot.title = element_text(size = 18, hjust = 0.5, face = "bold"),
        axis.text.x = element_text(hjust = 0.5),
        axis.text.y = element_text(hjust = 0.5),
        legend.title = element_text(hjust = 0.5, face = "bold"),
        legend.position = "bottom",
        legend.key.width = unit(2, "cm"),
        legend.key.height = unit(0.25, "cm"),
        legend.title.position = "top"
      )
  })

  output$plot02 <- renderWordcloud2({
    data() %>%
      as.data.frame() %>%
      dplyr::count(peptide) %>%
      dplyr::mutate(frequency = round(n / sum(n) * 100, 2)) %>%
      dplyr::select(-n) %>%
      wordcloud2::wordcloud2(
        color = rep_len(color_blue_seq, nrow(.)),
        backgroundColor = "white",
        size = 1,
        shuffle = TRUE,
        minRotation = -pi / 6,
        maxRotation = pi / 6,
        widgetsize = "100%"
      )
  })

  output$plot03 <- renderPlot({
    data() %>%
      as.data.frame() %>%
      dplyr::select(fingerprint_Nterm) %>%
      na.omit() %>%
      ggseqlogo::ggseqlogo(
        method = "bits",
        seq_type = "AA"
      ) +
      geom_hline(yintercept = 0, color = "black", linetype = "dashed") +
      geom_vline(xintercept = 4.5, color = "black", linetype = "dashed") +
      scale_x_continuous(
        breaks = c(1, 2, 3, 4, 5, 6, 7, 8),
        labels = c(
          "1" = "P4",
          "2" = "P3",
          "3" = "P2",
          "4" = "P1",
          "5" = "P1'",
          "6" = "P2'",
          "7" = "P3'",
          "8" = "P4'"
        )
      ) +
      labs(
        title = "SeqLogo of the N-termini fingerprint",
        x = "Amino acid position",
        y = "Bits"
      ) +
      theme_bw() +
      theme(
        plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
        text = element_text(size = 15, color = "black"),
        axis.title = element_text(face = "bold"),
        legend.position = "bottom",
        legend.title.position = "top",
        legend.title = element_text(size = 15, hjust = 0.5)
      )
  })

  output$plot04 <- renderPlot({
    data() %>%
      as.data.frame() %>%
      dplyr::select(fingerprint_Cterm) %>%
      na.omit() %>%
      ggseqlogo::ggseqlogo(
        method = "bits",
        seq_type = "AA"
      ) +
      geom_hline(yintercept = 0, color = "black", linetype = "dashed") +
      geom_vline(xintercept = 4.5, color = "black", linetype = "dashed") +
      scale_x_continuous(
        breaks = c(1, 2, 3, 4, 5, 6, 7, 8),
        labels = c(
          "1" = "P4",
          "2" = "P3",
          "3" = "P2",
          "4" = "P1",
          "5" = "P1'",
          "6" = "P2'",
          "7" = "P3'",
          "8" = "P4'"
        )
      ) +
      labs(
        title = "SeqLogo of the C-termini fingerprint",
        x = "Amino acid position",
        y = "Bits"
      ) +
      theme_bw() +
      theme(
        plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
        text = element_text(size = 15, color = "black"),
        axis.title = element_text(face = "bold"),
        legend.position = "bottom",
        legend.title.position = "top",
        legend.title = element_text(size = 15, hjust = 0.5)
      )
  })

  output$plot05 <- renderPlot({
    data() %>%
      as.data.frame() %>%
      ggplot(aes(x = retention / 60, y = observed_m_z)) +
      ggpointdensity::geom_pointdensity(size = 0.25) +
      viridis::scale_color_viridis(option = "plasma") +
      labs(
        x = "Retention time (min)",
        y = "Scan range (m/z)",
        color = "Number of Neighborhoods"
      ) +
      theme(
        legend.position = "bottom",
        legend.key.width = unit(2, "cm"),
        legend.key.height = unit(0.25, "cm")
      )
  })

  output$plot06 <- renderPlot({
    data() %>%
      as.data.frame() %>%
      dplyr::filter(abs(delta_mass_ppm) < 100) %>%
      ggplot(aes(x = retention / 60, y = delta_mass_ppm)) +
      geom_point(alpha = 0.1, color = "black", size = 1) +
      geom_hline(
        yintercept = c(10, 0, -10),
        color = "red",
        linetype = "dashed",
        linewidth = 0.2
      ) +
      labs(
        x = "Retention time (min)",
        y = "Mass error (ppm)",
        caption = "ppm error is calculated as:\n∆m/z over theoretical m/z * 1e6"
      )
  })

  output$plot07 <- renderPlot({
    data() %>%
      as.data.frame() %>%
      ggplot() +
      geom_density(aes(x = peptide_length), fill = input$plot_color) +
      labs(x = "Peptide Length", y = "Frequency (%)") +
      theme(
        text = element_text(size = 15, color = "black"),
        axis.text.x = element_text(hjust = 0.5),
        plot.title = element_text(hjust = 0.5),
      )
  })

  output$plot08 <- renderPlot({
    data() %>%
      as.data.frame() %>%
      ggplot(aes(x = gravy, fill = stat(x))) +
      geom_histogram(color = "black") +
      labs(
        x = NULL,
        y = "Count",
        caption = "GRAVY is a measure of the hydropathic character of a sequence"
      ) +
      scale_fill_viridis_c(name = "GRAVY index", option = "C") +
      theme(
        text = element_text(size = 15, color = "black"),
        axis.text.x = element_text(hjust = 0.5),
        legend.position = "bottom",
        legend.key.width = unit(2.5, "cm"),
        legend.key.height = unit(0.25, "cm")
      )
  })

  output$plot09 <- renderPlot({
    data() %>%
      as.data.frame() %>%
      ggplot(aes(x = isoelectric_point, fill = stat(x))) +
      geom_histogram(color = "black") +
      labs(x = NULL, y = "Count") +
      scale_fill_viridis_c(name = "Isoelectric Point (pI)", option = "C") +
      theme(
        text = element_text(size = 15, color = "black"),
        axis.text.x = element_text(hjust = 0.5),
        legend.position = "bottom",
        legend.key.width = unit(2.5, "cm"),
        legend.key.height = unit(0.25, "cm")
      )
  })

  output$plot10 <- renderPlot({
    data() %>%
      as.data.frame() %>%
      ggplot() +
      geom_bar(aes(x = charge), fill = input$plot_color, color = "black") +
      labs(x = "Charge state", y = "Count")
  })

  output$plot11 <- renderPlot({
    data() %>%
      as.data.frame() %>%
      dplyr::count(number_of_missed_cleavages) %>%
      dplyr::mutate(
        number_of_missed_cleavages = factor(number_of_missed_cleavages)
      ) %>%
      ggplot(aes(x = number_of_missed_cleavages, y = n)) +
      geom_bar(
        stat = "identity",
        position = "dodge",
        show.legend = FALSE,
        fill = input$plot_color,
        color = "black"
      ) +
      geom_text(aes(label = n), vjust = -0.5, size = 5) +
      labs(x = "Number of Missed Cleavages", y = "Count")
  })

  output$plot12 <- renderPlot({
    data() %>%
      dplyr::count(is_unique) %>%
      as.data.frame() %>%
      dplyr::mutate(
        uniqueness = case_when(
          is_unique == TRUE ~ "Unique",
          TRUE ~ "Shared"
        )
      ) %>%
      ggplot(aes(x = uniqueness, y = n)) +
      geom_bar(
        stat = "identity",
        position = "dodge",
        show.legend = FALSE,
        fill = input$plot_color,
        color = "black"
      ) +
      geom_text(aes(label = n), vjust = -0.5, size = 5) +
      labs(x = "Unique peptides", y = "Count")
  })

  output$plot13 <- renderPlot({
    data() %>%
      as.data.frame() %>%
      ggplot() +
      geom_histogram(
        aes(x = hyperscore),
        fill = input$plot_color,
        color = "black"
      ) +
      labs(
        x = "Hyperscore",
        y = "Count",
        caption = "Similarity score between observed and theoretical spectra.\nHigher values indicate greater similarity"
      )
  })

  output$plot14 <- renderPlot({
    data() %>%
      as.data.frame() %>%
      ggplot() +
      geom_histogram(
        aes(x = nextscore),
        fill = input$plot_color,
        color = "black"
      ) +
      labs(
        x = "Nextscore",
        y = "Count",
        caption = "Second-highest scoring match for the spectrum"
      )
  })

  output$plot15 <- renderPlot({
    data() %>%
      as.data.frame() %>%
      ggplot() +
      geom_histogram(
        aes(x = probability),
        fill = input$plot_color,
        color = "black"
      ) +
      labs(
        x = "PeptideProphet Probability",
        y = "Count",
        caption = "Confidence score determined by PeptideProphet.\nHigher values indicate greater confidence"
      )
  })

  output$plot16 <- renderPlot({
    data() %>%
      as.data.frame() %>%
      ggplot() +
      geom_histogram(
        aes(x = expectation),
        fill = input$plot_color,
        color = "black"
      ) +
      labs(
        x = "Expectation value",
        y = "Count",
        caption = "Expectation value from statistical modeling with PeptideProphet.\nLower values indicate higher likelihood"
      )
  })

  output$plot17 <- renderPlot({
    data() %>%
      as.data.frame() %>%
      tidyr::separate_rows(assigned_modifications, sep = ",") %>%
      dplyr::mutate(
        assigned_modifications = str_remove_all(
          assigned_modifications,
          ".*\\(|\\)"
        ),
        assigned_modifications = ifelse(
          is.na(assigned_modifications),
          "Unassigned modifications",
          assigned_modifications
        )
      ) %>%
      dplyr::count(assigned_modifications) %>%
      ggplot(aes(y = assigned_modifications, x = n)) +
      geom_col(fill = input$plot_color, color = "black") +
      geom_text(aes(label = n), hjust = -0.1, size = 5) +
      scale_x_continuous(expand = expansion(mult = c(0, 0.15))) +
      labs(
        y = "Assigned Modifications",
        x = "Count",
        caption = "Number of modifications assigned to the peptide sequences"
      ) +
      theme(
        axis.text.x = element_text(angle = 90)
      )
  })

  output$plot18 <- renderPlot({
    data() %>%
      as.data.frame() %>%
      dplyr::group_by(entry_name) %>%
      dplyr::summarize(n_psm = n()) %>%
      dplyr::arrange(desc(n_psm)) %>%
      dplyr::mutate(entry_name = factor(entry_name, levels = entry_name)) %>%
      head(20) %>%
      ggplot() +
      geom_bar(
        aes(x = n_psm, y = reorder(entry_name, n_psm)),
        fill = input$plot_color,
        color = "black",
        stat = "identity"
      ) +
      labs(x = "Number of PSMs", y = "Protein") +
      theme(
        axis.text.x = element_text(angle = 90)
      )
  })

  # Download handler for all PSM plots
  output$download_all_plots <- downloadHandler(
    filename = function() {
      paste0("PSM_plots_", Sys.Date(), ".zip")
    },
    content = function(file) {
      temp_dir <- tempdir()
      plots_to_save <- list()

      # Plot 01 - Protease fingerprint
      plots_to_save[["plot01_protease_fingerprint.png"]] <- function() {
        frequency_matrix_of_aa() %>%
          as.data.frame() %>%
          rownames_to_column(var = "residue") %>%
          pivot_longer(
            cols = -residue,
            names_to = "position",
            values_to = "frequency"
          ) %>%
          dplyr::mutate(
            position = factor(
              position,
              c("P4", "P3", "P2", "P1", "P1'", "P2'", "P3'", "P4'")
            ),
            residue = factor(
              residue,
              c(
                "A",
                "C",
                "D",
                "E",
                "F",
                "G",
                "H",
                "I",
                "K",
                "L",
                "M",
                "N",
                "P",
                "Q",
                "R",
                "S",
                "T",
                "V",
                "W",
                "Y"
              )
            )
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
          theme(
            plot.title = element_text(size = 18, hjust = 0.5, face = "bold"),
            axis.text.x = element_text(hjust = 0.5),
            axis.text.y = element_text(hjust = 0.5),
            legend.title = element_text(hjust = 0.5, face = "bold"),
            legend.position = "bottom",
            legend.key.width = unit(2, "cm"),
            legend.key.height = unit(0.25, "cm"),
            legend.title.position = "top"
          )
      }

      # Plot 03 - N-termini SeqLogo
      plots_to_save[["plot03_Nterm_seqlogo.png"]] <- function() {
        data() %>%
          as.data.frame() %>%
          dplyr::select(fingerprint_Nterm) %>%
          na.omit() %>%
          ggseqlogo::ggseqlogo(method = "bits", seq_type = "AA") +
          geom_hline(yintercept = 0, color = "black", linetype = "dashed") +
          geom_vline(xintercept = 4.5, color = "black", linetype = "dashed") +
          scale_x_continuous(
            breaks = c(1, 2, 3, 4, 5, 6, 7, 8),
            labels = c(
              "1" = "P4",
              "2" = "P3",
              "3" = "P2",
              "4" = "P1",
              "5" = "P1'",
              "6" = "P2'",
              "7" = "P3'",
              "8" = "P4'"
            )
          ) +
          labs(
            title = "SeqLogo of the N-termini fingerprint",
            x = "Amino acid position",
            y = "Bits"
          ) +
          theme_bw() +
          theme(
            plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
            text = element_text(size = 15, color = "black"),
            axis.title = element_text(face = "bold"),
            legend.position = "bottom",
            legend.title.position = "top",
            legend.title = element_text(size = 15, hjust = 0.5)
          )
      }

      # Plot 04 - C-termini SeqLogo
      plots_to_save[["plot04_Cterm_seqlogo.png"]] <- function() {
        data() %>%
          as.data.frame() %>%
          dplyr::select(fingerprint_Cterm) %>%
          na.omit() %>%
          ggseqlogo::ggseqlogo(method = "bits", seq_type = "AA") +
          geom_hline(yintercept = 0, color = "black", linetype = "dashed") +
          geom_vline(xintercept = 4.5, color = "black", linetype = "dashed") +
          scale_x_continuous(
            breaks = c(1, 2, 3, 4, 5, 6, 7, 8),
            labels = c(
              "1" = "P4",
              "2" = "P3",
              "3" = "P2",
              "4" = "P1",
              "5" = "P1'",
              "6" = "P2'",
              "7" = "P3'",
              "8" = "P4'"
            )
          ) +
          labs(
            title = "SeqLogo of the C-termini fingerprint",
            x = "Amino acid position",
            y = "Bits"
          ) +
          theme_bw() +
          theme(
            plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
            text = element_text(size = 15, color = "black"),
            axis.title = element_text(face = "bold"),
            legend.position = "bottom",
            legend.title.position = "top",
            legend.title = element_text(size = 15, hjust = 0.5)
          )
      }

      # Plot 05 - m/z over retention time
      plots_to_save[["plot05_mz_retention.png"]] <- function() {
        data() %>%
          as.data.frame() %>%
          ggplot(aes(x = retention / 60, y = observed_m_z)) +
          ggpointdensity::geom_pointdensity(size = 0.25) +
          viridis::scale_color_viridis(option = "plasma") +
          labs(
            x = "Retention time (min)",
            y = "Scan range (m/z)",
            color = "Number of Neighborhoods"
          ) +
          theme(
            legend.position = "bottom",
            legend.key.width = unit(2, "cm"),
            legend.key.height = unit(0.25, "cm")
          )
      }

      # Plot 06 - Mass error
      plots_to_save[["plot06_mass_error.png"]] <- function() {
        data() %>%
          as.data.frame() %>%
          dplyr::filter(abs(delta_mass_ppm) < 100) %>%
          ggplot(aes(x = retention / 60, y = delta_mass_ppm)) +
          geom_point(alpha = 0.1, color = "black", size = 1) +
          geom_hline(
            yintercept = c(10, 0, -10),
            color = "red",
            linetype = "dashed",
            linewidth = 0.2
          ) +
          labs(
            title = "Mass error in ppm",
            x = "Retention time (min)",
            y = "Mass error (ppm)",
            caption = "ppm error is calculated as:\n∆m/z over theoretical m/z * 1e6"
          )
      }

      # Plot 07 - Peptide length
      plots_to_save[["plot07_peptide_length.png"]] <- function() {
        data() %>%
          as.data.frame() %>%
          ggplot() +
          geom_density(aes(x = peptide_length), fill = input$plot_color) +
          labs(x = "Peptide Length", y = "Frequency (%)") +
          theme(
            text = element_text(size = 15, color = "black"),
            axis.text.x = element_text(hjust = 0.5),
            plot.title = element_text(hjust = 0.5)
          )
      }

      # Plot 08 - GRAVY
      plots_to_save[["plot08_gravy.png"]] <- function() {
        data() %>%
          as.data.frame() %>%
          ggplot(aes(x = gravy, fill = stat(x))) +
          geom_histogram(color = "black") +
          labs(
            x = NULL,
            y = "Count",
            caption = "GRAVY is a measure of the hydropathic character of a sequence"
          ) +
          scale_fill_viridis_c(name = "GRAVY index", option = "C") +
          theme(
            text = element_text(size = 15, color = "black"),
            axis.text.x = element_text(hjust = 0.5),
            legend.position = "bottom",
            legend.key.width = unit(2.5, "cm"),
            legend.key.height = unit(0.25, "cm")
          )
      }

      # Plot 09 - Isoelectric Point
      plots_to_save[["plot09_isoelectric_point.png"]] <- function() {
        data() %>%
          as.data.frame() %>%
          ggplot(aes(x = isoelectric_point, fill = stat(x))) +
          geom_histogram(color = "black") +
          labs(x = NULL, y = "Count") +
          scale_fill_viridis_c(name = "Isoelectric Point (pI)", option = "C") +
          theme(
            text = element_text(size = 15, color = "black"),
            axis.text.x = element_text(hjust = 0.5),
            legend.position = "bottom",
            legend.key.width = unit(2.5, "cm"),
            legend.key.height = unit(0.25, "cm")
          )
      }

      # Plot 10 - Charge state
      plots_to_save[["plot10_charge_state.png"]] <- function() {
        data() %>%
          as.data.frame() %>%
          ggplot() +
          geom_bar(aes(x = charge), fill = input$plot_color, color = "black") +
          labs(x = "Charge state", y = "Count")
      }

      # Plot 11 - Missed cleavages
      plots_to_save[["plot11_missed_cleavages.png"]] <- function() {
        data() %>%
          as.data.frame() %>%
          dplyr::count(number_of_missed_cleavages) %>%
          dplyr::mutate(
            number_of_missed_cleavages = factor(number_of_missed_cleavages)
          ) %>%
          ggplot(aes(x = number_of_missed_cleavages, y = n)) +
          geom_bar(
            stat = "identity",
            position = "dodge",
            show.legend = FALSE,
            fill = input$plot_color,
            color = "black"
          ) +
          geom_text(aes(label = n), vjust = -0.5, size = 5) +
          labs(x = "Number of Missed Cleavages", y = "Count")
      }

      # Plot 12 - Uniqueness
      plots_to_save[["plot12_uniqueness.png"]] <- function() {
        data() %>%
          dplyr::count(is_unique) %>%
          as.data.frame() %>%
          dplyr::mutate(
            uniqueness = case_when(
              is_unique == TRUE ~ "Unique",
              TRUE ~ "Shared"
            )
          ) %>%
          ggplot(aes(x = uniqueness, y = n)) +
          geom_bar(
            stat = "identity",
            position = "dodge",
            show.legend = FALSE,
            fill = input$plot_color,
            color = "black"
          ) +
          geom_text(aes(label = n), vjust = -0.5, size = 5) +
          labs(x = "Unique peptides", y = "Count")
      }

      # Plot 13 - Hyperscore
      plots_to_save[["plot13_hyperscore.png"]] <- function() {
        data() %>%
          as.data.frame() %>%
          ggplot() +
          geom_histogram(
            aes(x = hyperscore),
            fill = input$plot_color,
            color = "black"
          ) +
          labs(
            x = "Hyperscore",
            y = "Count",
            caption = "Similarity score between observed and theoretical spectra.\nHigher values indicate greater similarity"
          )
      }

      # Plot 14 - Nextscore
      plots_to_save[["plot14_nextscore.png"]] <- function() {
        data() %>%
          as.data.frame() %>%
          ggplot() +
          geom_histogram(
            aes(x = nextscore),
            fill = input$plot_color,
            color = "black"
          ) +
          labs(
            x = "Nextscore",
            y = "Count",
            caption = "Second-highest scoring match for the spectrum"
          )
      }

      # Plot 15 - PeptideProphet probability
      plots_to_save[["plot15_peptideprophet_probability.png"]] <- function() {
        data() %>%
          as.data.frame() %>%
          ggplot() +
          geom_histogram(
            aes(x = probability),
            fill = input$plot_color,
            color = "black"
          ) +
          labs(
            x = "PeptideProphet Probability",
            y = "Count",
            caption = "Confidence score determined by PeptideProphet.\nHigher values indicate greater confidence"
          )
      }

      # Plot 16 - Expectation
      plots_to_save[["plot16_expectation.png"]] <- function() {
        data() %>%
          as.data.frame() %>%
          ggplot() +
          geom_histogram(
            aes(x = expectation),
            fill = input$plot_color,
            color = "black"
          ) +
          labs(
            x = "Expectation value",
            y = "Count",
            caption = "Expectation value from statistical modeling with PeptideProphet.\nLower values indicate higher likelihood"
          )
      }

      # Plot 17 - Assigned modifications
      plots_to_save[["plot17_assigned_modifications.png"]] <- function() {
        data() %>%
          as.data.frame() %>%
          tidyr::separate_rows(assigned_modifications, sep = ",") %>%
          dplyr::mutate(
            assigned_modifications = str_remove_all(
              assigned_modifications,
              ".*\\(|\\)"
            ),
            assigned_modifications = ifelse(
              is.na(assigned_modifications),
              "Unassigned modifications",
              assigned_modifications
            )
          ) %>%
          dplyr::count(assigned_modifications) %>%
          ggplot(aes(y = assigned_modifications, x = n)) +
          geom_col(fill = input$plot_color, color = "black") +
          geom_text(aes(label = n), hjust = -0.1, size = 5) +
          scale_x_continuous(expand = expansion(mult = c(0, 0.15))) +
          labs(
            y = "Assigned Modifications",
            x = "Count",
            caption = "Number of modifications assigned to the peptide sequences"
          ) +
          theme(axis.text.x = element_text(angle = 90))
      }

      # Plot 18 - Top 20 proteins
      plots_to_save[["plot18_top20_proteins.png"]] <- function() {
        data() %>%
          as.data.frame() %>%
          dplyr::group_by(entry_name) %>%
          dplyr::summarize(n_psm = n()) %>%
          dplyr::arrange(desc(n_psm)) %>%
          dplyr::mutate(
            entry_name = factor(entry_name, levels = entry_name)
          ) %>%
          head(20) %>%
          ggplot() +
          geom_bar(
            aes(x = n_psm, y = reorder(entry_name, n_psm)),
            fill = input$plot_color,
            color = "black",
            stat = "identity"
          ) +
          labs(x = "Number of PSMs", y = "Protein") +
          theme(axis.text.x = element_text(angle = 90))
      }

      # Save all plots as PNG files with high resolution (300 DPI)
      file_paths <- c()
      for (plot_name in names(plots_to_save)) {
        file_path <- file.path(temp_dir, plot_name)
        file_paths <- c(file_paths, file_path)

        ggsave(
          filename = file_path,
          plot = plots_to_save[[plot_name]](),
          width = 12,
          height = 8,
          dpi = 300,
          units = "in",
          bg = "white"
        )
      }

      utils::zip(file, files = file_paths, flags = "-j")
    },
    contentType = "application/zip"
  )

  # Import and pre-process the uploaded protein.tsv file
  protein_data <- reactive({
    req(input$protein)
    protein_file <- readr::read_tsv(input$protein$datapath) %>%
      janitor::clean_names()
  })

  output$info_box2 <- renderInfoBox({
    infoBox(
      "protein.tsv files contain FDR-filtered protein results, where each row is an identified protein group",
      icon = icon("info"),
      color = "black"
    )
  })

  # Render plots for the protein viewer
  output$plot19 <- renderPlot({
    protein_data() %>%
      as.data.frame() %>%
      ggplot() +
      geom_histogram(
        aes(x = coverage),
        fill = input$plot_color,
        color = "black"
      ) +
      labs(x = "Protein coverage (%)", y = "Count")
  })

  output$plot20 <- renderPlot({
    protein_data() %>%
      as.data.frame() %>%
      dplyr::count(organism) %>%
      ggplot() +
      geom_bar(
        aes(x = n, y = reorder(organism, n)),
        fill = input$plot_color,
        color = "black",
        stat = "identity"
      ) +
      geom_text(
        aes(x = n, y = reorder(organism, n), label = n),
        hjust = -0.1,
        size = 5
      ) +
      labs(x = "Number of proteins", y = NULL) +
      theme(axis.text.y = element_text(face = "italic"))
  })

  output$plot21 <- renderPlot({
    protein_data() %>%
      as.data.frame() %>%
      dplyr::mutate(
        protein_existence = str_remove(protein_existence, ".*\\:"),
        protein_existence = factor(
          protein_existence,
          levels = c(
            "Experimental evidence at protein level",
            "Experimental evidence at transcript level",
            "Protein inferred from homology",
            "Protein predicted"
          )
        )
      ) %>%
      ggplot() +
      geom_bar(
        aes(y = protein_existence, fill = protein_existence),
        color = "black",
        show.legend = FALSE
      ) +
      scale_fill_manual(
        values = c("#5499c7", "#7fb3d5", "#a9cce3", "#d4e6f1")
      ) +
      geom_text(
        aes(y = protein_existence, label = ..count..),
        show.legend = FALSE,
        stat = "count",
        vjust = -0.5,
        size = 7,
        fontface = "bold"
      ) +
      labs(y = NULL, x = "Count", fill = NULL)
  })

  output$plot22 <- renderPlot({
    protein_data() %>%
      as.data.frame() %>%
      ggplot() +
      geom_histogram(
        aes(x = protein_probability),
        fill = input$plot_color,
        color = "black"
      ) +
      labs(x = "Protein Probability", y = "Count")
  })

  output$plot23 <- renderPlot({
    protein_data() %>%
      as.data.frame() %>%
      ggplot() +
      geom_histogram(
        aes(x = top_peptide_probability),
        fill = input$plot_color,
        color = "black"
      ) +
      labs(
        x = "Peptide Probability",
        y = "Count",
        caption = "Best peptide probability of supporting peptides"
      )
  })

  output$plot24 <- renderPlot({
    protein_data() %>%
      as.data.frame() %>%
      ggplot() +
      geom_histogram(
        aes(x = total_peptides),
        fill = input$plot_color,
        color = "black"
      ) +
      labs(x = "Total peptides mapped to proteins", y = "Count")
  })

  output$plot25 <- renderPlot({
    protein_data() %>%
      as.data.frame() %>%
      ggplot() +
      geom_histogram(
        aes(x = razor_spectral_count),
        fill = input$plot_color,
        color = "black"
      ) +
      labs(
        x = "Razor Spectral Count",
        y = "Count",
        caption = "Number of PSMs corresponding to the razor peptides"
      )
  })

  output$plot26 <- renderPlot({
    protein_data() %>%
      as.data.frame() %>%
      ggplot() +
      geom_histogram(
        aes(x = razor_intensity),
        fill = input$plot_color,
        color = "black"
      ) +
      labs(
        x = "Razor Intensity",
        y = "Count",
        caption = "Protein intensity calculated using the unique peptides (from the top-N algorithm)"
      )
  })

  output$plot27 <- renderPlot({
    protein_data() %>%
      as.data.frame() %>%
      dplyr::arrange(desc(razor_intensity)) %>%
      head(20) %>%
      ggplot() +
      geom_bar(
        aes(
          x = log2(razor_intensity),
          y = reorder(entry_name, razor_intensity)
        ),
        fill = input$plot_color,
        color = "black",
        stat = "identity"
      ) +
      labs(x = "log2 of Razor Intensity", y = "Protein")
  })

  # Import and pre-process the uploaded combined_protein.tsv file
  combined_protein_data <- reactive({
    req(input$combined_protein)
    combined_protein_file <- readr::read_tsv(
      input$combined_protein$datapath
    ) %>%
      janitor::clean_names() %>%
      dplyr::select(protein_id, ends_with("max_lfq_intensity")) %>%
      column_to_rownames("protein_id") %>%
      dplyr::rename_all(~ str_remove(., "_max_lfq_intensity")) %>%
      log2()
  })

  # Observe the uploaded file and update selectInput choices
  observe({
    req(combined_protein_data())
    colnames <- colnames(combined_protein_data())
    updateSelectInput(session, "xcol", choices = colnames)
    updateSelectInput(session, "ycol", choices = colnames)
  })

  output$plot28 <- renderPlot({
    combined_protein_data() %>%
      as.data.frame() %>%
      rownames_to_column(var = "protein_id") %>%
      tidyr::pivot_longer(
        cols = -protein_id,
        names_to = "sample",
        values_to = "maxlfq_intensity"
      ) %>%
      ggplot() +
      geom_violin(
        aes(x = sample, y = maxlfq_intensity),
        fill = input$plot_color,
        alpha = 0.7,
        color = "black"
      ) +
      geom_boxplot(
        aes(x = sample, y = maxlfq_intensity),
        fill = "white",
        outliers = FALSE,
        color = "black",
        width = 0.1,
        show.legend = FALSE
      ) +
      labs(x = NULL, y = "log2(MaxLFQ intensity)") +
      theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
  })

  output$plot29 <- renderPlotly({
    combined_protein_data() %>%
      as.data.frame() %>%
      ggplot(aes(x = !!sym(input$xcol), y = !!sym(input$ycol))) +
      geom_point(alpha = 0.7, show.legend = FALSE) +
      geom_smooth(method = "lm", se = FALSE, color = input$plot_color) +
      labs(
        x = paste0("log2(", input$xcol, ")"),
        y = paste0("log2(", input$ycol, ")")
      )
  })

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
      theme(
        text = element_text(size = 15),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        axis.text.y = element_text(angle = 0, hjust = 1, vjust = 0.5),
        legend.position = "bottom",
        legend.key.width = unit(2.5, "cm"),
        legend.key.height = unit(0.25, "cm")
      ) +
      labs(x = NULL, y = NULL, fill = "Cosine similarity")
  })

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
      theme(
        text = element_text(size = 15),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        axis.text.y = element_text(angle = 0, hjust = 1, vjust = 0.5),
        legend.position = "bottom",
        legend.key.width = unit(2.5, "cm"),
        legend.key.height = unit(0.25, "cm")
      ) +
      labs(x = NULL, y = NULL, fill = "Euclidean distance")
  })

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
      theme(
        text = element_text(size = 15),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        axis.text.y = element_text(angle = 0, hjust = 1, vjust = 0.5),
        legend.position = "bottom",
        legend.key.width = unit(2.5, "cm"),
        legend.key.height = unit(0.25, "cm")
      ) +
      labs(x = NULL, y = NULL, fill = "Jaccard similarity")
  })
}

shinyApp(ui, server)
