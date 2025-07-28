## PSM analyst dashboard for FragPipe search results
## The files psm.tsv, protein.tsv and combined_protein.tsv are the inputs for the PSManalyst dashboard
## It is possible to filter the PSMs by the hyperscore and PeptideProphet probability
## You can customize the color of the plots using the Select plot color option
## you can remove one organism from the analysis if you think it is a contaminant

# Check if the required R libraries are installed and install them if necessary.
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
  "colourpicker",
  "ggpointdensity"
)
not_installed_CRAN <- CRAN_packages[
  !(CRAN_packages %in% installed.packages()[, "Package"])
]
if (length(not_installed_CRAN)) {
  install.packages(not_installed_CRAN)
}

GitHub_packages <- c("wordcloud2")
not_installed_GitHub <- GitHub_packages[
  !(GitHub_packages %in% installed.packages()[, "Package"])
]
if (length(not_installed_GitHub)) {
  install.packages(not_installed_GitHub)
}

# Load required libraries
library(shiny) # from CRAN
library(shinydashboard) # from CRAN
library(tidyverse) # from CRAN
library(janitor) # from CRAN
library(ggseqlogo) # from CRAN
library(ggpointdensity) # from CRAN
library(wordcloud2) # from GitHub
library(ggtext) # from CRAN
library(lsa) # from CRAN
library(plotly) # from CRAN
library(viridis) # from CRAN
library(ggfortify) # from CRAN
library(vegan) # from CRAN
library(ggwordcloud) # from CRAN
library(colourpicker) # from CRAN

# Increase the maximum file size to 1000 MB
options(shiny.maxRequestSize = 1000 * 1024^2)

theme_set(theme_bw())
theme_update(
  text = element_text(color = "black", size = 12),
  axis.text = element_text(color = "black"),
  axis.title = element_text(color = "black", face = "bold"),
  strip.background = element_blank(),
  strip.text = element_text(face = "bold"),
  legend.title = element_text(face = "bold", hjust = 0.5),
  legend.title.position = "top",
  panel.grid = element_blank()
)

color_blue_seq <- c(
  "#d4e6f1",
  "#a9cce3",
  "#7fb3d5",
  "#5499c7",
  "#2980b9",
  "#1f618d",
  "#154360"
)

# GRAVY (Grand Average of Hydropathy)
# Kyte J, Doolittle RF. A simple method for displaying the hydropathic character of a protein. J Mol Biol. 1982 May 5;157(1):105-32. doi: 10.1016/0022-2836(82)90515-0.
GRAVY <- function(sequence) {
  # Define the hydropathy index for each amino acid
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
  # Calculate the GRAVY index for the peptide sequences
  scores <- sapply(strsplit(sequence, NULL)[[1]], function(aa) {
    hydropathy_index[aa]
  })
  return(mean(scores, na.rm = TRUE))
}

# Calculate the isoelectric point (pI) of peptide sequences
calculate_pI <- function(sequence) {
  # Define the pKa values for the amino acids
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

  # Calculate the pI based on the sequence
  pI <- mean(
    sapply(strsplit(sequence, NULL)[[1]], function(aa) pKa_values[aa]),
    na.rm = TRUE
  )

  return(pI)
}

# Calculate the Margalef’s index
# Margalef R. Information theory in ecology. General Systems 3. 1958;36-71.
calculate_margalef <- function(data) {
  S <- length(unique(data$peptide))
  N <- nrow(data)
  # Calculate Margalef's index: (S - 1) / ln(N)
  if (N > 1) {
    margalef_index <- (S - 1) / log(N)
  } else {
    margalef_index <- 0
  }
  
  return(margalef_index)
}

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
      textInput(
        "data_directory",
        label = "Path to search for PSM files",
        value = getwd(),
        placeholder = "/path/to/your/data"
      ),
      div(style = "text-align: center; margin: 10px 0; display: flex; justify-content: center;",
          actionButton("load_psm_files", "Load PSM Files", class = "btn-primary")
      ),
      br(),
      verbatimTextOutput("psm_files_status"),
      br(),
      sliderInput(
        "hyperscore",
        label = "PSM hyperscore filter",
        min = 0,
        max = 1000,
        value = 10,
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
      textInput(
        "protein_data_directory",
        label = "Path to search for protein.tsv",
        value = getwd(),
        placeholder = "/path/to/your/data"
      ),
      div(style = "text-align: center; margin: 10px 0; display: flex; justify-content: center;",
        actionButton(
        "load_protein_files",
        "Load protein files",
        class = "btn-primary"
      )
    ),
      br(),
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
      )
    )
  ),

  dashboardBody(
    tabItems(
      tabItem(
        tabName = "psm",
        fluidRow(
          infoBoxOutput("info_box1", width = 12),
          box(
            title = "m/z over retention time",
            status = "primary",
            solidHeader = TRUE,
            width = 12,
            plotOutput("plot1"),
            collapsible = TRUE
          ),
          box(
            title = "Mass error (ppm)",
            status = "primary",
            solidHeader = TRUE,
            width = 12,
            plotOutput("plot2"),
            collapsible = TRUE
          ),
          box(
            title = "Peptide length",
            status = "primary",
            solidHeader = TRUE,
            width = 12,
            plotOutput("plot3"),
            collapsible = TRUE
          ),
          box(
            title = "Charge state distribution",
            status = "primary",
            solidHeader = TRUE,
            width = 12,
            plotOutput("plot4"),
            collapsible = TRUE
          ),
          box(
            title = "Number of missed cleavages",
            status = "primary",
            solidHeader = TRUE,
            width = 12,
            plotOutput("plot5"),
            collapsible = TRUE
          ),
          box(
            title = "Uniqueness",
            status = "primary",
            solidHeader = TRUE,
            width = 12,
            plotOutput("plot6"),
            collapsible = TRUE
          ),
          box(
            title = "Hyperscore distribution",
            status = "primary",
            solidHeader = TRUE,
            width = 12,
            plotOutput("plot7"),
            collapsible = TRUE
          ),
          box(
            title = "Next Score distribution",
            status = "primary",
            solidHeader = TRUE,
            width = 12,
            plotOutput("plot8"),
            collapsible = TRUE
          ),
          box(
            title = "PeptideProphet probability",
            status = "primary",
            solidHeader = TRUE,
            width = 12,
            plotOutput("plot9"),
            collapsible = TRUE
          ),
          box(
            title = "Expectation (PeptideProphet)",
            status = "primary",
            solidHeader = TRUE,
            width = 12,
            plotOutput("plot10"),
            collapsible = TRUE
          ),
          box(
            title = "Assigned modifications",
            status = "primary",
            solidHeader = TRUE,
            width = 12,
            plotOutput("plot11"),
            collapsible = TRUE
          ),
          box(
            title = "Top 20 proteins with more PSMs",
            status = "primary",
            solidHeader = TRUE,
            width = 12,
            plotOutput("plot12"),
            collapsible = TRUE
          ),
          box(
            title = "PSMs by Source Folder",
            status = "primary",
            solidHeader = TRUE,
            width = 12,
            plotOutput("plot13"),
            collapsible = TRUE
          ),
          box(
            title = "GRAVY (Grand Average of Hydropathy)",
            status = "primary",
            solidHeader = TRUE,
            width = 12,
            plotOutput("plot_gravy"),
            collapsible = TRUE
          ),
          box(
            title = "Isoelectric point (pI)",
            status = "primary",
            solidHeader = TRUE,
            width = 12,
            plotOutput("plot_pI"),
            collapsible = TRUE
          ),
          box(
            title = "Margalef's Diversity Index",
            status = "primary",
            solidHeader = TRUE,
            width = 12,
            plotOutput("plot_margalef"),
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
            width = 12,
            plotOutput("plot14"),
            collapsible = TRUE
          ),
          box(
            title = "Number of proteins by organism",
            status = "primary",
            solidHeader = TRUE,
            width = 12,
            height = 700,
            plotOutput("plot15", height = "600px"),
            collapsible = TRUE
          ),
          box(
            title = "Protein existence evidence",
            status = "primary",
            solidHeader = TRUE,
            width = 12,
            plotOutput("plot16"),
            collapsible = TRUE
          ),
          box(
            title = "Protein probability (ProteinProphet)",
            status = "primary",
            solidHeader = TRUE,
            width = 12,
            plotOutput("plot17"),
            collapsible = TRUE
          ),
          box(
            title = "Top Peptide Probability",
            status = "primary",
            solidHeader = TRUE,
            width = 12,
            plotOutput("plot18"),
            collapsible = TRUE
          ),
          box(
            title = "Total peptides mapped to the proteins",
            status = "primary",
            solidHeader = TRUE,
            width = 12,
            plotOutput("plot19"),
            collapsible = TRUE
          ),
          box(
            title = "Razor spectral count",
            status = "primary",
            solidHeader = TRUE,
            width = 12,
            plotOutput("plot20"),
            collapsible = TRUE
          ),
          box(
            title = "Razor intensity",
            status = "primary",
            solidHeader = TRUE,
            width = 12,
            plotOutput("plot21"),
            collapsible = TRUE
          ),
          box(
            title = "Top 20 proteins with higher razor intensity",
            status = "primary",
            solidHeader = TRUE,
            width = 12,
            plotOutput("plot22"),
            collapsible = TRUE
          ),
          box(
            title = "MaxLFQ intensity distribution",
            status = "primary",
            solidHeader = TRUE,
            width = 12,
            plotOutput("plot23"),
            collapsible = TRUE
          ),
          box(
            title = "Sample correlation - Non-normalized log2(Intensity)",
            status = "primary",
            height = 600,
            solidHeader = TRUE,
            width = 12,
            plotlyOutput("plot24"),
            collapsible = FALSE
          ),
          tabBox(
            title = "Similarity metrics",
            side = "right",
            height = 750,
            width = 12,
            tabPanel("Cosine similarity", plotOutput("cosine_similarity", height = "700px")),
            tabPanel("Euclidean distance", plotOutput("euclidean_distance", height = "700px")),
            tabPanel("Jaccard similarity", plotOutput("jaccard_similarity", height = "700px"))
          )
        )
      )
    )
  )
)

server <- function(input, output, session) {
  multiple_psm_data <- eventReactive(input$load_psm_files, {
    req(input$data_directory)
    if (!dir.exists(input$data_directory)) {
      showNotification("Directory does not exist!", type = "error")
      return(NULL)
    }

    psm_files <- list.files(
      path = input$data_directory,
      pattern = "^psm\\.tsv$",
      recursive = TRUE,
      full.names = TRUE
    )
    if (length(psm_files) == 0) {
      showNotification("No psm.tsv files found!", type = "warning")
      return(NULL)
    }
    showNotification(
      paste("Found", length(psm_files), "PSM files"),
      type = "message"
    )
    all_data <- NULL
    for (file_path in psm_files) {
      tryCatch(
        {
          subfolder_name <- basename(dirname(file_path))
          psm_data <- readr::read_tsv(file_path, show_col_types = FALSE)
          psm_data <- janitor::clean_names(psm_data)
          psm_data$source_folder <- subfolder_name
          if (is.null(all_data)) {
            all_data <- psm_data
          } else {
            all_data <- dplyr::bind_rows(all_data, psm_data)
          }
        },
        error = function(e) {
          showNotification(
            paste("Error reading", basename(file_path), ":", e$message),
            type = "error"
          )
        }
      )
    }

    return(all_data)
  })
  
  output$psm_files_status <- renderText({
    if (input$load_psm_files == 0) {
      return("Click 'Load PSM Files' to search for files")
    }

    multi_data <- multiple_psm_data()
    if (!is.null(multi_data)) {
      n_folders <- length(unique(multi_data$source_folder))
      folders <- paste(unique(multi_data$source_folder), collapse = ", ")
      paste(
        "Loaded",
        nrow(multi_data),
        "PSMs from",
        n_folders,
        "folders:\n",
        folders
      )
    } else {
      "No PSM files loaded"
    }
  })

  data <- reactive({
    psm_file <- NULL
    if (input$load_psm_files > 0) {
      multi_data <- multiple_psm_data()
      if (!is.null(multi_data)) {
        psm_file <- multi_data
      }
    }
    if (is.null(psm_file) && !is.null(input$psm)) {
      psm_file <- readr::read_tsv(
        input$psm$datapath,
        show_col_types = FALSE
      ) %>%
        janitor::clean_names() %>%
        mutate(source_folder = "single_upload")
    }
    if (is.null(psm_file)) {
      return(NULL)
    }
    
    psm_file <- psm_file %>%
      dplyr::filter(
        hyperscore >= input$hyperscore & probability >= input$probability
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
        isoelectric_point = sapply(peptide, calculate_pI)
      ) %>%
      dplyr::group_by(source_folder) %>%
      dplyr::mutate(
        margalef_index = calculate_margalef(cur_data())
  ) %>%
      dplyr::ungroup() %>%
      dplyr::relocate(source_folder, .after = extended_peptide)

    return(psm_file)
  })

  output$info_box1 <- renderInfoBox({
    filter_text <- paste(
      "Showing PSMs with Hyperscore ≥",
      input$hyperscore,
      "and PeptideProphet probability ≥",
      input$probability
    )

    if (!is.null(data()) && "source_folder" %in% colnames(data())) {
      n_folders <- length(unique(data()$source_folder))
      if (n_folders > 1) {
        filter_text <- paste(filter_text, "from", n_folders, "folders")
      }
    }

    if (!is.null(input$protein_pattern) && input$protein_pattern != "") {
      pattern_text <- if (input$case_sensitive) {
        paste(
          "and organism entry name matching:",
          input$protein_pattern,
          "(case sensitive)"
        )
      } else {
        paste(
          "and organism entry name matching:",
          input$protein_pattern,
          "(case insensitive)"
        )
      }
      filter_text <- paste(filter_text, pattern_text)
    }

    infoBox(
      "Filter settings",
      filter_text,
      icon = icon("info"),
      color = "black"
    )
  })

  output$plot1 <- renderPlot({
    data() %>%
      as.data.frame() %>%
      ggplot(aes(x = retention / 60, y = observed_m_z)) +
      stat_density_2d(
        geom = "raster",
        aes(
          fill = after_stat(density)
        ),
        contour = FALSE
      ) +
      scale_fill_viridis_c() +
      labs(
        x = "Retention time (min)",
        y = "Scan range (m/z)",
        fill = "Density"
      ) +
      facet_wrap(~source_folder) +
      theme(
        strip.text = element_text(size = 12),
        legend.position = "bottom",
        legend.key.width = unit(2, "cm"),
        legend.key.height = unit(0.25, "cm")
      )
  })

  output$plot2 <- renderPlot({
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
        caption = "ppm error is calculated as ∆m/z over theoretical m/z * 1e6"
      ) +
      facet_wrap(~source_folder) +
      theme(
        strip.text = element_text(size = 12),
        plot.title = element_text(hjust = 0.5, face = "bold")
      )
  })

  output$plot3 <- renderPlot({
    data() %>%
      as.data.frame() %>%
      ggplot() +
      geom_density(aes(x = peptide_length), fill = input$plot_color) +
      labs(x = "Peptide Length", y = "Frequency (%)") +
      facet_wrap(~source_folder) +
      theme(
        text = element_text(size = 15, color = "black"),
        axis.text.x = element_text(hjust = 0.5),
        strip.text = element_text(size = 12),
        plot.title = element_text(hjust = 0.5),
      )
  })

  output$plot4 <- renderPlot({
    data() %>%
      as.data.frame() %>%
      ggplot() +
      geom_bar(aes(x = charge), fill = input$plot_color, color = "black") +
      labs(x = "Charge state", y = "Count") +
      facet_wrap(~source_folder) +
      theme(
        text = element_text(size = 15, color = "black"),
        axis.text.x = element_text(hjust = 0.5),
        strip.text = element_text(size = 12),
        plot.title = element_text(hjust = 0.5)
      )
  })

  output$plot5 <- renderPlot({
    data() %>%
      as.data.frame() %>%
      dplyr::group_by(source_folder, number_of_missed_cleavages) %>%
      dplyr::summarise(n = n(), .groups = "drop") %>%
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
      geom_text(aes(label = n), vjust = -0.5, size = 4) +
      labs(x = "Number of Missed Cleavages", y = "Count") +
      facet_wrap(~source_folder, scales = "free_y") +
      theme(
        text = element_text(size = 15, color = "black"),
        axis.text.x = element_text(hjust = 0.5),
        strip.text = element_text(size = 12),
        plot.title = element_text(hjust = 0.5)
      )
  })

  output$plot6 <- renderPlot({
    data() %>%
      dplyr::group_by(source_folder, is_unique) %>%
      dplyr::summarise(n = n(), .groups = "drop") %>%
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
      geom_text(aes(label = n), vjust = -0.5, size = 4) +
      labs(x = "Unique peptides", y = "Count") +
      facet_wrap(~source_folder) +
      theme(
        text = element_text(size = 15, color = "black"),
        axis.text.x = element_text(hjust = 0.5),
        strip.text = element_text(size = 12),
        plot.title = element_text(hjust = 0.5)
      )
  })

  output$plot7 <- renderPlot({
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
        caption = "Similarity score between observed and theoretical spectra, higher values indicate greater similarity"
      ) +
      facet_wrap(~source_folder) +
      theme(
        text = element_text(size = 15, color = "black"),
        axis.text.x = element_text(hjust = 0.5),
        strip.text = element_text(size = 12),
        plot.title = element_text(hjust = 0.5)
      )
  })

  output$plot8 <- renderPlot({
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
        caption = "Similarity score (hyperscore) of the second-highest scoring match for the spectrum"
      ) +
      facet_wrap(~source_folder) +
      theme(
        text = element_text(size = 15, color = "black"),
        axis.text.x = element_text(hjust = 0.5),
        strip.text = element_text(size = 12),
        plot.title = element_text(hjust = 0.5)
      )
  })

  output$plot9 <- renderPlot({
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
        caption = "Confidence score determined by PeptideProphet, higher values indicate greater confidence"
      ) +
      facet_wrap(~source_folder) +
      theme(
        text = element_text(size = 15, color = "black"),
        axis.text.x = element_text(hjust = 1, vjust = 0.5, angle = 90),
        strip.text = element_text(size = 12),
        plot.title = element_text(hjust = 0.5)
      )
  })

  output$plot10 <- renderPlot({
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
        caption = "Expectation value from statistical modeling with PeptideProphet, lower values indicate higher likelihood"
      ) +
      facet_wrap(~source_folder) +
      theme(
        text = element_text(size = 15, color = "black"),
        axis.text.x = element_text(hjust = 0.5),
        strip.text = element_text(size = 12),
        plot.title = element_text(hjust = 0.5)
      )
  })

  output$plot11 <- renderPlot({
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
          "Unmodified",
          assigned_modifications
        )
      ) %>%
      dplyr::group_by(source_folder, assigned_modifications) %>%
      dplyr::summarize(n = n(), .groups = "drop") %>%
      ggplot(aes(y = assigned_modifications, x = n)) +
      geom_col(fill = input$plot_color, color = "black") +
      geom_text(aes(label = n), hjust = -0.1, size = 4) +
      labs(
        y = "Assigned Modifications",
        x = "Count",
        caption = "Number of modifications assigned to the peptide sequences"
      ) +
      facet_wrap(~source_folder) +
      theme(
        text = element_text(size = 15, color = "black"),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        strip.text = element_text(size = 12),
        plot.title = element_text(hjust = 0.5)
      )
  })

  output$plot12 <- renderPlot({
    data() %>%
      as.data.frame() %>%
      dplyr::group_by(source_folder, entry_name) %>%
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
      facet_wrap(~source_folder, scales = "free_y") +
      theme(
        text = element_text(size = 15, color = "black"),
        axis.text.x = element_text(hjust = 1, angle = 90, vjust = 0.5),
        strip.text = element_text(size = 12),
        plot.title = element_text(hjust = 0.5)
      )
  })

  output$plot13 <- renderPlot({
    req(data())
    if ("source_folder" %in% colnames(data())) {
      data() %>%
        count(source_folder, sort = TRUE) %>%
        ggplot(aes(y = reorder(source_folder, n), x = n)) +
        geom_col(fill = input$plot_color, color = "black") +
        geom_text(aes(label = scales::comma(n)), hjust = -0.1, size = 4) +
        labs(
          y = "Source Folder",
          x = "Number of PSMs",
          title = "PSM Distribution by File"
        ) +
        theme(plot.title = element_text(hjust = 0.5, face = "bold"))
    } else {
      ggplot() +
        geom_text(
          aes(x = 1, y = 1, label = "No source folder data available"),
          size = 6
        ) +
        theme_void()
    }
  })

  output$plot_gravy <- renderPlot({
    data() %>%
      as.data.frame() %>%
      ggplot(aes(x = gravy, fill = stat(x))) +
      geom_histogram(olor = "black") +
      scale_fill_gradient2(
        name = "GRAVY Index",
        low = "#5499c7",
        mid = "grey90",
        high = "firebrick",
        midpoint = 0
      ) +
      labs(
        x = NULL,
        y = "Count",
        caption = "GRAVY is a measure of the hydropathic character of a sequence"
      ) +
      facet_wrap(~source_folder, scales = "free_y") +
      theme(
        text = element_text(size = 15, color = "black"),
        axis.text.x = element_text(hjust = 0.5),
        strip.text = element_text(size = 12),
        legend.position = "bottom",
        legend.key.width = unit(2.5, "cm"),
        legend.key.height = unit(0.25, "cm")
      )
  })

  output$plot_pI <- renderPlot({
    data() %>%
      as.data.frame() %>%
      ggplot(aes(x = isoelectric_point, fill = stat(x))) +
      geom_histogram(color = "black") +
      labs(x = NULL, y = "Count") +
      scale_fill_viridis_c(name = "Isoelectric Point (pI)", option = "C") +
      facet_wrap(~source_folder, scales = "free_y") +
      theme(
        text = element_text(size = 15, color = "black"),
        axis.text.x = element_text(hjust = 0.5),
        strip.text = element_text(size = 12),
        legend.position = "bottom",
        legend.key.width = unit(2.5, "cm"),
        legend.key.height = unit(0.25, "cm")
      )
  })

  output$plot_margalef <- renderPlot({
    margalef_data <- data() %>%
      as.data.frame() %>%
      dplyr::group_by(source_folder) %>%
      dplyr::summarise(
        margalef_index = calculate_margalef(cur_data()),
        unique_peptides = length(unique(peptide)),
        total_psms = n(),
        .groups = "drop"
      )
  
  margalef_data %>%
    ggplot(aes(x = source_folder, y = margalef_index)) +
    geom_col(fill = input$plot_color, color = "black") +
    geom_text(
      aes(label = paste0("D = ", round(margalef_index, 2))),
      vjust = -0.5,
      size = 4,
      fontface = "bold"
    ) +
    labs(
      x = NULL,
      y = "Margalef's Diversity Index",
      title = "Peptide Diversity by Sample",
      caption = "Margalef's index: D = (S-1)/ln(N), where S = unique peptides, N = total PSMs"
    ) +
    theme(
      text = element_text(size = 15, color = "black"),
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
      axis.text.y = element_text(hjust = 0.5, size = 10),
      plot.title = element_text(hjust = 0.5, face = "bold"),
      plot.caption = element_text(hjust = 0.5, size = 10)
    )
  })

  multiple_protein_data <- eventReactive(input$load_protein_files, {
    req(input$protein_data_directory)
    if (!dir.exists(input$protein_data_directory)) {
      showNotification("Directory does not exist!", type = "error")
      return(NULL)
    }
    protein_files <- list.files(
      path = input$protein_data_directory,
      pattern = "^protein\\.tsv$",
      recursive = TRUE,
      full.names = TRUE
    )

    if (length(protein_files) == 0) {
      showNotification("No protein files found!", type = "warning")
      return(NULL)
    }

    showNotification(
      paste("Found", length(protein_files), "protein files"),
      type = "message"
    )
    all_data <- NULL

    for (file_path in protein_files) {
      tryCatch(
        {
          subfolder_name <- basename(dirname(file_path))
          protein_data <- readr::read_tsv(file_path, show_col_types = FALSE)
          protein_data <- janitor::clean_names(protein_data)
          protein_data$source_folder <- subfolder_name
          if (is.null(all_data)) {
            all_data <- protein_data
          } else {
            all_data <- dplyr::bind_rows(all_data, protein_data)
          }
        },
        error = function(e) {
          showNotification(
            paste("Error reading", basename(file_path), ":", e$message),
            type = "error"
          )
        }
      )
    }

    return(all_data)
  })
  output$protein_files_status <- renderText({
    if (input$load_protein_files == 0) {
      return("Click 'Load protein Files' to search for files")
    }

    multi_data <- multiple_protein_data()
    if (!is.null(multi_data)) {
      n_folders <- length(unique(multi_data$source_folder))
      folders <- paste(unique(multi_data$source_folder), collapse = ", ")
      paste(
        "Loaded",
        nrow(multi_data),
        "proteins from",
        n_folders,
        "folders:\n",
        folders
      )
    } else {
      "No protein files loaded"
    }
  })
  data_protein <- reactive({
    protein_file <- NULL
    if (input$load_protein_files > 0) {
      multi_data <- multiple_protein_data()
      if (!is.null(multi_data)) {
        protein_file <- multi_data
      }
    }
    if (is.null(protein_file) && !is.null(input$protein)) {
      protein_file <- readr::read_tsv(
        input$protein$datapath,
        show_col_types = FALSE
      ) %>%
        janitor::clean_names() %>%
        mutate(source_folder = "single_upload")
    }
    if (is.null(protein_file)) {
      return(NULL)
    }

    return(protein_file)
  })

  output$info_box2 <- renderInfoBox({
    infoBox(
      "protein.tsv files contain FDR-filtered protein results, where each row is an identified protein group",
      icon = icon("info"),
      color = "black"
    )
  })

  output$plot14 <- renderPlot({
    data_protein() %>%
      as.data.frame() %>%
      ggplot() +
      geom_histogram(
        aes(x = coverage),
        fill = input$plot_color,
        color = "black"
      ) +
      labs(x = "Protein coverage (%)", y = "Count") +
      facet_wrap(~source_folder, scales = "free_y")
  })

  output$plot15 <- renderPlot({
    data_protein() %>%
      as.data.frame() %>%
      dplyr::count(source_folder, organism) %>%
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
      facet_wrap(~source_folder) +
      theme(
        axis.text.y = element_text(face = "italic"),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
      )
  })

  output$plot16 <- renderPlot({
    data_protein() %>%
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
        vjust = 0.5,
        size = 6,
        fontface = "bold"
      ) +
      labs(y = NULL, x = "Count", fill = NULL) +
      facet_wrap(~source_folder) +
      theme(
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
      )
  })

  output$plot17 <- renderPlot({
    data_protein() %>%
      as.data.frame() %>%
      ggplot() +
      geom_histogram(
        aes(x = protein_probability),
        fill = input$plot_color,
        color = "black"
      ) +
      labs(x = "Protein Probability", y = "Count") +
      facet_wrap(~source_folder, scales = "free_y") +
      theme(
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
      )
  })

  output$plot18 <- renderPlot({
    data_protein() %>%
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
      ) +
      facet_wrap(~source_folder, scales = "free_y") +
      theme(
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
      )
  })

  output$plot19 <- renderPlot({
    data_protein() %>%
      as.data.frame() %>%
      ggplot() +
      geom_histogram(
        aes(x = total_peptides),
        fill = input$plot_color,
        color = "black"
      ) +
      labs(x = "Total peptides mapped to proteins", y = "Count") +
      facet_wrap(~source_folder, scales = "free_y") +
      theme(
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
      )
  })

  output$plot20 <- renderPlot({
    data_protein() %>%
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
      ) +
      facet_wrap(~source_folder, scales = "free_y") +
      theme(
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
      )
  })

  output$plot21 <- renderPlot({
    data_protein() %>%
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
      ) +
      facet_wrap(~source_folder, scales = "free_y") +
      theme(
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
      )
  })

  output$plot22 <- renderPlot({
    data_protein() %>%
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
      labs(x = "log2 of Razor Intensity", y = "Protein") +
      facet_wrap(~source_folder, scales = "free_y") +
      theme(
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
      )
  })

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

  output$plot27 <- renderPlotly({
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
