## The psm.tsv, protein.tsv and combined_protein.tsv files are the inputs for the PSManalyst dashboard
## The user can filter PSMs by the hyperscore and PeptideProphet probability, as well as by enzymatic specificity
## The user can remove a contaminant organism as well
## Chaves AFA. PSManalyst: A Dashboard for Visual Quality Control of FragPipe Results. J Proteome Res. 2025 Sep 5;24(9):4344-4346. doi: 10.1021/acs.jproteome.5c00557. Epub 2025 Aug 15. PMID: 40815682.

CRAN_packages <- c(
  "shiny", "shinydashboard", "tidyverse", "janitor", "ggseqlogo", "ggtext", "lsa", "vegan", "plotly", "viridis", "ggfortify", "colourpicker"
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
    'A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y'
  )

  for (amino_acid in twenty_amino_acids) {
    if (!amino_acid %in% names(element)) {
      element[[amino_acid]] <- 0
    }
  }
  element <- element[match(twenty_amino_acids, names(element))]

  return(element)
}

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
    "P4", "P3", "P2", "P1", "P1'", "P2'", "P3'", "P4'"
  )
  rownames(final_matrix) <- c(
    "A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y"
  )

  return(final_matrix)
}

# GRAVY (Grand Average of Hydropathy)
# Ref: Kyte J, Doolittle RF. A simple method for displaying the hydropathic character of a protein. J Mol Biol. 1982 May 5;157(1):105-32. doi: 10.1016/0022-2836(82)90515-0
GRAVY <- function(sequence) {
  hydropathy_index <- c(
    A = 1.8, R = -4.5, N = -3.5, D = -3.5, C = 2.5, Q = -3.5, E = -3.5, G = -0.4, H = -3.2, I = 4.5, L = 3.8, K = -3.9, M = 1.9, F = 2.8, P = -1.6, S = -0.8, T = -0.7, W = -0.9, Y = -1.3, V = 4.2
  )
  scores <- sapply(strsplit(sequence, NULL)[[1]], function(aa) {
    hydropathy_index[aa]
  })
  return(mean(scores, na.rm = TRUE))
}

calculate_pI <- function(sequence) {
  pKa_values <- c(A = 2.34, R = 12.48, N = 10.76, D = 3.86, C = 8.33, Q = 10.76, E = 4.25, G = 2.34, H = 6.00, I = 6.04, L = 6.04, K = 9.74, M = 5.74, F = 5.48, P = 1.99, S = 2.21, T = 2.15, W = 9.39, Y = 10.07, V = 6.02)

  pI <- mean(
    sapply(strsplit(sequence, NULL)[[1]], function(aa) pKa_values[aa]),
    na.rm = TRUE
  )

  return(pI)
}

color_blue_seq <- c("#d4e6f1", "#a9cce3", "#7fb3d5", "#5499c7", "#2980b9", "#1f618d", "#154360")

# Function to calculate amino acid co-occurrence matrix
analyze_terminus_cooccurrence <- function(
  df,
  peptide_col,
  show_values = FALSE
) {
  amino_acids <- c("A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y")
  terminus_data <- df %>%
    dplyr::filter(
      !is.na(.data[[peptide_col]]) & nchar(.data[[peptide_col]]) > 0
    ) %>%
    dplyr::mutate(
      N_terminus = substr(.data[[peptide_col]], 1, 1),
      C_terminus = substr(
        .data[[peptide_col]],
        nchar(.data[[peptide_col]]),
        nchar(.data[[peptide_col]])
      )
    ) %>%
    dplyr::filter(N_terminus %in% amino_acids & C_terminus %in% amino_acids)
  contingency_table <- table(terminus_data$N_terminus, terminus_data$C_terminus)
  total_peptides <- sum(contingency_table)
  prob_matrix <- contingency_table / total_peptides
  full_matrix <- matrix(
    0,
    nrow = length(amino_acids),
    ncol = length(amino_acids),
    dimnames = list(amino_acids, amino_acids)
  )

  for (i in rownames(contingency_table)) {
    for (j in colnames(contingency_table)) {
      full_matrix[i, j] <- prob_matrix[i, j]
    }
  }

  heatmap_data <- expand.grid(
    N_terminus = rownames(full_matrix),
    C_terminus = colnames(full_matrix)
  ) %>%
    mutate(
      Probability = as.vector(full_matrix),
      N_terminus = factor(N_terminus, levels = rownames(full_matrix)),
      C_terminus = factor(C_terminus, levels = colnames(full_matrix))
    )

  p <- ggplot(
    heatmap_data,
    aes(x = C_terminus, y = N_terminus, fill = Probability)
  ) +
    geom_tile(color = "white", size = 0.1) +
    scale_fill_gradient2(
      low = "white",
      mid = "#B8E6D9",
      high = "#1BB99A",
      midpoint = max(heatmap_data$Probability) / 2,
      name = "Probability"
    ) +
    labs(
      title = "N:C-terminus co-occurrence",
      x = "C-terminus AA",
      y = "N-terminus AA"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 14, hjust = 0.5, face = "bold"),
      axis.text.x = element_text(
        hjust = 0.5,
        size = 10,
        face = "bold",
        color = "black"
      ),
      axis.text.y = element_text(size = 10, face = "bold", color = "black"),
      axis.title = element_text(size = 12, face = "bold"),
      legend.title = element_text(size = 10, face = "bold", hjust = 0.5),
      legend.title.position = "top",
      legend.text = element_text(size = 10, face = "bold"),
      legend.key.height = unit(0.25, "cm"),
      legend.key.width = unit(1, "cm"),
      legend.position = "bottom",
      panel.border = element_rect(color = "black", fill = NA)
    ) +
    coord_fixed()

  if (show_values) {
    p <- p +
      geom_text(
        aes(label = sprintf("%.3f", Probability)),
        size = 2,
        color = "black"
      )
  }

  return(list(
    matrix = full_matrix,
    plot = p,
    total_peptides = total_peptides,
    summary_stats = list(
      min_prob = min(full_matrix[full_matrix > 0]),
      max_prob = max(full_matrix),
      mean_prob = mean(full_matrix[full_matrix > 0]),
      n_observed_pairs = sum(full_matrix > 0)
    )
  ))
}

analyze_modification_type <- function(
  data,
  rt_threshold_insource = 0.5,
  rt_threshold_real = 1.5,
  min_observations = 3
) {
  mod_patterns <- list(
    H2O_loss = "-18\\.0106",
    NH3_loss = "-17\\.0265",
    Oxidation = "15\\.9949",
    Deamidation = "0\\.9840"
  )
  modified_psms <- data %>%
    dplyr::filter(str_detect(
      assigned_modifications,
      paste(unlist(mod_patterns), collapse = "|")
    )) %>%
    dplyr::mutate(
      mod_type = case_when(
        str_detect(assigned_modifications, mod_patterns$H2O_loss) ~ "H2O_loss",
        str_detect(assigned_modifications, mod_patterns$NH3_loss) ~ "NH3_loss",
        str_detect(
          assigned_modifications,
          mod_patterns$Oxidation
        ) ~ "Oxidation",
        str_detect(
          assigned_modifications,
          mod_patterns$Deamidation
        ) ~ "Deamidation",
        TRUE ~ "other"
      ),
      retention_min = retention / 60,
      base_peptide = str_remove_all(peptide, "\\[.*?\\]")
    ) %>%
    dplyr::filter(mod_type != "other")

  unmodified_psms <- data %>%
    dplyr::filter(
      assigned_modifications == "" | is.na(assigned_modifications)
    ) %>%
    dplyr::mutate(
      retention_min = retention / 60,
      base_peptide = str_remove_all(peptide, "\\[.*?\\]")
    ) %>%
    dplyr::select(base_peptide, retention_min, intensity) %>%
    rename(unmod_rt = retention_min, unmod_intensity = intensity)

  modified_stats <- modified_psms %>%
    group_by(base_peptide, mod_type) %>%
    summarise(
      n_modified = n(),
      median_mod_rt = median(retention_min),
      mean_mod_rt = mean(retention_min),
      sd_mod_rt = sd(retention_min),
      rt_range_mod = max(retention_min) - min(retention_min),
      cv_mod_rt = sd(retention_min) / mean(retention_min) * 100,
      mean_intensity_mod = mean(intensity, na.rm = TRUE),
      .groups = "drop"
    )

  unmodified_stats <- unmodified_psms %>%
    group_by(base_peptide) %>%
    summarise(
      n_unmodified = n(),
      median_unmod_rt = median(unmod_rt),
      mean_unmod_rt = mean(unmod_rt),
      sd_unmod_rt = sd(unmod_rt),
      rt_range_unmod = max(unmod_rt) - min(unmod_rt),
      cv_unmod_rt = sd(unmod_rt) / mean(unmod_rt) * 100,
      mean_intensity_unmod = mean(unmod_intensity, na.rm = TRUE),
      .groups = "drop"
    )

  rt_shift_analysis <- modified_stats %>%
    inner_join(unmodified_stats, by = c("base_peptide")) %>%
    dplyr::filter(
      n_modified >= min_observations,
      n_unmodified >= min_observations
    ) %>%
    dplyr::mutate(
      rt_shift = median_mod_rt - median_unmod_rt,
      abs_rt_shift = abs(rt_shift),
      rel_rt_shift = (rt_shift / median_unmod_rt) * 100,
      total_rt_variability = rt_range_mod + rt_range_unmod,
      shift_to_variability_ratio = abs_rt_shift / (total_rt_variability + 0.01),
      intensity_ratio = mean_intensity_mod / mean_intensity_unmod,
      classification = case_when(
        abs_rt_shift <= rt_threshold_insource &
          shift_to_variability_ratio < 2 ~ "Likely in-source loss",
        abs_rt_shift >= rt_threshold_real &
          shift_to_variability_ratio > 3 ~ "Likely real modification",
        TRUE ~ "Ambiguous"
      ),
      expected_behavior = case_when(
        mod_type %in% c("H2O_loss", "NH3_loss") ~ "In-source loss",
        mod_type %in% c("Oxidation", "Deamidation") ~ "Real modification",
        TRUE ~ "Unknown"
      ),
      classification_agreement = classification ==
        paste("Likely", expected_behavior)
    )

  summary_by_mod <- rt_shift_analysis %>%
    group_by(mod_type, classification) %>%
    summarise(
      n_peptides = n(),
      median_abs_shift = median(abs_rt_shift),
      mean_abs_shift = mean(abs_rt_shift),
      median_rel_shift = median(abs(rel_rt_shift)),
      median_shift_var_ratio = median(shift_to_variability_ratio),
      median_intensity_ratio = median(intensity_ratio),
      .groups = "drop"
    ) %>%
    arrange(mod_type, classification)

  p1 <- ggplot(
    rt_shift_analysis,
    aes(x = mod_type, y = abs_rt_shift, fill = classification)
  ) +
    geom_boxplot(outlier.alpha = 0.3) +
    geom_hline(
      yintercept = rt_threshold_insource,
      linetype = "dashed",
      color = "blue",
      linewidth = 0.8
    ) +
    geom_hline(
      yintercept = rt_threshold_real,
      linetype = "dashed",
      color = "red",
      linewidth = 0.8
    ) +
    scale_fill_manual(
      values = c(
        "Likely in-source loss" = "#4CAF50",
        "Likely real modification" = "#F44336",
        "Ambiguous" = "#FFC107"
      )
    ) +
    labs(
      title = "RT shifts: modified vs unmodified peptides",
      caption = paste(
        "Blue line: in-source threshold (",
        rt_threshold_insource,
        " min), ",
        "Red line: real modification threshold (",
        rt_threshold_real,
        " min)",
        sep = ""
      ),
      x = "Modification type",
      y = "Absolute RT shift (min)",
      fill = "Classification"
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(size = 14, hjust = 0.5, face = "bold"),
      axis.text.x = element_text(
        size = 10,
        face = "bold",
        color = "black",
        angle = 45,
        hjust = 1
      ),
      axis.text.y = element_text(size = 10, face = "bold", color = "black"),
      axis.title = element_text(size = 12, face = "bold"),
      legend.title = element_text(size = 10, face = "bold", hjust = 0.5),
      legend.title.position = "top",
      legend.text = element_text(size = 10, face = "bold"),
      legend.position = "bottom",
      panel.border = element_rect(color = "black", fill = NA)
    )

  p2 <- ggplot(
    rt_shift_analysis,
    aes(
      x = total_rt_variability,
      y = abs_rt_shift,
      color = mod_type,
      shape = classification
    )
  ) +
    geom_point(size = 3, alpha = 0.6) +
    geom_abline(slope = 2, intercept = 0, linetype = "dashed", color = "blue") +
    geom_abline(slope = 3, intercept = 0, linetype = "dashed", color = "red") +
    scale_shape_manual(
      values = c(
        "Likely in-source loss" = 16,
        "Likely real modification" = 17,
        "Ambiguous" = 15
      )
    ) +
    scale_color_manual(
      values = c(
        "H2O_loss" = "#1f77b4",
        "NH3_loss" = "#ff7f0e",
        "Oxidation" = "#2ca02c",
        "Deamidation" = "#d62728"
      )
    ) +
    labs(
      title = "RT shift vs total RT variability",
      caption = "Ratios > 3 suggest real modifications; < 2 suggest in-source losses",
      x = "Total RT variability (min)",
      y = "Absolute RT shift (min)",
      color = "Modification type",
      shape = "Classification"
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(size = 14, hjust = 0.5, face = "bold"),
      axis.text = element_text(size = 10, face = "bold", color = "black"),
      axis.title = element_text(size = 12, face = "bold"),
      legend.title = element_text(size = 10, face = "bold", hjust = 0.5),
      legend.title.position = "top",
      legend.text = element_text(size = 10, face = "bold"),
      legend.position = "bottom",
      panel.border = element_rect(color = "black", fill = NA)
    )

  p3 <- ggplot(
    rt_shift_analysis,
    aes(
      x = round(log10(intensity_ratio), 1),
      y = abs_rt_shift,
      color = classification
    )
  ) +
    geom_point(size = 2, alpha = 0.6) +
    facet_wrap(~mod_type, scales = "free") +
    geom_hline(
      yintercept = rt_threshold_insource,
      linetype = "dashed",
      alpha = 0.5
    ) +
    geom_hline(
      yintercept = rt_threshold_real,
      linetype = "dashed",
      alpha = 0.5
    ) +
    scale_color_manual(
      values = c(
        "Likely in-source Loss" = "#4CAF50",
        "Likely real modification" = "#F44336",
        "Ambiguous" = "#FFC107"
      )
    ) +
    labs(
      title = "RT shift vs intensity ratio by modification type",
      x = "log<sub>10</sub>(intensity ratio: modified/unmodified)",
      y = "Absolute RT shift (minutes)",
      color = "Classification"
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(size = 14, hjust = 0.5, face = "bold"),
      axis.text = element_text(size = 10, face = "bold", color = "black"),
      axis.title = element_markdown(size = 12, face = "bold"),
      legend.title = element_text(size = 10, face = "bold", hjust = 0.5),
      legend.title.position = "top",
      legend.text = element_text(size = 10, face = "bold"),
      legend.position = "bottom",
      panel.border = element_rect(color = "black", fill = NA),
      strip.background = element_blank(),
      strip.text = element_text(size = 12, face = "bold")
    )

  classification_summary <- rt_shift_analysis %>%
    group_by(mod_type, classification, expected_behavior) %>%
    summarise(count = n(), .groups = "drop")

  p4 <- ggplot(
    classification_summary,
    aes(x = mod_type, y = count, fill = classification)
  ) +
    geom_bar(stat = "identity", position = "dodge") +
    geom_text(
      aes(label = count),
      position = position_dodge(width = 0.9),
      vjust = -0.5
    ) +
    scale_fill_manual(
      values = c(
        "Likely in-source loss" = "#4CAF50",
        "Likely real modification" = "#F44336",
        "Ambiguous" = "#FFC107"
      )
    ) +
    labs(
      title = "Classification summary by modification type",
      x = "Modification type",
      y = "Number of peptides",
      fill = "Classification"
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(size = 14, hjust = 0.5, face = "bold"),
      axis.text = element_text(size = 10, face = "bold", color = "black"),
      axis.title = element_markdown(size = 12, face = "bold"),
      legend.title = element_text(size = 10, face = "bold", hjust = 0.5),
      legend.title.position = "top",
      legend.text = element_text(size = 10, face = "bold"),
      legend.position = "bottom",
      panel.border = element_rect(color = "black", fill = NA),
      axis.text.x = element_text(angle = 45, hjust = 1)
    )

  list(
    rt_shift_data = rt_shift_analysis,
    summary_by_modification = summary_by_mod,
    plots = list(
      rt_shift_boxplot = p1,
      shift_vs_variability = p2,
      shift_vs_intensity = p3,
      classification_summary = p4
    )
  )
}

ui <- dashboardPage(
  dashboardHeader(
    title = "PSM Analyst for FragPipe",
    titleWidth = "250",
    dropdownMenu(
      type = "messages",
      messageItem(
        from = "Communication",
        message = tags$a(
          "doi: 10.1021/acs.jproteome.5c00557",
          href = "https://pubs.acs.org/doi/10.1021/acs.jproteome.5c00557",
          target = "_blank"
        ),
        icon = icon("file-lines")
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
          ),
          box(
            title = "Co-occurrence probability matrix of N- and C-terminus AA",
            status = "primary",
            solidHeader = TRUE,
            plotOutput("plot19"),
            collapsible = TRUE
          ),
          box(
            title = "Cysteine counts in peptides",
            status = "primary",
            solidHeader = TRUE,
            plotOutput("plot20"),
            collapsible = TRUE
          ),
          box(
            title = "RT shift distribution by modification type",
            status = "primary",
            solidHeader = TRUE,
            plotOutput("plot21"),
            collapsible = TRUE
          ),
          box(
            title = "Shift-to-variability ratio",
            status = "primary",
            solidHeader = TRUE,
            plotOutput("plot22"),
            collapsible = TRUE
          ),
          box(
            title = "RT shift vs intensity ratio",
            status = "primary",
            solidHeader = TRUE,
            plotOutput("plot23"),
            collapsible = TRUE
          ),
          box(
            title = "Classification summary of modification types",
            status = "primary",
            solidHeader = TRUE,
            plotOutput("plot24"),
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
            plotOutput("plot01p"),
            collapsible = TRUE
          ),
          box(
            title = "Number of proteins by organim",
            status = "primary",
            solidHeader = TRUE,
            plotOutput("plot02p"),
            collapsible = TRUE
          ),
          box(
            title = "Protein existence evidence",
            status = "primary",
            solidHeader = TRUE,
            plotOutput("plot03p"),
            collapsible = TRUE
          ),
          box(
            title = "Protein probability (ProteinProphet)",
            status = "primary",
            solidHeader = TRUE,
            plotOutput("plot04p"),
            collapsible = TRUE
          ),
          box(
            title = "Top Peptide Probability",
            status = "primary",
            solidHeader = TRUE,
            plotOutput("plot05p"),
            collapsible = TRUE
          ),
          box(
            title = "Total peptides mapped to the proteins",
            status = "primary",
            solidHeader = TRUE,
            plotOutput("plot06p"),
            collapsible = TRUE
          ),
          box(
            title = "Razor spectral count",
            status = "primary",
            solidHeader = TRUE,
            plotOutput("plot07p"),
            collapsible = TRUE
          ),
          box(
            title = "Razor intensity",
            status = "primary",
            solidHeader = TRUE,
            plotOutput("plot08p"),
            collapsible = TRUE
          ),
          box(
            title = "Top 20 proteins with higher razor intensity",
            status = "primary",
            solidHeader = TRUE,
            plotOutput("plot09p"),
            collapsible = TRUE
          ),
          box(
            title = "MaxLFQ intensity distribution",
            status = "primary",
            solidHeader = TRUE,
            plotOutput("plot10p"),
            collapsible = TRUE
          ),
          box(
            title = "Sample correlation - Non-normalized log2(Intensity)",
            status = "primary",
            height = 600,
            solidHeader = TRUE,
            plotlyOutput("plot11p"),
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

  cooccurrence_data <- reactive({
    req(data())
    analyze_terminus_cooccurrence(data(), "peptide", show_values = TRUE)
  })

  rt_loss_analysis <- reactive({
    req(data())
    data() %>%
      analyze_modification_type(
        rt_threshold_insource = 0.5,
        rt_threshold_real = 1.5,
        min_observations = 3
      )
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
          c("A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y")
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
        labels = c("1" = "P4", "2" = "P3", "3" = "P2", "4" = "P1", "5" = "P1'", "6" = "P2'", "7" = "P3'", "8" = "P4'")
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
        labels = c("1" = "P4", "2" = "P3", "3" = "P2", "4" = "P1", "5" = "P1'", "6" = "P2'", "7" = "P3'", "8" = "P4'")
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

  output$plot19 <- renderPlot({
    cooccurrence_data()[["plot"]]
  })

  output$plot20 <- renderPlot({
    data() %>%
      as.data.frame() %>%
      dplyr::mutate(cysteine_count = str_count(peptide, "C")) %>%
      dplyr::count(cysteine_count) %>%
      ggplot(aes(x = cysteine_count, y = n)) +
      geom_bar(
        stat = "identity",
        position = "dodge",
        show.legend = FALSE,
        fill = input$plot_color,
        color = "black"
      ) +
      geom_text(aes(label = n), vjust = -0.5, size = 5) +
      labs(x = "Cysteine counts in peptides", y = "Cys count")
  })

  output$plot21 <- renderPlot({
    rt_loss_analysis()[["plots"]][["rt_shift_boxplot"]]
  })

  output$plot22 <- renderPlot({
    rt_loss_analysis()[["plots"]][["shift_vs_variability"]]
  })

  output$plot23 <- renderPlot({
    rt_loss_analysis()[["plots"]][["shift_vs_intensity"]]
  })

  output$plot24 <- renderPlot({
    rt_loss_analysis()[["plots"]][["classification_summary"]]
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
              c("A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y")
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
            labels = c("1" = "P4", "2" = "P3", "3" = "P2", "4" = "P1", "5" = "P1'", "6" = "P2'", "7" = "P3'", "8" = "P4'")
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
            labels = c("1" = "P4", "2" = "P3", "3" = "P2", "4" = "P1", "5" = "P1'", "6" = "P2'", "7" = "P3'", "8" = "P4'")
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

      # Plot 19 - Co-occurrence heatmap
      plots_to_save[["plot19_cooccurrence_heatmap.png"]] <- function() {
        cooccurrence_data()[["plot"]]
      }

      # Plot 20 - Cysteine counts
      plots_to_save[["plot20_cysteine_counts.png"]] <- function() {
        data() %>%
          as.data.frame() %>%
          dplyr::mutate(cysteine_count = str_count(peptide, "C")) %>%
          dplyr::count(cysteine_count) %>%
          ggplot(aes(x = cysteine_count, y = n)) +
          geom_bar(
            stat = "identity",
            position = "dodge",
            show.legend = FALSE,
            fill = input$plot_color,
            color = "black"
          ) +
          geom_text(aes(label = n), vjust = -0.5, size = 5) +
          labs(x = "Cysteine counts in peptides", y = "Cys count")
      }

      # Plot 21 - RT shift boxplot
      plots_to_save[["plot21_rt_shift_boxplot.png"]] <- function() {
        rt_loss_analysis()[["plots"]][["rt_shift_boxplot"]]
      }

      # Plot 22 - Shift vs variability
      plots_to_save[["plot22_shift_vs_variability.png"]] <- function() {
        rt_loss_analysis()[["plots"]][["shift_vs_variability"]]
      }

      # Plot 23 - Shift vs intensity
      plots_to_save[["plot23_shift_vs_intensity.png"]] <- function() {
        rt_loss_analysis()[["plots"]][["shift_vs_intensity"]]
      }

      # Plot 24 - Classification summary
      plots_to_save[["plot24_classification_summary.png"]] <- function() {
        rt_loss_analysis()[["plots"]][["classification_summary"]]
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
  output$plot01p <- renderPlot({
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

  output$plot02p <- renderPlot({
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

  output$plot03p <- renderPlot({
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

  output$plot04p <- renderPlot({
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

  output$plot05p <- renderPlot({
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

  output$plot06p <- renderPlot({
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

  output$plot07p <- renderPlot({
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

  output$plot08p <- renderPlot({
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

  output$plot09p <- renderPlot({
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
      labs(x = "log<sub>2</sub> of Razor Intensity", y = "Protein") +
      theme(
        axis.text = element_markdown()
      )
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

  output$plot10p <- renderPlot({
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

  output$plot11p <- renderPlotly({
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
