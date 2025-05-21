# PSM Analyst for FragPipe (PSManalyst)

## Dashboard for PSM and protein information visualization from FragPipe search

This is a shiny application that takes the psm.tsv, protein.tsv and combined_protein.tsv files from the FragPipe search and renders them into graphs for quick visualization of your results. To switch between PSM and Protein visualization, you just need to click on the **PSM viewer** or **Protein viewer** buttons. In the **PSM viewer** panel you have the option to filter the PSMs by hyperscore. Please, see the [FragPipe](https://fragpipe.nesvilab.org/docs/tutorial_fragpipe_outputs.html) documentation to have a complete understanding of the outputs.

Important: You need the following libraries in order to run the App:

```r
library(shiny)            # from CRAN
library(shinydashboard)   # from CRAN
library(tidyverse)        # from CRAN
library(janitor)          # from CRAN
library(ggseqlogo)        # from CRAN
library(ggpointdensity)   # from GitHub
library(wordcloud2)       # from GitHub
library(ggtext)           # from CRAN
```

Each visualization focuses on a specific aspect of the proteomics data:

**Protease fingerprint:** Heatmap of amino acid frequencies at cleavage sites
**Peptide length distribution:** Density plot of peptide lengths
**N/C-termini SeqLogo:** Sequence logo visualizations
**Charge state distribution:** Bar chart of charge states
**m/z over retention time:** Point density plot
**Missed cleavages:** Bar chart of missed cleavage counts
**Mass error:** Scatter plot of mass errors
**Score distributions:** Histograms of various scoring metrics
**Protein coverage:** Histogram of protein coverage percentages
**Sample correlation:** Scatter plot of intensity correlations between samples
**Similarity metrics:** Heatmaps of various similarity metrics between samples

<p align = "center">
<img src = "https://github.com/41ison/PSManalyst/blob/main/Screenshot%20PSManalyst.png" width = "1000">
</p>
