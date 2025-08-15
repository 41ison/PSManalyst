# PSM Analyst for FragPipe (PSManalyst)

If you use PSManalyst, please consider citing the following publication:
> Chaves, AFA. PSManalyst: A Dashboard for Visual Quality Control of FragPipe Results. Journal of Proteome Research. 2025 [https://doi.org/10.1021/acs.jproteome.5c00557](https://pubs.acs.org/doi/10.1021/acs.jproteome.5c00557).

## Dashboard for PSM and protein information visualization from FragPipe search

This is a shiny application that takes the psm.tsv, protein.tsv and combined_protein.tsv files from the FragPipe search and renders them into graphs for quick visualization of your results. To switch between PSM and Protein visualization, you just need to click on the **PSM viewer** or **Protein viewer** buttons. In the **PSM viewer** panel you have the option to filter the PSMs by hyperscore and peptideProphet probability scores.

There are two versions of PSManalyst shiny app:
- PSManalyst: this is the single psm.tsv evaluation version (some analysis are only available in this version, e.g. PICS fingerprint).
- PSManalyst_MB: this is the multi-batch psm.tsv evaluation version (you only need the path to the subfolders e all the psm and protein files will be automatically loaded).

Please, see the [FragPipe](https://fragpipe.nesvilab.org/docs/tutorial_fragpipe_outputs.html) documentation to have a complete understanding of the outputs.

PSManalyst have the instruction to install all required libraries at the first use:

```r
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
```

Each visualization focuses on a specific aspect of the proteomics data:

- **Protease fingerprint:** Heatmap of amino acid frequencies at cleavage sites
- **Peptide length distribution:** Density plot of peptide lengths
- **N/C-termini SeqLogo:** Sequence logo visualizations
- **GRAVY index:** Grand average of hydropathicity for peptide sequences
- **Isoelectric point:** Isoelectric point distribution for peptide sequences
- **Charge state distribution:** Bar chart of charge states
- **m/z over retention time:** Point density plot
- **Missed cleavages:** Bar chart of missed cleavage counts
- **Mass error:** Scatter plot of mass errors
- **Score distributions:** Histograms of various scoring metrics
- **Protein coverage:** Histogram of protein coverage percentages
- **Sample correlation:** Scatter plot of intensity correlations between samples
- **Similarity metrics:** Heatmaps of various similarity metrics between samples

<p align = "center">
<img src = "https://github.com/41ison/PSManalyst/blob/main/Screenshot_PSManalyst.png" width = "1000">
</p>
