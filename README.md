# PSM Analyst for FragPipe (PSManalyst)

If you use PSManalyst, please consider citing the following publication:
> Chaves AFA. PSManalyst: A Dashboard for Visual Quality Control of FragPipe Results. J Proteome Res. 2025 Sep 5;24(9):4344-4346. [doi: 10.1021/acs.jproteome.5c00557.](https://pubs.acs.org/doi/10.1021/acs.jproteome.5c00557) Epub 2025 Aug 15. PMID: 40815682.

## Dashboard for PSM and protein information visualization from FragPipe search

This is a shiny application that takes the psm.tsv, protein.tsv and combined_protein.tsv files from the FragPipe search and renders them into graphs for quick visualization of your results. To switch between PSM and Protein visualization, you just need to click on the **PSM viewer** or **Protein viewer** buttons. In the **PSM viewer** panel you have the option to filter the PSMs by hyperscore and peptideProphet probability scores.

There are two versions of PSManalyst shiny app:
- PSManalyst: this is the single psm.tsv evaluation version (some analysis are only available in this version, e.g. PICS fingerprint).
- PSManalyst_MB: this is the multi-batch psm.tsv evaluation version (you only need the path to the subfolders e all the psm and protein files will be automatically loaded).

Please, see the [FragPipe](https://fragpipe.nesvilab.org/docs/tutorial_fragpipe_outputs.html) documentation to have a complete understanding of the outputs.

PSManalyst have the instruction to install all the required libraries at the first use:

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

FragPipe is recognized as one of the fastest computational platforms in proteomics, making it a practical solution for the rapid quality control of high-throughput sample analyses. Starting with version 23.0, FragPipe introduced the “Generate Summary Report” feature, offering .pdf reports with essential quality control metrics to address the challenge of intuitively assessing large-scale proteomics data. While traditional spreadsheet formats (e.g., tsv files) are accessible, the complexity of the data often limits user-friendly interpretation. To further enhance accessibility, PSManalyst, a Shiny-based R application, was developed to process FragPipe output files (psm.tsv, protein.tsv, and combined_protein.tsv) and provide interactive, code-free data visualization. Users can filter peptide-spectrum matches (PSMs) by quality scores, visualize protease cleavage fingerprints as heatmaps and SeqLogos, and access a range of quality control metrics and representations such as peptide length distributions, ion densities, mass errors, and wordclouds for overrepresented peptides. The tool facilitates seamless switching between PSM and protein data visualization, offering insights into protein abundance discrepancies, samplewise similarity metrics, protein coverage, and contaminants evaluation. PSManalyst leverages several R libraries (lsa, vegan, ggfortify, ggseqlogo, wordcloud2, tidyverse, ggpointdensity, and plotly) and runs on Windows, MacOS, and Linux, requiring only a local R setup and an IDE.

<p align = "center">
<img src = "https://github.com/41ison/PSManalyst/blob/main/Screenshot_PSManalyst.png" width = "1000">
</p>
