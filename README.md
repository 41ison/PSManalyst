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

<p align = "center">
<img src = "https://github.com/41ison/PSManalyst/blob/main/Screenshot%20of%20the%20App.png" width = "1000">
</p>
