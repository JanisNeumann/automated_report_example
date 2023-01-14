# Install CRAN packages
install.packages("tidyverse", "caret", "viridis", "glmnet", "knitr","kableExtra", "htmltools")

# Install Bioconductor packages
bioc_packages <- c("GEOquery", "minfi",
                   "IlluminaHumanMethylation450kmanifest",
                   "IlluminaHumanMethylation450kanno.ilmn12.hg19")

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(bioc_packages)
