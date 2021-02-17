print("**CAUTION** this will install Seurat V3.2.3, replacing any other version of Seurat in your R library")

install.packages(c("tidyverse",
                   "Hmisc",
                   "reshape2",
                   "pheatmap",
                   "RColorBrewer",
                   "patchwork",
                   "FSA",
                   "ggpubr",
                   "rstatix",
                   "deldir"), 
                 repos = "https://repo.miserver.it.umich.edu/cran/")

remotes::install_version("Seurat",
                         version = "3.2.3", 
                         repos = "https://repo.miserver.it.umich.edu/cran/")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("biomaRt")