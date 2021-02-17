## Integrating growth and architecture: Cell cycle expression of polarity genes
Sandhya Payankaulam , Stephanie L. Hickey, and David N. Arnosti
##

This repository contains the results generated in this study as well as the code necessary for reproducing all of the analyses in figure 1 and the supplemental figures.

### Data
The WT Drosophila Wing imaginal disc 10x scRNA-seq data generated and described in Deng *et. al.* 2019 can be fond [here](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM3902311). The early Drosophila embryo Drop-seq data generated and described in Karaiskos et. al. 2017 can be found [here](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi). 

**Downloading raw data**

The specific data used in this study can be downloaded from GEO by running:

`$ sh ./run/download_data.sh`

`download_data.sh` creates a `Raw_data` to contain the downloaded files.

### Analysis

The following analyses were completed using R version 4.0.0. To install the necessary packages [install R](https://www.r-project.org/) and run:

**CAUTION** this will install Seurat V3.2.3, replacing any other version of Seurat in your R library
`$ Rscript ./run/install_packages.R`

**Preprocessing**

The data was preprocessed and normalized using:

`$ Rscript ./run/preprocess.R`

`preprocess.R` creates a `normalized_data` directory which contains R objects with normalized gene expression matrices.

**Phase assignment**

Each cell was assigned a putative cell-cycle phase label using [Seurat V3's](https://satijalab.org/seurat/) *CellCycleScoring* function. We used a custom list of cell-cycle phase marker genes, and the citations linking these genes to specific cell-cycle phases can be found in supplemental table 1 included with the manuscript. To recreated our labels run:

`$ Rscript ./run/phase_assignment.R`

`phase_assignment.R` creates a results folder and outputs:
- `results/20200922_genes_of_interest_detected_in_dataset.csv`, a table that lists the labels for polarity and cell-cycle genes and whether or not they are expressed in a particular single-cell data set. 
- Seurat objects that contain the phase assignments for each cell in the @meta.data slot

**Figures and Tables**

All of the tables and figures from this manuscript are available in `results` and can be recreated by running:

`$ Rscript ./run/figure1_plots.R`

- `figure1_plots.R` outputs the contents of `results/figure1`
- these plots were used in figure 1 (wing disc data) and figure S5 (embryo data)

and 

`$ Rscript ./run/suplemental_figures.R`

`suplemental_figures.R` outputs the contents of `results/supplemental`
- `*_fisher_pval.csv` and `*_fisher_FDR.csv` list the p-values or FDRs associated with the stacked bar plots (Figures 1C and E, S1, S3, S5C and E, S6, S8)
- `*_kw_dunns.csv` lists the p-values and FDRs associated with the expression vs. phase dot plots and heatmaps (Figures 1B, D and F, S2, S4, S5B, D and F, S7, S9)

