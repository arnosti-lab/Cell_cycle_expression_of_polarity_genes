################
library(tidyr)
library(dplyr)
################

setwd("./")

# Karaiskos et al. 2017, Drosophila embryo  -------------------------------

embryo_raw = read.delim("Raw_data/GSE95025_high_quality_cells_digital_expression.txt", header = TRUE)

# CMP normalzize 
# code modified from Zheng PBMC 
CPM_normalize <-function(x) {
  cs <- colSums(x)
  x_use_genes <- which(cs >= 1)
  x_filt<-x[,x_use_genes]
  scale_factor = 1e6 # like in Seurat's LogNormalize
  x_norm <- t(t(x)/cs) # divide by total UMI per cell
  x_norm * scale_factor
}

embryo_norm = CPM_normalize(embryo_raw)
# hyperbolic arcsine transform

embryo_asinh = asinh(embryo_norm)
dir.create("normalized_data")
save(embryo_asinh, file = "normalized_data/scRNA_dros_embyro_asinh.Rdata")

# Deng et al. 2019, Drosophila wing disc -----------------------------------

wing_orig = read.csv("Raw_data/GSM3902311_frt82b_normalization_data.csv", row.names = 1)

# asinh instead of ln +1 normalize
wing_asinh = asinh(exp(wing_orig)-1)
save(wing_asinh, file = "normalized_data/scRNA_dros_wing_disc_asinh.Rdata")