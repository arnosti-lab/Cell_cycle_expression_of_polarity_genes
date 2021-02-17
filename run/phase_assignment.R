################
library(Hmisc)
library(dplyr)
library(tidyr)
library(biomaRt)
library(Seurat)
################

setwd("./")
load("normalized_data/scRNA_dros_embyro_asinh.Rdata")
load("normalized_data/scRNA_dros_wing_disc_asinh.Rdata")
dir.create("results")

polarity_genes = read.delim("gene_files/polarity_genes_from_sandhya.txt", header=F)
polarity_genes_20200209 = read.delim("gene_files/20200209_polarity_genes_from_sandhya.txt", header=F)
polarity_genes = unique(rbind(polarity_genes, polarity_genes_20200209))

cc.sandhya = read.delim("gene_files/cell_cycle_genes_from_sandhya.txt", header =T)

# Organize cc genes and polarity genes ------------------------------------
# change Seurat cc.genes to fly homologs
s.genes = cc.genes[[1]]
g2m.genes = cc.genes[[2]]

listMarts()
fly.ensembl = useMart("ENSEMBL_MART_ENSEMBL", dataset="dmelanogaster_gene_ensembl")
human.ensembl = useMart("ENSEMBL_MART_ENSEMBL", dataset= "hsapiens_gene_ensembl")
head(listAttributes(fly.ensembl, page="feature_page"))

s.fly = getLDS(attributes = "hgnc_symbol", 
               filters =  "hgnc_symbol",
               values = s.genes,
               mart = human.ensembl,
               attributesL = "external_gene_name",
               martL = fly.ensembl)

g2m.fly = getLDS(attributes = "hgnc_symbol", 
                 filters =  "hgnc_symbol",
                 values = g2m.genes,
                 mart = human.ensembl,
                 attributesL = "external_gene_name",
                 martL = fly.ensembl)

save(s.fly, g2m.fly, file = "gene_files/Seurat_ccgenes_fly_homologs.rdata")

g2m.sandhya = as.character(cc.sandhya$Gene[which(cc.sandhya$cc.phase == "G2M.phase")])
g2m.both = unique(c(g2m.sandhya, g2m.fly$Gene.name))

s.sandhya = as.character(cc.sandhya$Gene[which(cc.sandhya$cc.phase == "S.phase")])
s.both = unique(c(s.sandhya, s.fly$Gene.name))

# make a df with all of the genes of interest
s.both = as.data.frame(s.both)
colnames(s.both) = "Gene"
s.both$Definition = rep("s.phase", nrow(s.both))

g2m.both = as.data.frame(g2m.both)
colnames(g2m.both) = "Gene"
g2m.both$Definition = rep("g2m.phase", nrow(g2m.both))

g1.arrest = cc.sandhya[which(cc.sandhya$cc.phase == "G1.arrest"),]
g1.arrest = g1.arrest[c(1)]
g1.arrest$Definition = rep("g1.arrest", nrow(g1.arrest))

colnames(polarity_genes) = "Gene"
polarity_genes$Definition = rep("polarity", nrow(polarity_genes))

goi_df = rbind(s.both, g2m.both, g1.arrest, polarity_genes) 

#Check which genes are expressed in which sc data sets

goi_df$Embryo = goi_df$Gene %in% rownames(embryo_asinh)
goi_df$Wing_Disc = goi_df$Gene %in% rownames(wing_asinh)

# "CycB" is is listed as an s.phase and a G2M phase gene. Remove it.
goi_df = goi_df[!grepl("CycB", goi_df$Gene),]
write.csv(goi_df, "results/20200922_genes_of_interest_detected_in_dataset.csv")


# Label cells with cell-cycle phase ---------------------------------------

scData_list = list(embryo_asinh, wing_asinh)
scData_list = lapply(scData_list, as.matrix)
names(scData_list) = c("embryo", "wing_disc")

for (i in 1:length(scData_list)){
  seurat_dat = CreateSeuratObject(scData_list[[i]], 
                                  project= paste(names(scData_list[i])))
  
  seurat_dat[["RNA"]]@data = scData_list[[i]]
  
  meta = goi_df[which(goi_df[,i+2] == "TRUE"),]
  s.fly = as.character(meta$Gene[which(meta$Definition == "s.phase")])
  g2m.fly = as.character(meta$Gene[which(meta$Definition == "g2m.phase")])
  seurat_dat = CellCycleScoring(seurat_dat, s.fly, g2m.fly, set.ident = TRUE)
  save(seurat_dat, file = paste("results/drosophia", 
                                names(scData_list[i]), 
                                "seurat_object.rdata",
                                sep = "_"))
}


  



