library(Seurat)
library(dplyr)
library(tidyr)
library(ggplot2)
library(reshape2)
library(pheatmap)
library(RColorBrewer)
library(patchwork)

setwd("./results")
load("drosophia_wing_disc_seurat_object.rdata")
seurat_wing = seurat_dat
load("drosophia_embryo_seurat_object.rdata")
seurat_embryo = seurat_dat
goi_df = read.csv("20200922_genes_of_interest_detected_in_dataset.csv",  row.names = 1)
dir.create("figure1")


seurat_list = list(seurat_embryo, seurat_wing)
names(seurat_list) = c("embryo", "wing_disc")

marker = goi_df[,2]
names(marker) = goi_df$Gene

## heatmaps scaled 

for (i in 1:length(seurat_list)){
  goi_df = goi_df[which(goi_df[,i+2] == TRUE),]
  goi = goi_df$Gene
  meta = seurat_list[[i]]@meta.data
  data = t(GetAssayData(seurat_list[[i]]))
  goi_exp= data[, as.character(goi)]
  phase_df = meta[c(6)]
  for_test = cbind(phase_df, goi_exp)
  
  marker = goi_df[,2]
  names(marker) = goi_df$Gene
  
  for_test$Cell_ID = rownames(for_test)
  
  tidy = melt(for_test, id = c("Cell_ID", "Phase"))
  colnames(tidy) = c("Cell_ID", "Phase", "Gene", "Expression")
  
  for_plot = tidy %>% group_by(Phase,Gene) %>% summarise(mean_exp = mean(Expression))
  for_plot_wide = as.data.frame(pivot_wider(for_plot, names_from = Phase, values_from = mean_exp))
  rownames(for_plot_wide) = for_plot_wide$Gene
  for_plot_wide$Gene = NULL
  for_plot_wide = for_plot_wide[,c(1,3,2)]
  for_plot_wide = for_plot_wide[which(rowSums(for_plot_wide) > 0),]
  marker2 = marker[names(marker) %in% rownames(for_plot_wide)]
  marker_df = as.data.frame(marker2[!is.na(marker2)])
  colnames(marker_df) = "marker"
  
  if("polarity" %in% marker2){
    
    anno = as.data.frame(marker2[which(marker2 == "polarity")])
    if(nrow(anno) > 1) {
      clust = TRUE
    } else {clust = FALSE}
    colnames(anno) = "marker"
    
    pdf(file = paste0("figure1/",
                      names(seurat_list)[i],
                      "_scaled_mean_expression_scaled_polarity_genes_heatmap.pdf"
    ),
    width = 3,
    height = 12)
    
    pheatmap(for_plot_wide[rownames(for_plot_wide) %in% rownames(anno),],
             color = colorRampPalette(rev(brewer.pal(n = 11, name ="RdBu")))(100),
             border_color = NA,
             cellheight = 12,
             cluster_cols = FALSE,
             cluster_rows = clust,
             scale = "row")
    
    dev.off()
    
    print(paste0(names(seurat_list)[i], " polarity genes heatmap finished"))
    
    
  } else {print(paste0("mean expression of all polarity genes is zero in ", names(seurat_list)[i]))}
  
  anno = as.data.frame(marker2[which(!marker2 == "polarity")])
  colnames(anno) = "marker"
  
  ann_colors = list(
    marker = c(g1.arrest = "#58D68D", 
                  s.phase = "#00838F", 
                  g2m.phase = "#FF9800")
  )
  
  for_markers = for_plot_wide[rownames(for_plot_wide) %in% rownames(anno),]
  for_markers = t(for_markers)
  
  pdf(file = paste0("figure1/",
                    names(seurat_list)[i],
                    "_scaled_mean_expression_cell_cycle_genes_pheatmap.pdf"
  ),
  width = 14,
  height = 5)
  
  pheatmap(for_markers,
           color = colorRampPalette(rev(brewer.pal(n = 11, name ="RdBu")))(100),
           border_color = NA,
           cellheight = 12,
           cellwidth = 9,
           fontsize = 10,
           cluster_cols = TRUE,
           cluster_rows = FALSE,
           scale = "column",
           annotation_col = anno,
           annotation_colors = ann_colors
  )
  
  dev.off()
  print(paste0(names(seurat_list)[i], " cc phase genes heatmap finished"))
}



for (i in 1:length(seurat_list)){
  
  goi = c("Vang", "Mcm5")
  meta = seurat_list[[i]]@meta.data
  phase_df = meta[c(6)]
  data = t(GetAssayData(seurat_list[[i]]))
  goi_exp = data[, as.character(goi)]
  for_test = cbind(phase_df, goi_exp)
  for_test$Cell_ID = rownames(for_test)
  
  for_plot = melt(for_test, id = c("Cell_ID", "Phase"))
  colnames(for_plot) = c("Cell_ID", "Phase", "Gene", "Expression")
  for_plot$Phase <- factor(for_plot$Phase,levels = c("G1", "S", "G2M"))
  for_plot$Detected = ifelse(for_plot$Expression > 0,
                             "above detection level",
                             "below detection level")
  
  # proportion bar graphs
  for_stack = 
    for_plot %>% 
    group_by(Phase,Gene,Detected) %>% 
    summarise(nCells = length(Phase))
  
  fisherp = function(x){fisher.test(matrix(x,2,2), alternative = "greater")$p.value} 
  
  fisher_gene = list()
  for(gene in unique(for_stack$Gene)){
    pergene = filter(for_stack, Gene == gene)
    pergene = pergene %>% pivot_wider(names_from = Detected, values_from = nCells)
    colnames(pergene) = c("Phase", "Gene", "phase_above", "phase_below")
    all_above = sum(pergene$phase_above)
    all_below = sum(pergene$phase_below)
    pergene$rest_above = all_above - pergene$phase_above
    pergene$rest_below = all_below - pergene$phase_below
    for_fisher = as.matrix(pergene[,3:6])
    fisher_gene[[gene]] = apply(for_fisher, 1, fisherp)
    names(fisher_gene[[gene]]) = c("G1", "S", "G2M")
  }
  fisher_pval = do.call(rbind,fisher_gene)
  rownames(fisher_pval) = unique(for_stack$Gene)
  
  write.csv(fisher_pval, file = paste0("figure1/",
                                      names(seurat_list)[i],
                                      "_Vang_Mcm5_fisher.csv"))
 
  p1 = ggplot(for_stack, 
              aes(fill=Detected, y=nCells, x=Phase)) + 
    geom_bar(position="stack", 
               stat="identity") +
    scale_fill_manual(values=c("#1565C0", "#90A4AE")) +
    geom_text(aes(label = nCells), 
              size = 3, 
              position = position_stack(vjust = 0.5)) +
    theme_classic() +
    ylab("number of cells") +
    facet_wrap(~Gene, 
               ncol = 2)
    
  # dot plots 
  p2 = ggplot(for_plot, 
              aes(x = Phase, y = Expression)) + 
    geom_jitter(shape=16, 
                colour = "#90A4AE",
                position=position_jitter(0.2)) +
    geom_boxplot(color = c("#C62828"),
                  fill = "NA",
                  outlier.size=0, 
                  outlier.shape=NA) +
    theme_classic() +
    ylab("Normalized Expression") +
    facet_wrap(~Gene, 
               ncol = 2)
  
  ggsave(p1/p2, 
         file = paste0("figure1/",
                       names(seurat_list)[i],
                       "_Vang_Mcm5.pdf"),
         width = 7,
         height = 5)
  }
  
  
  
  
  