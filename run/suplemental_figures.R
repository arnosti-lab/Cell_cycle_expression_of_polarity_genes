library(Seurat)
library(dplyr)
library(tidyr)
library(ggplot2)
library(reshape2)
library(FSA)
library(ggpubr)
library(rstatix)

setwd("./results")
load("drosophia_wing_disc_seurat_object.rdata")
seurat_wing = seurat_dat
load("drosophia_embryo_seurat_object.rdata")
seurat_embryo = seurat_dat
dir.create("supplemental")

seurat_list = list(seurat_embryo, seurat_wing)
names(seurat_list) = c("embryo", "wing_disc")

# Stacked bar plots + fisher -------------------------------------------------------

goi_df = read.csv("20200922_genes_of_interest_detected_in_dataset.csv",  row.names = 1)
for (i in 1:length(seurat_list)){
  goi_df = goi_df[which(goi_df[,i+2] == TRUE),]
  marker = goi_df[,2]
  names(marker) = goi_df$Gene
  goi = goi_df$Gene
  meta = seurat_list[[i]]@meta.data
  data = t(GetAssayData(seurat_list[[i]]))
  goi_exp= data[, as.character(goi)]
  phase_df = meta[c(6)]
  for_test = cbind(phase_df, goi_exp)
  for_test$Cell_ID = rownames(for_test)
  
  tidy = melt(for_test, id = c("Cell_ID", "Phase"))
  colnames(tidy) = c("Cell_ID", "Phase", "Gene", "Expression")
  
  tidy$Detected = ifelse(tidy$Expression > 0,
                          "above detection level",
                          "below detection level")
  
  tidy$Detected = factor(tidy$Detected, levels=c(unique(tidy$Detected)))
                                               
  for_stack = 
    tidy %>% 
    group_by(Phase,Gene,Detected, .drop=FALSE) %>% 
    summarise(nCells = length(Phase))
  
  for_stack$Phase = factor(for_stack$Phase, levels=c("G1", "S", "G2M"))
  
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
    names(fisher_gene[[gene]]) = c("G1", "G2M", "S")
  }
  fisher_pval = do.call(rbind,fisher_gene)
  rownames(fisher_pval) = unique(for_stack$Gene)
  
  write.csv(fisher_pval, file = paste0("supplemental/",
                                       names(seurat_list)[i],
                                       "_fisher_pval.csv"))
  
  adj_cols = apply(fisher_pval, 2, p.adjust, method = "BH")
  adj_rows = as.data.frame(t(apply(adj_cols, 1, p.adjust, method = "BH")))
  
  write.csv(adj_rows, file = paste0("supplemental/",
                                       names(seurat_list)[i],
                                       "_fisher_FDR.csv"))
  
  adj_rows$Gene = rownames(adj_rows)
  
  adjp = pivot_longer(adj_rows, c("G1", "G2M", "S"), names_to = c("Phase"))
  colnames(adjp)[3] = "FDR"
  adjp$FDR = signif(adjp$FDR, 2)
  adjp$FDRtext = ifelse(adjp$FDR < .001, "***",
                             ifelse(adjp$FDR < .01, "**",
                                    ifelse(adjp$FDR < .05, "*", "ns")))
  
  
  below = for_stack[which(for_stack$Detected == "below detection level"),]
  below$FDRtext = ""
  below$FDRpos = 0
  
  above = for_stack[which(for_stack$Detected == "above detection level"),]
  above$gene_phase = paste0(above$Gene, "_", above$Phase)
  above = above[order(above$gene_phase),]
  
  FDRpos = 
    tidy %>%
    group_by(Gene, Phase) %>%
    tally
  
  FDRpos$gene_phase = paste0(FDRpos$Gene, "_", FDRpos$Phase)
  FDRpos = FDRpos[order(FDRpos$gene_phase),]
  
  FDRpos$FDRpos = FDRpos$n + (.05 * max(FDRpos$n))
  above$FDRpos = FDRpos$FDRpos
  
  adjp$gene_phase = paste0(adjp$Gene, "_", adjp$Phase)
  adjp = adjp[order(adjp$gene_phase),]
  above$FDRtext = adjp$FDRtext
  
  for_stack = rbind(above, below)
  
  cycle = names(marker)[which(!marker == "polarity")]
  plot_cycle = for_stack[for_stack$Gene %in% cycle,]
  
  p2 = ggplot(plot_cycle, 
              aes(x=Phase, y=nCells, fill=Detected)) + 
    geom_bar(position="stack", 
             stat="identity") +
    scale_fill_manual(values=c("#1565C0", "#90A4AE")) +
    geom_text(aes(label = nCells),
              size = 2, 
              position = position_stack(vjust = 0.5)) +
    geom_text(aes(label = FDRtext,
                  x = Phase,
                  y = FDRpos), 
              size = 2) +
    theme_classic() +
    ylab("number of cells") +
    facet_wrap(~Gene, 
               ncol = 4)
    ggsave(p2, 
          file = paste0("supplemental/",
                        names(seurat_list)[i],
                        "_stacked_bar_cell_cycle_genes.pdf"),
          width = 10,
          height = 32.4)
    
    polarity = names(marker)[which(marker == "polarity")]
    plot_polarity = for_stack[for_stack$Gene %in% polarity,]
    
    p1 = ggplot(plot_polarity, 
                aes(x=Phase, y=nCells, fill=Detected)) + 
      geom_bar(position="stack", 
               stat="identity") +
      scale_fill_manual(values=c("#1565C0", "#90A4AE")) +
      geom_text(aes(label = nCells),
                size = 2, 
                position = position_stack(vjust = 0.5)) +
      geom_text(aes(label = FDRtext,
                    x = Phase,
                    y = FDRpos), 
                size = 2) +
      theme_classic() +
      ylab("number of cells") +
      facet_wrap(~Gene, 
                 ncol = 4)
    ggsave(p1, 
           file = paste0("supplemental/",
                         names(seurat_list)[i],
                         "_stacked_bar_polarity_genes.pdf"),
           width = 10,
           height = 9)
} 

# box/dot plots + kw dunns ------------------------------------------------
goi_df = read.csv("20200922_genes_of_interest_detected_in_dataset.csv",  row.names = 1)

for(i in 1:length(seurat_list)) {
  
  goi_df = goi_df[which(goi_df[,i+2] == TRUE),]
  marker = goi_df[,2]
  names(marker) = goi_df$Gene
  goi = goi_df$Gene
  meta = seurat_list[[i]]@meta.data
  exp = t(GetAssayData(seurat_list[[i]]))
  goi_exp= as.data.frame(exp[, goi])
  goi_exp$Phase = meta$Phase 
  goi_exp$Cell_ID = rownames(goi_exp)
  
  for_plot = melt(goi_exp, id = c("Cell_ID", "Phase"))
  colnames(for_plot) = c("Cell_ID", "Phase", "Gene", "Expression")
  for_plot$Phase <- factor(for_plot$Phase,levels = c("G1", "S", "G2M"))
  
  kw = for_plot %>% 
    group_by(Gene) %>%
    summarise(kw_pval = kruskal.test(Expression ~ Phase)$p.value)
  
  kw$kw_adj = p.adjust(kw$kw_pval, method = "BH") 
  
  dunn = for_plot %>% 
    group_by(Gene) %>% 
    summarise(dunnTest(Expression ~ Phase, method = "bh")[[2]]$P.adj)
  
  colnames(dunn) = c("Gene", "FDR")
  
  dunn$Phase_compare = rep(c("G1-G2M Dunns FDR",
                             "G1-S Dunns FDR", 
                             "G2M-S Dunns FDR"),
                           length(unique(dunn$Gene)))
  
  dunn = pivot_wider(dunn, names_from = Phase_compare, values_from = FDR)
  
  type_median = for_plot %>% 
    group_by(Gene,Phase) %>%
    summarise(median = median(Expression))
  
  type_median$namesfrom = rep(c("G1_median",
                                "G2M_median",
                                "S_median"),
                              length(unique(dunn$Gene)))
  
  type_median$Phase = NULL
  
  type_median = pivot_wider(type_median, names_from = namesfrom, values_from = median)
  
  type_mean = for_plot %>% 
    group_by(Gene,Phase) %>%
    summarise(mean = mean(Expression))
  
  type_mean$namesfrom = rep(c("G1_mean",
                              "G2M_mean",
                              "S_mean"),
                              length(unique(dunn$Gene)))
  
  type_mean$Phase = NULL
  
  type_mean = pivot_wider(type_mean, names_from = namesfrom, values_from = mean)
  
  mean_med = full_join(type_median, type_mean, by = "Gene")
  kw_dunn = full_join(mean_med, kw, by = "Gene")
  kw_dunn = full_join(kw_dunn, dunn, by = "Gene")
  
  write.csv(kw_dunn, file = paste0("supplemental/",
                                   names(seurat_list)[i],
                                   "_kw_dunns.csv"))
  
  kw_dunn$kw.signif = ifelse(kw_dunn$kw_adj < .001, "FDR < .001",
                                ifelse(kw_dunn$kw_adj < .01, "FDR < .01",
                                       ifelse(kw_dunn$kw_adj < .05, "FDR < .05", "ns")))
  
  kw_dunn = as.data.frame(kw_dunn)
  rownames(kw_dunn) = kw_dunn$Gene
  
  G1_S = kw_dunn[,c(1,11)]
  colnames(G1_S)[2] = "FDR"
  G1_S$group1 = rep("G1", nrow(G1_S))
  G1_S$group2 = rep("S", nrow(G1_S))
  
  G1_G2M = kw_dunn[,c(1,10)]
  colnames(G1_G2M)[2] = "FDR"
  G1_G2M$group1 = rep("G1", nrow(G1_G2M))
  G1_G2M$group2 = rep("G2M", nrow(G1_G2M))
  
  G2M_S = kw_dunn[,c(1,12)]
  colnames(G2M_S)[2] = "FDR"
  G2M_S$group1 = rep("S", nrow(G2M_S))
  G2M_S$group2 = rep("G2M", nrow(G2M_S))
  
  dunn_plot = rbind(G1_S, G1_G2M, G2M_S)
  dunn_plot$.y. = rep("Expression", nrow(dunn_plot))
  dunn_plot$FDR.signif = ifelse(dunn_plot$FDR < .001, "***",
                                ifelse(dunn_plot$FDR < .01, "**",
                                       ifelse(dunn_plot$FDR < .05, "*", "ns")))
  
  dunn_plot$term = "Phase"
  
  FDR_labeller = function(Gene) {paste(Gene, 
                                       "|", 
                                       kw_dunn[Gene, 13], 
                                       "|", 
                                       marker[Gene],
                                       "gene")}
  
  cycle = names(marker)[which(!marker == "polarity")]
  plot_cycle = for_plot[for_plot$Gene %in% cycle,]
  dunn_cycle = dunn_plot[dunn_plot$Gene %in% cycle,]
  
  dunn_cycle = dunn_cycle %>%
    add_y_position(formula = Expression ~ Phase, fun = "max", data = plot_cycle, step.increase = 3)
  
  P1 = ggplot(plot_cycle, aes(x = Phase, y = Expression)) + 
    geom_jitter(shape=16, 
                position=position_jitter(0.2)) +
    geom_boxplot(color = c("red"),
                 fill = "NA",
                 outlier.size=0, 
                 outlier.shape=NA) +
    theme_classic() +
    ylab("normalized expression") +
    facet_wrap(~Gene, 
               ncol = 5, 
               scales = "free_y",
               labeller = labeller(Gene = FDR_labeller)) +
    stat_pvalue_manual(dunn_cycle, label = "FDR.signif", tip.length = 0.01) 
    
  
  ggsave(file = paste0("supplemental/",
                       names(seurat_list)[i],
                       "_expression_vs_phase_dotplot_cell_cycle_genes.pdf"),
  width = 15,
  height = 40)
  
  polarity = names(marker)[which(marker == "polarity")]
  plot_polarity = for_plot[for_plot$Gene %in% polarity,]
  dunn_polarity = dunn_plot[dunn_plot$Gene %in% polarity,]
  
  dunn_polarity = dunn_polarity %>%
    add_y_position(formula = Expression ~ Phase, fun = "max", data = plot_polarity, step.increase = 5)
  
  P1 = ggplot(plot_polarity, aes(x = Phase, y = Expression)) + 
    geom_jitter(shape=16, 
                position=position_jitter(0.2)) +
    geom_boxplot(color = c("red"),
                 fill = "NA",
                 outlier.size=0, 
                 outlier.shape=NA) +
    theme_classic() +
    ylab("normalized expression") +
    facet_wrap(~Gene, 
               ncol = 4, 
               scales = "free_y",
               labeller = labeller(Gene = FDR_labeller)) +
    stat_pvalue_manual(dunn_polarity, 
                       label = "FDR.signif", 
                       tip.length = 0.01,
                       size = 2)
  
  
  ggsave(file = paste0("supplemental/",
                       names(seurat_list)[i],
                       "_expression_vs_phase_dotplot_polarity_genes.pdf"),
         width = 12,
         height = 12)
  
  
            }

