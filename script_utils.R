## Author : Ilango Guy
## Contact : guy.ilango@univ-tours.fr
## Doc : usefull function to source for R analysis


library(dplyr)
library(stringr)

prep_GO <- function(df){
original_gene_list <- df$log2FoldChange

# name the vector
names(original_gene_list) <- rownames(df)

# omit any NA values
gene_list<- original_gene_list

# sort the list in decreasing order (required for clusterProfiler)
gene_list = sort(gene_list, decreasing = TRUE)
gse <- gseGO(geneList=gene_list,
             ont ="ALL",
             keyType = "SYMBOL",
             minGSSize = 3,
             maxGSSize = 800,
             pvalueCutoff = 0.05,
             verbose = TRUE,
             OrgDb = organism,
             pAdjustMethod = "none")
return(gse)
}

baranal <- function(seurat_object){
  library(tidyverse)
  library(scCustomize)
  library(ggalluvial)
  cluster_stats <- as.data.frame(Cluster_Stats_All_Samples(seurat_object = seurat_object))
  cluster_stats <- cluster_stats %>% filter(row_number() <= n()-1)
  cluster_stats <- cluster_stats %>% select(Cluster , ends_with("%"))
  tab <- cluster_stats %>% gather(key = "keys" , value = "values" , -Cluster)
  tab %>% ggplot(aes(y = values, x = keys, fill = Cluster)) +
    geom_flow(aes(alluvium = Cluster), alpha= .5, color = "white",
              curve_type = "sigmoid",
              width = .5) +
    geom_col(width = .5, color = "white") +
    scale_y_continuous(NULL, expand = c(0,0)) +
    cowplot::theme_minimal_hgrid() +
    theme(panel.grid.major = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank())}

heatmap_bulk_zscore <- function(dds , title){
library(ClassDiscovery)
ntd <- normTransform(dds)
t <- assay(ntd)
zscores <- apply(t, 1, function(gene_expr) {
  (gene_expr - mean(gene_expr)) / sd(gene_expr)
})
zscores <- zscores[ , colSums(is.na(zscores)) == 0]
hc.features <- t(zscores) %>% t() %>%
  distanceMatrix(metric="pearson") %>%
  hclust(method="average")

hc.samples <- t(zscores) %>%
  distanceMatrix(metric="pearson") %>%
  hclust(method="average")

return(pheatmap::pheatmap(t(zscores),
                   scale = "row",

                   show_colnames  = T
                   ,main = title,
                   cluster_rows   = hc.features,
                   cluster_cols   = hc.samples) )
}


save_pheatmap_pdf <- function(x, filename, width=7, height=7) {
   stopifnot(!missing(x))
   stopifnot(!missing(filename))
   pdf(filename, width=width, height=height)
   grid::grid.newpage()
   grid::grid.draw(x$gtable)
   dev.off()
}





sc_annotate_de <- function(df) {
    df$gene <- rownames(df)
  df <- mutate(df, annot = case_when(
      str_detect(gene, '^Gm') & avg_log2FC > 0 & p_val_adj <= 0.05 ~ 'pas_interessant_up_sig',
      str_detect(gene, '^Gm') & avg_log2FC < 0 & p_val_adj <= 0.05 ~ 'pas_interessant_down_sig',
      str_detect(gene, '^Gm') &   p_val_adj > 0.05 ~ 'pas_interessant_unsig',
      str_detect(gene, '^mt') & avg_log2FC > 0 & p_val_adj <= 0.05 ~ 'pas_interessant_up_sig',
      str_detect(gene, '^mt') & avg_log2FC < 0 & p_val_adj <= 0.05 ~ 'pas_interessant_down_sig',
      str_detect(gene, '^mt') &   p_val_adj > 0.05 ~ 'pas_interessant_unsig',
      str_detect(gene, '^Rp[l|s]') & avg_log2FC > 0 & p_val_adj <= 0.05 ~ 'pas_interessant_up_sig',
      str_detect(gene, '^Rp[l|s]') & avg_log2FC < 0 & p_val_adj <= 0.05 ~ 'pas_interessant_down_sig',
      str_detect(gene, '^Rp[l|s]') &  p_val_adj > 0.05 ~ 'pas_interessant_unsig',
      str_detect(gene, '^Mir') & avg_log2FC > 0 & p_val_adj <= 0.05 ~ 'pas_interessant_up_sig',
      str_detect(gene, '^Mir') & avg_log2FC < 0 & p_val_adj <= 0.05 ~ 'pas_interessant_down_sig',
      str_detect(gene, '^Mir') &   p_val_adj > 0.05 ~ 'pas_interessant_unsig',
     avg_log2FC > 0 & p_val_adj <= 0.05 ~ 'up_sig',
      avg_log2FC < 0 & p_val_adj <= 0.05 ~ "down_sig",
      p_val_adj > 0.05 ~ "unsig",




))
}

sc_volcano_plotly <- function(df,title) {
  plot_ly(data = df, x = df$avg_log2FC, y = -log10(df$p_val_adj), text = rownames(df), mode = "markers", color = df$annot) %>%   layout(title = title)


}
sc_volcano <- function(df){
ggplot(df , aes( x = df$avg_log2FC , y = -log10(df$p_val_adj) , color = df$annot)) + geom_point()}



bulk_volcano_plotly <- function(df,title) {
  plot_ly(data = df, x = df$log2FoldChange, y = -log10(df$pvalue), text = rownames(df), mode = "markers", color = df$annot) %>%   layout(title = title)


}
bulk_volcano <- function(df){
ggplot(df , aes( x = df$log2FoldChange , y = -log10(df$pvalue) , color = df$annot)) + geom_point()}


make_it_bulk <- function(seurat_obj , sampleID){
    seurat_obj@meta.data$bulk <- seurat_obj@meta.data$sampleID
    Idents(seurat_obj) <- seurat_obj@meta.data$bulk
    markers <- FindAllMarkers(seurat_obj)
    return(markers)
    }

make_readable <- function(df){
df <- dplyr::filter(df , !grepl('Rik', rownames(df)))
df <- dplyr::filter(df , !grepl('^mt', rownames(df)))
df <- dplyr::filter(df , !grepl('^Rp[l|s]', rownames(df)))
df <- dplyr::filter(df , !grepl('^Mir', rownames(df)))
    df <- dplyr::filter(df , !grepl('^Gm[0-9]+', rownames(df)))
}

bulk_annotate_de <- function(df) {
  df <- mutate(df, annot = case_when(
    log2FoldChange > 0 ~ 'up',
    log2FoldChange < 0 ~ "down"

))
}


basic_pipe <- function(seurat_obj , dim , res ){
    seurat_obj <- NormalizeData(seurat_obj)
seurat_obj <- FindVariableFeatures(seurat_obj)
all.genes <- rownames(seurat_obj)
seurat_obj <- ScaleData(seurat_obj, features = all.genes)
seurat_obj <- RunPCA(seurat_obj, features = VariableFeatures(object = seurat_obj))

seurat_obj <- RunUMAP(seurat_obj, dims = 1:dim)
seurat_obj <- FindNeighbors(seurat_obj, dims = 1:dim)
seurat_obj <- FindClusters(seurat_obj , resolution = res)
return(seurat_obj)
}


doublet_finder <- function(seurat_obj){
    seurat_obj <- NormalizeData(seurat_obj)
    seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 2000)
    all.genes <- rownames(seurat_obj)
seurat_obj <- ScaleData(seurat_obj, features = all.genes)
    seurat_obj <- RunPCA(seurat_obj, features = VariableFeatures(object = seurat_obj))

# Find significant PCs
  stdv <- seurat_obj[["pca"]]@stdev
  sum.stdv <- sum(seurat_obj[["pca"]]@stdev)
  percent.stdv <- (stdv / sum.stdv) * 100
  cumulative <- cumsum(percent.stdv)
  co1 <- which(cumulative > 90 & percent.stdv < 5)[1]
  co2 <- sort(which((percent.stdv[1:length(percent.stdv) - 1] -
                       percent.stdv[2:length(percent.stdv)]) > 0.1),
              decreasing = T)[1] + 1
  min.pc <- min(co1, co2)
  min.pc

  # finish pre-processing
  seurat_obj <- RunUMAP(seurat_obj, dims = 1:min.pc)
  seurat_obj <- FindNeighbors(object = seurat_obj, dims = 1:min.pc)
  seurat_obj <- FindClusters(object = seurat_obj, resolution = 0.2)
sweep.list <- paramSweep(seurat_obj, PCs = 1:min.pc, num.cores = detectCores() - 1)
sweep.stats <- summarizeSweep(sweep.list)
bcmvn <- find.pK(sweep.stats)

# Optimal pK is the max of the bomodality coefficent (BCmvn) distribution
bcmvn.max <- bcmvn[which.max(bcmvn$BCmetric),]
optimal.pk <- bcmvn.max$pK
optimal.pk <- as.numeric(levels(optimal.pk))[optimal.pk]

## Homotypic doublet proportion estimate
annotations <- seurat_obj@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
nExp.poi <- round(optimal.pk * nrow(seurat_obj@meta.data)) ## Assuming 7.5% doublet formation rate - tailor for your dataset
nExp.poi.adj <- round(nExp.poi * (1 - homotypic.prop))

# run DoubletFinder
seurat_obj <- doubletFinder(seu = seurat_obj,
                                   PCs = 1:min.pc,
                                   pK = optimal.pk,
                                   nExp = nExp.poi.adj)
colnames(seurat_obj@meta.data)[8] <- "Doublet_Finder"
DimPlot(seurat_obj , group.by = "Doublet_Finder")
seurat_obj.singlets <- subset(seurat_obj, Doublet_Finder == "Singlet")
seurat_obj.singlets
seurat_obj.singlets
return(seurat_obj.singlets)
}

plot_cinetic <- function(table1, cluster, table2){
  p<-ggplot(table1 %>% filter(clusters == cluster) , aes(x = timepoint  , y = freq   , group = 1) ) +geom_line( linewidth = 1)

  p<- p + annotate("text", x = "J4", y = as.numeric(table1 %>% filter(clusters == cluster , timepoint == "J4") %>% select(freq)) -5, label = head(rownames(arrange(table2 %>% filter(condition == "NI_J4") , avg_log2FC) , 5))[5] , color ="blue")
  p<- p + annotate("text", x = "J4", y = as.numeric(table1 %>% filter(clusters == cluster , timepoint == "J4") %>% select(freq)) -4, label = head(rownames(arrange(table2 %>% filter(condition == "NI_J4") , avg_log2FC) , 5))[4] , color ="blue")
  p <- p + annotate("text", x = "J4", y = as.numeric(all_cell %>% filter(clusters == cluster , timepoint == "J4") %>% select(freq)) -3, label = head(rownames(arrange(table2 %>% filter(condition == "NI_J4") , avg_log2FC) , 5))[3] , color ="blue")
  p<-p + annotate("text", x = "J4", y = as.numeric(all_cell %>% filter(clusters == cluster , timepoint == "J4") %>% select(freq)) -2, label = head(rownames(arrange(table2 %>% filter(condition == "NI_J4") , avg_log2FC) , 5))[2] , color ="blue")
  p<- p + annotate("text", x = "J4", y = as.numeric(all_cell %>% filter(clusters == cluster , timepoint == "J4") %>% select(freq)) -1, label =head(rownames(arrange(table2 %>% filter(condition == "NI_J4") , avg_log2FC) , 5))[1] , color ="blue")
  p<- p + annotate("text", x = "J4", y = as.numeric(all_cell %>% filter(clusters == cluster , timepoint == "J4") %>% select(freq)) +5, label = tail(rownames(arrange(table2 %>% filter(condition == "NI_J4") , avg_log2FC) , 5))[1] , color ="red")
  p<- p + annotate("text", x = "J4", y = as.numeric(all_cell %>% filter(clusters ==cluster , timepoint == "J4") %>% select(freq)) +4, label = tail(rownames(arrange(table2 %>% filter(condition == "NI_J4") , avg_log2FC) , 5))[2] , color ="red")
  p<- p + annotate("text", x = "J4", y = as.numeric(all_cell %>% filter(clusters == cluster , timepoint == "J4") %>% select(freq)) +3, label = tail(rownames(arrange(table2 %>% filter(condition == "NI_J4") , avg_log2FC) , 5))[3] , color ="red")
  p<- p + annotate("text", x = "J4", y = as.numeric(all_cell %>% filter(clusters == cluster , timepoint == "J4") %>% select(freq)) +2, label = tail(rownames(arrange(table2 %>% filter(condition == "NI_J4") , avg_log2FC) , 5))[4] , color ="red")
  p<- p + annotate("text", x = "J4", y = as.numeric(all_cell %>% filter(clusters == cluster , timepoint == "J4") %>% select(freq)) +1, label = tail(rownames(arrange(table2 %>% filter(condition == "NI_J4") , avg_log2FC) , 5))[5] , color ="red")

  p<- p + annotate("text", x = "J7", y = as.numeric(table1 %>% filter(clusters == cluster , timepoint == "J7") %>% select(freq)) -5, label = head(rownames(arrange(table2 %>% filter(condition == "J4_J7") , avg_log2FC) , 5))[5] , color ="blue")
  p<- p + annotate("text", x = "J7", y = as.numeric(table1 %>% filter(clusters == cluster , timepoint == "J7") %>% select(freq)) -4, label = head(rownames(arrange(table2 %>% filter(condition == "J4_J7") , avg_log2FC) , 5))[4] , color ="blue")
  p <- p + annotate("text", x = "J7", y = as.numeric(all_cell %>% filter(clusters == cluster , timepoint == "J7") %>% select(freq)) -3, label = head(rownames(arrange(table2 %>% filter(condition == "J4_J7") , avg_log2FC) , 5))[3] , color ="blue")
  p<-p + annotate("text", x = "J7", y = as.numeric(all_cell %>% filter(clusters == cluster , timepoint == "J7") %>% select(freq)) -2, label = head(rownames(arrange(table2 %>% filter(condition == "J4_J7") , avg_log2FC) , 5))[2] , color ="blue")
  p<- p + annotate("text", x = "J7", y = as.numeric(all_cell %>% filter(clusters == cluster , timepoint == "J7") %>% select(freq)) -1, label =head(rownames(arrange(table2 %>% filter(condition == "J4_J7") , avg_log2FC) , 5))[1] , color ="blue")
  p<- p + annotate("text", x = "J7", y = as.numeric(all_cell %>% filter(clusters == cluster , timepoint == "J7") %>% select(freq)) +5, label = tail(rownames(arrange(table2 %>% filter(condition == "J4_J7") , avg_log2FC) , 5))[1] , color ="red")
  p<- p + annotate("text", x = "J7", y = as.numeric(all_cell %>% filter(clusters ==cluster , timepoint == "J7") %>% select(freq)) +4, label = tail(rownames(arrange(table2 %>% filter(condition == "J4_J7") , avg_log2FC) , 5))[2] , color ="red")
  p<- p + annotate("text", x = "J7", y = as.numeric(all_cell %>% filter(clusters == cluster , timepoint == "J7") %>% select(freq)) +3, label = tail(rownames(arrange(table2 %>% filter(condition == "J4_J7") , avg_log2FC) , 5))[3] , color ="red")
  p<- p + annotate("text", x = "J7", y = as.numeric(all_cell %>% filter(clusters == cluster , timepoint == "J7") %>% select(freq)) +2, label = tail(rownames(arrange(table2 %>% filter(condition == "J4_J7") , avg_log2FC) , 5))[4] , color ="red")
  p<- p + annotate("text", x = "J7", y = as.numeric(all_cell %>% filter(clusters == cluster , timepoint == "J7") %>% select(freq)) +1, label = tail(rownames(arrange(table2 %>% filter(condition == "J4_J7") , avg_log2FC) , 5))[5] , color ="red")

  p<- p + annotate("text", x = "J12", y = as.numeric(table1 %>% filter(clusters == cluster , timepoint == "J12") %>% select(freq)) -5, label = head(rownames(arrange(table2 %>% filter(condition == "J7_J12") , avg_log2FC) , 5))[5] , color ="blue")
  p<- p + annotate("text", x = "J12", y = as.numeric(table1 %>% filter(clusters == cluster , timepoint == "J12") %>% select(freq)) -4, label = head(rownames(arrange(table2 %>% filter(condition == "J7_J12") , avg_log2FC) , 5))[4] , color ="blue")
  p <- p + annotate("text", x = "J12", y = as.numeric(all_cell %>% filter(clusters == cluster , timepoint == "J12") %>% select(freq)) -3, label = head(rownames(arrange(table2 %>% filter(condition == "J7_J12") , avg_log2FC) , 5))[3] , color ="blue")
  p<-p + annotate("text", x = "J12", y = as.numeric(all_cell %>% filter(clusters == cluster , timepoint == "J12") %>% select(freq)) -2, label = head(rownames(arrange(table2 %>% filter(condition == "J7_J12") , avg_log2FC) , 5))[2] , color ="blue")
  p<- p + annotate("text", x = "J12", y = as.numeric(all_cell %>% filter(clusters == cluster , timepoint == "J12") %>% select(freq)) -1, label =head(rownames(arrange(table2 %>% filter(condition == "J7_J12") , avg_log2FC) , 5))[1] , color ="blue")
  p<- p + annotate("text", x = "J12", y = as.numeric(all_cell %>% filter(clusters == cluster , timepoint == "J12") %>% select(freq)) +5, label = tail(rownames(arrange(table2 %>% filter(condition == "J7_J12") , avg_log2FC) , 5))[1] , color ="red")
  p<- p + annotate("text", x = "J12", y = as.numeric(all_cell %>% filter(clusters ==cluster , timepoint == "J12") %>% select(freq)) +4, label = tail(rownames(arrange(table2 %>% filter(condition == "J7_J12") , avg_log2FC) , 5))[2] , color ="red")
  p<- p + annotate("text", x = "J12", y = as.numeric(all_cell %>% filter(clusters == cluster , timepoint == "J12") %>% select(freq)) +3, label = tail(rownames(arrange(table2 %>% filter(condition == "J7_J12") , avg_log2FC) , 5))[3] , color ="red")
  p<- p + annotate("text", x = "J12", y = as.numeric(all_cell %>% filter(clusters == cluster , timepoint == "J12") %>% select(freq)) +2, label = tail(rownames(arrange(table2 %>% filter(condition == "J7_J12") , avg_log2FC) , 5))[4] , color ="red")
  p<- p + annotate("text", x = "J12", y = as.numeric(all_cell %>% filter(clusters == cluster , timepoint == "J12") %>% select(freq)) +1, label = tail(rownames(arrange(table2 %>% filter(condition == "J7_J12") , avg_log2FC) , 5))[5] , color ="red")
  return(p)

}


dotplot_go_spe_mouse <- function(table , cluster){
cluster_down <- enrichGO(gene         = rownames(table %>% filter(avg_log2FC < 0)),
                OrgDb         = org.Mm.eg.db,
                keyType       = 'SYMBOL',
                ont           = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.05)

cluster_up <- enrichGO(gene         = rownames(table %>% filter( avg_log2FC > 0)),
                OrgDb         = org.Mm.eg.db,
                keyType       = 'SYMBOL',
                ont           = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.05)

a<-dotplot(cluster_down) + ggtitle(paste0("Cluster ", as.character(cluster) ," downregulated gene"))

b <- dotplot(cluster_up) + ggtitle(paste0("Cluster ",as.character(cluster) ," upregulated gene"))
return(list(a , b))

                                   }

dotplot_go_mouse <- function(table , cluster){
cluster_down <- enrichGO(gene         = rownames(table %>% filter(avg_log2FC < 0 , cluster == cluster) ),
                OrgDb         = org.Mm.eg.db,
                keyType       = 'SYMBOL',
                ont           = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.05)

cluster_up <- enrichGO(gene         = rownames(table %>% filter( avg_log2FC > 0, cluster == cluster)),
                OrgDb         = org.Mm.eg.db,
                keyType       = 'SYMBOL',
                ont           = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.05)

a<-dotplot(cluster_down) + ggtitle(paste0("Cluster ", as.character(cluster) ," downregulated gene"))

b <- dotplot(cluster_up) + ggtitle(paste0("Cluster ",as.character(cluster) ," upregulated gene"))
return(list(a , b))

                                   }



##### signatures

GOBP_REGULATION_OF_VASCULAR_WOUND_HEALING <- str_to_title(c("Tafa5","Alox5","Cxcr4","Foxc2","Hif1a","Serpine1","Slc12a2","Tnf","Vegfb","Xbp1","Smoc2"
))
GOBP_REGULATION_OF_WOUND_HEALING <- str_to_title(c("Tpsab1","Ano6","Tafa5","Kank1","Mylk","Adtrp","F11","Fgb","Actg1","Acta2","Adra2a","Adrb2","Alox12","Alox5","Anxa5","Apoe","Apoh","Serpinc1","Serping1","Anxa2","Cask","Cav1","Cd34","Cd36","Cd9","Cldn1","Cldn3","Cldn4","Cxcr4","Crk","Dmtn","Ephb2","F2","F2r","F2rl1","F3","F7","Ptk2","Fga","Fgf2","Foxc2","Gja1","Gp1ba","Gp5","Hbegf","Hif1a","Hmgb1","Hmgcr","Foxa2","Hpse","Hras","Insl3","Itgb1","Klkb1","Kng1","Sh2b3","Anxa1","Smad3","Nfe2l2","Ddr2","Pdgfa","Pdgfb","Pdgfra","Prkcd","Prkce","Serpine1","Plat","Plau","Plaur","Plg","Serpinf2","Prkg1","Proc","Pros1","Pten","Ptger3","Ptger4","Reg3a","Reg3g","S100a9","Ccl2","St3gal4","Slc12a2","Serpine2","Phldb2","Adamts18","Abcc8","Arfgef1","Tbxa2r","Duox2","Tspan8","Prdx2","Tfpi","Thbd","Thbs1","Fermt2","Tnf","Vegfb","Vil1","Eppk1","Vtn","Ccn4","Wnt4","Xbp1","Enpp4","Ajap1","Klrh1","Cd109","Fermt1","Cldn19","Emilin2","Cadm4","Ceacam1","Abat","Tnfrsf12a","Vkorc1","Kng2","Cpb2","Clec7a","Mtor","C1qtnf1","Psg23","Cldn13","F12","Myoz1","Smoc2","Wfdc1","Srsf6","Rreb1","Mmrn1","Ubash3b","Muc16","Clasp2","Clasp1","Hrg","Duox1","Fgg"

))

GOBP_VASCULAR_WOUND_HEALING <- str_to_title(c("Tafa5","Alox5","Cxcr4","Foxc2","Gata2","Hif1a","Hpse","Kdr","Serpine1","Slc12a2","Tnf","Vegfa","Vegfb","Xbp1","Npr2","Smoc2","Ndnf","Adipor2"


))
GOBP_REGULATION_OF_WOUND_HEALING <- str_to_title(c("Tpsab1","Ano6","Tafa5","Kank1","Mylk","Adtrp","F11","Fgb","Actg1","Acta2","Adra2a","Adrb2","Alox12","Alox5","Anxa5","Apoe","Apoh","Serpinc1","Serping1","Anxa2","Cask","Cav1","Cd34","Cd36","Cd9","Cldn1","Cldn3","Cldn4","Cxcr4","Crk","Dmtn","Ephb2","F2","F2r","F2rl1","F3","F7","Ptk2","Fga","Fgf2","Foxc2","Gja1","Gp1ba","Gp5","Hbegf","Hif1a","Hmgb1","Hmgcr","Foxa2","Hpse","Hras","Insl3","Itgb1","Klkb1","Kng1","Sh2b3","Anxa1","Smad3","Nfe2l2","Ddr2","Pdgfa","Pdgfb","Pdgfra","Prkcd","Prkce","Serpine1","Plat","Plau","Plaur","Plg","Serpinf2","Prkg1","Proc","Pros1","Pten","Ptger3","Ptger4","Reg3a","Reg3g","S100a9","Ccl2","St3gal4","Slc12a2","Serpine2","Phldb2","Adamts18","Abcc8","Arfgef1","Tbxa2r","Duox2","Tspan8","Prdx2","Tfpi","Thbd","Thbs1","Fermt2","Tnf","Vegfb","Vil1","Eppk1","Vtn","Ccn4","Wnt4","Xbp1","Enpp4","Ajap1","Klrh1","Cd109","Fermt1","Cldn19","Emilin2","Cadm4","Ceacam1","Abat","Tnfrsf12a","Vkorc1","Kng2","Cpb2","Clec7a","Mtor","C1qtnf1","Psg23","Cldn13","F12","Myoz1","Smoc2","Wfdc1","Srsf6","Rreb1","Mmrn1","Ubash3b","Muc16","Clasp2","Clasp1","Hrg","Duox1","Fgg"
))

GOBP_WOUND_HEALING <- str_to_title(c("Tpsab1","Fer1l5","Nlrp6","Chmp7","Ano6","Scrib","Tafa5","Mpig6b","Kank1","Mylk","Fermt3","Tspan9","Adtrp","Dsp","F11","Fgb","Fntb","Macf1","Actg1","Acta2","Acvrl1","Adra2a","Adra2b","Adra2c","Adrb1","Adrb2","Vps4a","Alox12","Alox15","Alox5","Bloc1s4","Anxa5","Anxa6","Anxa8","Ap3b1","Apoe","Apoh","Aqp1","Rab27a","Serpinc1","Bnc1","Serping1","C3","C9","Pdia4","Ddr1","Anxa2","Cask","Casp7","Cav1","Cav3","Cd151","Cd34","Cd36","Entpd1","Entpd2","Cd44","Cd9","Cdkn1a","Celsr1","Cfh","Cflar","Cldn1","Cldn3","Cldn4","Cxcr4","Ccr2","Cnn2","Col3a1","Col5a1","Col1a1","Comp","Crk","Crp","Ctsg","Cxadr","Drd5","Dst","Egfr","Elk3","Eng","Dmtn","Ephb2","Erbb2","Evl","Evpl","Ext1","F10","F13b","F2","F2r","F2rl1","F2rl2","F2rl3","F3","F5","F7","F8","F9","Ptk2","Fbln1","Fcer1g","Fga","Fgf1","Fgf10","Fgf2","Fgf7","Fkbp10","Foxc2","Fn1","Fzd6","Fzd7","Gas6","Gata1","Gata2","Gata4","B4galt1","Gja1","Gna13","Gnaq","Gnas","Gp1ba","Gp1bb","Pdpn","Lilrb4a","Gp5","Gpx1","Pdia3","Serpind1","Ptpn6","Hbegf","Hif1a","Hmgb1","Hmgcr","Hmox1","Foxa2","Hnf4a","Hpse","Hras","Igf1","Ccn1","Il1a","Il6","Il6ra","Ins1","Ins2","Insl3","Itga2b","Itgb1","Itgb3","Itgb6","Jak2","F11r","Ajuba","Kdr","Klkb1","Kng1","Krt6a","Sh2b3","Lox","Anxa1","Lyn","Lyst","Smad3","Smad4","Mertk","Clec10a","Mmp12","Msx2","Myh2","Myh9","Naca","Nf1","Nfe2l2","Nog","Notch2","Slc11a1","Ddr2","P2rx1","P2ry1","Bloc1s6","Pak1","Pdgfa","Pdgfb","Pdgfra","Pecam1","Pip5k1c","Prkca","Prkcd","Prkce","Prkcq","Pla2g4a","Serpine1","Serpinb2","Plat","Plau","Plaur","Plec","Plg","Serpinf2","Pou2f3","Ppara","Ppard","Pparg","Ppl","Prkg1","Proc","Procr","Pros1","Pten","Flna","Ptger3","Ptger4","Hps4","Hps1","Ptprj","Rab3a","Reg3a","Reg3g","Hps6","S100a10","S100a9","Scnn1b","Scnn1g","Ccl2","Cx3cl1","Selp","Shh","St3gal4","Vps4b","Slc12a2","Slc4a1","Snai2","Smpd1","Serpine2","Srf","Chmp6","Phldb2","Adamts18","Stxbp1","Stxbp3","Abcc8","Syk","Sdc1","Sdc4","Nrg1","Arfgef1","Arhgef19","Tbxa2r","Duox2","Tspan8","Prdx2","Tec","Serpina10","Tfpi","Tfpi2","Tgfb1","Tgfb2","Thbd","Thbs1","Timp1","Fermt2","Tlr4","Tnf","Cd40lg","Tpm1","Txk","Tyro3","Vegfa","Vegfb","Vil1","Eppk1","Vtn","Mrtfa","Vwf","Ccn4","Wnt3a","Wnt4","Wnt5a","Wnt7a","Xbp1","Enpp4","Yap1","Myof","Ccm2l","Syt11","Npr2","Grhl3","Ajap1","Arhgap24","Klrh1","Arhgap35","Bloc1s3","Vps33b","Fgl1","Carmil2","Chmp1a","Cd109","Nbeal2","Coro1b","Papss2","Fermt1","Cldn19","Stard13","Gp6","Tsku","Hps5","Emilin2","Cadm4","Axl","Ceacam1","Map3k1","Slc7a11","Ppia","Abat","Dysf","Jaml","Tspan32","Tnfrsf12a","Sytl4","Adamts13","Vkorc1","Tor1a","Gpr4","Abi3bp","Ubash3a","Mia3","Kng2","Mfsd2b","Trim72","Epb41l4b","Gp9","Hgfac","Syt7","Tubb1","Plek","Tmeff2","Cpb2","Pdcd10","Clec7a","Mtor","Pf4","C1qtnf1","Psg23","Cldn13","F12","Myoz1","C1galt1c1","Smoc2","Mustn1","Chmp4c","Ahnak","Proz","Chmp1b","Arl8b","Fundc2","Fgfr1op2","Wfdc1","Plpp3","Srsf6","Ndnf","Adipor2","Rreb1","Chmp2b","Chmp2a","Pdia2","Prss56","Lnpk","P2ry12","Mmrn1","Treml1","Ptk7","Prcp","Ubash3b","Pear1","Dcbld2","Muc16","Rap2b","F13a1","Pik3cb","Chmp4b","Clasp2","Plet1","Clasp1","Lrg1","Chmp5","Myh10","BC004004","Pard3","Vangl2","Hrg","Dtnbp1","Duox1","Fgg"


))

GOBP_POSITIVE_REGULATION_OF_RESPONSE_TO_WOUNDING <- str_to_title(c("Tpsab1","Fer1l5","Nlrp6","Chmp7","Ano6","Scrib","Tafa5","Mpig6b","Kank1","Mylk","Fermt3","Tspan9","Adtrp","Dsp","F11","Fgb","Fntb","Macf1","Actg1","Acta2","Acvrl1","Adra2a","Adra2b","Adra2c","Adrb1","Adrb2","Vps4a","Alox12","Alox15","Alox5","Bloc1s4","Anxa5","Anxa6","Anxa8","Ap3b1","Apoe","Apoh","Aqp1","Rab27a","Serpinc1","Bnc1","Serping1","C3","C9","Pdia4","Ddr1","Anxa2","Cask","Casp7","Cav1","Cav3","Cd151","Cd34","Cd36","Entpd1","Entpd2","Cd44","Cd9","Cdkn1a","Celsr1","Cfh","Cflar","Cldn1","Cldn3","Cldn4","Cxcr4","Ccr2","Cnn2","Col3a1","Col5a1","Col1a1","Comp","Crk","Crp","Ctsg","Cxadr","Drd5","Dst","Egfr","Elk3","Eng","Dmtn","Ephb2","Erbb2","Evl","Evpl","Ext1","F10","F13b","F2","F2r","F2rl1","F2rl2","F2rl3","F3","F5","F7","F8","F9","Ptk2","Fbln1","Fcer1g","Fga","Fgf1","Fgf10","Fgf2","Fgf7","Fkbp10","Foxc2","Fn1","Fzd6","Fzd7","Gas6","Gata1","Gata2","Gata4","B4galt1","Gja1","Gna13","Gnaq","Gnas","Gp1ba","Gp1bb","Pdpn","Lilrb4a","Gp5","Gpx1","Pdia3","Serpind1","Ptpn6","Hbegf","Hif1a","Hmgb1","Hmgcr","Hmox1","Foxa2","Hnf4a","Hpse","Hras","Igf1","Ccn1","Il1a","Il6","Il6ra","Ins1","Ins2","Insl3","Itga2b","Itgb1","Itgb3","Itgb6","Jak2","F11r","Ajuba","Kdr","Klkb1","Kng1","Krt6a","Sh2b3","Lox","Anxa1","Lyn","Lyst","Smad3","Smad4","Mertk","Clec10a","Mmp12","Msx2","Myh2","Myh9","Naca","Nf1","Nfe2l2","Nog","Notch2","Slc11a1","Ddr2","P2rx1","P2ry1","Bloc1s6","Pak1","Pdgfa","Pdgfb","Pdgfra","Pecam1","Pip5k1c","Prkca","Prkcd","Prkce","Prkcq","Pla2g4a","Serpine1","Serpinb2","Plat","Plau","Plaur","Plec","Plg","Serpinf2","Pou2f3","Ppara","Ppard","Pparg","Ppl","Prkg1","Proc","Procr","Pros1","Pten","Flna","Ptger3","Ptger4","Hps4","Hps1","Ptprj","Rab3a","Reg3a","Reg3g","Hps6","S100a10","S100a9","Scnn1b","Scnn1g","Ccl2","Cx3cl1","Selp","Shh","St3gal4","Vps4b","Slc12a2","Slc4a1","Snai2","Smpd1","Serpine2","Srf","Chmp6","Phldb2","Adamts18","Stxbp1","Stxbp3","Abcc8","Syk","Sdc1","Sdc4","Nrg1","Arfgef1","Arhgef19","Tbxa2r","Duox2","Tspan8","Prdx2","Tec","Serpina10","Tfpi","Tfpi2","Tgfb1","Tgfb2","Thbd","Thbs1","Timp1","Fermt2","Tlr4","Tnf","Cd40lg","Tpm1","Txk","Tyro3","Vegfa","Vegfb","Vil1","Eppk1","Vtn","Mrtfa","Vwf","Ccn4","Wnt3a","Wnt4","Wnt5a","Wnt7a","Xbp1","Enpp4","Yap1","Myof","Ccm2l","Syt11","Npr2","Grhl3","Ajap1","Arhgap24","Klrh1","Arhgap35","Bloc1s3","Vps33b","Fgl1","Carmil2","Chmp1a","Cd109","Nbeal2","Coro1b","Papss2","Fermt1","Cldn19","Stard13","Gp6","Tsku","Hps5","Emilin2","Cadm4","Axl","Ceacam1","Map3k1","Slc7a11","Ppia","Abat","Dysf","Jaml","Tspan32","Tnfrsf12a","Sytl4","Adamts13","Vkorc1","Tor1a","Gpr4","Abi3bp","Ubash3a","Mia3","Kng2","Mfsd2b","Trim72","Epb41l4b","Gp9","Hgfac","Syt7","Tubb1","Plek","Tmeff2","Cpb2","Pdcd10","Clec7a","Mtor","Pf4","C1qtnf1","Psg23","Cldn13","F12","Myoz1","C1galt1c1","Smoc2","Mustn1","Chmp4c","Ahnak","Proz","Chmp1b","Arl8b","Fundc2","Fgfr1op2","Wfdc1","Plpp3","Srsf6","Ndnf","Adipor2","Rreb1","Chmp2b","Chmp2a","Pdia2","Prss56","Lnpk","P2ry12","Mmrn1","Treml1","Ptk7","Prcp","Ubash3b","Pear1","Dcbld2","Muc16","Rap2b","F13a1","Pik3cb","Chmp4b","Clasp2","Plet1","Clasp1","Lrg1","Chmp5","Myh10","BC004004","Pard3","Vangl2","Hrg","Dtnbp1","Duox1","Fgg"


))

GOBP_NEGATIVE_REGULATION_OF_TISSUE_REMODELING <- str_to_title(c("Gpr137","Abr","Bcr","Agt","Cd24a","Cd38","Csk","Cst3","Cyp19a1","Fshr","Gata4","Iapp","Il6","Inpp5d","Tnfrsf11b","P2rx7","Pparg","Sfrp1","Vegfa","Tmem119","Grem1","Ypel4","Ceacam1","Cartpt","Abi3bp","Cldn18","Ubash3b","Gpr137b","Hamp"


))

GOBP_ACUTE_INFLAMMATORY_RESPONSE_TO_ANTIGENIC_STIMULUS <- c("Nlrp6","Ano6","Npy","Selenos","H2-Q6","H2-M10.3","Acvr1","Adam8","Adcyap1","Adora1","Ahsg","Alox5ap","Aoc3","Bdkrb2","Btk","C3","Casp6","Cd6","Cxcr2","Ccr5","Ccr7","Cnr1","Crp","Dnase1","Dnase1l3","Ednrb","Ephb6","Ext1","F2","F3","F8","Fcer1a","Fcer1g","Fcgr1","Fcgr2b","Fcgr3","Fn1","Fut7","B4galt1","Gstp1","H2-Bl","H2-D1","H2-K1","","H2-M10.1","H2-M2","H2-M9","H2-Q1","H2-Q10","H2-Q2","H2-Q4","H2-Q7","H2-T10","H2-T22","H2-T23","H2-T3","Hp","Eif2ak1","Icam1","Ighg2b","Ighg1","","","Il1a","Il1b","Il1rn","Il4","Il6","Ins1","Ins2","Itih4","Klk1b1","Lbp","Lipa","Npy5r","Orm1","Orm2","Reg3b","Pla2g2d","Serpinf2","Pparg","Ptger3","Ash1l","Ptgs2","Trpv1","Reg3a","Reg3g","Saa1","Saa2","Saa3","Saa4","Ccl5","Serpina1a","Serpina1b","Serpina3n","Serpinb9","Spn","Stat3","Stat5b","Il20rb","Tac1","Mylk3","Nlrp3","Il31ra","Tnf","Tnfrsf11a","Tnfsf11","Plscr1","Tnfsf4","Vnn1","Vwf","H2-M10.4","H2-M11","H2-M1","H2-M10.5","Zp3","Ffar2","Ffar3","Mrgpra3","H2-M5","Scn11a","Sigirr","C2cd4a","Rhbdd3","Pik3cg","H2-M10.2","Ighe","","H2-M10.6","Elane","Nupr1","Park7","F12","Ptges","H2-T-ps","Ctnnbip1","Prcp","C2cd4b","Cd163")
GOBP_ACTIVATION_OF_IMMUNE_RESPONSE <- c("Ifi208","Btnl1","Mndal","Lrrc19","Smpdl3b","Epg5","Zfp683","Oas1d","Slc15a4","Plekha1","Prkd2","Nlrp6","Mapkapk3","Ifi206","Btnl12","Pvrig","Cd276","Dgkz","Dusp22","Nckap1l","Ticam1","Slc39a6","Lpxn","Lgr4","Nod1","Rapgef1","Eif2b3","Foxp1","Card11","Irak2","C7","Braf","Rap1a","C8b","Rps6ka3","D1Pas1","Abl1","Oas1c","Pawr","Ada","Cfd","Tirap","Ap3b1","Cd5l","App","Arf6","Bax","Phb2","Bcl10","Bcl2","Bcl2a1d","Blk","Bmx","Btk","Btn1a1","Serping1","C1qa","C1qb","C1qbp","C1qc","C2","C3","C3ar1","C4b","C4bp","C5ar1","C6","C9","Cacnb3","Cacnb4","Casp1","Casp4","Casp6","Cav1","Cd14","Ctla4","Cd19","Ms4a1","Cd22","Cd28","Cd36","Cd38","Cd3e","Cd247","Cd59a","Cd79a","Cd81","Cd86","Cd8a","Cfh","Cfi","Chuk","Ccr7","Cr2","Bcar1","Crkl","Crp","Cr1l","Csk","Cyba","Cd55","Cd55b","Ddx3x","Eif2b4","Elf1","Esr1","F2rl1","Colec12","Fcer1g","Fcna","Fcnb","Flot1","Fosl1","Fosl2","Fpr2","Fpr-rs3","Fpr-rs4","Fpr1","Fpr3","Tlr3","Fyn","Xrcc6","Gata3","Gbp2b","Gbp2","Usp15","Gcsam","Gfi1","Lilrb4b","Lilrb4a","Cmklr1","Gpld1","Gpr33","Grb2","Cfb","Hc","Ptpn6","Hmgb1","Hras","Hspa8","Hspd1","Hsp90aa1","Irgm1","Ifi203","Ifi204","Ifng","Cd79b","Ighg2b","Ighg1","Ighm","","Igtp","Ikbkg","Il1b","Irak1","Irf1","Irf2","Irf4","Acod1","Cd47","Itk","Kcnj8","Kcnn4","Klrc1","Klrc2","Klrd1","Krt1","Lamp2","Laptm5","Lat","Lbp","Lck","Lcp2","Lgals3","Lipa","Ltf","Blnk","Znrf1","Tlr7","Tlr8","Ly96","Lyn","Havcr2","Mapkapk2","Masp1","Masp2","Matr3","Mbl1","Mbl2","Cd46","Mef2c","Mog","Myd88","Nfatc2","Nfkb1","Nfkbia","Nfkbil1","Ninj1","Nr4a3","Nras","Pde4b","Pdpk1","Cfp","Phb","Pik3r1","Prkcb","Prkce","Prkch","Prkd1","Plcg1","Plscr2","Ppp2ca","Prkdc","Eif2ak2","Prnp","Psen1","Psen2","Pten","Btnl10","Hexim1","Ptpn2","Ptpn22","Ptprc","Ripk2","Ptprj","Ptprs","Nectin2","Nlrp1a","Skint3","Reg3g","Rela","Trim30a","Nr1h4","Khdrbs1","Foxp3","Sh2d1a","Sin3a","Clpb","Sos1","Src","Znrf4","Cblb","Stk11","Eif2b1","Trim30d","Syk","Themis","Brcc3","Lacc1","Tifa","Tarbp2","Cfhr4","Cgas","","Appl2","Nlrp3","Tec","Arrb2","Dhx33","Trim25","Nr1d1","Cd300a","Nploc4","Eif2b2","Thy1","Tlr1","Tlr4","Tlr6","Tnf","Tnfaip3","Cd40","Traf3","Traf6","Plscr1","Trex1","Txk","Tyro3","Tyrobp","Usp12","Ufd1","Nr1h3","Usp9x","Ezr","Lrrc14","Gramd4","Eif2b5","Bag6","Trim31","Treml4","Plcl2","Pja2","Ticam2","Cd226","Tkfc","Xrcc5","Zap70","Rab7b","Rab29","Ifi207","Ifi205","Slc39a10","Zp3r","Zdhhc5","Mavs","Bpifb1","Gbp5","Gbp7","Rigi","Shb","C8a","Skint10","Skint6","Skint11","Zc3h12a","Themis2","Tnip2","Oasl1","Oas1e","Lrch4","Wnk1","A2m","C1rl","Mark4","Ffar2","Ddx60","Carmil2","Plcg2","Ifi209","Btnl9","Sarm1","Igha","Btn2a2","Fyb","Pde4d","Tlr11","Sh2b2","Colec10","Oas1g","Oas1b","Klhl6","Malt1","Peli3","Lax1","Tlr2","Rbck1","Vtcn1","Fyb2","Skint5","Gpatch3","Oas1f","Klre1","Nfkbid","Trim30b","Nlrp10","Myo1g","Oas3","Oas1h","Oas1a","Cd300lf","Nod2","Map3k7","Mapk1","Mapk8","Rnf31","Nlrc3","Nlrc4","Ecsit","Klrk1","Nlrx1","Ermap","Rps3","Naglu","Vsig4","Tlr13","Tomm70a","Slc46a2","C1s2","Trim12c","C5ar2","Rc3h2","Ptgs2os","Klri2","Skint4","Fpr-rs6","Fpr-rs7","Scimp","Skint7","Ubash3a","Dennd1b","Fcrl5","Skint9","Skint2","Cd59b","Brcc3dc","Pram1","Ighe","","","Ighg3","Rc3h1","Ifi211","Cd300ld3","Aim2","Tlr12","Ighg2c","Tnip3","Otulin","Trim30c","Klri1","Cfhr1","Icosl","C1s1","C1ra","Mfhas1","Nono","Tlr5","Irf7","Irf3","Cd160","Irgm2","Unc93b1","Mefv","Cfhr2","Ifi214","Pqbp1","Btnl2","Gbp3","Trim3","Ubqln1","Nagk","Gps2","Tbk1","Tspan6","Clec4e","Clec4n","Clec7a","Rabgef1","Lat2","Stap1","Vav3","Slc15a2","Tnip1","Klrc3","Rsad2","Zbp1","Erbin","Nek7","Ifi213","Btnl6","C4a","Gm12250","Btnl4","Trav7-2","Nlrp1b","Skint8","Skint1","Sirt2","Nmi","Svep1","Slc15a3","S100a14","Rgcc","Zdhhc12","Tmem126a","Stoml2","Pspc1","Nepn","C1rb","Ncr3-ps","Trim5","Pycard","Tril","Riok3","Lsm14a","Peli1","Inava","Tespa1","Rnf125","Stmp1","Zcchc3","Hcfc2","Nop53","Cmtm3","Rtn4","Trim15","Wdfy1","C8g","Dab2ip","Usp46","Ifi35","Cactin","Ipo5","Lrrfip2","Tasl","Alpk1","Sfpq","Ifih1","Colec11","Rnf135","Fbxl2","Dusp3","Sting1","Lime1","Appl1","Atat1","Irak3","Otud4","Fcho1","Nfam1","Sec14l1","Cyld","Themis3","Sppl3","Rab11fip2","Usp50","Phpt1","Rftn1","Trim12a","Trat1","Sla2","Gpr108","Skap1","Zc3hav1","Cptp","Nfkbiz","Dhx58","Pum1","Pum2","Ankrd17","Tlr9","Trem2","Pik3ap1","Clec2i","Csnk1a1","Ube2n","Trim11","Tnfrsf21","Susd4","Znfx1")
GOBP_NEUTROPHIL_ACTIVATION <- c("Traf3ip2","Abr","Bcr","Anxa3","Cxcr2","Camp","Ctsg","Dnase1","Dnase1l3","F2rl1","Fcer1g","Grn","Il15","Il16","Il18","Il18rap","Itgam","Itgb2","Itgb2l","Anxa1","Myd88","Myo1f","Pikfyve","Prkcd","Pla2g2a","Ptafr","Scnn1b","Ccl5","Spi1","Syk","Lypd11","Cd300a","Tnf","Tyrobp","Lypd10","Fcgr4","Pram1","Prg3","Cd177","Kmt2e","Plpp6","Stx11")
GOBP_NEUTROPHIL_ACTIVATION_INVOLVED_IN_IMMUNE_RESPONSE <- c("Abr","Bcr","Anxa3","Dnase1","Dnase1l3","Fcer1g","Itgam","Itgb2","Itgb2l","Myd88","Myo1f","Pikfyve","Ptafr","Scnn1b","Spi1","Syk","Lypd11","Tyrobp","Lypd10","Pram1","Cd177","Stx11")
NEUTROPHIL_MEDIATED_CYTOTOXICITY <- c("Nlrp6","Arg1","Ctsg","Dnase1","Dnase1l3","F2","F2rl1","Cxcl1","Myd88","Ncf1","Pomc","Scnn1b","Cxcl5","Elane","Trem1","Trem3","Tusc2")
PHAGOCYTOSIS_ENGULFMENT <- c("Ano6","Nckap1l","Rab31","Adgrb1","Abca1","Ager","Aif1","Alox15","C3","Cd36","Cdc42","Clcn3","Elmo1","F2rl1","Fcer1g","Fcgr1","Fcgr2b","Fcgr3","Gata2","Ighg2b","Ighg1","","Itga2","Itgam","Itgb2","Lbp","Rac3","Havcr1","Timd2","Marco","Myh9","Pparg","Sirpa","Rac1","Msr1","Sh3bp1","","Gm12169","Appl2","Cd300a","Thbs1","Xkr6","Treml4","Gsn","Xkr7","Arhgap25","Plcg2","Abca7","Timd4","Xkr5","Dppa1","","Xkr9","Xkr8","BC053393","Xkr4","Becn1","Clec7a","Stap1","Bin2","Megf10","Gulp1","Arhgap12","Siglece","Trem2")
neutro_degra <- c("Abr","Bcr","Anxa3","Itgam","Itgb2","Itgb2l","Myo1f","Pikfyve","Ptafr","Spi1","","Syk","Lypd11","Lypd10","Pram1","Cd177","Stx11")
neutro_kill_bact <- c("Nlrp6","Ctsg","F2","F2rl1","Myd88","Ncf1","Scnn1b","Cxcl5","Elane","Trem1","Trem3","Tusc2")
pos_reg_neutro_acti <- c("Il16","Itgam","Itgb2","Itgb2l","Ptafr","Lypd11","Tnf","Lypd10","Cd177","Plpp6")
neutro_reg_deg <- c("Abr","Bcr","Itgam","Itgb2","Itgb2l","Ptafr","Spi1","Syk","Lypd11","Lypd10","Pram1","Cd177")
cytotox <- c("Cd5l","Rab27a","C3","Cd59a","Cfh","Cr1l","Cd55","Hsp90ab1","Il13","Il4","Cd59b")
alt_compl <- c("C8b","Cfd","C3","C9","Cfh","Cr2","Cfb","Hc","Cfp","Cfhr4","C8a","Vsig4","C8g","Susd4")
HALLMARK_COMP <- c("Cdk5r1","Usp15","Mmp14","Gng2","Klkb1","Klk1","Scg3","F10","Fdx1","Calm3","Gnb4","Was","Dgkg","Adam9","Cfb","Anxa5","Dock9","Ctsh","Ctss","Notch4","Phex","Casp4","Casp3","Cd36","Kif2a","Tfpi2","Hnf4a","Actn2","Gzma","Gzmb","Mmp15","F7","Casp7","Ctsc","S100a13","C9","Pla2g4a","Tnfaip3","Rnf4","Mmp8","Psen1","Cd46","Pik3ca","Olr1","Prep","Casp9","Gzmk","Rasgrp1","Ppp2cb","Lcp2","Dyrk2","Lgmn","Gp1ba","Atox1","Zfpm2","Rce1","S100a9","Mmp13","Fcnb","Ehd1","Zeb1","Psmb9","Pclo","Pla2g7","Msrb1","Pik3cg","Sirt6","C1s1","Irf7","Gp9","Plek","Cpq","Akap10","Spock2","Ppp4c","Lap3","Dusp6","Cpm","Vcpip1","Dock4","Kynu","Gca","Tmprss6","Cda","Prcp","Prdm4","Usp16","Prss36","Rbsn","Gmfb","Rhog","Usp14","Kcnip3","Usp8","Kcnip2","Rabif","Ctso","Pcsk9","Dock10","Cblb","Brpf3","Hpcal4","Xpnpep1","Pik3r5","Dgkh","L3mbtl4","Dusp5","Pdp1","Adra2b","Ang","Apoa4","Apoc1","Serpinc1","C1qa","C1qc","C2","C3","Calm1","Car2","Cd40lg","Cebpb","F2","F3","F5","F8","Cfh","Clu","Col4a2","Cp","Cr2","Csrp1","Ctsb","Ctsd","Ctsl","Plscr1","Stx4a","Gngt2","Sh2b3","Serping1","Dpp4","Fcer1g","Fn1","Fyn","Gata3","Gnai2","Gnai3","Gnb2","Grb2","Hspa5","Hspa1a","Casp1","Il6","Irf1","Irf2","Itgam","Itih1","Jak2","Lamp2","Lck","Lgals3","Lipa","Lrp1","Lta4h","Ltf","Lyn","Maff","Mmp12","Me1","Mt3","Pdgfb","Pfn1","Pim1","Prkcd","Serpine1","Serpinb2","Plat","Plaur","Plg","Raf1","Ccl5","Src","Timp1","Timp2","Cdh13","Gpd2")
GOBP_COMPLEMENT_ACTIVATION <- c("C7","C8b","Cfd","Cd5l","Serping1","C1qa","C1qb","C1qbp","C1qc","C2","C3","C4b","C4bp","C6","C9","Cd59a","Cfh","Cfi","Cr2","Crp","Cr1l","Cd55","Cd55b","Fcna","Fcnb","Cfb","Hc","Ighg2b","Ighg1","Ighm","Il1b","Krt1","Masp1","Masp2","Mbl1","Mbl2","Cd46","Cfp","Phb","Cfhr4","Zp3r","C8a","A2m","C1rl","Igha","Colec10","Vsig4","C1s2","Cd59b","Ighe","Ighg3","Ighg2c","Cfhr1","C1s1","C1ra","Cfhr2","C4a","Svep1","Rgcc","C1rb","C8g","Colec11","Trem2","Susd4")
clas_compl <- c("C8b","Serping1","C1qa","C1qb","C1qbp","C1qc","C2","C3","C4bp","C9","Cfi","Cr2","Crp","Cr1l","Cd55","Cd55b","Hc","Ighg2b","Ighg1","Ighm","","Masp2","Mbl1","Mbl2","Cd46","Zp3r","C8a","C1rl","Igha","C1s2","Ighe","","","Ighg3","Ighg2c","C1s1","C1ra","Svep1","C1rb","C8g","Trem2","Susd4")
GOBP_CELLULAR_RESPONSE_TO_OXIDATIVE_STRESS  <- c("Adprs","Plekha1","Pex12","Etv5","Ppif","Prkaa1","Prkaa2","Foxp1","Selenos","Msra","Abl1","Pawr","Parp1","Aif1","Akt1","Abcd1","Akr1b3","Alox5","Ambp","Prdx3","Apex1","Apoa4","Aqp1","Rhob","Arnt","Arntl","Atf2","Atf4","Atm","Atp2a2","Atp7a","Bmp4","Bnip3","Btk","Capn1","Cat","Ccs","Cd36","Cdkn2a","Cfl1","Chuk","Coq7","Crygd","Cryge","Crygf","Cst3","Cyp1b1","Dhfr","Ect2","Edn1","Ednra","Egfr","Eif2s1","Endog","Epas1","Stx2","Chchd2","Ezh2","Fabp1","Fancc","Fer","Fos","Fxn","Fyn","G6pd2","G6pdx","Gch1","Rack1","Gpr37","Gpx5","Gsr","Hdac2","Hdac6","Hgf","Hif1a","Hmox1","Foxa1","Hsf1","Hspb1","Il6","Jun","Lcn2","Anxa1","Oxr1","Ppargc1b","Gpr37l1","Met","Mgat3","Mmp2","Mmp3","Mmp9","Mpo","Mpv17","Mt3","Nfe2l1","Nfe2l2","Ngfr","Nqo1","Nos3","Slc11a2","Ddr2","Nr4a2","Sqstm1","Prdx1","Reg3b","Pax2","Pdgfra","Pdgfrb","Pdk2","Prkcd","Prkd1","Pla2r1","Pml","Ppargc1a","Psap","Sirpa","Ptprk","Pex2","Pex5","Rad52","Rela","Ripk1","Sin3a","Slc1a1","Slc25a14","Slc8a1","Snca","Sod1","Sod2","Sod3","Sphk1","Src","Stat6","Pycr1","Stx4a","Ncoa7","Fancd2","Hk3","Txndc2","Agap3","Prdx2","Kdm6b","Tnfaip3","Trex1","Trp53","Trpc6","Txn1","Ucp1","Xbp1","Rbm11","Tbc1d24","Ermp1","Pcgf2","Mapkap1","Slc25a24","Zc3h12a","Sesn2","Cpeb2","Pyroxd1","Tmem161a","Mpv17l2","Fbln5","Mapk7","Prkra","Zfp277","Axl","Map2k4","Map3k5","Mapk1","Mapk13","Mapk3","Mapk8","Mapk9","Slc7a11","Ppia","Setx","Slc4a11","Rps3","Naglu","Trpa1","Trpm2","Ep300","Pjvk","Tldc2","Chchd2-ps","Keap1","Prkn","Scly","Ankzf1","Fut8","Prdx5","Pex14","Net1","Pdcd10","Foxo1","Foxo3","Ripk3","Mgst1","Ankrd2","Park7","Smpd3","Stk25","Trp53inp1","Ngb","Sirt2","Htra2","Tsc1","Arl6ip5","Rwdd1","Tmigd1","Brf2","Oser1","Lrrk2","Pex10","Romo1","Pnpla8","Aldh3b1","Trap1","Rnf146","Pink1","Zfp580","Pycr2","Vkorc1l1","Dapk1","Vrk2","Aifm2","Dhrs2","Pnpt1","Pex13","Chchd4","Prr5l","Nme8","Ddias","Lonp1","Atg7","Meak7","Atp13a2","Selenon","Ercc6l2","Srxn1","Mpv17l","Wnt16","Sirt1")
GOBP_GRANULOCYTE_ACTIVATION <- c("Traf3ip2","Abr","Bcr","Anxa3","Cxcr2","Ccr2","Camp","Ctsg","Dnase1","Dnase1l3","F2rl1","Fcer1a","Fcer1g","Grn","Il15","Il16","Il18","Il18rap","Itgam","Itgb2","Itgb2l","Anxa1","Myd88","Myo1f","Pikfyve","Prkcd","Pla2g2a","Ptafr","Scnn1b","Ccl5","Spi1","Stx4a","","Enpp3","Syk","Lypd11","Cd300a","Tnf","Tyrobp","Vamp2","Lypd10","Fcgr4","Pram1","Ighe","Prg3","Tex101","Cd177","Kmt2e","Plpp6","Stx11","Pi4k2a","Kars")
GOBP_INFLAMMASOME_MEDIATED_SIGNALING_PATHWAY <- c("Nlrp6","Btk","Casp4","Cd36","Ddx3x","Gbp2","Hspa8","Irgm1","Igtp","Kcnj8","Lamp2","Myd88","Prkd1","Ppp2ca","Eif2ak2","Ptpn22","Nlrp1a","Trim30a","Brcc3","Nlrp3","Dhx33","Tlr4","Tlr6","Trim31","Mavs","Gbp5","Mark4","Plcg2","Mapk8","Nlrc3","Brcc3dc","Irgm2","Mefv","Nek7","Gm12250","Nlrp1b","Sirt2","Zdhhc12","Pycard","Stmp1","Fbxl2","Atat1","Usp50","Cptp","Trem2","Csnk1a1","Trim11")
GOBP_INFLAMMATORY_RESPONSE  <- c("Msmp","Ccl21b","Cd200l2","Lrrc19","Mir883b","Akna","Smpdl3b","Nlrp6","Mir7116","Mir7578","Ttc39aos1","Cd276","Traf3ip2","Chil4","Bap1","Pld4","Ano6","Sharpin","Gpsm3","Ticam1","Ttbk1","Ffar4","Il1rl2","Lgals2","Olr1","Foxp1","Npy","Hyal3","Selenos","Cela1","Abr","Bcr","H2-Q6","H2-M10.3","Acp5","Chrna7","Adipoq","Acvr1","Ada","Adam17","Adam8","Adcy7","Adcy8","Adcyap1","Adora1","Adora2a","Adora2b","Adora3","Ager","Agt","Agtr1a","Agtr1b","Agtr2","Ahcyl","Ahr","Ahsg","Aif1","Akt1","Abcd1","Alox15","Alox5","Alox5ap","Tirap","Aoc3","Fabp4","Ap3b1","Cd5l","Apoa1","Apod","Apoe","App","Arnt","Atm","Slc7a2","Atrn","Bcl6","Bdkrb1","Bdkrb2","Bmp2","Bmp6","Bmpr1b","Bst1","Bpgm","Btk","C1qa","C3","C3ar1","C5ar1","Calca","Camk4","Casp1","Casp4","Casp6","Ccr6","Cd14","Cd24a","Cd28","Cd36","Cd44","Cd6","Cd68","Cd81","Cdh5","Cebpa","Cebpb","Cfh","Chil1","Chil3","Socs3","Clock","Clu","Cxcr2","Cxcr3","Ccr1","Ccr1l1","Ccr3","Ccr2","Ccr4","Ccr5","Ccr7","Camp","Cnr1","Cnr2","Cntf","Cr2","Crh","Crhbp","Crp","Csf1","Csf1r","Csrp3","Cst7","Ctla2a","Ctsc","Ctss","Celf1","Cx3cr1","Cyba","Cybb","Cyp19a1","Ddt","Ddx3x","Dhx9","Ackr1","Dnase1","Dnase1l3","Dpep1","Ecm1","S1pr3","Ednrb","Egfr","Elf3","Aimp1","Epha2","Ephb2","Ephb6","Ephx2","Esr1","Drosha","Mecom","Ext1","Ezh2","F2","F2r","F2rl1","F3","F8","Il25","Fanca","Fasn","Fcer1a","Fcer1g","Fcgr1","Fcgr2b","Fcgr3","Fem1a","Fn1","Fosl1","Fosl2","Fpr2","Fpr-rs3","Fpr-rs4","Fpr1","Fpr3","Tlr3","Fut4","Fut7","Fxr1","Gata3","Gbp2","B4galt1","Ggt1","Gja1","Gnat2","Lilrb4b","Lilrb4a","Cmklr1","Gpr33","Gpx1","Gpx2","Lpcat3","Nr3c1","Grn","Cxcl1","Gstp1","H2-Bl","H2-D1","H2-K1","","H2-M10.1","H2-M2","H2-M9","H2-Q1","H2-Q10","H2-Q2","H2-Q4","H2-Q7","H2-T10","H2-T22","H2-T23","H2-T3","Hc","Hck","Ptpn6","Hdac5","Foxf1","Hgf","Hif1a","Hk1","Hmgb1","Hmga1","Hmox1","Hp","Eif2ak1","Hspa8","Hspa4","Ndst1","Hyal1","Hyal2","Icam1","Ido1","Ier3","Irgm1","Cxcl10","Ifng","Ifngr1","Ifngr2","Igf1","Ighg2b","Ighg1","","","Igtp","Il10","Il10ra","Il12b","Il13","Il16","Il17a","Il17ra","Il18","Il18rap","Il1a","Il1b","Il1r1","Il1r2","Il1rap","Il1rn","Il18r1","Il2","Il2ra","Il4","Il4ra","Il5ra","Il6","Ins1","Ins2","Acod1","Isl1","Itgam","Itgav","Itgb1","Itgb2","Itgb2l","Itgb6","Cd47","Itih4","Jak2","Jun","Kcnj8","Kit","Klkb1","Klk1b1","Kng1","Kpna6","Krt16","Krt1","Lamp2","Large1","Lat","Lbp","Ldlr","Lep","Lgals1","Lgals9","Lipa","Loxl3","Anxa1","Lpl","Xcl1","Lrp1","Lta","Lxn","Tlr7","Tlr8","Cd180","Il1rl1","Ly86","Ly96","Lyn","Il17rc","Smad1","Smad3","Havcr2","Mapkapk2","Mas1","Cma1","Tpsb2","Mdk","Abcc1","Mep1b","Clec10a","Mif","Cxcl9","Mmp8","Cd200","Mvk","Myd88","Myo5a","Naip1","Naip2","Naip5","Naip6","Ncf1","Ndp","Ndufs4","Nfe2l1","Nfe2l2","Nfix","Nfkb1","Nfkb2","Nfkbia","Nfkbib","Ninj1","Nos2","Notch1","Notch2","Ccn3","Npy5r","Slc11a1","Orm1","Orm2","Osm","P2rx1","P2rx7","Reg3b","Pdcd4","Enpp1","Per1","Pik3cd","Prkca","Prkd1","Prkcq","Prkcz","Pla2g2d","Pla2g5","Plaa","Serpine1","Pld3","Serpinf2","Plp1","Ccl21a","Pmp22","Polb","Ppara","Ppard","Pparg","Ppp2ca","Eif2ak2","Proc","Psen1","Psen2","Psmb4","Pstpip1","Ptafr","Dicer1","Ptgdr","Ptger1","Ptger2","Ptger3","Stab1","Ptger4","Ash1l","Ptgfr","Ptgir","Ptgis","Hps1","Ptgs2","Ptn","Ptpn2","Ptpn22","Sirpa","Trpv1","Rasgrp1","Nlrp1a","Rb1","Rbpj","Reg3a","Reg3g","Rel","Rela","Relb","Ripk1","Rora","Rps19","Trim30a","Nr1h4","S100a8","S100a9","Saa1","Saa2","Saa3","Saa4","Scn9a","Scnn1b","Ccl1","Ccl11","Ccl12","Ccl17","Ccl2","Ccl20","Ccl22","Ccl25","Ccl3","Ccl4","Ccl5","Ccl6","Ccl7","Ccl8","Ccl9","Cxcl15","Cxcl2","Cxcl5","Cx3cl1","Selp","Sema7a","Foxp3","Snca","Snap23","Sod1","Sphk1","Serpina1a","Serpina1b","Serpina3n","Serpinb9","Spn","Pde2a","Cd200l1","Ghsr","Stat3","Stat5a","Stat5b","Ulk4","","Enpp3","Il23r","Syk","Metrnl","Nlrp4b","Brcc3","Lacc1","Fancd2","Il20rb","Tac1","Mylk3","Duoxa1","Tbxa2r","Slc18a2","Il36g","Il1f10","Sbno2","Appl2","Fem1al","Prdx2","Nlrp3","Dhrs7b","Kdm6b","Dhx33","Git1","Nr1d1","Cd300a","Pomt2","Tff2","Tgfb1","Thbs1","Timp1","Il31ra","Pxk","Tlr1","Tlr4","Tlr6","Tnc","Tnf","Tnfaip3","Tnfaip6","Tnfrsf11a","Tnfrsf1a","Tnfrsf1b","Tnfsf11","Pglyrp1","Cd40lg","Plscr1","Trex1","Trp73","Tnfrsf4","Tnfsf4","Tyro3","Tyrobp","Umod","Nr1h3","Scgb1a1","Vamp8","Vcam1","Vnn1","Vwf","Ccn4","Nrros","Wnt5a","Setd4","H2-M10.4","H2-M11","H2-M1","H2-M10.5","Trim31","Pja2","Hrh4","Ticam2","Afap1l2","Ywhaz","Zfp35","Zfp36","Camk1d","Zp3","Mavs","Syt11","Pbxip1","S100a7l2","Trim45","Chil5","Chil6","Rap1gds1","Gbp5","Zc3h12a","Themis2","Il22ra1","Nppa","Tnip2","Daglb","Cyp26b1","Klrh1","Lilra5","Nlrp2","Mark4","Cxcl17","Nlrp9a","Ffar2","Ffar3","Mrgpra3","Plcg2","Il17c","Il22ra2","Pla2g3","Agr2","Lrfn5","Ets1","Ggt5","Tlr11","Il17d","Adamts12","Mgll","Nt5e","Muc19","Nlrp5","Cd200r4","H2-M5","Scn11a","Ccl19","Sigirr","Tnfsf18","Tlr2","Scyl3","Zeb2","Cers6","Tspan18","Pde5a","Napepld","Fkrp","Nlrp9b","Nlrp4a","Nfkbid","Siglecg","Nlrp10","Mcph1","C2cd4a","Tarm1","Vps54","Il27","Il17f","Nod2","Axl","Spata2","Map2k3","Map3k7","Mapk14","Mapk8","Mapk9","Psma1","Cul3","Pla2g10","Abcd2","Nlrc3","Chst4","Nlrc4","Dagla","Ahcy","Pla2g2e","Nlrx1","Aoah","Irf5","Tcirg1","Cd200r2","Slamf1","Pla2g7","Naglu","Tlr13","Rhbdd3","Pik3cg","H2bc1","Gpr4","C5ar2","Smo","","Abi3bp","Fpr-rs6","Fpr-rs7","Parp4","Tafa3","Cxcl3","Nlrp9c","Card9","H2-M10.2","Trim65","Nr1d2","Brcc3dc","Nlrp12","Ighe","","Nr1h5","Trim55","S100a7a","Aim2","Tlr12","Mir21a","Mir147","Mir155","H2-M10.6","Spink7","Adcy1","Otulin","Gpr31b","Rab44","Nlrp4e","Elane","Il17rb","Il22","Pbk","Mfhas1","Macir","Stk39","Tslp","Tlr5","","Ccl26","Irf3","Ccrl2","Chst2","Irgm2","Il36a","Il36rn","Tollip","Mefv","Gm5849","Extl3","Cxcl13","Cxcl11","Il17b","Ccl24","Hdac7","Gps2","Nupr1","Tmed2","Socs5","Elf4","Rps6ka4","Clec7a","Clcf1","Rabgef1","Mtor","Pf4","Stap1","Brd4","Park7","Ppbp","Gpr17","Pglyrp2","Cd200r1","Tnip1","Il17re","Crlf2","Zbp1","Cysltr1","Ghrl","F12","Nampt","Nek7","Ackr2","Gm12250","Trav7-2","Nlrp1b","Trpv4","Dusp10","Ptges","Gprc5b","Fndc4","Sirt2","Nmi","Ndfip1","Vps35","Setd6","Acer3","Zdhhc12","Serpinb1a","Camk2n1","Lrrk2","H2-T-ps","Cntnap2","Duoxa2","Pycard","Tril","Mkrn2","Ctnnbip1","Cuedc2","C1qtnf12","Armh4","Slc39a8","Tbc1d23","Stmp1","Plgrkt","Wdr83","Wfdc1","Chid1","Pomgnt1","Arel1","Rtn4","Letmd1","Tab2","Zfp580","Tmem258","Gsdmd","Snx4","Odam","Dab2ip","Il36b","Tnfaip8l2","Wnk4","Ifi35","Unc13d","Pnma1","Tspan2","Tradd","Rarres2","Nkiras2","Seh1l","Fbxl2","Nkg7","Mir301","Mir324","Prcp","Sting1","Uaca","Appl1","Rps6ka5","Atat1","Ankrd42","Cyld","Cd200r3","Shpk","Slamf8","Usp50","C2cd4b","Abhd12","Il34","Gper1","Nlrp14","Chst1","Il33","Hnrnpa0","Rictor","Ak7","Scyl1","Hdac9","Lias","Cptp","Tusc2","Nfkbiz","Cxcr6","Chia1","C1qtnf3","Tlr9","Siglece","Il23a","Trem2","Pik3ap1","Nlrp4c","Jam3","Sucnr1","Hamp","Cd96","Kars","Tac4","Cd163","Csnk1a1","Trim11","Hmgb2","Nlrp4f","Stard7")
GOBP_INNATE_IMMUNE_RESPONSE  <- c("Ifi208","Wfdc17","Isg15","Defa30","Gm10358","Mndal","Defa35","Defa41","Defa40","Defa27","Defa37","Defa34","Gapdh-ps15","Gm3839","Ccl21b","Klrb1","Lrrc19","Mir511","Defa2","Smpdl3b","Epg5","Zfp683","AY761185","Mptx2","Oas1d","Slc15a4","Spon2","Gbp6","Nlrp6","Trim68","Mapkapk3","Ifi206","Cdc42ep2","Pik3r6","Ddx1","Pld4","Ticam1","Eprs","Lgr4","Il1rl2","Nod1","Adgrb1","Shmt2","Ankhd1","Irak2","Sp110","Npy","Fgb","C8b","H2-Q6","H2-Q7","Rps6ka3","H2-M10.3","D1Pas1","Oas1c","Actg1","Adam15","Adam8","Cfd","Parp1","Akap1","Akt1","Gfer","Tirap","Ang","Ang5","Ang2","Ap1g1","Ap3b1","Atg5","Apoa4","Apoe","App","Aqp4","Arf6","Arg1","Arg2","Rab27a","B2m","Phb2","Bcl10","Prdm1","Blk","Btk","Serping1","C1qa","C1qb","C1qbp","C1qc","C2","Ciita","C3","C4bp","C9","Calm1","Camk2a","Capg","Casp1","Casp4","Casp6","Casp8","Cav1","Cd14","Cd1d1","Cd2","Cd24a","Cd36","Cd84","Cd86","Cdc37","Cdc42","Cebpg","Cfh","Cfi","Ch25h","Chuk","Cited1","Coro1a","Ccr1","Camp","Cr2","Crcp","Crk","Crp","Cr1l","Csf1","Csf1r","Cx3cr1","Cyba","Cybb","Cyp27b1","Cd55","Cd55b","Dapk3","Daxx","Dhx15","Ddx3x","Dhx9","Defb1","Defb2","","Defa29","","","","Defa31","","","","","","","Defa25","Defa3","","Defa5","","","","","Dpp4","Drd2","Ear1","Ear2","Ereg","Esr1","Evl","Wfdc18","F2rl1","Colec12","Fadd","Fau","Fcer1g","Fcgr1","Fcna","Fcnb","Fga","Fgr","Flot1","Fosl1","Fosl2","Fpr2","Fpr-rs3","Fpr-rs4","Fpr3","Tlr3","Frk","Fyn","Xrcc6","Gapdh","Gata3","Gbp2b","Gbp2","Usp15","Gch1","Gfi1","Grn","Gzmb","Gzmc","Gzmd","Gzme","Gzmf","Gzmg","H2-Aa","H2-Ab1","Cfb","H2-Bl","H2-D1","H2-Eb1","H2-K1","","H2-M10.1","H2-M2","H2-M3","H2-M9","H2-Q1","H2-Q10","H2-Q2","H2-Q4","H2-Q7","H2-T10","H2-T22","H2-T23","H2-T3","Mr1","","Hc","Hck","Ptpn6","Hk1","Hmgb1","Hmgb3","Hpx","Hspa8","Hspd1","Hsp90aa1","Irf8","Irgm1","Ifi203","Ifi204","Ifit1","Ifit2","Ifit3","Ifna1","Ifnar1","Ifnar2","Ifnb1","Ifng","Ifngr1","Ifngr2","Igf2","Ighm","Jchain","Igtp","Cd74","Il12a","Il12b","Il12rb1","Il17a","Il17ra","Il18","Il18rap","Irak1","Il1rap","Il4","Ins1","Ins2","Irf1","Irf2","Irf4","Acod1","Itch","Cd47","Jak1","Jak2","Jak3","Kcnj8","Kif16b","Kif5b","Klrc1","Klrc2","Klrd1","Krt16","Krt1","Lag3","Lamp1","Lamp2","Arhgef2","Lbp","Lck","Lcn2","Lep","Lgals3","Lgals9","Anxa1","Xcl1","Lrp8","Ltf","Klrb1a","Klrb1c","Znrf1","Sertad3","Tlr7","Tlr8","Cd180","Ly86","Ly9","Ly96","Lyar","Lyn","Lyst","Havcr2","Slc26a6","Mapkapk2","Marco","Masp1","Masp2","Matr3","Mbl1","Mbl2","Mif","Mmp12","Mmp3","Gbp4","Clec4d","Mpeg1","Mrc1","Mx1","Mx2","Myd88","Myo1c","Myo1f","Naip2","Naip5","Naip6","Ncf1","Nck1","Ndufs4","Nedd4","Nfe2l2","Nfkb1","Nfkb2","Nfkbia","Nfkbil1","Ninj1","Nmbr","Nqo1","Cd244a","Nos2","Slc11a1","Prdx1","Pcbp2","Padi4","Pdpk1","Cfp","Phb","Pik3cd","Pik3r1","Prkce","Prkd1","Pla2g1b","Pla2g5","Pld3","Plscr2","Ccl21a","Bpifa1","Pml","Ppp1r14b","Cnot7","Med1","Pparg","Ppp2ca","Prkdc","Eif2ak2","Pstpip1","Wfdc12","Wfdc15b","Hexim1","Ptpn2","Ptpn22","Sirpa","Ripk2","Rfpl4","Ptprs","Ptx3","Nectin2","Rab12","Rab20","Rnf185","Raet1d","Rasgrp1","Nlrp1a","Trim40","Reg3g","Rel","Rela","Relb","Trim27","Trim10","Mst1r","Rps19","Trim30a","Nr1h4","S100a8","S100a9","Apcs","Ccl1","Ccl11","Ccl12","Ccl17","Ccl2","Ccl20","Ccl22","Ccl25","Ccl3","Ccl4","Ccl5","Ccl6","Ccl7","Ccl8","Ccl9","Cx3cl1","Spi1","Sftpd","Sh2d1a","Ptk6","Sin3a","Clpb","Slc22a5","Slpi","Smpd1","Snca","Sp100","Serpinb9b","Serpinb9c","Serpinb9f","Serpinb9e","Serpinb9","Serpinb9d","Dtx4","Src","Srebf1","Srms","Trim21","Znrf4","Stat1","Stat2","Stat5a","Stat5b","Hdac4","Slk","Aurkb","Stx4a","Stxbp1","Stxbp2","Stxbp3","Stxbp4","Dtx3l","Wfdc5","Trim30d","","Il23r","Syk","Nlrp4b","Brcc3","Lacc1","Trim41","Tifa","Pde12","Trim52","N4bp3","Ifitm6","Tac1","Tap1","Ythdf2","Tap2","Tarbp2","Trim38","Tfeb","Cfhr4","Cgas","Neurl3","Arid5a","Traf3ip3","Il36g","Ccdc92","9230019H11Rik","Rfpl4b","Appl2","Trim58","Nlrp3","Arrb2","Dhx33","Trim25","Nr1d1","Tmem106a","Cd300a","Nploc4","Cybc1","Ifi27l2b","Tgfb1","Tgtp1","Serinc5","Trim28","Polr3a","Tlr1","Tlr4","Tlr6","Ang4","Otop1","Tnf","Tnfaip3","Cd40","Pglyrp1","Traf2","Traf3","Traf4","Traf6","Plscr1","Trex1","Trp53","Rpl13a","Tnfsf4","Txk","Tyro3","Tyrobp","Ufd1","Umod","Nr1h3","Vamp3","Vamp8","Vav1","Vim","Vip","Vnn1","Lrrc14","Ttll12","Wap","Gramd4","Was","Nrros","Wnt5a","Ilrun","Marchf2","H2-M10.4","H2-M11","H2-M1","H2-M10.5","Trim31","Treml4","Pja2","Ticam2","Cd226","Tkfc","Xrcc5","Yes1","Ywhaz","Zap70","Rab7b","Ifi207","Ifi205","Trim26","Gigyf2","Gsn","Zyx","Zdhhc5","Mavs","Bpifb1","Ythdf3","Gbp5","Gbp7","Rigi","C8a","Tnfrsf14","Tnip2","Oasl1","Oas1e","Trafd1","Lrch4","A2m","C1rl","Klrb1f","Mark4","Nlrp9a","Ffar2","Ipo7","Ddx60","Plcg2","Zfp809","Tmem45b","Setd2","Ifi209","Gbp9","Sarm1","Banf1","Clec5a","Defa17","G3bp2","Tlr11","Colec10","Mid2","Oas1g","Oas1b","Oasl2","H2-M5","Malt1","Ccl19","Peli3","Tlr2","Ubd","Pglyrp3","Gpatch3","Oas1f","Klre1","Nlrp9b","Nlrp4a","Siglecg","Trim30b","Nlrp10","Triml1","Tarm1","Gzmn","Atg9a","Rbm47","Defb9","Defb11","Defb15","Defb35","Defb10","Fcgr4","Defb19","Oas3","Oas2","Oas1h","Oas1a","Cd300lf","Il27","Il17f","Nod2","Axl","Ceacam1","Map3k5","Map3k7","Map4k2","Mapk8","Defb36","Irak4","Zbtb1","Nlrc3","Clec4a2","Nlrc4","Sh2d1b1","Polr3e","Ecsit","Serinc3","Pla2g2f","Clec4a1","Ssc5d","Eif4e2","Klrk1","Nlrx1","G3bp1","Irf5","Bpifc","Slamf1","Clnk","Msrb1","Skp2","Naglu","Wfdc16","Vsig4","Tlr13","Cep63","Tomm70a","Flnb","Slamf6","Slc46a2","C1s2","Gimap5","H2bc12","H2bc21","Trim12c","Shfl","Defb20","Ptgs2os","Klri2","Fpr-rs6","Fpr-rs7","Usp44","Scimp","Bpi","Nlrp9c","Ifnl2","Card9","H2-M10.2","Trim65","Ifnl3","Morc3","Defb37","Myo18a","Defb34","Defb38","Defb39","Defb40","Brcc3dc","Ifi211","Clec4b2","Defa39","Defa22","Cd300ld3","Aim2","Tlr12","Trim56","Pglyrp4","H2-M10.6","Defb21","Wfdc13","Tnip3","Otulin","Akirin2","Trim34b","Trim30c","Nlrc5","Defb22","Nlrp4e","Clec4a4","Defa23","Defa24","Klri1","Zdhhc18","Ang5","Ear14","C1s1","C1ra","Klrg1","Rab11fip5","Mfhas1","Pvr","Ube2k","Vamp4","Nono","Tlr5","","Ear14","Irf7","Ccl26","Irf3","Rnase2b","Cd160","Irgm2","Unc93b1","Il36a","Il36rn","Tollip","Mefv","Serpinb9h","Sh2d1b2","Ifi214","Pqbp1","Usp27x","Crtam","Trim43a","Tyk2","Cadm1","Parp14","Gbp3","Stx8","Trim3","Samhd1","Ubqln1","Nagk","Ddx21","Ccl24","Rbm14","Gps2","Mettl3","Akap8","Syncrip","Adar","Tbk1","Ikbke","Tspan6","Slc22a21","Raet1d","Clec4e","Clec4n","Clec7a","Cdc42ep4","Ube2l6","U90926","Isg20","Slc15a2","Usp29","Tnip1","Klrc3","Rsad2","Zbp1","Trem1","Trem3","Usp14","Erbin","Nek7","Iigp1","Il21","Defb42","Triml2","Ifi213","","Gbp10","Defa28","Defa26","Wfdc9","Wfdc10","Wfdc11","Ang6","Gm12250","Gbp11","Defa38","Nlrp1b","Dusp10","Sirt2","Nmi","Slc15a3","Defb43","Defb25","Defb18","Cxcl16","Wfdc21","Ifitm3","S100a14","Zdhhc12","Serpinb1a","Tmem126a","Mptx1","Fam3a","Defa21","Senp7","Dus2","Defa42","Trim13","Pspc1","Nepn","Trim43c","Trim43b","Actr2","C1rb","H60b","Ifit3b","Ppp1r14bl","Ncr3-ps","H2-T-ps","Trim5","Ube2w","Pycard","Trim35","Tril","Riok3","Trim59","Polr3k","H60c","Polr3d","Lsm14a","Gp2","Bpifa5","Herc6","Arl8b","Peli1","Rpl39","Inava","Polr3g","Trim62","Atg12","Rnf125","Wfdc2","Stmp1","Ppp6c","Tmem33","Zcchc3","Hcfc2","Sdhaf4","Defa20","Chid1","Nmb","Nop53","Wfdc15a","Mcoln2","Mul1","Rtn4","Ifitm1","Rnf166","Cd177","Vps26b","Trim15","Gsdmd","Dhx16","Wdfy1","C8g","Ubl7","Bst2","Plekhm2","Dab2ip","Dapk1","Il36b","Tnfaip8l2","Trim32","Clec4b1","Rab43","Ifi35","Cactin","Polr3f","Polr3b","Unc13d","Ipo5","Kynu","Zdhhc1","Zdhhc11","Clec12b","Lrrfip2","Tasl","Alpk1","Sfpq","Fbxo9","Ifih1","Optn","Rarres2","Colec11","Wfdc3","Endod1","Rnf135","Cnpy3","Dhx36","Trim29","Fbxl2","Nkg7","Ttc4","Rps6kb1","Sting1","Appl1","Clec4a3","Tbkbp1","Atat1","Defb30","Irak3","Otud4","Actr3","Tmem43","Sec14l1","Uba7","Cyld","Usp20","Plpp6","Polr3c","Ifitm7","Stx11","Trim14","Slamf8","Usp38","Rab11fip2","Usp50","Rnf19b","Slamf7","Defb29","Garin5a","Gbp8","Rab2b","Rftn1","Il34","Trim12a","Calcoco2","Gper1","Defb41","Defb12","Dcst1","Ulbp1","Gpr108","Rnase6","Zc3hav1","Wrnip1","Polr3h","Cptp","Parp9","Apobec3","N4bp1","Klrb1b","Nfkbiz","Dhx58","Ifitm2","Pum1","Pum2","Kat5","Ankrd17","Tlr9","Gimap3","Il23a","Trem2","Pik3ap1","Nlrp4c","Dnaja3","Cd96","Trim8","Csnk1a1","Clec2d","Ear6","Ear10","Rnase2a","Ube2n","Serpinb9g","Trim6","Trim7","Trim11","Trim34a","Susd4","Hmgb2","Nlrp4f","Znfx1","Fgg")




th17 <- c("Mir873a","Ccr6","Il2","Il4","Il6","Il6ra","Irf4","Loxl3","Ly9","Smad7","Ascl2","Rora","Rorc","Ccl20","Foxp3","Stat3","Nlrp3","Tgfb1","Zbtb7b","Zc3h12a","Malt1","Tnfsf18","Nfkbid","Slamf6","Rc3h2","Rc3h1","Batf","Otud5","Tbx21","Mir301","Mir326","Nfkbiz","Il23a","Entpd7")

Tissue_residency <- (c("2900026A02Rik","AA467197","Abcb1b","Adam19","Atf3","B4galnt4","Bag3","Bhlhe40","Btg2","Btg3","Ccl4","Cd244","Cd69","Cdh1","Chd7","Chn2","Cish","Coq10b","Cpd","Crem","Ctla4","Dennd4a","Dgat1","Dnaja1","Dnaja4","Dnajb1","Dnajb4","Dnajb6","Dusp1","Dusp4","Dusp6","Egr1","Ehd1","Ell2","Fasl","Fbxo30","Fgl2","Fmnl3","Fos","Fosl2","Gabarapl1","Gadd45b","Gem","Glrx","Gpr171","Gzmb","H3f3b","Havcr2","Hip1","Hpgds","Hsp90aa1","Hspa5","Icos","Ifng","Ifrd1","Il21r","Il4ra","Inpp4b","Irf4","Isg20","Isy1","Itga1","Jun","Junb","Kdm6b","Litaf","Lmnb1","Ly6g5b","Mxd1","Mxi1","Neurl3","Nfil3","Nfkbie","Nr4a1","Nr4a2","Nr4a3","Ntan1","Odc1","P2ry10","P4hb","Pdcd1","Per1","Plk3","Plscr1","Pmepa1","Ppp1r15a","Prdx6","Ptp4a1","Qpct","Rbpj","Rgs1","Rgs10","Rgs16","Rhob","Rnf149","Rrad","Sik1","Skil","Slc16a6","Slc7a5","Spty2d1","Ssbp2","Stk17b","Tgif1","Tigit","Tiparp","Tjp1","Tnf","Tnfaip3","Tnfrsf1b","Tnfrsf9","Tra2a","Trib1","Trp53inp2","Ube2s","Vps37b","Wsb1","Xcl1","Zfand5","Zfp36l1"))

Tissue_repair1 <- (c("Flg","Flg2","Furin","Areg","Bmp2","Bmp4","Bmp7","Ngf","Tgfa","Tgfb1","Tgfb2","Tgfb3","Pdgfa","Pdgfb","Csf2","Fgf1","Fgf7","Fgf10","Fgf22","Igf1","Igf2","Vegfa","Ctgf","Mmp3","Mmp10","Mmp13","Mmp15","Mmp25","Mmp28","Defb1","Defb6","Epgn","Hbegf","Egf","Mmp2","Tnf","Il1a","Il1b","Csf1","Ccl2","Ccl3","Cxcl1","Cxcl2","Cxcl10","Cxcl15","Il6","Pgf","Angpt1","Angpt2","Angpt4","Hgf","Cyr61","Inhba","Inhbb","Btc","Ereg","Nrg1","Nrg2","Nrg4","Nrg3","Lep","Vegfb","Vegfc","Vegfd","Pdgfc","Hmgb1","Ptges2","Cxcl12","Shh","Dhh","Ihh","Disp1","Dll1","Dll3","Dll4","Jag1","Jag2","Wnt1","Wnt2","Wnt2b","Wnt3","Wnt3a","Wnt4","Wnt5a","Wnt5b","Wnt6","Wnt7a","Wnt7b","Wnt8a","Wnt8b","Wnt9a","Wnt9b","Wnt10a","Wnt10b","Wnt11","Wnt16","Chat","Thbs1","Hif1a"))

Tissue_repair2 <- (c("Pcyt1a","Pcsk1","Ccr1","Areg","Frmd5","Havcr2","Slc15a3","Ccr3","Il10rb","Lyn","Npnt","Ctsh","Arnt2","Il23r","Tubb6","Neb","Snx9","Sik1","Rgs2","Dgat2","Fam46a","Trf","Kcna4","Arl5a","Tnfrsf10b","Padi2","Gpr55","Sdc4","Cd200r1","Pparg","Camk2n1","Ehd4","Tmbim1","Ifrd1","Serpinb6a","Plxnd1","Rrad","Rasgef1b","Plin2","Msrb3","Per1","Sgms1","Dusp14","Spty2d1","Myo1d","Coq10b","Il9r","Evi5","Cables1","Cd80","Cpm","Ehd1","Nr1d1","Alcam","Epcam","Arl5b","Ets2","Acot11","Nebl","Klf4","Klf8","Slc25a19","Cnnm2","L1cam","Gprc5a","Sytl3","Tagln2","Axl","Zdhhc23","Ky","Ttbk2","Fam129b","Trim46","Gla","Gm2a","Cd74","Clic4","Slc16a3","Dusp3","Pdgfb","Itgav","Tdpoz4","Rnf125","Lgmn","Crem","Tnfrsf9","Rassf6","Bmpr1a","N4bp1","Rab11fip1","Zc3h12c","Got1","Hpse","Dot1l","Zfp36","Icam1","Rab4a","Plxdc1","Rorc","Asns","Irak3","Kdm6b","Ero1l","Fos"))

Tissue_repair3 <- (c("Il18r1","Il17re","Itgb1","Ramp3","Tnfsf11","Cd40lg","Plcb4","B3gnt5","Cep290","Tmem176b","Plxnc1","Tmem176a","Lif","Nrp1","Itga7","Colq","Sema4f"))

Tissue_repair4 <- (c("Adrb2","Aimp1","Aqp3","Arhgap1","Arnt","B4galt1","Bcl11b","Bub3","Cap2","Casp1","Cd151","Cd9","Cdkn1b","Cebpg","Csf2","Ctnnb1","Epb41","Fam129b","Fgfr1op2","Fhl2","Flna","Hmox2","Icam1","Icos","Il22","Il33","Ilk","Itga3","Itga9","Itgb1","Jun","Junb","Klf10","Klf4","Klk8","Lgals1","Map3k1","Map3k5","Mfge8","Mgat5","Mif","Myc","Myd88","Ncoa3","Nf1","Nfe2l2","Nfic","Nr4a1","Panx1","Pkd1","Ppard","Prkca","Prkce","Prkch","Ptk2","Ptk2b","Rac1","Rac2","Raf1","Ramp1","Ripk3","Sdc1","Sdc4","Siah2","Smad4","Smad7","Sod1","Stat3","Tcf3","Tgfb1","Thra","Thy1","Tnf","Tnfsf14","Trp53","Vav1","Vav3","Vegfa","Vim","Slc6a6","Itgb2","Hif1a","Clic4","Akt1","Anxa1","Dst","Egr1","Kit","Plaur","Pten","Il4","Ncoa6","Cd44"))

Circulating <- (c("1700025G04Rik","A430078G23Rik","Acss1","Aff3","Apobec2","Arhgap26","Arhgef18","Armc7","As3mt","Aven","B3gnt5","Bin2","Car5b","Cdc14b","Cmah","D1Ertd622e","Dapl1","Dennd2d","Dgka","Dock5","Dusp7","Ehd3","Elovl7","Eomes","F2rl2","Fam49a","Fam53b","Fgf13","Flna","Gab3","Gramd4","Gzmm","Haao","Icam2","Il17ra","Il18rap","Kbtbd11","Kcnj8","Klf2","Klhdc1","Klrb1c","Klrg1","Lair1","Lfng","Limd1","Lpin1","Ly9","Mboat1","Ms4a4c","Nck2","Nfatc3","Nod1","Nsg2","Pde2a","Prss12","Prune1","Pycard","Qprt","Rap1gap2","Rasa3","Rasgrp2","Rgs19","S1pr1","S1pr4","Samd3","Sell","Sfxn3","Sh2d1a","Sh3bp5","Sidt1","Sike1","Slamf6","Smpdl3b","Spn","St3gal1","Stard10","Tcf7","Tlr1","Tmem71","Tnfrsf26","Traf3ip3","Ubxn2b","Usp33","Ypel1","Zeb2","Zfp595","Zfp760"))

iT_17 <- (c("Rorc","Itgae","Il23r","Serpinb1a","Mmp25","1700040D17Rik","Pxdc1","Il1r1","Lingo4","Ccr8"))

iT_cycling_S <- (c("Actb","Ptma","Slc25a5","Ppia","Hmgn2","Lgals1","Cfl1","Gapdh","Ran","H2afz","Ldha","Eif5a","Hmgb1","Npm1","Actg1","Ybx1","Hmgb1","Ranbp1","Set","Mif"))

iT_cycling_G2M <- (c("Tuba1b","Stmn1","Tubb5","Ptma","Pclaf","H2afz","Hmgb2","Hmgn2","Hmgb1","Cks1b","Top2a","Dut","Hnrnpab","Lmnb1","Dek","Ran","Ppia","Smsc2","Hist1h2ap","Ranbp1"))

TCR_signaling <- str_to_title(c("AKT1","AKT2","AKT3","ARAF","BRAF","CD247","CD28","CD3D","CD3E","CD3G","CD80","CD86","CDC42","CHUK","CSK","FOS","GRAP2","HRAS","IKBKB","ITPR1","JUN","LAT","LCK","LCP2","MAP2K1","MAP2K2","MAP3K1","MAPK1","MAPK3","MAPK8","MAPK9","NCK1","NCK2","NFKBIA","NRAS","PAK1","PAK2","PAK3","PIK3C3","PIK3CA","PIK3CB","PIK3CD","PIK3CG","PIK3R1","PIK3R2","PIK3R3","PLCG1","PPP3CA","PPP3CB","PPP3CC","PRKCQ","PTPRC","RAC1","RAF1","SOS1","SOS2","VAV1","VAV2","VAV3","WAS","ZAP70","WAS","ZAP70"))

Fas_pathway <- (c("Arhgdib","Casp3","Casp6","Casp7","Casp8","Cflar","Daxx","Dffa","Dffb","Fadd","Faf1","Fas","Fasl","Jun","Lmna","Lmnb1","Lmnb2","Map2k4","Map3k1","Map3k7","Mapk8","Pak1","Pak2","Parp1","Prkdc","Ptpn13","Rb1","Sptan1"))

Il1R_signaling <- (c("Chuk","Ecsit","Ifna1","Ifnb1","Ikbkb","Il1a","Il1b","Il1r1","Il1rap","Il1rn","Il6","Irak1","Irak2","Irak3","Jun","Map2k3","Map2k6","Map3k1","Map3k14","Map3k7","Mapk14","Mapk8","Myd88","Nfkb1","Nfkbia","Rela","Tgfb1","Tgfb2","Tgfb3","Tnf","Traf6"))

ISG <- (c("Gbp5","I830012O16Rik","Ifit3","Ifi205","Ifit2","Cxcl10","Ifit1","Rsad2","Ifitm3","E030037K03Rik","Ifi204","Gbp3","E430029J22Rik","Iigp1","Tgtp1","Mx1","Slfn5","Tgtp1","Tnfsf10","Cmpk2","Slfn9","Mx2","Usp18","BC023105","Gbp2","Ifi44","Herc6","Ifi47","Isg20","AI607873","Ppm1k","Parp12","Gbp4","Fcgr1","Cd69","Gbp6","Pnp","Slfn1","Nlrc5","Mmp13","Nt5c3","Pnpt1","Irg1","Ms4a4c","Pyhin1","Oas3","Ifih1","Setdb2","BC094916","Oasl2","Nlrc5","Oas2","Oasl1","Pnp2","Eif2ak2","2010002M12Rik","D14Ertd668e","Ddx60","Tmem67","Irf7","Tagap","Ifi203","Rtp4","Nlrc5","Acer3","A530064D06Rik","Gpr141","Zufsp","Ms4a4a","Gdap10","Ly6i","Il18"))
