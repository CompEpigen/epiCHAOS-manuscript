
#--- test epiCHAOS on scTAM-seq from mouse hematopoiesis (Scherer et al. 2024)

library(Seurat)
library(pheatmap)

set.seed(10)

source("/omics/groups/OE0219/internal/KatherineK/ATACseq/epiCHAOS-analysis-scripts/epiCHAOS-functions/compute_epiCHAOS.R")

#--- load scTAM data as seurat object downloaded from figshare
scTAM <- readRDS("/omics/groups/OE0219/internal/KatherineK/data/scTAM/seurat_for_figshare.rds")

#--- extract single-cell methylation data as matrix
mat <- scTAM@assays$DNAm@counts %>% as.matrix()
mat[1:10,1:10]

#--- create matrix per cell type for epiCHAOS calculation
datasets <- list()
for (i in unique(scTAM@meta.data$CellType)) {
  ids <- scTAM@meta.data[scTAM@meta.data$CellType==i,] %>% rownames()
  if (length(ids)>200) { 
  ids <- ids %>% sample(200, replace = F)
  }
  datasets[[i]] <- mat[,ids]
}

lapply(datasets, dim)


#--- repeat with subsampling approach, takin pseudoreplicates of 5 x 100 cells for clusters which have more cells
datasets <- list()
for (i in unique(scTAM@meta.data$CellType)) {
  ids <- scTAM@meta.data[scTAM@meta.data$CellType==i,] %>% rownames()
  print(length(ids))
  
  if (length(ids)>100) { 
    for (n in 1:5) {
      datasets[[paste0(i, "-", n)]] <- mat[,sample(ids, 100)]
    }
  } else {datasets[[i]] <- mat[,ids]}
}

lapply(datasets, dim)


#--- compute epiCHAOS scores
het <- compute.eITH(datasets)
het$group <- het$state %>% str_split("-") %>% lapply("[", 1) %>% unlist()
het$group[het$group=="pre/pro"] <- "pre/pro-B"

#--- assign mean of subsamples
for (i in unique(het$group)) {
  het$mean[het$group==i] <- mean(het$mean.het[het$group==i])
}


# save scores
#saveRDS(het, "epiCHAOS_scores.Rds")
#het %>% arrange(desc(het)) %>% write.csv("/omics/groups/OE0219/internal/KatherineK/ATACseq/epiCHAOS-supplementary-data/epiCHAOS_scTAMseq_hemato.csv")


#--- add annotation for epiCHAOS scores to metadata
for (i in unique(het$state)) {
  scTAM@meta.data$epiCHAOS[scTAM@meta.data$CellType==i] <- het$mean[het$group==i]
}

#--- plot cell types and heterogeneity scores in umaps & barplots
p1 <- DimPlot(scTAM, label=T, label.size = 3) + scale_color_brewer(palette = "Set3") + labs(x="UMAP1", y="UMAP2")
p2 <- FeaturePlot(scTAM, features = "epiCHAOS") + scale_color_distiller(palette = "Blues", direction = 1, values = c(-0.1, 1)) + labs(x="UMAP1", y="UMAP2")

p3 <- ggplot(het, aes(x = reorder(state, mean.het), y = mean.het)) +
  geom_bar(stat="identity", position = "dodge", alpha=0.7, width = 0.5, fill="black")+
  labs(x="", y="epiCHAOS")+
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

pdf("umaps_hemato_scTAM_epiCHAOS.pdf", 9, 3.5)
ggarrange(p1,p2, widths = c(1,0.9), align = "h")
dev.off()
