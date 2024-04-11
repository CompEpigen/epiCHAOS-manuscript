#--- test epiCHAOS on scTAM-seq

setwd("/omics/groups/OE0219/internal/KatherineK/ATACseq/scTAM")

library(Seurat)
library(pheatmap)

source("/omics/groups/OE0219/internal/KatherineK/ATACseq/eITH-test-scripts/jaccard.R")

scTAM <- readRDS("/omics/groups/OE0219/internal/KatherineK/data/scTAM/seurat_for_figshare.rds")

mat <- scTAM@assays$DNAm@counts %>% as.matrix()
mat[1:10,1:10]

datasets <- list()
for (i in unique(scTAM@meta.data$CellType)) {
  ids <- scTAM@meta.data[scTAM@meta.data$CellType==i,] %>% rownames()
  #if (length(ids)>500) { 
  ids <- ids %>% sample(100, replace = F)
  #}
  datasets[[i]] <- mat[,ids]
}

lapply(datasets, dim)

het <- compute.eITH(datasets)

ggplot(het, aes(y=mean.het, x=reorder(state, mean.het))) +
  geom_dotplot(color="black")+
  labs( x="", y="heterogeneity score") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

for (i in unique(het$state)) {
  scTAM@meta.data$epiCHAOS[scTAM@meta.data$CellType==i] <- het$mean.het[het$state==i]
}


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
