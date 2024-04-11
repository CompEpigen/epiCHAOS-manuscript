
#--- epiCHAOS analysis of scATAC-seq data from ependyoma in 5 primary and 2 metastases from https://doi.org/10.1038/s41467-022-31683-9

library(Seurat)
library(Signac)
library(EnsDb.Hsapiens.v86)
library(magrittr)
library(stringr)


set.seed(10)

#--- scATAC from ependymoma
epn <- readRDS("/omics/groups/OE0219/internal/KatherineK/data/scATAC/Ependymoma/GSE206579_snATAC_PF_EPN_Seurat.rds")

# check umap
DimPlot(epn)
DimPlot(epn, split.by = "biological_replicates")

# compute eITH per cluster
epn@meta.data$celltype <- Idents(epn)

epn@meta.data$group <- paste0(epn@meta.data$celltype, "-", epn@meta.data$biological_replicates)
epn@meta.data$group <- epn@meta.data$celltype

epn@meta.data$group %>% table()

atac <- epn@assays$ATAC@counts %>% na.omit()

dim(atac)

# populate list of datasets for each cell-type-sample pair
datasets <- list()
for (i in unique(epn$group)) {
  ids <- epn@meta.data[epn$group==i,] %>% rownames()
  if (length(ids)<20) { next }
  temp <- atac[,ids]
  if (ncol(temp)>200) { temp <- temp[,sample(ncol(temp), 200)]}
  rownames(temp) <- rownames(temp) %>% str_replace("-", ":")
  temp[temp>1] <- 1
  datasets[[i]] <- temp
}

lapply(datasets, dim)


# compute eITH
het <- compute.eITH(datasets)
het <- compute.eITH.cancer(datasets)
#het$state <- het$state %>% str_replace("_", " ")

#--- check associations with metastasis
epn@meta.data$origin <- ifelse(epn@meta.data$biological_replicates %in% c("M7", "M8"), "metastasis", "primary")

DimPlot(epn, group.by = "origin")
DimPlot(epn, group.by = "celltype")

# calculate average eITH score per population
for (i in unique(het$state)) {
  epn@meta.data$epiCHAOS[epn@meta.data$celltype==i] <- het$mean.het[het$state==i]
}


het$malignant <- ifelse(het$state %in% c("Undifferentiated_cells", "NPCs", "Ependymal_cells", "MLCs", "Astrocytes"), "malignant", "non-malignant")
het$state <- het$state %>% str_replace_all("_", " ")
Idents(epn) <- epn@meta.data$celltype %>% str_replace_all("_", " ")

library(RColorBrewer)
library(ggpubr)
p1 <- DimPlot(epn) + scale_color_brewer(palette = "Set3") +labs(x="UMAP1", y="UMAP2")
p2 <- FeaturePlot(epn, features="epiCHAOS")  + labs(x="UMAP1", y="UMAP2") + scale_color_distiller(palette = "Blues", direction = 1)
p3 <- ggplot(temp, aes(x=reorder(state, mean.het), y=mean.het, fill=malignant)) +
  geom_bar(position="dodge", stat = "identity", alpha=0.9, width = 0.6) +
  scale_fill_manual(values = c("steelblue4", "grey20"))+
  labs(x="", y="epiCHAOS")+
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90))

# pdf("umaps_and_barplots_eICH_vs_celltypes.pdf", 13, 3)
# ggarrange(p1,p2,p3, ncol=3, widths = c(1, 0.8, 0.8))
# dev.off()

pdf("/omics/groups/OE0219/internal/KatherineK/ATACseq/epiCHAOS-Figures/Figure 3/umap_ependymoma_celltypes.pdf", 6, 4)
p2
dev.off()

svg("/omics/groups/OE0219/internal/KatherineK/ATACseq/epiCHAOS-Figures/Figure 3/barplot_ependymoma_celltypes.svg", 6, 3.5)
p3
dev.off()


#--- try clustering the cells and performing eITH calculation on the clusters

epn <- FindNeighbors(object = epn, reduction = 'lsi', dims = 2:20)
epn <- FindClusters(object = epn, verbose = FALSE, algorithm = 3)
DimPlot(object = epn, label = TRUE) + NoLegend()
epn@meta.data$seurat_clusters

datasets <- list()
for (i in unique(epn$seurat_clusters)) {
  ids <- epn@meta.data[epn$seurat_clusters==i,] %>% rownames()
  if (length(ids)<20) { next }
  temp <- atac[,ids]
  if (ncol(temp)>100) { temp <- temp[,sample(ncol(temp), 100)]}
  datasets[[i]] <- temp
}

lapply(datasets, dim)
names(datasets) <- paste0("cluster ", names(datasets))

# compute eITH
het <- compute.eITH(datasets)


# add eITH score to sample annotation
for (i in unique(epn$seurat_clusters)) {
  epn@meta.data$eICH[epn@meta.data$seurat_clusters==i] <- het$mean.het[str_remove(het$state, "cluster ")==i]
}

p1 <- DimPlot(epn, group.by = "celltype", label = T) + NoLegend() + labs(x="UMAP 1", y = "UMAP 2")
p2 <- FeaturePlot(epn, features = "eICH") + 
  scale_colour_gradientn(colours = rev(brewer.pal(n = 9, name = "RdBu"))) + labs(x="UMAP 1", y = "UMAP 2")

pdf("ependymoma/umaps_epn_celltype_vs_eICH.pdf", 9, 4)
ggarrange(p1,p2, widths = c(4,5))
dev.off()


Idents(epn) <- epn$celltype
FeaturePlot(epn, features = "eICH", split.by = "origin", label=T) &
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))

