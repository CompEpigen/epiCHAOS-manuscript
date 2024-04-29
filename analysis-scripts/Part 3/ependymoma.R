
#--- epiCHAOS analysis of scATAC-seq data from ependyoma in 5 primary and 2 metastases from https://doi.org/10.1038/s41467-022-31683-9

library(Seurat)
library(Signac)
library(EnsDb.Hsapiens.v86)
library(magrittr)
library(stringr)


set.seed(10)

data.dir <- "/omics/groups/OE0219/internal/KatherineK/data/scATAC/Ependymoma/"

#--- scATAC from ependymoma
epn <- readRDS(file.path(data.dir, "GSE206579_snATAC_PF_EPN_Seurat.rds"))

#--- check umap
DimPlot(epn)

#--- set group for epiCHAOS calculation
epn@meta.data$group <- epn@meta.data$celltype <- Idents(epn)
epn@meta.data$group %>% table()

#--- peaks-by-cells matrix
atac <- epn@assays$ATAC@counts %>% na.omit()
dim(atac)

#--- populate list of peaks-by-cells matrices for each cell type
datasets <- list()
for (i in unique(epn$group)) {
  
  #--- cell ids for the selected cell type
  ids <- epn@meta.data[epn$group==i,] %>% rownames()
  
  #--- skip if too few cells, subset to 200 if many cells
  if (length(ids)<20) { next }
  temp <- atac[,ids]
  if (ncol(temp)>200) { temp <- temp[,sample(ncol(temp), 200)]}
  
  #--- adjust rownames
  rownames(temp) <- rownames(temp) %>% str_replace("-", ":")
  
  #--- binarise
  temp[temp>1] <- 1
  
  datasets[[i]] <- temp
}

lapply(datasets, dim)


#--- compute epiCHAOS scores per celltype
het <- compute.eITH(datasets)
het <- compute.eITH.cancer(datasets)
#het %>% write.csv("/omics/groups/OE0219/internal/KatherineK/ATACseq/epiCHAOS-supplementary-data/epiCHAOS_ependymoma.csv")

#--- add epiCHAOS scores to epn metadata
for (i in unique(het$state)) {
  epn@meta.data$epiCHAOS[epn@meta.data$celltype==i] <- het$mean.het[het$state==i]
}

#--- adjust labels for plotting
het$malignant <- ifelse(het$state %in% c("Undifferentiated_cells", "NPCs", "Ependymal_cells", "MLCs", "Astrocytes"), "malignant", "non-malignant")
het$state <- het$state %>% str_replace_all("_", " ")
Idents(epn) <- epn@meta.data$celltype %>% str_replace_all("_", " ")

library(RColorBrewer)
library(ggpubr)

#--- umaps and barplots of epiCHAOS scores per-celltype
p1 <- DimPlot(epn) + scale_color_brewer(palette = "Set3") +labs(x="UMAP1", y="UMAP2")

p2 <- FeaturePlot(epn, features="epiCHAOS")  + labs(x="UMAP1", y="UMAP2") + scale_color_distiller(palette = "Blues", direction = 1)

p3 <- ggplot(temp, aes(x=reorder(state, mean.het), y=mean.het, fill=malignant)) +
  geom_bar(position="dodge", stat = "identity", alpha=0.9, width = 0.6) +
  scale_fill_manual(values = c("steelblue4", "grey20"))+
  labs(x="", y="epiCHAOS")+
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90))


ggarrange(p1,p2,p3, ncol=3, widths = c(1, 0.8, 0.8))

