

#--- testing epiCHAOS score on increasing numbers of cell types from Granja 2019 scATAC-seq data from bone marrow (relating to Figure 1E)

library(Seurat)
library(Signac)
library(stringr)
library(magrittr)

#--- load scATAC from normal hematopoiesis from Granja et al. 2019

atac <- readRDS("/omics/groups/OE0219/internal/KatherineK/data/scATAC/Granja2019/scATAC-Healthy-Hematopoiesis-191120.rds")
atac <- atac[,grepl("BMMC", atac$Group)]

# #--- identify cell type markers
# 
# #--- create as a signac object and find marker genes of each cell type
# counts <- atac@assays$data$counts
# rownames(counts) <- names(atac@rowRanges)
# 
# chrom_assay <- CreateChromatinAssay(
#   counts = counts,
#   sep = c("_", "-"),
#   min.cells = 10,
#   min.features = 200
# )
# 
# obj <- CreateSeuratObject(counts = chrom_assay, assay="peaks", meta.data = data.frame(atac@colData))
# 
# #--- find differentially accessible peaks for each celltype
# Idents(obj) <- obj@meta.data$BioClassification
# markers <- FindAllMarkers(obj)
# saveRDS(markers, file = "hemato_celltype_markers.Rds")
# 

#--- repeat to find marker peaks for all clusters (intensive, run as bsub)
markers <- readRDS("hemato_celltype_markers.Rds")

#--- select cell types to use for creation of "mixtures"
select.celltypes <- c("11_CD14.Mono.1", "17_B", "01_HSC", "24_CD8.CM", "09_pDC")

#--- select top 500 markers per selected cell type
top.markers <- list()
for (i in select.celltypes) {
  top <- markers[markers$cluster==i & abs(markers$avg_log2FC>2),] %>% arrange(p_val_adj) %>% head(500) %>% pull(gene) %>% str_replace_all("-", "_") 
  top.markers[[i]] <- gsub("\\..*","", top)
}

names(top.markers)
length(unlist(top.markers))
head(top.markers[[1]])

#--- create a list of datasets of cell type "mixtures"
datasets <- list()
cells.select <- c()

#--- individual cell types
for (i in 1:length(select.celltypes)) {
  hema <- atac@assays$data$counts[,atac$BioClassification == select.celltypes[i]]
  rownames(hema) <- names(atac@rowRanges)
  cells <- hema[rownames(hema) %in% unlist(top.markers),] %>% colSums() %>% sort() %>% tail(100) %>% names()
  cells.select <- c(cells.select, cells)
  hema <- hema[rownames(hema) %in% unlist(top.markers),cells]
  hema[hema>1] <- 1
  datasets[[select.celltypes[i]]] <- hema
}

lapply(datasets, dim)

#--- mix 2 celltypes
comparisons <- combn(c(1,2,3,4,5), 2)
for (i in 1:ncol(comparisons)) {
  hema <- atac@assays$data$counts[,atac$BioClassification %in% select.celltypes[comparisons[,i]]]
  hema <- hema[,colnames(hema) %in% cells.select]
  rownames(hema) <- names(atac@rowRanges)
  hema <- hema[rownames(hema) %in% unlist(top.markers),sample(ncol(hema), 100)]
  hema[hema>1] <- 1
  datasets[[paste0("mix2-", i)]] <- hema
}

#--- mix 3 celltypes
comparisons <- combn(c(1,2,3,4,5), 3)
for (i in 1:ncol(comparisons)) {
  hema <- atac@assays$data$counts[,atac$BioClassification %in% select.celltypes[comparisons[,i]]]
  hema <- hema[,colnames(hema) %in% cells.select]
  rownames(hema) <- names(atac@rowRanges)
  hema <- hema[rownames(hema) %in% unlist(top.markers),sample(ncol(hema), 100)]
  hema[hema>1] <- 1
  datasets[[paste0("mix3-", i)]] <- hema
}

#--- mix 4 celltypes
comparisons <- combn(c(1,2,3,4,5), 4)
for (i in 1:ncol(comparisons)) {
  hema <- atac@assays$data$counts[,atac$BioClassification %in% select.celltypes[comparisons[,i]]]
  hema <- hema[,colnames(hema) %in% cells.select]
  rownames(hema) <- names(atac@rowRanges)
  hema <- hema[rownames(hema) %in% unlist(top.markers),sample(ncol(hema), 100)]
  hema[hema>1] <- 1
  datasets[[paste0("mix4-", i)]] <- hema
}

#--- mix 5 celltypes
hema <- atac@assays$data$counts[,atac$BioClassification %in% select.celltypes[1:5]]
hema <- hema[,colnames(hema) %in% cells.select]
rownames(hema) <- names(atac@rowRanges)
hema <- hema[rownames(hema) %in% unlist(top.markers),sample(ncol(hema), 100)]
hema[hema>1] <- 1
datasets[["mix-5"]] <- hema

lapply(datasets, dim)
names(datasets)[1:5] <- c("Monocyte", "B cell", "HSC", "CD8T", "DC")

#--- compute epiCHAOS scores
het <- compute.eITH(datasets)

#--- add column "n" indicating number of cell types
het$n[het$state %in% c("Monocyte", "CD8T", "HSC", "B cell", "DC")] <- 1
het$n[grepl("mix2", het$state)] <- 2
het$n[grepl("mix3", het$state)] <- 3
het$n[grepl("mix4", het$state)] <- 4
het$n[grepl("mix-5", het$state)] <- 5


#--- Figure 1E. boxplot comparing epiCHAOS scores vs number of celltypes mixed
ggplot(unique(het[,c("mean.het", "state", "n")]), aes(y=mean.het, x=as.factor(n))) +
  geom_boxplot(fill="lightsteelblue")+
  geom_jitter(size=0.5)+
  labs( x="number of celltypes", y="epiCHAOS") +
  theme_classic()


#--- create a umap on the selected celltypes
atac <- atac[,atac$BioClassification %in% c("11_CD14.Mono.1", "17_B", "01_HSC", "24_CD8.CM", "09_pDC")]

#--- create a signac object of the selected cell types & regions
counts <-  atac@assays@.xData$data$counts
rownames(counts) <- names(atac@rowRanges) %>% str_replace_all("_", "-")
chrom_assay <- CreateChromatinAssay(
  counts = counts,
  sep = c(":", "-"),
  min.cells = 10,
  min.features = 200
)

obj <- CreateSeuratObject(
  counts = chrom_assay,
  assay = "peaks",
  meta.data = atac@colData
)


#--- run TF-IDF
obj <- RunTFIDF(obj)
obj <- FindTopFeatures(obj, min.cutoff = 'q0')
obj <- RunSVD(obj)
obj <- RunUMAP(object = obj, reduction = 'lsi', dims = 2:30)
obj <- FindNeighbors(object = obj, reduction = 'lsi', dims = 2:30)
obj <- FindClusters(object = obj, verbose = FALSE, algorithm = 3)
DimPlot(object = obj, label = TRUE) + NoLegend()

temp <- obj@reductions$umap@cell.embeddings %>% merge(atac@colData, by=0)
head(temp)

#--- Figure 1E. plot UMAP from subsetted data
ggplot(temp, aes(x = umap_1, y = umap_2, col=BioClassification)) +
  geom_point(size=0.5, alpha=0.7)+
  labs(x="UMAP1", y="UMAP2") +
  theme_classic() +
  scale_color_manual(values = c("thistle", "darkseagreen", "salmon1", "lightsteelblue", "lightgoldenrod2"))

#--- or with original umap coordinates
ggplot(temp, aes(x = UMAP1, y = UMAP1, col=BioClassification)) +
  geom_point(size=0.5, alpha=0.7)+
  labs(x="UMAP1", y="UMAP2") +
  theme_classic() +
  scale_color_manual(values = c("thistle", "darkseagreen", "salmon1", "lightsteelblue", "lightgoldenrod2"))

