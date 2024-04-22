

#--- (iii) testing eITH score on increasing numbers of cell types from Granja 2019 hematopoiesis

#--- scATAC from normal hematopoiesis from Granja et al. 2019

atac <- readRDS("/omics/groups/OE0219/internal/KatherineK/data/scATAC/Granja2019/scATAC-Healthy-Hematopoiesis-191120.rds")
atac <- atac[,grepl("BMMC", atac$Group)]

# select sites corresponding to hg19 promoters
require("TxDb.Hsapiens.UCSC.hg19.knownGene")
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
promoters <- promoters(genes(txdb), upstream = 1500, downstream = 500)
promoters.atac <- subsetByOverlaps(atac@rowRanges, promoters) %>% names()

#--- repeat to find marker peaks for all clusters (intensive, run as bsub)
markers <- readRDS("/omics/groups/OE0219/internal/KatherineK/ATACseq/cell-mixtures/hemato_celltype_markers.Rds")

unique(markers$cluster)

# select celltypes to include
select.celltypes <- c("11_CD14.Mono.1", "17_B", "01_HSC", "24_CD8.CM", "09_pDC")

# select top 1000 markers per selected cell type
top.markers <- list()
for (i in select.celltypes) {
  top <- markers[markers$cluster==i & abs(markers$avg_log2FC>2),] %>% arrange(p_val_adj) %>% head(500) %>% pull(gene) %>% str_replace_all("-", "_") 
  top.markers[[i]] <- gsub("\\..*","", top)
}

names(top.markers)
length(unlist(top.markers))
head(top.markers[[1]])

# create datasets of cell mixtures
datasets <- list()
cells.select <- c()

# individual cell types
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

# mix 2 celltypes
comparisons <- combn(c(1,2,3,4,5), 2)
for (i in 1:ncol(comparisons)) {
  hema <- atac@assays$data$counts[,atac$BioClassification %in% select.celltypes[comparisons[,i]]]
  hema <- hema[,colnames(hema) %in% cells.select]
  rownames(hema) <- names(atac@rowRanges)
  hema <- hema[rownames(hema) %in% unlist(top.markers),sample(ncol(hema), 100)]
  hema[hema>1] <- 1
  datasets[[paste0("mix2-", i)]] <- hema
}

# mix 3 celltypes
comparisons <- combn(c(1,2,3,4,5), 3)
for (i in 1:ncol(comparisons)) {
  hema <- atac@assays$data$counts[,atac$BioClassification %in% select.celltypes[comparisons[,i]]]
  hema <- hema[,colnames(hema) %in% cells.select]
  rownames(hema) <- names(atac@rowRanges)
  hema <- hema[rownames(hema) %in% unlist(top.markers),sample(ncol(hema), 100)]
  hema[hema>1] <- 1
  datasets[[paste0("mix3-", i)]] <- hema
}

# mix 4 celltypes
comparisons <- combn(c(1,2,3,4,5), 4)
for (i in 1:ncol(comparisons)) {
  hema <- atac@assays$data$counts[,atac$BioClassification %in% select.celltypes[comparisons[,i]]]
  hema <- hema[,colnames(hema) %in% cells.select]
  rownames(hema) <- names(atac@rowRanges)
  hema <- hema[rownames(hema) %in% unlist(top.markers),sample(ncol(hema), 100)]
  hema[hema>1] <- 1
  datasets[[paste0("mix4-", i)]] <- hema
}

# mix 5 celltypes
hema <- atac@assays$data$counts[,atac$BioClassification %in% select.celltypes[1:5]]
hema <- hema[,colnames(hema) %in% cells.select]
rownames(hema) <- names(atac@rowRanges)
hema <- hema[rownames(hema) %in% unlist(top.markers),sample(ncol(hema), 100)]
hema[hema>1] <- 1
datasets[["mix-5"]] <- hema

lapply(datasets, dim)
names(datasets)[1:5] <- c("Monocyte", "B cell", "HSC", "CD8T", "DC")

# compute eITH
het <- compute.eITH(datasets)

het$n[het$state %in% c("Monocyte", "CD8T", "HSC", "B cell", "DC")] <- 1
het$n[grepl("mix2", het$state)] <- 2
het$n[grepl("mix3", het$state)] <- 3
het$n[grepl("mix4", het$state)] <- 4
het$n[grepl("mix-5", het$state)] <- 5



# add scatterplot for number of celltypes mixed
het$color <- "Mixture"
het$color[het$state=="Monocyte"] <- "Monocyte"
het$color[het$state=="DC"] <- "DC"
het$color[het$state=="HSC"] <- "HSC"
het$color[het$state=="B cell"] <- "B cell"
het$color[het$state=="CD8T"] <- "CD8T"


dotplot <- ggplot(unique(het[,c("mean.het", "state", "n", "color")]), aes(y=mean.het, x=n, color=color)) +
  geom_point( size=2)+
  labs( x="number of celltypes", y="epiCHAOS") +
  scale_color_manual(values = c( "#1B9E77",  "#D95F02",  "#7570B3", "#E7298A", "grey20", "#66A61E"))+
  #stat_cor(method = "spearman")+
  #geom_smooth(method='lm', color="grey20")+
  lims(y=c(0, 1), x=c(0.9, 5.1))+
  theme_classic()

boxplot <- ggplot(unique(het[,c("mean.het", "state", "n", "color")]), aes(y=mean.het, x=as.factor(n))) +
  geom_boxplot(fill="lightsteelblue")+
  geom_jitter(size=0.5)+
  labs( x="number of celltypes", y="epiCHAOS") +
  theme_classic()

pdf("Figure 1/boxplot_epiCHAOS_celltype_mixtures.pdf", 4, 3)
boxplot
dev.off()

svg("Figure 1/boxplot_epiCHAOS_celltype_mixtures.svg", 4, 3)
boxplot
dev.off()

# plot heatmap of selected celltype markers
ca <- data.frame(row.names = cells.select, celltype=rep(c("Monocytes", "B cells", "HSCs", "CD8T", "Erythroid"), each=100))
temp <- atac@assays$data$counts[rownames(atac) %in% unlist(top.markers),atac$BioClassification %in% select.celltypes]
temp[temp>1] <- 1

heatmap <- pheatmap::pheatmap(temp[,cells.select], show_rownames = F, show_colnames = F, treeheight_row = 0, color = c("lightsteelblue", "salmon"), legend=F, cutree_cols = 5, annotation_col = ca)

pdf("Figure 1/correlation_epiCHAOS_vs_ncelltypes.pdf", 4, 2.5)
dotplot
dev.off()

pdf("Figure 1/heatmap_marker_peaks_per_celltype.pdf", 5, 5)
heatmap
dev.off()


#--- create a umap on the selected celltypes
atac <- atac[,atac$BioClassification %in% c("11_CD14.Mono.1", "17_B", "01_HSC", "24_CD8.CM", "09_pDC")]
atac$BioClassification
atac@metadata

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

obj <- RunTFIDF(obj)
obj <- FindTopFeatures(obj, min.cutoff = 'q0')
obj <- RunSVD(obj)
obj <- RunUMAP(object = obj, reduction = 'lsi', dims = 2:30)
obj <- FindNeighbors(object = obj, reduction = 'lsi', dims = 2:30)
obj <- FindClusters(object = obj, verbose = FALSE, algorithm = 3)
DimPlot(object = obj, label = TRUE) + NoLegend()

temp <- obj@reductions$umap@cell.embeddings %>% merge(atac@colData, by=0)
head(temp)

pdf("Figure 1/umap_selected_celltypes_to_mix.pdf", 5, 3)
ggplot(temp, aes(x = UMAP1, y = UMAP2, col=BioClassification)) +
  geom_point(size=0.5, alpha=0.7)+
  theme_classic() +
  scale_color_manual(values = c("thistle", "darkseagreen", "salmon1", "lightsteelblue", "lightgoldenrod2"))
dev.off()

svg("Figure 1/umap_selected_celltypes_to_mix.svg", 5, 3)
ggplot(temp, aes(x = UMAP1, y = UMAP2, col=BioClassification)) +
  geom_point(size=0.5, alpha=0.7)+
  theme_classic() +
  scale_color_manual(values = c("thistle", "darkseagreen", "salmon1", "lightsteelblue", "lightgoldenrod2"))
dev.off()