

#--- use the scATAC-seq data from one of the analysed datasets - e.g. data from breast cancer - to explore how the epiCHAOS score is affected with different clustering resolution

#--- scATAC data from breast cancer samples from : https://pubmed.ncbi.nlm.nih.gov/35676392/

data.dir <- "/omics/groups/OE0219/internal/KatherineK/data/scATAC/Kumegawa_Brca"
analysis.dir <- "/omics/groups/OE0219/internal/KatherineK/ATACseq/breast-cancer"

setwd(analysis.dir)

#--- load ArchR packages
library(ArchR)
library(BSgenome.Hsapiens.UCSC.hg19)


addArchRThreads(threads = 16)

addArchRGenome("hg19")

brca <- loadArchRProject("epithelial/")

#--- plot clusters
p1 <- plotEmbedding(ArchRProj = brca, colorBy = "cellColData", name = "Clusters.epith", embedding = "UMAP.epith")
p1

#--- run clustering with different resolutions

# 0.1
brca <- addClusters(
  input = brca,
  reducedDims = "IterativeLSI",
  method = "Seurat",
  name = "Clusters.01",
  resolution = 0.1
)

#--- plot clusters
p1 <- plotEmbedding(ArchRProj = brca, colorBy = "cellColData", name = "Clusters.01", embedding = "UMAP.epith")
p1 <- p1 + theme_classic() + labs(x="UMAP1", y="UMAP2", title = "resolution: 0.1") + NoLegend()

# 0.2
brca <- addClusters(
  input = brca,
  reducedDims = "IterativeLSI",
  method = "Seurat",
  name = "Clusters.02",
  resolution = 0.2
)

#--- plot clusters
p2 <- plotEmbedding(ArchRProj = brca, colorBy = "cellColData", name = "Clusters.02", embedding = "UMAP.epith")
p2 <- p2 + theme_classic() + labs(x="UMAP1", y="UMAP2", title = "resolution: 0.2") + NoLegend()

# 0.3
brca <- addClusters(
  input = brca,
  reducedDims = "IterativeLSI",
  method = "Seurat",
  name = "Clusters.03",
  resolution = 0.3
)

#--- plot clusters
p3 <- plotEmbedding(ArchRProj = brca, colorBy = "cellColData", name = "Clusters.03", embedding = "UMAP.epith")
p3 <- p3 + theme_classic() + labs(x="UMAP1", y="UMAP2", title = "resolution: 0.3") + NoLegend()

# 0.4
brca <- addClusters(
  input = brca,
  reducedDims = "IterativeLSI",
  method = "Seurat",
  name = "Clusters.04",
  resolution = 0.4
)

#--- plot clusters
p4 <- plotEmbedding(ArchRProj = brca, colorBy = "cellColData", name = "Clusters.04", embedding = "UMAP.epith")
p4 <- p4 + theme_classic() + labs(x="UMAP1", y="UMAP2", title = "resolution: 0.4") + NoLegend()

# 0.5
brca <- addClusters(
  input = brca,
  reducedDims = "IterativeLSI",
  method = "Seurat",
  name = "Clusters.05",
  resolution = 0.5
)

#--- plot clusters
p5 <- plotEmbedding(ArchRProj = brca, colorBy = "cellColData", name = "Clusters.05", embedding = "UMAP.epith")
p5 <- p5 + theme_classic() + labs(x="UMAP1", y="UMAP2", title = "resolution: 0.5") + NoLegend()

# 0.6
brca <- addClusters(
  input = brca,
  reducedDims = "IterativeLSI",
  method = "Seurat",
  name = "Clusters.06",
  resolution = 0.6
)

#--- plot clusters
p6 <- plotEmbedding(ArchRProj = brca, colorBy = "cellColData", name = "Clusters.06", embedding = "UMAP.epith")
p6 <- p6 + theme_classic() + labs(x="UMAP1", y="UMAP2", title = "resolution: 0.6") + NoLegend()

# 0.7
brca <- addClusters(
  input = brca,
  reducedDims = "IterativeLSI",
  method = "Seurat",
  name = "Clusters.07",
  resolution = 0.7
)

#--- plot clusters
p7 <- plotEmbedding(ArchRProj = brca, colorBy = "cellColData", name = "Clusters.07", embedding = "UMAP.epith")
p7 <- p7 + theme_classic() + labs(x="UMAP1", y="UMAP2", title = "resolution: 0.7") + NoLegend()

# 0.8
brca <- addClusters(
  input = brca,
  reducedDims = "IterativeLSI",
  method = "Seurat",
  name = "Clusters.08",
  resolution = 0.8
)

#--- plot clusters
p8 <- plotEmbedding(ArchRProj = brca, colorBy = "cellColData", name = "Clusters.08", embedding = "UMAP.epith")
p8 <- p8 + theme_classic() + labs(x="UMAP1", y="UMAP2", title = "resolution: 0.8") + NoLegend()

# 0.9
brca <- addClusters(
  input = brca,
  reducedDims = "IterativeLSI",
  method = "Seurat",
  name = "Clusters.09",
  resolution = 0.9
)

#--- plot clusters
p9 <- plotEmbedding(ArchRProj = brca, colorBy = "cellColData", name = "Clusters.09", embedding = "UMAP.epith")
p9 <- p9 + theme_classic() + labs(x="UMAP1", y="UMAP2", title = "resolution: 0.9") + NoLegend()

pdf("epithelial/umaps_clustering_different_resolutions.pdf", 20, 5)
ggpubr::ggarrange(p1, p2, p3, p4, p5, p6, p7, p8, p9, legend = F, nrow = 1)
dev.off()


saveArchRProject(brca, "epithelial/")

#--- run epiCHAOS on the different resolutions
mat <- readRDS("epithelial/peaks_matrix.Rds")
counts <- mat@assays@data$PeakMatrix

res <- "Clusters.09"

clusters <- brca@cellColData[,res] %>% unique()


sample.rows <- sample(nrow(mat), 20000)
datasets <- list()
for (i in clusters) {
  for (n in 1:5) {
    ids <- brca$cellNames[brca@cellColData[,res]==i]
    if (length(ids)>100) { ids <- sample(ids, 100)}
    datasets[[paste0(i, "-", n)]] <- counts[sample.rows,ids]
  }
}

lapply(datasets, dim)

het.09 <- compute_eITH(datasets)
het.09

#--- save scores for each clustering resolution
saveRDS(het.01, file = "epithelial/epiCHAOS_clusters01.Rds")
saveRDS(het.02, file = "epithelial/epiCHAOS_clusters02.Rds")
saveRDS(het.03, file = "epithelial/epiCHAOS_clusters03.Rds")
saveRDS(het.04, file = "epithelial/epiCHAOS_clusters04.Rds")
saveRDS(het.05, file = "epithelial/epiCHAOS_clusters05.Rds")
saveRDS(het.06, file = "epithelial/epiCHAOS_clusters06.Rds")
saveRDS(het.07, file = "epithelial/epiCHAOS_clusters07.Rds")
saveRDS(het.08, file = "epithelial/epiCHAOS_clusters08.Rds")
saveRDS(het.09, file = "epithelial/epiCHAOS_clusters09.Rds")


#--- for each of het.01 to het.09, get the average epiCHAOS score per cluster, add it as a column to the project CellColData, and draw a UMAP colored by that

het.res <- list(het.01=het.01, het.02=het.02, het.03=het.03,
                het.04=het.04, het.05=het.05, het.06=het.06,
                het.07=het.07, het.08=het.08, het.09=het.09)

#--- make umap plots colored by epiCHAOS score for different clustering resolutions
gg <- list(p1,p2,p3,p4,p5,p6,p7,p8,p9)
for (i in names(het.res)) {
  print(i)
  
  # column name corresponing to clusters on that resolution
  res <- i %>% str_replace("het", "Clusters")
  
  # mean heterogeneity score from subsamples per cluster
  het <- het.res[[i]] %>% group_by(state) %>% summarise(mean=mean(het.adj))
  
  # add epiCHAOS column to project metadata
  for (cluster in het$state) {
    brca@cellColData$epiCHAOS[brca@cellColData[,res]==cluster] <- het$mean[het$state==cluster]
  }
  
  # create plot
  gg[[res]] <- plotEmbedding(ArchRProj = brca, colorBy = "cellColData", name = "epiCHAOS", embedding = "UMAP.epith") +
    theme_classic() + labs(x="UMAP1", y="UMAP2", title = i)
  
}

#--- print umaps to pds for clusters and epiCHAOS scores
pdf("epithelial/umaps_epiCHAOS_score_different_resolutions.pdf", 20, 5)
ggpubr::ggarrange(plotlist = gg[10:18], nrow = 1, ncol=9, common.legend = T)
dev.off()

#--- print umaps to pds for clusters and epiCHAOS scores
pdf("epithelial/umaps_clustering_and_epiCHAOS_score_different_resolutions.pdf", 20, 10)
ggpubr::ggarrange(plotlist = gg, nrow = 2, ncol=9, common.legend = T)
dev.off()
