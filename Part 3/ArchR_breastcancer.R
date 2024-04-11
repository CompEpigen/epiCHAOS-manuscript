
#--- scATAC data from breast cancer samples from : https://pubmed.ncbi.nlm.nih.gov/35676392/

setwd("/omics/groups/OE0219/internal/KatherineK/ATACseq/breast-cancer/epithelial/")

library(ArchR)
library(BSgenome.Hsapiens.UCSC.hg19)

# 
# addArchRThreads(threads = 16)
# 
# addArchRGenome("hg19")
# fragments <- paste0("/omics/groups/OE0219/internal/KatherineK/data/scATAC/Kumegawa_Brca/",
#                     list.files("/omics/groups/OE0219/internal/KatherineK/data/scATAC/Kumegawa_Brca/", pattern = "fragments.tsv.gz$"))
# names(fragments) <- list.files("/omics/groups/OE0219/internal/KatherineK/data/scATAC/Kumegawa_Brca/", pattern = "fragments.tsv.gz$") %>% str_split("\\.") %>% lapply("[", 1) %>% unlist()
# 
# # create arrow files
# ArrowFiles <- createArrowFiles(
#   inputFiles = fragments,
#   sampleNames = names(fragments),
#   minTSS = 4, #Dont set this too high because you can always increase later
#   minFrags = 1000,
#   addTileMat = TRUE,
#   addGeneScoreMat = TRUE
# )
# 
# ArrowFiles <- list.files(pattern = "arrow")
# 
# # doublet prediction
# doubScores <- addDoubletScores(
#   input = ArrowFiles,
#   k = 10, #Refers to how many cells near a "pseudo-doublet" to count.
#   knnMethod = "UMAP", #Refers to the embedding to use for nearest neighbor search with doublet projection.
#   LSIMethod = 1
# )
# 
# # create an ArchR project
# brca <- ArchRProject(
#   ArrowFiles = ArrowFiles,
#   outputDirectory = "/omics/groups/OE0219/internal/KatherineK/ATACseq/breast-cancer/",
#   copyArrows = TRUE #This is recommened so that if you modify the Arrow files you have an original copy for later usage.
# )
# 
# # save ArchR project
# saveArchRProject(ArchRProj = brca, load = FALSE)
# 
# brca <- readRDS("/omics/groups/OE0219/internal/KatherineK/ATACseq/breast-cancer/Save-ArchR-Project.rds")
# 
# # filter doublets
# brca <- filterDoublets(brca)
# 
# brca <- addIterativeLSI(
#   ArchRProj = brca,
#   useMatrix = "TileMatrix",
#   name = "IterativeLSI",
#   iterations = 2,
#   clusterParams = list( #See Seurat::FindClusters
#     resolution = c(0.2),
#     sampleCells = 10000,
#     n.start = 10),
#   varFeatures = 25000,
#   dimsToUse = 1:30)
# 
# 
# # add clusters
# brca <- addClusters(
#   input = brca,
#   reducedDims = "IterativeLSI",
#   method = "Seurat",
#   name = "Clusters",
#   resolution = 0.8
# )
# 
# # umap
# brca <- addUMAP(
#   ArchRProj = brca,
#   reducedDims = "IterativeLSI",
#   name = "UMAP",
#   nNeighbors = 30,
#   minDist = 0.5,
#   metric = "cosine"
# )
# 
# p1 <- plotEmbedding(ArchRProj = brca, colorBy = "cellColData", name = "Sample", embedding = "UMAP")
# p2 <- plotEmbedding(ArchRProj = brca, colorBy = "cellColData", name = "Sample", embedding = "UMAP")
# 
# ggAlignPlots(p1, p2, type = "h")
# 
# plotPDF(p1,p2, name = "Plot-UMAP-Sample-Clusters.pdf", ArchRProj = brca, addDOC = FALSE, width = 5, height = 5)
# 
# saveArchRProject(ArchRProj = brca, load = FALSE)
# 

# 
# brca <- addGroupCoverages(brca)
# 
# # calling peaks
# pathToMacs2 <- findMacs2()
# 
# brca <- addReproduciblePeakSet(
#   ArchRProj = brca, 
#   pathToMacs2 = pathToMacs2
# )
# 
# saveArchRProject(ArchRProj = brca, load = FALSE)
# 
# # add a peaks matrix
# brca <- addPeakMatrix(brca)
# 
# mat <- getMatrixFromProject(
#   ArchRProj = brca,
#   useMatrix = "PeakMatrix",
#   useSeqnames = NULL,
#   verbose = TRUE,
#   binarize = T,
#   threads = getArchRThreads(),
#   logFile = createLogFile("getMatrixFromProject")
# )
# 

# # save gene score matrix unbinarised
# mat <- getMatrixFromProject(
#   ArchRProj = brca,
#   useMatrix = "GeneScoreMatrix",
#   useSeqnames = NULL,
#   verbose = TRUE,
#   binarize = F,
#   threads = getArchRThreads(),
#   logFile = createLogFile("getMatrixFromProject")
# )
# 
# saveRDS(mat, "matrices/gene_score_matrix.Rds")


#--- identify which clusters are epithelial cells based on gene scores for EPCAM and KRT
# mat.ge <- readRDS("matrices/gene_score_matrix.Rds")
# temp <- mat.ge@assays@data$GeneScoreMatrix
# rownames(temp) <- mat.ge@elementMetadata$name
# colnames(temp) <- rownames(mat.ge@colData)
# temp[1:10,1:10]
# 
# temp <- temp[c("EPCAM", "KRT1", "CD3D", "KDR", "ENG", "PECAM1", "PAX5", "ITGAM", "ITGAX", "CXCL12", "THY1", "CD38", "IGLL5"),] %>% t() %>% merge(brca@embeddings$UMAP$df, by=0) %>% keepRow() %>% merge(brca@cellColData[,c("Sample", "Clusters")], by=0) %>% data.frame()
# temp[,2:14][temp[,2:14]>0] <- 1
# colnames(temp)[grepl("UMAP", colnames(temp))] <- c("UMAP1", "UMAP2")
# 
# # clusters
# p1 <- ggplot(temp, aes(y=UMAP1, x=UMAP2, color=Clusters)) +
#   geom_point(size=0.5)+
#   theme_classic()
# 
# p2 <- ggplot(temp, aes(y=UMAP1, x=UMAP2, color=Sample)) +
#   geom_point(size=0.5)+
#   theme_classic()
# 
# p3 <- ggplot(temp, aes(y=UMAP1, x=UMAP2, color=CD3D)) +
#   geom_point(size=0.5)+
#   labs(subtitle="T cells") +
#   theme_classic()
# 
# p4 <- ggplot(temp, aes(y=UMAP1, x=UMAP2, color=EPCAM)) +
#   geom_point(size=0.5)+
#   labs(subtitle="Epithelial cells") +
#   theme_classic()
# 
# p5 <- ggplot(temp, aes(y=UMAP1, x=UMAP2, color=KDR+ENG)) +
#   geom_point(size=0.5)+
#   labs(subtitle="Endothelial cells") +
#   theme_classic()
# 
# p6 <- ggplot(temp, aes(y=UMAP1, x=UMAP2, color=PAX5)) +
#   geom_point(size=0.5)+
#   labs(subtitle="B cells") +
#   theme_classic()
# 
# p7 <- ggplot(temp, aes(y=UMAP1, x=UMAP2, color=ITGAX+ITGAM)) +
#   geom_point(size=0.5)+
#   labs(subtitle="Myeloid cells") +
#   theme_classic()
# 
# p8 <- ggplot(temp, aes(y=UMAP1, x=UMAP2, color=CD38+IGLL5)) +
#   geom_point(size=0.5)+
#   labs(subtitle="Plasma cells") +
#   theme_classic()
# 
# p9 <- ggplot(temp, aes(y=UMAP1, x=UMAP2, color=CXCL12+THY1)) +
#   geom_point(size=0.5)+
#   labs(subtitle="Fibroblasts") +
#   theme_classic()
# 
# pdf("umaps_celltype_annotating_by_marker_genes.pdf", 15,10)
# ggarrange(p1,p2,p3,p4,p5,p6,p7,p8,p9, ncol=3, nrow=3)
# dev.off()
# 
# plotEmbedding(brca, colorBy = "geneScoreMatrix", name = "CXCL12")
# 
# 
#--- epigenetic heterogeneity analysis
brca <- loadArchRProject()
mat <- readRDS("peak_matrix.Rds")
promoters <- brca@peakSet$idx[brca@peakSet$peakType=="Distal"]
mat <- mat[promoters,]
dim(mat)

colnames(mat) <- paste0(brca$Sample, "_", brca$cellNames)
#clusters <- colnames(mat) %>% str_split("_") %>% lapply("[", 2) %>% unlist() %>% unique() %>% paste0("_")
clusters <- unique(brca$Clusters.epith)
table(brca$Clusters.epith)

# clusters 14:24 are likely epithelial cell clusters abd likely containing the malignant cells
#clusters <- paste0("C", 13:24)

# after subsetting for epitehlial cells
#mat <- mat[,brca$cellNames]

datasets <- list()
for (i in clusters) {
  #ids <- colnames(mat)[grepl(i, colnames(mat))]
  ids <- colnames(mat)[brca$Clusters.epith==i]
  print(length(ids))
  datasets[[i]] <- mat[,ids]@assays@data$PeakMatrix
  if (ncol(datasets[[i]])>200) { datasets[[i]] <- datasets[[i]][,sample(ncol(datasets[[i]]), 200)]  }
}

lapply(datasets, dim)

# compute heterogeneity scores
het <- compute.eITH(datasets)

saveRDS(het, "heterogeneity_scores_per_cluster_epithelial_after_reclustering.Rds")

# #--- plot clusters and epiCHAOS scores in UMAP
# temp <- brca@embeddings$UMAP$df %>% merge(brca@cellColData[,c("Sample", "Clusters")], by=0) %>% data.frame()
# colnames(temp)[grepl("UMAP", colnames(temp))] <- c("UMAP1", "UMAP2")
# temp <- merge(temp, het, by.x="Clusters", by.y="state", all.x=T) %>% keepRow()
# temp$Clusters[temp$Clusters %in% het$state] <- "Epithelial cells"
# temp$Clusters[temp$Clusters %in% c("C1", "C4", "C5", "C6", "C7", "C8", "C9")] <- "Lymphoid cells"
# temp$Clusters[temp$Clusters=="C12"] <- "Plasma cells"
# temp$Clusters[temp$Clusters=="C2"] <- "Endothelial cells"
# temp$Clusters[temp$Clusters=="C10"] <- "Myeloid cells"
# temp$Clusters[temp$Clusters=="C11"] <- "B cells"
# temp$Clusters[temp$Clusters=="C3"] <- "Unknown"
# 
# # clusters
# p1 <- ggplot(temp, aes(x=UMAP1, y=UMAP2, color=Clusters)) +
#   geom_point(size=0.1)+
#   scale_color_brewer(palette = "Set3", direction = -1)+
#   theme_classic()
# 
# p2 <- ggplot(temp, aes(x=UMAP1, y=UMAP2, color=mean.het)) +
#   geom_point(size=0.1)+
#   scale_color_distiller(palette = "Blues", direction = 1, na.value = "grey70")+
#   theme_classic()
#   
# pdf("umaps_clusters_epiCHAOS_brca.pdf", 8, 2.5)
# ggarrange(p1,p2,ncol=2, widths = c(1, 0.8))
# dev.off()
# 
# 
# # #--- correlate with gene score matrix
# # mat <- getMatrixFromProject(
# #   ArchRProj = brca,
# #   useMatrix  = "GeneScoreMatrix",
# #   useSeqnames = NULL,
# #   verbose = TRUE,
# #   binarize = T,
# #   threads = getArchRThreads(),
# #   logFile = createLogFile("getMatrixFromProject")
# # )
# # 
# # saveRDS(mat, "matrices/gene_score_matrix.Rds")
# 
# #--- correlate per-cluster eICH with gene activity matrix
# mat.ge <- readRDS("matrices/gene_score_matrix.Rds")
# temp <- mat.ge@assays@data$GeneScoreMatrix
# rownames(temp) <- mat.ge@elementMetadata$name
# colnames(temp) <- rownames(mat.ge@colData)
# temp <- t(temp)         
# temp[1:10,1:10]
# 
# gex.per.cluster <- data.frame()
# for (i in unique(clusters)) {
#   print(i)
#   ids <- brca$cellNames[brca$Clusters==i]
#   gex.per.cluster <- rbind(gex.per.cluster, colMeans(temp[ids,]))
# }
# 
# rownames(gex.per.cluster) <- unique(clusters)
# colnames(gex.per.cluster) <- colnames(temp)
# gex.per.cluster[1:10,1:10]
# gex.per.cluster <- merge(gex.per.cluster , het, by.x=0, by.y="state") %>% keepRow()
# 
# cor <- cor(gex.per.cluster, gex.per.cluster$mean.het) %>% data.frame()
# cor$gene <- rownames(cor)
# cor <- cor[!is.na(cor$.),]
# 
# cor[cor$.>0.8,]
# cor[cor$.<(-0.80),]
# 
# #--- test gene ontology enrichment of the most highly correlated genes
# library(clusterProfiler)
# #top.genes <- cor[cor$.<(-0.7),] %>% rownames()
# top.genes <- cor[cor$.>0.8,] %>% rownames()
# enriched <- enrichGO(top.genes, OrgDb="org.Hs.eg.db", keyType = "SYMBOL", ont = "BP", pvalueCutoff = 0.05, pAdjustMethod = "BH", universe=rownames(mat.ge), qvalueCutoff = 0.1)
# enriched@result %>% filter(qvalue<0.1)
# 
# 
# 
# # plot some of interest
# p1 <- ggplot(gex.per.cluster, aes(y=mean.het, x=TEAD1)) +
#   geom_point(size=3, alpha=0.7, color="steelblue4")+
#   labs(x="TEAD1", y="epiCHAOS") +
#   stat_cor()+
#   theme_bw()
# 
# p2 <- ggplot(gex.per.cluster, aes(y=mean.het, x=CTNNB1)) +
#   geom_point(size=3, alpha=0.7, color="steelblue4")+
#   labs(x="CTNNB1", y="epiCHAOS") +
#   stat_cor()+
#   theme_bw()
# 
# ggarrange(p1,p2)



# #--- subset the object for epithelial cells and perform clustering again on the subseted object
# brca <- readRDS("/omics/groups/OE0219/internal/KatherineK/ATACseq/breast-cancer/Save-ArchR-Project.rds")
# epith.cells <- brca$cellNames[brca$Clusters %in% paste0("C", 13:24)]
# 
# epith <- subsetCells(ArchRProj = brca, cellNames = epith.cells)
# 
# epith <- addIterativeLSI(
#   ArchRProj = epith,
#   useMatrix = "TileMatrix",
#   name = "IterativeLSI.epith",
#   iterations = 2,
#   clusterParams = list( #See Seurat::FindClusters
#     resolution = c(0.2),
#     sampleCells = 10000,
#     n.start = 10),
#   varFeatures = 25000,
#   dimsToUse = 1:30)
# 
# 
# # add clusters
# epith <- addClusters(
#   input = epith,
#   reducedDims = "IterativeLSI.epith",
#   method = "Seurat",
#   name = "Clusters.epith",
#   resolution = 0.8
# )
# 
# # umap
# epith <- addUMAP(
#   ArchRProj = epith,
#   reducedDims = "IterativeLSI.epith",
#   name = "UMAP.epith",
#   nNeighbors = 30,
#   minDist = 0.5,
#   metric = "cosine"
# )
# 
# p1 <- plotEmbedding(ArchRProj = epith, colorBy = "cellColData", name = "Sample", embedding = "UMAP.epith")
# p2 <- plotEmbedding(ArchRProj = epith, colorBy = "cellColData", name = "Clusters.epith", embedding = "UMAP.epith")
# 
# ggAlignPlots(p1, p2, type = "h")
# 
# plotPDF(p1,p2, name = "Plot-UMAP-Sample-Clusters_epithelial.pdf", ArchRProj = epith, addDOC = FALSE, width = 5, height = 5)

brca <- loadArchRProject("/omics/groups/OE0219/internal/KatherineK/ATACseq/breast-cancer/")

brca <- addGroupCoverages(brca, force = T)

# calling peaks
pathToMacs2 <- findMacs2()

brca <- addReproduciblePeakSet(
  ArchRProj = brca,
  pathToMacs2 = pathToMacs2
)


# add a peaks matrix
brca <- addPeakMatrix(brca)

mat <- getMatrixFromProject(
  ArchRProj = brca,
  useMatrix = "PeakMatrix",
  useSeqnames = NULL,
  verbose = TRUE,
  binarize = T,
  threads = getArchRThreads(),
  logFile = createLogFile("getMatrixFromProject")
)

saveRDS(mat, "/omics/groups/OE0219/internal/KatherineK/ATACseq/breast-cancer/epithelial/peaks_matrix.Rds")

# save gene score matrix unbinarised
mat <- getMatrixFromProject(
  ArchRProj = brca,
  useMatrix = "GeneScoreMatrix",
  useSeqnames = NULL,
  verbose = TRUE,
  binarize = F,
  threads = getArchRThreads(),
  logFile = createLogFile("getMatrixFromProject")
)

saveRDS(mat, "/omics/groups/OE0219/internal/KatherineK/ATACseq/breast-cancer/epithelial/gene_score_matrix.Rds")
saveArchRProject(ArchRProj = brca, load = FALSE)

# 
# 
# #--- downstream analysis on epithelial subsetted data
# 
# #--- plot clusters and epiCHAOS scores in UMAP
# temp <- brca@embeddings$UMAP.epith$df %>% merge(brca@cellColData[,c("Sample", "Clusters.epith")], by=0) %>% data.frame()
# colnames(temp)[grepl("UMAP", colnames(temp))] <- c("UMAP1", "UMAP2")
# temp <- merge(temp, het, by.x="Clusters.epith", by.y="state", all.x=T) %>% keepRow()
# temp$Clusters <- temp$Clusters.epith
# 
# colors <- c("1"="#D51F26","2"="#272E6A","3"="#208A42","4"="#89288F","5"="#F47D2B", "6"="#FEE500","7"="#8A9FD1","8"="#C06CAB","19"="#E6C2DC",
#             "10"="#90D5E4", "11"="#89C75F","12"="#F37B7D","13"="#9983BD","14"="#D24B27","15"="#3BBCA8", "16"="#6E4B9E","17"="#0C727C", "18"="#7E1416","9"="#D8A767","20"="#3D3D3D") %>% as.vector()
# 
# # clusters
# p1 <- ggplot(temp, aes(x=UMAP1, y=UMAP2, color=Clusters)) +
#   geom_point(size=0.1)+
#   scale_color_manual(values = colors)+
#   theme_classic()
# 
# p2 <- ggplot(temp, aes(x=UMAP1, y=UMAP2, color=mean.het)) +
#   geom_point(size=0.1)+
#   scale_color_distiller(palette = "Blues", direction = 1, na.value = "grey70")+
#   theme_classic()
# 
# pdf("umaps_clusters_epiCHAOS_brca_epithelial.pdf", 8, 2.5)
# ggarrange(p1,p2,ncol=2, widths = c(1, 1))
# dev.off()
# 
# #--- correlate per-cluster eICH with gene activity matrix
# mat.ge <- readRDS("matrices/gene_score_matrix.Rds")
# temp <- mat.ge@assays@data$GeneScoreMatrix
# rownames(temp) <- mat.ge@elementMetadata$name
# colnames(temp) <- rownames(mat.ge@colData)
# temp <- temp[,brca$cellNames]
# temp <- t(temp)
# temp[1:10,1:10]
# 
# gex.per.cluster <- data.frame()
# for (i in unique(clusters)) {
#   print(i)
#   ids <- brca$cellNames[brca$Clusters.epith==i]
#   gex.per.cluster <- rbind(gex.per.cluster, colMeans(temp[ids,]))
# }
# 
# rownames(gex.per.cluster) <- unique(clusters)
# colnames(gex.per.cluster) <- colnames(temp)
# gex.per.cluster[1:10,1:10]
# gex.per.cluster <- merge(gex.per.cluster , het, by.x=0, by.y="state") %>% keepRow()
# 
# cor <- cor(gex.per.cluster, gex.per.cluster$mean.het) %>% data.frame()
# cor$gene <- rownames(cor)
# cor <- cor[!is.na(cor$.),]
# 
# cor[cor$.>0.75,]
# cor[cor$.<(-0.70),]

#--- compare epiCHAOS in peaks vs windows matrix
mat <- getMatrixFromProject(
  ArchRProj = brca,
  useMatrix = "TileMatrix",
  useSeqnames = NULL,
  verbose = TRUE,
  binarize = T,
  threads = getArchRThreads(),
  logFile = createLogFile("getMatrixFromProject")
)

saveRDS(mat, file = "tile_matrix.Rds")

#--- compute epiCHAOS and compare to peaks
mat <- readRDS("tile_matrix.Rds")
mat <- mat@assays@data$TileMatrix
mat <- mat[rowSums(mat)>5,]
dim(mat)
mat[mat>1] <- 1

sample.sites <- sample(nrow(mat), 100000)
mat <- mat[sample.sites,]
dim(mat)

colnames(mat) <- paste0(brca$Sample, "_", brca$cellNames)
clusters <- unique(brca$Clusters.epith)

datasets <- list()
for (i in clusters) {
  ids <- colnames(mat)[brca$Clusters.epith==i]
  print(length(ids))
  datasets[[i]] <- mat[,ids]
  if (ncol(datasets[[i]])>200) { datasets[[i]] <- datasets[[i]][,sample(ncol(datasets[[i]]), 200)]  }
}

lapply(datasets, dim)

# compute heterogeneity scores
het2 <- compute.eITH(datasets)
saveRDS(het2, "heterogeneity_scores_tilematrix.Rds")

temp <- data.frame(TileMatrix=het2$mean.het, PeakMatrix=het$mean.het)

pdf("correlate_epiCHAOS_tiles_vs_peaks_brca.pdf", 5, 4)
ggplot(temp, aes(y=TileMatrix, x=PeakMatrix)) +
  geom_point(size=3, alpha=0.7,  color="steelblue4")+
  stat_cor() +
  labs(x="epiCHAOS PeakMatrix", y="epiCHAOS TileMatrix") +
  theme_bw()
dev.off()


#-- umap with color for heterogeneity scores
for (i in het$state) {
  brca$epiCHAOS[brca$Clusters.epith==i] <- het$mean.het[het$state==i]
}

p1 <- plotEmbedding(brca, embedding = "UMAP.epith", name = "Clusters.epith")
p2 <- plotEmbedding(brca, embedding = "UMAP.epith", name = "epiCHAOS", continuousSet = "whiteBlue")

plotPDF(p1,p2, name = "UMAP_epiCHAOS_epithelial.pdf", ArchRProj = brca, addDOC = FALSE, width = 5, height = 5)


pdf("barplot_epiCHAOS_percluster.pdf", 5, 2)
ggplot(het, aes(x = reorder(state, mean.het), y = mean.het)) +
  geom_bar(stat="identity", position = "dodge", alpha=0.8, width = 0.6)+
  labs(x="", y="epiCHAOS")+
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()




#--- use ArchR to find enrichments of TFBS activity per cluster

#--- find marker peaks
markersPeaks <- getMarkerFeatures(
  ArchRProj = brca,
  useMatrix = "PeakMatrix",
  groupBy = "Clusters.epith",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)

# add peak annotations for encode TFBS
brca <- addArchRAnnotations(ArchRProj = brca, collection = "EncodeTFBS")

# find enrichments
enrichEncode <- peakAnnoEnrichment(
  seMarker = markersPeaks,
  ArchRProj = brca,
  peakAnnotation = "EncodeTFBS",
  cutOff = "FDR <= 0.1 & Log2FC >= 0.5"
)

heatmapEncode <- plotEnrichHeatmap(enrichEncode, n = 7, transpose = TRUE)
ComplexHeatmap::draw(heatmapEncode, heatmap_legend_side = "bot", annotation_legend_side = "bot")

# add codex annotations
brca <- addArchRAnnotations(ArchRProj = brca, collection = "Codex")

# find enrichments
enrichCodex <- peakAnnoEnrichment(
  seMarker = markersPeaks,
  ArchRProj = brca,
  peakAnnotation = "Codex",
  cutOff = "FDR <= 0.1 & Log2FC >= 0.5"
)

heatmapCodex <- plotEnrichHeatmap(enrichCodex, n = 7, transpose = TRUE)
ComplexHeatmap::draw(heatmapCodex, heatmap_legend_side = "bot", annotation_legend_side = "bot")


#--- motif enrichment in differential peaks
brca <- addMotifAnnotations(ArchRProj = brca, motifSet = "cisbp", name = "Motif")
brca$Clusters2 <- ifelse(brca$Clusters.epith=="C19", "C19", "Other")

# find differential markers between C19 (the highest epiCHAOS cluster) and all other clusters
markerTest <- getMarkerFeatures(
  ArchRProj = brca, 
  useMatrix = "PeakMatrix",
  groupBy = "Clusters2",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "C19",
  bgdGroups = "Other"
)

motifsUp <- peakAnnoEnrichment(
  seMarker = markerTest,
  ArchRProj = brca,
  peakAnnotation = "Motif",
  cutOff = "FDR <= 0.1 & Log2FC >= 0.5"
)

df <- data.frame(TF = rownames(motifsUp), mlog10Padj = assay(motifsUp)[,1])
df <- df[order(df$mlog10Padj, decreasing = TRUE),]
df$rank <- seq_len(nrow(df))
head(df)

ggUp <- ggplot(df, aes(rank, mlog10Padj, color = mlog10Padj)) + 
  geom_point(size = 1) +
  ggrepel::geom_label_repel(
    data = df[rev(seq_len(30)), ], aes(x = rank, y = mlog10Padj, label = TF), 
    size = 1.5,
    nudge_x = 2,
    color = "black"
  ) + theme_ArchR() + 
  ylab("-log10(P-adj) Motif Enrichment") + 
  xlab("Rank Sorted TFs Enriched") +
  scale_color_gradientn(colors = paletteContinuous(set = "comet"))

ggUp

plotPDF(ggUp, name = "C19_up-Markers-Motifs-Enriched", width = 5, height = 5, ArchRProj = brca, addDOC = FALSE)


#--- motif deviations
if("Motif" %ni% names(brca@peakAnnotation)){
  brca <- addMotifAnnotations(ArchRProj = brca, motifSet = "cisbp", name = "Motif")
}

brca <- addBgdPeaks(brca)

brca <- addDeviationsMatrix(
  ArchRProj = brca, 
  peakAnnotation = "Motif",
  force = TRUE
)

# access the deviations
plotVarDev <- getVarDeviations(brca, name = "MotifMatrix", plot = TRUE)

saveArchRProject(brca, load=F)


#--- try computing epiCHAOS on all PRC2 targetted regions and compare to the global

mat <- readRDS("/omics/groups/OE0219/internal/KatherineK/ATACseq/breast-cancer/epithelial/peaks_matrix.Rds")

#--- encode TFBS
lola <- get(load("/omics/groups/OE0219/internal/KatherineK/Genomic_Annotations/LOLA/hg19/encode_tfbs/encode_tfbs.RData"))
index <- read.table("/omics/groups/OE0219/internal/KatherineK/Genomic_Annotations/LOLA/hg19/encode_tfbs/index.txt", header = T)
names(lola) <- index$filename
ezh2 <- lola$wgEncodeAwgTfbsBroadH1hescEzh239875UniPk.narrowPeak
names(mat@rowRanges) <- 1:length(mat@rowRanges)
ezh2 <- subsetByOverlaps(mat@rowRanges, ezh2) %>% names() %>% as.numeric()

row.names <- mat@rowRanges %>% paste0()
mat <- mat@assays@data$PeakMatrix
rownames(mat) <- row.names
mat <- mat[ezh2,]
dim(mat)
mat[1:10,1:10]


clusters <- brca$Clusters.epith %>% unique()

datasets <- list()
for (i in  clusters) {
  print(i)
  ids <- brca$cellNames[brca$Clusters.epith==i]
  if (length(ids)>100) {
    datasets[[i]] <- mat[,sample(ids, 100)]
  } else {
    datasets[[i]] <- mat[,ids]
  }
}

lapply(datasets, dim)

source("/omics/groups/OE0219/internal/KatherineK/ATACseq/eITH-test-scripts/jaccard.R")

het <- compute.eITH(datasets)
