
#--- scATAC data from breast cancer samples from : https://pubmed.ncbi.nlm.nih.gov/35676392/

data.dir <- "/omics/groups/OE0219/internal/KatherineK/data/scATAC/Kumegawa_Brca"
analysis.dir <- "/omics/groups/OE0219/internal/KatherineK/ATACseq/breast-cancer"

setwd(analysis.dir)

#--- load ArchR packages
library(ArchR)
library(BSgenome.Hsapiens.UCSC.hg19)


addArchRThreads(threads = 16)

addArchRGenome("hg19")

#--- path to fragments files
fragments <- paste0(data.dir, list.files(data.dir, pattern = "fragments.tsv.gz$"))
names(fragments) <- list.files(data.dir, pattern = "fragments.tsv.gz$") %>% str_split("\\.") %>% lapply("[", 1) %>% unlist()

#--- create arrow files
ArrowFiles <- createArrowFiles(
  inputFiles = fragments,
  sampleNames = names(fragments),
  minTSS = 4, 
  minFrags = 1000,
  addTileMat = TRUE,
  addGeneScoreMat = TRUE
)

ArrowFiles <- list.files(pattern = "arrow")

#--- doublet prediction
doubScores <- addDoubletScores(
  input = ArrowFiles,
  k = 10, #Refers to how many cells near a "pseudo-doublet" to count.
  knnMethod = "UMAP", #Refers to the embedding to use for nearest neighbor search with doublet projection.
  LSIMethod = 1
)

#--- create an ArchR project
brca <- ArchRProject(
  ArrowFiles = ArrowFiles,
  outputDirectory = analysis.dir,
  copyArrows = TRUE #This is recommened so that if you modify the Arrow files you have an original copy for later usage.
)

#--- save ArchR project
saveArchRProject(ArchRProj = brca, load = FALSE)

#--- filter doublets
brca <- filterDoublets(brca)

#--- add iterative LSI
brca <- addIterativeLSI(
  ArchRProj = brca,
  useMatrix = "TileMatrix",
  name = "IterativeLSI",
  iterations = 2,
  clusterParams = list( #See Seurat::FindClusters
    resolution = c(0.2),
    sampleCells = 10000,
    n.start = 10),
  varFeatures = 25000,
  dimsToUse = 1:30)


#--- add clusters
brca <- addClusters(
  input = brca,
  reducedDims = "IterativeLSI",
  method = "Seurat",
  name = "Clusters",
  resolution = 0.8
)

#--- add umap
brca <- addUMAP(
  ArchRProj = brca,
  reducedDims = "IterativeLSI",
  name = "UMAP",
  nNeighbors = 30,
  minDist = 0.5,
  metric = "cosine"
)

p1 <- plotEmbedding(ArchRProj = brca, colorBy = "cellColData", name = "Sample", embedding = "UMAP")
p2 <- plotEmbedding(ArchRProj = brca, colorBy = "cellColData", name = "Sample", embedding = "UMAP")

ggAlignPlots(p1, p2, type = "h")

plotPDF(p1,p2, name = "Plot-UMAP-Sample-Clusters.pdf", ArchRProj = brca, addDOC = FALSE, width = 5, height = 5)

saveArchRProject(ArchRProj = brca, load = FALSE)


#--- add group coverages for peak calling
brca <- addGroupCoverages(brca)

#--- calling peaks using macs2
pathToMacs2 <- findMacs2()

brca <- addReproduciblePeakSet(
  ArchRProj = brca,
  pathToMacs2 = pathToMacs2
)

saveArchRProject(ArchRProj = brca, load = FALSE)

#--- add a peaks matrix
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


#--- save gene score matrix unbinarised
mat <- getMatrixFromProject(
  ArchRProj = brca,
  useMatrix = "GeneScoreMatrix",
  useSeqnames = NULL,
  verbose = TRUE,
  binarize = F,
  threads = getArchRThreads(),
  logFile = createLogFile("getMatrixFromProject")
)

saveRDS(mat, "matrices/gene_score_matrix.Rds")


#--- identify which clusters are epithelial cells based on gene scores for EPCAM and KRT
mat.ge <- readRDS("matrices/gene_score_matrix.Rds")
temp <- mat.ge@assays@data$GeneScoreMatrix
rownames(temp) <- mat.ge@elementMetadata$name
colnames(temp) <- rownames(mat.ge@colData)
temp[1:10,1:10]

temp <- temp[c("EPCAM", "KRT1", "CD3D", "KDR", "ENG", "PECAM1", "PAX5", "ITGAM", "ITGAX", "CXCL12", "THY1", "CD38", "IGLL5"),] %>% t() %>% merge(brca@embeddings$UMAP$df, by=0) %>% keepRow() %>% merge(brca@cellColData[,c("Sample", "Clusters")], by=0) %>% data.frame()
temp[,2:14][temp[,2:14]>0] <- 1
colnames(temp)[grepl("UMAP", colnames(temp))] <- c("UMAP1", "UMAP2")

#--- clusters
p1 <- ggplot(temp, aes(y=UMAP1, x=UMAP2, color=Clusters)) +
  geom_point(size=0.5)+
  theme_classic()

p2 <- ggplot(temp, aes(y=UMAP1, x=UMAP2, color=Sample)) +
  geom_point(size=0.5)+
  theme_classic()

p3 <- ggplot(temp, aes(y=UMAP1, x=UMAP2, color=CD3D)) +
  geom_point(size=0.5)+
  labs(subtitle="T cells") +
  theme_classic()

p4 <- ggplot(temp, aes(y=UMAP1, x=UMAP2, color=EPCAM)) +
  geom_point(size=0.5)+
  labs(subtitle="Epithelial cells") +
  theme_classic()

p5 <- ggplot(temp, aes(y=UMAP1, x=UMAP2, color=KDR+ENG)) +
  geom_point(size=0.5)+
  labs(subtitle="Endothelial cells") +
  theme_classic()

p6 <- ggplot(temp, aes(y=UMAP1, x=UMAP2, color=PAX5)) +
  geom_point(size=0.5)+
  labs(subtitle="B cells") +
  theme_classic()

p7 <- ggplot(temp, aes(y=UMAP1, x=UMAP2, color=ITGAX+ITGAM)) +
  geom_point(size=0.5)+
  labs(subtitle="Myeloid cells") +
  theme_classic()

p8 <- ggplot(temp, aes(y=UMAP1, x=UMAP2, color=CD38+IGLL5)) +
  geom_point(size=0.5)+
  labs(subtitle="Plasma cells") +
  theme_classic()

p9 <- ggplot(temp, aes(y=UMAP1, x=UMAP2, color=CXCL12+THY1)) +
  geom_point(size=0.5)+
  labs(subtitle="Fibroblasts") +
  theme_classic()

pdf("umaps_celltype_annotating_by_marker_genes.pdf", 15,10)
ggarrange(p1,p2,p3,p4,p5,p6,p7,p8,p9, ncol=3, nrow=3)
dev.off()


#--- subset the object for epithelial cells (cluster 13:24) and perform clustering again on the subseted object
brca <- loadArchRProject(analysis.dir)
epith.cells <- brca$cellNames[brca$Clusters %in% paste0("C", 13:24)]

epith <- subsetCells(ArchRProj = brca, cellNames = epith.cells)

#--- add new iterative LSI
epith <- addIterativeLSI(
  ArchRProj = epith,
  useMatrix = "TileMatrix",
  name = "IterativeLSI.epith",
  iterations = 2,
  clusterParams = list( #See Seurat::FindClusters
    resolution = c(0.2),
    sampleCells = 10000,
    n.start = 10),
  varFeatures = 25000,
  dimsToUse = 1:30)


#--- add new clusters
epith <- addClusters(
  input = epith,
  reducedDims = "IterativeLSI.epith",
  method = "Seurat",
  name = "Clusters.epith",
  resolution = 0.8
)

#--- add new umap
epith <- addUMAP(
  ArchRProj = epith,
  reducedDims = "IterativeLSI.epith",
  name = "UMAP.epith",
  nNeighbors = 30,
  minDist = 0.5,
  metric = "cosine"
)

p1 <- plotEmbedding(ArchRProj = epith, colorBy = "cellColData", name = "Sample", embedding = "UMAP.epith")
p2 <- plotEmbedding(ArchRProj = epith, colorBy = "cellColData", name = "Clusters.epith", embedding = "UMAP.epith")

ggAlignPlots(p1, p2, type = "h")

plotPDF(p1,p2, name = "Plot-UMAP-Sample-Clusters_epithelial.pdf", ArchRProj = epith, addDOC = FALSE, width = 5, height = 5)

#--- load the subsetted object with only epithelial cells
brca <- loadArchRProject(file.path(analysis.dir,"epithelial/"))

#--- add group coverages for peak calling
brca <- addGroupCoverages(brca, force = T)

#--- calling peaks
pathToMacs2 <- findMacs2()

brca <- addReproduciblePeakSet(
  ArchRProj = brca,
  pathToMacs2 = pathToMacs2
)


#--- add a peaks matrix
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

saveRDS(mat, file.path(analysis.dir, "epithelial/peaks_matrix.Rds"))

#--- save gene score matrix unbinarised
mat <- getMatrixFromProject(
  ArchRProj = brca,
  useMatrix = "GeneScoreMatrix",
  useSeqnames = NULL,
  verbose = TRUE,
  binarize = F,
  threads = getArchRThreads(),
  logFile = createLogFile("getMatrixFromProject")
)

saveRDS(mat, file.path(analysis.dir, "epithelial/gene_score_matrix.Rds"))
saveArchRProject(ArchRProj = brca, load = FALSE)


#--- downstream analysis on epithelial subsetted data
setwd("epithelial/")

#--- epigenetic heterogeneity analysis
brca <- loadArchRProject()
mat <- readRDS("peaks_matrix.Rds")
dim(mat)

colnames(mat) <- paste0(brca$Sample, "_", brca$cellNames)
clusters <- unique(brca$Clusters.epith)
table(brca$Clusters.epith)

#--- list of peaks-by-cells matrices for each epithelial cluster
datasets <- list()
for (i in clusters) {
  ids <- colnames(mat)[brca$Clusters.epith==i]
  print(length(ids))
  datasets[[i]] <- mat[,ids]@assays@data$PeakMatrix
  if (ncol(datasets[[i]])>200) { datasets[[i]] <- datasets[[i]][,sample(ncol(datasets[[i]]), 200)]  }
}

lapply(datasets, dim)

#--- compute epiCHAOS scores
het <- compute_eITH(datasets)


#--- plot clusters and epiCHAOS scores in UMAP
temp <- brca@embeddings$UMAP.epith$df %>% merge(brca@cellColData[,c("Sample", "Clusters.epith")], by=0) %>% data.frame()
colnames(temp)[grepl("UMAP", colnames(temp))] <- c("UMAP1", "UMAP2")
temp <- merge(temp, het, by.x="Clusters.epith", by.y="state", all.x=T) %>% keepRow()
temp$Clusters <- temp$Clusters.epith

colors <- c("1"="#D51F26","2"="#272E6A","3"="#208A42","4"="#89288F","5"="#F47D2B", "6"="#FEE500","7"="#8A9FD1","8"="#C06CAB","19"="#E6C2DC",
            "10"="#90D5E4", "11"="#89C75F","12"="#F37B7D","13"="#9983BD","14"="#D24B27","15"="#3BBCA8", "16"="#6E4B9E","17"="#0C727C", "18"="#7E1416","9"="#D8A767","20"="#3D3D3D") %>% as.vector()

# clusters
p1 <- ggplot(temp, aes(x=UMAP1, y=UMAP2, color=Clusters)) +
  geom_point(size=0.1)+
  scale_color_manual(values = colors)+
  theme_classic()

p2 <- ggplot(temp, aes(x=UMAP1, y=UMAP2, color=het.adj)) +
  geom_point(size=0.1)+
  scale_color_distiller(palette = "Blues", direction = 1, na.value = "grey70")+
  theme_classic()

ggarrange(p1,p2,ncol=2, widths = c(1, 1))


#--- compute per-chromosome count-corrected epiCHAOS scores and perform subsampling of 5 x 100 cells per cluster
setwd("/omics/groups/OE0219/internal/KatherineK/ATACseq/breast-cancer/epithelial/")

#--- load data
brca <- loadArchRProject("/omics/groups/OE0219/internal/KatherineK/ATACseq/breast-cancer/epithelial")
mat <- readRDS("/omics/groups/OE0219/internal/KatherineK/ATACseq/breast-cancer/epithelial/peaks_matrix.Rds")
row.names <- mat@rowRanges %>% paste0()
mat <- mat@assays@data$PeakMatrix
rownames(mat) <- row.names
dim(mat)
mat[1:10,1:10]

#--- clusters
clusters <- brca$Clusters.epith %>% unique()

set.seed(10)

#--- create list of datasets
datasets <- list()
for (i in  clusters) {
  print(i)
  ids <- brca$cellNames[brca$Clusters.epith==i]
  if (length(ids)<30) { next }
  
  #--- number of cells for subsampling
  ncells <- min(c(100, length(ids)))
  
  #--- select 5 subsamples
  for (n in 1:5) {
    datasets[[paste0(i, "-", n)]] <- mat[,sample(ids, ncells)]
  }
}

lapply(datasets, dim)

#--- compute heterogeneity scores
het <- compute_eITH.cancer(datasets)

saveRDS(het, "epiCHAOS_scores_count_corrected_allpeaks_subsampling.Rds")


het$state <- het$state %>% str_split("-") %>% laply("[", 1) %>% unlist()

ggplot(het, aes(y=het, x=state)) +
  geom_boxplot()+
  theme_classic()
