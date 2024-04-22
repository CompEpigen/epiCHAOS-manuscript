


#--- liver cancer analysis with excluding non-malignant cells

#--- scATAC data from liver cancer: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE227265

setwd("/omics/groups/OE0219/internal/KatherineK/ATACseq/Liver-Cancer/malignant-subset/")


library(ArchR)
library(BSgenome.Hsapiens.UCSC.hg38)

addArchRThreads(threads = 16)

addArchRGenome("hg38")
fragments <- paste0("/omics/groups/OE0219/internal/KatherineK/data/scATAC/Liver-Cancer/",
                    list.files("/omics/groups/OE0219/internal/KatherineK/data/scATAC/Liver-Cancer/", pattern = "Samples.tsv.gz$"))

names(fragments) <- list.files("/omics/groups/OE0219/internal/KatherineK/data/scATAC/Liver-Cancer/", pattern = "Samples.tsv.gz$") %>% str_split("\\.") %>% lapply("[", 1) %>% unlist()

# create arrow files
ArrowFiles <- createArrowFiles(
  inputFiles = fragments,
  sampleNames = names(fragments),
  minTSS = 6, 
  minFrags = 3000,
  addTileMat = TRUE,
  addGeneScoreMat = TRUE
)

ArrowFiles <- list.files(pattern = "arrow")

# doublet prediction
doubScores <- addDoubletScores(
  input = ArrowFiles,
  k = 10, #Refers to how many cells near a "pseudo-doublet" to count.
  knnMethod = "UMAP", #Refers to the embedding to use for nearest neighbor search with doublet projection.
  LSIMethod = 1
)

# create an ArchR project
lica <- ArchRProject(
  ArrowFiles = ArrowFiles,
  outputDirectory = "/omics/groups/OE0219/internal/KatherineK/ATACseq/Liver-Cancer/malignant-subset/",
  copyArrows = TRUE #This is recommened so that if you modify the Arrow files you have an original copy for later usage.
)

# save ArchR project
saveArchRProject(ArchRProj = lica, load = FALSE)

# filter doublets
lica <- filterDoublets(lica)

# exclude cells from non-malignant clusters
celltypes <- read.table("/omics/groups/OE0219/internal/KatherineK/data/scATAC/Liver-Cancer/celltypes.txt", header = T, sep="\t")
exclude.cells <- celltypes$Cell_Barcode[celltypes$Predicted_Cell_Type!="Malignant cells"]
exclude.cells <- paste0("GSE227265_fragments_AllSamples#", exclude.cells) %>% intersect(lica$cellNames)
malignant.cells <- setdiff(lica$cellNames, exclude.cells)

lica <- subsetCells(ArchRProj = lica, cellNames = malignant.cells)

lica <- addIterativeLSI(
  ArchRProj = lica,
  useMatrix = "TileMatrix",
  name = "IterativeLSI",
  iterations = 2,
  clusterParams = list( #See Seurat::FindClusters
    resolution = c(0.2),
    sampleCells = 10000,
    n.start = 10),
  varFeatures = 25000,
  dimsToUse = 1:30)


# add clusters
lica <- addClusters(
  input = lica,
  reducedDims = "IterativeLSI",
  method = "Seurat",
  name = "Clusters",
  resolution = 0.8
)

# umap
lica <- addUMAP(
  ArchRProj = lica,
  reducedDims = "IterativeLSI",
  name = "UMAP",
  nNeighbors = 30,
  minDist = 0.5,
  metric = "cosine"
)

p1 <- plotEmbedding(ArchRProj = lica, colorBy = "cellColData", name = "Sample", embedding = "UMAP")
p2 <- plotEmbedding(ArchRProj = lica, colorBy = "cellColData", name = "Clusters", embedding = "UMAP")

ggAlignPlots(p1, p2, type = "h")

plotPDF(p1,p2, name = "Plot-UMAP-Sample-Clusters.pdf", ArchRProj = lica, addDOC = FALSE, width = 5, height = 5)

saveArchRProject(ArchRProj = lica, load = FALSE)


lica <- addGroupCoverages(lica)

# calling peaks
pathToMacs2 <- findMacs2()

lica <- addReproduciblePeakSet(
  ArchRProj = lica,
  pathToMacs2 = pathToMacs2
)

saveArchRProject(ArchRProj = lica, load = FALSE)

# add a peaks matrix
lica <- addPeakMatrix(lica)

mat <- getMatrixFromProject(
  ArchRProj = lica,
  useMatrix = "PeakMatrix",
  useSeqnames = NULL,
  verbose = TRUE,
  binarize = T,
  threads = getArchRThreads(),
  logFile = createLogFile("getMatrixFromProject")
)

saveRDS(mat, file = "peaks_matrix.Rds")

# gene score matrix
mat <- getMatrixFromProject(
  ArchRProj = lica,
  useMatrix  = "GeneScoreMatrix",
  useSeqnames = NULL,
  verbose = TRUE,
  binarize = T,
  threads = getArchRThreads(),
  logFile = createLogFile("getMatrixFromProject")
)

saveRDS(mat, "gene_score_matrix.Rds")



#--- downstream analysis on peaks matrix
lica <- loadArchRProject()
mat <- readRDS("peaks_matrix.Rds")
names(mat@rowRanges) <- 1:length(mat@rowRanges)

# select sites corresponding to hg38 promoters
# require("TxDb.Hsapiens.UCSC.hg38.knownGene")
# txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
# promoters <- promoters(genes(txdb), upstream = 1500, downstream = 500)
# promoter.peaks <- subsetByOverlaps(mat@rowRanges, promoters) %>% names() %>% as.numeric()
#mat <- mat[promoter.peaks,]
row.names <- mat@rowRanges %>% paste0()
dim(mat)
mat <- mat@assays@data$PeakMatrix
rownames(mat) <- row.names
mat[1:10,1:10]
clusters <- unique(lica$Clusters)

datasets <- list()
for (i in clusters) {
  #ids <- colnames(mat)[grepl(i, colnames(mat))]
  ids <- colnames(mat)[lica$Clusters==i]
  if (length(ids)<50) { next }
  datasets[[i]] <- mat[,ids]
  if (ncol(datasets[[i]])>100) { datasets[[i]] <- datasets[[i]][,sample(ncol(datasets[[i]]), 100)]  }
}

lapply(datasets, dim)

# compute epiCHAOS
het <- compute.eITH(datasets)

saveRDS(het, "epiCHAOS_scores_copy_corrected_all_peaks.Rds")


#--- compare epiCHAOS in peaks vs windows matrix
mat <- getMatrixFromProject(
  ArchRProj = lica,
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
mat <- mat[rowSums(mat)>30,]
dim(mat)
mat[mat>1] <- 1

sample.sites <- sample(nrow(mat), 500000)
mat <- mat[sample.sites,]
dim(mat)

#colnames(mat) <- paste0(lica$Sample, "_", lica$cellNames)
clusters <- unique(lica$Clusters)

datasets <- list()
for (i in clusters) {
  ids <- colnames(mat)[lica$Clusters==i]
  print(length(ids))
  datasets[[i]] <- mat[,ids]
  if (ncol(datasets[[i]])>100) { datasets[[i]] <- datasets[[i]][,sample(ncol(datasets[[i]]), 100)]  }
}

lapply(datasets, dim)

# compute heterogeneity scores
het2 <- compute.eITH(datasets)
saveRDS(het2, "epiCHAOS_scores_tilematrix.Rds")

temp <- data.frame(TileMatrix=het2$mean.het, PeakMatrix=het$mean.het)

pdf("correlate_epiCHAOS_tiles_vs_peaks_brca.pdf", 5, 4)
ggplot(temp, aes(y=TileMatrix, x=PeakMatrix)) +
  geom_point(size=3, alpha=0.7,  color="steelblue4")+
  stat_cor() +
  labs(x="epiCHAOS PeakMatrix", y="epiCHAOS TileMatrix") +
  theme_bw()
dev.off()

saveRDS(het2, file = "epiCHAOS_tilematrix.Rds")

#--- barplot for heterogeneity scores
pdf("barplot_epiCHAOS_percluster.pdf", 5, 2)
ggplot(het, aes(x = reorder(state, mean.het), y = mean.het)) +
  geom_bar(stat="identity", position = "dodge", alpha=0.8, width = 0.6)+
  labs(x="", y="epiCHAOS")+
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()


#--- calculate heterogeneity scores with count adjustement per chromosome

# tiles matrix
mat <- readRDS("tile_matrix.Rds")
row.names <- paste(mat@elementMetadata$seqnames, mat@elementMetadata$start, sep=":")
mat <- mat@assays@data$TileMatrix
rownames(mat) <- row.names
mat <- mat[rowSums(mat)>30,]
mat[mat>1] <- 1
dim(mat)
sample.sites <- sample(nrow(mat), 500000)
mat <- mat[sample.sites,]
dim(mat)

clusters <- unique(lica$Clusters)

datasets <- list()
for (i in clusters) {
  ids <- colnames(mat)[lica$Clusters==i]
  print(length(ids))
  datasets[[i]] <- mat[,ids]
  if (ncol(datasets[[i]])>100) { datasets[[i]] <- datasets[[i]][,sample(ncol(datasets[[i]]), 100)]  }
}

lapply(datasets, dim)


### source corrected eper chromosome function

# compute heterogeneity scores
het <- compute.eITH(datasets)
