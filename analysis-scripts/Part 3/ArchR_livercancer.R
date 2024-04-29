

#--- liver cancer scATAC-seq analysis with excluding non-malignant cells

#--- scATAC data from liver cancer from Craig et al 2023: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE227265

setwd("/omics/groups/OE0219/internal/KatherineK/ATACseq/Liver-Cancer/malignant-subset/")

#--- directory containing fragments files, and for analysis
data.dir <- "/omics/groups/OE0219/internal/KatherineK/data/scATAC/Liver-Cancer/"
analysis.dir <- "/omics/groups/OE0219/internal/KatherineK/ATACseq/Liver-Cancer/malignant-subset/"

#--- load ArchR packages
library(ArchR)
library(BSgenome.Hsapiens.UCSC.hg38)

addArchRThreads(threads = 16)
addArchRGenome("hg38")

#--- fragments files
fragments <- paste0(file.path(data.dir, "scATAC/Liver-Cancer/"), list.files(file.path(data.dir, "scATAC/Liver-Cancer/"), pattern = "Samples.tsv.gz$"))
names(fragments) <- list.files(file.path(data.dir, "scATAC/Liver-Cancer/"), pattern = "Samples.tsv.gz$") %>% str_split("\\.") %>% lapply("[", 1) %>% unlist()

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
lica <- ArchRProject(
  ArrowFiles = ArrowFiles,
  outputDirectory = file.path(analysis.dir, "malignant-subset/"),
  copyArrows = TRUE #This is recommened so that if you modify the Arrow files you have an original copy for later usage.
)

#--- save ArchR project
saveArchRProject(ArchRProj = lica, load = FALSE)

#--- filter doublets
lica <- filterDoublets(lica)


#--- exclude cells annotated to non-malignant clusters in the Craig et al. publication
celltypes <- read.table(file.path(data.dir, "celltypes.txt"), header = T, sep="\t")
exclude.cells <- celltypes$Cell_Barcode[celltypes$Predicted_Cell_Type!="Malignant cells"]
exclude.cells <- paste0("GSE227265_fragments_AllSamples#", exclude.cells) %>% intersect(lica$cellNames)
malignant.cells <- setdiff(lica$cellNames, exclude.cells)

lica <- subsetCells(ArchRProj = lica, cellNames = malignant.cells)

#--- add iterative LSI
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


#--- add clusters
lica <- addClusters(
  input = lica,
  reducedDims = "IterativeLSI",
  method = "Seurat",
  name = "Clusters",
  resolution = 0.8
)

#--- umap
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

#--- add group coverages
lica <- addGroupCoverages(lica)

#--- calling peaks with macs2
pathToMacs2 <- findMacs2()

lica <- addReproduciblePeakSet(
  ArchRProj = lica,
  pathToMacs2 = pathToMacs2
)

saveArchRProject(ArchRProj = lica, load = FALSE)

#--- add a peaks matrix
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

#--- gene score matrix
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

row.names <- mat@rowRanges %>% paste0()
dim(mat)
mat <- mat@assays@data$PeakMatrix
rownames(mat) <- row.names
mat[1:10,1:10]
clusters <- unique(lica$Clusters)

#--- peaks-by-cells matrices for each cluster, subset for 100 cells for larger clusters
datasets <- list()
for (i in clusters) {
  ids <- colnames(mat)[lica$Clusters==i]
  if (length(ids)<50) { next }
  datasets[[i]] <- mat[,ids]
  if (ncol(datasets[[i]])>100) { datasets[[i]] <- datasets[[i]][,sample(ncol(datasets[[i]]), 100)]  }
}

lapply(datasets, dim)

#--- compute and save epiCHAOS scores, use "cancer" version to correct for per-chromosome counts
het <- compute.eITH(datasets)
het <- compute.eITH.cancer(datasets)

saveRDS(het, "epiCHAOS_scores_copy_corrected_all_peaks.Rds")

