
#--- this script contains quality control tests perfurmed on the epiCHAOS score, e.g. 
#--- (i) to check correlation between epiCHAOS score and quality metrics to see if it reflects signal to noise ratio
#--- (ii) to show that is is not influenced by the number of cells per cluster
#--- (iii) to show that it was correlated with total counts before coverage adjustment, and that this effect is removed after adjustment

# epiCHAOS function
source("/omics/groups/OE0219/internal/KatherineK/ATACseq/eITH-test-scripts/jaccard.R")

library(Seurat)
library(Signac)
library(EnsDb.Hsapiens.v86)
library(magrittr)

#--- use liver cancer cell line data, separate data into bins ordered by QC metrics and compare their heterogeneity scores
hcc <- readRDS("/omics/groups/OE0219/internal/KatherineK/ATACseq/HCC-celllines/Save-ArchR-Project.rds")
mat <- readRDS("/omics/groups/OE0219/internal/KatherineK/ATACseq/HCC-celllines/matrices/peak_matrix.Rds")

# select ids for one cell line
ids <- rownames(hcc@cellColData[hcc@cellColData$Sample=="Hep-1",])

#promoter.peaks <- hcc@peakSet$idx[hcc@peakSet$peakType=="Promoter"]
#mat <- mat@assays@data$PeakMatrix[as.numeric(promoter.peaks),ids]

# select a random 50,000 peaks
mat <- mat@assays@data$PeakMatrix[sample(nrow(mat), 50000),ids]
dim(mat)

qc <- data.frame(id= rownames(hcc@cellColData)[hcc$Sample=="Hep-1"], qc=hcc$TSSEnrichment[hcc$Sample=="Hep-1"])
qc$qc <- ntile(qc$qc, n = 20)

datasets <- list()
for (i in unique(qc$qc)) {
  ids <- qc$id[qc$qc==i] %>% sample(100)
  temp <- mat[,ids]
  datasets[[paste0("group-",i)]] <- temp
}

lapply(datasets, dim)

het <- compute.eITH(datasets)
het$qc <- het$state %>% str_remove("group-") %>% as.numeric()

barplot1 <- ggplot(het, aes(y=mean.het+0.01, x=reorder(state, qc))) +
  geom_bar(stat="identity", width = 0.3, fill="salmon")+
  labs( x="cells binned by TSS enrichment score", y="epiCHAOS") +
  lims(y=c(0,1.01))+
  theme_classic() +
  theme(axis.text.x = element_blank())

qc <- data.frame(id= rownames(hcc@cellColData)[hcc$Sample=="Hep-1"], qc=hcc$NucleosomeRatio[hcc$Sample=="Hep-1"])
qc$qc <- ntile(qc$qc, n = 20)

datasets <- list()
for (i in unique(qc$qc)) {
  ids <- qc$id[qc$qc==i] %>% sample(100)
  temp <- mat[,ids]
  datasets[[paste0("group-",i)]] <- temp
}

lapply(datasets, dim)

het <- compute.eITH(datasets)
het$qc <- het$state %>% str_remove("group-") %>% as.numeric()

barplot2 <- ggplot(het, aes(y=mean.het+0.01, x=reorder(state, qc))) +
  geom_bar(stat="identity", width = 0.3, fill="salmon")+
  labs( x="cells binned by nucleosome ratio", y="epiCHAOS") +
  lims(y=c(0,1.01))+
  theme_classic() +
  theme(axis.text.x = element_blank())

qc <- data.frame(id= rownames(hcc@cellColData)[hcc$Sample=="Hep-1"], qc=hcc$nFrags[hcc$Sample=="Hep-1"])
qc$qc <- ntile(qc$qc, n = 20)

datasets <- list()
for (i in unique(qc$qc)) {
  ids <- qc$id[qc$qc==i] %>% sample(100)
  temp <- mat[,ids]
  datasets[[paste0("group-",i)]] <- temp
}

lapply(datasets, dim)

het <- compute.eITH(datasets)
het$qc <- het$state %>% str_remove("group-") %>% as.numeric()

barplot3 <- ggplot(het, aes(y=mean.het+0.01, x=reorder(state, qc))) +
  geom_bar(stat="identity", width = 0.3, fill="salmon")+
  labs( x="cells binned by Doublet score", y="epiCHAOS") +
  lims(y=c(0,1.01))+
  theme_classic() +
  theme(axis.text.x = element_blank())

# FRIP
hcc$FRIP <- hcc$ReadsInPromoter/hcc$nFrags
qc <- data.frame(id= rownames(hcc@cellColData)[hcc$Sample=="Hep-1"], qc=hcc$FRIP[hcc$Sample=="Hep-1"])
qc$qc <- ntile(qc$qc, n = 20)

datasets <- list()
for (i in unique(qc$qc)) {
  ids <- qc$id[qc$qc==i] %>% sample(100)
  temp <- mat[,ids]
  datasets[[paste0("group-",i)]] <- temp
}

lapply(datasets, dim)

het <- compute.eITH(datasets)
het$qc <- het$state %>% str_remove("group-") %>% as.numeric()

barplot3 <- ggplot(het, aes(y=mean.het+0.01, x=reorder(state, qc))) +
  geom_bar(stat="identity", width = 0.3, fill="salmon")+
  labs( x="cells binned by FRIP score", y="epiCHAOS") +
  lims(y=c(0,1.01))+
  theme_classic() +
  theme(axis.text.x = element_blank())


#--- check if epiCHAOS is affected by the number of cells per cluster
ncells <- 25
datasets <- list()
for (i in 1:20) {
  ids <- qc$id %>% sample(ncells)
  temp <- mat[,ids]
  datasets[[paste0("group-",i)]] <- temp
  ncells <- ncells + 25
  print(ncells)
}

lapply(datasets, dim)
het <- compute.eITH(datasets)
het$clustersize <- seq(25,500,by=25)
het$n <- 1:20
plot(het$mean.het, het$clustersize)

barplot4 <- ggplot(het, aes(y=mean.het+0.01, x=reorder(state, n))) +
  geom_bar(stat="identity", width = 0.3, fill="salmon")+
  labs( x="cells binned by Cluster size", y="epiCHAOS") +
  lims(y=c(0,1.01))+
  theme_classic() +
  theme(axis.text.x = element_blank())


svg("/omics/groups/OE0219/internal/KatherineK/ATACseq/epiCHAOS-Figures/Supplementary/S1/QC_barplots_hep1.svg", 4, 7)
ggarrange(barplot3,barplot1, barplot2, barplot4, nrow=4)
dev.off()





#--- use the Granja et al. data to show that epiCHAOS correlated with the total counts before adjustment, and does not after adjustment

#--- load scATAC from normal hematopoiesis from Granja et al. 2019, subset for BM-derived
atac <- readRDS("/omics/groups/OE0219/internal/KatherineK/data/scATAC/Granja2019/scATAC-Healthy-Hematopoiesis-191120.rds")
atac <- atac[,grepl("BMMC", atac@colData$Group)]

hema <- atac@assays$data$counts # %>% as.matrix()

# create list of datasets per cell type for eITH calculation
datasets <- list()
n <- 100 # select n cells form each celltype
for (i in unique(atac$BioClassification[!grepl("Unk", atac$BioClassification)])) {
  cellnames <- colnames(atac)[atac$BioClassification==i]
  if (length(cellnames)>n) { cellnames <- sample(cellnames, n, replace = F) }
  datasets[[i]] <- hema[,colnames(hema) %in% cellnames]
  datasets[[i]][datasets[[i]]>1] <- 1
}

lapply(datasets, dim)

# epiCHAOS function before adjusting for total counts


# compute heterogeneity
het <- compute.eITH(datasets)



