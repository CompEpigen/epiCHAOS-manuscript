
#--- this script contains quality control tests performed on the epiCHAOS score, e.g. (related to Figure S1)
#--- - to check correlation between epiCHAOS score and quality metrics to see if it reflects signal to noise ratio
#--- - to show that is is not influenced by the number of cells per cluster

#--- load packages
library(EnsDb.Hsapiens.v86)
library(magrittr)
library(ArchR)

#--- use liver cancer cell line data, separate data into bins ordered by QC metrics and compare their heterogeneity scores
hcc <- loadArchRProject("/omics/groups/OE0219/internal/KatherineK/ATACseq/HCC-celllines/")
mat <- readRDS("/omics/groups/OE0219/internal/KatherineK/ATACseq/HCC-celllines/matrices/peak_matrix.Rds")

#--- select ids for one cell line so as to have a group of cells that we expect to be relatively homogeneous
ids <- rownames(hcc@cellColData[hcc@cellColData$Sample=="Hep-1",])

#--- subset for a random 50,000 peaks to ease computation
mat <- mat@assays@data$PeakMatrix[sample(nrow(mat), 50000),ids]
dim(mat)

#--- TSS enrichment score
qc <- data.frame(id= rownames(hcc@cellColData)[hcc$Sample=="Hep-1"], qc=hcc$TSSEnrichment[hcc$Sample=="Hep-1"])
qc$qc <- ntile(qc$qc, n = 20)

#--- create 20 bins of cells ordered by TSS enrichment score
datasets <- list()
for (i in unique(qc$qc)) {
  ids <- qc$id[qc$qc==i] %>% sample(100)
  temp <- mat[,ids]
  datasets[[paste0("group-",i)]] <- temp
}

lapply(datasets, dim)

#--- compute epiCHAOS score
het <- compute.eITH(datasets)
het$qc <- het$state %>% str_remove("group-") %>% as.numeric()

barplot1 <- ggplot(het, aes(y=mean.het+0.01, x=reorder(state, qc))) +
  geom_bar(stat="identity", width = 0.3, fill="salmon")+
  labs( x="cells binned by TSS enrichment score", y="epiCHAOS") +
  lims(y=c(0,1.01))+
  theme_classic() +
  theme(axis.text.x = element_blank())

#--- nucleosome ratio
qc <- data.frame(id= rownames(hcc@cellColData)[hcc$Sample=="Hep-1"], qc=hcc$NucleosomeRatio[hcc$Sample=="Hep-1"])
qc$qc <- ntile(qc$qc, n = 20)

#--- create 20 bins of cells ordered by nucleosome ratio
datasets <- list()
for (i in unique(qc$qc)) {
  ids <- qc$id[qc$qc==i] %>% sample(100)
  temp <- mat[,ids]
  datasets[[paste0("group-",i)]] <- temp
}

lapply(datasets, dim)

#--- compute epiCHAOS score
het <- compute.eITH(datasets)
het$qc <- het$state %>% str_remove("group-") %>% as.numeric()

barplot2 <- ggplot(het, aes(y=mean.het+0.01, x=reorder(state, qc))) +
  geom_bar(stat="identity", width = 0.3, fill="salmon")+
  labs( x="cells binned by nucleosome ratio", y="epiCHAOS") +
  lims(y=c(0,1.01))+
  theme_classic() +
  theme(axis.text.x = element_blank())


#--- FRIP score
hcc$FRIP <- hcc$ReadsInPromoter/hcc$nFrags
qc <- data.frame(id= rownames(hcc@cellColData)[hcc$Sample=="Hep-1"], qc=hcc$FRIP[hcc$Sample=="Hep-1"])
qc$qc <- ntile(qc$qc, n = 20)

#--- create 20 bins of cells ordered by FRIP score
datasets <- list()
for (i in unique(qc$qc)) {
  ids <- qc$id[qc$qc==i] %>% sample(100)
  temp <- mat[,ids]
  datasets[[paste0("group-",i)]] <- temp
}

lapply(datasets, dim)

#--- compute epiCHAOS scores
het <- compute.eITH(datasets)
het$qc <- het$state %>% str_remove("group-") %>% as.numeric()

barplot3 <- ggplot(het, aes(y=mean.het+0.01, x=reorder(state, qc))) +
  geom_bar(stat="identity", width = 0.3, fill="salmon")+
  labs( x="cells binned by FRIP score", y="epiCHAOS") +
  lims(y=c(0,1.01))+
  theme_classic() +
  theme(axis.text.x = element_blank())


#--- check if epiCHAOS is affected by the number of cells per group/cluster
#--- create 20 groups of cells with increasing number of cells from 25 to 500
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

#--- compute epiCHAOS scores
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

ggarrange(barplot3,barplot1, barplot2, barplot4, nrow=4)


