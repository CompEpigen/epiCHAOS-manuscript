
#--- testing epiCHAOS in single cell ChIP data from breast cancer cells sensitive / resistant to Capecitabine (Grosselin et al. 2019)

# *** the below code for analysis of the scCHIP-seq data were adapted from https://github.com/vallotlab/scChIPseq ***

##Importing packages
library(scater)

## Initializing fixed parameters ## 

set.seed(10)

# Capacitabine resistance
path.1 <- "/omics/groups/OE0219/internal/KatherineK/data/scCHIP/Breastcancer/GSM3290889_CountTable_HBCx-95_scChIP_H3K27me3_hg38.txt"
path.2 <- "/omics/groups/OE0219/internal/KatherineK/data/scCHIP/Breastcancer/GSM3290890_CountTable_HBCx-95-CapaR_scChIP_H3K27me3_hg38.txt"


datamatrix=NULL
annot_raw = NULL

#QC & Filtering
min_coverage_cell = 1600
min_cells_window = 1
quant_removal = 95

#--- first counts matrix
datamatrix_single <- read.table(path.1, header=TRUE, stringsAsFactors=FALSE)
datamatrix_single <- datamatrix_single[!duplicated(rownames(datamatrix_single)),] #put IN for new format

# restructure to set rownames
rownames(datamatrix_single) = as.character(datamatrix_single[,1])
datamatrix_single = datamatrix_single [,-1]

# number of cells
total_cell <- length(datamatrix_single[1,])
sample_name <- "sensitive"
annot_single <- data.frame(barcode=colnames(datamatrix_single), cell_id=paste0(sample_name, "_c", 1:total_cell), sample_id=rep(sample_name, total_cell))

colnames(datamatrix_single) <- annot_single$cell_id

datamatrix <- datamatrix_single
annot_raw <- annot_single

#--- second counts matrix
datamatrix_single <- read.table(path.2, header=TRUE, stringsAsFactors=FALSE)
datamatrix_single <- datamatrix_single[!duplicated(rownames(datamatrix_single)),] #put IN for new format

# restructure to set rownames
rownames(datamatrix_single) = as.character(datamatrix_single[,1])
datamatrix_single = datamatrix_single [,-1]

# number of cells
total_cell <- length(datamatrix_single[1,])
sample_name <- "resistant"
annot_single <- data.frame(barcode=colnames(datamatrix_single), cell_id=paste0(sample_name, "_c", 1:total_cell), sample_id=rep(sample_name, total_cell))

colnames(datamatrix_single) <- annot_single$cell_id

# merge matrices
common_regions <- intersect(rownames(datamatrix), rownames(datamatrix_single))
datamatrix <- cbind(datamatrix[common_regions,], datamatrix_single[common_regions,])
annot_raw <- rbind(annot_raw, annot_single)
dim(datamatrix)

#Removing non-standard chromosomes
splitID <- sapply(rownames(datamatrix), function(x) strsplit(as.character(x), split="_"))
normalChr <- which(sapply(splitID, length) <= 3) # weird chromosomes contain underscores in the name
datamatrix <- datamatrix[normalChr,]
dim(datamatrix)

#Remove chrM and chrY from mattrix
if(length(grep("chrM",rownames(datamatrix)))>0)  datamatrix <- datamatrix[-grep("chrM",rownames(datamatrix)),]
if(length(grep("chrY",rownames(datamatrix)))>0)  datamatrix <- datamatrix[-grep("chrY",rownames(datamatrix)),]


## Filtering and QC

# create a single cell experiment
assays = list(counts = as.matrix(datamatrix)) #, colData =  annot_raw
umi <- SingleCellExperiment(assays = assays, colData = annot_raw)
umi <- umi[rowSums(counts(umi)>0)>0, ] # remove windows that do not have any read in any cells
qc <- scater::perCellQCMetrics(umi)
thresh <- quantile(colSums(counts(umi)), probs=seq(0, 1, 0.01))

## Filtering & Window selection

#Cell selection based on number of total counts, between min_cov_cell and upper 5%
sel1000 =  (colSums(counts(umi))>1000 & colSums(counts(umi))< thresh[input$quant_removal+1])

sel <- ( colSums(counts(umi))>min_coverage_cell & colSums(counts(umi)) < thresh[quant_removal+1] )

SelMatCov1000 <- counts(umi)[,sel1000]
bina_counts <- SelMatCov1000
bina_counts[bina_counts<2] <-0
bina_counts[bina_counts>1] <-1
fixedWin <- names(which((rowSums(bina_counts) > ( (min_cells_window/100.0)*(dim(bina_counts)[2])) ))) # window selection for regions covered in > 1% of cells
length(fixedWin)

SelMatCov <- counts(umi)[,sel]
SelMatCov <- SelMatCov[fixedWin,]

annot <- colData(umi)
annot <- as.data.frame(annot[sel,])
#annot = cbind(annot,annot_raw[which(annot_raw$cell_id %in% rownames(annot)),])


#--- matrices for epiCHAOS
mat <- SelMatCov
dim(SelMatCov)

mat[mat>1] <- 1

mat[,grepl("resistant", colnames(mat))] %>% dim()
mat[,grepl("sensitive", colnames(mat))] %>% dim()

set.seed(11)

#--- subsample groups of 100 cells from each condition for epiCHAOS calculation
datasets <- list()
for (i in 1:10) {
  samples <- sample(ncol(mat[,grepl("resistant", colnames(mat))]), 100) %>% as.numeric()
  datasets[[paste0("resistant-", i)]] <- mat[,samples]
  
  samples <- sample(ncol(mat[,grepl("sensitive", colnames(mat))]), 100) %>% as.numeric()
  datasets[[paste0("sensitive-", i)]] <- mat[,samples]
  
}

lapply(datasets, dim)
lapply(datasets, sum)


#--- compute epiCHAOS scores
het <- compute.eITH(datasets)
het$group <- "resistant"
het$group[grepl("sensitive", het$state)] <- "sensitive"
het$state <- het$state %>% str_replace("resistant", "R") %>% str_replace("sensitive", "S")

#--- create boxplot
p2 <- ggplot(het, aes(y=mean.het+0.01, x=reorder(group, mean.het))) +
  geom_boxplot(fill="lightsteelblue3")+
  geom_jitter(size=0.5, width = 0.1)+
  ggpubr::stat_compare_means(comparisons = list(c(1,2)))+
  labs( x="", y="epiCHAOS")+
  theme_classic() +
  theme(axis.text.x = element_text(angle=90))


# plot 
svg("/omics/groups/OE0219/internal/KatherineK/ATACseq/epiCHAOS-Figures/Figure 5/barplot_epiCHAOS_scCHIP_CapR.svg", 7, 2.5)
ggpubr::ggarrange(p1,p2, widths = c(5, 2), align = "h")
dev.off()
