
#--- comparison of epigenetic heterogeneity in old vs young HSCs from Itowaka et al. https://www.nature.com/articles/s41467-022-30440-2#Sec11

library(Signac)
library(Seurat)
library(EnsDb.Hsapiens.v75)
library(ggplot2)
library(patchwork)
library(plyr)
library(ggbeeswarm)

#--- create seurat object from young HSCs
counts <- Read10X_h5(filename="/omics/groups/OE0219/internal/KatherineK/data/scATAC/HspcAging/Itokawa-2022/GSM5723631_Young_HSC_filtered_peak_bc_matrix.h5")
meta.data <- read.csv("/omics/groups/OE0219/internal/KatherineK/data/scATAC/HspcAging/Itokawa-2022/GSM5723631_Young_HSC_singlecell.csv", header = T, row.names = 1)
chrom_assay <- CreateChromatinAssay(counts=counts, sep=c(":", "-"), fragments = "/omics/groups/OE0219/internal/KatherineK/data/scATAC/HspcAging/Itokawa-2022/GSM5723631_Young_HSC_fragments.tsv.gz", min.cells = 10, min.features = 200)
young <- CreateSeuratObject(counts = chrom_assay, assay="peaks", meta.data = meta.data)

young[["peaks"]]
granges(young)

# extract gene annotations from EnsDb
#BiocManager::install("EnsDb.Mmusculus.v79")
require("EnsDb.Mmusculus.v79")
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
# change to UCSC style since the data was mapped to hg19
seqlevels(annotations) <- paste0('chr', seqlevels(annotations))
genome(annotations) <- "mm10"

# add the gene information to the object
Annotation(young) <- annotations

# computing QC metrics
young <- NucleosomeSignal(object = young)
young <- TSSEnrichment(object=young, fast = F)

# add blacklist ratio and fraction of reads in peaks
young$pct_reads_in_peaks <- young$peak_region_fragments / young$passed_filters * 100
young$blacklist_ratio <- young$blacklist_region_fragments / young$peak_region_fragments

DensityScatter(young, x = 'nCount_peaks', y = 'TSS.enrichment', log_x = TRUE, quantiles = TRUE)

young <- subset(
  x = young,
  subset = nCount_peaks > 3000 &
    nCount_peaks < 20000 &
    pct_reads_in_peaks > 15 &
    blacklist_ratio < 0.05 &
    nucleosome_signal < 4 #&
    #TSS.enrichment > 3
)

# normalisation and linear dimensionality reduction
young <- RunTFIDF(young)
young <- FindTopFeatures(young, min.cutoff = 'q0')
young <- RunSVD(young)

# check if the first LSI component captures sequencing depth, if so remove it
DepthCor(young)

# non-linear dimension reduction
young <- RunUMAP(object = young, reduction = 'lsi', dims = 2:10)
young <- FindNeighbors(object = young, reduction = 'lsi', dims = 2:10)
young <- FindClusters(object = young, verbose = FALSE, algorithm = 3)
DimPlot(object = young, label = F) + NoLegend()

# save signac object
saveRDS(young, "/omics/groups/OE0219/internal/KatherineK/data/scATAC/HspcAging/Itokawa-2022/young_hscs_signac.Rds")


#--- create seurat object from old HSCs
counts <- Read10X_h5(filename="/omics/groups/OE0219/internal/KatherineK/data/scATAC/HspcAging/Itokawa-2022/GSM5723632_Aged_HSC_filtered_peak_bc_matrix.h5")
meta.data <- read.csv("/omics/groups/OE0219/internal/KatherineK/data/scATAC/HspcAging/Itokawa-2022/GSM5723632_Aged_HSC_singlecell.csv", header = T, row.names = 1)
chrom_assay <- CreateChromatinAssay(counts=counts, sep=c(":", "-"), fragments = "/omics/groups/OE0219/internal/KatherineK/data/scATAC/HspcAging/Itokawa-2022/GSM5723632_Aged_HSC_fragments.tsv", min.cells = 10, min.features = 200)
aged <- CreateSeuratObject(counts = chrom_assay, assay="peaks", meta.data = meta.data)

aged[["peaks"]]
granges(aged)

# already done above
# # extract gene annotations from EnsDb
# BiocManager::install("EnsDb.Mmusculus.v79")
# annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
# # change to UCSC style since the data was mapped to hg19
# seqlevels(annotations) <- paste0('chr', seqlevels(annotations))
# genome(annotations) <- "mm10"

# add the gene information to the object
Annotation(aged) <- annotations

# computing QC metrics
aged <- NucleosomeSignal(object = aged)
aged <- TSSEnrichment(object=aged, fast = F)

# add blacklist ratio and fraction of reads in peaks
aged$pct_reads_in_peaks <- aged$peak_region_fragments / aged$passed_filters * 100
aged$blacklist_ratio <- aged$blacklist_region_fragments / aged$peak_region_fragments

DensityScatter(aged, x = 'nCount_peaks', y = 'TSS.enrichment', log_x = TRUE, quantiles = TRUE)

aged <- subset(
  x = aged,
  subset = nCount_peaks > 3000 &
    nCount_peaks < 20000 &
    pct_reads_in_peaks > 15 &
    blacklist_ratio < 0.05 &
    nucleosome_signal < 4 #&
  #TSS.enrichment > 3
)

# normalisation and linear dimensionality reduction
aged <- RunTFIDF(aged)
aged <- FindTopFeatures(aged, min.cutoff = 'q0')
aged <- RunSVD(aged)

# check if the first LSI component captures sequencing depth, if so remove it
DepthCor(aged)

# non-linear dimension reduction
aged <- RunUMAP(object = aged, reduction = 'lsi', dims = 2:10)
aged <- FindNeighbors(object = aged, reduction = 'lsi', dims = 2:10)
aged <- FindClusters(object = aged, verbose = FALSE, algorithm = 3)
DimPlot(object = aged, label = F) + NoLegend()

# save signac object
saveRDS(aged, "/omics/groups/OE0219/internal/KatherineK/data/scATAC/HspcAging/Itokawa-2022/aged_hscs_signac.Rds")



#---
#--- downstream analyses of eITH
#---

young <- readRDS("/omics/groups/OE0219/internal/KatherineK/data/scATAC/HspcAging/Itokawa-2022/young_hscs_signac.Rds")
aged <- readRDS("/omics/groups/OE0219/internal/KatherineK/data/scATAC/HspcAging/Itokawa-2022/aged_hscs_signac.Rds")

# create list of counts matrices for eITH analysis
datasets <- list()

# young HSCs
counts <- young@assays$peaks@counts
counts[counts>1] <- 1
datasets$young <- counts

# aged HSCs
counts <- aged@assays$peaks@counts
counts[counts>1] <- 1
datasets$aged <- counts
lapply(datasets, dim)

#--- overlap peaks from old and young
young.gr <- granges(young)
aged.gr <- granges(aged)
names(young.gr) <- 1:length(young.gr)
names(aged.gr) <- 1:length(aged.gr)
y.in.a <- subsetByOverlaps(young.gr, aged.gr) %>% names() %>% as.numeric() # index of peaks in young which overlap in aged
a.in.y <- subsetByOverlaps(aged.gr, young.gr) %>% names() %>% as.numeric() # index of peaks in aged which overlap in young

datasets$young <- datasets$young[y.in.a, ]
datasets$aged <- datasets$aged[a.in.y, ]

lapply(datasets, dim)

datasets.small <- list()

for (i in 1:5) {
  datasets.small[[paste0("young-", i)]] <- datasets$young[,sample(ncol(datasets$young), 100)] %>% as.matrix()
  datasets.small[[paste0("aged-", i)]] <- datasets$aged[,sample(ncol(datasets$aged), 100)] %>% as.matrix()
}

lapply(datasets.small, dim)

het <- compute_eITH(datasets.small)
het$state <- het$state %>% str_remove("[[:digit:]]")

pdf("Aging/boxplot_epiCHAOS_young_vs_old_HSCs_Itokawa.pdf", 3, 4)
ggplot(het, aes(y=het.adj, x=reorder(state, het.adj))) +
  geom_boxplot(fill="honeydew3", color="black", linewidth=0.5)+
  labs( x="", y="epiCHAOS") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()
