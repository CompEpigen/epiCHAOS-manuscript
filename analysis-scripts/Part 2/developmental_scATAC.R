
#--- A. test epiCHAOS in four developmental contexts and create sorted dotplots
#--- (i) hematopoietic system
#--- (ii) mouse gastrulation
#--- (iii) drosophila embryogenesis

#------------------------------------------------------------------------------

set.seed(10)

#--- load packages
library(magrittr)
library(stringr)
library(ggplot2)
library(ggbeeswarm)
library(ggpubr)

source("/omics/groups/OE0219/internal/KatherineK/ATACseq/eITH-test-scripts/jaccard.R")
source("~/Scripts/miscellaneous/myfunctions.R")

#--- (i) hematopoietic system

#--- load scATAC from normal hematopoiesis from Granja et al. 2019, subset for BM-derived
atac <- readRDS("/omics/groups/OE0219/internal/KatherineK/data/scATAC/Granja2019/scATAC-Healthy-Hematopoiesis-191120.rds")
atac <- atac[,grepl("BMMC", atac@colData$Group)]

hema <- atac@assays$data$counts # %>% as.matrix()

# # select sites corresponding to hg19 promoters
# require("TxDb.Hsapiens.UCSC.hg19.knownGene")
# txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
# promoters <- promoters(genes(txdb), upstream = 1500, downstream = 500)
# promoters.atac <- subsetByOverlaps(atac@rowRanges, promoters) %>% names()
# rownames(hema) <- names(atac@rowRanges)
# hema <- hema[promoters.atac,]
#hema[hema>1] <- 1

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

# compute heterogeneity
het <- compute.eITH(datasets)
saveRDS(het, "/omics/groups/OE0219/internal/KatherineK/ATACseq/Hematopoietic/epiCHAOS_scores_allpeaks.Rds")

het <- readRDS("/omics/groups/OE0219/internal/KatherineK/ATACseq/Hematopoietic/epiCHAOS_scores_allpeaks.Rds")

het$state <- het$state %>% str_split("_") %>% lapply("[", 2) %>% str_replace_all("\\.", " ")
het$progenitor <- ifelse(het$state %in% c("HSC", "GMP", "GMP Neut", "CLP 1", "CLP 2", "CMP LMPP", "Early Eryth"), T, F)

# create dotplot
dotplot1 <- ggplot(het, aes(x = reorder(state, mean.het), y = mean.het, fill=progenitor, color=progenitor)) +
  geom_dotplot(binaxis = "y", stackdir = "center", position = "dodge", alpha=0.7)+
  labs(x="", y="epiCHAOS")+
  scale_color_manual(values = rev(c("steelblue4", "black")))+
  scale_fill_manual(values = rev(c("steelblue4", "black")))+
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

barplot1 <- ggplot(het, aes(x = reorder(state, mean.het), y = mean.het, fill=progenitor, color=progenitor)) +
  geom_bar(stat="identity", position = "dodge", alpha=0.8, width = 0.7)+
  labs(x="", y="epiCHAOS")+
  scale_color_manual(values = rev(c("steelblue4", "black")))+
  scale_fill_manual(values = rev(c("steelblue4", "black")))+
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

# create UMAP
atac <- readRDS("/omics/groups/OE0219/internal/KatherineK/data/scATAC/Granja2019/scATAC-Healthy-Hematopoiesis-191120.rds")
atac@colData$celltype <- atac@colData$BioClassification %>% str_split("_") %>% lapply("[", 2) %>% str_replace_all("\\.", " ")
temp <- data.frame(umap1=atac@colData$UMAP1, umap2=atac@colData$UMAP2, celltype=atac@colData$celltype) %>% merge(unique(het[,c("state", "mean.het")]), by.x="celltype", by.y="state")
temp$epiCHAOS <- temp$mean.het
umap1 <- ggplot(temp, aes(y=umap2, x=umap1, color=epiCHAOS)) +
  geom_point(size=0.1, alpha=0.7)+
  scale_color_distiller(palette = "Blues", direction = 1)+
  labs( x="UMAP1", y="UMAP2") +
  theme_classic()

#--- correlate with cytoTRACE
cytotrace <- readRDS("/omics/groups/OE0219/internal/KatherineK/ATACseq/eICH-vs-tICH/cytotrace_Granja_2019_hematopoietic_BM.Rds")
cytotrace$CytoTRACErank

rna <- readRDS("/omics/groups/OE0219/internal/KatherineK/data/scATAC/Granja2019/scRNA-Healthy-Hematopoiesis-191120.rds")
rna <- rna[,grepl("BMMC", rna$Group)]
rna@colData$cytotrace <- cytotrace$CytoTRACErank
cyto <- data.frame(rna@colData)

for (i in unique(cyto$BioClassification)) {
  cyto$mean.cytotrace[cyto$BioClassification==i] <- median(cyto$cytotrace[cyto$BioClassification==i])
}

cyto <- cyto[,c("BioClassification", "mean.cytotrace")] %>% unique()
cyto$BioClassification <- cyto$BioClassification %>% str_split("_") %>% lapply("[", 2) %>% unlist() %>% str_replace_all("\\.", " ")
cyto <- merge(cyto, het, by.x="BioClassification", by.y="state")
cyto$progenitor
cyto$label <- cyto$BioClassification
fit <- lm(cyto$mean.het~cyto$mean.cytotrace)
cyto$label[cyto$progenitor==F] <- ""

pseudotime1 <- ggplot(cyto, aes(x=mean.cytotrace, y=mean.het, label=label)) +
  geom_point(size=3, alpha=0.7, aes(color=progenitor))+
  scale_color_manual(values = c("black", "steelblue4"))+
  #stat_cor()+
  lims(y=c(0, 1)) +
  geom_smooth(color="grey30", method = "lm", se=F)+
  #ggrepel::geom_text_repel(aes(color=progenitor), size=3, max.overlaps = 10)+
  labs( x="CytoTRACE", y="epiCHAOS", subtitle = paste0("R^2 = ", signif(summary(fit)$r.squared, 3))) +
  theme_classic()


#--- (ii) mouse gastrulation

#--- read atac data from mouse gastrulation
atac <- readRDS("/omics/groups/OE0219/internal/KatherineK/data/scATAC/mouse-gastrulation/Save-ArchR-Project.rds")
mat <- readRDS("/omics/groups/OE0219/internal/KatherineK/data/scATAC/mouse-gastrulation/Matrices/PeakMatrix_summarized_experiment.rds")

atac.gr <- mat@rowRanges
names(atac.gr) <- 1:length(atac.gr)

mat <- mat@assays@data$PeakMatrix

# select sites corresponding to mm10 promoters
#require("TxDb.Mmusculus.UCSC.mm10.knownGene")
#txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
#promoters <- promoters(genes(txdb), upstream = 1500, downstream = 500)
#rownames <- paste0(data.frame(subsetByOverlaps(atac.gr, promoters))$seqnames, "-", data.frame(subsetByOverlaps(atac.gr, promoters))$start, "-", data.frame(subsetByOverlaps(atac.gr, promoters))$end)
#promoters <- subsetByOverlaps(atac.gr, promoters) %>% names() %>% as.numeric()
#mat <- mat[promoters,]
#rownames(mat) <- rownames

atac$group <- atac$celltype.mapped

datasets <- list()
for (celltype in unique(na.omit(atac$group))) {
  print(celltype)
  ids <- atac$cellNames[!is.na(atac$group) & atac$group==celltype & atac$genotype=="WT"] %>% intersect(colnames(mat))
  if (length(ids)<30) { next }
  if (length(ids)>100) { ids <- sample(ids, 100)}
  datasets[[celltype]] <- mat[,colnames(mat) %in% ids]
  datasets[[celltype]][datasets[[celltype]]>1] <- 1
}

lapply(datasets, dim)
length(datasets)

# compute eICH
het <- compute.eITH(datasets)
saveRDS(het, "/omics/groups/OE0219/internal/KatherineK/ATACseq/embryogenesis/epiCHAOS_scores_allpeaks.Rds")


het <- readRDS("/omics/groups/OE0219/internal/KatherineK/ATACseq/embryogenesis/epiCHAOS_scores_allpeaks.Rds")
het$state <- het$state %>% str_replace_all("_", " ")
het$group <- ifelse(grepl("Primitive|Epiblast|PGC", het$state), "primitive", "other")


# create dotplot
dotplot2 <- ggplot(het, aes(x = reorder(state, mean.het), y = mean.het, fill=group, color=group)) +
  geom_dotplot(binaxis = "y", stackdir = "center", position = "dodge", alpha=0.6)+
  labs(x="", y="epiCHAOS")+
  scale_color_manual(values = rev(c("steelblue4", "black")))+
  scale_fill_manual(values = rev(c("steelblue4", "black")))+
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


barplot2 <- ggplot(het, aes(x = reorder(state, mean.het), y = mean.het, fill=group, color=group)) +
  geom_bar(stat="identity", position = "dodge", alpha=0.8, width = 0.7)+
  labs(x="", y="epiCHAOS")+
  scale_color_manual(values = rev(c("steelblue4", "black")))+
  scale_fill_manual(values = rev(c("steelblue4", "black")))+
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

#--- plot umap using reference cells from https://www.nature.com/articles/s41586-019-0933-9#MOESM3
ref <- read.csv("/omics/groups/OE0219/internal/KatherineK/data/scATAC/mouse-gastrulation/reference_cells_metadata.csv")

ref$celltype <- ref$celltype %>% str_replace_all("/", " ")

for (i in unique(ref$celltype)) {
  if (i %notin% het$state) { print(i); next }
  ref$het[ref$celltype==i] <- unique(het$mean.het[het$state==i])
}

ref <- ref[!is.na(ref$het),]
umap2 <- ggplot(ref[!is.na(ref$het),], aes(x=umapX, y=umapY, color=het)) +
  geom_point(size=0.1, alpha=0.7)+
  scale_color_distiller(palette = "Blues", direction = 1)+
  labs( x="UMAP1", y="UMAP2") +
  theme_classic()


#--- correlate epiCHAOS with developmental stage
temp <- ref[,c("stage", "het", "celltype")] %>% unique()
temp$stage <- temp$stage %>% str_remove("E") %>% as.numeric()
temp <- na.omit(temp)

for (i in unique(temp$celltype)) {
  temp$average.stage[temp$celltype==i] <- mean(temp$stage[temp$celltype==i])
}
temp$stage <- NULL
temp <- unique(temp)

temp$highlight <- ifelse(temp$celltype %in% c("Anterior Primitive Streak", "Primitive Streak", "Epiblast", "Nascent mesoderm"), T, F)
temp$label <- ifelse(temp$highlight==T|temp$celltype %in% c("Nascent mesoderm"), temp$celltype, "")
fit <- lm(temp$het~temp$average.stage)

pseudotime2 <- ggplot(temp, aes(x=average.stage, y=het, label=label)) +
  geom_point(size=3, alpha=0.5, aes(color=highlight))+
  scale_color_manual(values = c("black", "steelblue4"))+
  #ggrepel::geom_text_repel(size=3, aes(color=highlight))+
  geom_smooth(color="grey30", method="lm", se = F)+
  #stat_cor()+
  lims(y=c(0, 1))+
  labs( x="Developmental stage", y="epiCHAOS", subtitle = paste0("R^2 = ", signif(summary(fit)$r.squared, 3))) +
  theme_classic()

#--- (iii) drosophila embryogenesis

# read in atac seurat object
atac <- readRDS("/omics/groups/OE0219/internal/KatherineK/data/scATAC/embryogenesis/fly.all.downsampled_seurat_filtered_processed.rds")

# metadata including cellty/tissue annotation
temp <- readRDS("/omics/groups/OE0219/internal/KatherineK/data/scATAC/embryogenesis/atac_meta.rds")
temp <- temp[,c("NNv1_age", "NNv1_time.new", "refined_annotation", "seurat_clusters")]

# samples at time windows from 0-2, 2-4, 4-6 hours etc. until 20 hours
atac@meta.data <- atac@meta.data %>% merge(temp, by=0) %>% keepRow()
atac$group <- atac$refined_annotation
samples <- atac$group %>% unique()

## select regions corresponding to dm6 promoters
# require("TxDb.Dmelanogaster.UCSC.dm6.ensGene")
# txdb <- TxDb.Dmelanogaster.UCSC.dm6.ensGene
# promoters <- promoters(genes(txdb), upstream = 1500, downstream = 500)
# row.ranges <- rownames(atac) %>% str_replace("-", ":") %>% GRanges()
# names(row.ranges) <- 1:nrow(atac)
# promoters <- subsetByOverlaps(row.ranges, promoters) %>% names() %>% as.numeric()

atac@meta.data$refined_annotation <- as.character(atac@meta.data$refined_annotation) # 
atac@meta.data$refined_annotation[atac@meta.data$refined_annotation=="Unknown"] <- "Undifferentiated"

# create a peaks x cells dataset for each celltype/timepoint
datasets <- list()
for (i in unique(atac@meta.data$group)) {
  print(i)
  cells.select <- atac@meta.data[atac@meta.data$group==i,] %>% rownames() %>% sample(100)
  #temp <- atac@assays$RNA@data[promoters,cells.select]
  temp <- atac@assays$RNA@data[,cells.select]
  temp[temp>1] <- 1
  if (ncol(temp)<100) { next }
  datasets[[i]] <- temp[,sample(ncol(temp), 100)]
}

lapply(datasets, dim)

# compute heterogeneity
het <- compute.eITH(datasets)

saveRDS(het, "/omics/groups/OE0219/internal/KatherineK/ATACseq/mouse-gastrulation/epiCHAOS_scores_allpeaks.Rds")
het <- readRDS("/omics/groups/OE0219/internal/KatherineK/ATACseq/mouse-gastrulation/epiCHAOS_scores_allpeaks.Rds")

het$state[het$state=="Unknown"] <- atac@meta.data$refined_annotation[atac@meta.data$refined_annotation=="Unknown"] <- "Undifferentiated"
het$highlight <- ifelse(het$state %in% c("Undifferentiated", "Germ cell", "Blastoderm"), T, F)

# create dotplot
dotplot3 <- ggplot(het, aes(x = reorder(state, mean.het), y = mean.het, fill=highlight, color=highlight)) +
  geom_dotplot(binaxis = "y", stackdir = "center", position = "dodge", alpha=0.6)+
  labs(x="", y="epiCHAOS")+
  scale_color_manual(values = rev(c("steelblue4", "black")))+
  scale_fill_manual(values = rev(c("steelblue4", "black")))+
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

barplot3 <- ggplot(het, aes(x = reorder(state, mean.het), y = mean.het, fill=highlight, color=highlight)) +
  geom_bar(stat="identity", position = "dodge", alpha=0.8, width = 0.7)+
  labs(x="", y="epiCHAOS")+
  scale_color_manual(values = rev(c("steelblue4", "black")))+
  scale_fill_manual(values = rev(c("steelblue4", "black")))+
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


for (i in unique(het$state)) {
  atac@meta.data$epiCHAOS[atac@meta.data$refined_annotation==i] <- het$mean.het[het$state==i]
}

temp <- data.frame(umap2=atac@reductions$umap@cell.embeddings) %>% merge(atac@meta.data, by=0)
colnames(temp)[2:3] <- c("umap1", "umap2")

# create umap
umap3 <- ggplot(temp, aes(x=umap1, y=umap2, color=epiCHAOS)) +
  geom_point(size=0.1, alpha=0.7)+
  scale_color_distiller(palette = "Blues", direction = 1)+
  labs( x="UMAP1", y="UMAP2") +
  theme_classic()

#--- plot epiCHAOS vs pseudotime
for (i in unique(temp$refined_annotation)) {
  temp$pseudotime.avg[temp$refined_annotation==i] <- mean(temp$NNv1_age[temp$refined_annotation==i])
}

temp$highlight <- ifelse(temp$refined_annotation %in% c("Undifferentiated", "Germ cell", "Blastoderm", "Unknown"), T, F)
temp$label <- ifelse(temp$highlight==T|temp$refined_annotation %in% c("Neural", "PNS & sense"), temp$refined_annotation, "")
temp <- unique(temp[,c("epiCHAOS", "pseudotime.avg", "refined_annotation", "highlight", "label")])
fit <- lm(temp$epiCHAOS~temp$pseudotime.avg)

peudotime3 <- ggplot(temp, aes(x=pseudotime.avg, y=epiCHAOS, label=label)) +
  geom_point(size=3, alpha=0.5, aes(color=highlight))+
  scale_color_manual(values = c("black", "steelblue4"))+
  #ggrepel::geom_text_repel(aes(color=highlight), size=3)+
  geom_smooth(color="grey30", method="lm", se=F)+
  lims(y=c(0, 1))+
  #stat_cor()+
  labs( x="Pseudotime", y="epiCHAOS", subtitle = paste0("R^2 = ", signif(summary(fit)$r.squared, 3))) +
  theme_classic()


setwd("/omics/groups/OE0219/internal/KatherineK/ATACseq/epiCHAOS-Figures/Figure 2/")

pdf("dotplots_combined.pdf", 5, 10)
ggpubr::ggarrange(dotplot3, dotplot2, dotplot1, legend = F,  align = "hv", nrow=3)
dev.off()

pdf("barplots_combined.pdf", 5.5, 10)
ggpubr::ggarrange(barplot3, barplot2, barplot1, legend = F,  align = "hv", nrow=3)
dev.off()

pdf("umaps_combined.pdf", 10, 3)
ggpubr::ggarrange(umap1, umap2, umap3, legend = F,  align = "h", ncol=3)
dev.off()


pdf("pseudotime_combined_unlabelled.pdf", 12, 3.5)
ggpubr::ggarrange(pseudotime1,pseudotime2, peudotime3, legend = F, ncol=3)
dev.off()

#--- .svg images

svg("epiCHAOS_vs_devtime.svg", 18,4)
ggarrange(pseudotime1, pseudotime2, peudotime3, ncol = 3)
dev.off()

svg("epiCHAOS_barplots_development.svg", 20,4)
ggarrange(barplot1, barplot2, barplot3, ncol = 3, align = "h", widths = c(0.75, 1, 1), legend=F)
dev.off()

svg("umaps_combined.svg", 10, 3)
ggarrange(umap1, umap2, umap3, legend = F,  align = "h", ncol=3)
dev.off()



