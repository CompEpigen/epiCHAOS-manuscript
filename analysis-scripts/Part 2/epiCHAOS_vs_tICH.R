

#--- compare epiCHAOS with transcriptional heterogeneity scores in three developmental datasets (relating to Figure S3)

#--- load epiCHAOS and transcriptional heterogeneity functions
source("compute_epiCHAOS.R")
source("transcriptional_heterogeneity.R")

#--- required packages
library(Seurat)
library(magrittr)
library(dplyr)
library(stringr)
library(ggplot2)
library(msigdbr)

#--- (i) Granja et al. hematopoiesis

#--- rna
rna <- readRDS("/omics/groups/OE0219/internal/KatherineK/data/scATAC/Granja2019/scRNA-Healthy-Hematopoiesis-191120.rds")
rna <- rna[,grepl("BMMC", rna$Group)]
table(rna$Group)

#--- create a seurat object, then normalise and scale the data
gex <- CreateSeuratObject(counts = rna@assays@.xData$data$counts, meta.data = data.frame(rna@colData))
gex@meta.data

gex <- NormalizeData(gex)
gex <- ScaleData(gex, features = rownames(gex))
gex <- gex@assays$RNA$scale.data
dim(gex)

#--- create matrices per celltype, subsetting 200 cells for larger groups
datasets <- list()
for (i in unique(rna$BioClassification)) {
  
  ids <- rna@colData[rna@colData$BioClassification==i,] %>% rownames()
  datasets[[i]] <- gex[,ids] %>% as.matrix()
  if (length(ids)>200) {
    datasets[[i]] <- gex[,sample(ids, 200)] %>% as.matrix()
  }
}

lapply(datasets, dim)

#--- compute transcriptional heterogeneity
tICH <- compute.tITH(datasets)

#--- load previously computed epiCHAOS scores
eICH <- readRDS("/omics/groups/OE0219/internal/KatherineK/ATACseq/Hematopoietic/epiCHAOS_scores_allpeaks.Rds")


#--- correlate epiCHAOS and tICH
t.vs.eICH <- merge(tICH, eICH, by="state")
colnames(t.vs.eICH)[2:3] <- c("tICH", "eICH")
t.vs.eICH$HSPC <- ifelse(t.vs.eICH$state %in% c("01_HSC", "05_CMP.LMPP", "06_CLP.1" ,  "07_GMP",  "08_GMP.Neut", "15_CLP.2", "02_Early.Eryth", "04_Early.Baso", "16_Pre.B"), "HSPC", "other")
t.vs.eICH$label <- t.vs.eICH$state %>% strsplit("_") %>% lapply("[", 2)  %>% str_replace_all("\\.", " ") %>% unlist()
cor.test(t.vs.eICH$tICH, t.vs.eICH$eICH, method="spearman")
fit <- lm(t.vs.eICH$tICH~t.vs.eICH$eICH)

summary(fit)$r.squared

#--- scatterplot
scatterplot1 <- ggplot(t.vs.eICH, aes(y=eICH, x=tICH, label=label)) +
  labs( x="transcriptional heterogeneity", y="epiCHAOS") +
  geom_point(alpha=0.6, size=3, aes(color=HSPC))+
  ggpubr::stat_cor()+
  scale_color_manual(values = rev(c("grey10", "steelblue4")))+
  ggrepel::geom_text_repel(size=3, aes(color=HSPC))+
  geom_smooth(color="grey10", se = F, method = "lm")+
  theme_classic()

svg("/omics/groups/OE0219/internal/KatherineK/ATACseq/epiCHAOS-Figures/Supplementary/S2/scatterplot_epiCHAOS_vs_tICH_hematopoiesis.svg", 5, 3.5)
scatterplot1
dev.off()


#--- (ii) mouse gastrulation data from Argelaguet et al.

rna <- readRDS("/omics/groups/OE0219/internal/KatherineK/data/scATAC/mouse-gastrulation/seurat.rds")

head(rna@meta.data)

#--- merge metadata with celltype annotation
metadata <- read.table("/omics/groups/OE0219/internal/KatherineK/data/scATAC/mouse-gastrulation/sample_metadata_after_archR.txt", header = T, sep="\t", comment.char = "")
metadata <- merge(metadata, rna@meta.data, by.x="cell", by.y=0)
rownames(metadata) <- metadata$cell
metadata <- metadata[rownames(metadata) %in% colnames(rna),]
rna <- rna[,rownames(metadata)]
rna@meta.data <- metadata[colnames(rna),]


#--- create matrices per celltype
datasets <- list()
for (i in na.omit(unique(rna$celltype.mapped.x))) {

  ids <- rna@meta.data[!is.na(rna@meta.data$celltype.mapped.x) & rna@meta.data$celltype.mapped.x==i,] %>% rownames() # %>% str_replace_all("#|-", ".")
  datasets[[i]] <- rna@assays$RNA@scale.data[,ids] %>% as.matrix()
  if (length(ids)>200) {
    datasets[[i]] <- rna@assays$RNA@scale.data[,sample(ids, 200)] %>% as.matrix()
  }
}

lapply(datasets, dim)

#--- calculate trancriptional heterogeneity
tICH <- compute.tITH(datasets)


#--- scale data for tICH calculation
ids <- lapply(datasets, colnames) %>% unlist()
rna <- rna[,ids]

all.genes <- rownames(rna)
rna <- ScaleData(rna, features = all.genes)

gex <- rna@assays$RNA@scale.data
dim(gex)

#--- create matrices per celltype
datasets <- list()
for (i in na.omit(unique(rna$celltype.mapped.x))) {
  
  ids <- rna@meta.data[!is.na(rna@meta.data$celltype.mapped.x) & rna@meta.data$celltype.mapped.x==i,] %>% rownames() # %>% str_replace_all("#|-", ".")
  datasets[[i]] <- rna@assays$RNA@scale.data[,ids] %>% as.matrix()
  if (length(ids)>200) {
    datasets[[i]] <- rna@assays$RNA@scale.data[,sample(ids, 200)] %>% as.matrix()
  }
}

lapply(datasets, dim)

#--- calculate trancriptional heterogeneity
tICH <- compute.tITH(datasets)

#--- correlate with previously computed epiCHAOS scores
eICH <- readRDS("/omics/groups/OE0219/internal/KatherineK/ATACseq/embryogenesis/epiCHAOS_scores_allpeaks.Rds")


eICH.vs.tICH <- tICH[,c("mean.het", "state")] %>% unique() %>% merge(eICH, by="state")
colnames(eICH.vs.tICH)[2:3] <- c("transcriptional", "epigenetic")

#--- scatterplot
scatterplot2 <- ggplot(eICH.vs.tICH, aes(y=epigenetic, x=transcriptional, label=state)) +
  labs( x="transcriptional heterogeneity", y="epiCHAOS") +
  geom_point(alpha=0.6, size=3)+
  #ggpubr::stat_cor()+
  ggrepel::geom_text_repel(size=3, max.overlaps = 10)+
  geom_smooth(color="grey10", se = F, method = "lm")+
  theme_classic()

 
svg("/omics/groups/OE0219/internal/KatherineK/ATACseq/epiCHAOS-Figures/Supplementary/S2/scatterplot_epiCHAOS_vs_tICH_gastrulation.svg", 5, 3.5)
scatterplot2
dev.off()


#--- Calderon et al. drosophila embryogenesis

rna <- readRDS("/omics/groups/OE0219/internal/KatherineK/data/scATAC/embryogenesis/main.Rds")
rna <- rna[,sample(ncol(rna), 50000)] # subset cells because data is large

#--- scale data 
rna <- ScaleData(rna, features = rownames(rna))


#--- load metadata
meta <- readRDS("/omics/groups/OE0219/internal/KatherineK/data/scATAC/embryogenesis/rna_meta.rds")
meta <- merge(meta, rna@meta.data, by="cell")
head(meta)
head(rna@meta.data)

#--- normalised and scaled gene expression matrix
gex <- rna@assays$RNA$scale.data


#--- create matrices per celltype
datasets <- list()
for (i in unique(meta$manual_annot)) {
  ids <- meta$cell[meta$manual_annot==i]
  datasets[[i]] <- gex[,colnames(gex) %in% ids] %>% as.matrix()
  if (length(ids)>200) { datasets[[i]] <- datasets[[i]][,sample(colnames(datasets[[i]]), 100)] %>% as.matrix() }
}

lapply(datasets, dim)

#--- compute transcriptional heterogeneity
het <- compute.tITH(datasets)


#--- compare to epiCHAOS scores
eICH <- readRDS("/omics/groups/OE0219/internal/KatherineK/ATACseq/mouse-gastrulation/epiCHAOS_scores_allpeaks.Rds")
match.anno <- read.table("/omics/groups/OE0219/internal/KatherineK/ATACseq/embryogenesis/epiCHAOS_scores.tsv", sep="\t", header = T, fill = T)
match.anno$mean.het <- NULL
match.anno$tissue <- match.anno$tissue %>% str_replace_all("_", " ")
eICH <- merge(eICH, match.anno, by="state")


eICH$tissue <- eICH$tissue %>% str_replace_all("_", " ")

het <- het %>% merge(eICH[eICH$tissue %in% het$state,], by.x="state", by.y="tissue")

#--- scatterplot
scatterplot3 <- ggplot(het, aes(y=mean.het.x, x=mean.het.y, label=state)) +
  labs( x="transcriptional heterogeneity", y="epiCHAOS") +
  geom_point(alpha=0.6, size=3)+
  ggpubr::stat_cor()+
  ggrepel::geom_text_repel(size=3, max.overlaps = 10)+
  geom_smooth(color="grey10", se = F, method = "lm")+
  theme_classic()

svg("/omics/groups/OE0219/internal/KatherineK/ATACseq/epiCHAOS-Figures/Supplementary/S2/scatterplot_epiCHAOS_vs_tICH_embryogenesis.svg", 5, 3.5)
scatterplot3
dev.off()
