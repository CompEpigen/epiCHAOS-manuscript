
#--- test epiCHAOS in three different developmental contexts and create sorted barplots (related to Figure 2)
#--- (i) hematopoietic system (Granja et al.)
#--- (ii) mouse gastrulation (Argelaguet et al.)
#--- (iii) drosophila embryogenesis (Calderon et al.)

#------------------------------------------------------------------------------

set.seed(10)

#--- load packages
library(magrittr)
library(stringr)
library(ggplot2)
library(ggbeeswarm)
library(ggpubr)
library(epiCHAOS)

source("~/Scripts/miscellaneous/myfunctions.R")

# 
# #--- (i) hematopoietic system
# 
# #--- load scATAC from normal hematopoiesis from Granja et al. 2019, subset for bone marrow (BM)-derived
# atac <- readRDS("/omics/groups/OE0219/internal/KatherineK/data/scATAC/Granja2019/scATAC-Healthy-Hematopoiesis-191120.rds")
# atac <- atac[,grepl("BMMC", atac@colData$Group)]
# 
# #--- counts matrix
# hema <- atac@assays$data$counts # %>% as.matrix()
# 
# #--- create matrices per cell type for epiCHAOS calculation
# datasets <- list()
# for (i in unique(atac$BioClassification[!grepl("Unk", atac$BioClassification)])) { # exclude cells annotated as "Unknown"
#   
#   cellnames <- colnames(atac)[atac$BioClassification==i] %>% intersect(colnames(hema))
#   ncells <- min(c(100, length(cellnames))) # select n cells form each celltype
#                 
#   #--- subsample five times x 100 cells from each celltype
#   for (n in 1:5) {
# 
#     temp <- hema[,sample(cellnames, ncells)]
#     temp[temp>1] <- 1 # binarise data
#     datasets[[paste0(i, "-", n)]] <- temp
#   }
# }
# 
# lapply(datasets, dim)
# 
# #--- compute and save epiCHAOS scores
# het <- compute_eITH(datasets)
#  
# saveRDS(het, "/omics/groups/OE0219/internal/KatherineK/ATACseq/Hematopoietic/epiCHAOS_scores_allpeaks_subsampling.Rds")
# 
#--- if loading downstream
het <- readRDS("/omics/groups/OE0219/internal/KatherineK/ATACseq/Hematopoietic/epiCHAOS_scores_allpeaks_subsampling.Rds")

#--- adjust names and highlight progenitor cells for plotting
het$state <- het$state %>% str_split("_") %>% lapply("[", 2) %>% str_replace_all("\\.", " ") %>% str_split("-") %>% lapply("[", 1) %>% unlist()
het$progenitor <- ifelse(het$state %in% c("HSC", "GMP", "GMP Neut", "CLP 1", "CLP 2", "CMP LMPP", "Early Eryth"), T, F)

violin1 <- ggplot(het, aes(x = reorder(state, het.adj), y = het.adj, fill=progenitor, color=progenitor)) +
  geom_violin()+
  labs(x="", y="epiCHAOS")+
  scale_color_manual(values = rev(c("steelblue4", "black")))+
  scale_fill_manual(values = rev(c("steelblue4", "black")))+
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

# #--- barplot (if only one measurement per condition)
# barplot1 <- ggplot(het, aes(x = reorder(state, mean.het), y = mean.het, fill=progenitor, color=progenitor)) +
#   geom_bar(stat="identity", position = "dodge", alpha=0.8, width = 0.7)+
#   labs(x="", y="epiCHAOS")+
#   scale_color_manual(values = rev(c("steelblue4", "black")))+
#   scale_fill_manual(values = rev(c("steelblue4", "black")))+
#   theme_classic() +
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

#--- reload atac data to get umap coordinates
atac <- readRDS("/omics/groups/OE0219/internal/KatherineK/data/scATAC/Granja2019/scATAC-Healthy-Hematopoiesis-191120.rds")
atac@colData$celltype <- atac@colData$BioClassification %>% str_split("_") %>% lapply("[", 2) %>% str_replace_all("\\.", " ")

#--- add "mean.het" column to epiCHAOS result, computing mean measurements after subsampling
het <- het %>% group_by(state) %>% summarise(mean.het=mean(het.adj)) %>% data.frame()

#--- merge previous cell annotations including umap coordinates with epiCHAOS scores
temp <- data.frame(umap1=atac@colData$UMAP1, umap2=atac@colData$UMAP2, celltype=atac@colData$celltype) %>% merge(het[,c("state", "mean.het")], by.x="celltype", by.y="state")
temp$epiCHAOS <- temp$mean.het
umap1 <- ggplot(temp, aes(y=umap2, x=umap1, color=epiCHAOS)) +
  geom_point(size=0.1, alpha=0.7)+
  scale_color_distiller(palette = "Blues", direction = 1)+
  labs( x="UMAP1", y="UMAP2") +
  theme_classic()

#--- correlate with cytoTRACE scores previously computed
cytotrace <- readRDS("/omics/groups/OE0219/internal/KatherineK/ATACseq/eICH-vs-tICH/cytotrace_Granja_2019_hematopoietic_BM.Rds")
cytotrace$CytoTRACErank

#--- load scRNA-seq data from which cytoTRACE scores were computed, then assign a mean cytoTRACE score per celltype so they can be compared with epiCHAOS scores computed on cell-type level
rna <- readRDS("/omics/groups/OE0219/internal/KatherineK/data/scATAC/Granja2019/scRNA-Healthy-Hematopoiesis-191120.rds")
rna <- rna[,grepl("BMMC", rna$Group)]
rna@colData$cytotrace <- cytotrace$CytoTRACErank
cyto <- data.frame(rna@colData)

#--- mean cytoTRACE score per cell-type
for (i in unique(cyto$BioClassification)) {
  cyto$mean.cytotrace[cyto$BioClassification==i] <- median(cyto$cytotrace[cyto$BioClassification==i])
}

#--- adjust dataframe labels etc. for plotting
cyto <- cyto[,c("BioClassification", "mean.cytotrace")] %>% unique()
cyto$BioClassification <- cyto$BioClassification %>% str_split("_") %>% lapply("[", 2) %>% unlist() %>% str_replace_all("\\.", " ")
cyto <- merge(cyto, het, by.x="BioClassification", by.y="state")
cyto$progenitor <- ifelse(cyto$BioClassification %in% c("HSC", "GMP", "GMP Neut", "CLP 1", "CLP 2", "CMP LMPP", "Early Eryth"), T, F)
cyto$label <- cyto$BioClassification
fit <- lm(cyto$mean.het~cyto$mean.cytotrace)
cyto$label[cyto$progenitor==F] <- ""

#--- scatter plot of CytoTRACE scores vs epiCHAOS scores
pseudotime1 <- ggplot(cyto, aes(x=mean.cytotrace, y=mean.het, label=label)) +
  geom_point(size=3, alpha=0.7, aes(color=progenitor))+
  scale_color_manual(values = c("black", "steelblue4"))+
  stat_cor()+
  lims(y=c(0, 1)) +
  geom_smooth(color="grey30", method = "lm", se=F)+
  #ggrepel::geom_text_repel(aes(color=progenitor), size=3, max.overlaps = 10)+
  labs( x="CytoTRACE", y="epiCHAOS", subtitle = paste0("R^2 = ", signif(summary(fit)$r.squared, 3))) +
  theme_classic()

# 
# # --- (ii) mouse gastrulation
#
# #--- read scATAC data from mouse gastrulation
# atac <- readRDS("/omics/groups/OE0219/internal/KatherineK/data/scATAC/mouse-gastrulation/Save-ArchR-Project.rds")
# mat <- readRDS("/omics/groups/OE0219/internal/KatherineK/data/scATAC/mouse-gastrulation/Matrices/PeakMatrix_summarized_experiment.rds")
# 
# #--- peaks-by-cells matrix
# mat <- mat@assays@data$PeakMatrix
# 
# atac$group <- atac$celltype.mapped
# 
# #--- create peaks-cy-cells matrices for each cell/tissue type
# datasets <- list()
# for (celltype in unique(na.omit(atac$group))) {
#   print(celltype)
#   
#   ids <- atac$cellNames[!is.na(atac$group) & atac$group==celltype & atac$genotype=="WT"] %>% intersect(colnames(mat))
#   if (length(ids)<30) { next } # in case there are very small groups, skip them
#   
#   #--- number of cells to subsample
#   ncells <- min(c(100, length(ids)))
#   
#   for (n in 1:5) {
#     temp <- mat[,sample(ids, ncells)]
#     temp[temp>1] <- 1
#     datasets[[paste0(celltype, "-", n)]] <- temp
#   }
# 
# }
# 
# lapply(datasets, dim)
# length(datasets)
# 
# #--- compute epiCHAOS scores and save
# het <- compute_eITH(datasets)
# 
# saveRDS(het, "/omics/groups/OE0219/internal/KatherineK/ATACseq/embryogenesis/epiCHAOS_scores_allpeaks_subsampling.Rds")
# 
het <- readRDS("/omics/groups/OE0219/internal/KatherineK/ATACseq/embryogenesis/epiCHAOS_scores_allpeaks_subsampling.Rds")

#--- adjust labels and highlight undifferentiated celltypes for plotting
het$state <- het$state %>% str_replace_all("_", " ") %>% str_split("-") %>% lapply("[", 1) %>% unlist()
het$group <- ifelse(grepl("Primitive|Epiblast|PGC", het$state), "primitive", "other")

violin2 <- ggplot(het, aes(x = reorder(state, het.adj), y = het.adj, fill=group, color=group)) +
  geom_violin()+
  labs(x="", y="epiCHAOS")+
  scale_color_manual(values = rev(c("steelblue4", "black")))+
  scale_fill_manual(values = rev(c("steelblue4", "black")))+
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

# #--- create barplot
# barplot2 <- ggplot(het, aes(x = reorder(state, mean.het), y = mean.het, fill=group, color=group)) +
#   geom_bar(stat="identity", position = "dodge", alpha=0.8, width = 0.7)+
#   labs(x="", y="epiCHAOS")+
#   scale_color_manual(values = rev(c("steelblue4", "black")))+
#   scale_fill_manual(values = rev(c("steelblue4", "black")))+
#   theme_classic() +
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

#--- get mean of epiCHAOS scores per group
het <- het %>% group_by(state) %>% summarise(mean.het=mean(het.adj)) %>% data.frame()
het$group <- ifelse(grepl("Primitive|Epiblast|PGC", het$state), "primitive", "other")


#--- plot umap using reference cells from https://www.nature.com/articles/s41586-019-0933-9#MOESM3
ref <- read.csv("/omics/groups/OE0219/internal/KatherineK/data/scATAC/mouse-gastrulation/reference_cells_metadata.csv")

ref$celltype <- ref$celltype %>% str_replace_all("/", " ")

for (i in unique(ref$celltype)) {
  if (!i %in% het$state) { print(i); next }
  ref$het[ref$celltype==i] <- unique(het$mean.het[het$state==i])
}

ref <- ref[!is.na(ref$het),]

#--- plot umap
umap2 <- ggplot(ref[!is.na(ref$het),], aes(x=umapX, y=umapY, color=het)) +
  geom_point(size=0.1, alpha=0.7)+
  scale_color_distiller(palette = "Blues", direction = 1)+
  labs( x="UMAP1", y="UMAP2") +
  theme_classic()


#--- correlate epiCHAOS with developmental stage
temp <- ref[,c("stage", "het", "celltype")] %>% unique()
temp$stage <- temp$stage %>% str_remove("E") %>% as.numeric()
temp <- na.omit(temp)

#--- assign average developmental time per cell/tissue type so taht it can be correlated with epiCHAOS scores
for (i in unique(temp$celltype)) {
  temp$average.stage[temp$celltype==i] <- mean(temp$stage[temp$celltype==i])
}

temp$stage <- NULL
temp <- unique(temp)

#--- highlight least differentiated tissues
temp$highlight <- ifelse(temp$celltype %in% c("Anterior Primitive Streak", "Primitive Streak", "Epiblast", "Nascent mesoderm"), T, F)
temp$label <- ifelse(temp$highlight==T|temp$celltype %in% c("Nascent mesoderm"), temp$celltype, "")
fit <- lm(temp$het~temp$average.stage)

#--- plot epiCHAOS scores vs developmental time
pseudotime2 <- ggplot(temp, aes(x=average.stage, y=het, label=label)) +
  geom_point(size=3, alpha=0.5, aes(color=highlight))+
  scale_color_manual(values = c("black", "steelblue4"))+
  #ggrepel::geom_text_repel(size=3, aes(color=highlight))+
  geom_smooth(color="grey30", method="lm", se = F)+
  stat_cor()+
  lims(y=c(0, 1))+
  labs( x="Developmental stage", y="epiCHAOS", subtitle = paste0("R^2 = ", signif(summary(fit)$r.squared, 3))) +
  theme_classic()

# #--- (iii) drosophila embryogenesis
# 
# #--- read in atac seurat object
# atac <- readRDS("/omics/groups/OE0219/internal/KatherineK/data/scATAC/embryogenesis/fly.all.downsampled_seurat_filtered_processed.rds")
# 
# #--- load metadata including celltype/tissue annotation
# temp <- readRDS("/omics/groups/OE0219/internal/KatherineK/data/scATAC/embryogenesis/atac_meta.rds")
# temp <- temp[,c("NNv1_age", "NNv1_time.new", "refined_annotation", "seurat_clusters")]
# 
# #--- merge metadata columns
# atac@meta.data <- atac@meta.data %>% merge(temp, by=0) %>% keepRow()
# atac$group <- atac$refined_annotation
# samples <- atac$group %>% unique()
# 
# #--- update annotations
# atac@meta.data$refined_annotation <- as.character(atac@meta.data$refined_annotation) #
# atac@meta.data$refined_annotation[atac@meta.data$refined_annotation=="Unknown"] <- "Undifferentiated"
# 
# # create a peaks-by-cells matrices for each celltype
# datasets <- list()
# for (i in unique(atac@meta.data$group)) {
#   print(i)
#   cells.select <- atac@meta.data[atac@meta.data$group==i,] %>% rownames() %>% intersect(colnames(atac))
#   
#   if (length(cells.select)<30) { next }
#   
#   # number of cells for subsampling
#   ncells <- min(c(100, length(cells.select)))
#   
#   for (n in 1:5) {
#     temp <- atac@assays$RNA@data[,sample(cells.select, ncells)]
#     temp[temp>1] <- 1
#     datasets[[paste0(i, "-", n)]] <- temp
#   }
#   
# }
# 
# lapply(datasets, dim)
# 
# #--- compute epiCHAOS scores and save
# het <- compute_eITH(datasets)
# 
# saveRDS(het, "/omics/groups/OE0219/internal/KatherineK/ATACseq/mouse-gastrulation/epiCHAOS_scores_allpeaks_subsampling.Rds")
# 

het <- readRDS("/omics/groups/OE0219/internal/KatherineK/ATACseq/mouse-gastrulation/epiCHAOS_scores_allpeaks_subsampling.Rds")
het$state <- het$state %>% str_split("-") %>% lapply("[", 1) %>% unlist()

#--- noted in Caldeeron et al. manuscript that this "unknown" cluster comprises undifferentiated cells and maternal cells, so rename it
het$state[het$state=="Unknown"] <- atac@meta.data$refined_annotation[atac@meta.data$refined_annotation=="Unknown"] <- "Undifferentiated"

#--- highlight undifferentiated etc. cells for plotting
het$highlight <- ifelse(het$state %in% c("Undifferentiated", "Germ cell", "Blastoderm"), T, F)

violin3 <- ggplot(het, aes(x = reorder(state, het.adj), y = het.adj, fill=highlight, color=highlight)) +
  geom_violin()+
  labs(x="", y="epiCHAOS")+
  scale_color_manual(values = rev(c("steelblue4", "black")))+
  scale_fill_manual(values = rev(c("steelblue4", "black")))+
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

# #--- create barplot
# barplot3 <- ggplot(het, aes(x = reorder(state, mean.het), y = mean.het, fill=highlight, color=highlight)) +
#   geom_bar(stat="identity", position = "dodge", alpha=0.8, width = 0.7)+
#   labs(x="", y="epiCHAOS")+
#   scale_color_manual(values = rev(c("steelblue4", "black")))+
#   scale_fill_manual(values = rev(c("steelblue4", "black")))+
#   theme_classic() +
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


het <- het %>% group_by(state) %>% summarise(mean.het=mean(het.adj)) %>% data.frame()
het$highlight <- ifelse(het$state %in% c("Undifferentiated", "Germ cell", "Blastoderm"), T, F)

#--- add epiCHAOS scores to metadata
for (i in unique(het$state)) {
  atac@meta.data$epiCHAOS[atac@meta.data$refined_annotation==i] <- het$mean.het[het$state==i]
}

temp <- data.frame(umap2=atac@reductions$umap@cell.embeddings) %>% merge(atac@meta.data, by=0)
colnames(temp)[2:3] <- c("umap1", "umap2")

#--- create umap
umap3 <- ggplot(temp, aes(x=umap1, y=umap2, color=epiCHAOS)) +
  geom_point(size=0.1, alpha=0.7)+
  scale_color_distiller(palette = "Blues", direction = 1)+
  labs( x="UMAP1", y="UMAP2") +
  theme_classic()

#--- get average developmental time per cell/tissue type in order to correlate with epiCHAOS scores
for (i in unique(temp$refined_annotation)) {
  temp$pseudotime.avg[temp$refined_annotation==i] <- mean(temp$NNv1_age[temp$refined_annotation==i])
}

#--- highlight undifferentiated cells etc. for plotting
temp$highlight <- ifelse(temp$refined_annotation %in% c("Undifferentiated", "Germ cell", "Blastoderm", "Unknown"), T, F)
temp$label <- ifelse(temp$highlight==T|temp$refined_annotation %in% c("Neural", "PNS & sense"), temp$refined_annotation, "")
temp <- unique(temp[,c("epiCHAOS", "pseudotime.avg", "refined_annotation", "highlight", "label")])

#--- fit regression model
fit <- lm(temp$epiCHAOS~temp$pseudotime.avg)

#--- scatterplot of epiCHAOS scores vs developmental time
peudotime3 <- ggplot(temp, aes(x=pseudotime.avg, y=epiCHAOS, label=label)) +
  geom_point(size=3, alpha=0.5, aes(color=highlight))+
  scale_color_manual(values = c("black", "steelblue4"))+
  #ggrepel::geom_text_repel(aes(color=highlight), size=3)+
  geom_smooth(color="grey30", method="lm", se=F)+
  lims(y=c(0, 1))+
  stat_cor()+
  labs( x="Pseudotime", y="epiCHAOS", subtitle = paste0("R^2 = ", signif(summary(fit)$r.squared, 3))) +
  theme_classic()


#--- save plots related to Figure 2
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

svg("epiCHAOS_violinplots_development.svg", 20,4)
ggarrange(violin1, violin2, violin3, ncol = 3, align = "h", widths = c(0.7, 1, 1), legend=F)
dev.off()

svg("umaps_combined.svg", 10, 3)
ggarrange(umap1, umap2, umap3, legend = F,  align = "h", ncol=3)
dev.off()


#--- write to csv for suplementary information
het <- readRDS("/omics/groups/OE0219/internal/KatherineK/ATACseq/Hematopoietic/epiCHAOS_scores_allpeaks_subsampling.Rds")
write.csv(het, "/omics/groups/OE0219/internal/KatherineK/ATACseq/epiCHAOS-supplementary-data/Revision/epiCHAOS_hematopoiesis_subsampling.csv")

het <- readRDS("/omics/groups/OE0219/internal/KatherineK/ATACseq/embryogenesis/epiCHAOS_scores_allpeaks_subsampling.Rds")
het$state <- het$state %>% str_split("-") %>% lapply("[", 1) %>% unlist()
write.csv(het, "/omics/groups/OE0219/internal/KatherineK/ATACseq/epiCHAOS-supplementary-data/Revision/epiCHAOS_mouse_gastrulation_subsampling.csv")

het <- readRDS("/omics/groups/OE0219/internal/KatherineK/ATACseq/mouse-gastrulation/epiCHAOS_scores_allpeaks_subsampling.Rds")
het$state <- het$state %>% str_split("-") %>% lapply("[", 1) %>% unlist()
write.csv(het, "/omics/groups/OE0219/internal/KatherineK/ATACseq/epiCHAOS-supplementary-data/Revision/epiCHAOS_drosophila_embryo_subsampling.csv")
