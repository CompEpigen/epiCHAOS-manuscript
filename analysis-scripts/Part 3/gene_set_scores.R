

#--- compute gene set/pathway scores per cluster using the scATAC gene score matrix generated from ArchR

library(ArchR)
library(msigdbr)
library(magrittr)
library(dplyr)
library(stringr)
library(ggpubr)
library(corrplot)


analysis.dir <- "/omics/groups/OE0219/internal/KatherineK/ATACseq/Liver-Cancer/malignant-subset/"

setwd(analysis.dir)

#--- load ArchR project
lica <- readRDS("Save-ArchR-Project.rds")

#--- load gene activity matrix from the above archr project
mat.ge <- readRDS("gene_score_matrix.Rds")
row.names <- mat.ge@elementMetadata$name
mat.ge <- mat.ge@assays@data$GeneScoreMatrix
rownames(mat.ge) <- row.names
dim(mat.ge)
mat.ge[1:10,1:10]

# get cluster annotation from metadata
meta <- lica@cellColData[,c("Sample", "Clusters")] %>% data.frame()

# msigdbr gene sets
gene_sets = msigdbr(species = "Homo sapiens")
head(gene_sets)

set.names <- gene_sets$gs_name %>% unique()
length(set.names)

# malignant cell clusters
clusters <- unique(lica$Clusters)

#--- compute per cluster scores per gene set
gene.set.scores <- data.frame(row.names = clusters)

#--- for each gene set, retrieve the list of genes within that set, 
for (set in set.names) {
  print(set)
  genes <- gene_sets$gene_symbol[gene_sets$gs_name==set] %>% intersect(rownames(mat.ge)) # retrieve genes in the selected gene set, of which correspond to rownames of gene score matrix
  if (length(genes)<10) { next } # exclude sets with fewer than 10 genes
  mat.set <- mat.ge[rownames(mat.ge) %in% genes,]
  temp <- t(mat.set) %>% as.matrix() %>% data.frame() # transpose subset gene scores matrix and coerce to dataframe
  temp <- merge(temp, meta, by=0) %>% keepRow() # merge gene scores and cluster annotation
  
  # compute average gene scores for all cells within each
  gex.per.cluster <- data.frame()
  for (i in clusters) {
    gex.per.cluster <- rbind(gex.per.cluster, colMeans(temp[temp$Clusters==i,colnames(temp) %notin% c("Sample", "Clusters")]))
  }
  
  # set rownames to cluster names and colnames to gene names
  rownames(gex.per.cluster) <- clusters
  colnames(gex.per.cluster) <- colnames(temp[,colnames(temp) %notin% c("Sample", "Clusters")])
  
  # compute average of gene scores for each cluster to produce the gene set scores per cluster
  gene.set.scores[,set] <- rowMeans(gex.per.cluster)
  
}

#--- load epigenetic heterogeneity scores
het <- readRDS("epiCHAOS_scores_copy_corrected_all_peaks.Rds")

temp <- merge(gene.set.scores, het[,c("state", "mean.het")], by.x=0, by.y="state") %>% unique()
rownames(temp) <- temp$Row.names
temp$Row.names <- NULL
cor <- cor(as.matrix(temp))
cor["mean.het",] %>% sort()

# save list contianing the pathway scores and their correlations with epiCHAOS
obj <- list(scores=temp, correlation=cor)

saveRDS(obj, file = "gene_set_scores_vs_het_allgenesets.Rds")

#--- downstream analysis
obj <- readRDS("gene_set_scores_vs_het_allgenesets.Rds")

scores <- obj$scores
cor <- obj$correlation
cor <- cor[,"mean.het"]

temp <- data.frame(cor)
colnames(temp) <- "pearson"
temp$n <- 1:nrow(temp)
temp <- temp[grepl("HALLMARK", rownames(temp)),]
labels <- c("HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION")
temp$label <- ifelse(rownames(temp) %in% labels, rownames(temp), "")
temp$label <- temp$label %>% str_replace_all("_", " ")

# plot of correlations
ggplot(temp, aes(y=pearson, x=reorder(n, pearson), label=label)) +
  geom_point(size=1, alpha=0.7)+
  labs(x="", y="Correlation vs epiCHAOS") +
  #ggrepel::geom_label_repel(max.overlaps = 10000, size=3)+
  theme_classic() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

# check top 5
temp %>% arrange(pearson) %>% tail() %>% rownames() %>% rev()


# correlation plot of interesting gene sets in liver cancer
interest <- c("MEBARKI_HCC_PROGENITOR_WNT_UP", "CHIANG_LIVER_CANCER_SUBCLASS_CTNNB1_UP", "HOSHIDA_LIVER_CANCER_SURVIVAL_UP")

cor <- obj$correlation
M <- cor[interest,"mean.het"] %>% as.matrix()
colnames(M) <- ""
rownames(M) <- c("HCC Progenitor WNT Up (Mebarki et al.)", "Liver Cancer Beta-Catenin subclass (Chiang et al.)", "Liver Cancer Survival (Hoshida et al.)")
corrplot(M, method = 'ellipse',  addCoef.col = 'black',col = rev(COL2('RdBu', 10)), cl.pos="n", tl.col="black", tl.pos = "right", title = "pearson R")


