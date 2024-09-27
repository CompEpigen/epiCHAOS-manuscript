

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

#--- get cluster annotation from metadata
meta <- lica@cellColData[,c("Sample", "Clusters")] %>% data.frame()

#--- msigdbr gene sets
gene_sets = msigdbr(species = "Homo sapiens")
head(gene_sets)

set.names <- gene_sets$gs_name %>% unique()
length(set.names)

#--- malignant cell clusters
clusters <- unique(lica$Clusters)

#--- compute per cluster scores per gene set
gene.set.scores <- data.frame(row.names = clusters)

#--- for each gene set, retrieve the list of genes within that set, and compute the average gene scores and gene-set scores for all cells within each cluster
for (set in set.names) {
  
  print(set)
  
  #--- retrieve genes in the selected gene set, of which correspond to rownames of gene score matrix
  genes <- gene_sets$gene_symbol[gene_sets$gs_name==set] %>% intersect(rownames(mat.ge)) 
  
  #--- exclude sets with fewer than 10 genes
  if (length(genes)<10) { next } 
  
  #--- subset for genes in set of interest
  mat.set <- mat.ge[rownames(mat.ge) %in% genes,]
  
  #--- transpose subset gene scores matrix and coerce to dataframe
  temp <- t(mat.set) %>% as.matrix() %>% data.frame() 
  
  #--- merge gene scores and cluster annotation
  temp <- merge(temp, meta, by=0) %>% keepRow() 
  
  #--- compute average gene scores for all cells within each
  gex.per.cluster <- data.frame()
  for (i in clusters) {
    gex.per.cluster <- rbind(gex.per.cluster, colMeans(temp[temp$Clusters==i,colnames(temp) %notin% c("Sample", "Clusters")]))
  }
  
  #--- set rownames to cluster names and colnames to gene names
  rownames(gex.per.cluster) <- clusters
  colnames(gex.per.cluster) <- colnames(temp[,colnames(temp) %notin% c("Sample", "Clusters")])
  
  #--- compute average of gene scores for each cluster to produce the gene set scores per cluster
  gene.set.scores[,set] <- rowMeans(gex.per.cluster)
  
}

#--- load previously computed epiCHAOS scores
het <- readRDS("epiCHAOS_scores_copy_corrected_all_peaks_subsampling.Rds")
het$state <- het$state %>% str_split("-") %>% lapply("[", 1) %>% unlist()
het <- het %>% group_by(state) %>% summarise(mean.het=mean(het)) %>% data.frame()

gene.set.scores <- readRDS("gene_set_scores_vs_het_allgenesets.Rds")
gene.set.scores <- gene.set.scores$scores
gene.set.scores$mean.het <- NULL

#--- merge gene set scores and epiCHAOS scores
temp <- merge(gene.set.scores, het[,c("state", "mean.het")], by.x=0, by.y="state") %>% unique()
rownames(temp) <- temp$Row.names
temp$Row.names <- NULL

#--- compute correlations between each gene set score and epiCHAOS score
require(Hmisc) # use rcorr so that you can also get the p-values
cor <- rcorr(as.matrix(temp))
pvals <- cor$P
cor <- cor$r
cor["mean.het",] %>% sort()

#--- write to supplementary info
data.frame(r=cor[,"mean.het"], pval=pvals[,"mean.het"]) %>% arrange(desc(r)) %>% write.csv("/omics/groups/OE0219/internal/KatherineK/ATACseq/epiCHAOS-supplementary-data/Revision/correlations_and_pvals_lica_epiCHAOS_vs_gene_sets.csv")


#--- save list contianing the gene set scores and their correlations with epiCHAOS (now saved result after subsampling)
obj <- list(scores=temp, correlation=cor)
saveRDS(obj, file = "gene_set_scores_vs_het_allgenesets.Rds")

#--- downstream analysis
obj <- readRDS("gene_set_scores_vs_het_allgenesets.Rds")

#--- scores and correlations
scores <- obj$scores
cor <- obj$correlation
cor <- cor[,"mean.het"]

#--- for plotting
temp <- data.frame(cor)
colnames(temp) <- "pearson"
temp$n <- 1:nrow(temp)
temp <- temp[grepl("HALLMARK", rownames(temp)),]
labels <- c("HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION")
temp$label <- ifelse(rownames(temp) %in% labels, rownames(temp), "")
temp$label <- temp$label %>% str_replace_all("_", " ")

#--- Figure 3B-C. plot of correlations
elbow.lica <- ggplot(temp, aes(y=pearson, x=reorder(n, pearson), label=label)) +
  geom_point(size=1, alpha=0.7)+
  labs(x="", y="Correlation vs epiCHAOS") +
  #ggrepel::geom_label_repel(max.overlaps = 10000, size=3)+
  theme_classic() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

#--- check top 5
temp %>% arrange(pearson) %>% tail(5) %>% rownames() %>% rev()
# [1] "HALLMARK_KRAS_SIGNALING_UP"                 "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION"
# [3] "HALLMARK_UV_RESPONSE_DN"                    "HALLMARK_COAGULATION"                      
# [5] "HALLMARK_PANCREAS_BETA_CELLS" 

#--- correlation plot of interesting gene sets in liver cancer
interest <- c("MEBARKI_HCC_PROGENITOR_WNT_UP", "CHIANG_LIVER_CANCER_SUBCLASS_CTNNB1_UP", "HOSHIDA_LIVER_CANCER_SURVIVAL_UP")

cor <- obj$correlation
M <- cor[interest,"mean.het"] %>% as.matrix()
colnames(M) <- ""
rownames(M) <- c("HCC Progenitor WNT Up (Mebarki et al.)", "Liver Cancer Beta-Catenin subclass (Chiang et al.)", "Liver Cancer Survival (Hoshida et al.)")

#--- save corrplot for genes of interest (supplementary figure)
pdf("/omics/groups/OE0219/internal/KatherineK/ATACseq/epiCHAOS-Figures/Figure 3/corrplot_lica.pdf", 9, 2.5)
corrplot(M, method = 'ellipse',  addCoef.col = 'black',col = rev(COL2('RdBu', 10)), cl.pos="n", tl.col="black", tl.pos = "right", title = "pearson R")
dev.off()


#--- brca

setwd("/omics/groups/OE0219/internal/KatherineK/ATACseq/breast-cancer/epithelial/")


#--- load previously computed epiCHAOS scores
het <- readRDS("epiCHAOS_scores_count_corrected_allpeaks_subsampling.Rds")
het$state <- het$state %>% str_split("-") %>% lapply("[", 1) %>% unlist()
het <- het %>% group_by(state) %>% summarise(mean.het=mean(het)) %>% data.frame()

gene.set.scores <- readRDS("gene_set_scores_vs_het_allgenesets_epithelial_clusters.Rds")
gene.set.scores <- gene.set.scores$scores
gene.set.scores$mean.het <- NULL

#--- merge gene set scores and epiCHAOS scores
temp <- merge(gene.set.scores, het[,c("state", "mean.het")], by.x=0, by.y="state") %>% unique()
rownames(temp) <- temp$Row.names
temp$Row.names <- NULL

#--- compute correlations between each gene set score and epiCHAOS score
require(Hmisc) # use rcorr so that you can also get the p-values
cor <- rcorr(as.matrix(temp))
pvals <- cor$P
cor <- cor$r
cor["mean.het",] %>% sort()

#--- write to supplementary info
data.frame(r=cor[,"mean.het"], pval=pvals[,"mean.het"]) %>% arrange(desc(r)) %>% head() 
data.frame(r=cor[,"mean.het"], pval=pvals[,"mean.het"]) %>% arrange(desc(r)) %>% write.csv("/omics/groups/OE0219/internal/KatherineK/ATACseq/epiCHAOS-supplementary-data/Revision/correlations_and_pvals_brca_epiCHAOS_vs_gene_sets.csv")


#--- save list contianing the gene set scores and their correlations with epiCHAOS (now saved result after subsampling)
obj <- list(scores=temp, correlation=cor)
saveRDS(obj, file = "gene_set_scores_vs_het_allgenesets.Rds")

#--- downstream analysis
obj <- readRDS("gene_set_scores_vs_het_allgenesets.Rds")

#--- scores and correlations
scores <- obj$scores
cor <- obj$correlation
cor <- cor[,"mean.het"]

#--- for plotting
temp <- data.frame(cor)
colnames(temp) <- "pearson"
temp$n <- 1:nrow(temp)
temp <- temp[grepl("HALLMARK", rownames(temp)),]
labels <- c("HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION")
temp$label <- ifelse(rownames(temp) %in% labels, rownames(temp), "")
temp$label <- temp$label %>% str_replace_all("_", " ")

#--- Figure 3B-C. plot of correlations
elbow.brca <- ggplot(temp, aes(y=pearson, x=reorder(n, pearson), label=label)) +
  geom_point(size=1, alpha=0.7)+
  labs(x="", y="Correlation vs epiCHAOS") +
  #ggrepel::geom_label_repel(max.overlaps = 10000, size=3)+
  theme_classic() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

#--- check top 5
temp %>% arrange(pearson) %>% tail(5) %>% rownames() %>% rev()
# "HALLMARK_ALLOGRAFT_REJECTION"               "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION" "HALLMARK_INFLAMMATORY_RESPONSE"            
# [4] "HALLMARK_COMPLEMENT"                        "HALLMARK_UV_RESPONSE_DN"

#--- correlation plot of interesting gene sets in liver cancer
interest <- c("NIKOLSKY_OVERCONNECTED_IN_BREAST_CANCER", "SCHUETZ_BREAST_CANCER_DUCTAL_INVASIVE_UP", "LIEN_BREAST_CARCINOMA_METAPLASTIC", "JECHLINGER_EPITHELIAL_TO_MESENCHYMAL_TRANSITION_UP", "GOBP_TRANSFORMING_GROWTH_FACTOR_BETA_ACTIVATION", "GOBP_CANONICAL_WNT_SIGNALING_PATHWAY")


cor <- obj$correlation
M <- cor[interest,"mean.het"] %>% as.matrix()
colnames(M) <- ""
rownames(M) <- c("Overconnected in Breast Cancer (Nikolsky et al.)", "Ductal Invasive Breast Cancer (Schuetz et al.)", "Metaplastic Breast Cancer (Lien et al.)",  "EMT in Breast Cancer (Jechlinger et al.)",  "TGF-beta activation", "Canonical WNT Signaling")

#--- save corrplot for genes of interest (supplementary figure)
pdf("/omics/groups/OE0219/internal/KatherineK/ATACseq/epiCHAOS-Figures/Figure 3/corrplot_brca.pdf", 10, 4)
corrplot(M, method = 'ellipse',  addCoef.col = 'black',col = rev(COL2('RdBu', 10)), cl.pos="n", tl.col="black", tl.pos = "right", title = "pearson R")
dev.off()


#--- elbow plots for each tumor type hallmark gene set correlations
svg("/omics/groups/OE0219/internal/KatherineK/ATACseq/epiCHAOS-Figures/Figure 3/elbowplots_brca_lica.svg", 4, 3)
ggarrange(elbow.brca, elbow.lica, ncol=2)
dev.off()

