

#--- test epiCHAOS for comparison of old vs young HSC

setwd("/omics/groups/OE0219/internal/KatherineK/ATACseq")

# library(devtools)
# install_github("Katherine-Kelly/epiCHAOS-R")

#--- load libraries
library(magrittr)
library(stringr)
library(ggplot2)
library(ggpubr)
library(epiCHAOS)
library(TxDb.Mmusculus.UCSC.mm9.knownGene)


#--- load data from young (2 independent experiments) vs old (3 independent experiments) vwf+ HSCc
atac <- read.table("/omics/groups/OE0219/internal/KatherineK/data/scATAC/Meng2023/GSE219096_SingleCell_ATACseq_Young_vs_Old_HSC_counts.txt", header = T, sep=" ")
atac[1:4,1:6] 
dim(atac)
atac[atac > 0] <- 1  # make binary

#--- create metadata table
grouping <- colnames(atac) %>% str_split("\\.") %>% lapply("[", 1) %>% unlist()
meta <- data.frame(row.names = colnames(atac), grouping=grouping)

set.seed(10)

#--- compute epiCHAOS scores
het <- epiCHAOS(counts = atac, meta = meta, n = 100, subsample = 5)
het$state <- rownames(het) %>% str_remove("group-") %>% str_split("-") %>% lapply("[", 1) %>% unlist()
het$group <- ifelse(grepl("Young", het$state), "Young", "Old")

saveRDS(het, file = "Aging/epiCHAOS_scores_subsampling.Rds")

pdf("Aging/boxplot_epiCHAOS_young_vs_old_HSCs.pdf", 3, 4)
svg("epiCHAOS-Figures/Figure 3/violinplot_epiCHAOS_young_vs_old_HSCs.svg", 4, 4)
ggplot(het, aes(y=het.adj, x=reorder(state, het.adj), color=group, fill=group)) +
  geom_violin(linewidth=0.5)+
  labs( x="", y="epiCHAOS") +
  scale_fill_manual(values = c("steelblue4", "grey20"))+
  scale_color_manual(values = c("steelblue4", "grey20"))+
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()

#--- compute epiCHAOS scores on promoters

#--- mm9 promoters
atac.gr <- rownames(atac) %>% str_replace("-", ":") %>% GRanges()
names(atac.gr) <- 1:nrow(atac)
txdb <- TxDb.Mmusculus.UCSC.mm9.knownGene
promoters <- GenomicFeatures::promoters(genes(txdb), upstream = 500, downstream = 1500)
rows.select <- subsetByOverlaps(atac.gr, promoters) %>% names() %>% as.numeric()

#--- compute epiCHAOS scores at promoters
het <- epiCHAOS(counts = atac[rows.select,], meta = meta, n = 100, subsample = 5)
het$state <- rownames(het) %>% str_remove("group-") %>% str_split("-") %>% lapply("[", 1) %>% unlist()

pdf("Aging/boxplot_epiCHAOS_young_vs_old_HSCs_promoters.pdf", 3, 4)
ggplot(het, aes(y=het.adj, x=reorder(state, het.adj))) +
  geom_boxplot(fill="honeydew3", color="black", linewidth=0.5)+
  labs( x="", y="epiCHAOS") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()
