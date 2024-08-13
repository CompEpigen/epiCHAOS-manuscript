
#--- demonstrate for example using the HSPC aging data how epiCHAOS can be used to compute statistical significance when comparing epiCHAOS scores between two groups

#--- load libraries
library(magrittr)
library(stringr)
library(ggplot2)
library(ggpubr)
library(epiCHAOS)
library(TxDb.Mmusculus.UCSC.mm9.knownGene)
library(epiCHAOS)

#--- the compute_diff_eICH() function can be used for simple comparison of two groups on e.g. all peaks, or for comparison of multiple region types between groups
# epiCHAOS::compute_diff_eICH()


#--- load data from young (2 independent experiments) vs old (3 independent experiments) vwf+ HSCc
atac <- read.table("/omics/groups/OE0219/internal/KatherineK/data/scATAC/Meng2023/GSE219096_SingleCell_ATACseq_Young_vs_Old_HSC_counts.txt", header = T, sep=" ")
atac[1:4,1:6] 
dim(atac)
atac[atac > 0] <- 1  # make binary

#--- create metadata table
grouping <- colnames(atac) %>% str_split("\\.") %>% lapply("[", 1) %>% unlist()
meta <- data.frame(row.names = colnames(atac), grouping=grouping)


group1 <- atac[,sample(grepl("Young1", colnames(atac)), 100)]
group2 <- atac[,sample(grepl("Old1", colnames(atac)), 100)]


comparison <- compute_diff_eICH(group1 = group1, group2 = group2, region.type = "all peaks", niter = 1000)

saveRDS(comparison, file = "/omics/groups/OE0219/internal/KatherineK/ATACseq/Aging/old1_vs_young1_differential.Rds")
