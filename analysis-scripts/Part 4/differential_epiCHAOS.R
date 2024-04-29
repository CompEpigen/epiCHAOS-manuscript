

#--- compute differential epigenetic heterogeneity using epiCHAOS scores at a set of genomic regions using a permutation approach

set.seed(10)

library(magrittr)
library(dplyr)
library(stringr)
library(ggplot2)
library(GenomicRanges)
library(ggpubr)
library(msigdbr)

setwd("/omics/groups/OE0219/internal/KatherineK/ATACseq/differential-eICH")


#--- scATAC from normal hematopoiesis from Granja et al. 2019, subset for bone-marrow (BM) derived cells
atac <- readRDS("/omics/groups/OE0219/internal/KatherineK/data/scATAC/Granja2019/scATAC-Healthy-Hematopoiesis-191120.rds")
atac <- atac[,atac@colData$Group %in% c("BMMC_D5T1", "BMMC_D6T1")]
counts <- atac@assays$data$counts 

#--- select sites corresponding to hg19 promoters (or silence this block, use all peaks)
require("TxDb.Hsapiens.UCSC.hg19.knownGene")
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
promoters <- promoters(genes(txdb), upstream = 1500, downstream = 500)
promoters.atac <- subsetByOverlaps(atac@rowRanges, promoters) %>% names()
dim(counts)
rownames(counts) <- names(atac@rowRanges)
counts <- counts[promoters.atac,]

#--- binarise counts
counts[counts>1] <- 1
dim(counts)

#--- select 100 cells from two groups: e.g. HSCs and monocytes, B cells or T cells
ids1 <- colnames(atac)[atac$BioClassification=="01_HSC"] %>% sample(100)
ids2 <- colnames(atac)[atac$BioClassification=="19_CD8.N"] %>% sample(100)  # 11_CD14_Mono.1, 17_B, 19_CD8.N

counts1 <- counts[,ids1]
counts2 <- counts[,ids2]

dim(counts1)
dim(counts2)

#--- differential heterogeneity at a set of genomic regions


#--- apply differentail epiCHAOSfunction to compute differential epiCHAOS per TFBS between HSCs and monocytes

#--- load encode TFBS, subset for K562 cells since we focus on the hematopoietic lineage
lola <- get(load("/omics/groups/OE0219/internal/KatherineK/Genomic_Annotations/LOLA/hg19/encode_tfbs/encode_tfbs.RData"))
index <- read.table("/omics/groups/OE0219/internal/KatherineK/Genomic_Annotations/LOLA/hg19/encode_tfbs/index.txt", header = T)
names(lola) <- index$filename
lola <- lola[grepl("K562", names(lola))]

#--- ENCODE TFBSs from K562
regions.interest <- names(lola)

#--- genomic ranges of each peak in the peaks-by-cells matrix, which we will overlap with those of the TFBSs
data.gr <- rownames(counts) %>% str_replace("_", ":") %>% str_replace("_", "-") %>% GRanges()
names(data.gr) <- 1:length(data.gr)

#--- compute differential heterogenetiy for each selected region 
gg <- list()
results <- list()
for (region in regions.interest) {
  
  state <- paste0(index$cellType[index$filename==region], " - ", index$antibody[index$filename==region])
  print(state)
  
  #--- GRanges object containing selected regions
  temp <- lola[[region]]
  
  #--- find overlaps with ATAC peaks, subset to 1000 peaks if too many are overlapping, exclude the region if too few
  select.sites <- subsetByOverlaps(data.gr, temp) %>% names() %>% as.numeric()
  if (length(select.sites)>1000) { select.sites <- sample(select.sites, 1000)}
  if (length(select.sites)<10) { next }
  
  #--- counts matrices, subset for region of interest
  group1 <- counts1[select.sites, ]
  group2 <- counts2[select.sites, ]
  
  #--- apply differential epiCHAOS function
  print("computing differentials...")
  plot.dif <- compute.diff.eICH(group1, group2, state, niter=1000)
  
  #--- save plot and results
  gg[[state]] <- plot.dif$plot
  results[[state]] <- list(mix=plot.dif$result.mix, test=plot.dif$result.test)
}

#--- save result
#saveRDS(results, "HSC-vs-CD8T/differential_epiCHAOS_TFBSs.Rds")


#--- finally calculate heterogeneity scores per region for the two groups in order to establish differences for plotting against permutation p-values
datasets <- list()
for (region in regions.interest) {
  
  state <- paste0(index$cellType[index$filename==region], " - ", index$antibody[index$filename==region])
  
  #--- GRanges object containing selected regions
  temp <- lola[[region]]
  
  #--- find overlaps with ATAC peaks, subset to 1000 peaks if too many are overlapping, exclude the region if too few
  select.sites <- subsetByOverlaps(data.gr, temp) %>% names() %>% as.numeric()
  if (length(select.sites)>1000) { select.sites <- sample(select.sites, 1000)}
  if (length(select.sites)<10) { next }
  
  datasets[[paste0("group1-", state)]] <- counts1[select.sites, ]
  datasets[[paste0("group2-", state)]] <- counts2[select.sites, ]
  
}

lapply(datasets, dim)

#--- compute epiCHAOS scores
het <- compute.eITH(datasets)

het$group <- het$state %>% str_split("-") %>% lapply("[", 1) %>% unlist()

#saveRDS(het, "HSC-vs-CD8T/epiCHAOS_scores_TFBSs.Rds")


#--- load the above saved results for pvalues computation and plotting
results <- readRDS("/omics/groups/OE0219/internal/KatherineK/ATACseq/differential-eICH/HSC-vs-Mono/differential_epiCHAOS_TFBSs.Rds")
het <- readRDS("/omics/groups/OE0219/internal/KatherineK/ATACseq/differential-eICH/HSC-vs-Mono/epiCHAOS_scores_TFBSs.Rds")
het$state <- het$state %>% str_remove("group1-|group2-")

#--- check
length(results)
names(results)
results[[10]]


#--- compute p-values and differences per region
temp <- data.frame(row.names=names(results))
for (i in 1:length(results)) {

  state <- names(results)[i]
  
  #--- the proportion of the e.g. 1000 permutations (plus the true comparison) in which the difference in epiCHAOS score between two groups is greater than that of the test comparison
  pval <- length(which(c(results[[i]]$test, results[[i]]$mix$dif.dist)>=results[[i]]$test))/length(c(results[[i]]$test, results[[i]]$mix$dif.dist))
  
  #--- the difference in epiCHAOS score between the two test groups
  dif <- het$mean.het[het$state==state & het$group=="group1"] - het$mean.het[het$state==state & het$group=="group2"]

  temp[i, "pval"] <- pval
  temp[i, "dif"] <- dif

}

#--- adjust labels for plotting
temp$TF <- rownames(temp) %>% str_split(" - ") %>% lapply("[", 2) %>% str_split("_") %>% lapply("[", 1) %>% unlist()
temp$label <- temp$TF
temp$label[temp$pval>0.01|temp$dif<0.42] <- ""

#--- Figure 4F
ggplot(temp, aes(y = -log10(pval), x = dif, label=label)) +
  geom_point(size=1, alpha=0.6)+
  labs(x="EpiCHAOS Diff.", subtitle = "HSCs vs CD8T")+
  lims(y=c(0,3.1), x=c(-.75, 0.75))+
  ggrepel::geom_text_repel(max.overlaps = 1000, size=3)+
  theme_classic()
