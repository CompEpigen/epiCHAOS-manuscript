

#--- testing performance of epiCHAOS using synthetic datasets of controlled heterogeneity (relating to Figure 1B-D)

#--- load packages
library(magrittr)
library(stringr)
library(ggplot2)
library(ggbeeswarm)
library(GenomicFeatures)
library(Signac)
library(Seurat)
library(dplyr)
library(RColorBrewer)
library(pheatmap)


#--- Figure 1A. perturbing a random dataset towards increased homogeneity while maintaing a constant number of counts

#--- create 100 synthetic datasets with equal total number of counts, which begin random, and gain increasing homogeneity, as we remove 1's from selected n rows, and add them to a different selected nrows

#--- first create a random binary matrix with 10000 rows/"loci" and 100 columns/"cells", randomly populate it with 200000 1's
x <- matrix(nrow=10000,ncol=100)
x[is.na(x)] <- 0
x[sample(length(x), 200000)] <- 1
rownames(x) <- 1:nrow(x)
colnames(x) <- paste0("c", 1:ncol(x))
head(x)

#--- add to the series of matrices, each time randomly selecting 50% of the rows, removing 100, 200, 300, 400 etc. 1's from those selected rows, and adding the same number of 1's to a different 50% of rows
#--- heterogenetiy should be highest in the first (randomly populated) matrix, and should gradually decrease in successive matrices
datasets <- list()
datasets$baseline <- x

for (i in seq(100,10000, 100)) {
  temp <- x
  select.rows <- sample(nrow(temp), nrow(temp)/2)
  to.remove <- sample(which(temp[select.rows, ]==1), i)
  temp[select.rows,][to.remove] <- 0
  to.add <- sample(which(temp[setdiff(rownames(temp), select.rows), ]==0), i)
  temp[setdiff(rownames(temp), select.rows),][to.add] <- 1
  
  datasets[[paste0("t",i)]] <- temp
}


#--- check that each matrix has the same total count
lapply(datasets, sum)

#--- compute epiCHAOS score
het <- compute_eITH(datasets)

#--- "n" represents the expected level of homogeneity, being zero in the baseline (random) dataset, so "-n" should correlate with the epiCHAOS score
het$n <- het$state %>% str_remove("t") %>% as.numeric()
het$n[het$state=="baseline"] <- 0

cor.test(het$het.adj, het$n, method = "spearman")

#--- plot the correlation between epiCHAOS scores and controlled heterogeneity ("-n") (Figure 1B)
ggplot(het, aes(x = -n, y = het.adj)) +
  geom_point(size=1, alpha=0.6)+
  labs(y="epiCHAOS", x="Controlled heterogeneity")+
  stat_cor(method = "spearman")+
  lims(x=c(-10000, 200))+
  theme_classic() +
  theme(axis.text.x = element_blank())


#--- Figure 1C. testing the increased heterogeneity is detected both in cases of increasing and decreasing counts

#--- load datasets created by perturbing scATAC-seq data from human monocyte by increasing/decreasing 50% of 1's
datasets <- readRDS("synthetic-data-tests/perturbed_data_monocytes.Rds")
names(datasets)

lapply(datasets, dim)

#--- compute epiCHAOS scores on datasets where 10-50% of counts are added and where 10-50% of counts are removed
het <- compute.eITH(datasets[c("baseline", "rem10", "rem20", "rem30", "rem40", "rem50","add10", "add20", "add30", "add40", "add50")])
het$n <- het$state %>% str_remove("rem|add") %>% str_remove("\\.") %>% str_remove("hom") %>% as.numeric()
het$n[het$state=="baseline"] <- 0

p1 <- ggplot(het[het$state %in% c("baseline", "rem10", "rem20", "rem30", "rem40", "rem50"),], aes(y=mean.het+0.01, x=reorder(state, mean.het))) +
  geom_bar(stat="identity", width = 0.3, fill="steelblue4", alpha=0.5)+
  labs( x="", y="epiCHAOS") +
  #lims(y=c(0,1.01))+
  coord_cartesian(ylim=c(0,1.01))+
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

p2 <- ggplot(het[het$state %in% c("baseline", "add10", "add20", "add30", "add40", "add50"),], aes(y=mean.het+0.01, x=reorder(state, mean.het))) +
  geom_bar(stat="identity", width = 0.3, fill="steelblue4", alpha=0.5)+
  labs( x="", y="epiCHAOS") +
  #lims(y=c(0,1.01))+
  coord_cartesian(ylim=c(0,1.01))+
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

#--- plot to show that epiCHAOS scores increase with increasing perturbation regardless of whether counts are added or removed
ggarrange(p1,p2, ncol=2, align = "h")




#--- Figure 1D. create a series of synthetic random datasets with incrementarlly increasing coverage to see how epiCHAOS score behaves in relation to sparsity

#--- each matrix in the series will have 10000 rows/"loci" and 100 columns/"cells"
x <- matrix(nrow=10000,ncol=100)
x[is.na(x)] <- 0

#--- add 10000 random counts to the first matrix in the series. This will be the lowest count.
x[sample(length(x), 10000)] <- 1
rownames(x) <- 1:nrow(x)
colnames(x) <- paste0("c", 1:ncol(x))
x[1:4,1:10]


datasets <- list()
lapply(datasets, dim)

#--- let the first matrix in the series be "x" created above which contains 10000 1's
#--- each subsequent matrix in the series will have 1000 more 1's than the previous
for (i in 1:100) {
  
  temp <- x
  temp[temp==0][sample(length(temp[temp==0]), i*1000)] <- 1
  datasets[[paste0("t", i)]] <- temp 
}

#--- check that the counts in each dataset is incrementally increasing by 1000
lapply(datasets, sum)

#--- compute epiCHAOS scores
het <- compute.eITH(datasets)

#--- add column to "het" to be the total number of 1's in each matrix
#het$ncount <- unlist(lapply(datasets, sum))

#--- otherwise, just note the rank of their counts from 1:nrow
temp <- het[,c("mean.het", "state")] %>% unique()
temp$n <- 1:nrow(temp)

#--- plot to confirm no correlation between epiCHAOS score and tota counts (Figure 1D)
ggplot(temp, aes(x = n, y = mean.het)) +
  geom_point(size=1, alpha=0.6)+
  labs(y="epiCHAOS", x="Total count")+
  ggpubr::stat_cor()+
  #lims(y=c(0,1)) +
  theme_classic() +
  theme(axis.text.x = element_blank())


