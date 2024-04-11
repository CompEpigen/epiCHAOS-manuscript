

#--- (i) test epiCHAOS in synthetic datasets of controlled heterogeneity
#--- (ii) test epiCHAOS in mixtures of celltypes

# load packages
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

source("/omics/groups/OE0219/internal/KatherineK/ATACseq/eITH-test-scripts/jaccard.R")

setwd("/omics/groups/OE0219/internal/KatherineK/ATACseq/epiCHAOS-Figures/")


#--- (i) testing epiCHAOS in synthetic datasets of controlled heterogeneity
datasets <- readRDS("/omics/groups/OE0219/internal/KatherineK/ATACseq/variability_methods/perturbed_datasets_signac_monocytes.Rds")
names(datasets)
#datasets$perturb.50 <- datasets$random.add50 <- datasets$random.rem50 <- NULL

lapply(datasets, dim)

# het <- compute.eITH(datasets)
# het
# 
# het <- compute.eITH(datasets[c("baseline", "rem10.hom", "rem20.hom", "rem30.hom", "rem40.hom", "rem50.hom","add10.hom", "add20.hom", "add30.hom", "add40.hom", "add50.hom")])
# het$n <- het$state %>% str_remove("rem|add") %>% str_remove("\\.") %>% str_remove("hom") %>% as.numeric()
# het$n[het$state=="baseline"] <- 0
# p1 <- ggplot(het[het$state %in% c("baseline", "rem10.hom", "rem20.hom", "rem30.hom", "rem40.hom", "rem50.hom"),], aes(y=mean.het+0.01, x=reorder(state, n))) +
#   geom_bar(stat="identity", width = 0.3, fill="steelblue4", alpha=0.5)+
#   labs( x="", y="epiCHAOS") +
#   #lims(y=c(0,1.01))+
#   coord_cartesian(ylim=c(0,1.01))+
#   theme_classic() +
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
# 
# p2 <- ggplot(het[het$state %in% c("baseline", "add10.hom", "add20.hom", "add30.hom", "add40.hom", "add50.hom"),], aes(y=mean.het+0.01, x=reorder(state, n))) +
#   geom_bar(stat="identity", width = 0.3, fill="steelblue4", alpha=0.5)+
#   labs( x="", y="epiCHAOS") +
#   #lims(y=c(0,1.01))+
#   coord_cartesian(ylim=c(0,1.01))+
#   theme_classic() +
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


het <- compute.eITH(datasets[c("baseline", "rem10", "rem20", "rem30", "rem40", "rem50","add10", "add20", "add30", "add40", "add50")])
het$n <- het$state %>% str_remove("rem|add") %>% str_remove("\\.") %>% str_remove("hom") %>% as.numeric()
het$n[het$state=="baseline"] <- 0

p3 <- ggplot(het[het$state %in% c("baseline", "rem10", "rem20", "rem30", "rem40", "rem50"),], aes(y=mean.het+0.01, x=reorder(state, mean.het))) +
  geom_bar(stat="identity", width = 0.3, fill="steelblue4", alpha=0.5)+
  labs( x="", y="epiCHAOS") +
  #lims(y=c(0,1.01))+
  coord_cartesian(ylim=c(0,1.01))+
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

p4 <- ggplot(het[het$state %in% c("baseline", "add10", "add20", "add30", "add40", "add50"),], aes(y=mean.het+0.01, x=reorder(state, mean.het))) +
  geom_bar(stat="identity", width = 0.3, fill="steelblue4", alpha=0.5)+
  labs( x="", y="epiCHAOS") +
  #lims(y=c(0,1.01))+
  coord_cartesian(ylim=c(0,1.01))+
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

pdf("Figure 1/barplots_synthetic_datasets.pdf", 7, 2.5)
ggarrange(p1,p2,p3,p4, ncol=4, align = "h")
dev.off()

# plot all in a single plot
het$group[grepl("rem",het$state)] <- "remove counts"
het$group[grepl("add",het$state)] <- "add counts"
het$n <- as.factor(het$n)

pdf("Figure 1/barplots_synthetic_datasets_oneplot.pdf", 5, 2.5)
ggplot(het, aes(y=mean.het+0.01, x=n, fill=group)) +
  geom_bar(stat="identity", alpha=0.5, position=position_dodge2(), width=0.5)+
  labs( x="", y="epiCHAOS") +
  scale_fill_manual(values = c("grey40", "red3"))+
  #lims(y=c(0,1.01))+
  coord_cartesian(ylim=c(0,1.01))+
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()

svg("Figure 1/barplots_synthetic_datasets_oneplot.svg", 5, 2.5)
ggplot(het, aes(y=mean.het+0.01, x=n, fill=group)) +
  geom_bar(stat="identity", alpha=0.5, position=position_dodge2(), width=0.5)+
  labs( x="", y="epiCHAOS") +
  scale_fill_manual(values = c("grey40", "red3"))+
  #lims(y=c(0,1.01))+
  coord_cartesian(ylim=c(0,1.01))+
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()



#--- (ii) perturbing a random dataset towards increased homogeneity while maintaing a constant number of counts

# create 100 synthetic datasets with equal total number of counts, which begin random, and gain increasing homogeneity, as we remove 1's from selected n rows, and add them to a different selected nrows

x <- matrix(nrow=10000,ncol=100)
x[is.na(x)] <- 0
x[sample(length(x), 200000)] <- 1
rownames(x) <- 1:nrow(x)
colnames(x) <- paste0("c", 1:ncol(x))
head(x)

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

lapply(datasets, dim)
lapply(datasets, sum)

het <- compute.eITH(datasets)
het$n <- het$state %>% str_remove("t") %>% as.numeric()
het$n[het$state=="baseline"] <- 0

pdf("Figure 1/scatterplot_epiCHAOS_vs_controlled_heterogeneity.pdf", 4, 3)
ggplot(het, aes(x = -n, y = mean.het)) +
  geom_point(size=1, alpha=0.6)+
  labs(y="epiCHAOS", x="Controlled heterogeneity")+
  stat_cor(method = "spearman")+
  lims(x=c(-10000, 200))+
  theme_classic() +
  theme(axis.text.x = element_blank())
dev.off()

svg("Figure 1/scatterplot_epiCHAOS_vs_controlled_heterogeneity.svg", 4, 3)
ggplot(het, aes(x = -n, y = mean.het)) +
  geom_point(size=1, alpha=0.6)+
  labs(y="epiCHAOS", x="Controlled heterogeneity")+
  stat_cor(method = "spearman")+
  lims(x=c(-10000, 200))+
  theme_classic() +
  theme(axis.text.x = element_blank())
dev.off()

#--- (iii) create a series of synthetic random datasets with incrementarlly increasing coverage to see how eITH score behaves in relation to sparsity

x <- matrix(nrow=10000,ncol=100)
x[is.na(x)] <- 0
x[sample(length(x), 10000)] <- 1
rownames(x) <- 1:nrow(x)
colnames(x) <- paste0("c", 1:ncol(x))

datasets <- list()
lapply(datasets, dim)

for (i in 1:100) {
  temp <- x
  temp[sample(length(temp), i*1000)] <- 1
  datasets[[paste0("t", i)]] <- temp 
}

het <- compute.eITH(datasets)
het
het$ncount <- unlist(lapply(datasets, sum))

temp <- het[,c("mean.het", "state")] %>% unique()
temp$n <- 1:nrow(temp)

pdf("Figure 1/scatterplot_epiCHAOS_vs_ncount_random_datasets.pdf", 4, 3)
ggplot(temp, aes(x = n, y = mean.het)) +
  geom_point(size=1, alpha=0.6)+
  labs(y="epiCHAOS", x="Total count")+
  ggpubr::stat_cor()+
  #lims(y=c(0,1)) +
  theme_classic() +
  theme(axis.text.x = element_blank())
dev.off()

svg("Figure 1/scatterplot_epiCHAOS_vs_ncount_random_datasets.svg", 4, 3)
ggplot(temp, aes(x = n, y = mean.het)) +
  geom_point(size=1, alpha=0.6)+
  labs(y="epiCHAOS", x="Total count")+
  ggpubr::stat_cor()+
  #lims(y=c(0,1)) +
  theme_classic() +
  theme(axis.text.x = element_blank())
dev.off()


