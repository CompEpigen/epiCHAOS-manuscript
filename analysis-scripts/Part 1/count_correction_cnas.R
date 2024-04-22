
#---
#--- check effect of CNVs on heterogeneity score using CNAs caled by epiaceufinder in Hep-1 cell line
#---

set.seed(10)

library(ggplot2)
library(ggpubr)


#--- deletion chromosome 13

#--- CNV data for the Hep1 cell line inferred using epianeufinder
cnvs <- read.table("/omics/groups/OE0219/internal/KatherineK/ATACseq/epiAneufinder/epiAneufinder_results_Hep1/epiAneufinder_results/results_table.tsv", header = T)

# subgroup cells based on apparent subclonal deletion on chromosome 13
colnames(cnvs)[1] <- "chr"
cnvs <- cnvs[cnvs$chr=="chr13",]
coords <- cnvs[2]

gr <- GRanges(cnvs[,1:3])
cnvs <- cnvs[,4:ncol(cnvs)]
dim(cnvs)

# identify the frequently amplified region on chromosome 5 and subgroup cells based on that regions copy number
temp <- data.frame(coords=coords$start, copy=rowMeans(cnvs))

copyplot1 <- ggplot(data=temp, aes(x=coords, y=copy)) +
  geom_line(color="grey30") +
  geom_vline(xintercept = c(53000000, 107000000), lty="dashed", color="red3") +
  labs(x="position on chr13", y="mean copy number") +
  theme_classic()

# select approximate coordinates of subclonal deletion
select.rows <- cnvs[coords$start>53000000 & coords$start< 107000000,] %>% rownames()

# density plot
temp <- data.frame(copy=colMeans(cnvs[select.rows,]))

density1 <- ggplot(temp, aes(x=copy)) +
  geom_density(fill="grey90", color="grey30", alpha=0.8) +
  lims(x=c(-0.5, 2)) +
  labs(x="mean copy number") +
  theme_classic()


# note which cells have deletions and which are diploid at the selected region
dip <- names(cnvs)[colMeans(cnvs[select.rows,]) > 0.5 & colMeans(cnvs[select.rows,]) < 1.5] %>% str_remove("cell.SK.")
del <- names(cnvs)[colMeans(cnvs[select.rows,]) < 0.5] %>% str_remove("cell.SK.")

# load atac data for HCC cell lines
atac <- readRDS("/omics/groups/OE0219/internal/KatherineK/ATACseq/HCC-celllines/matrices/peak_matrix.Rds")
atac.gr <- atac@rowRanges
names(atac.gr) <- 1:length(atac.gr)

# subset for 
atac <- atac[,colnames(atac) %in% paste0("Hep-1#SK-", c(del, dip))]
colnames(atac)
atac <- atac@assays@data$PeakMatrix
dim(atac)

# subset regions corresponding to the altered region
select.sites <- subsetByOverlaps(atac.gr, gr[select.rows]) %>% names()

atac <- atac[as.numeric(select.sites),]

datasets <- list()
datasets$del <- atac[,colnames(atac) %in%  paste0("Hep-1#SK-", del)]                 
datasets$dip <- atac[,colnames(atac) %in%  paste0("Hep-1#SK-", dip)]                 


# random samples of 100 "diploid cells
datasets$del <- datasets$del[,1:100]
datasets$dip1 <- datasets$dip[,sample(ncol(datasets$dip), 100)]
datasets$dip2 <- datasets$dip[,sample(ncol(datasets$dip), 100)]
datasets$dip3 <- datasets$dip[,sample(ncol(datasets$dip), 100)]
datasets$dip4 <- datasets$dip[,sample(ncol(datasets$dip), 100)]
datasets$dip5 <- datasets$dip[,sample(ncol(datasets$dip), 100)]
datasets$dip <- NULL

lapply(datasets, dim)

# compute eICH
het <- compute.eITH.raw(datasets)
het.adj <- compute.eITH(datasets)

het$group <- ifelse(grepl("del", het$state), "deletion", "diploid")
het.adj$group <- ifelse(grepl("del", het.adj$state), "deletion", "diploid")

#--- plot correlations of epiCHAOS vs total counts before and after count adjustment
temp <- lapply(datasets, colSums) %>% lapply(mean) %>% unlist() %>% data.frame() %>% merge(het, by.x=0, by.y="state")
scatter1 <- ggplot(temp, aes(x=., y=mean.het, color=group)) +
  geom_point(size=3, alpha=0.7, aes(color=group))+
  scale_color_manual(values = c("steelblue4", "black"))+
  lims(y=c(0, 1)) +
  geom_smooth(color="grey30", method = "lm", se=F)+
  labs( x="Count", y="epiCHAOS") +
  theme_classic()

temp <- lapply(datasets, colSums) %>% lapply(mean) %>% unlist() %>% data.frame() %>% merge(het.adj, by.x=0, by.y="state")
scatter2 <- ggplot(temp, aes(x=., y=mean.het, color=group)) +
  geom_point(size=3, alpha=0.7, aes(color=group))+
  scale_color_manual(values = c("steelblue4", "black"))+
  lims(y=c(0, 1)) +
  labs( x="Count", y="epiCHAOS") +
  theme_classic()

ggpubr::ggarrange(copyplot1, density1, scatter1, scatter2, ncol=4, common.legend = T, widths = c(1, 1, 0.8, 0.8))



#--- repeat for an example of copy number gain, e.g. on chromosome 5

#--- CNV data for the Hep1 cell line inferred using epianeufinder
cnvs <- read.table("/omics/groups/OE0219/internal/KatherineK/ATACseq/epiAneufinder/epiAneufinder_results_Hep1/epiAneufinder_results/results_table.tsv", header = T)

# subgroup cells based on chromoosme 5 CNVs, which appear subclonal
colnames(cnvs)[1] <- "chr"
cnvs <- cnvs[cnvs$chr=="chr5",]
coords <- cnvs[2]

gr <- GRanges(cnvs[,1:3])
cnvs <- cnvs[,4:ncol(cnvs)]
dim(cnvs)

# identify the frequently amplified region on chromosome 5 and subgroup cells based on that regions copy number
temp <- data.frame(coords=coords$start, copy=rowMeans(cnvs))

copyplot2 <- ggplot(data=temp, aes(x=coords, y=copy)) +
  geom_line(color="grey30") +
  geom_vline(xintercept = c(30000000, 45000000), lty="dashed", color="red3") +
  labs(x="position on chr5", y="mean copy number") +
  theme_classic()

select.rows <- coords$start>30000000 & coords$start< 45000000

cnvs[select.rows,] %>% colMeans() %>% density() %>% plot()

temp <- data.frame(copy=colMeans(cnvs[select.rows,]))

density2 <- ggplot(temp, aes(x=copy)) +
  geom_density(fill="grey90", color="grey30", alpha=0.8) +
  lims(x=c(0.5, 2.5)) +
  labs(x="mean copy number") +
  theme_classic()

# note which cells have amplificationa nd which are diploid at the selected region
amp <- names(cnvs)[colMeans(cnvs[select.rows,]) > 1.5] %>% str_remove("cell.SK.")
dip <- names(cnvs)[colMeans(cnvs[select.rows,]) > 0.5 & colMeans(cnvs[select.rows,]) < 1.5] %>% str_remove("cell.SK.")

# atac data for HCC cell lines
atac <- readRDS("/omics/groups/OE0219/internal/KatherineK/ATACseq/HCC-celllines/matrices/peak_matrix.Rds")
atac.gr <- atac@rowRanges
names(atac.gr) <- 1:length(atac.gr)

# subset for 
atac <- atac[,colnames(atac) %in% paste0("Hep-1#SK-", c(amp, dip))]
colnames(atac)
atac <- atac@assays@data$PeakMatrix
dim(atac)

# subset regions corresponding to the chromosome 5 commonly amplified region
select.sites <- subsetByOverlaps(atac.gr, gr[290:420]) %>% names()

atac <- atac[as.numeric(select.sites),]

datasets <- list()
datasets$amp <- atac[,colnames(atac) %in%  paste0("Hep-1#SK-", amp)]                 
datasets$dip <- atac[,colnames(atac) %in%  paste0("Hep-1#SK-", dip)]                 

lapply(datasets, dim)
datasets$amp <- datasets$amp[,1:150]
datasets$dip1 <- datasets$dip[,sample(ncol(datasets$dip), 150)]
datasets$dip2 <- datasets$dip[,sample(ncol(datasets$dip), 150)]
datasets$dip3 <- datasets$dip[,sample(ncol(datasets$dip), 150)]
datasets$dip4 <- datasets$dip[,sample(ncol(datasets$dip), 150)]
datasets$dip5 <- datasets$dip[,sample(ncol(datasets$dip), 150)]
datasets$dip <- NULL

lapply(datasets, dim)
lapply(datasets, sum)

het <- compute.eITH.raw(datasets)
het.adj <- compute.eITH(datasets)

het$group <- ifelse(grepl("amp", het$state), "gain", "diploid")
het.adj$group <- ifelse(grepl("amp", het.adj$state), "gain", "diploid")

#--- plot correlations of epiCHAOS vs total counts before and after count adjustment
temp <- lapply(datasets, colSums) %>% lapply(mean) %>% unlist() %>% data.frame() %>% merge(het, by.x=0, by.y="state")
scatter3 <- ggplot(temp, aes(x=., y=mean.het, color=group)) +
  geom_point(size=3, alpha=0.7, aes(color=group))+
  scale_color_manual(values = c("black", "red3"))+
  lims(y=c(0, 1)) +
  geom_smooth(color="grey30", method = "lm", se=F)+
  labs( x="Count", y="epiCHAOS") +
  theme_classic()

temp <- lapply(datasets, colSums) %>% lapply(mean) %>% unlist() %>% data.frame() %>% merge(het.adj, by.x=0, by.y="state")
scatter4 <- ggplot(temp, aes(x=., y=mean.het, color=group)) +
  geom_point(size=3, alpha=0.7, aes(color=group))+
  scale_color_manual(values = c("black", "red3"))+
  lims(y=c(0, 1)) +
  labs( x="Count", y="epiCHAOS") +
  theme_classic()

# plot svgs to supplementary figure 1
svg("/omics/groups/OE0219/internal/KatherineK/ATACseq/epiCHAOS-Figures/Supplementary/S1/correct_count_test_cnas_deletion.svg", 12, 3)
ggarrange(copyplot1, density1, scatter1, scatter2, ncol=4, common.legend = T, widths = c(1.2, 1, 0.9, 0.9))
dev.off()

svg("/omics/groups/OE0219/internal/KatherineK/ATACseq/epiCHAOS-Figures/Supplementary/S1/correct_count_test_cnas_gain.svg", 12, 3)
ggarrange(copyplot2, density2, scatter3, scatter4, ncol=4, common.legend = T, widths = c(1.2, 1, 0.9, 0.9))
dev.off()

