
#--- download scNMT-seq DNA methylation data from Argelaguet et al. 2019
#--- processed data are available within the SingleCellMultiModal R package

set.seed(10)

#BiocManager::install("SingleCellMultiModal")
library(SingleCellMultiModal)
library(MultiAssayExperiment)

scNMT("mouse_gastrulation", version = "2.0.0")

#--- select data for genomic regions of interest
nmt <- scNMT("mouse_gastrulation", mode = c("met_cgi", "met_genebody", "met_promoter",  "acc_promoter", "acc_genebody", "acc_cgi"), version = "1.0.0", dry.run = FALSE)

#--- check data
# nmt
# nmt@ExperimentList@listData$met_promoter[1:10,1:10]
# nmt@ExperimentList@listData$met_promoter %>% hist()
# nmt@colData$lineage10x_2 %>% table()
# nmt@colData$stage %>% table()
# nmt@colData$stage_lineage %>% table()

#--- compute epiCHAOS scores on single cell methylome data for promoters, gene bodies, CpG islands
gg <- list()
het.list <- list()
for (i in c("met_promoter", "met_genebody", "met_cgi")) {
  temp <- nmt@ExperimentList@listData[[i]]
  temp[is.na(temp)] <- 0
  temp[temp<0.2] <- 0
  temp[temp>0] <- 1
  
  #--- create matrices for each cell type
  datasets <- list()
  for (celltype in na.omit(unique(nmt@colData$lineage10x_2))) {
    ids <- nmt@colData$cellID[nmt@colData$lineage10x_2==celltype]
    if (length(ids)<40) { next }
    ncells <- min(c(50, length(ids)))
    
    
    datasets[[celltype]] <- temp[,colnames(temp) %in% ids]
  }
  
  #--- compute epiCHAOS scores
  het <- compute_eITH(datasets)
  het$group <- ifelse(het$state=="Epiblast", "epiblast", "other")
  het$state <- het$state %>% str_replace("_", " ")
  
  #--- save barplots
  het.list[[i]] <- het
  gg[[i]] <- ggplot(het, aes(x = reorder(state, het.adj), y = het.adj+0.01, fill=group)) +
    geom_bar(stat="identity", position = "dodge", alpha=0.8, width = 0.7)+
    scale_fill_manual(values = rev(c("black", "steelblue4")))+
    labs(x="", y="epiCHAOS")+
    theme_classic() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  
}

library(ggpubr)


#--- write scres to files
# het.list$met_promoter %>% arrange(desc(mean.het)) %>% write.csv("/omics/groups/OE0219/internal/KatherineK/ATACseq/epiCHAOS-supplementary-data/epiCHAOS_scNMTseq_meth_promoters.csv")
# het.list$met_genebody %>% arrange(desc(mean.het)) %>% write.csv("/omics/groups/OE0219/internal/KatherineK/ATACseq/epiCHAOS-supplementary-data/epiCHAOS_scNMTseq_meth_genebody.csv")
# het.list$met_cgi %>% arrange(desc(mean.het)) %>% write.csv("/omics/groups/OE0219/internal/KatherineK/ATACseq/epiCHAOS-supplementary-data/epiCHAOS_scNMTseq_meth_cgi.csv")


#--- save plot results
pdf("/omics/groups/OE0219/internal/KatherineK/ATACseq/epiCHAOS-Figures/Figure 5/barplot_gastrulation_scNMTseq.pdf", 6, 3)
plot <- ggpubr::ggarrange(plotlist = gg, ncol=3, legend=F)
annotate_figure(plot, top=text_grob("promoters, gene bodies, CpG Islands"))
dev.off()

svg("/omics/groups/OE0219/internal/KatherineK/ATACseq/epiCHAOS-Figures/Figure 5/barplot_gastrulation_scNMTseq.svg", 6, 3)
annotate_figure(plot, top=text_grob("promoters, gene bodies, CpG Islands"))
dev.off()


#--- compare epiCHAOS on methylation vs ATAC from same cells

#--- previously saved methylation scores for promoters
het <- read.csv("/omics/groups/OE0219/internal/KatherineK/ATACseq/epiCHAOS-supplementary-data/epiCHAOS_scNMTseq_meth_promoters.csv")

#--- corresponding atac data for promoters
temp <- nmt@ExperimentList@listData$acc_promoter
temp[is.na(temp)] <- 0
temp[temp<0.2] <- 0
temp[temp>0] <- 1

#--- create matrices for each celltype
datasets <- list()
for (celltype in na.omit(unique(nmt@colData$lineage10x_2))) {
  ids <- nmt@colData$cellID[nmt@colData$lineage10x_2==celltype]
  if (length(ids)<40) { next }
  datasets[[celltype]] <- temp[,colnames(temp) %in% ids]
}

#--- compute and compare epiCHAOS scores from both modalities
het2 <- compute.eITH(datasets)
het2$state <- het2$state %>% str_replace("_", " ")

het <- merge(het, het2, by="state")
colnames(het)[c(3, 5)] <- c("DNAm", "ATAC")

#--- fit linear model to get r-squared
fit <- lm(het$DNAm~het$ATAC)

#--- scatterplot
svg("/omics/groups/OE0219/internal/KatherineK/ATACseq/epiCHAOS-Figures/Figure 5/scatterplot_DNAm_vs_ATAC_epiCHAOS_gastrulation.svg", 4, 3)
ggplot(het, aes(x=DNAm, y=ATAC, label=state)) +
  geom_point(size=3, alpha=0.5)+
  geom_smooth(color="lightsteelblue3", method="lm", se=F)+
  ggrepel::geom_text_repel(size=3)+
  lims(y=c(0, 1))+
  stat_cor()+
  labs( x="epiCHAOS (DNAm)", y="epiCHAOS (ATAC)", subtitle = paste0("R^2 = ", signif(summary(fit)$r.squared, 3))) +
  theme_classic()
dev.off()
