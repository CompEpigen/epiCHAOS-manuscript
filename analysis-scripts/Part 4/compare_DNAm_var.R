

#--- compare scATAC-based and DNA methylation heterogeneity at different genomic regions

#--- DNA methylation variation per-CpG site from WGBS of HSCs (hg38)
meth.var <- read.csv("/omics/groups/OE0219/internal/KatherineK/Genomic_Annotations/variably_methylated_CpGs_hematopoesis_hg38.csv")

#--- liftover using a GRanges object (hg19)
library(liftOver)
grObject <- meth.var[,1:3] %>% GRanges()
grObject$variance <- meth.var$variance
chainObject <- import.chain("/omics/groups/OE0219/internal/KatherineK/Genomic_Annotations/liftover_chain_files/hg38ToHg19.over.chain")
results <- as.data.frame(liftOver(grObject, chainObject)) %>% GRanges()
results

#--- calculate variance at each ENCODE TFBS

#--- load encode TFBS
lola <- get(load("/omics/groups/OE0219/internal/KatherineK/Genomic_Annotations/LOLA/hg19/encode_tfbs/encode_tfbs.RData"))
index <- read.table("/omics/groups/OE0219/internal/KatherineK/Genomic_Annotations/LOLA/hg19/encode_tfbs/index.txt", header = T)
names(lola) <- index$filename
length(lola)

#--- compute average of the per-CpG variances for each TFBS
var.per.region <- data.frame(row.names = names(lola))
for (i in names(lola)) {
  
  print(i)
  regions <- lola[[i]]
  
  #--- CpGs overlapping with the selected TFBS
  ovl <- subsetByOverlaps(results, regions)
  
  #--- mean of variances of overlapping CpG
  var.per.region[i, "mean.var"] <- mean(ovl$variance)
}

#--- check top
var.per.region$region <- rownames(var.per.region)
var.per.region %>% arrange(mean.var) %>% tail(20)

#--- set labels for plotting, highlighting PRC2, CTCF and Cohesin binding sites which had higest epiCHAOS scores
var.per.region$label <- ""
var.per.region$label[grepl("Ezh2|Suz12", rownames(var.per.region))] <- "PRC2"
var.per.region$label[grepl("Ctcf", rownames(var.per.region))] <- "CTCF"
var.per.region$label[grepl("Rad21|Smc3", rownames(var.per.region))] <- "Cohesin"


#--- load epiCHAOS scores previously computed at each ENCODE TFBS in HSCs
het <- readRDS("/omics/groups/OE0219/internal/KatherineK/ATACseq/per-region-eITH/epiCHAOS_encodeTFBS_HSCs.Rds")

temp <- het
temp$state <- temp$state %>% str_remove(paste0("-01_HSC"))
temp <- merge(temp, var.per.region, by.x="state", by.y=0)
temp$label <- ifelse(grepl("Ezh2|Suz12",temp$state), T, F)

#--- highlights for plotting
temp[grepl("Ctcf", temp$state),"state"] <- "CTCF"
temp[grepl("Rad21|Smc3", temp$state),"state"] <- "Cohesin"
temp[grepl("Ezh2|Suz12", temp$state),"state"] <- "PRC2"

temp <- temp[,c("state", "mean.het", "mean.var")] %>% unique()
temp$label <- ifelse(temp$state %in% c("CTCF", "PRC2", "Cohesin"), temp$state, "")

#saveRDS(temp, file = "/omics/groups/OE0219/internal/KatherineK/ATACseq/epiCHAOS-vs-DNAmvar/epiCHAOS_vs_mvar_scores.Rds")

#--- set to range of 0-1 for easier interpretation
temp$mean.var <- (temp$mean.var - min(temp$mean.var)) / (max(temp$mean.var) - min(temp$mean.var)) 

#--- Figure 4B. scatter plot comparing DNA methylation variation & epiCHAOS scores per- ENCODE TFBS
ggplot(temp, aes(y=mean.het, x=mean.var, color=label)) +
  geom_point( alpha=0.5, size=1)+
  scale_color_manual(values = c("grey10", "royalblue", "orange2", "red3")) +
  labs( x="DNAm variation", y="epiCHAOS") +
  geom_hline(yintercept=mean(temp$mean.het), lty="dashed", color="red3")+
  geom_vline(xintercept=mean(temp$mean.var), lty="dashed", color="red3")+
  theme_classic() +
  NoLegend()