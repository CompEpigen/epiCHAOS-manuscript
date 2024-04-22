

#--- compare scATAC-based and DNA methylation heterogeneity at different genomic regions

meth.var <- read.csv("/omics/groups/OE0219/internal/KatherineK/Genomic_Annotations/variably_methylated_CpGs_hematopoesis_hg38.csv")

#--- liftover using a GRanges object
library(liftOver)
grObject <- meth.var[,1:3] %>% GRanges()
grObject$variance <- meth.var$variance
chainObject <- import.chain("/omics/groups/OE0219/internal/KatherineK/Genomic_Annotations/liftover_chain_files/hg38ToHg19.over.chain")
results <- as.data.frame(liftOver(grObject, chainObject)) %>% GRanges()
results

#--- calculate variance at each LOLA regions

# encode TFBS
lola <- get(load("/omics/groups/OE0219/internal/KatherineK/Genomic_Annotations/LOLA/hg19/encode_tfbs/encode_tfbs.RData"))
index <- read.table("/omics/groups/OE0219/internal/KatherineK/Genomic_Annotations/LOLA/hg19/encode_tfbs/index.txt", header = T)
names(lola) <- index$filename
#names(lola) <- paste0(index$cellType, "-", index$antibody, "-", index$treatment)
#lola <- lola[grepl("K562", names(lola))]
length(lola)

var.per.region <- data.frame(row.names = names(lola))
for (i in names(lola)) {
  print(i)
  regions <- lola[[i]]
  ovl <- subsetByOverlaps(results, regions)
  var.per.region[i, "mean.var"] <- mean(ovl$variance)
}

var.per.region %>% arrange(mean.var) %>% tail(20)
var.per.region$region <- rownames(var.per.region)
var.per.region$label <- ""
var.per.region$label[grepl("Ezh2|Suz12", rownames(var.per.region))] <- "PRC2"
var.per.region$label[grepl("Ctcf", rownames(var.per.region))] <- "CTCF"
var.per.region$label[grepl("Rad21|Smc3", rownames(var.per.region))] <- "Cohesin"


pdf("dotplot_DNAm_var_hematopoietic_highlight_PRC2_cohesin_CTCF.pdf", 6, 3)
ggplot(var.per.region[grepl("K562", var.per.region$region),], aes(y=mean.var, x=reorder(region, mean.var), label=label)) +
  geom_point(size=1)+
  labs( x="", y="DNAm var")+
  ggrepel::geom_text_repel(max.overlaps = 100)+
  theme_classic() +
  expand_limits(x= c(-2, 5))+
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
dev.off()

# ggplot(var.per.region, aes(y=mean.var, x=reorder(region, mean.var), color=label)) +
#   geom_point(size=1, alpha=0.7)+
#   labs( x="", y="DNAm var")+
#   theme_classic() +
#   scale_color_manual(values=c("grey", "orange", "steelblue4", "red3")) +
#   expand_limits(x= c(-20, 5))+
#   theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
# dev.off()

#--- compute same with epiCHAOS in HSCs
het <- readRDS("/omics/groups/OE0219/internal/KatherineK/ATACseq/per-region-eITH/epiCHAOS_encodeTFBS_HSCs.Rds")
#temp <- het[het$celltype=="01_HSC",]
temp <- het
temp$state <- temp$state %>% str_remove(paste0("-01_HSC"))
temp <- merge(temp, var.per.region, by.x="state", by.y=0)
temp$label <- ifelse(grepl("Ezh2|Suz12",temp$state), T, F)
p1 <- ggplot(temp, aes(y=mean.het, x=mean.var, color=label)) +
  geom_point( alpha=0.8, size=1)+
  scale_color_manual(values = c("grey20", "salmon"))+
  labs( x="DNAm variation", y="epiCHAOS", subtitle="PRC2") +
  theme_classic() +
  NoLegend()

temp$label <- ifelse(grepl("Ctcf",temp$state), T, F)
p2 <- ggplot(temp, aes(y=mean.het, x=mean.var, color=label)) +
  geom_point( alpha=0.8, size=1)+
  scale_color_manual(values = c("grey20", "salmon"))+
  labs( x="DNAm variation", y="epiCHAOS", subtitle="CTCF") +
  theme_classic() +
  NoLegend()

temp$label <- ifelse(grepl("Rad21|Smc3",temp$state), T, F)
p3 <- ggplot(temp, aes(y=mean.het, x=mean.var, color=label)) +
  geom_point( alpha=0.8, size=1)+
  scale_color_manual(values = c("grey20", "salmon"))+
  labs( x="DNAm variation", y="epiCHAOS", subtitle="Cohesin") +
  theme_classic() +
  NoLegend()

pdf("/omics/groups/OE0219/internal/KatherineK/ATACseq/epiCHAOS-vs-DNAmvar/scatterplots_epiCHAOS_vs_DNAm_var_TFBS_HSCs.pdf", 8, 2)
ggarrange(p1,p2,p3, ncol=3)
dev.off()


temp[grepl("Ctcf", temp$state),"state"] <- "CTCF"
temp[grepl("Rad21|Smc3", temp$state),"state"] <- "Cohesin"
temp[grepl("Ezh2|Suz12", temp$state),"state"] <- "PRC2"

# for (i in c("CTCF", "PRC2", "Cohesin")) {
#   temp$mean.het[temp$state==i] <- mean(temp$mean.het[temp$state==i])
#   temp$mean.var[temp$state==i] <- mean(temp$mean.var[temp$state==i])
# }

temp <- temp[,c("state", "mean.het", "mean.var")] %>% unique()
temp$label <- ifelse(temp$state %in% c("CTCF", "PRC2", "Cohesin"), temp$state, "")

saveRDS(temp, file = "/omics/groups/OE0219/internal/KatherineK/ATACseq/epiCHAOS-vs-DNAmvar/epiCHAOS_vs_mvar_scores.Rds")

pdf("/omics/groups/OE0219/internal/KatherineK/ATACseq/epiCHAOS-vs-DNAmvar/scatterplot_epiCHAOS_vs_DNAm_var_TFBS_HSCs_average_groups.pdf", 4, 3)
ggplot(temp, aes(y=mean.het, x=mean.var, label=label)) +
  geom_point( alpha=0.5, size=1)+
  ggrepel::geom_text_repel()+
  labs( x="DNAm variation", y="epiCHAOS") +
  geom_hline(yintercept=mean(temp$mean.het), lty="dashed", color="red3")+
  geom_vline(xintercept=mean(temp$mean.var), lty="dashed", color="red3")+
  theme_bw() +
  NoLegend()
dev.off()