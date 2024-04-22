#--- investigate CV, EVA and scde methods for computing expression variability in ordar to correlate per-geneset expression variability with epiCHAOS

library(GSReg)
library(msigdbr)

#--- scRNAseq data for HSCs
rna <- readRDS("/omics/groups/OE0219/internal/KatherineK/data/scATAC/Granja2019/scRNA-Healthy-Hematopoiesis-191120.rds")
rna <- rna[,grepl("BMMC", rna$Group)]

hscs <- rna@colData[rna@colData$BioClassification=="01_HSC",] %>% rownames() %>% sample(100)
monos <- rna@colData[rna@colData$BioClassification=="11_CD14.Mono.1",] %>% rownames() %>% sample(100)
gex <- rna@assays$data$counts

# create a seurat object, then normalise and scale the data
gex <- CreateSeuratObject(counts = rna@assays@.xData$data$counts, meta.data = data.frame(rna@colData))

gex <- NormalizeData(gex)
counts <- gex@assays$RNA$data[,c(hscs, monos)] %>% as.matrix()
dim(counts)


# # create a list of gene sets
gene_sets = msigdbr(species = "Homo sapiens")
gene_sets <- gene_sets[gene_sets$gene_symbol %in% rownames(counts),]


#--- try coefficient of variation
celltype <- "01_HSC"
ge <- gex@assays$RNA$data[,gex@meta.data$BioClassification==celltype & grepl("BMMC", gex@meta.data$Group)]
dim(ge)
ge <- ge[,sample(ncol(ge), 200)]

# calculate coefficient of variation for every gene
cov <- sapply(as.data.frame(t(ge)), function(x) sd(x) / mean(x) * 100)
cov
sum <- apply(ge, 1, sum)
sum <- sum[sum>0]

# compare geneSummary.numeric_version()# compare gene sets
prc2 <- gene_sets$gene_symbol[gene_sets$gs_name == "BENPORATH_PRC2_TARGETS"]
cov.prc2 <- cov[names(cov) %in% prc2] %>% na.omit() %>% as.numeric()
sum.prc2 <- sum[names(sum) %in% prc2] %>% na.omit() %>% as.numeric()
plot(cov.prc2, sum.prc2)


# hk genes
hk <- c("Actb", "Atp5f1", "Atp5pb", "B2m", "Gapdh", "Hprt1", "Hprt", "Pgk1", "Rer1", "Rpl13a", "Rpl27", "Sdha", "Tbp", "Ubc") %>% str_to_upper()
cov.hk <- cov[names(cov) %in% hk]

# bivalent genes
bivalent <- read.table("/omics/groups/OE0219/internal/KatherineK/Genomic_Annotations/gene_annotations/bivalent_genes_hg19_Court_2017.txt", header = T, sep = "\t")
bivalent <- bivalent$Gene
cov.biv <- cov[names(cov) %in% bivalent] %>% na.omit() %>% as.numeric()

random <- names(cov) %>% setdiff(c(bivalent, hk, emt, prc2)) %>% sample(1000)
cov.random <- cov[names(cov) %in% random]  %>% na.omit() %>% as.numeric()
sum.random <- sum[names(sum) %in% random]  %>% na.omit() %>% as.numeric()
plot(cov.random, sum.random)


boxplot(cov.hk, cov.random, cov.biv, cov.prc2)
boxplot(cov.random, cov.prc2)
boxplot(sum[random], sum[prc2])
wilcox.test(cov.random, cov.biv)

temp <- data.frame(cv=c(cov.hk, cov.random, cov.biv, cov.prc2), type=c(rep("HK", length(cov.hk)), rep("Other", length(cov.random)), rep("Bivalent", length(cov.biv)), rep("PRC2", length(cov.prc2))))

saveRDS(temp, file = "/omics/groups/OE0219/internal/KatherineK/ATACseq/eICH-vs-tICH/CV_pergene_HSCs_categorised.Rds")

ggplot(temp, aes(y=cv, x=reorder(type, cv))) +
  geom_boxplot(linewidth=0.5, fill="honeydew3")+
  stat_compare_means(comparisons = list(c(1,2), c(2,3))) +
  #geom_jitter(color="grey30", size=0.1, alpha=0.3, width = 0.3)+
  labs( x="", y="transcriptional noise") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

pdf("boxplot_tx_noise_cv_prc2_targets_vs_other_genes_HSCs.pdf", 2.5, 4)
ggplot(temp[temp$type %in% c("Other", "PRC2"),], aes(y=log(cv), x=reorder(type, cv))) +
  geom_boxplot(linewidth=0.5, fill="steelblue4", alpha=0.5)+
  stat_compare_means(comparisons = list(c(1,2))) +
  #geom_jitter(color="grey30", size=0.1, alpha=0.3, width = 0.3)+
  labs( x="", y="transcriptional noise") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, color="black"))
dev.off()