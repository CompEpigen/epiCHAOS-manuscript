
#--- compare transcriptional noise (measured using the coefficient of variation) between PRC2 targets and non-PRC2 target genes (related to Figure S6)

library(msigdbr)

#--- scRNAseq data for hematopoietic cells from Granja et al., subset for bone-marrow-derived (BM), subset fo HSCs
rna <- readRDS("/omics/groups/OE0219/internal/KatherineK/data/scATAC/Granja2019/scRNA-Healthy-Hematopoiesis-191120.rds")
rna <- rna[,grepl("BMMC", rna$Group)]
gex <- rna@assays$data$counts

#--- create a seurat object, then normalise
gex <- CreateSeuratObject(counts = rna@assays@.xData$data$counts, meta.data = data.frame(rna@colData))
gex <- NormalizeData(gex)


#--- load msigdb gene sets
gene_sets = msigdbr(species = "Homo sapiens")
gene_sets <- gene_sets[gene_sets$gene_symbol %in% rownames(gex),]

#--- subset data for HSCs
celltype <- "01_HSC"
ge <- gex@assays$RNA$data[,gex@meta.data$BioClassification==celltype & grepl("BMMC", gex@meta.data$Group)]
ge <- ge[,sample(ncol(ge), 200)]

#--- calculate coefficient of variation for every gene
cov <- sapply(as.data.frame(t(ge)), function(x) sd(x) / mean(x) * 100)
cov

#--- note also that CV correlates with the gene expression levels, which has been previously described
sum <- apply(ge, 1, sum)
sum <- sum[sum>0]

#--- list of PRC2 target genes
prc2 <- gene_sets$gene_symbol[gene_sets$gs_name == "BENPORATH_PRC2_TARGETS"]
cov.prc2 <- cov[names(cov) %in% prc2] %>% na.omit() %>% as.numeric()
sum.prc2 <- sum[names(sum) %in% prc2] %>% na.omit() %>% as.numeric()

#--- check relationship between CV and gene expression level
plot(cov.prc2, sum.prc2)


#--- housekeeping genes
hk <- c("Actb", "Atp5f1", "Atp5pb", "B2m", "Gapdh", "Hprt1", "Hprt", "Pgk1", "Rer1", "Rpl13a", "Rpl27", "Sdha", "Tbp", "Ubc") %>% str_to_upper()
cov.hk <- cov[names(cov) %in% hk]

#--- bivalent genes
bivalent <- read.table("/omics/groups/OE0219/internal/KatherineK/Genomic_Annotations/gene_annotations/bivalent_genes_hg19_Court_2017.txt", header = T, sep = "\t")
bivalent <- bivalent$Gene
cov.biv <- cov[names(cov) %in% bivalent] %>% na.omit() %>% as.numeric()

#--- randomly select 1000 genes which are not PRC2 targets, bivalent or hk genes
random <- names(cov) %>% setdiff(c(bivalent, hk, prc2)) %>% sample(1000)
cov.random <- cov[names(cov) %in% random]  %>% na.omit() %>% as.numeric()
sum.random <- sum[names(sum) %in% random]  %>% na.omit() %>% as.numeric()

#--- check relationship between CV and gene expression level
plot(cov.random, sum.random)

#--- compare coefficient of variaton for the different gene sets
temp <- data.frame(cv=c(cov.hk, cov.random, cov.biv, cov.prc2), type=c(rep("HK", length(cov.hk)), rep("Other", length(cov.random)), rep("Bivalent", length(cov.biv)), rep("PRC2", length(cov.prc2))))

#saveRDS(temp, file = "CV_pergene_HSCs_categorised.Rds")

#--- Figure S6 boxplot
ggplot(temp[temp$type %in% c("Other", "PRC2"),], aes(y=log(cv), x=reorder(type, cv))) +
  geom_boxplot(linewidth=0.5, fill="steelblue4", alpha=0.5)+
  stat_compare_means(comparisons = list(c(1,2))) +
  labs( x="", y="transcriptional noise") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, color="black"))
