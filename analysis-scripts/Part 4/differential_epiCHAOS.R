

#---
#--- test permutation approach to compute differential epigenetic heterogeneity at a given region
#--- 

setwd("/omics/groups/OE0219/internal/KatherineK/ATACseq/differential-eICH")

#--- take for example the set of bivalent genes
bivalent <- read.table("/omics/groups/OE0219/internal/KatherineK/Genomic_Annotations/gene_annotations/bivalent_genes_hg19_Court_2017.txt", header = T, sep = "\t")$Gene
bivalent <- bivalent %>% str_split(",") %>% unlist() %>% unique()

#--- select hg19 promoters
require("TxDb.Hsapiens.UCSC.hg19.knownGene")
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
promoters <- promoters(genes(txdb), upstream = 1500, downstream = 500)

# use gene sets table to translate bivalent gene ids to entrez ids
gene_sets = msigdbr(species = "Homo sapiens")
bivalent <- gene_sets$entrez_gene[gene_sets$gene_symbol %in% bivalent] %>% unique()

#--- scATAC-seq example dataset - data from normal hematopoiesis from Granja et al. 2019 - subset for HSCs and monocytes
atac <- readRDS("/omics/groups/OE0219/internal/KatherineK/data/scATAC/Granja2019/scATAC-Healthy-Hematopoiesis-191120.rds")
atac <- atac[,grepl("BMMC", atac$Group)]
hscs <- atac[,atac$BioClassification=="01_HSC"] %>% colnames() %>% sample(100)
hscs <- atac@assays$data$counts[,hscs]
hscs[hscs>1] <- 1
rownames(hscs) <- names(atac@rowRanges)

mono <- atac[,atac$BioClassification== "05_CMP.LMPP" ] %>% colnames() %>% sample(100)
mono <- atac@assays$data$counts[,mono]
mono[mono>1] <- 1
rownames(mono) <- names(atac@rowRanges)

dim(hscs)
dim(mono)

# peaks overlapping with bivalent genes
bivalent <- subsetByOverlaps(atac@rowRanges, promoters[names(promoters) %in% bivalent,]) %>% names()
hscs <- hscs[bivalent,]
mono <- mono[bivalent,]
merged <- cbind(hscs, mono)

select.sites <-  subsetByOverlaps(atac@rowRanges, promoters) %>% names() %>% sample(10000)
hscs <- hscs[select.sites,]
mono <- mono[select.sites,]
merged <- cbind(hscs, mono)


datasets <- list(hscs=hscs, mono=mono)
for (i in 1:100) {
  sample.cells <- c(colnames(merged)) %>% sample(100)
  datasets[[paste0("mix", i)]] <- merged[,sample.cells]
}

lapply(datasets, dim)

het <- compute.eITH(datasets)
het
#barplot(het$mean.het)

het.dist <- het[het$state %notin% c("mono", "hscs"),]

dif.dist <- c()
for (i in 1:nrow(het.dist)/2) {
  dif <- het.dist$mean.het[i]-het.dist$mean.het[i+1]
  dif.dist <- c(dif.dist, dif)
}

plot(density(dif.dist), xlim=c(-1,1), xlab="", main="")
hist(dif.dist, xlim=c(-1,1), xlab="", main="")
abline(v = (het$mean.het[het$state=="hscs"]-het$mean.het[het$state=="mono"]), lty="dashed", col="red3")


#---
#--- differential heterogeneity at a set of genomic regions
#---

#--- function to compute differential epigenetic heterogeneity between two groups
compute.diff.eICH <- function(group1, group2, region.type) {
  
  merged <- cbind(group1, group2)
  
  datasets <- list(group1=group1, group2=group2)
  
  for (i in 1:50) {
    sample.cells <- c(colnames(merged)) %>% sample(100)
    datasets[[paste0("mix", i)]] <- merged[,sample.cells]
  }
  
  het <- compute.eITH(datasets)
  het.dist <- het[het$state %notin% c("group1", "group2"),]
  
  dif.dist <- c()
  for (i in 1:nrow(het.dist)/2) {
    dif <- het.dist$mean.het[i]-het.dist$mean.het[i+1]
    dif.dist <- c(dif.dist, dif)
  }
  
  df <- data.frame(dif.dist=dif.dist, condition="mix")
  gg <- ggplot(df, aes(x=dif.dist)) + 
    geom_density() +
    lims(x=c(-1,1)) +
    labs(x="", subtitle = region.type) +
    geom_vline(aes(xintercept=(het$mean.het[het$state=="group1"]-het$mean.het[het$state=="group2"])), color="red3", linetype="dashed") +
    theme_bw()
  
  return(gg)
  
}


#--- apply above function to compute eICH per region between HSCs and monocytes
dim(hscs)
dim(mono)


# encode TFBS
lola <- get(load("/omics/groups/OE0219/internal/KatherineK/Genomic_Annotations/LOLA/hg19/encode_tfbs/encode_tfbs.RData"))
index <- read.table("/omics/groups/OE0219/internal/KatherineK/Genomic_Annotations/LOLA/hg19/encode_tfbs/index.txt", header = T)
names(lola) <- index$filename
lola <- lola[grepl("K562", names(lola))]

regions.interest <- names(lola)

data.gr <- atac@rowRanges
gg <- list()
for (region in regions.interest) {
  
  state <- paste0(index$cellType[index$filename==region], " - ", index$antibody[index$filename==region])
  print(state)
  
  # GRanges object containing selected regions
  temp <- lola[[region]]
  
  select.sites <- subsetByOverlaps(data.gr, temp) %>% names()
  if (length(select.sites)>1000) { select.sites <- sample(select.sites, 1000)}
  group1 <- hscs[select.sites, ]
  group2 <- mono[select.sites, ]
  
  print("computing differentials...")
  plot.dif <- compute.diff.eICH(group1, group2, state)
  
  gg[[state]] <- plot.dif
  
}

pdf("densityplots_per_region_diff_eICH_HSCs_monocytes.pdf", 10, 10)
ggarrange(plotlist = gg, ncol=4, nrow=4)
dev.off()
