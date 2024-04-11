

#---
#--- compute eITH on different gene sets
#---

source("/omics/groups/OE0219/internal/KatherineK/ATACseq/eITH-test-scripts/jaccard.R")

library(msigdbr)
library(magrittr)
library(dplyr)
library(stringr)


setwd("/omics/groups/OE0219/internal/KatherineK/ATACseq/per-region-eITH")

#--- scATAC-seq example dataset - data from normal hematopoiesis from Granja et al. 2019 - subset for HSCs
atac <- readRDS("/omics/groups/OE0219/internal/KatherineK/data/scATAC/Granja2019/scATAC-Healthy-Hematopoiesis-191120.rds")
atac <- atac[,grepl("BMMC", atac$Group)]

#--- select hg19 promoters
require("TxDb.Hsapiens.UCSC.hg19.knownGene")
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
promoters <- promoters(genes(txdb), upstream = 1500, downstream = 500)

#--- gene sets

# hallmark gene ontologies for human
gene_sets = msigdbr(species = "Homo sapiens")
gene_sets <- gene_sets[gene_sets$gs_subcat=="GO:BP",]

set.names <- gene_sets$gs_name %>% unique()

gene_sets_all = msigdbr(species = "Homo sapiens")

# hk genes
hk <- c("Actb", "Atp5f1", "Atp5pb", "B2m", "Gapdh", "Hprt1", "Hprt", "Pgk1", "Rer1", "Rpl13a", "Rpl27", "Sdha", "Tbp", "Ubc") %>% str_to_upper()
hk <- gene_sets_all$entrez_gene[gene_sets_all$gene_symbol %in% hk] %>% unique()

#--- bivalent genes
bivalent <- read.table("/omics/groups/OE0219/internal/KatherineK/Genomic_Annotations/gene_annotations/bivalent_genes_hg19_Court_2017.txt", header = T, sep = "\t")$Gene
bivalent <- bivalent %>% str_split(",") %>% unlist() %>% unique()

# use gene sets table to translate bivalent gene ids to entrez ids
bivalent <- gene_sets_all$entrez_gene[gene_sets_all$gene_symbol %in% bivalent] %>% unique()
bivalent <- promoters[names(promoters) %in% bivalent,] %>% names()


#--- selected cell type
datasets <- list()
for (celltype in unique(atac$BioClassification)) {
  
  if (length(colnames(atac[,atac$BioClassification==celltype]))<100) { next }
  print(celltype)
  
  hscs <- atac[,atac$BioClassification==celltype] %>% colnames() %>% sample(100)
  hema <- atac@assays$data$counts[,hscs]
  hema[hema>1] <- 1
  rownames(hema) <- names(atac@rowRanges)
  
  for (set in set.names) {
    
    print(set)
    
    # select promoters for gene set of interest
    genes <- gene_sets$entrez_gene[gene_sets$gs_name==set] %>% intersect(names(promoters))
    select.promoters <- promoters[genes,]
    select.promoters <- subsetByOverlaps(atac@rowRanges, select.promoters) %>% names()
    if (length(select.promoters)<20) { next }
    
    datasets[[paste0(celltype, "-", set)]] <- hema[select.promoters,]
    
  }
  
  select.promoters <- promoters[bivalent,]
  select.promoters <- subsetByOverlaps(atac@rowRanges, select.promoters) %>% names()
  datasets[[paste0(celltype, "-BIV")]] <- hema[select.promoters,]
  
  select.promoters <- promoters[names(promoters) %in% hk,]
  select.promoters <- subsetByOverlaps(atac@rowRanges, select.promoters) %>% names()
  datasets[[paste0(celltype, "-HK")]] <- hema[select.promoters,]
  
}

lapply(datasets, dim)


# compute eITH
het <- compute.eITH(datasets)

saveRDS(het, "epiCHAOS_GOBP_Hematopoiesis.Rds")


#--- downstream
het <- readRDS("/omics/groups/OE0219/internal/KatherineK/ATACseq/per-region-eITH/epiCHAOS_GOBP_Hematopoiesis.Rds")

het$celltype <- het$state %>% str_split("-") %>% lapply("[", 1) %>% unlist()
het$geneset <- het$state %>% str_split("-") %>% lapply("[", 2) %>% unlist()
unique(het$celltype)
plot(het$mean.het[het$celltype=="01_HSC"], het$mean.het[het$celltype=="25_NK"])

het[het$celltype=="01_HSC",1:2] %>% arrange(mean.het) %>% tail(20)

# temp <- data.frame(row.names = unique(het$celltype))
# for (i in unique(het$celltype)) {
#   for (g in unique(het$geneset)) {
#     temp[i,g] <- het$mean.het[het$celltype==i & het$geneset==g]
#   }
# }

het$label <- ifelse(het$geneset %in% c("BIV"), "Bivalent genes", "")
het$label[het$geneset=="GOBP_POSITIVE_REGULATION_OF_CELL_FATE_COMMITMENT"] <- "Positive regulation of cell fate commitment"
het$label[het$geneset=="GOBP_CELL_FATE_SPECIFICATION"] <- "Cell fate specification"

ggplot(het[het$celltype=="17_B",], aes(y=mean.het^2, x=reorder(state, mean.het), label=label)) +
  labs(y="epiCHAOS", x="Gene set") +
  geom_point(size=0.5, color="grey30")+
  theme_classic() +
  expand_limits(x= c(-100, 100))+
  ggrepel::geom_text_repel(size=4, max.overlaps=6000)+
  lims(y=c(min(het$mean.het),max(het$mean.het)))+
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

a <- het[grepl("Bivalent", het$label),c("celltype", "mean.het")]
b <- het[grepl("Positive", het$label),c("celltype", "mean.het")]
c <- het[grepl("specification", het$label),c("celltype", "mean.het")]

a$celltype==b$celltype
temp <- data.frame(row.names = a$celltype, A=a$mean.het, B=b$mean.het, C=c$mean.het)
pheatmap(temp, scale = "column")
#pheatmap(temp, show_colnames = F, show_rownames = F, scale = "column")


# pdf("hockeyplot_epiCHAOS_GOBP_HSCs.pdf", 5, 3.5)
# ggplot(het, aes(y=mean.het^2, x=reorder(state, mean.het), label=label)) +
#   labs(y="epiCHAOS", x="Gene set") +
#   geom_point(size=0.5, color="grey30")+
#   theme_classic() +
#   expand_limits(x= c(-100, 100))+
#   ggrepel::geom_text_repel(size=4, max.overlaps=600)+
#   lims(y=c(min(het$mean.het),max(het$mean.het)))+
#   theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
# dev.off()


#--- for a given gene set, for each celltype, find what is the rank of that gene set among all gene sets for that celltype

geneset <- "GOBP_POSITIVE_REGULATION_OF_CELL_FATE_COMMITMENT"
geneset <- "GOBP_CELL_FATE_SPECIFICATION"
geneset <- "GOBP_SOMATIC_STEM_CELL_DIVISION"

ranks <- data.frame()
for (celltype in unique(het$celltype)) {
  temp <- het[het$celltype==celltype,c("mean.het", "geneset")] %>% arrange(mean.het) %>% pull(geneset) %>% rev()
  ranks[celltype, geneset] <- which(temp==geneset)/length(unique(het$geneset))
  
}
# ranks[1:10,1:10]
# 
# ranks[is.na(ranks)] <- 0
# pheatmap(ranks)

ranks$celltype <- rownames(ranks) %>% str_split("_") %>% lapply("[", 2) %>% str_replace_all("\\.", " ")
ranks <- ranks[!grepl("Unk", ranks$celltype),]
p1 <- ggplot(ranks, aes(x=-log10(GOBP_CELL_FATE_SPECIFICATION), y=celltype, color=-log10(GOBP_CELL_FATE_SPECIFICATION))) +
  geom_point(alpha=0.5, shape=18, size=2) +
  labs(x="Rank", y="",subtitle="Cell fate specification") +
  scale_color_gradient(low="darkblue", high="red3") +
  theme_bw() +
  NoLegend()

p2 <- ggplot(ranks, aes(x=-log10(GOBP_POSITIVE_REGULATION_OF_CELL_FATE_COMMITMENT), y=celltype, color=-log10(GOBP_POSITIVE_REGULATION_OF_CELL_FATE_COMMITMENT))) +
  geom_point(alpha=0.5, shape=18, size=2) +
  labs(x="Rank", y="",subtitle="Positive regulation of cell fate commitment") +
  scale_color_gradient(low="darkblue", high="red3") +
  theme_bw() +
  NoLegend()

p3 <- ggplot(ranks, aes(x=-log10(GOBP_SOMATIC_STEM_CELL_DIVISION), y=celltype, color=-log10(GOBP_SOMATIC_STEM_CELL_DIVISION))) +
  geom_point(alpha=0.5, shape=18, size=2) +
  labs(x="Rank", y="",subtitle="Somatic stem cell division") +
  scale_color_gradient(low="darkblue", high="red3") +
  theme_bw() +
  NoLegend()

pdf("dotplots_selected_GOBP_compare_celltypes.pdf", 7, 4)
ggpubr::ggarrange(p1,p2, p3, ncol=3)
dev.off()

# order so that progenitors are first
order <- c("01_HSC",  "05_CMP.LMPP" , "06_CLP.1" , "15_CLP.2" , "07_GMP" ,  "02_Early.Eryth" , "16_Pre.B", "04_Early.Baso","08_GMP.Neut" ,  "11_CD14.Mono.1" , "12_CD14.Mono.2" ,"10_cDC",  "09_pDC" , "03_Late.Eryth",  "17_B" , "25_NK" ,   "20_CD4.N1",  "19_CD8.N" , "23_CD8.EM", "22_CD4.M" , "24_CD8.CM" )       
ranks <- ranks[order,]                             
ranks$order <- 1:nrow(ranks)

#--- or same as barplots
p1 <- ggplot(ranks, aes(x=-log10(GOBP_CELL_FATE_SPECIFICATION), y=reorder(celltype, -GOBP_CELL_FATE_SPECIFICATION), fill=-log10(GOBP_CELL_FATE_SPECIFICATION))) +
  geom_bar(alpha=0.5, stat = "identity", width=0.4) +
  labs(x="Rank", y="",subtitle="Cell fate specification") +
  scale_fill_gradient(low="darkblue", high="red3") +
  theme_bw() +
  NoLegend()

p2 <- ggplot(ranks, aes(x=-log10(GOBP_POSITIVE_REGULATION_OF_CELL_FATE_COMMITMENT), y=reorder(celltype, -GOBP_POSITIVE_REGULATION_OF_CELL_FATE_COMMITMENT), fill=-log10(GOBP_POSITIVE_REGULATION_OF_CELL_FATE_COMMITMENT))) +
  geom_bar(alpha=0.5, stat = "identity", width=0.4) +
  labs(x="Rank", y="",subtitle="Positive regulation of cell fate commitment") +
  scale_fill_gradient(low="darkblue", high="red3") +
  theme_bw() +
  NoLegend()

p3 <- ggplot(ranks, aes(x=-log10(GOBP_SOMATIC_STEM_CELL_DIVISION), y=reorder(celltype, -GOBP_SOMATIC_STEM_CELL_DIVISION), fill=-log10(GOBP_SOMATIC_STEM_CELL_DIVISION))) +
  geom_bar(alpha=0.5, stat = "identity", width=0.4) +
  labs(x="Rank", y="",subtitle="Somatic stem cell division") +
  scale_fill_gradient(low="darkblue", high="red3") +
  theme_bw() +
  NoLegend()

pdf("barplots_selected_GOBP_compare_celltypes.pdf", 7, 4)
ggpubr::ggarrange(p1,p2, p3, ncol=3)
dev.off()


#--- correlate for each celltype pair
temp <- het[!grepl("Unk", het$celltype),]
cor.cells <- data.frame(row.names = unique(temp$celltype))
for (state in unique(temp$geneset)) {
  for (celltype in unique(temp$celltype)) {
    cor.cells[celltype,state] <- temp$mean.het[temp$geneset==state & temp$celltype==celltype]
  }
}

dim(cor.cells)

pdf("correlation_heatmap_compare_per-geneset_eICH_per_celltype.pdf", 10,10)
cor(t(cor.cells)) %>% pheatmap(color = brewer.pal(9, "Blues"), cellheight = 17, cellwidth = 17, display_numbers = T, number_color = "black", fontsize_number = 7, border_color = "grey90")
dev.off()
