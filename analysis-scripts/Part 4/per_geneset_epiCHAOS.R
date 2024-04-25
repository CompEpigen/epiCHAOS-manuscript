

#--- compare epiCHAOS scores on different gene sets

library(msigdbr)
library(magrittr)
library(dplyr)
library(stringr)

#--- scATAC-seq data from normal hematopoiesis from Granja et al. 2019
atac <- readRDS("/omics/groups/OE0219/internal/KatherineK/data/scATAC/Granja2019/scATAC-Healthy-Hematopoiesis-191120.rds")
atac <- atac[,grepl("BMMC", atac$Group)]

#--- select hg19 promoters
require("TxDb.Hsapiens.UCSC.hg19.knownGene")
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
promoters <- promoters(genes(txdb), upstream = 1500, downstream = 500)

#--- load hallmark gene ontologies for human from msigdb
gene_sets = msigdbr(species = "Homo sapiens")
gene_sets <- gene_sets[gene_sets$gs_subcat=="GO:BP",]

set.names <- gene_sets$gs_name %>% unique()

gene_sets_all = msigdbr(species = "Homo sapiens")

#--- housekeeping genes
hk <- c("Actb", "Atp5f1", "Atp5pb", "B2m", "Gapdh", "Hprt1", "Hprt", "Pgk1", "Rer1", "Rpl13a", "Rpl27", "Sdha", "Tbp", "Ubc") %>% str_to_upper()
hk <- gene_sets_all$entrez_gene[gene_sets_all$gene_symbol %in% hk] %>% unique()

#--- bivalent genes
bivalent <- read.table("/omics/groups/OE0219/internal/KatherineK/Genomic_Annotations/gene_annotations/bivalent_genes_hg19_Court_2017.txt", header = T, sep = "\t")$Gene
bivalent <- bivalent %>% str_split(",") %>% unlist() %>% unique()

#--- use gene sets table to translate bivalent gene ids to entrez ids
bivalent <- gene_sets_all$entrez_gene[gene_sets_all$gene_symbol %in% bivalent] %>% unique()
bivalent <- promoters[names(promoters) %in% bivalent,] %>% names()

#--- to hold peaks matrices for each selected cell type/gene set pair
datasets <- list()

#--- loop through cell types
for (celltype in unique(atac$BioClassification)) {
  
  #--- if want to exclude cell types with few cells, otherwise silence this
  if (length(colnames(atac[,atac$BioClassification==celltype]))<100) { next }
  print(celltype)
  
  #--- select ids corresponding to the cell type of interest, subset 100 cells, subset the peaks matrix for those cells
  cell.ids <- atac[,atac$BioClassification==celltype] %>% colnames() %>% sample(100)
  hema <- atac@assays$data$counts[,cell.ids]
  
  #--- binarise & update rownames
  hema[hema>1] <- 1
  rownames(hema) <- names(atac@rowRanges)
  
  #--- loop through gene sets
  for (set in set.names) {
    
    print(set)
    
    #--- select promoters for gene set of interest
    genes <- gene_sets$entrez_gene[gene_sets$gs_name==set] %>% intersect(names(promoters))
    select.promoters <- promoters[genes,]
    select.promoters <- subsetByOverlaps(atac@rowRanges, select.promoters) %>% names()
    
    #--- exclude gene sets with small number of associated genes/peaks
    if (length(select.promoters)<20) { next }
    
    datasets[[paste0(celltype, "-", set)]] <- hema[select.promoters,]
    
  }
  
  #--- if wanted, add gene set for bivalent genes
  select.promoters <- promoters[bivalent,]
  select.promoters <- subsetByOverlaps(atac@rowRanges, select.promoters) %>% names()
  datasets[[paste0(celltype, "-BIV")]] <- hema[select.promoters,]
  
  #--- if wanted, add gene set for housekeeping genes
  select.promoters <- promoters[names(promoters) %in% hk,]
  select.promoters <- subsetByOverlaps(atac@rowRanges, select.promoters) %>% names()
  datasets[[paste0(celltype, "-HK")]] <- hema[select.promoters,]
  
}

lapply(datasets, dim)


#--- compute and save epiCHAOS scores
het <- compute.eITH(datasets)

#saveRDS(het, "epiCHAOS_GOBP_Hematopoiesis.Rds")

###

#--- load epiCHAOS scores for downstream analysis
het <- readRDS("/omics/groups/OE0219/internal/KatherineK/ATACseq/per-region-eITH/epiCHAOS_GOBP_Hematopoiesis.Rds")

#--- adjust labels for plotting
het$celltype <- het$state %>% str_split("-") %>% lapply("[", 1) %>% unlist()
het$geneset <- het$state %>% str_split("-") %>% lapply("[", 2) %>% unlist()
het$label <- ifelse(het$geneset %in% c("BIV"), "Bivalent genes", "")
het$label[het$geneset=="GOBP_POSITIVE_REGULATION_OF_CELL_FATE_COMMITMENT"] <- "Positive regulation of cell fate commitment"
het$label[het$geneset=="GOBP_CELL_FATE_SPECIFICATION"] <- "Cell fate specification"

#--- plot of ordered per-gene set epiCHAOS scores in HSCs
ggplot(het[het$celltype=="01_HSC",], aes(y=mean.het^2, x=reorder(state, mean.het), label=label)) +
  labs(y="epiCHAOS", x="Gene set") +
  geom_point(size=0.5, color="grey30")+
  theme_classic() +
  expand_limits(x= c(-100, 100))+
  ggrepel::geom_text_repel(size=4, max.overlaps=6000)+
  lims(y=c(min(het$mean.het),max(het$mean.het)))+
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())


#--- investigate differences in epiCHAOS ranks between celltypes for gene sets of interest (Figure 4 and Figure S6)
het <- readRDS("/omics/groups/OE0219/internal/KatherineK/ATACseq/per-region-eITH/epiCHAOS_GOBP_Hematopoiesis.Rds")
het$celltype <- het$state %>% str_split("-") %>% lapply("[", 1) %>% unlist()
het$geneset <- het$state %>% str_split("-") %>% lapply("[", 2) %>% unlist()

#--- for a given gene set, for each celltype, find what is the rank of that gene set among all gene sets for that celltype
genesets <- c("GOBP_CELL_FATE_SPECIFICATION", "BIV")

ranks <- data.frame()
for (geneset in genesets) {
  for (celltype in unique(het$celltype)) {
    temp <- het[het$celltype==celltype,c("mean.het", "geneset")] %>% arrange(mean.het) %>% pull(geneset) %>% rev()
    ranks[celltype, geneset] <- which(temp==geneset)/length(unique(het$geneset))
    
  }
}

#--- adjust named for plotting
ranks$celltype <- rownames(ranks) %>% str_split("_") %>% lapply("[", 2) %>% str_replace_all("\\.", " ")
ranks <- ranks[!grepl("Unk", ranks$celltype),]

#--- order so that progenitors are first in the barplot
order <- c("01_HSC",  "05_CMP.LMPP" , "06_CLP.1" , "15_CLP.2" , "07_GMP" ,  "02_Early.Eryth" , "16_Pre.B", "04_Early.Baso","08_GMP.Neut" ,  "11_CD14.Mono.1" , "12_CD14.Mono.2" ,"10_cDC",  "09_pDC" , "03_Late.Eryth",  "17_B" , "25_NK" ,   "20_CD4.N1",  "19_CD8.N" , "23_CD8.EM", "22_CD4.M" , "24_CD8.CM" )       
ranks <- ranks[order,]                             
ranks$order <- 1:nrow(ranks)

#--- barplots
barplot1 <- ggplot(ranks, aes(x=-log10(GOBP_CELL_FATE_SPECIFICATION), y=reorder(celltype, -GOBP_CELL_FATE_SPECIFICATION), fill=-log10(GOBP_CELL_FATE_SPECIFICATION))) +
  geom_bar(alpha=0.5, stat = "identity", width=0.4) +
  labs(x="Rank", y="",subtitle="Cell fate specification") +
  scale_fill_gradient(low="darkblue", high="red3") +
  theme_classic() +
  NoLegend()

barplot2 <- ggplot(ranks, aes(x=-log10(BIV), y=reorder(celltype, -BIV), fill=-log10(BIV))) +
  geom_bar(alpha=0.5, stat = "identity", width=0.4) +
  labs(x="Rank", y="",subtitle="Bivalent genes") +
  scale_fill_gradient(low="darkblue", high="red3") +
  theme_classic() +
  NoLegend()

ggarrange(barplot1, barplot2)



#--- correlate per-gene set epiCHAOS scores for each celltype pair (for Figure S5)
temp <- het[!grepl("Unk", het$celltype),]
cor.cells <- data.frame(row.names = unique(temp$celltype))
for (state in unique(temp$geneset)) {
  for (celltype in unique(temp$celltype)) {
    cor.cells[celltype,state] <- temp$mean.het[temp$geneset==state & temp$celltype==celltype]
  }
}

dim(cor.cells)

#--- heatmap for Figure S5
pdf("correlation_heatmap_compare_per-geneset_eICH_per_celltype.pdf", 10,10)
cor(t(cor.cells)) %>% pheatmap(color = brewer.pal(9, "Blues"), cellheight = 17, cellwidth = 17, display_numbers = T, number_color = "black", fontsize_number = 7, border_color = "grey90")
dev.off()
