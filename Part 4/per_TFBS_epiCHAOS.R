

#---
#--- compute eITH on different genomic regions in each hematopoietic celltype
#---

library(magrittr)
library(dplyr)
library(stringr)
library(GenomicRanges)


setwd("/omics/groups/OE0219/internal/KatherineK/ATACseq/per-region-eITH")
source("/omics/groups/OE0219/internal/KatherineK/ATACseq/eITH-test-scripts/jaccard.R")

#--- encode TFBS
lola <- get(load("/omics/groups/OE0219/internal/KatherineK/Genomic_Annotations/LOLA/hg19/encode_tfbs/encode_tfbs.RData"))
index <- read.table("/omics/groups/OE0219/internal/KatherineK/Genomic_Annotations/LOLA/hg19/encode_tfbs/index.txt", header = T)
names(lola) <- index$filename
regions.interest <- index$filename
chrom.states <- lola

#--- scATAC-seq example dataset - data from normal hematopoiesis from Granja et al. 2019 - subset for HSCs
atac <- readRDS("/omics/groups/OE0219/internal/KatherineK/data/scATAC/Granja2019/scATAC-Healthy-Hematopoiesis-191120.rds")

datasets <- list()
data.gr <- atac@rowRanges
n <- 2000

for (celltype in unique(atac$BioClassification)) {
  
  cells <- atac[,atac$BioClassification==celltype] %>% colnames() %>% sample(100)
  hema <- atac@assays$data$counts[,cells]
  dim(hema)
  hema[hema>1] <- 1
  rownames(hema) <- names(atac@rowRanges)
  #hema <- hema[rowSums(hema)>0 & rowSums(hema)<30,]
  
  for (region in regions.interest) {
    
    state <- paste0(index$cellType[index$filename==region], " - ", index$antibody[index$filename==region])
    print(region)
    
    # GRanges object containing selected regions
    temp <- chrom.states[[region]]
    
    select.sites <- subsetByOverlaps(data.gr, temp) %>% names()
    datasets[[paste0(region, "-", celltype)]] <- hema[select.sites, ]
    
  }
  
  
  # downsample because larger nrow datasets are percieved as being more homogenous
  for (i in names(datasets)) {
    temp <- datasets[[i]]
    if (nrow(temp)<n) {
      datasets[[i]] <- NULL
    } else {
      datasets[[i]] <- temp[sample(nrow(temp), n, replace = F), ]
      
    }
  }
  
  
}

# compute eITH
het <- compute.eITH(datasets)

saveRDS(het, "epiCHAOS_encodeTFBS_HSCs.Rds")

# calculate average eITH score per population
for (i in names(datasets)) {
  het$mean.het[het$state==i] <- mean(het$het[het$state==i])
}

het$celltype <- het$state %>% str_split("-") %>% lapply("[", 2) %>% unlist()
temp <- unique(het[,c("mean.het", "state", "celltype")])
head(temp)

saveRDS(temp, file = "per_region_eICH_per_hemato_celltype.Rds")


#--- downstream analysis
temp <- readRDS("/omics/groups/OE0219/internal/KatherineK/ATACseq/per-region-eITH/per_region_eICH_per_hemato_celltype.Rds")
head(temp)
temp$regionset <- temp$state %>% str_split("-") %>% lapply("[", 1) %>% unlist()
temp$label <- NULL
temp$label[grepl("Ezh2", temp$regionset)|grepl("Suz12", temp$regionset)] <- "PRC2"

gg <- list()
for (celltype in unique(temp$celltype)) {
  
  if (grepl("Unk", celltype)) { next }
  gg[[celltype]] <- ggplot(temp[temp$celltype==celltype,], aes(y=mean.het, x=reorder(state, mean.het), color=label)) +
    labs( x="", y="", subtitle = celltype) +
    geom_point(size=0.1)+
    theme_classic() +
    #geom_text_repel(size=1, max.overlaps=400)+
    lims(y=c(min(temp$mean.het),max(temp$mean.het)))+
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
}

length(gg)

pdf("eICH_per_TFBS_per_hemato_highlight_PRC2.pdf", 15, 2)
ggarrange(plotlist = gg, ncol=12, legend = F)
dev.off()


#--- correlate for each celltype pair
temp <- temp[!grepl("Unk", temp$celltype),]
cor.cells <- data.frame(row.names = unique(temp$celltype))
for (state in unique(temp$regionset)) {
  for (celltype in unique(temp$celltype)) {
    cor.cells[celltype,state] <- temp$mean.het[temp$regionset==state & temp$celltype==celltype]
  }
}

dim(cor.cells)

pdf("correlation_heatmap_compare_per-region_eICH_per_celltype.pdf", 10,10)
cor(t(cor.cells)) %>% pheatmap(color = brewer.pal(9, "Blues"), cellheight = 17, cellwidth = 17, display_numbers = T, number_color = "black", fontsize_number = 7, border_color = "grey90")
dev.off()
