
#--- compute epiCHAOS scores on different genomic regions


#--- the lola encode tfbs should be included as example data..
lola <- get(load("/omics/groups/OE0219/internal/KatherineK/Genomic_Annotations/LOLA/hg19/encode_tfbs/encode_tfbs.RData"))
index <- read.table("/omics/groups/OE0219/internal/KatherineK/Genomic_Annotations/LOLA/hg19/encode_tfbs/index.txt", header = T)
names(lola) <- index$filename
regions.interest <- index$filename
chrom.states <- lola

#--- load scATAC-seq data from normal hematopoiesis from Granja et al. 2019
atac <- readRDS("/omics/groups/OE0219/internal/KatherineK/data/scATAC/Granja2019/scATAC-Healthy-Hematopoiesis-191120.rds")

#--- to hold scATAC matrices for each celltype/TFBS pair
datasets <- list()

#--- Genomic Ranges for the scATAC peaks
data.gr <- atac@rowRanges

#--- so that the same number of peaks is selected for each region type
n <- 2000

#--- loop through cell types
for (celltype in unique(atac$BioClassification)) {
  
  #--- ids for the selected celltype, subset 100 cells
  cells <- atac[,atac$BioClassification==celltype] %>% colnames() %>% sample(100)
  
  #--- subset the peaks matrix for the selected cells
  hema <- atac@assays$data$counts[,cells]

  #--- binarise
  hema[hema>1] <- 1
  
  #--- assign rownames
  rownames(hema) <- names(atac@rowRanges)
  
  #--- loop through each encode TFBS
  for (region in regions.interest) {
    
    state <- paste0(index$cellType[index$filename==region], " - ", index$antibody[index$filename==region])
    print(region)
    
    #--- get GRanges object containing selected regions and subset the peaks matrix for the selected regions
    temp <- chrom.states[[region]]
    select.sites <- subsetByOverlaps(data.gr, temp) %>% names()
    datasets[[paste0(region, "-", celltype)]] <- hema[select.sites, ]
    
  }
  
  
  #--- downsample so that the same number of peaks are kept for each region type, e.g. 2000
  for (i in names(datasets)) {
    temp <- datasets[[i]]
    if (nrow(temp)<n) {
      
      #--- if wanted to remove region types with small number of peaks, otherwise let datasets[[i]] <- temp
      datasets[[i]] <- NULL  
    } else {
      datasets[[i]] <- temp[sample(nrow(temp), n, replace = F), ]
      
    }
  }
  
}

#--- compute and save epiCHAOS scores
het <- compute.eITH(datasets)

#--- adjust names for plotting
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

#--- plot ordered per-region epiCHAOS scores, highlighting PRC2 targets
# gg <- list()
# for (celltype in unique(temp$celltype)) {
#   
#   #--- exclude celltypes annotated as "Unknown"
#   if (grepl("Unk", celltype)) { next }
#   gg[[celltype]] <- ggplot(temp[temp$celltype==celltype,], aes(y=mean.het, x=reorder(state, mean.het), color=label)) +
#     labs( x="", y="", subtitle = celltype) +
#     geom_point(size=0.1)+
#     theme_classic() +
#     lims(y=c(min(temp$mean.het),max(temp$mean.het)))+
#     theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
# }
# 
# length(gg)
# ggarrange(plotlist = gg, ncol=12, legend = F)


#--- correlate for each celltype pair
temp <- temp[!grepl("Unk", temp$celltype),]
cor.cells <- data.frame(row.names = unique(temp$celltype))
for (state in unique(temp$regionset)) {
  for (celltype in unique(temp$celltype)) {
    cor.cells[celltype,state] <- temp$mean.het[temp$regionset==state & temp$celltype==celltype]
  }
}

dim(cor.cells)

#--- heatmap for figure S5
pdf("correlation_heatmap_compare_per-region_eICH_per_celltype.pdf", 10,10)
cor(t(cor.cells)) %>% pheatmap(color = brewer.pal(9, "Blues"), cellheight = 17, cellwidth = 17, display_numbers = T, number_color = "black", fontsize_number = 7, border_color = "grey90")
dev.off()
