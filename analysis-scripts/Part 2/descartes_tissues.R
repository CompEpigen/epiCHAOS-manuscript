
#--- script for computation of epiCHAOS scores per tissue and per celltype using the human pan-tissue scATAC atlas from https://descartes.brotmanbaty.org/bbi/human-chromatin-during-development/ (relating to figure S2)

#--- required packages
library(Seurat)
library(magrittr)

#--- data were downloaded in .rds files, one per tissue. Since the data are large, below we collect a list where each element contains 100 cells of each cell type and save those for downstream analysis
# datasets <- list()
# for (i in list.files("/omics/groups/OE0219/internal/KatherineK/data/scATAC/pan-tissue-descartes/", pattern = "seurat")) {
# 
#   obj <- readRDS(paste0("/omics/groups/OE0219/internal/KatherineK/data/scATAC/pan-tissue-descartes/", i))
#   obj <- UpdateSeuratObject(obj)
#   meta <- obj@meta.data
# 
#   for (celltype in unique(meta$cell_type)) {
#     ids <- meta[meta$cell_type==celltype, ] %>% rownames()
#     print(celltype)
#     datasets[[paste0(i, "-", celltype)]] <- obj@assays$peaks$counts[,ids]
#     #if (length(ids)>1000) {}
#   }
# 
# }
# 
# saveRDS(datasets, "/omics/groups/OE0219/internal/KatherineK/data/scATAC/pan-tissue-descartes/atac_list_celltypes.Rds")

#--- load the datasets from above - a list of scATAC-seq matrices, one for each cell/tissue type
datasets <- readRDS("/omics/groups/OE0219/internal/KatherineK/data/scATAC/pan-tissue-descartes/atac_list_celltypes.Rds")

#--- subset 50,000 rows (peaks) for ease of computation
select.rows <- sample(nrow(datasets[[1]]), 50000)
for (i in names(datasets)) {
  temp <- datasets[[i]]
  temp <- temp[select.rows,1:min(c(ncol(temp), 100))]
  
  #--- binarise the data
  temp[temp>1] <- 1
  datasets[[i]] <- temp
}

length(datasets)
lapply(datasets, dim)

#--- compute epiCHAOS scores
het <- compute.eITH(datasets)

#--- check result
het %>% arrange(mean.het)

#--- add labels and adjust names for plotting
het$label <- ""
het$label[grepl("placenta", het$state)] <-  "placenta"
het$label[grepl("cerebellum", het$state)] <-  "neural"
het$label[grepl("cerebrum", het$state)] <-  "neural"
het <- het[!grepl("Unknown", het$state),]
het$state <- het$state %>% str_remove("_filtered.seurat.for_website.RDS")

#--- save scores
saveRDS(het, file = "epiCHAOS_scores_per_celltype.Rds")

#--- save plots
svg("/omics/groups/OE0219/internal/KatherineK/ATACseq/epiCHAOS-Figures/Supplementary/S2/barplot_epiCHAOS_descartes_percelltype.svg", 10, 4)
ggplot(het, aes(y=mean.het, x=reorder(state, mean.het), color=label, fill=label)) +
  geom_bar(size=1,  stat="identity", width=0.3)+
  labs( x="", y="epiCHAOS")+
  scale_color_manual(values = c("grey70", "steelblue4", "red3")) +
  scale_fill_manual(values = c("grey70", "steelblue4", "red3")) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=5))
dev.off()

svg("/omics/groups/OE0219/internal/KatherineK/ATACseq/epiCHAOS-Figures/Supplementary/S2/barplot_epiCHAOS_descartes_pertissue.svg", 6, 4)
ggplot(het, aes(y=mean.het, x=reorder(state, mean.het), color=label, fill=label)) +
  geom_bar(size=1,  stat="identity", width=0.3)+
  labs( x="", y="epiCHAOS")+
  scale_color_manual(values = c("grey70", "steelblue4", "red3")) +
  scale_fill_manual(values = c("grey70", "steelblue4", "red3")) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=10))
dev.off()
