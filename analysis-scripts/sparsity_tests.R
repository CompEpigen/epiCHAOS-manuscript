
#--- try using the scTAM-seq data to test the effect of increasing sparsity on epiCHAOS score

library(Seurat)
library(pheatmap)

set.seed(10)

#--- load scTAM data as seurat object downloaded from figshare
scTAM <- readRDS("/omics/groups/OE0219/internal/KatherineK/data/scTAM/seurat_for_figshare.rds")

#--- extract single-cell methylation data as matrix
mat <- scTAM@assays$DNAm@counts %>% as.matrix()
mat[1:10,1:10]

celltypes <- c("MPP4", "GMP", "HSC/MPP1",  "MPP3", "MPP2")
celltypes <- scTAM@meta.data$CellType %>%

#--- create matrix per cell type for epiCHAOS calculation
datasets <- list()
for (i in celltypes) {
  ids <- scTAM@meta.data[scTAM@meta.data$CellType==i,] %>% rownames()
  if (length(ids)>200) { 
    ids <- ids %>% sample(200, replace = F)
  }
  datasets[[i]] <- mat[,ids]
}

lapply(datasets, dim)

#--- compute epiCHAOS scores
het <- compute.eITH(datasets)

#--- function to introduce sparsity into 
sparsify <- function(x, percent) {
  x.sparse <- list()
  
  for (i in names(x)) {
    temp <- x[[i]]
    n <- round(nrow(temp)*ncol(temp)*percent)
    temp[sample(n)] <- 0
    x.sparse[[i]] <- temp
  }
  return(x.sparse)
}

#--- 20%
datasets.10 <- sparsify(datasets, 0.1)
lapply(datasets.10, sum)  
het.10 <- compute.eITH(datasets.10)
plot(het$het.raw, het.10$het.raw)

#--- 20%
datasets.20 <- sparsify(datasets, 0.2)
lapply(datasets.20, sum)
het.20 <- compute.eITH(datasets.20)
plot(het$het.raw, het.20$het.raw)

#--- 30%
datasets.30 <- sparsify(datasets, 0.3)
lapply(datasets.30, sum)  
het.30 <- compute.eITH(datasets.30)
plot(het$het.raw, het.30$het.raw)

#--- 40%
datasets.40 <- sparsify(datasets, 0.4)
lapply(datasets.40, sum)
het.40 <- compute.eITH(datasets.40)
plot(het$het.raw, het.40$het.raw)

#--- 50%
datasets.50 <- sparsify(datasets, 0.5)
lapply(datasets.50, sum)
het.50 <- compute.eITH(datasets.50)
plot(het$het.raw, het.50$het.raw)

# cbind(het$het.raw, het.10$het.raw, het.20$het.raw, het.30$het.raw, het.40$het.raw, het.50$het.raw) %>% cor() %>% 
#   pheatmap(cluster_rows = F, cluster_cols = F, display_numbers = T, number_color = "red3", 
#            fontsize_number = 15, number_format = "%.3f", border_color = "grey90")

temp <- cbind(het$mean.het, het.10$mean.het, het.20$mean.het, het.30$mean.het, het.40$mean.het, het.50$mean.het) 
rownames(temp) <- het$state
colnames(temp) <- c("0%", "10%", "20%", "30%", "40%", "50%")
order <- het %>% arrange(mean.het) %>% pull(state)
p1 <- temp[order,] %>% pheatmap(cluster_rows = F, cluster_cols = F,  border_color = "grey90", cellwidth = 15, cellheight = 15, color = rev(brewer.pal(10, "RdBu")))

pdf("/omics/groups/OE0219/internal/KatherineK/ATACseq/QualityControl/Sparsity tests/heatmap_scTAM_sparsity_test.pdf")
p1
dev.off()


#--- try in another dataset --> scATAC-seq hematopoietic data

#--- load scATAC from normal hematopoiesis from Granja et al. 2019, subset for bone marrow (BM)-derived
atac <- readRDS("/omics/groups/OE0219/internal/KatherineK/data/scATAC/Granja2019/scATAC-Healthy-Hematopoiesis-191120.rds")
atac <- atac[,grepl("BMMC", atac@colData$Group)]

#--- counts matrix
hema <- atac@assays$data$counts # %>% as.matrix()

#--- create matrices per cell type for epiCHAOS calculation
datasets <- list()
n <- 100 # select n cells form each celltype
n.row <- 50000
for (i in unique(atac$BioClassification[!grepl("Unk", atac$BioClassification)])) { # exclude cells annotated as "Unknown"
  cellnames <- colnames(atac)[atac$BioClassification==i]
  if (length(cellnames)>n) { cellnames <- sample(cellnames, n, replace = F) }
  datasets[[i]] <- hema[sample(nrow(hema), n.row),colnames(hema) %in% cellnames]
  datasets[[i]][datasets[[i]]>1] <- 1 # binarise data
}

lapply(datasets, dim)

#--- compute and save epiCHAOS scores
het <- compute.eITH(datasets)

#--- introduce sparsity
datasets.10 <- sparsify(datasets, 0.1)
datasets.20 <- sparsify(datasets, 0.2)
datasets.30 <- sparsify(datasets, 0.3)
datasets.40 <- sparsify(datasets, 0.4)
datasets.50 <- sparsify(datasets, 0.5)

#--- compute epiCHAOS scores in data of varying sparsity
het.10 <- compute.eITH(datasets.10)
het.20 <- compute.eITH(datasets.20)
het.30 <- compute.eITH(datasets.30)
het.40 <- compute.eITH(datasets.40)
het.50 <- compute.eITH(datasets.50)

#--- compare the resulting heterogeneity scores from each test
temp <- cbind(het$mean.het, het.10$mean.het, het.20$mean.het, het.30$mean.het, het.40$mean.het, het.50$mean.het) 
rownames(temp) <- het$state
colnames(temp) <- c("0%", "10%", "20%", "30%", "40%", "50%")
order <- het %>% arrange(mean.het) %>% pull(state)
p1 <- temp[order,] %>% pheatmap(cluster_rows = F, cluster_cols = F,  border_color = "grey90", cellwidth = 15, cellheight = 15, color = rev(brewer.pal(10, "RdBu")))

pdf("/omics/groups/OE0219/internal/KatherineK/ATACseq/QualityControl/Sparsity tests/heatmap_scATAC_hemato_sparsity_test.pdf")
p1
dev.off()

plot(temp[,1], temp[,6])

temp <- data.frame(temp)
pdf("/omics/groups/OE0219/internal/KatherineK/ATACseq/QualityControl/Sparsity tests/scatterplot_scATAC_hemato_original_score_vs_added_sparsity.pdf", 5, 4)
ggplot(temp, aes(x=X0., y=X50.)) +
  geom_point(size=3, alpha=0.7)+
  lims(y=c(0, 1)) +
  geom_smooth(color="grey30", method = "lm", se=F)+
  stat_cor()+
  labs( x="EpiCHAOS (original)", y="epiCHAOS (+50% sparsity)") +
  theme_classic()
dev.off()
