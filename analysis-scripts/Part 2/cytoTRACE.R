

#--- run cytoTRACE on scRNAseq data from granja et al. hematopoiesis

#--- rna
rna <- readRDS("/omics/groups/OE0219/internal/KatherineK/data/scATAC/Granja2019/scRNA-Healthy-Hematopoiesis-191120.rds")

#--- subset bone marrow derived cells
rna <- rna[,grepl("BMMC", rna$Group)]
table(rna$Group)

results <- CytoTRACE(as.matrix(rna@assays@.xData$data$counts))