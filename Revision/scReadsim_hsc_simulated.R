
library(epiCHAOS)
library(dplyr)
library(stringr)

#--- read in HSC scATAC subset data and write to txt file for scReadSim
# setwd("/Users/kayke/scReadSim")
# hsc <- readRDS("subset_hscs_for_simulation.Rds")
# hsc <- hsc[sample(nrow(hsc), 10000),] # subset for ease of processing
# dim(hsc)
# hsc %>% as.matrix() %>% write.table(file = "subset_hscs_for_simulation.txt", sep = "\t")

#--- path to files containing simulated scATAC-seq data
files <- list.files("/omics/groups/OE0219/internal/KatherineK/data/scATAC/scReadSim_HSCs")

set.seed(11)
datasets <- list()
for (i in 1:length(files)) {
  temp <- read.table(file.path("/omics/groups/OE0219/internal/KatherineK/data/scATAC/scReadSim_HSCs", files[[i]]), row.names = 1)
  temp[temp>1] <- 1
  datasets[[paste0("group-", i)]] <- temp
  #--- take subsamples of 100 cells for each sequencing depth
  # for (rep in 1:5) {
  #   datasets[[paste0(files[[i]], "-", rep)]] <- temp[,sample(ncol(temp), 100)]
  # }
}

lapply(datasets, dim)
counts <- lapply(datasets, sum) %>% unlist()
length(datasets)

#--- compute epiCHAOS scores
het <- compute_eITH(datasets)
het$count <- counts

#--- check how adjustment for sparsity affected the scores
plot(het$het.adj, het$het.raw)

#--- adjust labels for plotting
het$n <- het$state %>% str_remove("subset_hscs_for_simulation.scDesign2Simulated_") %>% str_remove("\\.txt") %>%
  str_split("-") %>% lapply("[", 1) %>% unlist() %>% as.numeric()

het$state <- as.factor(het$n)

#--- check relationship with simulated sequencing depth before and after plotting
plot(het$n, het$het.adj)
plot(het$n, het$het.raw)

#--- boxplot
p1 <- ggplot(het, aes(x=state, y=het.raw)) +
  geom_boxplot(fill="steelblue") +
  geom_beeswarm()+
  labs(x="", y="epiCHAOS")+
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90))

p2 <- ggplot(het, aes(x=state, y=het.adj)) +
  geom_boxplot(fill="steelblue") +
  geom_beeswarm()+
  labs(x="", y="epiCHAOS")+
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90))

svg("/omics/groups/OE0219/internal/KatherineK/ATACseq/epiCHAOS-Figures/Supplementary/Revision/boxplots_epiCHAOS_scores_simulate_seq_depth.svg", 6, 2.5)
ggpubr::ggarrange(p1,p2, ncol=2)
dev.off()
