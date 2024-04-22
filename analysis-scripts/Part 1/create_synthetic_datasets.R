
#---
#--- create perturbed datasets using the monocytes from cluster 0 from the Signac tutorial dataset
#---

counts <- readRDS("/omics/groups/OE0219/internal/KatherineK/data/scATAC/SignacTutorial_PBMCs/peak_counts_monocytes.Rds")

dim(counts)

# first convert to binary
counts[counts>1] <- 1

# filter for loci with high coverage across cells
select.loci <- which(rowSums(counts)>500)
counts <- counts[select.loci,]

# subset for fewer cells since computation is slow
counts <- counts[,which(colSums(counts)>2500)]
dim(counts)

# list to hold datasets for each perturbation
datasets <- list()
datasets$baseline <- counts

# we will perturb the cells, first try a 50% perturbation of 1's
base <- counts 

# remove 10% of 1's
base <- counts
len <- length(which(base>0))
select.sites <- len %>% sample(0.1*length(which(base>0))) # select at random 50% of non-zero elements to perturb
base[which(base>0)][select.sites] <- 0 # remove counts at randomly selected sites
datasets$rem10 <- base

# remove 20% of 1's
base <- counts
len <- length(which(base>0))
select.sites <- len %>% sample(0.2*length(which(base>0))) # select at random 50% of non-zero elements to perturb
base[which(base>0)][select.sites] <- 0 # remove counts at randomly selected sites
datasets$rem20 <- base

# remove 30% of 1's
base <- counts
len <- length(which(base>0))
select.sites <- len %>% sample(0.3*length(which(base>0))) # select at random 50% of non-zero elements to perturb
base[which(base>0)][select.sites] <- 0 # remove counts at randomly selected sites
datasets$rem30 <- base

# remove 40% of 1's
base <- counts
len <- length(which(base>0))
select.sites <- len %>% sample(0.4*length(which(base>0))) # select at random 50% of non-zero elements to perturb
base[which(base>0)][select.sites] <- 0 # remove counts at randomly selected sites
datasets$rem40 <- base

# remove 50% of 1's
len <- length(which(base>0))
select.sites <- len %>% sample(0.5*length(which(base>0))) # select at random 50% of non-zero elements to perturb
base[which(base>0)][select.sites] <- 0 # remove counts at randomly selected sites
datasets$rem50 <- base # replace the original portion of the data with perturbed version

lapply(datasets, sum)

# add 10% of 1's
base <- counts
len <- length(which(base==0))
select.sites <- len %>% sample(0.1*length(which(base>0))) # select at random 50% of non-zero elements to perturb
base[which(base==0)][select.sites] <- 1 # remove counts at randomly selected sites
datasets$add10 <- base

# add 20% of 1's
base <- counts
len <- length(which(base==0))
select.sites <- len %>% sample(0.3*length(which(base>0))) # select at random 50% of non-zero elements to perturb
base[which(base==0)][select.sites] <- 1 # remove counts at randomly selected sites
datasets$add30 <- base

# add 30% of 1's
base <- counts
len <- length(which(base==0))
select.sites <- len %>% sample(0.4*length(which(base>0))) # select at random 50% of non-zero elements to perturb
base[which(base==0)][select.sites] <- 1 # remove counts at randomly selected sites
datasets$add40 <- base

# add 40% of 1's
base <- counts
len <- length(which(base==0))
select.sites <- len %>% sample(0.2*length(which(base>0))) # select at random 50% of non-zero elements to perturb
base[which(base==0)][select.sites] <- 1 # remove counts at randomly selected sites
datasets$add20 <- base

# add 50% of 1's
base <- counts
len <- length(which(base==0))
select.sites <- len %>% sample(0.5*length(which(base>0))) # select at random 50% of non-zero elements to perturb
base[which(base==0)][select.sites] <- 1 # remove counts at randomly selected sites
datasets$add50 <- base

#--- datasets with increasing homogeneity

# here, cells are perturbed non-randomly, i.e. by selecting loci (at random) and removing 1's at those same loci from all cells
ctrl <- counts 
ctrl[sample(nrow(ctrl), nrow(counts)/2),] <- 0
datasets$nr.ctrl <- ctrl

lapply(datasets, sum)


#--- select 50% of rows/loci, then remove counts in increments of 10-50 % at the selected loci
ctrl <- counts 
sample.rows <- sample(nrow(ctrl), nrow(counts)/2)

# 10%
n.perturb <- ctrl[sample.rows,] %>% sum() %>% multiply_by(0.1) %>% round()
ctrl[sample.rows,][ctrl[sample.rows,]==1][sample(n.perturb)] <- 0
datasets$rem10.hom <- ctrl

# 20%
ctrl <- counts
n.perturb <- ctrl[sample.rows,] %>% sum() %>% multiply_by(0.2) %>% round()
ctrl[sample.rows,][ctrl[sample.rows,]==1][sample(n.perturb)] <- 0
datasets$rem20.hom <- ctrl

# 30%
ctrl <- counts
n.perturb <- ctrl[sample.rows,] %>% sum() %>% multiply_by(0.3) %>% round()
ctrl[sample.rows,][ctrl[sample.rows,]==1][sample(n.perturb)] <- 0
datasets$rem30.hom <- ctrl

# 40%
ctrl <- counts
n.perturb <- ctrl[sample.rows,] %>% sum() %>% multiply_by(0.4) %>% round()
ctrl[sample.rows,][ctrl[sample.rows,]==1][sample(n.perturb)] <- 0
datasets$rem40.hom <- ctrl

# 50%
ctrl <- counts
n.perturb <- ctrl[sample.rows,] %>% sum() %>% multiply_by(0.5) %>% round()
ctrl[sample.rows,][ctrl[sample.rows,]==1][sample(n.perturb)] <- 0
datasets$rem50.hom <- ctrl

lapply(datasets, sum)



#--- now cases which 50% are added non-randomly, for comparison to that in which I randomly added 50%

# 10 %
n <- sum(datasets$add10)-(sum(datasets$baseline)) # n is the number of 1s I want to add
base <- datasets$baseline
n.row <- sample(nrow(base), nrow(base)/2) # select, arbitrarily, 25% of rows (I don't want to create rows which have all 1's because those never occur in reality, so I will randomly add 1's to non-random selection of rows)
select.sites <- length(which(base[n.row,]==0)) %>% sample(n)
base[n.row,][base[n.row,]==0][select.sites] <- 1
datasets$add10.hom <- base

# 20%
n <- sum(datasets$add20)-(sum(datasets$baseline)) # n is the number of 1s I want to add
base <- datasets$baseline
n.row <- sample(nrow(base), nrow(base)/2) # select, arbitrarily, 25% of rows (I don't want to create rows which have all 1's because those never occur in reality, so I will randomly add 1's to non-random selection of rows)
select.sites <- length(which(base[n.row,]==0)) %>% sample(n)
base[n.row,][base[n.row,]==0][select.sites] <- 1
datasets$add20.hom <- base

# 30%
n <- sum(datasets$add30)-(sum(datasets$baseline)) # n is the number of 1s I want to add
base <- datasets$baseline
n.row <- sample(nrow(base), nrow(base)/2) # select, arbitrarily, 25% of rows (I don't want to create rows which have all 1's because those never occur in reality, so I will randomly add 1's to non-random selection of rows)
select.sites <- length(which(base[n.row,]==0)) %>% sample(n)
base[n.row,][base[n.row,]==0][select.sites] <- 1
datasets$add30.hom <- base

# 40%
n <- sum(datasets$add40)-(sum(datasets$baseline)) # n is the number of 1s I want to add
base <- datasets$baseline
n.row <- sample(nrow(base), nrow(base)/2) # select, arbitrarily, 25% of rows (I don't want to create rows which have all 1's because those never occur in reality, so I will randomly add 1's to non-random selection of rows)
select.sites <- length(which(base[n.row,]==0)) %>% sample(n)
base[n.row,][base[n.row,]==0][select.sites] <- 1
datasets$add40.hom <- base

# 50%
n <- sum(datasets$add50)-(sum(datasets$baseline)) # n is the number of 1s I want to add
base <- datasets$baseline
n.row <- sample(nrow(base), nrow(base)/2) # select, arbitrarily, 25% of rows (I don't want to create rows which have all 1's because those never occur in reality, so I will randomly add 1's to non-random selection of rows)
select.sites <- length(which(base[n.row,]==0)) %>% sample(n)
base[n.row,][base[n.row,]==0][select.sites] <- 1
datasets$add50.hom <- base

lapply(datasets, sum)


saveRDS(datasets, "/omics/groups/OE0219/internal/KatherineK/ATACseq/variability_methods/perturbed_datasets_signac_monocytes.Rds")
