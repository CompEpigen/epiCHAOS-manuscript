

#--- use data for sorted monocytes as an example to create perturbed datasets
counts <- readRDS("SignacTutorial_PBMCs/peak_counts_monocytes.Rds")

dim(counts)

#--- first binarise the scATAC matrix
counts[counts>1] <- 1

#--- select loci with relatively high coverage across cells
select.loci <- which(rowSums(counts)>500)
counts <- counts[select.loci,]

#--- subset for fewer cells for ease of computation
counts <- counts[,which(colSums(counts)>2500)]
dim(counts)

#--- "datasets" list will hold datasets for the baseline and for each perturbation
datasets <- list()
datasets$baseline <- counts

#--- begin perturbations from the baseline counts
base <- counts 

#--- remove 10% of 1's
base <- counts
len <- length(which(base>0))
select.sites <- len %>% sample(0.1*length(which(base>0))) # select at random 10% of non-zero elements to perturb
base[which(base>0)][select.sites] <- 0 # remove counts at randomly selected sites
datasets$rem10 <- base

#--- remove 20% of 1's
base <- counts
len <- length(which(base>0))
select.sites <- len %>% sample(0.2*length(which(base>0))) # select at random 20% of non-zero elements to perturb
base[which(base>0)][select.sites] <- 0 # remove counts at randomly selected sites
datasets$rem20 <- base

#--- remove 30% of 1's
base <- counts
len <- length(which(base>0))
select.sites <- len %>% sample(0.3*length(which(base>0))) # select at random 30% of non-zero elements to perturb
base[which(base>0)][select.sites] <- 0 # remove counts at randomly selected sites
datasets$rem30 <- base

#--- remove 40% of 1's
base <- counts
len <- length(which(base>0))
select.sites <- len %>% sample(0.4*length(which(base>0))) # select at random 40% of non-zero elements to perturb
base[which(base>0)][select.sites] <- 0 # remove counts at randomly selected sites
datasets$rem40 <- base

#--- remove 50% of 1's
len <- length(which(base>0))
select.sites <- len %>% sample(0.5*length(which(base>0))) # select at random 50% of non-zero elements to perturb
base[which(base>0)][select.sites] <- 0 # remove counts at randomly selected sites
datasets$rem50 <- base # replace the original portion of the data with perturbed version

lapply(datasets, sum)

#--- add 10% of 1's
base <- counts
len <- length(which(base==0))
select.sites <- len %>% sample(0.1*length(which(base>0))) # select at random 10% of non-zero elements to perturb
base[which(base==0)][select.sites] <- 1 # add counts at randomly selected sites
datasets$add10 <- base

#--- add 20% of 1's
base <- counts
len <- length(which(base==0))
select.sites <- len %>% sample(0.2*length(which(base>0))) # select at random 20% of non-zero elements to perturb
base[which(base==0)][select.sites] <- 1 # add counts at randomly selected sites
datasets$add20 <- base

#--- add 30% of 1's
base <- counts
len <- length(which(base==0))
select.sites <- len %>% sample(0.3*length(which(base>0))) # select at random 30% of non-zero elements to perturb
base[which(base==0)][select.sites] <- 1 # add counts at randomly selected sites
datasets$add30 <- base

#--- add 40% of 1's
base <- counts
len <- length(which(base==0))
select.sites <- len %>% sample(0.4*length(which(base>0))) # select at random 40% of non-zero elements to perturb
base[which(base==0)][select.sites] <- 1 # add counts at randomly selected sites
datasets$add40 <- base

#--- add 50% of 1's
base <- counts
len <- length(which(base==0))
select.sites <- len %>% sample(0.5*length(which(base>0))) # select at random 50% of non-zero elements to perturb
base[which(base==0)][select.sites] <- 1 # add counts at randomly selected sites
datasets$add50 <- base

saveRDS(datasets, "synthetic-data-tests/perturbed_data_monocytes.Rds")
