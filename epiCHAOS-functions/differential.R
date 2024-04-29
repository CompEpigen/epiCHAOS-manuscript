
library(dplyr)
library(magrittr)
library(ggplot2)

#--- function to compute differential epigenetic heterogeneity between two groups
#-   group1 and group2 are peaks-by-cells matrices for the two groups of interest
#-   region.type is the name of the selected region type e.g. the TFBS
#-   n.iter is the number of permutations, defaults to 1000
compute.diff.eICH <- function(group1, group2, region.type, niter=1000) {
  
  #--- bind the two peaks-by-cells matrices
  merged <- cbind(group1, group2)
  
  #--- list to hold the individual peaks-by-cells matrices
  datasets <- list(group1=group1, group2=group2)
  
  #--- randomly sample 100 cells from the pool of cell from the two groups, repeat in 100 (or selected "n.iter") iterations, then add the permuted matrices to the list of peaks-by-cells matrices
  for (i in 1:niter) {
    sample.cells <- c(colnames(merged)) %>% sample(100)
    datasets[[paste0("mix", i)]] <- merged[,sample.cells]
  }
  
  #--- compute epiCHAOS scores on the permuted and true matrices
  het <- compute.eITH(datasets)
  
  #--- epiCHAOS scores for permuted groups, excluding those of the true groupings
  het.dist <- het[het$state %notin% c("group1", "group2"),]
  
  #--- find the differences in epiCHAOS scores between pairs of randomly permuted matrices - that will indicate the expected difference in heterogeneity between two randomly selected groups
  dif.dist <- c()
  for (i in 1:nrow(het.dist)/2) {
    dif <- het.dist$mean.het[i]-het.dist$mean.het[i+1]
    dif.dist <- c(dif.dist, dif)
  }
  
  #--- density plot of the differences in epiCHAOS scores from pairs of randomly permuted matrices, with red line indicating the difference in epiCHAOS scores between the true groupings
  df <- data.frame(dif.dist=dif.dist, condition="mix")
  gg <- ggplot(df, aes(x=dif.dist)) + 
    geom_density() +
    lims(x=c(-1,1)) +
    labs(x="", subtitle = region.type) +
    geom_vline(aes(xintercept=(het$mean.het[het$state=="group1"]-het$mean.het[het$state=="group2"])), color="red3", linetype="dashed") +
    theme_bw()
  
  #--- return the results from the permuted mixtures, from the test comparison, and the density plot
  return(list(plot=gg, result.mix=df, result.test=(het$mean.het[het$state=="group1"]-het$mean.het[het$state=="group2"])))
  
}