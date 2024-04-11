
library(jaccard)


#--- trying chance corrected Jaccard with regressing out count --> adjusted method so that 
compute.eITH <- function(x) {  # x should be a list of scATAC datasets with equal dimensions (n cells/cols and n loci/rows)
  
  # dists is a list which will hold the raw obs-exp deviations before scaling and negation
  dists <-  list()
  ids <- list()
  for (d in names(x)) { 
    
    print(paste0("calculating eITH for ", d))
    temp <- x[[d]] %>% as.matrix() %>% as.data.frame()
    
    if (sum(temp)==0) { next }
    temp <- temp[,colSums(temp)>0]
    x[[d]] <- temp
    
    jac <- corrr::colpair_map(temp, jaccard, center=T)
    jac <- apply(jac, 1, as.numeric)
    
    dists[[d]] <- colMeans(jac, na.rm = T)  #*colMeans(sums, na.rm=T)
    ids[[d]] <- colnames(temp)
    
    
  }
  
  #--- create dataframe to hold heterogeneity scores
  het <- data.frame(het=unlist(dists), id=unlist(ids))
  state <- c()
  for (cond in names(dists)) {
    state <- c(state, rep(cond, length(dists[[cond]])))
  }
  
  het$state <- state
  
  for (state in unique(het$state)) {
    het$mean.het[het$state==state] <- mean(het$het[het$state==state])
  }
  
  het <- unique(het[,c("mean.het", "state")])
  
  #--- using residuals to correct score for total count
  het.d <- het$mean.het
  
  sum.whole <- c()
  for (i in names(x)) {
    sum.whole <- c(sum.whole, mean(colMeans(x[[i]])))  # changed colSums to colMeans 31.01.2024
  }
  
  fit <- lm(het.d~sum.whole)
  het$mean.het <- residuals(fit)
  
  #--- convert values to a range of 0,1
  het$mean.het <- (het$mean.het - min(het$mean.het)) / (max(het$mean.het) - min(het$mean.het))  # Scale to a range of 0-1
  het$mean.het <- 1-het$mean.het
  
  #--- return a dataframe with colnames "het" and "state" representing the heterogeneity score and condition respectively
  return(het[,c("mean.het", "state")])
  
}




#--- epiCHAOS calculation without adjustment for total counts
compute.eITH.raw <- function(x) {  # x should be a list of scATAC datasets with equal dimensions (n cells/cols and n loci/rows)
  
  # dists is a list which will hold the raw obs-exp deviations before scaling and negation
  dists <-  list()
  ids <- list()
  for (d in names(x)) { 
    
    print(paste0("calculating eITH for ", d))
    temp <- x[[d]] %>% as.matrix() %>% as.data.frame()
    
    if (sum(temp)==0) { next }
    temp <- temp[,colSums(temp)>0]
    x[[d]] <- temp
    
    jac <- corrr::colpair_map(temp, jaccard, center=T)
    jac <- apply(jac, 1, as.numeric)
    
    dists[[d]] <- colMeans(jac, na.rm = T)  #*colMeans(sums, na.rm=T)
    ids[[d]] <- colnames(temp)
    
    
  }
  
  #--- create dataframe to hold heterogeneity scores
  het <- data.frame(het=unlist(dists), id=unlist(ids))
  state <- c()
  for (cond in names(dists)) {
    state <- c(state, rep(cond, length(dists[[cond]])))
  }
  
  het$state <- state
  
  for (state in unique(het$state)) {
    het$mean.het[het$state==state] <- mean(het$het[het$state==state])
  }
  
  het <- unique(het[,c("mean.het", "state")])
  
  #--- convert values to a range of 0,1
  het$mean.het <- (het$mean.het - min(het$mean.het)) / (max(het$mean.het) - min(het$mean.het))  # Scale to a range of 0-1
  het$mean.het <- 1-het$mean.het
  
  #--- return a dataframe with colnames "het" and "state" representing the heterogeneity score and condition respectively
  return(het[,c("mean.het", "state")])
  
}

# counts: a counts matrix representing e.g. a peaks by cells matrix in the case of scATAC data. The counts matrix will be binarised if not already.
# meta: metadata including the column on which the data should be grouped, e.g. cluster or celltype.
# colname: the name of the column in "meta" on which to group the data e.g. "cluster" or "celltype". If not specified this defaults the the first column in "meta".
# n: the number of cells to subset for each group/cluster. This defults to 100.
# index: the rows in counts on which to subset the counts matrix. If not provided the whole counts matrix will be used by default. Otherwise "index" may be specified as either a (i) vector of numerical indices, (ii) a vector or names corresponding to the rownames of interest in "counts"
# plot: if true, a barplot will be returned as well as the epiCHAOS scores. Defaults to FALSE.

create.group.matrices <- function(counts, meta, colname, n=100, index=NULL) {

  meta$group <- meta[,colname]
  counts <- na.omit(counts)
  
  if (index) { counts <- counts[index, ] }
  
  # a list to hold counts matrices for each group/cluster
  matrices <- list()
  
  for (group in unique(meta$group)) {
    #print(group)
    ids <- meta[meta$group==group, ] %>% rownames()
    matrices[[group]] <- counts[,ids] %>% as.matrix()
    
    # if more cells than selected n, downsample for n cells for that group
    if (length(ids)>n) { matrices[[group]] <- counts[,sample(ids, n)] %>% as.matrix() }
    
    # binarise counts if not already
    matrices[[group]][matrices[[group]]>1] <- 1
  }
  
  return(matrices)
}


#--- epiCHAOS with count correction per chromosome
compute.eITH.cancer <- function(x) {  # x should be a list of scATAC datasets with equal dimensions (n cells/cols and n loci/rows)
  
  chromosomes <- rownames(x[[1]]) %>% str_split("-|_|:") %>% lapply("[", 1) %>% unlist() %>% unique()
  
  # dists will hold jaccard distances, het.chr will hold per-chromoosme heterogeneity scores, coverage will hold per chromosome coverage
  het.chr <- dists <- coverage <- list()
  for (chr in paste0(chromosomes, ":")) {
    
    
    for (d in names(x)) { 
      
      print(paste0("calculating eITH for ", d, " ", chr))
      temp <- x[[d]] %>% as.matrix() %>% as.data.frame()
      
      # subset atac matrix for selected chromosome
      temp <- temp[grepl(chr, rownames(temp)),]
      
      #if (sum(temp)==0) { next }
      #temp <- temp[,colSums(temp)>0]
      
      # compute pairwise Jaccard distances
      jac <- corrr::colpair_map(temp, jaccard, center=T)
      jac <- apply(jac, 1, as.numeric)
      
      # get mean of pairwise jaccard indices
      dists[[d]] <- colMeans(jac, na.rm = T) %>% mean()
      coverage[[d]] <- colMeans(temp) %>% mean()
      
    }
    
    # regress out counts
    fit <- lm(unlist(dists)~unlist(coverage))
    
    
    #--- create dataframe to hold heterogeneity scores
    het.chr[[chr]] <- (fit$residuals - min(fit$residuals)) / (max(fit$residuals) - min(fit$residuals)) 
    
  }
  
  require(plyr)
  het <-  rbind.fill(lapply(het.chr,function(y){as.data.frame(t(y),stringsAsFactors=FALSE)})) %>% colMeans(na.rm = T)
  het <- data.frame(het=het, state=names(het))
  
  #--- convert values to a range of 0,1
  het$het <- (het$het - min(het$het)) / (max(het$het) - min(het$het))  # Scale to a range of 0-1
  het$het <- 1-het$het
  
  
  return(het)
  
}


# compute epiCHAOS if given a counts matrix and metadata
epiCHAOS <- function(counts, meta, colname=colnames(meta)[1], n=100, index=NULL, plot=F, cancer=F) {
  
  print("creating group matrices")
  matrices <- create.group.matrices(counts, meta, colname, n, index)
  print("computing epiCHAOS scores")
  het <- ifelse (cancer==T, compute.eITH.cancermatrices, compute.eITH(matrices))
  
  if (plot) {
    p <- ggplot2::ggplot(het, aes(x = reorder(state, mean.het), y = mean.het+0.01)) +
      geom_bar(stat="identity", position = "dodge", alpha=0.8, width = 0.6)+
      labs(x="", y="epiCHAOS")+
      theme_classic() +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
    return(list(het=het, barplot=p))
  } else { 
    return(het) 
    }
  
}



