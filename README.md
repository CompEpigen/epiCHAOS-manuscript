
# EpiCHAOS

A metric for calculating cell-to-cell heterogeneity in single-cell epigenomics data


#### Calculation of epiCHAOS scores:  

Required: (i) a single cell epigenomics dataset in binarised matrix form, e.g. a peaks-by-cells or tiles-by-cells matrix in the case of scATAC data, (ii) cell annotation matching column names to clusters, cell types or other groups on which heterogenetiy scores should be computed

1. Pairwise chance-corrected Jaccard indices are calculated between all cells in each group/cluster
2. The mean of all pairwise indices is computed per cluster as a raw heterogeneity score. Scores are negated so that a higher score refects higher heterogeneity.
3. Raw heterogeneity scores are adjusted for differences in total counts between groups by fitting a linear regression model and taking the residuals of the model as an adjusted score
5. Scores are normalised to a range of 0-1 and subtracted from one so that a higher score indicates higher cell-to-cell heterogeneity
