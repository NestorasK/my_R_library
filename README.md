# my_R_library
Useful generic R functions

### File: calculate_correlations.R
<code>calculate_cors_fast()</code> function calculating correlations, t-statistic, pvalues and standard errors. 

- Outputs a flat data.table, 
- Faster in calculating pvalues than what I have found in base R and some available packages such as Hmisc, psych. 
- It handles one or two different tables

### File: Permutation_analysis/stats2pvalues.R
Functions converting a vector of statistics to p.values. Three functions have been implemented with small differences. My suggestion is to use <code>my.t2p()</code> for reasons analyzed here: 
https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0165545
Citations: omicsNPC: Applying the Non-Parametric Combination Methodology to the Integrative Analysis of Heterogeneous Omics Data.PLoS One 11, e0165545 (2016)

### File: Permutation_analysis/stats2fdr.R
<code>FDR_calculation()</code>: function converting a matrix of statistic to qvalues. 
Citation: A nonparametric approach for identifying differential expression in RNA-Seq data, Tibshirani


