# Short description #
# Correlations calculations using R

# Function: calculate_cors_fast()
# - Outputs a flat data.table
# - Calculates correlations, t-statistic, pvalues, standard error
# - It is faster in calculating pvalues than what I have found in base R and some available packages such as Hmisc, psych.
# - It can handle one or two different tables

library(data.table)
cor2pvalue = function(r, n) {
    
    # Code adopted from 
    # https://stackoverflow.com/questions/13112238/a-matrix-version-of-cor-test
    
    t <- (r*sqrt(n-2))/sqrt(1-r^2)
    p <- 2*(1 - pt(abs(t),(n-2)))
    se <- sqrt((1-r*r)/(n-2))
    out <- list(#r = r, n = n, 
        t = t,
        p = p,
        se = se
    )
    return(out)
}
flattenCorrMatrix <- function(cormat) {
    
    # Code adopted from 
    # https://rdrr.io/github/heuselm/mocode/src/R/flattenCorrMatrix.R
    
    ut <- upper.tri(cormat)
    data.table(
        row = rownames(cormat)[row(cormat)[ut]],
        column = rownames(cormat)[col(cormat)[ut]],
        cor  = (cormat)[ut]
    )
}
calculate_cors_fast <- function(table_i, table_j = NULL, ...) {
    
    # Input #
    # table_i: data.table, Column "Sample" is expected
    # table_j: data.table, Column "Sample" is expected
    # ... : arguments for cor() function of base R
    
    # Output # 
    # Flat data.table (not a square matrix) with the following columns
    # - row: columns from table_i
    # - column: columns from table_j
    # - tstat: t-statistic as calculated from cor2pvalue
    # - pvalue: respective p.value 
    # - se: standard error of the correlation
    
    if (is.null(table_j)) {
        cors <- cor(x = as.matrix(table_i[,-"Sample"]), ...)
        corMat_flat_i <- flattenCorrMatrix(cormat = cors)
        
    }else{
        setkey(x = table_i, Sample)
        setkey(x = table_j, Sample)
        cors <- cor(x = as.matrix(table_i[,-"Sample"]), y = as.matrix(table_j[, -"Sample"]), ...)
        corMat_flat_i <- melt.data.table(data = data.table(cors, keep.rownames = TRUE), id.vars = "rn", 
                                         variable.name = "column", value.name = "cor")
        colnames(corMat_flat_i)[1] <- "row"
    }
    
    corMat_flat_i[, c("tstat", "pvalue", "se") := cor2pvalue(r = cor, n = length(table_i[, Sample]))]
    
    # Output
    return(corMat_flat_i)
    
}

