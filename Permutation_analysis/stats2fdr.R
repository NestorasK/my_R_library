# Permutation based FDR #
# Citation: A nonparametric approach for identifying differential expression in RNA-Seq data, Tibshirani

library(data.table)

# Information of function #
# Calculate FDR based on statistics
# Function 
# Small modifications in relation to publication
# - No 0s are allowed in Rs 
# - The one 0 that is produced we turn it to 1
FDR_calculation <- function(statistics){
    # statistics: matrix
    # columns: statistics, first column initial statistics, 
    #          other columns statistics created during permutations
    # rows: genes  
    
    statistics <- abs(statistics)
    num_perm <- ncol(statistics) - 1
    
    # Calculating Rs
    init_stats = statistics[,1]
    
    # Calculating testing Rs
    init_stats_rank <- (frank(-init_stats, ties.method = 'min') - 1)
    Rs <- (frank(-init_stats, ties.method = 'min') - 1)
    Rs[Rs == 0] <- 1
    
    # Calculate Vs
    perm_stats <- as.numeric(statistics[,2:(num_perm + 1)])
    
    # Calculate testing Vs
    all_stat_rank <- (frank(-c(init_stats, perm_stats), ties.method = 'min') - 1)
    Vs <- (all_stat_rank[1:length(init_stats_rank)] - init_stats_rank) / num_perm
    
    # Calculate pi0
    q <- median(x = perm_stats, na.rm = TRUE)
    pi0 <- (2*sum(init_stats <= q))/ length(init_stats)
    
    # FDRs
    FDR <- pi0*Vs/Rs
    
    # add names
    names(FDR) <- rownames(statistics)
    
    # return object
    return(FDR)
}
