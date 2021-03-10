library(data.table)
# Statistics to pvalues
# Citations: 
# Pesarin, F. & Salmaso, L. Permutation tests for complex data: theory, applications, and software. (Wiley, 2010).
t2p <- function(stats){    
    # Pesarin
    if(is.null(dim(stats))){stats<-array(stats,dim=c(length(stats),1))}
    oth<-seq(1:length(dim(stats)))[-1]    
    B<-dim(stats)[1]-1          		
    p<-dim(stats)[2]	
    if(length(dim(stats))==3){C<-dim(stats)[3]}
    rango<-function(x){
        r=1-rank(x[-1],ties.method="min")/B+1/B
        return(c(mean(x[-1]>=x[1]),r))
    }  
    P=apply(stats,oth,rango)
    return(P)
}

# omicsNPC: Applying the Non-Parametric Combination Methodology to the Integrative Analysis of Heterogeneous Omics Data.PLoS One 11, e0165545 (2016)
# Similar to pesarin but pvalues are never 0
my.t2p <- function(statsVec){    
    # Number of permutations
    B <- sum(!is.na(statsVec))
    
    # p.values
    p.vec <- (1 - rank(x = statsVec, na.last = "keep", ties.method = "min")/B) + 1/B
    return(p.vec)
}

# Code adopted from omicsNPC online package
# Pvalues are never 0 nor 1
computePvaluesVect <- function(statsVect){        
    
    #how many values are present in the vector of statistics?
    numStatistics <- sum(!is.na(statsVect))
    
    #the rank function with ties.method = 'min' indicates how many statistics are 
    #larger or equal than each value in the vector
    #Note: we assume that the greater the statistics, the more significant the finding
    #This is the case with, for example, chi-square statistics
    #This is not the case for on-tailed t-test where more negative statistics are actually more extreme
    #The functions computing the statistics should be programmed accordingly; 
    #for example, by changing the sign of negative statistics and putting to zero positive once, 
    #or by returning the -log(desiredOneTailedPvalue)
    numLargerValues <- numStatistics - rank(statsVect, na.last = "keep", ties.method = "min");
    
    #correction for avoiding zero p-values
    pvaluesVect <- (numLargerValues + 1) / (numStatistics + 1);
    
    #returning the p-values
    return(pvaluesVect)
    
}  
