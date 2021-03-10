# Simulation analysis
rm(list = ls())
source(file = "~/scratch/my_R_library/Permutation_analysis/stats2fdr.R")
source(file = "~/scratch/my_R_library/Permutation_analysis/stats2pvalues.R")
set.seed(124)
numP <- 1000
numPerm <- 2001
data_NULL <- matrix(data = rnorm(n = numP * numPerm, mean = 10), nrow = numP, ncol = numPerm)
data_Sig <- data_NULL
data_Sig[1:(nrow(data_Sig) * 0.1), 1] <- 20


# Pvalues NULL
pvalues_t2p_NULL <- apply(X = data_NULL, MARGIN = 1, FUN = t2p)
# hist(pvalues_t2p_NULL[1,])

pvalues_Vin_NULL <- apply(X = data_NULL, MARGIN = 1, FUN = computePvaluesVect)
# hist(pvalues_Vin_NULL[1,])

pvalues_myt2p_NULL <- apply(X = data_NULL, MARGIN = 1, FUN = my.t2p)
# hist(pvalues_myt2p_NULL[1,])

pvalues_NULL <- cbind(pvalues_t2p_NULL[1,], pvalues_Vin_NULL[1,], pvalues_myt2p_NULL[1,])
colnames(pvalues_NULL) <- c("t2p", "Vin", "myt2p")

cor(pvalues_NULL)
boxplot(pvalues_NULL, main = "pvalues - NULL")
apply(X = pvalues_NULL, MARGIN = 2, FUN = summary)


# Pvalues sig
pvalues_t2p_Sig <- apply(X = data_Sig, MARGIN = 1, FUN = t2p)
# hist(pvalues_t2p_Sig[1,])

pvalues_Vin_Sig <- apply(X = data_Sig, MARGIN = 1, FUN = computePvaluesVect)
# hist(pvalues_Vin_Sig[1,])

pvalues_myt2p_Sig <- apply(X = data_Sig, MARGIN = 1, FUN = my.t2p)
# hist(pvalues_myt2p_Sig[1,])

pvalues_Sig <- cbind(pvalues_t2p_Sig[1,], pvalues_Vin_Sig[1,], pvalues_myt2p_Sig[1,])
colnames(pvalues_Sig) <- c("t2p", "Vin", "myt2p")

cor(pvalues_Sig)
boxplot(pvalues_Sig, main = "pvalues - Sig")
apply(X = pvalues_Sig, MARGIN = 2, FUN = summary)


# FDR - p.adjust
# - NULL
PermFDR <- FDR_calculation(statistics = data_NULL)
qvalues_NULL <- apply(X = pvalues_NULL, MARGIN = 2, FUN = p.adjust, method = "fdr")
qvalues_NULL <- cbind(qvalues_NULL, PermFDR)

boxplot(qvalues_NULL, main = "qvalues - NULL")
cor(qvalues_NULL)
apply(X = qvalues_NULL, MARGIN = 2, FUN = summary)

apply(X = qvalues_NULL <= 0.005, MARGIN = 2, FUN = sum)


# - Sig
PermFDR <- FDR_calculation(statistics = data_Sig)
qvalues_Sig <- apply(X = pvalues_Sig, MARGIN = 2, FUN = p.adjust, method = "fdr")
qvalues_Sig <- cbind(qvalues_Sig, PermFDR)

boxplot(qvalues_Sig, main = "qvalues - Sig")
cor(qvalues_Sig)
apply(X = qvalues_Sig, MARGIN = 2, FUN = summary)


apply(X = qvalues_Sig <= 0.005, MARGIN = 2, FUN = sum)
# apply(X = qvalues_Sig <= 0.05, MARGIN = 2, FUN = which)




