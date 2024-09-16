# This is basically analyzing real data with all FixCov
# without covariate selection, where the number of FixCov will be
# 0, 1, 2, 6, 7, 8; corresponding to the model size of each set of
# covariates.

rm(list = ls())
pm1 <- proc.time()
library(csrnaseq)
FixCov0 <- covset["Line"]
VarCov0 <- covset[c("Diet", "RFI", "Concb", "RINb", "Conca", "RINa", "neut", "lymp",
                    "mono", "eosi", "baso", "Block", "Order")]

Config <- list(ModelSize_0 = covset["Line"],
               ModelSize_1 = covset[c("Line", "mono")],
               ModelSize_2 = covset[c("Line", "mono", "Concb")],
               ModelSize_3 = covset[c("Line", "mono", "Concb", "neut")],
               ModelSize_4 = covset[c("Line", "mono", "Concb", "neut", "Block")],
               ModelSize_5 = covset[c("Line", "mono", "Concb", "neut", "Block", "RINa")],
               ModelSize_6 = covset[c("Line", "mono", "Concb", "neut", "Block", "RINa", "baso")],
               ModelSize_7 = covset[c("Line", "mono", "Concb", "neut", "Block", "RINa", "baso", "lymp")],
               ModelSize_8 = covset[c("Line", "mono", "Concb", "neut", "Block", "RINa", "baso", "lymp", "RFI")])

RealDataOut <- NULL
dir.create(path = "RealDataOutBS", showWarnings = F, recursive = T)
for (i in names(Config)){
  AllCov <- Config[[i]]
  RealDataOut$VoomPvOutER <- csrnaseq:::VoomPv(counts, AllCov)
  RealDataOut$FixCov <- FixCov0
  RealDataOut$VarCov <- VarCov0
  RealDataOut$BestCovOut$BestER <- rep(2, length(AllCov)-1) # nuisance
  names(RealDataOut$BestCovOut$BestER) <- names(AllCov[-1]) # relevance
  saveRDS(RealDataOut, file = paste0("RealDataOutBS/", i, ".rds")) # use this if run on Wahba HPC
  }


