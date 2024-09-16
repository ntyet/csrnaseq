# This is basically analyzing real data with all FixCov
# without covariate selection, where the number of FixCov will be
# 0, 1; corresponding to the model size of each set of
# covariates.

rm(list = ls())
pm1 <- proc.time()
library(csrnaseq)

covsets <- readRDS("./extra-rna-seq-data/covsets2.rds")
counts <- readRDS("./extra-rna-seq-data/counts2.rds")
FixCov <- covsets["Treatment"]
VarCov <- covsets[c("Batch", "RIN")]

Config <- list(ModelSize_0 = covsets[c("Treatment")],
               ModelSize_1 = covsets[c("Treatment", "Batch")])

RealDataOut <- NULL
dir.create(path = "revise/RealDataOutBSzebrafish", showWarnings = F, recursive = T)
for (i in names(Config)){
  AllCov <- Config[[i]]
  RealDataOut$VoomPvOutER <- csrnaseq:::VoomPv(counts, AllCov)
  RealDataOut$FixCov <- FixCov
  RealDataOut$VarCov <- VarCov
  RealDataOut$BestCovOut$BestER <- rep(2, length(AllCov)-1) # just initial value, not important
  names(RealDataOut$BestCovOut$BestER) <- names(AllCov[-1]) # relevance
  saveRDS(RealDataOut, file = paste0("revise/RealDataOutBSzebrafish/", i, ".rds")) # use this if run on Wahba HPC
  }


