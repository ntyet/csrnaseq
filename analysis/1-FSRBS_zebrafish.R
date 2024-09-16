rm(list = ls())
pm1 <- proc.time()
# devtools::install_github("ntyet/csrnaseq")
library(tidyverse)
library(csrnaseq)
Bvec <- 100
mvec <- 4
alpha0vec <- 0.05
Bmalpha0 <- expand.grid(alpha0 = alpha0vec,m = mvec,  B = Bvec)
covsets <- readRDS("./extra-rna-seq-data/covsets2.rds")
counts <- readRDS("./extra-rna-seq-data/counts2.rds")
FixCov <- covsets["Treatment"]
VarCov <- covsets[c("Batch", "RIN")]

lambdamax <- 10
ncores <- 20 # ncores <- 4
print.progress <- F
saveall <- TRUE
dir.create(path = "revise/RealDataOutBSzebrafish", recursive = TRUE,  showWarnings = F)

### For HPC -----------------------------
pm1 <- proc.time()
env = as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID")) # env <- 1
print(env)

B <- Bmalpha0$B[env]; alpha0 <- Bmalpha0$alpha0[env]; m <- Bmalpha0$m[env]
res <- plyr::llply(c("WN", "OWN", "RX", "ORX"),
                     function(option)FSRAnalysisBS(counts,
                                                   FixCov, VarCov,
                                                   option, B, m,
                                                   lambdamax, alpha0, ncores,
                                                   print.progress, saveall))
names(res) <- c("WN", "OWN", "RX", "ORX")
saveRDS(res, file = paste0("revise/RealDataOutBSzebrafish/RealDataOutBSzebrafish_B_", B, "_m_", m, "_alpha0_", alpha0, ".rds"))
proc.time() - pm1

