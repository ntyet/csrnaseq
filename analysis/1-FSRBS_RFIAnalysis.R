rm(list = ls())
pm1 <- proc.time()
# devtools::install_github("ntyet/csrnaseq")
library("csrnaseq")
Bvec <- 100
mvec <- c(1, 3, 5, 7)
alpha0vec <- c(0.01, 0.05, 0.1, 0.2)
Bmalpha0 <- expand.grid(alpha0 = alpha0vec,m = mvec,  B = Bvec)
data(counts)
data(covset)
FixCov <- covset["Line"]
VarCov <- covset[c("Diet", "RFI", "Conca", "Concb", "RINa", "RINb", "Block", "Order","neut", "lymp", "mono", "eosi", "baso")]
lambdamax <- 10
ncores <- 20
print.progress <- F
saveall <- TRUE
dir.create(path = "RealDataOutBS", showWarnings = F)

### For HPC -----------------------------
pm1 <- proc.time()
env = as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
print(env)

B <- Bmalpha0$B[env]; alpha0 <- Bmalpha0$alpha0[env]; m <- Bmalpha0$m[env]
res <- plyr::llply(c("WN", "OWN", "RX", "ORX"),
                     function(option)FSRAnalysisBS(counts,
                                                   FixCov, VarCov,
                                                   option, B, m,
                                                   lambdamax, alpha0, ncores,
                                                   print.progress, saveall))
names(res) <- c("WN", "OWN", "RX", "ORX")
saveRDS(res, file = paste0("RealDataOutBS/RealDataOutBS_B_", B, "_m_", m, "_alpha0_", alpha0, ".rds"))
proc.time() - pm1



### For personal computer --------------------
# pm1 <- proc.time()
# for(i in 1:nrow(Bmalpha0)){
# B <- Bmalpha0$B[i]; alpha0 <- Bmalpha0$alpha0[i]; m <- Bmalpha0$m[i]
# res <- plyr::llply(c("WN", "OWN", "RX", "ORX"),
#                    function(option)FSRAnalysisBS(counts,
#                                                  FixCov, VarCov,
#                                                  option, B, m,
#                                                  lambdamax, alpha0, ncores,
#                                                  print.progress, saveall))
# names(res) <- c("WN", "OWN", "RX", "ORX")
# saveRDS(res, file = paste0("RealDataOutBS/RealDataOutBS_B_", B, "_m_", m, "_alpha0_", alpha0, ".rds"))
# }
# proc.time() - pm1
#


###-------------------------------------------
###-------------------------------------------


### Check computational time using 1 core-----
# pm <- proc.time()
# FSR1option <- FSRAnalysisBS(counts = counts, FixCov = FixCov, VarCov = VarCov,
#                         option = "WN", B = 100,
#                  m = 7, lambdamax = 10, alpha0 = 0.1, ncores = 1, print.progress = T, saveall = F)
# proc.time() - pm

###-------------------------------------------


