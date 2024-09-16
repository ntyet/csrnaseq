## This function run all simulation for total 24 scenarios of
## ModelSize (6: 0, 1, 2, 6, 7, 8)
## alpha0 (1:0.05)
## m (4: 1, 3, 5, 7)

rm(list = ls())
pm1 <- proc.time()
# devtools::install_github("ntyet/csrnaseq")
library("csrnaseq")
env <-  as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
print(env)
nGene <- 2000
nSim <- 1:100
B <- 100
lambdamax <- 10
ncores <- 20
print.progress <- F
ModelSizevec <- c(0, 1, 2, 6, 7, 8)
alpha0vec <- 0.05
mvec <- c(1, 3, 5, 7)
Scenarios <- expand.grid(alpha0 = alpha0vec, m = mvec, ModelSize = ModelSizevec)
saveall <- FALSE
savesim <- TRUE
foldername <- "revise/SimulationOutRFIprojectiononLine"
RealDataOut <- readRDS(file = paste0(foldername, "/ModelSize_",Scenarios$ModelSize[env], ".rds" ))
alpha0 <- Scenarios$alpha0[env]; m <- Scenarios$m[env]
FSRSimOut <- csrnaseq:::FSRnSimBS(RealDataOut, nGene, nSim, B, m,
                       lambdamax, alpha0, ncores, print.progress,
                       saveall, savesim, foldername)
proc.time() -pm1
