#' Simulate Count Data from RFI RNA-seq Data
#'
#'  This function simulates count datasets that contain both EE and DE genes
#'   with respect to the
#'  primary variable and each of the relevant covariates. For each simulation scenario,
#'  as true parameters to simulate new data, we used the precision weights,
#'  the scaled error variances and the partial regression coefficient estimates
#'  from the `voom-limma` fit \link{VoomPV} of the corresponding model to the RFI RNA-seq data,
#'  except that we set partial regression coefficients on each variable to zero
#'  for a subset of genes to allow simulation of EE genes. More specifically,
#'  for each variable \eqn{j} (either relevant covariates or the primary variable
#'  Line), the \eqn{\hat{G}_0^{(j)}} least significant partial regression
#' coefficients were set to zero, where  \eqn{\hat{G}_0^{(j)}}  is the estimated
#' number of the \eqn{j}-covariate partial regression coefficients equal to zero
#' when the histogram-based method of Nettleton et al. 2006 is applied to the
#' \eqn{j}-variable's \eqn{p}-values
#' from the \link{VoomPv} fit of the corresponding model to the RFI RNA-seq data.
#' This strategy yielded a parameter vector (consisting of a scaled error
#' variance, precision weights, and  partial regression coefficients) for
#' each of 12280 genes. To simulate any particular dataset for a given set
#' of truly relevant covariates (either 0, 1, 2, 6, 7, or 8 relevant covariates),
#' we randomly sampled \code{nGene} (default value is 2000) gene parameter vectors.
#' The selected parameters
#' and the respective  values of the primary variable and the relevant
#' covariates for the 31 pigs were used to simulate a \code{nGene} x 31 dataset
#' of read counts following the inverse steps of log2 count transformation and
#' the linear model used in the main paper.
#'
#' @param RealDataOut the output of voom-limma differential expression analysis
#' implemented in \link{VoomPv} for RFI RNA-seq data. This analysis is performed in the
#' script 2-RSRBS_ModelSize_Simulation.R  in analysis folder.
#' @param nGene number of genes in simulated count data
#' @param nrep replication index, used in set.seed within this function
#' @return A list of 3 components
#' \item{SimCnt}{simulated count data of nGene rows and 31 columns}
#' \item{IndSample}{a vector of original row ID of the RFI RNA-seq data for nGene simulated genes}
#' \item{EEGene}{a list of vectors of original row ID of the RFI RNA-seq data for  null genes of each covariate}
#' @author Yet Nguyen \email{tienyettoan@gmail.com}
#' @examples
#'\dontrun{
#' RealDataOut <- readRDS("./analysis/RealDataOutBS/ModelSize_Simulation_6.rds")
#' nGene <- 2000
#' nrep <- 1
#' SimCnt <- SimCounts(RealDataOut, nGene,nrep)
#' }

SimCounts <- function(RealDataOut, nGene = 2000, nrep){
  set.seed(2017+nrep)
  # res <- ((RealDataOut[[1]])[[option]])[[paste0("VoomPvOut",ErrType)]]
  res <- RealDataOut$VoomPvOutER
  EEGene <- list()
  AllCov <- cbind(RealDataOut$FixCov, RealDataOut$VarCov)
  ConCov <- AllCov[colnames(res$pvs)]
  lib.size <- res$lib.size; design <- res$design; Beta0 <- res$Beta;
  sigma <- res$sigma;  weights <- res$weights; VarMat0 <- sigma^2*1/weights
  m0 <- apply(res$pvs, 2, function(x)estimate.m0(x))
  for(i in 1:ncol(ConCov)){
    if(is.factor(ConCov[,i]) | is.character(ConCov[,i])) {
      ct <- grep(paste0(names(ConCov)[i]),  x = colnames(design), value = T)
    }else{
      ct <- grep(paste0(names(ConCov)[i], "$"),  x = colnames(design), value = T)
    }
    EEGene[[i]] <- (sort(res$pvs[,i], decreasing = T,  index.return = T))$ix[1:m0[i]]
    Beta0[EEGene[[i]], ct] <- 0
  }
  names(EEGene) <- colnames(res$pvs)
  IndSample <- sample(1:nrow(Beta0), size = nGene)
  Beta <- Beta0[IndSample,]; VarMat <- VarMat0[IndSample,]
  sim.counts <- vapply(1:nrow(Beta), function(i){
    y <- design%*%Beta[i,]+ MASS::mvrnorm(n = 1, mu = rep(0, nrow(design)), Sigma = diag(VarMat[i,])) # library(MASS)
    unname(round(2^y*(lib.size+1)/10^6 )[ , drop = T])
  }, FUN.VALUE = rep(1.0, nrow(design)))
  SimCnt <- t(sim.counts)
  list(SimCnt = SimCnt, IndSample = IndSample, EEGene = EEGene)
}

#' Calculate evaluation metrics on differential expression analysis
#'
#' This function calculate the number of true positive (NTP),
#' the number of declared differentially expressed genes (R),
#' number of false positive (V), False Discovery Proportion (FDP),
#' and the partial area under ROC
#' curve with false positive rate less than 0.05 (PAUC) from the p-values vector
#' of the main factor of interest when analyzing simulated count data.
#'
#'
#' @param  p vector of p-values
#' @param lab vector of 0 and 1. 0 means true EE, 1 means true DE
#' @param method the methods that handles available covariates
#' (WN_RE, WN_ER, RX_ER, RX_RE, OWN_RE, OWN_ER, ORX_ER, ORX_RE, Oracle, OnlyLine, All, OldBS)
#' @return A vector of 5 elements: NTP, R, V, FDP, PAUC
#' @examples
#' sim <- readRDS("./analysis/SimulationOut/ModelSize_6_nGene_2000_B_100_m_7_alpha0_0.05/nrep_2.rds")
#' AllCov <- cbind(sim$SimCntOut$FixCov, sim$SimCntOut$VarCov[,sim$SimSelCov.ER$ORX.ER])
#' counts <- sim$SimCnt$SimCnt
#' p <- VoomPv(counts, AllCov)$pvs[,"Line"]
#' lab <-as.factor(as.numeric(!(sim$SimCnt$IndSample%in%sim$SimCnt$EEGene$Line)))
#' method <- "ORX.ER"
#' DEAevalOut <- DEAeval(p, lab, method)
#' DEAeval

DEAeval <- function(p, lab,  method){
  roc.out <- AUC::roc(1-p, lab) # plot(roc.out)
  roc.ind <- sum(roc.out$fpr<=.05)
  roc.min <- roc.out$cutoffs[roc.ind]
  pauc <- AUC::auc(roc.out, min =roc.min)
  qv <- jabes.q(p)
  R<- sum(qv <=.05)
  V <- sum((qv <= .05)*(!as.numeric(as.character(lab))))
  FDP <- V/max(1, R)
  S = R - V
  out <- c(NTP = S, R = R, V = V,FDP = FDP, PAUC = pauc)
  names(out) <- paste0("DEA.", names(out), ".", method)
  out
}


#' Analyze one simulated dataset
#'
#' This function perform variable selection and differential expression analysis
#' for one simulated count data. The
#' output is saved in the directory `SimulationOut/...`. The result is
#' list of 4 components, each is a list of 2 vectors of selected covariates
#' under alpha_RE and alpha_ER FSR.
#' @inheritParams FSRAnalysisBS
#' @inheritParams SimCounts
#' @param savesim logical. If TRUE, save all the output including SimCounts to
#' `SimulationOut/ModelSize_*_nGene_*_B_*_m_*_alpha0_*/nrep_*.rds`
#' @examples
#'\dontrun{
#' RealDataOut <- readRDS("./analysis/RealDataOutBS/ModelSize_6.rds")
#' ModelSize <- 6
#' nGene <- 2000
#' nrep <- 1
#' B <- 4
#' m <- 4
#' lambdamax <- 5
#' alpha0 <- 0.05
#' ncores <- 4
#' print.progress <- FALSE
#' saveall <- FALSE
#' savesim <- TRUE
#' FSR1SimBSOut <- FSR1SimBS(RealDataOut, nGene, nrep, B, m, lambdamax, alpha0, ncores, print.progress, saveall, savesim)
#' }
FSR1SimBS <- function(RealDataOut, nGene, nrep, B, m, lambdamax, alpha0, ncores, print.progress, saveall = FALSE, savesim = TRUE){
  # cat("Sim = ", nrep, "\n")
  FixCov <- RealDataOut$FixCov
  VarCov <- RealDataOut$VarCov
  TrueSelCov <- names(RealDataOut$BestCovOut$BestER)
  SimCnt <- SimCounts(RealDataOut, nGene,nrep)

  option <- c("RX", "ORX", "WN", "OWN")

  SimCntOut <- plyr::llply(c("WN", "OWN", "RX", "ORX"),
                     function(option)FSRAnalysisBS(counts = SimCnt$SimCnt,
                                                   FixCov, VarCov,
                                                   option, B, m,
                                                   lambdamax, alpha0, ncores,
                                                   print.progress,
                                                   saveall))
  names(SimCntOut) <- c("WN", "OWN", "RX", "ORX")
#
#   SimCntOut <- FSRAnalysisBSAll(counts= SimCnt$SimCnt, FixCov, VarCov,
#                                 B, m, lambdamax, alpha0, ncores, print.progress, saveall)
  SimSelCov.ER <- plyr::llply(c("WN", "OWN", "RX", "ORX"),
                              function(option)names(SimCntOut$res[[option]]$BestCovOut$BestER))
  names(SimSelCov.ER) <- c("WN.ER", "OWN.ER", "RX.ER", "ORX.ER")
  SimSelCov.RE <- plyr::llply(c("WN", "OWN", "RX", "ORX"),
                              function(option)names(SimCntOut$res[[option]]$BestCovOut$BestRE))
  names(SimSelCov.RE) <- c("WN.RE", "OWN.RE", "RX.RE", "ORX.RE")

  #FSROldBS result
  SimSelCov.OldBS <- names(jabes.bs(counts = SimCnt$SimCnt, FixCov = FixCov, VarCov = VarCov, print.progress)$BestCov)

  # Calculate the empirical FSR, S, U, R for FSRBS

  S.ER <-  plyr::laply(SimSelCov.ER, function(x) length(intersect(TrueSelCov, x)))
  R.ER <-  plyr::laply(SimSelCov.ER, function(x)length(x))
  U.ER <- R.ER - S.ER
  FSP.ER <- U.ER/(R.ER+1) # i
  res1 <- c(S.ER, R.ER, U.ER, FSP.ER)
  names(res1) <- paste(rep(c("S", "R", "U", "FSP"), each = 4), rep(names(SimSelCov.ER), 4), sep = ".")
  S.RE <-  plyr::laply(SimSelCov.RE, function(x) length(intersect(TrueSelCov, x)))
  R.RE <-  plyr::laply(SimSelCov.RE, function(x)length(x))
  U.RE <- R.RE - S.RE
  FSP.RE <- U.RE/(R.RE+1) # i
  res2 <- c(S.RE, R.RE, U.RE, FSP.RE)
  names(res2) <- paste(rep(c("S", "R", "U", "FSP"), each = 4), rep(names(SimSelCov.RE), 4), sep = ".")

  # FSR for OldBS
  R <- length(SimSelCov.OldBS) # S in the Wu2007 paper
  S <- length(intersect(TrueSelCov, SimSelCov.OldBS))# I in the Wu2007 paper
  U <- R - S
  FSP <- U/(R+1)
  res3 <- c(S.OldBS=S, R.OldBS=R, U.OldBS=U, FSP.OldBS=FSP)
  FSROut <- c(res1, res2, res3)
  names(FSROut) <- paste0("FSR.", names(FSROut))

  # DEA result for 1Sim
  lab <-as.factor(as.numeric(!(SimCnt$IndSample%in%SimCnt$EEGene$Line)))
  # OnlyLine
  OnlyLine <- DEAeval(p = VoomPv(counts = SimCnt$SimCnt, AllCov = FixCov)$pvs[,"Line"], lab = lab,  method = "OnlyLine")
  # AllCov
  All <- DEAeval(p = VoomPv(counts = SimCnt$SimCnt, AllCov = cbind(FixCov, VarCov))$pvs[,"Line"], lab = lab,  method = "All")
  #Oracle
  TrueCov <- VarCov[TrueSelCov]
  Oracle <- DEAeval(p = VoomPv(counts = SimCnt$SimCnt, AllCov = cbind(FixCov,TrueCov))$pvs[,"Line"], lab = lab,  method = "Oracle")
  #jabes.bs
  BestCov <- jabes.bs(counts = SimCnt$SimCnt, FixCov = FixCov, VarCov = VarCov, print.progress = F)$BestCov
  OldBS <- DEAeval(p = VoomPv(counts = SimCnt$SimCnt, AllCov = cbind(FixCov,BestCov))$pvs[,"Line"], lab = lab, method = "OldBS")
  option <- c("RX", "ORX", "WN", "OWN")
  EstType <- c("ER", "RE")
  oE <- expand.grid(option = option, EstType = EstType)
  DEAevalOut <- lapply(1:nrow(oE), function(i){
    option <-as.character(oE$option[i]); EstType <- as.character(oE$EstType[i])
    BestCovOut <- VarCov[names(SimCntOut$res[[option]]$BestCovOut[[paste0("Best", EstType)]])]
    DEAeval(p = VoomPv(counts = SimCnt$SimCnt, AllCov = cbind(FixCov,BestCovOut))$pvs[,"Line"],lab = lab,  method = paste0( option, ".", EstType))
  })
  DEAevalOut <- do.call("c", DEAevalOut)
  DEAevalOut <- c(DEAevalOut, OnlyLine, All, Oracle, OldBS)

  SimOut <- list(SimCnt = SimCnt,
                 SimCntOut = SimCntOut,
                 TrueSelCov = TrueSelCov,
                 SimSelCov.ER = SimSelCov.ER,
                 SimSelCov.RE = SimSelCov.RE,
                 SimSelCov.OldBS = SimSelCov.OldBS,
                 FSROut = FSROut,
                 DEAevalOut = DEAevalOut,
                 FixCov = FixCov,
                 VarCov = VarCov)
  ModelSize <- length(RealDataOut$BestCovOut$BestER)
  if(savesim){
    SimPath <- paste0("SimulationOut/ModelSize_", ModelSize, "_nGene_", nGene, "_B_", B, "_m_", m,  "_alpha0_", alpha0)
    dir.create(path = SimPath, showWarnings = F, recursive = T)
    saveRDS(SimOut, file = paste0(SimPath, "/nrep_", nrep, ".rds"))
  }

  SimSelCov <- list(TrueSelCov = TrueSelCov,
                    SimSelCov.ER = SimSelCov.ER,
                    SimSelCov.RE = SimSelCov.RE,
                    SimSelCov.OldBS = SimSelCov.OldBS,
                    FSROut = FSROut,
                    DEAevalOut = DEAevalOut)
  SimSelCov
}


#' Analyze multiple simulated datasets
#'
#' This function performs the analysis \link{FSR1SimBS} for length(nSim) simulated
#' datasets. For each nrep in nSim, it will call \link{FSR1SimBS}, then extract the number of
#' selected covariates (R), the number of selected relevant covariates  (S),
#' and the number of selected irrelevant covariates (U). Here, R = S + U.
#' Then it calculates the empirical False Selection Rate using either \eqn{\alpha_{RE}} or
#' \eqn{\alpha{ER}}.
#'
#' @inheritParams FSR1SimBS
#' @param nSim a vector of  n replication indices, in the form  1:n
#' @return A list of 2 components
#' \item{TrueSelCov}{a vector of the name of the relevant covariates used in simulation}
#' \item{nSimSelCov}{a data frame of length(nSim) row, each coresponds to the set of FSP and FDP of the considered methods}
#' @examples
#'\dontrun{
#' RealDataOut <- readRDS("./analysis/RealDataOutBS/ModelSize_6.rds")
#' nGene <- 2000
#' nSim <- 1:3
#' B <- 4
#' m <- 2
#' lambdamax <- 5
#' alpha0 <- 0.05
#' ncores <- 4
#' print.progress <- FALSE
#' SimPath <- paste0("SimulationOut/ModelSize_", ModelSize, "_nGene_", nGene, "_B_", B, "_m_", m,  "_alpha0_", alpha0)
#' FSRnSimBSOut <- FSRnSimBS(RealDataOut, nGene, nSim, B, m, lambdamax, alpha0, ncores, print.progress, saveall = FALSE, savesim = TRUE)
#' }
FSRnSimBS <- function(RealDataOut, nGene, nSim, B, m, lambdamax, alpha0, ncores, print.progress,saveall = FALSE, savesim = TRUE){
  nSimSelCov <- plyr::ldply(nSim, function(nrep){
    FSR1SimBSOut <- FSR1SimBS(RealDataOut, nGene, nrep, B, m, lambdamax, alpha0, ncores, print.progress, saveall, savesim)
    c(FSR1SimBSOut$FSROut, FSR1SimBSOut$DEAevalOut)
  })
  res <- list(TrueSelCov = names(RealDataOut$BestCovOut$BestER),
              nSimSelCov = nSimSelCov)
  ModelSize <- length(RealDataOut$BestCovOut$BestER)
  SimPath <- paste0("SimulationOut/ModelSize_", ModelSize, "_nGene_", nGene, "_B_", B, "_m_", m,  "_alpha0_", alpha0)
  dir.create(path = SimPath, showWarnings = F, recursive = T)
  saveRDS(res, file = paste0(SimPath, "/nSim_", min(nSim), "_", max(nSim), ".rds"))
  saveRDS(res, file = paste0(SimPath, "_nSim_", min(nSim), "_", max(nSim), ".rds"))
  res
}


