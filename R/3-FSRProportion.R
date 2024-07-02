#' Estimate False Selection Rate Proportion
#'
#' This function  estimates  the False Selection Rate Proportion (FSRP)
#' via the pseudo-variable. \eqn{\alpha_{RE,P}} is based
#' on Ratio of Expectation, \eqn{\alpha_{ER,P}} is based on
#' Expectation of Ratio.
#'
#' \deqn{\alpha_{RE, P} = \frac{U}{1+SP}}
#' where \code{SP} is the average of the total number of selected covariates among both available covariates and
#' pseudo-variables, U is the average of the number of selected pseudo-variables over B replications.
#' \deqn{\alpha_{ER, P} = \frac{U}{1+S}}
#'where \code{S} is the number of covariates selected in the set
#'of the available covariates,  U is the average of the number of selected pseudo-variables over B replications.
#'
#' @param DataBSOut the output from function \code{\link{DataBS}}.
#' @param BBootBSOut the output from function \code{\link{BBootBS}}.
#'
#' @return A list of 6 components
#' \item{B}{a numerical value representing the number of
#' sets of the pseudo-variables.}
#' \item{nPhoCov}{the number of pseudo-covariates.}
#' \item{BS}{the \code{BS} component of the output DataBSOut from function
#' \code{\link{DataBS}}, which is a vector of \code{r} value of the eliminated covariates
#' during the selection process.}
#' \item{lambda}{the \code{lambda} component of the output DataBSOut from
#' the function \code{\link{DataBS}}, which is a numerical vector of the critical thresholds
#' \eqn{\lambda}.}
#' \item{alphahat.P.RE}{\eqn{\alpha_{RE,P}} defined as above. }
#' \item{alphahat.P.ER}{\eqn{\alpha_{ER,P}} defined as above. }
#' @author Yet Nguyen \email{tienyettoan@gmail.com}
#' @examples
#' lambdamax <- 5
#'lambda <- c(seq(1, lambdamax, length = 400))
#'data(FixCov)
#'data(VarCov)
#'data(counts)
#'print.progress = FALSE
#' m <- 5
#' option <- "OWN"
#' B <- 2
#' ncores <- 1
#' svamethod <- FALSE
#' BBootBSOut<-csrnaseq:::BBootBS(B, ncores, counts[1:100,],
#' FixCov, VarCov, m,  lambda, option, print.progress)
#'DataBSOut <- csrnaseq:::DataBS(counts[1:100,], FixCov, VarCov,
#'lambda, print.progress, svamethod)
#'FSRPBSOut <- csrnaseq:::FSRPBS(DataBSOut, BBootBSOut)
#'names(FSRPBSOut)
FSRPBS <- function(DataBSOut, BBootBSOut){
  S <- DataBSOut$S # vector of selected covariates in real data w.r.t. each threshold lambda after backward selection
  lambda <- DataBSOut$lambda # vector of threshold lambda
  nLambda <- length(lambda)
  nVarCov <- length(DataBSOut$BS)
  nPhoCov <- BBootBSOut$mPhoCov
  B <- nrow(BBootBSOut$USP)
  USP <- apply(BBootBSOut$USP, 2, sum)/B
  U <- USP[1:nLambda]
  SP <- USP[-c(1:nLambda)]
  alphahat.P.RE <- U/(1+SP)
  alphahat.P.ER <- U/(1+S)
  out <- list(B = B, nPhoCov = nPhoCov, BS = DataBSOut$BS, lambda = lambda, alphahat.P.RE = alphahat.P.RE, alphahat.P.ER = alphahat.P.ER)
  out
}
#' Calculate The Optimal Cut-off Threshold Lambda Star
#'
#' This function estimates the optimal threshold \code{lambda*}
#' and also the estimates of  the number of irrelevant covariates
#' using \eqn{\alpha_{RE}} and \eqn{\alpha_{ER}}.
#'
#' @param FSRPBSOut the output of \code{\link{FSRPBS}}.
#' @param alpha0 the nominal FSR  level.
#' @return a numerical vector consisting of
#' \item{lambdahat.RE}{ the optimal \code{lambda} when using \eqn{\alpha_{RE}}.}
#' \item{lambdahat.ER}{ the optimal \code{lambda} when using \eqn{\alpha_{ER}}.}
#' \item{kU.RE}{an estimate of the number of irrelevant covariates when using \eqn{\alpha_{RE}}.}
#' \item{kU.ER}{an estimate of the number of irrelevant covariates when using \eqn{\alpha_{ER}}.}
#' @author Yet Nguyen \email{tienyettoan@gmail.com}.
#' @examples
#' lambdamax <- 10
#'lambda <- c(seq(1, lambdamax, length = 800))
#'data(FixCov)
#'data(VarCov)
#'data(counts)
#'print.progress = FALSE
#' m <- 5
#' option <- "OWN"
#' B <- 2
#' ncores <- 1
#' svamethod <- FALSE
#' BBootBSOut<-csrnaseq:::BBootBS(B, ncores, counts[1:100,],
#' FixCov, VarCov, m,  lambda, option, print.progress)
#'DataBSOut <- csrnaseq:::DataBS(counts[1:100,], FixCov, VarCov,
#'lambda, print.progress)
#' FSRPBSOut <- csrnaseq:::FSRPBS(DataBSOut, BBootBSOut)
#' alpha0 <- 0.05
#' BestLambdaOut <- csrnaseq:::BestLambdaBS(FSRPBSOut, alpha0)
#' BestLambdaOut
BestLambdaBS <- function(FSRPBSOut, alpha0=0.05){
  lambda <- FSRPBSOut$lambda
  BS <- FSRPBSOut$BS
  kT <- length(BS)
  kP <- FSRPBSOut$nPhoCov
  alphahat.P.RE <- FSRPBSOut$alphahat.P.RE
  lambdamin.RE<-min(lambda[alphahat.P.RE==max(alphahat.P.RE)])  # lambda with largest alphahat.RE
  kU.RE <- kT
  alpha0.P <- kP*alpha0/(kU.RE +kP*alpha0)# c^(0) in step 3 of the FSR Algorithm
  repeat{
    ind.RE <- (cummin(alphahat.P.RE)<=alpha0.P& lambda>=lambdamin.RE)*1
    lambdahat.RE <- lambda[(which(ind.RE==1))[1]]       # RE est. of lambda
    # or using lambdahat.RE <- lambda[(ind.RE==1) & (cumsum(ind.RE==1) == 1)] which is the first element of lambda
    # that ind.RE = 1

    if(is.na(lambdahat.RE)){
      warning(paste0("All Estimated FSRs are greater than ", alpha0, ". This happens in BestLambda function."))
      lambdahat.RE <- 10000
    }
    temp <- kT - length(BS[cummax(BS)>=lambdahat.RE])
    if(temp==kU.RE) {break}
    kU.RE <- temp
    alpha0.P <- kP*alpha0/(kU.RE +kP*alpha0)
  }

  alphahat.P.ER <- FSRPBSOut$alphahat.P.ER
  lambdamin.ER<-min(lambda[alphahat.P.ER==max(alphahat.P.ER)])  # lambda with largest alphahat.RE
  kU.ER <- kT
  alpha0.P <- kP*alpha0/kU.ER
  repeat{
    ind.ER <- (cummin(alphahat.P.ER)<=alpha0.P& lambda>=lambdamin.ER)*1
    lambdahat.ER <- lambda[(which(ind.ER==1))[1]]       # ER est. of lambda
    if(is.na(lambdahat.ER)){
      warning(paste0("All Estimated FSRs are greater than ", alpha0, ". This happens in BestLambda function."))
      lambdahat.ER <- 10000
    }
    temp <- kT - length(BS[cummax(BS)>=lambdahat.ER])
    if(temp==kU.ER) {break}
    kU.ER <- temp
    alpha0.P <- kP*alpha0/(kU.ER)
  }
  c(lambdahat.RE  = lambdahat.RE, lambdahat.ER  = lambdahat.ER, kU.RE = kU.RE, kU.ER = kU.ER)
}

#' Estimate False Selection Rate FSR
#'
#' This function estimates FSR  in two cases
#' \eqn{\alpha_{RE}} and \eqn{\alpha_{RE}}.
#' @param FSRPBSOut the output from \code{\link{FSRPBS}}.
#' @param BestLambdaOut the output from \code{\link{BestLambdaBS}}.
#' @return A data frame with 3 columns
#' \item{lambda}{the numerical vector of the critical thresholds.}
#' \item{alphahat.RE}{estimate of \eqn{\alpha_{RE}}  at lambda.}
#' \item{alphahat.ER}{estimate of \eqn{\alpha_{ER}}  at lambda.}
#' @author Yet Nguyen \email{tienyettoan@gmail.com}.
#' @examples
#' lambdamax <- 10
#'lambda <- c(seq(1, lambdamax, length = 800))
#'data(FixCov)
#'data(VarCov)
#'data(counts)
#'print.progress = FALSE
#' m <- 5
#' option <- "OWN"
#' B <- 2
#' ncores <- 1
#' svamethod <- FALSE
#' BBootBSOut<-csrnaseq:::BBootBS(B, ncores, counts[1:100,],
#' FixCov, VarCov, m,  lambda, option, print.progress)
#'DataBSOut <- csrnaseq:::DataBS(counts[1:100,], FixCov, VarCov,
#'lambda, print.progress, svamethod)
#' FSRPBSOut <- csrnaseq:::FSRPBS(DataBSOut, BBootBSOut)
#' alpha0 <- 0.05
#' BestLambdaOut <- csrnaseq:::BestLambdaBS(FSRPBSOut, alpha0)
#' FSROut <- csrnaseq:::FSR(FSRPBSOut, BestLambdaOut)
#' dim(FSR)
#' head(FSR)
FSR <- function(FSRPBSOut, BestLambdaOut){
  # FSRPOut <- FSRP(DataFSOut, BBootFSOut)
  # BestLambdaOut <- BestLambda(FSRPOut, alpha0 = alpha0)
  alphahat.RE<-pmin(1, FSRPBSOut$alphahat.P.RE/(1-FSRPBSOut$alphahat.P.RE)*BestLambdaOut["kU.RE"]/FSRPBSOut$nPhoCov)      # alphahat_RE, SP-U plays role of S
  ghat<- pmin(1, FSRPBSOut$alphahat.P.ER*BestLambdaOut["kU.ER"]/FSRPBSOut$nPhoCov)             # alphahat_ER
  out <- cbind(lambda = FSRPBSOut$lambda, alphahat.RE = alphahat.RE, ghat = ghat) # ghat should be alphahat.ER
  rownames(out) <- NULL
  data.frame(out)
}


#' The Selected Covariates Using FSR Variable Selection Method
#'
#' This function returns the selected covariates after estimating the optimal
#' threshold lambda to control FSR.
#'
#' @param DataBSOut The output of \code{\link{DataBS}} function.
#' @param BestLambdaOut The output of \code{\link{BestLambdaBS}} function.
#' @return A list consists of two elements
#' \item{BestER}{A vector of  relevance measure r
#' for the selected covariates when using \eqn{\alpha_{ER}}.}
#' \item{BestRE}{A vector of  relevance measure  r
#' for the selected covariates when using\eqn{\alpha_{RE}}.}
#' @author Yet Nguyen \email{tienyettoan@gmail.com}.
#' @examples
#' lambdamax <- 10
#'lambda <- c(seq(1, lambdamax, length = 800))
#'data(FixCov)
#'data(VarCov)
#'data(counts)
#'print.progress = FALSE
#' m <- 5
#' option <- "OWN"
#' B <- 2
#' ncores <- 1
#' BBootBSOut<-csrnaseq:::BBootBS(B, ncores, counts[1:100,],
#' FixCov, VarCov, m,  lambda, option, print.progress)
#'DataBSOut <- csrnaseq:::DataBS(counts[1:100,], FixCov, VarCov,
#'lambda, print.progress)
#' FSRPBSOut <- csrnaseq:::FSRPBS(DataBSOut, BBootBSOut)
#' alpha0 <- 0.05
#' BestLambdaOut <- csrnaseq:::BestLambdaBS(FSRPBSOut, alpha0)
#' BestCovBSOut <- csrnaseq:::BestCovBS(DataBSOut, BestLambdaOut)
#' BestCovBSOut

BestCovBS <- function(DataBSOut, BestLambdaOut){
  list(BestER = DataBSOut$BS[cummax(DataBSOut$BS)>=BestLambdaOut["lambdahat.ER"]],
       BestRE = DataBSOut$BS[cummax(DataBSOut$BS)>=BestLambdaOut["lambdahat.RE"]])
}



#' Plot false selection rate FDR estimates as a function of lambda
#'
#' This function uses \code{\link[ggplot2]{ggplot}} to make a plot of
#' the estimators \eqn{\alpha_{RE}}
#' and \eqn{\alpha_{ER}} as functions of \eqn{\lambda}.
#' @param FSROut the output of \link{FSR} function
#' @inheritParams BBootBS
#' @inheritParams BestLambdaBS
#' @return A plot of estimators \eqn{\alpha_{RE}}
#' and \eqn{\alpha_{ER}} as functions of \eqn{\lambda}.
#' @examples
#'lambdamax <- 5
#'lambda <- c(seq(1, lambdamax, length = 400))
#'data(FixCov)
#'data(VarCov)
#'data(counts)
#'print.progress = FALSE
#' m <- 7
#' option <- "OWN"
#' B <- 20
#' ncores <- 10
#' svamethod <- FALSE
#' BBootBSOut<-csrnaseq:::BBootBS(B, ncores, counts[1:3000,],
#' FixCov, VarCov, m,  lambda, option, print.progress)
#'DataBSOut <- csrnaseq:::DataBS(counts[1:3000,], FixCov, VarCov,
#'lambda, print.progress, svamethod)
#' FSRPBSOut <- csrnaseq:::FSRPBS(DataBSOut, BBootBSOut)
#' alpha0 <- 0.05
#' BestLambdaOut <- csrnaseq:::BestLambdaBS(FSRPBSOut, alpha0)
#' FSROut <- csrnaseq:::FSR(FSRPBSOut, BestLambdaOut)
#' csrnaseq:::PlotFSR(FSROut, alpha0 = .05, option, B, m)

PlotFSR <- function(FSROut, alpha0, option, B, m){
  FSROut.melt <- reshape::melt(FSROut, id.vars ="lambda", measure.vars =  c("ghat", "alphahat.RE"))
  names(FSROut.melt) <- c("lambda", "Formula", "Est.FSR")
  levels(FSROut.melt$Formula) <- c("alpha_ER", "alpha_RE")
  p <- ggplot2::ggplot(FSROut.melt, ggplot2::aes(x=lambda, y=Est.FSR, group = Formula)) +
    ggplot2::geom_line(ggplot2::aes(linetype=Formula, color = Formula), linewidth = .5) +
    ggplot2::scale_linetype_manual(values=c("dotdash", "solid"))+
    ggplot2::scale_color_manual(values=c("red", "blue"))+
    ggplot2::theme(legend.title=ggplot2::element_blank(), axis.text=ggplot2::element_text(size=14),
          axis.title=ggplot2::element_text(size=14,face="bold"))+
    # geom_hline(yintercept=alpha0) +
    ggplot2::geom_hline(ggplot2::aes(yintercept=alpha0)) +
    ggplot2::geom_text(ggplot2::aes(0,alpha0,label = alpha0, vjust = -1))+
    ggplot2::labs(title = paste(paste0(option, ", B=", B, ", m=", m, ", alpha0=", alpha0)))+
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))+
    ggplot2::xlab( expression(lambda))+
    ggplot2::ylab(expression(hat(alpha)))
  p
}



#' Differential expression analysis using voom-limma
#'
#' This function performs voom limma differential expressin analysis and returns  selected
#' outputs of the \code{\link[limma]{voom}} result including a
#' vector of p-values, transformed counts,...
#' @param counts an RNA-seq read count data
#' @param AllCov a data frame of all covariates
#' @returns A list of the following elements
#' \item{y}{ A matrix of logcpm count data}
#' \item{pvs}{A matrix of p-value vectors, one p-value vector for each covariate}
#' \item{Beta}{A matrix of regression coefficent estimates}
#' \item{Yhat}{A matrix of fitted values}
#' \item{lib.size}{A vector of ibrary size  (75quartile)}
#' \item{weights}{A matrix of observational weights}
#' \item{sigma}{sigma from output of \code{\link[limma]{lmFit}}}.
#' \item{design}{the design matrix}
#' @examples
#' data(FixCov)
#' data(VarCov)
#' data(counts)
#' AllCov <- cbind(FixCov, VarCov)
#' VoomPvOut <- csrnaseq:::VoomPv(counts, AllCov)
#' names(VoomPvOut)
VoomPv <- function(counts, AllCov){
  dm <- model.matrix(formula(paste0("~", paste0(names(AllCov), collapse = "+"))), data = AllCov)
  colnames(dm)[1] <- "Intercept"
  vout <- limma::voom(counts = counts, design = dm, lib.size = apply(counts, 2, quantile, .75), plot = F)
  fit <- limma::lmFit(vout)
  pvs <- plyr::ldply(1:ncol(AllCov), function(i){
    if(is.factor(AllCov[,i]) | is.character(AllCov[,i])) {
      ct <- paste0(grep(paste0(names(AllCov)[i]),  x = colnames(dm), value = T), collapse = ",  ")
    }else{
      ct <- paste0(grep(paste0(names(AllCov)[i], "$"),  x = colnames(dm), value = T), collapse = ",  ")
    }
    C.matrix <- eval(parse(text=paste0("limma::makeContrasts(",  ct, ",levels = dm)")))
    fit1 <- limma::contrasts.fit(fit, contrasts =C.matrix)
    fit1 <- limma::eBayes(fit1)
    tt <- limma::topTable(fit1, sort ="none", n = Inf)
    pv <- tt$P.Value
  })
  pvs <- t(pvs)
  colnames(pvs) <- colnames(AllCov)
  res <- list(y = vout$E,
              pvs = pvs,
              Beta = fit$coef,
              Yhat = fit$coef%*%t(vout$design),
              lib.size = vout$targets$lib.size,
              weights = vout$weights,
              sigma = fit$sigma,
              design = vout$design)
  res
}

#' Histogram of p-values of covariates
#'
#' This function plots a histogram of the p-values for all covariates,
#' the primary ones and those subjected to variable selection.
#' @param VoomPvOut The output of \code{\link{VoomPv}}
#' @inheritParams BBootBS
#' @inheritParams BestLambdaBS
#' @param ErrType either "ER" (for \eqn{\alpha_{ER}})  or "RE" (for \eqn{\alpha_{RE}})
#' @param lambdamax the maximum value of  lambda vector
#' @return A plot of histogram of p-values for the covariates of the final model
#' @examples
#' lambdamax <- 5
#' lambda <- c(seq(1, lambdamax, length = 400))
#' data(FixCov)
#' data(VarCov)
#' data(counts)
#' print.progress = FALSE
#' m <- 5
#' option <- "OWN"
#' B <- 5
#' ncores <- 1
#' alpha0 <- 0.05
#' DataBSOut <- DataBS(counts = counts[1:1000,],
#' FixCov = FixCov, VarCov=VarCov, lambda = lambda,
#' print.progress = print.progress, svamethod = svamethod)
#' BBootBSOut <- BBootBS(B= B, ncores = ncores, counts = counts[1:1000,],
#' FixCov, VarCov = DataBSOut$VarCov, m= m, lambda = lambda,
#' option = option, print.progress = print.progress)
#' FSRPBSOut <- FSRPBS(DataBSOut=DataBSOut, BBootBSOut = BBootBSOut)
#' BestLambdaOut <- BestLambdaBS(FSRPBSOut, alpha0)
#' FSROut <- FSR(FSRPBSOut, BestLambdaOut)
#' PlotFSROut <- PlotFSR(FSROut = FSROut, alpha0 = alpha0, option = option, B = B, m = m)
#' BestCovOut <- BestCovBS(DataBSOut = DataBSOut, BestLambdaOut = BestLambdaOut)
#' VoomPvOutER <- VoomPv(counts = counts,
#'                       AllCov = cbind(FixCov, (DataBSOut$VarCov)[names(BestCovOut$BestER)]))
#'PlotVoomPvOutER <-  PlotVoomPv(VoomPvOut = VoomPvOutER,
#'option = option, B = B, m = m, ErrType = "ER", alpha0 = alpha0)
#'PlotVoomPvOutER
#'VoomPvOutRE <- VoomPv(counts = counts,
#'AllCov = cbind(FixCov, (DataBSOut$VarCov)[names(BestCovOut$BestRE)]))
#'PlotVoomPvOutRE <-  PlotVoomPv(VoomPvOut = VoomPvOutRE,
#'option = option, B = B, m = m, ErrType = "RE", alpha0 = alpha0)
#'PlotVoomPvOutRE

PlotVoomPv <- function(VoomPvOut, option= "OWN", B = 100, m = 3, ErrType= "ER", alpha0 = 0.05){
  VoomPvOutp <- VoomPvOut$pvs
  VoomQvOut <- apply(VoomPvOutp, 2, function(x)jabes.q(x))
  DEGs <- apply(VoomQvOut <= .05, 2, sum)
  VoomPvOut.melt <- reshape2::melt(VoomPvOutp, measure.vars = .)
  names(VoomPvOut.melt)[2:3] <- c("Covariate", "pvalue")
  levels(VoomPvOut.melt$Covariate)<- paste(levels(VoomPvOut.melt$Covariate),  DEGs[levels(VoomPvOut.melt$Covariate)], sep = ", q.05 = ")
  p <- ggplot2::ggplot(data = VoomPvOut.melt)+
    ggplot2::geom_histogram(mapping = aes(x = pvalue), breaks=seq(0,1,by=0.05))+
    ggplot2::facet_wrap(~ Covariate,  scales = "free_y") +
    ggplot2::labs(title = paste0("Selected Model, alpha.", ErrType, ", ",option,  ", B =", B, ", m = ", m, ", alpha0=", alpha0))+
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
  p
}

#' Variable Selection in RNA-seq Differential Expression
#' Analysis Using Pseudo-Variables
#'
#' This function performs the proposed backward variable selection
#' using pseudo-variables for RNA-seq differential expression analysis. The idea
#' is to select the most relevant covariates such that the false selection rate
#' is  below some pre-specified threshold.
#' @inheritParams DataBS
#' @inheritParams BBootBS
#' @inheritParams BestLambdaBS
#' @inheritParams PlotFSR
#' @inheritParams PlotVoomPv
#' @param saveall logical. If it is TRUE, then save all output, otherwise, just save
#' the best covariate BestCovOut, which is the output from \link{BestCovBS}
#' @return A list of the following components if \code{saveall = TRUE}
#' \item{DataBSOut}{ Output of DataBS function, see \link{DataBS}}
#' \item{BBootBSOut}{Output of BBootBS function, see \link{BBootBS}}
#' \item{FSROut}{Output of \link{FSR}}
#' \item{PlotFSROut}{Output plot of \link{PlotFSR}}
#' \item{BestLambdaOut}{Output of \link{BestLambdaBS}}
#' \item{BestCovOut}{Output of \link{BestCovBS}}
#' \item{VoomPvOutER}{Output of \link{PlotVoomPv} for ErrType ER}
#' \item{PlotVoomPvOutER}{Output plot of \link{PlotVoomPv} for ErrType ER}
#' \item{option}{Method generates pseudo-variables: WN, RX, OWN, ORX}
#' \item{B}{Number of sets of pseudo-variables}
#' \item{m}{Number of pseudo-variables}
#' @export
#' @examples
#' data(counts)
#' data(FixCov)
#' data(VarCov)
#' option <- "OWN"
#' B <- 5
#' m <- 5
#' lambdamax <- 5
#' alpha0 <- 0.05
#' ncores <- 1
#' print.progress <- FALSE
#' saveall <- TRUE
#' FSRAnalysisBSOut <- FSRAnalysisBS(counts, FixCov, VarCov, option, B, m, lambdamax, alpha0, ncores, print.progress, saveall)
#' names(FSRAnalysisBSOut)

FSRAnalysisBS <- function(counts, FixCov, VarCov, option = "OWN", B= 100, m = 3, lambdamax= 10, alpha0 = .05, ncores, print.progress, saveall = TRUE){
  cat("Option=", option, ", B=", B, ", m=",m, ", alpha0 = ", alpha0, "\n")
  pm <- proc.time()
  lambda <- c(seq(1, lambdamax, length = 800))
  DataBSOut <- DataBS(counts = counts, FixCov = FixCov, VarCov=VarCov, lambda = lambda, print.progress = print.progress)
  BBootBSOut <- BBootBS(B= B, ncores = ncores, counts = counts, FixCov, VarCov = DataBSOut$VarCov, m= m, lambda = lambda, option = option, print.progress = print.progress)
  FSRPBSOut <- FSRPBS(DataBSOut=DataBSOut, BBootBSOut = BBootBSOut)
  BestLambdaOut <- BestLambdaBS(FSRPBSOut, alpha0)
  FSROut <- FSR(FSRPBSOut, BestLambdaOut)
  if(saveall){
    BestCovOut <- BestCovBS(DataBSOut = DataBSOut, BestLambdaOut = BestLambdaOut)
    PlotFSROut <- PlotFSR(FSROut = FSROut, alpha0 = alpha0, option = option, B = B, m = m)
    VoomPvOutER <- VoomPv(counts = counts,
                        AllCov = cbind(FixCov, (DataBSOut$VarCov)[names(BestCovOut$BestER)]))
    PlotVoomPvOutER <-  PlotVoomPv(VoomPvOut = VoomPvOutER, option = option, B = B, m = m, ErrType = "ER", alpha0 = alpha0)
    VoomPvOutRE <- VoomPv(counts = counts,
                        AllCov = cbind(FixCov, (DataBSOut$VarCov)[names(BestCovOut$BestRE)]))
    PlotVoomPvOutRE <-  PlotVoomPv(VoomPvOut = VoomPvOutRE, option = option, B = B, m = m, ErrType = "RE", alpha0 = alpha0)
    res <- list(DataBSOut = DataBSOut, BBootBSOut = BBootBSOut,
              FSROut = FSROut, PlotFSROut = PlotFSROut,
              BestLambdaOut = BestLambdaOut,
              BestCovOut = BestCovOut,
              VoomPvOutER = VoomPvOutER, PlotVoomPvOutER = PlotVoomPvOutER,
              VoomPvOutRE = VoomPvOutRE, PlotVoomPvOutRE = PlotVoomPvOutRE,
              option = option, B = B, m = m)}else{
                BestCovOut <- BestCovBS(DataBSOut = DataBSOut, BestLambdaOut = BestLambdaOut)
                res <- list(BestCovOut = BestCovOut)
              }

  pm <- proc.time()-pm
  cat("Option=", option, ", B=", B, ", m=",m,  " runs in ", pm[3], " s\n")
  res
}
