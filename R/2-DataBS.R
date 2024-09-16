#' Backward Selection with Critical Thresholds
#'
#' This function performs backward selection algorithm.
#' The differences between \code{DataBS} and \code{bs} are that the
#' former returns also the number of selected covariates \code{S} with respect
#' to a critical threshold\code{lambda}, which are described as in
#' the \code{Return} section below.
#' @inheritParams bs
#' @param lambda a vector of critical thresholds for the relevance measure.
#' The typical \code{lambda} has range from 1 to 10.
#' @return A list of 6 components
#' \item{BS}{a vector of \code{r} value of the eliminated covariates
#' during the selection process.}
#' \item{S}{A vector of the number of selected covariates w.r.t  the
#' critical threshold vector \code{lambda}, i.e, S[i] is the largest number
#' of selected covariates whose their r values > lambda[i].}
#' \item{lambda}{ a numerical vector of the critical thresholds.}
#' \item{FixCov}{it is the same as the input argument.}
#'\item{VarCov}{it is the same as the input argument.}
#'
#' @examples
#'lambdamax <- 10
#'lambda <- c(seq(1, lambdamax, length = 400))
#'data(FixCov)
#'data(VarCov)
#'data(counts)
#'print.progress = FALSE
#'DataBSOut <- csrnaseq:::DataBS(counts, FixCov, VarCov, lambda, print.progress)
#'names(DataBSOut)
DataBS <- function(counts, FixCov, VarCov, lambda, print.progress = FALSE){
  BSout <- bs(counts = counts, FixCov = FixCov, VarCov = VarCov, print.progress = print.progress)
  BS <- BSout$WorstP5
  S <- vapply(1:length(lambda), function(i)sum(cummax(BS) > lambda[i]), FUN.VALUE = 1.0)
  return(list(BS = BS, S = S, lambda = lambda, FixCov = FixCov, VarCov = VarCov))
}


#' Pseudo Variable Generation
#'
#' This function generates pseudo-variables that are used in the
#' False Selection Rate variable selection method.
#' @param dm0  full model design matrix, the first column is the \code{intercept}
#' @param m  number of generated pseudo-variables
#' @param option the method to generate pseudo-variables
#' This is either one of these \code{"WN", "RX", "OWN", "ORX"}.
#' @return a data frame of \code{n} generated pseudo-variables.
#' @examples
#' data(FixCov)
#' data(VarCov)
#' AllCov <- cbind(FixCov, VarCov)
#' dm0 <- model.matrix(~., data = AllCov)
#' m <- 5
#' option <- "OWN"
#' PhoVarOut <- csrnaseq:::PhoVar(dm0, m, option)
#' head(PhoVarOut)
PhoVar <- function(dm0,  m, option){
  ProjMat <- function(dm0){
    Hx <- dm0%*%chol2inv(chol(crossprod(dm0, dm0)))%*%t(dm0)
    Ix <- diag(nrow(dm0))
    Ix - Hx
  }
  Px <- ProjMat(dm0)
  if(option =="WN"){
    Nois <- t(MASS::mvrnorm(n = m, mu = rep(0, nrow(dm0)), Sigma = diag(nrow(dm0))))
    if(m==1) Nois <- data.matrix(MASS::mvrnorm(n = m, mu = rep(0, nrow(dm0)), Sigma = diag(nrow(dm0))))
  }
  else if(option == "OWN"){
    Nois0 <- t(MASS::mvrnorm(n = m, mu = rep(0, nrow(dm0)), Sigma = diag(nrow(dm0))))
    if(m==1) Nois0 <- data.matrix(MASS::mvrnorm(n = m, mu = rep(0, nrow(dm0)), Sigma = diag(nrow(dm0))))
    Nois <- Px%*%Nois0
  } else if (option == "RX"){
    Nois <- data.matrix((data.frame(dm0)[-1])[sample(nrow(dm0)), sample(ncol(data.frame(dm0)[-1]), m)])
  } else {
    Nois0 <- data.matrix((data.frame(dm0)[-1])[sample(nrow(dm0)), sample(ncol(data.frame(dm0)[-1]), m)])
    Nois <- Px%*%Nois0
  }
  colnames(Nois) <- NULL
  Nois <- data.frame(Nois)
  Nois
}


#' Backward Selection for the Covariates Including Pseudo-Variables
#'
#' This function performs backward selection on the covariates
#' including available measured covariates and pseudo-variables and returns
#' a vector of the number of selected pseudo-variables U  and a vector
#' of the total number of selected covariates SP for the critical
#' threshold \code{lambda}.
#' @inheritParams DataBS
#' @inheritParams PhoVar
#' @param nrep replication index for each pseudo-variables set.
#' @return a vector of 2*length(\code{lambda}) elements.
#' Its first length(\code{lambda})
#' elements are U(\eqn{\lambda}), the number of pseudo-variables selected at
#' the critical threshold \eqn{\lambda}.
#' Its last length(\eqn{\lambda}) elements are SP(\eqn{\lambda}), the total number
#' of selected covariates including pseudo-variables  at the critical threshold \eqn{\lambda}.
#' @author Yet Nguyen \email{tienyettoan@gmail.com}.
#' @examples
#'lambdamax <- 10
#'lambda <- c(seq(1, lambdamax, length = 800))
#'data(FixCov)
#'data(VarCov)
#'data(counts)
#'print.progress = FALSE
#' m <- 5
#' nrep <- 1
#' option <- "OWN"
#' BootBSOut <- csrnaseq:::BootBS(counts[1:100,], FixCov, VarCov, m,  nrep, lambda, option, print.progress)
#' length(BootBSOut)
#' head(BootBSOut)

BootBS <- function(counts, FixCov, VarCov, m=3, nrep, lambda, option = "OWN", print.progress = T){
  set.seed(nrep)
   # cat("nrep = ", nrep, "\n")
  AllCov0 <- cbind(FixCov, VarCov)
  dm0 <- stats::model.matrix(stats::formula(paste0("~", paste0(names(AllCov0), collapse = "+"))), data = AllCov0)
  colnames(dm0)[1] <- "Intercept"
  if(Matrix::rankMatrix(dm0) + m > nrow(dm0)) stop("try smaller m because at least one degree of freedom is needed for estimating residual variance.")
  repeat{
    Nois <- PhoVar(dm0 = dm0, m = m, option = option)
    VarCov2 <- cbind(VarCov, Nois)
    AllCov1 <- cbind(FixCov, VarCov2)
    dm1 <- stats::model.matrix(stats::formula(paste0("~", paste0(names(AllCov1), collapse = "+"))), data = AllCov1)
    if(limma::is.fullrank(dm1)) break
  }
  bsout <- bs(counts = counts, FixCov = FixCov, VarCov = VarCov2, print.progress = print.progress)
  pv <- bsout$WorstP5
  pvmax <- cummax(pv)
  vor <- names(pv)
  U <- vapply(1:length(lambda), function(i)sum(vor%in%names(Nois) & pvmax>lambda[i]), FUN.VALUE = 1.0)
  SP <- vapply(1:length(lambda), function(i)sum( pvmax>lambda[i]), FUN.VALUE = 1.0)
  res <- c(U = U, SP = SP)
  # cat("nrep = ", nrep, " done!\n")
  res
}

#'  Repeatedly Run Backward Selection for the Covariates Including Pseudo-Variables
#'
#' This function performs \code{\link{BootBS}} \code{B} times to
#' obtain a matrix of B x 2*length(lambda), each row is a vector of U and S, which are the output of \link{BootBS} function
#' for each replication.
#' @inheritParams BootBS
#' @param B a numerical value representing the number of
#' sets of the pseudo-variables. Usually, \code{B = 100} is
#' good enough for this method.
#' @param ncores number of cores to use for parallel computing via
#' \code{\link[parallel]{mclapply}}.
#' @return A list of two components
#' \item{USP}{a matrix of B x 2*length(lambda) where each
#' row includes  \code{U} and \code{SP}, which is the output of \link{BootBS} for
#' each replication.}
#' \item{mPhoCov}{number of generated pseudo-variables}
#' @examples
#' lambdamax <- 5
#'lambda <- c(seq(1, lambdamax, length = 800))
#'data(FixCov)
#'data(VarCov)
#'data(counts)
#'print.progress <- FALSE
#' m <- 5
#' option <- "OWN"
#' B <- 2
#' ncores <- 1
#' BBootBSOut<-csrnaseq:::BBootBS(B, ncores, counts[1:100,],
#' FixCov, VarCov, m,  lambda, option, print.progress)
#' names(BBootBSOut)
#' dim(BBootBSOut$USP)
#' BBootBSOut$m
BBootBS <- function(B = 100, ncores = 1, counts, FixCov,
                    VarCov, m = 3, lambda, option = "OWN",
                    print.progress = FALSE){
  USP <- parallel::mclapply(1:B, function(nrep)BootBS(counts, FixCov,
                                            VarCov, m=m, nrep,
                                            lambda, option = option,
                                            print.progress),
                  mc.cores = ncores)
  USP <- do.call("rbind", USP)
  list(USP=USP, mPhoCov = m)
}
