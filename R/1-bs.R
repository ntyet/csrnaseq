#' #' Calculating surrogate variables on top of the existing covariates
#' #'
#' #' This function calculates the surrogate variables on top of the
#' #' existing the covariates that includes just fixed covariates or
#' #'  both the fix covariates
#' #' that are not subjected to variable selection and the nuisance
#' #' covariates that are subjected to variable selection using sva
#' #' method
#' #' @param counts a numerical  \code{matrix} of raw counts
#' #' Counts must be non-negative and NAs are not permitted.
#' #' @param FixCov a \code{dataframe} of all covariates that
#' #' are always included in the model, either they are the main
#' #' factors of interest or the covariates that researchers
#' #' want to include in the model.
#' #' @param VarCov a \code{dataframe} of all covariates that
#' #' are subjected to  selection.
#' #' @param combine logical. If TRUE (default), combining all FixCov
#' #' and VarCov, if FALSE, use only FixCov
#'
#' #' @return a \code{dataframe} which is the combination of \code{AllCov}
#' #'and the \code{sva} covariates if there are any.
#' #' @examples
#' #' data(counts)
#' #' data(FixCov)
#' #' data(VarCov)
#' #' svacovout <-svacov(counts[1:100,], FixCov, VarCov)
#' #' dim(svacovout)
#' #' head(svacovout)
#' svacov <- function(counts, FixCov, VarCov, combine = FALSE){
#'   set.seed(1)
#'   if(combine){
#'     AllCov <- cbind(FixCov, VarCov)
#'   }else{
#'       AllCov <- FixCov
#'     }
#'   dm <- stats::model.matrix(stats::formula(paste0("~", paste0(names(AllCov), collapse = "+"))), data = AllCov)
#'   # dm0 <-stats::model.matrix(stats::formula(paste0("~", paste0(names(FixCov), collapse = "+"))), data = FixCov)
#'   lib.size = apply(counts, 2, stats::quantile, .75)
#'   y <- t(log2(t(counts + 0.5)/(lib.size + 1) * 1e+06))
#'   svaout <- sva::sva(dat = y, mod = dm)
#'   n.sv <- svaout$n.sv
#'   if(n.sv ==0){
#'     sv <- data.frame(dm[,0])
#'   }else{
#'     sv  <- data.frame(svaout$sv)
#'     colnames(sv) <- paste0("sva", 1:ncol(sv))
#'   }
#'   VarCov_sva <- cbind(AllCov, sv)
#'   VarCov_sva
#' }

#' Backward Selection Algorithm Using the P-values Ratio Relevance Measurement
#'
#' This function performs a backward selection algorithm on the set of
#' available covariates in RNA-seq analysis using the covariate
#' relevance measurement \code{r} defined as the ratio of the number
#' of p-values less than 0.05 to a fifth of the number of p-values greater than
#' .75. It starts with all covariates included in the model, then using \code{\link{limma}{voom}}
#' to fit the model, calculate $p$-values vector and the relevance measurement
#' \code{r} for each covariate. It will iteratively removes the covariate
#'  with smallest \code{r} value.  The process is repeated until
#'   all covariates are eliminated.
#'
#' @param counts a numerical  \code{matrix} of raw counts
#' Counts must be non-negative and NAs are not permitted.
#' @param FixCov a \code{data.frame} of the covariates that
#' are always included in the model.
#' @param VarCov a \code{data.frame} of the covariates subjected to  selection.
#' @param print.progress logical. If \code{TRUE} the print out the progress of
#' backward selection.
#'
#' @return A list of 4 components
#' \item{WorstP5}{a vector of \code{r} value of the eliminated covariates
#' during the selection process. The first element of the vector is the first
#' covariate eliminated, and so on.}
#' \item{VarCov}{it is the same as the input argument.}
#' \item{FixCov}{it is the same as the input argument.}
#' @examples
#' data(counts)
#' data(FixCov)
#' data(VarCov)
#' bsout <- bs(counts = counts, FixCov, VarCov, print.progress = TRUE)
#' names(bsout)
#' bsout$WorstP5
bs <- function(counts, FixCov, VarCov, print.progress = FALSE){
  if(!is.data.frame(FixCov)|!is.data.frame(VarCov)) stop("FixCov and VarCov have to be dataframe with different column names")
  DelCov <- VarCov[,0] # Initiaize Set of Deleted Covariates
  ConCov <- VarCov # Initialize Set of Considered Covariates
  PvList <- list() # List of P-values for the covariate selected at each iteration
  WorstP5 <- list() # Ratio of Number of p-value less than 0.05 and( >.075)/5 for the covariate selected at each iteration
  if(ncol(ConCov)!=0){
    Iter <- 1
    repeat{
      if(print.progress){
        cat("----------------------------------------\n")
        cat("Iteration = ", Iter, "\n")
        cat("DelCov = ", names(DelCov), "\n")
      }
      WorstP5[[Iter]] <- Inf
      AllCov <- cbind(FixCov, ConCov)
      dm <- stats::model.matrix(stats::formula(paste0("~",
                                                      paste0(names(AllCov),
                                                             collapse = "+"))),
                                data = AllCov)
      colnames(dm)[1] <- "Intercept"
      # if(!is.fullrank(dm)|ncol(dm)==nrow(dm)){ Ind <- T;
      # warning("Sigularity design matrix! Check generated pseudo-variables!"); break}
      vout <- limma::voom(counts = counts, design = dm, lib.size = apply(counts, 2, stats::quantile, .75), plot = F)
      suppressWarnings(fit <- limma::lmFit(vout))
      for(i in 1:ncol(ConCov)){
        if(is.factor(unlist(ConCov[i])) | is.character(unlist(ConCov[i]))) {
          ct <- paste0(grep(paste0(names(ConCov)[i]),  x = colnames(dm), value = T), collapse = ",  ")
        }else{
          ct <- paste0(grep(paste0(names(ConCov)[i], "$"),  x = colnames(dm), value = T), collapse = ",  ")
        }
        C.matrix <- eval(parse(text=paste0("limma::makeContrasts(",  ct, ",levels = dm)")))
        fit1 <- limma::contrasts.fit(fit, contrasts =C.matrix)
        fit1 <- limma::eBayes(fit1)
        tt <- limma::topTable(fit1, sort ="none", n = Inf)
        pv <- tt$P.Value
        P5 <- max(sum(pv<=.05, na.rm = T), 1)/max(sum(pv>.75, na.rm = T)/5, 1) # sometime pv has missing value NA make P5 cannot calculable
        if(print.progress)cat("ConCov = ", names(ConCov)[i], ", P5 = ", P5,  "\n")
        if(WorstP5[[Iter]] > P5) {
          WorstCov <- i
          WorstP5[[Iter]] <- P5
          PvList[[Iter]] <- pv
          names(PvList)[Iter] <- names(WorstP5)[Iter] <- names(ConCov)[i]
        }
      }
      Iter <- Iter +1
      DelCov <- cbind(DelCov, ConCov[WorstCov])
      ConCov <- ConCov[-WorstCov]
      if(ncol(ConCov)==0)break
    }
    WorstP5 <- do.call("c", WorstP5)
  }else{WorstP5 <- numeric(0)}
  res <- list(WorstP5 = WorstP5, VarCov = VarCov, FixCov = FixCov)
  res
}

#' Backward Selection Algorithm  in Nguyen et al. 2015, JABES
#'
#' This function implements the backward selection algorithm
#' in Nguyen et al. 2015, JABES. The covariate relevance measure
#'  is the number of p-values less than
#' 0.05. The algorithm is applied to the data until the main
#' factor of interest is eliminated from the model. Among the
#'  sets of remaining covariates, one that possesses the highest number
#'  of DDE genes (FDR <= 0.05) w.r.t. the main factor of interest
#'  is the final set of covariates selected by the algorithm.
#' @inheritParams bs
#' @return A list of 5 components
#'\item{WorstP5}{a vector of the number of p-values less than 0.05 for each
#'eliminated covariate during the selection process. The first element of the vector is the first
#' covariate eliminated, and so on.}
#' \item{Q5}{a vector of the number of q-values less than
#' 0.05 for the test of the effect of the main variable of interest
#' at each iteration.}
#' \item{VarCov}{it is the same as the input argument.}
#' \item{FixCov}{it is the same as the input argument.}
#' \item{BestCov}{ the set of selected covariates from
#' the list of sequences of subset of covariates remaining
#' at each iteration of the backward selection. }
#' @references Yet Nguyen, Dan Nettleton, Haibo Liu, Christopher K. Tuggle.
#' Detecting Differentially Expressed Genes with
#' RNA-seq Data Using Backward Selection to Account
#' for the Effects of Relevant Covariates.
#' J Agric Biol Environ Stat. 2015; 20(4): 577–597.
#' @author Yet Nguyen \email{tienyettoan@gmail.com}
#' @examples
#' data(counts)
#' data(FixCov)
#' data(VarCov)
#' print.progress <- FALSE
#' jabes.bsout <- jabes.bs(counts, FixCov, VarCov, print.progress)
#' names(jabes.bsout)
#' jabes.bsout$Q5
#' jabes.bsout$BestCov
jabes.bs <- function(counts, FixCov, VarCov,  print.progress = FALSE){
  # if(!(criteria %in% c("ksg", "pv"))) stp("criteria have to be either ksg or pv")
  if(!is.data.frame(FixCov)|!is.data.frame(VarCov)) stop("FixCov and VarCov have to be data frame with different column names")
  DelCov <- VarCov[,0] # Initiaize Set of Deleted Covariates
  ConCov <- VarCov # Initialize Set of Considered Covariates
  PvList <- list() # List of P-values for the covariate selected at each iteration
  WorstP5 <- list() # Number of p-value less than 0.05 for the covariate selected at each iteration
  Q5 <- list()
  if(ncol(ConCov)!=0){
    Iter <- 1
    repeat{
      if(print.progress){
        cat("----------------------------------------\n")
        cat("Iteration = ", Iter, "\n")
        cat("DelCov = ", names(DelCov), "\n")
      }
      WorstP5[[Iter]] <- Inf
      AllCov <- cbind(FixCov, ConCov)
      dm <- stats::model.matrix(stats::formula(paste0("~", paste0(names(AllCov), collapse = "+"))), data = AllCov)
      colnames(dm)[1] <- "Intercept"
      # if(!is.fullrank(dm)|ncol(dm)==nrow(dm)){ Ind <- T;
      # warning("Sigularity design matrix! Check generated pseudo-variables!"); break}
      vout <- limma::voom(counts = counts, design = dm, lib.size = apply(counts, 2, stats::quantile, .75), plot = F)
      suppressWarnings(fit <- limma::lmFit(vout))
      for(i in 1:ncol(AllCov)){
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
        P5 <- max(sum(pv<=.05, na.rm = T), 1) # sometime pv has missing value NA make P5 cannot calculable
        if(i ==1) Q5[[Iter]] <- sum(jabes.q(pv)<=.05)
        if(print.progress)cat("ConCov = ", names(AllCov)[i], ", P5 = ", P5,  "\n")
        if(WorstP5[[Iter]] > P5) {
          WorstCov <- i
          WorstP5[[Iter]] <- P5
          PvList[[Iter]] <- pv
          names(PvList)[Iter] <- names(WorstP5)[Iter] <- names(AllCov)[i]
        }
      }
      Iter <- Iter +1
      DelCov <- cbind(DelCov, AllCov[WorstCov])
      ConCov <- ConCov[-(WorstCov-1)]
      if(WorstCov==1)break
    }
    WorstP5 <- do.call("c", WorstP5)
    Q5 <- do.call("c", Q5)
  }else{WorstP5 <- numeric(0)}

  BestCov <- VarCov[setdiff(names(VarCov),names(WorstP5[1:(which(Q5==max(Q5))[1]-1)]))]
  res <- list(WorstP5 = WorstP5, Q5 = Q5,
              VarCov = VarCov, FixCov = FixCov,
              BestCov = BestCov)
  res
}



# Function Calculates KSG distance -----------------------------------------

#' Function Calculates KSG distance
#'
#' Description: This function calculates Komogorov-Smirnov Grenander Distance of
#' the empirical cummulative distribution function of the p-value vector and uniform unif(0, 1)
#' @param z a vector of p-values
#' @return Kolmogorov-Smirnov distance of Grenander estimator of p-values's ECDF
#' and uniform distribution u(0,1)
#'
#'
kss <- function(z){
  e <- ecdf(z)
  g <- fdrtool::grenander(e)
  if(tail(g$x.knots, 1)!=1){
    gx <- c(0, g$x.knots, 1); gF <- c(0, g$F.knots, 1)
  }
  if (tail(g$x.knots, 1)==1){
    gx <- c(0, g$x.knots); gF <- c(0, g$F.knots)
  }
  out <- max(sqrt(length(z))*abs(gF - gx))
  return(out)
}
