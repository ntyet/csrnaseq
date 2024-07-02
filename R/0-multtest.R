#' Estimate Number of True Null Hypotheses Using Histogram-based Method
#' (Nettleton 2006)
#'
#' This function estimates the number of true null hypotheses given a vector of
#' p-values using the method of Nettleton et al. (2006) JABES 11, 337-356. The
#' estimate obtained is identical to the estimate obtained by the iterative
#' procedure described by Mosig et al. Genetics 157:1683-1698. The number of
#' p-values falling into B equally sized bins are counted. The count of each bin
#' is compared to the average of all the bin counts associated with the current
#' bins and all bins to its right.  Working from left to right, the first bin
#' that has a count less than or equal to the average is identified. That
#' average is multiplied by the total number of bins to obtain an estimate of
#' m0, the number of tests for which the null hypothesis is true.
#' @param p a numerical vector of p-value
#' @param B number of bin
#' @return The function returns an estimate of m0, the number of tests for which
#' the null hypothesis is true.
#' @author Dan Nettleton \email{dnett@iastate.edu}
#' @references Dan NETTLETON, J. T. Gene HWANG, Rico A. CALDO, and Roger P. WISE.
#' Estimating the Number of True Null Hypotheses From a Histogram of p Values.
#' Journal of Agricultural, Biological, and Environmental Statistics,
#' Volume 11, Number 3, Pages 337–356.
#' @examples
#' p <- runif(100)
#' m0 <- estimate.m0(p)
#' m0
estimate.m0 = function(p, B = 20) {
    m <- length(p)
    m0 <- m
    bin <- c(-0.1, (1:B)/B)
    bin.counts = rep(0, B)
    for (i in 1:B) {
        bin.counts[i] = sum((p > bin[i]) & (p <= bin[i + 1]))
    }
    tail.means <- rev(cumsum(rev(bin.counts))/(1:B))
    temp <- bin.counts - tail.means
    index <- min((1:B)[temp <= 0])
    m0 <- B * tail.means[index]
    return(m0)
}

#' Estimating The threshold lambda in Estimating FDR Using Storey's Formula
#'
#' This function calculate lambda threshold to estimate FDR
#'
#' @inheritParams estimate.m0
#' @return The function returns an estimate of lambdastar, the threshold
#' to calculate FDR based on Storey's formula.
#' @author Dan Nettleton \email{dnett@iastate.edu}
#' @references Dan NETTLETON, J. T. Gene HWANG, Rico A. CALDO, and Roger P. WISE.
#' Estimating the Number of True Null Hypotheses From a Histogram of p Values.
#' Journal of Agricultural, Biological, and Environmental Statistics,
#' Volume 11, Number 3, Pages 337–356.
#' @examples
#' p <- runif(100)
#' lambdastar <- estimate.lambda(p)
#' lambdastar
estimate.lambda = function(p, B = 20) {
  m <- length(p)
  m0 <- m
  bin <- c(-0.1, (1:B)/B)
  bin.counts = rep(0, B)
  for (i in 1:B) {
    bin.counts[i] = sum((p > bin[i]) & (p <= bin[i + 1]))
  }
  tail.means <- rev(cumsum(rev(bin.counts))/(1:B))
  temp <- bin.counts - tail.means
  index <- min((1:B)[temp <= 0])
  lambdastar <- index/B
  return(lambdastar)
}
#' Q-value Using Histogram-based Method (Nettleton 2006)
#'
#' This function computes q-values using the approach of Nettleton et al.
#'(2006) JABES 11, 337-356.  Author: Dan Nettleton
#'@return The function returns a q-value vector of the input p-value vector.
#' @inheritParams estimate.m0
#' @author Dan Nettleton \email{dnett@iastate.edu}
#' @references Dan NETTLETON, J. T. Gene HWANG, Rico A. CALDO, and Roger P. WISE.
#' Estimating the Number of True Null Hypotheses From a Histogram of p Values.
#' Journal of Agricultural, Biological, and Environmental Statistics,
#' Volume 11, Number 3, Pages 337–356.
#' @examples
#' p <- runif(100)
#' q <- jabes.q(p)
#' sum(q <= .05)
jabes.q = function(p, B = 20) {

    m = length(p)
    m0 = estimate.m0(p, B)
    k = 1:m
    ord = order(p)
    p[ord] = (p[ord] * m0)/(1:m)
    qval = p
    qval[ord] = rev(cummin(rev(qval[ord])))
    return(qval)
}

