% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/3-FSRProportion.R
\name{BestLambdaBS}
\alias{BestLambdaBS}
\title{Calculate The Optimal Cut-off Threshold Lambda Star}
\usage{
BestLambdaBS(FSRPBSOut, alpha0 = 0.05)
}
\arguments{
\item{FSRPBSOut}{the output of \code{\link{FSRPBS}}.}

\item{alpha0}{the nominal FSR  level.}
}
\value{
a numerical vector consisting of
\item{lambdahat.RE}{ the optimal \code{lambda} when using \eqn{\alpha_{RE}}.}
\item{lambdahat.ER}{ the optimal \code{lambda} when using \eqn{\alpha_{ER}}.}
\item{kU.RE}{an estimate of the number of irrelevant covariates when using \eqn{\alpha_{RE}}.}
\item{kU.ER}{an estimate of the number of irrelevant covariates when using \eqn{\alpha_{ER}}.}
}
\description{
This function estimates the optimal threshold \code{lambda*}
and also the estimates of  the number of irrelevant covariates
using \eqn{\alpha_{RE}} and \eqn{\alpha_{ER}}.
}
\examples{
lambdamax <- 10
lambda <- c(seq(1, lambdamax, length = 800))
data(FixCov)
data(VarCov)
data(counts)
print.progress = FALSE
m <- 5
option <- "OWN"
B <- 2
ncores <- 1
svamethod <- FALSE
BBootBSOut<-BBootBS(B, ncores, counts[1:100,],
FixCov, VarCov, m,  lambda, option, print.progress)
DataBSOut <- DataBS(counts[1:100,], FixCov, VarCov,
lambda, print.progress, svamethod)
FSRPBSOut <- FSRPBS(DataBSOut, BBootBSOut)
alpha0 <- 0.05
BestLambdaOut <- BestLambdaBS(FSRPBSOut, alpha0)
BestLambdaOut
}
\author{
Yet Nguyen \email{tienyettoan@gmail.com}.
}
