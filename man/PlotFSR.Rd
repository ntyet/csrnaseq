% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/3-FSRProportion.R
\name{PlotFSR}
\alias{PlotFSR}
\title{Plot false selection rate FDR estimates as a function of lambda}
\usage{
PlotFSR(FSROut, alpha0, option, B, m)
}
\arguments{
\item{FSROut}{the output of \link{FSR} function}

\item{alpha0}{the nominal FSR  level.}

\item{option}{the method to generate pseudo-variables
This is either one of these \code{"WN", "RX", "OWN", "ORX"}.}

\item{B}{a numerical value representing the number of
sets of the pseudo-variables. Usually, \code{B = 100} is
good enough for this method.}

\item{m}{number of generated pseudo-variables}
}
\value{
A plot of estimators \eqn{\alpha_{RE}}
and \eqn{\alpha_{ER}} as functions of \eqn{\lambda}.
}
\description{
This function uses \code{\link[ggplot2]{ggplot}} to make a plot of
the estimators \eqn{\alpha_{RE}}
and \eqn{\alpha_{ER}} as functions of \eqn{\lambda}.
}
\examples{
lambdamax <- 5
lambda <- c(seq(1, lambdamax, length = 400))
data(FixCov)
data(VarCov)
data(counts)
print.progress = FALSE
m <- 7
option <- "OWN"
B <- 20
ncores <- 10
svamethod <- FALSE
BBootBSOut<-BBootBS(B, ncores, counts[1:3000,],
FixCov, VarCov, m,  lambda, option, print.progress)
DataBSOut <- DataBS(counts[1:3000,], FixCov, VarCov,
lambda, print.progress, svamethod)
FSRPBSOut <- FSRPBS(DataBSOut, BBootBSOut)
alpha0 <- 0.05
BestLambdaOut <- BestLambdaBS(FSRPBSOut, alpha0)
FSROut <- FSR(FSRPBSOut, BestLambdaOut)
PlotFSR(FSROut, alpha0 = .05, option, B, m)
}
