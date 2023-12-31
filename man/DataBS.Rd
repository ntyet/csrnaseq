% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/2-DataBS.R
\name{DataBS}
\alias{DataBS}
\title{Backward Selection with Critical Thresholds}
\usage{
DataBS(counts, FixCov, VarCov, lambda, print.progress = FALSE)
}
\arguments{
\item{counts}{a numerical  \code{matrix} of raw counts
Counts must be non-negative and NAs are not permitted.}

\item{FixCov}{a \code{data.frame} of the covariates that
are always included in the model.}

\item{VarCov}{a \code{data.frame} of the covariates subjected to  selection.}

\item{lambda}{a vector of critical thresholds for the relevance measure.
The typical \code{lambda} has range from 1 to 10.}

\item{print.progress}{logical. If \code{TRUE} the print out the progress of
backward selection.}
}
\value{
A list of 6 components
\item{BS}{a vector of \code{r} value of the eliminated covariates
during the selection process.}
\item{S}{A vector of the number of selected covariates w.r.t  the
critical threshold vector \code{lambda}, i.e, S[i] is the largest number
of selected covariates whose their r values > lambda[i].}
\item{lambda}{ a numerical vector of the critical thresholds.}
\item{FixCov}{it is the same as the input argument.}
\item{VarCov}{it is the same as the input argument.}
}
\description{
This function performs backward selection algorithm.
The differences between \code{DataBS} and \code{bs} are that the
former returns also the number of selected covariates \code{S} with respect
to a critical threshold\code{lambda}, which are described as in
the \code{Return} section below.
}
\examples{
lambdamax <- 10
lambda <- c(seq(1, lambdamax, length = 400))
data(FixCov)
data(VarCov)
data(counts)
print.progress = FALSE
DataBSOut <- DataBS(counts, FixCov, VarCov, lambda, print.progress)
names(DataBSOut)
}
