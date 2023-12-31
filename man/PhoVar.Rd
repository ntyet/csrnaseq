% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/2-DataBS.R
\name{PhoVar}
\alias{PhoVar}
\title{Pseudo Variable Generation}
\usage{
PhoVar(dm0, m, option)
}
\arguments{
\item{dm0}{full model design matrix, the first column is the \code{intercept}}

\item{m}{number of generated pseudo-variables}

\item{option}{the method to generate pseudo-variables
This is either one of these \code{"WN", "RX", "OWN", "ORX"}.}
}
\value{
a data frame of \code{n} generated pseudo-variables.
}
\description{
This function generates pseudo-variables that are used in the
False Selection Rate variable selection method.
}
\examples{
data(FixCov)
data(VarCov)
AllCov <- cbind(FixCov, VarCov)
dm0 <- model.matrix(~., data = AllCov)
m <- 5
option <- "OWN"
PhoVarOut <- PhoVar(dm0, m, option)
head(PhoVarOut)
}
