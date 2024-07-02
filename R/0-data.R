#' The set of available covariates
#'
#' @format \code{covset} is a data frame with 31 rows and 14 columns
#' representing 14 covariates of 31 pig, including
#' \describe{
#' \item{Line}{is the categorical factor of primary scientific interest.
#' \code{Line} has two levels, high RFI (HRFI)  and low RFI (LRFI).
#' 15 pigs were from LRFI line and 16 were from  HRFI line.}
#' \item{RFI}{RFI is a continuous covariate that provides a measure
#' of the residual feed intake for each of the 31 pigs from which
#' blood samples were drawn for RNA-seq analysis. Pigs in the HRFI
#' line tend to have high RFI values, while pigs in the LRFI line
#' tend to have low RFI values.}
#' \item{Diet}{is a categorical factor with two levels corresponding to
#' the two diets (high fiber, low energy vs. low fiber, high energy)
#' that were fed to the pigs in this study. Approximately half the
#' pigs within each line were fed each diet. Because RNA-seq analysis
#' was performed on blood samples collected prior to the initiation
#' of the two diets, this factor is not expected to be associated with
#' the transcript abundance levels measured by RNA-seq. However,
#' other covariates in this study may be associated with Diet.
#' For example, feed consumption on the assigned diet played
#' a role in the calculation of RFI values, and pigs on the high fiber,
#' low energy diet tended to have higher RFI values than pigs on
#' the low fiber, high energy diet.}
#' \item{Baso}{is a continuous covariate that provides a measure of
#' the concentration of basophil cells in the blood sample drawn
#' from each pig.}
#' \item{Eosi}{is a continuous covariate that provides a measure
#' of the concentration of eosinophil cells in the blood sample
#' drawn from each pig.}
#' \item{Lymp}{is a continuous covariate that provides a measure of
#' the concentration of lymphocyte cells in the blood sample
#' drawn from each pig.}
#' \item{Mono}{is a continuous covariate that provides a measure
#' of the concentration of monocyte cells in the blood sample
#' drawn from each pig.}
#' \item{Neut}{is a continuous covariate that provides a measure
#' of the concentration of neutrophil cells in the blood sample
#' drawn from each pig.}
#' \item{Concb}{is a continuous measure of the RNA concentration
#' in each sample before globin depletion (a step that is
#' necessary to focus sequencing efforts on messenger RNA
#' molecules other than highly abundant globin messenger
#' RNA in each blood sample).}
#' \item{Conca}{is a continuous measure of the RNA concentration
#' in each sample after globin depletion.}
#'\item{RINb}{is a continuous measure of RNA integrity within
#'each sample before globin depletion.}
#'\item{RINa}{is a continuous measure of RNA integrity within
#'each sample after globin depletion.}
#'\item{Block}{is a categorical factor with four levels
#'corresponding to the four blocks used to organize sample
#'collection and processing. Initially, each block involved
#'eight samples, two for each combination of Line and Diet.
#'One LRFI sample from the first block was removed from the
#'study due to low-quality RNA.}
#'\item{Order}{is a categorical factor with eight levels
#'indicating the random order samples were processed
#'within each block.}
#' }
#'
#' @examples
#' data(covset)
#' colnames(covset)
#' head(covset)

"covset"


#' RFI RNA-seq Data
#'
#' This is the RFI RNA-seq count data used  as an example for this study.
#'
#' @format The dataset has 12280 rows and 31 columns corresponding
#' to 31 RNA-seq samples from 31 pigs.
#'
#' @examples
#' data(counts)
#' dim(counts)
#' head(counts)
"counts"

#' Dataframe Containing The Set of Fix Covariates
#'
#'This data set contains the covariates that are not subjected to variable
#'selection. In this case, it is \code{Line}, the main factor of
#'interest.
#'
#'@format \code{FixCov} is a data frame containing only \code{Line}
#' @examples
#' data(FixCov)
#' names(FixCov)
#' FixCov

"FixCov"

#' A Data frame of Covariates
#'
#'This data set contains the covariates subjected to variable
#'selection. In this case, it includes \code{RFI, Diet, Neut, Lymp, Mono,
#'Neut, Baso, Eosi, Conca, Concb, RINa, RINb, Block, Order.}
#'
#'@format \code{VarCov} is a data frame containing 13 covariates.
#' @examples
#' data(VarCov)
#' names(VarCov)
#' VarCov

"VarCov"

