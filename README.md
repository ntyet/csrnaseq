
<!-- README.md is generated from README.Rmd. Please edit that file -->

# Summary

<!-- badges: start -->
<!-- badges: end -->

The R package `csrnaseq` performs a backward variable selection
procedure using pseudo-variables for RNA-seq differential expression
analysis (Nguyen and Nettleton 2024+). The idea is to select the most
relevant covariates such that the false selection rate is below a
pre-specified threshold. The method is built upon the approach of Wu et
al. (2007). While Wu et al. (2007)’s method works for one response
variable, our method works for multiple response variables such as
RNA-seq data. The selected covariates are then included in differential
expression analysis using `voom-limma` pipeline (Law et al. 2014). The
proposed method is implemented in function `FSRAnalysisBS`.

## Installation

`csrnaseq` can be installed from
[GitHub](https://github.com/ntyet/csrnaseq):

``` r
# install.packages("devtools")
devtools::install_github("ntyet/csrnaseq")
```

## Example

This is a basic example that shows how to use our method:

``` r
library(csrnaseq)
data(counts)
data(FixCov)
data(VarCov)
option <- "OWN"
B <- 2
m <- 2
alphamax <- 5
alpha0 <- 0.05
ncores <- 1
print.progress <- FALSE
saveall <- TRUE
FSRAnalysisBSOut <- FSRAnalysisBS(counts, FixCov, VarCov, 
                                  option, B, m, alphamax, alpha0, 
                                  ncores, print.progress, saveall)
names(FSRAnalysisBSOut)
```

# Data Analysis and Simulation Based on the RFI RNA-seq Dataset

- Codes for RFI RNA-seq data analysis are
  [here](https://github.com/ntyet/csrnaseq/blob/main/analysis/1-FSRBS_RFIAnalysis.R)
  and
  [here](https://github.com/ntyet/csrnaseq/blob/main/analysis/1-FSRBS_RFIAnalysis.sh).

- Codes for generating six simulation scenarios are
  [here](https://github.com/ntyet/csrnaseq/blob/main/analysis/2-FSRBS_ModelSize.R)
  and
  [here](https://github.com/ntyet/csrnaseq/blob/main/analysis/2-FSRBS_ModelSize.sh).
  The outputs are
  [here](https://github.com/ntyet/csrnaseq/blob/main/analysis/RealDataOutBS).

- Codes for simulation are
  [here](https://github.com/ntyet/csrnaseq/blob/main/analysis/3-Simulation.R)
  and
  [here](https://github.com/ntyet/csrnaseq/blob/main/analysis/3-Simulation.sh).
  The outputs are
  [here](https://github.com/ntyet/csrnaseq/blob/main/analysis/SimulationOut).

- Codes for additional simulation to investigate covariates orthogonal
  to the primary variables that include:

  - Codes for generating six simulation scenarios:
    [here](https://github.com/ntyet/csrnaseq/blob/main/analysis/2-FSRBS_ModelSize-RFI_projection_on_Line.R).
    The outputs are
    [here](https://github.com/ntyet/csrnaseq/blob/main/analysis/revise/RealDataOutBSRFIprojectiononLine).
  - Codes for simulation:
    [here](https://github.com/ntyet/csrnaseq/blob/main/analysis/3-Simulation-RFI_projection_on_Line.R)
    and
    [here](https://github.com/ntyet/csrnaseq/blob/main/analysis/3-Simulation-RFI_projection_on_Line.sh).
    The outputs are
    [here](https://github.com/ntyet/csrnaseq/blob/main/analysis/revise/SimulationOutRFIprojectiononLine).

# Data Analysis and Simulation Based on the Zebrafish RNA-seq Dataset

- Zebrafish RNA-seq dataset are available at
  [here](https://github.com/ntyet/csrnaseq/tree/main/analysis/extra-rna-seq-data).

- Codes for Zebrafish RNA-seq data analysis are
  [here](https://github.com/ntyet/csrnaseq/blob/main/analysis/1-FSRBS_zebrafish.R)
  and
  [here](https://github.com/ntyet/csrnaseq/blob/main/analysis/1-FSRBS-zebrafish.sh).

- Codes for generating two simulation scenarios are
  [here](https://github.com/ntyet/csrnaseq/blob/main/analysis/2-FSRBS_ModelSize-zebrafish.R).
  The outputs are
  [here](https://github.com/ntyet/csrnaseq/blob/main/analysis/revise/RealDataOutBSzebrafish).

- Codes for simulation are
  [here](https://github.com/ntyet/csrnaseq/blob/main/analysis/3-Simulation-zebrafish.R)
  and
  [here](https://github.com/ntyet/csrnaseq/blob/main/analysis/3-Simulation-zebrafish.sh).
  The outputs are
  [here](https://github.com/ntyet/csrnaseq/blob/main/analysis/revise/SimulationOutzebrafish).

## References

<div id="refs" class="references csl-bib-body hanging-indent"
entry-spacing="0">

<div id="ref-law2014" class="csl-entry">

Law, C. W., Chen, Y., Shi, W., and Smyth, G. K. (2014), “Voom: Precision
weights unlock linear model analysis tools for
<span class="nocase">RNA-seq</span> read counts,” *Genome Biology*, 15,
R29. <https://doi.org/10.1186/gb-2014-15-2-r29>.

</div>

<div id="ref-nguyen2023" class="csl-entry">

Nguyen, Y., and Nettleton, D. (2024+), “Identifying relevant covariates
in <span class="nocase">RNA-seq</span> analysis by pseudo-variable
augmentation,” *Journal of Agricultural, Biological and Environmental
Statistics*, accepted.

</div>

<div id="ref-wu2007" class="csl-entry">

Wu, Y., Boos, D. D., and Stefanski, L. A. (2007), “Controlling variable
selection by the addition of pseudovariables,” *Journal of the American
Statistical Association*, Taylor & Francis, 102, 235–243.
<https://doi.org/10.1198/016214506000000843>.

</div>

</div>
