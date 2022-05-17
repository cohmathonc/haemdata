
<!-- README.md is generated from README.Rmd. Please edit that file -->

# haemdata <img src='man/figures/logo.png' align="right" height="139" />
<!-- [![R-CMD-check](https://github.com/drejom/haemdata/workflows/R-CMD-check/badge.svg)](https://github.com/drejom/haemdata/actions)  -->

<!-- badges: start -->
![](https://img.shields.io/badge/code-unstable-red) <br>
<!-- badges: end -->
A data package for RNAseq of AML and CML mouse models generated at COH.

This package provides harmonised expression matrices and limited metadata for mice collected between 2016 and 2022.

## Documentation

See the [package website](http://cgt.coh.org/haemdata) for usage with examples, including exporting CSV tables for matlab for instance. 

## Installation

You can install the default version of haemdata from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("drejom/haemdata")
```

The latest release provides `SummarizedExperiments` produced with version 3.5 of the nf-core/rnaseq pipeline using the GRCm38_HLT reference. 

## Publications

* Rockne et al 2020 *Cancer Research* <br>[![DOI](https://zenodo.org/badge/DOI/10.1158/0008-5472.CAN-20-0354.svg)](https://doi.org/10.1158/0008-5472.CAN-20-0354)

* Frankhouser et al 2022 *Science Advances* <br>[![DOI](https://zenodo.org/badge/DOI/10.1126/sciadv.abj1664.svg)](https://doi.org/10.1126/sciadv.abj1664)
