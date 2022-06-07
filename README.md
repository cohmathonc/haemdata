
<!-- README.md is generated from README.Rmd. Please edit that file -->

# haemdata <img src='man/figures/logo.png' align="right" height="139" />
<!-- [![R-CMD-check](https://github.com/drejom/haemdata/workflows/R-CMD-check/badge.svg)](https://github.com/drejom/haemdata/actions)  -->

<!-- badges: start -->
![](https://img.shields.io/badge/code-unstable-red) <br>
<!-- badges: end -->
An R package for accessing genomics data from the PSON AML State-transition project at COH.

This package provides harmonised `SummarisedExperiments` for mice collected between 2016 and 2022, MDS patients from EGA and AML patients from the COH Biobank.

Briefly, raw sequence reads are processed using nf-core pipelines, sample metadata is collated from a range of supplied excel sheets, and outlying samples are flagged using [OUTRIDER](https://doi.org/10.1016/j.ajhg.2018.10.025). Package functions assist in manipulating and subsetting the SUmmarisedExperiment objects that should enable most analyses. If something is lacking, please [contribute](http://cgt.coh.org/haemdata/CONTRIBUTING.html) a solution or raise an [issue](https://github.com/drejom/haemdata/issues).

## Documentation

See the [package website](http://cgt.coh.org/haemdata) for usage with examples, including exporting CSV tables for matlab, for instance. 

## Installation

You can install the default version of haemdata from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("drejom/haemdata")
```

The latest release provides `SummarizedExperiments` produced with version 3.5 of the nf-core/rnaseq pipeline using the [GRCm38_HLT](articles/genomes.html) reference. 

## Publications

* Rockne et al 2020 *Cancer Research* <br>[![DOI](https://zenodo.org/badge/DOI/10.1158/0008-5472.CAN-20-0354.svg)](https://doi.org/10.1158/0008-5472.CAN-20-0354)

* Frankhouser et al 2022 *Science Advances* <br>[![DOI](https://zenodo.org/badge/DOI/10.1126/sciadv.abj1664.svg)](https://doi.org/10.1126/sciadv.abj1664)
