# haemdata <img src='man/figures/logo.png' align="right" height="139" />
<!-- [![R-CMD-check](https://github.com/drejom/haemdata/workflows/R-CMD-check/badge.svg)](https://github.com/drejom/haemdata/actions)  -->

<!-- badges: start -->
![](https://img.shields.io/badge/code-unstable-red) <br>
<!-- badges: end -->
An R package for accessing genomics data from the PSON AML State-transition project at COH.

This package provides harmonised `SummarisedExperiments` for RNAseq of mice collected between 2016 and 2022, and from MDS & AML patients. `SummarisedExperiments` and sample metadata tables a are posted to [Haemdata Teams channel](https://cityofhope.sharepoint.com/:f:/r/sites/PSONAMLState-Transition/Shared%20Documents/haemdata?csf=1&web=1&e=Uh4VFb).

To build the package, raw sequence reads were processed using nf-core pipelines, poorly mapped samples were removed and outlying samples flagged using [OUTRIDER](https://doi.org/10.1016/j.ajhg.2018.10.025). Sample metadata was collated from a range of supplied excel sheets. Package [functions](http://cgt.coh.org/haemdata/reference/index.html) assist in visualising and sub-setting `SummarisedExperiment` objects and should facilitate most types of analyses. If something is incorrect or lacking, please [contribute](http://cgt.coh.org/haemdata/CONTRIBUTING.html) a solution or raise an [issue](https://github.com/drejom/haemdata/issues).

## Documentation
See [Getting started](http://cgt.coh.org/haemdata) for usage with examples, including exporting CSV tables for analyses with Matlab, for instance. 

## Installation

Installation is not *required* to use Haemdata, but installing the R package provides a number of useful helper functions that make working with these data more convenient within the R environment. Adding new data requires the package to be installed.

You can install the default version of Haemdata from the package website with:

``` r
install.packages("haemdata", repo = "http://cgt.coh.org/MHO")
```

The latest release provides `SummarizedExperiments` produced with version 3.5 of the nf-core/rnaseq pipeline using the [GRCm38_HLT](articles/genomes.html) reference. 

Can't find something that was previously available? Check the [release history](https://github.com/drejom/haemdata/releases).

## Publications
* Frankhouser et al 2022 *Science Advances* <br>[![DOI](https://zenodo.org/badge/DOI/10.1126/sciadv.abj1664.svg)](https://doi.org/10.1126/sciadv.abj1664)
* Rockne et al 2020 *Cancer Research* <br>[![DOI](https://zenodo.org/badge/DOI/10.1158/0008-5472.CAN-20-0354.svg)](https://doi.org/10.1158/0008-5472.CAN-20-0354)
