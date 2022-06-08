
# haemdata <img src='man/figures/logo.png' align="right" height="139" />
<!-- [![R-CMD-check](https://github.com/drejom/haemdata/workflows/R-CMD-check/badge.svg)](https://github.com/drejom/haemdata/actions)  -->

<!-- badges: start -->
![](https://img.shields.io/badge/code-unstable-red) <br>
<!-- badges: end -->
An R package for accessing genomics data from the PSON AML State-transition project at COH.

This package provides harmonised `SummarisedExperiments` for mice collected between 2016 and 2022, and R code to produce the same from MDS & AML patients. `SummarisedExperiments` for human data are posted to the [package website](http://cgt.coh.org/haemdata), [Teams channel](https://teams.microsoft.com/l/channel/19%3a210be89215cc4b2c878442a07b1580db%40thread.tacv2/haemdata?groupId=22521432-ac7e-43f8-be63-eb9f86a6f561&tenantId=972a3ea3-f979-4875-a2b8-cff001ab69e7) and to `MHO/haemdata` but cannot be installed with the package.

To build the package, raw sequence reads were processed using nf-core pipelines, outlying samples flagged using [OUTRIDER](https://doi.org/10.1016/j.ajhg.2018.10.025), and sample metadata collated from a range of supplied excel sheets. Package [functions](http://cgt.coh.org/haemdata/reference/index.html) assist in visualising and sub-setting `SummarisedExperiment` objects and should facilitate most types of analyses. If something is incorrect or lacking, please [contribute](http://cgt.coh.org/haemdata/CONTRIBUTING.html) a solution or raise an [issue](https://github.com/drejom/haemdata/issues).
## Documentation
See [Getting started](http://cgt.coh.org/haemdata) for usage with examples, including exporting CSV tables for analyses with Matlab, for instance. 

## Installation
You can install the default version of haemdata from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("drejom/haemdata")
```

The latest release provides `SummarizedExperiments` produced with version 3.5 of the nf-core/rnaseq pipeline using the [GRCm38_HLT](articles/genomes.html) reference. 

Can't find previously available? Check the [release history](https://github.com/drejom/haemdata/releases).
## Publications
* Frankhouser et al 2022 *Science Advances* <br>[![DOI](https://zenodo.org/badge/DOI/10.1126/sciadv.abj1664.svg)](https://doi.org/10.1126/sciadv.abj1664)
* Rockne et al 2020 *Cancer Research* <br>[![DOI](https://zenodo.org/badge/DOI/10.1158/0008-5472.CAN-20-0354.svg)](https://doi.org/10.1158/0008-5472.CAN-20-0354)
