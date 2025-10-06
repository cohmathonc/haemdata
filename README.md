# haemdata <img src='man/figures/logo.png' align="right" height="139" />
<!-- [![R-CMD-check](https://github.com/cohmathonchonchonc/haemdata/workflows/R-CMD-check/badge.svg)](https://github.com/cohmathonc/haemdata/actions)  -->

<!-- badges: start -->
![](https://img.shields.io/badge/code-unstable-red) <br>
[![R-CMD-check](https://github.com/cohmathonc/haemdata/workflows/R-CMD-check/badge.svg)](https://github.com/cohmathonc/haemdata/actions)
<!-- badges: end -->
An R package for processing and accessing genomics data from the PSON AML State-transition project at COH.

The package provides harmonised, preprocessed genomic data (mRNA, 10X, miRNA) from:
- mice collected between 2016 and 2023; 
- human samples from COH AML FLT3 patients; 
- diagnosis, remission and relapse AML patients from [Kim et al 2020](https://doi.org/10.1038/s41598-020-76933-2);
- MDS patients & cell lines from [Shiozawa et al 2017](http://doi.org/10.1182/blood-2017-05-783050).

Expression data and sample metadata are posted to the [Haemdatda Teams channel](https://cityofhope.sharepoint.com/:f:/r/sites/PSONAMLState-Transition/Shared%20Documents/haemdata?csf=1&web=1&e=Uh4VFb).

## Documentation
See [Getting started](https://cohmathonc.github.io/haemdata/articles/haemdata.html) for usage with examples, including exporting CSV tables for analysis with Matlab, for instance. See the [articles](articles) for more comprehensive analyses. 

## Installation
Installation is not *required* to use Haemdata, but installing the R package provides a number of useful helper functions that make working with these data more convenient within the R environment. 

You can install the default version of Haemdata from the package website (campus or VPN connection required):

``` r
install.packages("haemdata", repos = "http://cgt.coh.org/MHO")
```

Can't find something that was previously available? Check the [release history](https://github.com/cohmathonc/haemdata/releases).

### Windows users
Windows users require the `Rtools` toolchain bundle installed to build libraries, including some dependencies used by Haemdata. Find the appropriate installer for your version of R [here](https://cran.r-project.org/bin/windows/Rtools/).

## Publications
* Frankhouser et al 2022 *Science Advances* [![DOI](https://zenodo.org/badge/DOI/10.1126/sciadv.abj1664.svg)](https://doi.org/10.1126/sciadv.abj1664)
* Rockne et al 2020 *Cancer Research* [![DOI](https://zenodo.org/badge/DOI/10.1158/0008-5472.CAN-20-0354.svg)](https://doi.org/10.1158/0008-5472.CAN-20-0354)
