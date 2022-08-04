# haemdata <img src='man/figures/logo.png' align="right" height="139" />
<!-- [![R-CMD-check](https://github.com/drejom/haemdata/workflows/R-CMD-check/badge.svg)](https://github.com/drejom/haemdata/actions)  -->

<!-- badges: start -->
![](https://img.shields.io/badge/code-unstable-red) <br>
[![R-CMD-check](https://github.com/drejom/haemdata/workflows/R-CMD-check/badge.svg)](https://github.com/drejom/haemdata/actions)
<!-- badges: end -->
An R package for processing and accessing genomics data from the PSON AML State-transition project at COH.

The package provides harmonised, preprocessed genomic data (mRNA, 10X, miRNA) from mice collected between 2016 and 2022, COH AML FLT3 patients, and MDS patients & cell lines from [Shiozawa et al 2017](http://doi.org/10.1182/blood-2017-05-783050). Expression data and sample metadata are posted to the [Haemdatda Teams channel](https://cityofhope.sharepoint.com/:f:/r/sites/PSONAMLState-Transition/Shared%20Documents/haemdata?csf=1&web=1&e=Uh4VFb).

To build the `SummarizedExperiments` and expression matrices, raw sequence reads were processed using nf-core pipelines, poorly mapped samples were removed and outlying samples flagged using [`OUTRIDER`](https://doi.org/10.1016/j.ajhg.2018.10.025). Single-cell RNAseq was preprocessed with 10X Cellranger and annotated `Seurat` and `Scanpy` objects produced with Seurat. Sample metadata was collated from a range of supplied excel sheets. Package [functions](http://cgt.coh.org/haemdata/reference/index.html) assist in visualising and sub-setting `SummarizedExperiment` and `Seurat` objects and should facilitate most types of analyses. If something is incorrect or lacking, please [contribute](http://cgt.coh.org/haemdata/CONTRIBUTING.html) a solution or raise an [issue](https://github.com/drejom/haemdata/issues).

## Documentation
See [Getting started](http://cgt.coh.org/haemdata) for usage with examples, including exporting CSV tables for analysis with Matlab, for instance. 

## Installation

Installation is not *required* to use Haemdata, but installing the R package provides a number of useful helper functions that make working with these data more convenient within the R environment. Adding new data requires the package to be installed.

You can install the default version of Haemdata from the package website with:

``` r
install.packages("haemdata", repos = "http://cgt.coh.org/MHO", dependencies = TRUE)
```

The latest release uses v3.7 of the nf-core/rnaseq pipeline, Cellranger v6.1.1 & `Seurat` v4.1.1 and `sctranscform` v0.3.3.

Can't find something that was previously available? Check the [release history](https://github.com/drejom/haemdata/releases).

### Windows users

Windows users require the `Rtools` toolchain bundle installed to build libraries, including some dependencies used by Haemdata. Find the appropriate installer for your version of R [here](https://cran.r-project.org/bin/windows/Rtools/).

## Publications
* Frankhouser et al 2022 *Science Advances* <br>[![DOI](https://zenodo.org/badge/DOI/10.1126/sciadv.abj1664.svg)](https://doi.org/10.1126/sciadv.abj1664)
* Rockne et al 2020 *Cancer Research* <br>[![DOI](https://zenodo.org/badge/DOI/10.1158/0008-5472.CAN-20-0354.svg)](https://doi.org/10.1158/0008-5472.CAN-20-0354)
