#' haemdata: PSON omics data
#'
#' Harmonised data from mouse RNAseq projects processed with [nf-core](https://nf-co.re) pipelines
#' and published as the R package `haemdata` on GitHub.
#' @importFrom stats na.omit setNames
#' @importFrom utils read.csv read.table write.table
#' @importFrom magrittr %>%
#' @importFrom dplyr select case_when filter left_join mutate distinct 
#' @importFrom ggplot2 theme geom_vline aes 

#' @docType package
#' @name haemdata
NULL

# Some global variables
globalVariables(c(
    ".", "DOB", "Gender", "Mice", "bam", "cellAmt", "ega_run_accession_id",
    "htbIdPb", "nFeature_RNA", "name", "number", "patient_id", "percent_mt",
    "percent_ribo", "teams", "timepoint_project", "value", "weeks",
    "width"
))
