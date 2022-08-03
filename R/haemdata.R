#' @importFrom stats na.omit setNames
#' @importFrom utils read.csv read.table write.table
#' @importFrom magrittr %>%
#' @importFrom dplyr select case_when filter left_join mutate distinct
#' @importFrom ggplot2 theme geom_vline aes
NULL

# Some global variables
globalVariables(c(
    ".", "DOB", "Gender", "Mice", "bam", "cellAmt", "ega_run_accession_id",
    "htbIdPb", "nFeature_RNA", "name", "number", "patient_id", "percent_mt",
    "percent_ribo", "timepoint_project", "value", "weeks",
    "width", "csv_pin", "csv_version","Sample", "basepairs",
    "fastq_1", "fastq_2", "gene_id", "gene_name", "mouse_id",
    "strandedness"
))

# setup package environment
#' @export
haemdata_env <- new.env(parent = emptyenv())

# Message to display when no pinboard is set
haemdata_env$pin_board_msg <- "The pinboard is not set up. Please set the pinboard using the use_pinboard() function.
See ?use_pinboard() for more information."

# pinboard is set to NULL for building/installing the package
haemdata_env$pin_board <- NULL

# Package URL
haemdata_env$package_url <- "http://cgt.coh.org/haemdata"
