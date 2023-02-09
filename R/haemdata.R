#' @importFrom stats na.omit setNames
#' @importFrom utils read.csv read.table write.table
#' @importFrom magrittr %>%
#' @importFrom dplyr select case_when filter left_join mutate distinct
#' @import ggplot2
NULL

# Some global variables
globalVariables(c(
    ".", "DOB", "Gender", "Mice", "bam", "cellAmt", "ega_run_accession_id",
    "htbIdPb", "nFeature_RNA", "name", "number", "patient_id", "percent_mt",
    "percent_ribo", "value", "weeks",
    "width", "csv_pin", "csv_version", "Sample", "basepairs",
    "fastq_1", "fastq_2", "gene_id", "gene_name", "mouse_id",
    "sample_weeks", "strandedness", "cell_type", "cell_type_fine",
    "pruned.labels", "assay", "hdf5", "cohort", "library_id", "sample_id",
    "age_at_sample", "batch", "genotype", "percent_ckit",
    "sex", "tissue", "treatment", "filename", "dead", "published_metadata_mmu"
))
# setup package environment

#' @title Package environment
#' @details The `haemdata_env` package environment holds variables and objects used
#' by the package, for example the `pin_board`. See see the [package source](https://github.com/drejom/haemdata/blob/main/R/haemdata.R)
#' for more context.
#' @rdname haemdata_env
#' @export
haemdata_env <- new.env(parent = emptyenv())

# Message to display when no pinboard is set
haemdata_env$pin_board_msg <- "The pinboard is not set up. Please set the pinboard using the use_pinboard() function either in this R session or in the project \".Rprofile\" file.
See ?use_pinboard() for more information."

# pinboard is set to NULL for building/installing the package
haemdata_env$pin_board <- NULL

# Package URL
haemdata_env$package_url <- "http://cgt.coh.org/haemdata"


#### Utility functions for the target pipeline ----

# cross platform hostname
get_hostname <- function() {
    return(as.character(Sys.info()["nodename"]))
}

#' Write rda file to `data/`
#'
#' This function writes an `rda` file in `data`
#' and names it `schemeName`. Makes saving objects easier
#' from within a function. Inspired by [this](https://stackoverflow.com/questions/56293910/create-r-data-with-a-dynamic-variable-name-from-function-for-package)
#' StackOverflow post. Uses 'xz' to compress the file.
#'
#' @param schemeName the name of the rda file, without the `.rda` extension
#' @param data the object to save
#' @examples
#' \dontrun{
#' if (interactive()) {
#'     write_rda("an_rda_file_named_that", an_object_named_this)
#' }
#' }
#' @return the path to the rda file
#' @export
#' @rdname write_data
write_data <- function(schemeName, data) {
    rdaFile <- paste0(schemeName, ".rda")
    fileLocation <- file.path(".", "data", rdaFile)
    varName <- paste0(schemeName)

    assign(varName, data)
    eval(parse(text = sprintf("save(%s, file = '%s', compress = 'xz')", varName, fileLocation)))
    return(paste0("data/", rdaFile))
}
