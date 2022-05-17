#' @title write_tpm_matrix
#' @description Write an expression matrix to a CSV file in `extdata`
#' @param summarised_experiment a `SummarisedExperiment` object
#' @param drop samples (columns) to remove
#' @param tpm minimum TPM for a gene to be kept, Default: 1
#' @param samples minimum number of samples a gene must be present in to be kept, Default: 5
#' @return path to the CSV file
#' @details This function writes an expression matrix to a CSV file in `extdata`,
#' filtering out genes with fewer than `tpm` reads in at least `samples` samples.
#' @examples
#' \dontrun{
#' if (interactive()) {
#'     write_tpm_matrix(data(cml_mrna_GRCm38_HLT), c("F545.7", "X490.17"), tpm = 1, samples = 5)
#' }
#' }
#' @seealso
#'  \code{\link[SummarizedExperiment]{SummarizedExperiment-class}}
#'  \code{\link[utils]{write.table}}
#' @rdname write_tpm_matrix
#' @export
#' @importFrom SummarizedExperiment assay
#' @importFrom utils write.csv
write_tpm_matrix <- function(summarised_experiment, drop=NULL, tpm=1, samples=5) {
    out_name <- paste(tolower(summarised_experiment@metadata$project), summarised_experiment@metadata$reference_genome, sep = "_")
    # Subset to samples to drop
    summarised_experiment <- summarised_experiment[, !(summarised_experiment$names %in% drop)]
    # Make a CSV of TPMs, keeping genes with > 1 TPM in 5 samples, and transgenes
    mat <- SummarizedExperiment::assay(summarised_experiment, "abundance")
    filter <- rowSums(mat >= tpm) >= samples
    filter[grep("^HSA_", rownames(summarised_experiment))] <- TRUE
    filtered <- summarised_experiment[filter, ]
    tpm_matrix <- SummarizedExperiment::assay(filtered, "abundance")
    file_name <- paste0("inst/extdata/", out_name, "_", tpm, "tpm_in_", samples, "samples.csv")
    utils::write.csv(tpm_matrix, file_name)
    return(file_name)
}




#' Write rda file to `data/`
#'
#' This function writes an `rda` file in `data`
#' and names it `schemeName`. Makes saving objects easier
#' from within a function. Inspired by [this](https://stackoverflow.com/questions/56293910/create-r-data-with-a-dynamic-variable-name-from-function-for-package)
#' StackOverflow post.
#'
#' @param schemeName the name of the rda file
#' @param data the object to save
#'
#https://stackoverflow.com/questions/56293910/create-r-data-with-a-dynamic-variable-name-from-function-for-package
write_data <- function(schemeName, data) {
    rdaFile <- paste0(schemeName, ".rda")
    fileLocation <- file.path(".", "data", rdaFile)
    varName <- paste0(schemeName)

    assign(varName, data)
    eval(parse(text = sprintf("save(%s, file = '%s', compress = 'xz')", varName, fileLocation)))
    return(paste0("data/", rdaFile))
}
