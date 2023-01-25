#' Get the counts table from mirtop output
#'
#' Given the path of a pipeline multiqc report file, this function will read in the counts from the mirtop output and return a matrix of counts.
#'
#' @title nfcore_mirtop_counts
#' @param multiqc_path the full path to a multiqc pipeline report on the Isilon storage
#' @return a matrix of counts, with genes as rows, and samples as columns
#' @author Denis O'Meally
#' @export
nfcore_mirtop_counts <- function(multiqc_path) {

    #Test if the multiqc is from a smrnaseq run
    if(!stringr::str_detect(multiqc_path, "nfcore-smrnaseq")) stop("The multiqc is not from a nfcore/smrnaseq run")

    mirtop_counts_tsv <- gsub("multiqc/multiqc_report.html", "mirtop/mirna.tsv", multiqc_path)

    mirtop_counts <- readr::read_delim(mirtop_counts_tsv, show_col_types = FALSE) |>
        rename_all(~ gsub("_seqcluster", "", .))

    return(mirtop_counts)

}
