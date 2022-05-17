#' Get a `SummarisedExperiment` object from an nf-core/rnaseq run
#'
#' Pulls in a `SummarisedExperiment` from the nf-core/rnaseq Salmon folder, and annotates
#' experiment metadata derived from the path provided.
#'
#' @title get_rnaseq_se
#' @param isilon_path a path on the COH Isilon store for raw data. Should begin with `/net`.
#' @return a `SummarisedExperiment` object
#' @author Denis O'Meally
#' @export
get_rnaseq_se <- function(isilon_path) {

    summarised_experiment <- readRDS(paste0(isilon_path, "/salmon/salmon.merged.gene_counts_length_scaled.rds"))

    project <- isilon_path |>
        stringr::str_extract("(?<=MHO\\/)(.*?)(?=\\/)") |>
        stringr::str_replace_all("[.]", "_")

    reference_genome <- ifelse(grepl("HLT", isilon_path),
            dplyr::case_when(
                grepl("GRCm38", isilon_path) ~ "GRCm38_HLT",
                grepl("GENCODEv33", isilon_path) ~ "GENCODEv33_HLT"),
            stop("Reference genome doesn't include human leukemic transgenes (HLT)")
        )

    pipeline <- isilon_path |>
        stringr::str_extract("(?<=\\/results\\/)(.*?)(?=[\\._]G)")

    object_name <- paste0(tolower(project), "_", reference_genome)

    S4Vectors::metadata(summarised_experiment)$object_name <- object_name
    S4Vectors::metadata(summarised_experiment)$project <- project
    S4Vectors::metadata(summarised_experiment)$reference_genome <- reference_genome
    S4Vectors::metadata(summarised_experiment)$pipeline <- pipeline

    # Stash a copy of the raw data in the repo
    write_data(object_name, summarised_experiment)

    return(summarised_experiment)
}
