# Functions for amending metadata_mmu with QC information

#' Amend metadata_mmu with qc_pass_mapping column
#'
#' The function returns the `metadata_mmu` table amended with a
#' `qc_pass_mapping` column, where `TRUE` indicates the per cent of uniquely
#' mapped reads exceeds the `mapping_threshold`, or `FALSE` if not.
#'
#' @title flag_lowly_mapped_mmu_se
#' @param summarised_experiment a SummarisedExperiment from a "QC" run (`qc = TRUE`
#' set in [run_nf_core_rnaseq()]
#' @param mapping_threshold an integer. Threshold below which the `qc_pass_mapping`
#' will be `FALSE`. Default: 5 per cent.
#' @return a data.frame
#' @author Denis O'Meally
#' @export
flag_lowly_mapped_mmu_se <- function(summarised_experiment, mapping_threshold = 5) {

    # Check that the summarised_experiment is from a mouse run
    if (stringr::str_detect(S4Vectors::metadata(summarised_experiment)["object_name"], "mmu", negate = TRUE)) {
        stop("The summarised_experiment is not from a mmu run.")
    }

    # Check that the summarised_experiment is from a qc run
    if (!S4Vectors::metadata(summarised_experiment)["workflow"] == "qc") {
        stop("The summarised_experiment is not from a qc run and so no count of uniquely mapping reads is available.")
    }

    qc_pass_mapping <- SummarizedExperiment::colData(summarised_experiment) |>
        dplyr::as_tibble() |>
        dplyr::select(sample, star.uniquely_mapped_percent) |>
        dplyr::mutate(qc_pass_mapping = ifelse(star.uniquely_mapped_percent >= mapping_threshold, TRUE, FALSE)) |>
        dplyr::select(sample, qc_pass_mapping)

    ## Update metadata with qc_pass_mapping
    use_pinboard("onedrive")
    metadata_mmu <- get_pin("metadata_mmu.csv") |>
        dplyr::left_join(qc_pass_mapping, by = "sample")

    return(metadata_mmu)
}
