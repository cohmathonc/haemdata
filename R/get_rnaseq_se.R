#' Get a `SummarisedExperiment` object from an nf-core/rnaseq run
#'
#' Pulls in a `SummarisedExperiment` from the nf-core/rnaseq Salmon folder, and annotates
#' experiment metadata derived from the multiqc path provided. Used in conjunction with [run_nf_core_rnaseq()].
#'
#' @name get_rnaseq_se
#' @param multiqc_path path to the multiqc report of a successful nf-core/rnaseq pipeline run
#' @param gtf annotate the SummarisedExperiment using the `gtf` recorded in the multiqc report, Default: TRUE
#' @return a `SummarisedExperiment` object
#' @author Denis O'Meally
#' @seealso
#' \code{\link{run_nf_core_rnaseq}}
#' @rdname get_rnaseq_se
#' #' @examples
#' \dontrun{
#' if(interactive()) {
#'  get_rnaseq_se("/net/isi-dcnl/ifs/user_data/rrockne/MHO/haemdata-nf-core-cache/AML.mRNA.2016/nfcore-rnaseq-v3.7_GENCODEm28_HLT")
#'  }
#' }
#' @export
get_rnaseq_se <- function(multiqc_path, gtf = TRUE) {

    star_salmon <- ifelse(stringr::str_detect(multiqc_path, "star_salmon"), TRUE, FALSE)

    out_path <- ifelse(star_salmon,
        gsub("/multiqc/star_salmon/multiqc_report.html", "", multiqc_path),
        gsub("/multiqc/multiqc_report.html", "", multiqc_path)
    )

    # load SummarisedExperiment
    salmon_se <- readRDS(glue::glue("{out_path}/salmon/salmon.merged.gene_counts.rds"))

    # If gtf = true, extract the gtf file from the multiqc report and add gene annotations to rowData
    if (gtf) {
        gtf <- readr::read_file(multiqc_path) |>
            stringr::str_extract("\\/ref_genome\\/igenomes\\/.*\\.gtf")

        type <- ifelse(stringr::str_detect(gtf, "GENCODE"), "gene_type", "gene_biotype")

        row_data <- rtracklayer::import(gtf) |>
            tibble::as_tibble() |>
            dplyr::select(gene_id, gene_name, basepairs = width, gene_type = !!sym(type)) |>
            dplyr::group_by(gene_id) |>
            dplyr::slice(which.max(basepairs)) |>
            dplyr::ungroup()

        SummarizedExperiment::rowData(salmon_se) <- dplyr::left_join(
            data.frame(gene_id = rownames(salmon_se)),
            row_data,
            by = "gene_id"
        )
    }

    # Add some experiment and analysis metadata
    run_folder <- multiqc_path |>
        stringr::str_extract("(?<=haemdata-nf-core-cache\\/)(.*?)(?=\\/)") |>
        stringr::str_replace_all("[.]", "_")

    reference_genome <- dplyr::case_when(
        stringr::str_detect(multiqc_path, "GRCm38_HLT") ~ "GRCm38_HLT",
        stringr::str_detect(multiqc_path, "GENCODEm28_HLT") ~ "GENCODEm28_HLT",
        stringr::str_detect(multiqc_path, "GENCODEr40") ~ "GENCODEr40")

    rnaseq_release <- multiqc_path |>
        stringr::str_extract("(?<=\\/nfcore-rnaseq-)(.*?)(?=[_]G)")

    workflow <- ifelse(star_salmon, "qc", "salmon")

    object_name <- glue::glue("{run_folder}_{reference_genome}_{workflow}_se")

    S4Vectors::metadata(salmon_se) <- data.frame(
        object_name = object_name,
        run_folder = run_folder,
        reference_genome = reference_genome,
        rnaseq_release = rnaseq_release,
        workflow = workflow,
        multiqc_url = sub("/net/isi-dcnl/ifs/user_data/rrockne/", "http://cgt.coh.org/", multiqc_path)
    )

    return(salmon_se)
}
