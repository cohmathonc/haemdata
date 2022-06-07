#' Update SummarisedExperiment
#'
#' Updates or replaces the `colData` or `rowData` of a SummarisedExperiment.
#'
#' Genes and associated metadata come from the `gtf` logged in the multiqc report.
#'
#' Sample metadata is matched by the `sample` ID from the `sample_metadata` table and the
#' QC metrics from the pipeline run. If the multiqc path refers to a full QC run (`qc = TRUE`
#' set in [run_nf_core_rnaseq()], the QC metrics come from
#' `{out_path}/multiqc/star_salmon/multiqc_data/multiqc_general_stats.txt`. If a Salmon-only
#' run (`qc = FALSE` set in [run_nf_core_rnaseq()]), the QC metrics come from
#' `{out_path}/multiqc/multiqc_data/multiqc_salmon.txt`.
#'
#' @name annotate_se
#' @param summarised_experiment the SummarisedExperiment to update
#' @param sample_metadata a data.frame to replace or amend colData
#' @param multiqc_path the path to a multiqc report from which to extract the QC metrics or gtf
#' @param gtf a logical indicating whether to extract the gtf file from the multiqc report and replace gene annotations for rowData
#' @param qc a logical indicating whether to extract the QC metrics from the multiqc report and replace or amend gene annotations to rowData
#' @param mode one of c("replace", "amend"), should the provided data replace or amend the existing elements? Default: "replace"
#' @return a SummarisedExperiment
#' @author Denis O'Meally
#' @export
annotate_se <- function(summarised_experiment, sample_metadata, multiqc_path = NULL, gtf = FALSE, qc = TRUE, mode = "replace") {
    # multiqc_path <- AML.mRNA.2016_qc
    #   multiqc_path <- all_mice.mRNA_salmon_GRCm38_HLT
    # summarised_experiment <- AML.mRNA.2016_qc_se
    # sample_metadata <- minimal_metadata_mmu
    # sample_metadata <- AML.mRNA.2016_qc_se_outliers$

    # remove library metadata from sample_metadata
    sample_metadata <- sample_metadata |>
        dplyr::select(!dplyr::matches("fastq_1|fastq_2|strandedness")) |>
        dplyr::distinct()

    existing_colData <- SummarizedExperiment::colData(summarised_experiment) |>
        as.data.frame() |>
        tibble::rownames_to_column("sample")

    # get the QC metrics from the multiqc report; update gene annotation from the multiqc gtf
    if (!is.null(multiqc_path)) {
        if (gtf == FALSE && qc == FALSE) {
            stop("A multiqc path was provided but qc & gtf = FALSE - one must be TRUE")
        } else {
            star_salmon <- ifelse(stringr::str_detect(multiqc_path, "star_salmon"), TRUE, FALSE)

            out_path <- ifelse(star_salmon,
                gsub("/multiqc/star_salmon/multiqc_report.html", "", multiqc_path),
                gsub("/multiqc/multiqc_report.html", "", multiqc_path)
            )

            # extract qc data if qc = TRUE
            if (qc == TRUE) {
                qc_path <- ifelse(star_salmon,
                    glue::glue("{out_path}/multiqc/star_salmon/multiqc_data/multiqc_general_stats.txt"),
                    glue::glue("{out_path}/multiqc/multiqc_data/multiqc_salmon.txt")
                )

                qc_data <- utils::read.csv(qc_path, sep = "\t", header = TRUE, stringsAsFactors = FALSE) |>
                    dplyr::mutate(sample = gsub("_1$|_2$", "", Sample)) |> # for paired end data, take the mean
                    dplyr::group_by(sample) %>% # of the FastQC stats, which are reported for each read pair
                    dplyr::mutate(
                        dplyr::across(dplyr::starts_with("FastQC"), ~ replace_na(.x, mean(.x, na.rm = TRUE)))
                    ) |>
                    dplyr::filter(!grepl("_1$|_2$", Sample)) |>
                    dplyr::ungroup() |>
                    janitor::remove_constant(na.rm = TRUE) |> # remove invariant columns and tidy up column names
                    dplyr::select(-dplyr::contains(c("Sample", "start_time", "end_time"), ignore.case = FALSE)) %>%
                    setNames(gsub(".*\\.generalstats.", "", names(.)))
            }
            # Update gene annotations if gtf = TRUE
            if (gtf == TRUE) {
                gtf <- readr::read_file(multiqc_path) |>
                    stringr::str_extract("\\/ref_genome\\/igenomes\\/.*\\.gtf")

                type <- ifelse(stringr::str_detect(gtf, "GENCODE"), "gene_type", "gene_biotype")

                row_data <- rtracklayer::import(gtf) |>
                    tibble::as_tibble() |>
                    dplyr::select(gene_id, gene_name, basepairs = width, gene_type = !!sym(type)) |>
                    dplyr::group_by(gene_id) |>
                    dplyr::slice(which.max(basepairs)) |>
                    dplyr::ungroup()

                # Update rowData
                SummarizedExperiment::rowData(summarised_experiment) <- dplyr::left_join(
                    data.frame(gene_id = rownames(summarised_experiment)),
                    row_data,
                    by = "gene_id"
                )
            }
        }
    }
    # Update colData
    if (mode == "replace" && qc == TRUE) {
        SummarizedExperiment::colData(summarised_experiment) <-
            dplyr::left_join(
                dplyr::left_join(
                    data.frame(sample = colnames(summarised_experiment)),
                    sample_metadata,
                    by = "sample"
                ),
                qc_data,
                by = "sample"
            ) |>
            `rownames<-`(colnames(summarised_experiment)) |>
            S4Vectors::DataFrame()
            message("Replacing colData with sample_metadata, including QC metrics...")
    } else if (mode == "replace" && qc == FALSE) {
        SummarizedExperiment::colData(summarised_experiment) <-
            dplyr::left_join(
                data.frame(sample = colnames(summarised_experiment)),
                sample_metadata,
                by = "sample"
            ) |>
            `rownames<-`(colnames(summarised_experiment)) |>
            S4Vectors::DataFrame()
            message("Replacing colData with sample_metadata, no QC metrics...")
    } else if (mode == "amend" && qc == TRUE) {
        SummarizedExperiment::colData(summarised_experiment) <-
            dplyr::left_join(
                dplyr::left_join(
                    dplyr::left_join(
                        data.frame(sample = colnames(summarised_experiment)),
                        sample_metadata,
                        by = "sample"
                    ),
                    qc_data,
                    by = "sample"
                ),
                existing_colData,
                by = "sample"
            ) |>
            `rownames<-`(colnames(summarised_experiment)) |>
            S4Vectors::DataFrame()
            message("Amending colData with sample_metadata, including QC metrics...")
    } else if (mode == "amend" && qc == FALSE) {
        SummarizedExperiment::colData(summarised_experiment) <-
            dplyr::left_join(
                dplyr::left_join(
                    data.frame(sample = colnames(summarised_experiment)),
                    sample_metadata,
                    by = "sample"
                ),
                existing_colData,
                by = "sample"
            ) |>
            `rownames<-`(colnames(summarised_experiment)) |>
            S4Vectors::DataFrame()
            message("Amending colData with sample_metadata, no QC metrics...")
    } else {
        stop("Invalid mode: mut be one of 'replace' or 'amend'")
    }

    return(summarised_experiment)
}
