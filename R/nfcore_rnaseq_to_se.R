#' nf-core/rnaseq to SummarisedExperiment
#'
#' Makes a SummarisedExperiment from the pipeline run given by `multiqc_path`.
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
#' Assays slots for both `counts` and `abundance` come from the Salmon output `{out_path}/salmon/salmon.merged.gene_counts.rds`
#'
#' @name nfcore_rnaseq_to_se
#' @param multiqc_path path to the multiqc report of a successful nf-core/rnaseq pipeline run
#' @param sample_metadata a data.frame to combine with the QC metrics from the pipeline run as colData in the SummarisedExperiment
#' @return a SummarisedExperiment
#' @author Denis O'Meally
#' @export
nfcore_rnaseq_to_se <- function(multiqc_path, sample_metadata) {
    #  multiqc_path <- AML.mRNA.2016_qc
    #   multiqc_path <- all_mice.mRNA_salmon_GRCm38_HLT

    #   sample_metadata <- minimal_metadata_mmu

    # remove library metadata from sample_metadata
    sample_metadata <- sample_metadata |>
        select(-c(fastq_1, fastq_2, strandedness)) |>
        distinct()

    star_salmon <- ifelse(stringr::str_detect(multiqc_path, "star_salmon"), TRUE, FALSE)

    out_path <- ifelse(star_salmon,
        gsub("/multiqc/star_salmon/multiqc_report.html", "", multiqc_path),
        gsub("/multiqc/multiqc_report.html", "", multiqc_path)
    )

    qc_path <- ifelse(star_salmon,
        glue::glue("{out_path}/multiqc/star_salmon/multiqc_data/multiqc_general_stats.txt"),
        glue::glue("{out_path}/multiqc/multiqc_data/multiqc_salmon.txt"))

    qc_data <- utils::read.csv(qc_path, sep = "\t", header = TRUE, stringsAsFactors = FALSE)|>
            dplyr::mutate(sample = gsub("_1$|_2$", "", Sample)) |> # for paired end data, take the mean
            dplyr::group_by(sample) %>%                            # of the FastQC stats, which are reported for each read pair
            dplyr::mutate(
                dplyr::across(dplyr::starts_with("FastQC"), ~ replace_na(.x, mean(.x, na.rm = TRUE)))
            ) |>
            dplyr::filter(!grepl("_1$|_2$", Sample)) |>
                dplyr::ungroup() |>
                janitor::remove_constant(na.rm = TRUE) |>           # remove invariant columns and tidy up column names
                dplyr::select(-dplyr::contains(c("Sample", "start_time", "end_time"), ignore.case = FALSE)) %>%
                setNames(gsub(".*\\.generalstats.", "", names(.)))

    gtf <- readr::read_file(multiqc_path) |>
        stringr::str_extract("\\/ref_genome\\/igenomes\\/.*\\.gtf")

    type <- ifelse(stringr::str_detect(gtf, "GENCODE"), "gene_type", "gene_biotype")

    row_data <- rtracklayer::import(gtf) |>
        tibble::as_tibble() |>
        dplyr::select(gene_id, gene_name, basepairs = width, gene_type = !!sym(type)) |>
        dplyr::group_by(gene_id) |>
        dplyr::slice(which.max(basepairs)) |>
        dplyr::ungroup()

    # load SummarisedExperiment
    salmon_se <- readRDS(glue::glue("{out_path}/salmon/salmon.merged.gene_counts.rds"))

    # Update rowData
    SummarizedExperiment::rowData(salmon_se) <- dplyr::left_join(
        data.frame(gene_id = rownames(salmon_se)),
        row_data, by = "gene_id"
    )

    # Update colData
    SummarizedExperiment::colData(salmon_se) <-
        dplyr::left_join(
            dplyr::left_join(
                data.frame(sample = colnames(salmon_se)),
                sample_metadata,
                by = "sample"
            ),
            qc_data, by = "sample"
        ) |>
        `rownames<-`(colnames(salmon_se)) |>
            S4Vectors::DataFrame()

    # Add some metadata

    run_folder <- multiqc_path |>
        stringr::str_extract("(?<=haemdata-nf-core-cache\\/)(.*?)(?=\\/)") |>
        stringr::str_replace_all("[.]", "_")

    reference_genome <-  dplyr::case_when(
            grepl("GRCm38", multiqc_path) ~ "GRCm38_HLT",
            grepl("GENCODEv33", multiqc_path) ~ "GENCODEv33_HLT",
            grepl("GENCODEr40", multiqc_path) ~ "GENCODEr40"
        )

    pipeline <- multiqc_path |>
        stringr::str_extract("(?<=\\/nfcore-rnaseq-)(.*?)(?=[_]G)")

    workflow <- ifelse(star_salmon, "qc", "salmon")

    abundance_estimation <- ifelse(star_salmon, "STAR-Salmon", "Salmon")

    object_name <- glue::glue("{run_folder}_{reference_genome}_{workflow}")

    S4Vectors::metadata(salmon_se) <- data.frame(
        object_name = object_name,
        run_folder = run_folder,
        reference_genome = reference_genome,
        pipeline = pipeline,
        abundance_estimation = abundance_estimation,
        workflow = workflow
    )

    return(salmon_se)
}
