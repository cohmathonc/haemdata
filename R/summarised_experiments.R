# Functions for making and using SummarisedExperiments

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
#' @examples
#' \dontrun{
#' multiqc_path <- "path_to_a_nf-core_multiqc_report.html"
#' get_rnaseq_se(multiqc_path)
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
            dplyr::select(gene_id, gene_name, basepairs = width, gene_type = eval(rlang::sym(type))) |>
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
        stringr::str_detect(multiqc_path, "GENCODEr40") ~ "GENCODEr40"
    )

    rnaseq_release <- multiqc_path |>
        stringr::str_extract("(?<=\\/nfcore-rnaseq-)(.*?)(?=[_]G)")

    workflow <- ifelse(star_salmon, "qc", "salmon")

    object_name <- glue::glue("{run_folder}_{reference_genome}_{workflow}")

    S4Vectors::metadata(salmon_se) <- data.frame(
        object_name = object_name,
        run_folder = run_folder,
        reference_genome = reference_genome,
        rnaseq_release = rnaseq_release,
        workflow = workflow,
        multiqc_url = sub("/net/isi-dcnl/ifs/user_data/rrockne/", "http://cgt.coh.org/", multiqc_path),
        qc_removed = NA_character_
    )

    return(salmon_se)
}

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
    # sample_metadata <- metadata_mmu
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
                    dplyr::select(gene_id, gene_name, basepairs = width, gene_type = eval(rlang::sym(type))) |>
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

#' Run OUTRIDER pipeline on a SummarizedExperiment
#'
#' Runs the OUTRIDER pipeline on a summarised experiment, identifying outlying samples and genes.
#'
#' @name find_outliers_se
#' @param summarised_experiment a SummarisedExperiment to run the OUTRIDER pipeline on
#' @return a list of plots and tables
#' @author Denis O'Meally
#' @export
find_outliers_se <- function(summarised_experiment) {
    start_time <- Sys.time()
    # include sampleExclusionMask if duplicates are present

    # extract some objects
    mat <- SummarizedExperiment::assay(summarised_experiment, "counts") |>
        round() |>
        as.matrix()
    sample_sheet <- SummarizedExperiment::colData(summarised_experiment) |> tibble::as_tibble(rownames = "sampleID")
    # get annot_data
    annot_data <- S4Vectors::mcols(summarised_experiment) |> as.data.frame()
    results <- list()

    # Setup parallel processing
    ncores <- parallelly::availableCores()
    BiocParallel::register(BiocParallel::MulticoreParam(ncores, ncores * 2, progressbar = TRUE))

    # create OutriderDataSet object
    ods <- OUTRIDER::OutriderDataSet(countData = mat, colData = sample_sheet)

    S4Vectors::mcols(ods) <- annot_data

    # remove non-expressed genes
    ods <- OUTRIDER::filterExpression(ods, filterGenes = FALSE, savefpkm = TRUE)

    results$plotFPKM_raw <- OUTRIDER::plotFPKM(ods)
    results$plotExpressedGenes_raw <- OUTRIDER::plotExpressedGenes(ods)
    results$counts_raw <- OUTRIDER::counts(ods, normalized = FALSE)

    # normalise by libsize
    ods <- OUTRIDER::estimateSizeFactors(ods)

    # filter based on expression
    ods <- ods[S4Vectors::mcols(ods)$passedFilter, ]

    # counts corrected for library size
    results$counts_corr_libsize <- OUTRIDER::counts(ods, normalized = TRUE) #

    # remove samples wth very low mapping
    #    ods <- OUTRIDER::estimateSizeFactors(ods)
    ods <- ods[, OUTRIDER::sizeFactors(ods) > 0.1]
    results$removed_low_mapping <- setdiff(colnames(summarised_experiment), colnames(ods))

    # find number of encoding dimensions
    ods <- OUTRIDER::findEncodingDim(ods, BPPARAM = BiocParallel::bpparam(), params = seq(2, 20, by = 2))

    # counts corrected for confounds using PCA
    ods <- OUTRIDER::controlForConfounders(ods, implementation = "PCA", q = OUTRIDER::getBestQ(ods), BPPARAM = BiocParallel::bpparam())
    results$counts_corr_pca <- OUTRIDER::counts(ods, normalized = TRUE) #

    # counts corrected for confounds using autoencoder
    ods <- OUTRIDER::controlForConfounders(ods, q = OUTRIDER::getBestQ(ods), BPPARAM = BiocParallel::bpparam())
    results$counts_corr_ae <- OUTRIDER::counts(ods, normalized = TRUE) #

    # P-values & z-scores
    ods <- OUTRIDER::computePvalues(ods, alternative = "two.sided", method = "BH", BPPARAM = BiocParallel::bpparam())
    ods <- OUTRIDER::computeZscores(ods)
    results$plot_expressed_genes <- OUTRIDER::plotExpressedGenes(ods)
    results$plot_power <- OUTRIDER::plotPowerAnalysis(ods)
    results$plot_encoding_dimensions <- OUTRIDER::plotEncDimSearch(ods)
    results$plot_aberrant_per_sample <- OUTRIDER::plotAberrantPerSample(ods, padjCutoff = 0.05)
    results$aberrant_per_sample <- OUTRIDER::aberrant(ods, by = "sample")
    results$aberrant_per_gene <- OUTRIDER::aberrant(ods, by = "gene")
    results$results <- OUTRIDER::results(ods)
    results$normalizationFactors <- OUTRIDER::normalizationFactors(ods)
    results$counts_raw_flt <- OUTRIDER::counts(ods, normalized = FALSE) # filtered counts
    results$counts_raw_flt <- OUTRIDER::counts(ods, normalized = FALSE) # filtered counts
    results$annot_data <- annot_data
    results$ncores <- ncores
    end_time <- Sys.time()
    results$run_time <- (end_time - start_time)
    return(results)
}

#' Remove poorly mapped samples
#'
#' If the QC metrics contain a `uniquely_mapped_percent` column, remove samples below the specified `mapping_threshold`.
#'
#' @title qc_filter_se
#' @param summarised_experiment a SummarisedExperiment from a "QC" run (`qc = TRUE`
#' set in [run_nf_core_rnaseq()]
#' @param mapping_threshold threshold below which samples will be discarded (STAR uniquely-mapped reads). Default: 5 (%)
#' @return a SummarisedExperiment. The `metadata` slot of the SummarisedExperiment is updated with a record of the samples
#' removed by the filter (`qc_removed`).
#' @author Denis O'Meally
#' @export
qc_filter_se <- function(summarised_experiment, mapping_threshold = 5) {

    # Check that the summarised_experiment is from a qc run
    if (!S4Vectors::metadata(summarised_experiment)["workflow"] == "qc") {
        stop("The summarised_experiment is not from a qc run")
    }
    col_pc_mapped <- grepl("uniquely_mapped_percent", colnames(SummarizedExperiment::colData(summarised_experiment)))
    object_name <- S4Vectors::metadata(summarised_experiment)["object_name"]
    n_removed <- 0

    ## Filter out samples with low mapping rate
    if (any(col_pc_mapped)) {
        percent_mapped <- SummarizedExperiment::colData(summarised_experiment)[, col_pc_mapped]

        filtered_se <- subset(
            summarised_experiment,
            select = (percent_mapped > mapping_threshold)
        )
    } else {
        stop("'uniquely_mapped_percent' not found in summarised_experiment colData...")
    }

    ## Update metadata with removed sample if needs be
    if (ncol(summarised_experiment) != ncol(filtered_se)) {
        metadata <- S4Vectors::metadata(filtered_se)

        metadata["qc_removed"] <- paste(
            c(
                na.omit(metadata["qc_removed"]),
                setdiff(colnames(summarised_experiment), colnames(filtered_se))
            ),
            collapse = ", "
        )
        S4Vectors::metadata(filtered_se) <- metadata

        n_removed <- abs(ncol(filtered_se) - ncol(summarised_experiment))
    }

    # ## Describe what was done

    message <- ifelse(n_removed > 0,
        paste0(object_name, ": removed ", n_removed, " samples with low mapping rate (", mapping_threshold, "%)"),
        paste0(object_name, ": no samples removed due to poor mapping threshold (", mapping_threshold, "%)")
    )

    message(message)

    return(filtered_se)
}

#' Plot human transgene expression
#'
#' Plots expression of genes matching the pattern `HSA_.*_gene` from a SummarisedExpression object.
#' Log2TPM is values are plotted from the `abundance` assay.
#'
#' @name plot_transgenes
#' @param summarised_experiment a SummarisedExperiment object
#' @param assay the matrix to plot, Default: "abundance"
#' @param group_by a string indicating the column to group by
#' @return a ggplot2 plot
#' @author Denis O'Meally
#' @export
plot_transgenes_se <- function(summarised_experiment, assay = "abundance", group_by = "treatment") {
    mat <- SummarizedExperiment::assay(summarised_experiment, assay)

    col_data <- SummarizedExperiment::colData(summarised_experiment) |>
        as.data.frame() |>
        droplevels()

    mat |>
        tibble::rownames_to_column("gene_id") |>
        dplyr::filter(stringr::str_detect(gene_id, "HSA_.*_gene")) |>
        tidyr::pivot_longer(!gene_id, names_to = "sample") |>
        dplyr::left_join(col_data, by = "sample") |>
        ggplot2::ggplot(ggplot2::aes(y = log2(value + 0.1), x = weeks, group = mouse_id, col = eval(rlang::sym(group_by)))) +
        ggplot2::geom_line(alpha = 0.3) +
        ggplot2::geom_point(alpha = 0.3) +
        ggplot2::geom_smooth(ggplot2::aes(group = eval(rlang::sym(group_by))), span = 0.4) +
        ggplot2::facet_wrap("gene_id") +
        ggplot2::ggtitle("Abundance of human transgenes") +
        ggplot2::ylab("log2(TPM)") +
        ggplot2::xlab("time point") +
        ggsci::scale_color_d3("category20") +
        ggplot2::theme_minimal() +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
}

#' @title Perform PCA on a SummarisedExperiment
#' @description Perform PCA on the `assay` of a SummarisedExperiment object and make a biplot, pairs plot and
#' correlation of metadata columns with PCs.
#' @param summarised_experiment a `SummarisedExperiment` object
#' @param assay the assay to perform PCA on, Default: "abundance"
#' @param col_by column name of colData to use for colouring samples in PCA plots. If none is given,
#' samples will be coloured by the "condition", "treatment" or "sample" column if any exists, in that order. Default: NULL
#' @return a list of ggplot2 plots
#' @details Could be rewritten, but this is a quick and dirty way to get something working using \code{vignette("PCAtools", package = "PCAtools")}.
#' @examples
#' \dontrun{
#' if (interactive()) {
#'     se_pca(AML.mRNA.2016_qc_se)
#' }
#' }
#' @seealso
#'  \code{\link[SummarizedExperiment]{SummarizedExperiment-class}}
#'  \code{\link[dplyr]{case_when}}, \code{\link[dplyr]{select}}, \code{\link[dplyr]{reexports}}
#'  \code{\link[PCAtools]{pca}}, \code{\link[PCAtools]{findElbowPoint}}, \code{\link[PCAtools]{biplot}}, \code{\link[PCAtools]{pairsplot}}, \code{\link[PCAtools]{eigencorplot}}, \code{\link[PCAtools]{getComponents}}
#'  \code{\link[janitor]{remove_constant}}
#' @rdname se_pca
#' @author Denis O'Meally
#' @export
pca_se <- function(summarised_experiment, assay = "abundance", col_by = NULL) {
    # summarised_experiment <- AML.mRNA.2016_qc_se
    # col_by <- "treatment"

    mat <- SummarizedExperiment::assay(summarised_experiment, assay)
    mat <- mat[rowSums(mat) != 0, ] # remove rows with all zeros

    sample_sheet <- SummarizedExperiment::colData(summarised_experiment) |> as.data.frame()

    # colour the samples by `col_by`; If not given, colour by condition; if there is no condition,
    # colour by treatment; if there is no treatment, colour by batch
    # colour by batch; if there is no treatment, colour by sample
    if (is.null(col_by)) {
        condition <- dplyr::case_when(
            "condition" %in% colnames(sample_sheet) ~ "condition",
            "treatment" %in% colnames(sample_sheet) ~ "treatment",
            "batch" %in% colnames(sample_sheet) ~ "treatment",
            TRUE ~ "sample"
        )
    } else {
        condition <- col_by
    }

    p <- PCAtools::pca(as.matrix(mat),
        metadata = sample_sheet,
        scale = TRUE,
        removeVar = 0.1
    )

    # Use the elbow point to determine the number of components to show, if fewer than 10
    num_components <- min(PCAtools::findElbowPoint(p$variance), 10)

    corr_metadata <- sample_sheet |>
        dplyr::select(-dplyr::contains(c("sample", "fastq_1", "fastq_2"))) |>
        janitor::remove_constant(na.rm = TRUE) |>
        colnames()

    results <- list()
    results$biplot <- PCAtools::biplot(p,
        colby = condition,
        lab = sample_sheet$batch,
        drawConnectors = F,
        labSize = 2,
        legendPosition = "right",
        legendLabSize = 9,
        axisLabSize = 8
    ) + ggplot2::labs(colour = condition)

    results$pairs_plot <- PCAtools::pairsplot(p,
        colby = condition,
        trianglelabSize = 8,
        plotaxes = F,
        margingaps = grid::unit(c(0.0001, 0.0001, 0.0001, 0.0001), "cm")
    )

    results$correlation_plot <- PCAtools::eigencorplot(
        p,
        components = PCAtools::getComponents(p, 1:num_components),
        metavars = corr_metadata,
        col = c("white", "cornsilk1", "gold", "forestgreen", "darkgreen"),
        cexCorval = 0.6,
        main = bquote("Pearson" ~ r^2),
        cexMain = 0.9,
        fontCorval = 1,
        fontLabY = 3,
        cexLabY = 0.8,
        posLab = "bottomleft",
        rotLabX = 45,
        scale = TRUE,
        plotRsquared = TRUE,
        corFUN = "pearson",
        corUSE = "pairwise.complete.obs",
        corMultipleTestCorrection = "BH",
        signifSymbols = c( '**', '*', ''),
        signifCutpoints = c(0, 0.01, 0.05, 1)
    )

    return(results)
}

#' @title Write an expression matrix to file
#' @name make_tpm_matrix
#' @description Write an TPM expression matrix to a CSV file in `extdata`
#' @param summarised_experiment a `SummarisedExperiment` object
#' @param drop samples (columns) to remove
#' @param tpm minimum TPM for a gene to be kept, Default: 1
#' @param samples minimum number of samples a gene must be present in to be kept, Default: 5
#' @return a list
#' * `tpm_matrix` an expression matrix, with genes as rows nad samples as columns
#' * `name` a name indicating the `tpm` and `sample` thresholds used
#' @details This function makes an expression matrix, filtering out genes with
#' fewer than `tpm` reads in at least `samples` samples.
#' @examples
#' \dontrun{
#' if (interactive()) {
#'     make_tpm_matrix(cml_mrna_GRCm38_HLT, c("F545.7", "X490.17"), tpm = 1, samples = 5)
#' }
#' }
#' @seealso
#'  \code{\link[SummarizedExperiment]{SummarizedExperiment-class}}
#'  \code{\link[utils]{write.table}}
#' @rdname make_tpm_matrix
#' @export
make_tpm_matrix <- function(summarised_experiment, drop = NULL, tpm = 1, samples = 5) {
    object_name <- summarised_experiment@metadata$object_name
    # Subset to samples to drop
    summarised_experiment <- summarised_experiment[, !(colnames(summarised_experiment) %in% drop)]
    # Make a CSV of TPMs, keeping genes with > 1 TPM in 5 samples, and transgenes
    mat <- SummarizedExperiment::assay(summarised_experiment, "abundance")
    filter <- rowSums(mat >= tpm) >= samples
    filter[grep("^HSA_", rownames(summarised_experiment))] <- TRUE
    filtered <- summarised_experiment[filter, ]
    tpm_matrix <- SummarizedExperiment::assay(filtered, "abundance") |> tibble::rownames_to_column("gene_id")
    name <- paste0(object_name, "_", tpm, "tpm_in_", samples, "samples")

    return(list(tpm_matrix = tpm_matrix, name = name))
}

#' Merge a pair of `SummarisedExperiments`
#'
#' Takes a pair of `SummarisedExperiments` produced by [get_rnaseq_se()],
#' and merges them into a single `SummarisedExperiment`
#'
#' @name merge_mrna
#' @param summarised_experiment1, name of a `SummarisedExperiment` object to merge
#' @param summarised_experiment2, name of the other `SummarisedExperiment` object to merge
#' @param drop a vector of sample names to remove from the merged `SummarisedExperiment`
#' @return a `SummarisedExperiment` object
#' @author Denis O'Meally
#' @export
merge_mrna <- function(summarised_experiment1, summarised_experiment2, drop = NULL) {

    # metadata1 <- S4Vectors::metadata(cml_mrna_2021_m38)
    # metadata2 <- S4Vectors::metadata(cml_mrna_2022_m38)
    # drop <- c("F545.7", "X490.17")

    metadata1 <- S4Vectors::metadata(summarised_experiment1)
    metadata2 <- S4Vectors::metadata(summarised_experiment2)

    # Check the $reference_genome version matches
    identical(
        metadata1$reference_genome,
        metadata2$reference_genome
    )
    reference_genome <- metadata1$reference_genome

    # Check the pipeline version matches
    #    identical(metadata1$pipeline, metadata2$pipeline)
    pipeline <- metadata1$pipeline

    # Check the project name matches
    identical(
        stringr::str_match(metadata1$project, "$*(.*_.*?)_")[, 2],
        stringr::str_match(metadata2$project, "$*(.*_.*?)_")[, 2]
    )
    project <- stringr::str_match(metadata1$project, "$*(.*_.*?)_")[, 2]

    # update the metadata
    metadata <- list(
        project = project,
        reference_genome = reference_genome,
        pipeline = pipeline,
        date = date()
    )

    # merge the datasets
    summarised_experiment <- cbind(summarised_experiment1, summarised_experiment2)
    S4Vectors::metadata(summarised_experiment) <- metadata

    # remove samples that failed QC
    summarised_experiment[, !(colnames(summarised_experiment) %in% drop)]

    # name the output
    out_name <- paste(tolower(project), reference_genome, sep = "_")

    # write to the package data folder
    helpeRs::write_data(out_name, summarised_experiment)

    # write a CSV to extdata
    tpm_matrix_csv <- helpeRs::write_tpm_matrix(summarised_experiment, tpm = 1, samples = 5)

    # publish to Teams
    team <- Microsoft365R::get_team("PSON AML State-Transition")
    channel <- team$get_channel("haemdata")
    channel$send_message(paste0("New dataset available for ", out_name, ": ", basename(tpm_matrix_csv)),
        attachments = tpm_matrix_csv
    )

    # PINS https://pins.rstudio.com/reference/board_ms365.html
    # # A board in a SharePoint Online document library
    # sp <- Microsoft365R::get_sharepoint_site("PSON AML State-Transition")
    # chan_folder <- chan$get_folder()
    # board <- pins::board_ms365(chan_folder, "haemdata")

    return(summarised_experiment)
}

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

    #   sample_metadata <- metadata_mmu

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

    gtf <- readr::read_file(multiqc_path) |>
        stringr::str_extract("\\/ref_genome\\/igenomes\\/.*\\.gtf")

    type <- ifelse(stringr::str_detect(gtf, "GENCODE"), "gene_type", "gene_biotype")

    row_data <- rtracklayer::import(gtf) |>
        tibble::as_tibble() |>
        dplyr::select(gene_id, gene_name, basepairs = width, gene_type = eval(rlang::sym(type))) |>
        dplyr::group_by(gene_id) |>
        dplyr::slice(which.max(basepairs)) |>
        dplyr::ungroup()

    # load SummarisedExperiment
    salmon_se <- readRDS(glue::glue("{out_path}/salmon/salmon.merged.gene_counts.rds"))

    # Update rowData
    SummarizedExperiment::rowData(salmon_se) <- dplyr::left_join(
        data.frame(gene_id = rownames(salmon_se)),
        row_data,
        by = "gene_id"
    )

    # Update colData
    SummarizedExperiment::colData(salmon_se) <-
        dplyr::left_join(
            dplyr::left_join(
                data.frame(sample = colnames(salmon_se)),
                sample_metadata,
                by = "sample"
            ),
            qc_data,
            by = "sample"
        ) |>
        `rownames<-`(colnames(salmon_se)) |>
        S4Vectors::DataFrame()

    # Add some metadata

    run_folder <- multiqc_path |>
        stringr::str_extract("(?<=haemdata-nf-core-cache\\/)(.*?)(?=\\/)") |>
        stringr::str_replace_all("[.]", "_")

    reference_genome <- dplyr::case_when(
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
        workflow = workflow,
        qc_removed = NA_character_
    )

    return(salmon_se)
}