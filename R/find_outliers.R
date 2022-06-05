#' Run OUTRIDER pipeline on a SummarizedExperiment
#'
#' Runs the OUTRIDER pipeline on a summarised experiment, identifying outlying samples and genes.
#'
#' @name find_outliers
#' @param summarised_experiment a SummarisedExperiment to run the OUTRIDER pipeline on
#' @return a list of plots and tables
#' @author Denis O'Meally
#' @export
find_outliers <- function(summarised_experiment) {
    start_time <- Sys.time()
    # include sampleExclusionMask if duplicates are present

    # extract some objects ---------------------------------------------------------
    mat <- SummarizedExperiment::assay(summarised_experiment, "counts") |> round() |> as.matrix()
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
