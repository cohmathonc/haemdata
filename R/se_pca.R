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
#' @importFrom SummarizedExperiment assay colData
#' @importFrom dplyr case_when select contains
#' @importFrom PCAtools pca findElbowPoint biplot pairsplot eigencorplot getComponents
#' @importFrom janitor remove_constant
se_pca <- function(summarised_experiment, assay = "abundance", col_by = NULL) {
    # summarised_experiment <- AML.mRNA.2016_qc_se
    # col_by <- "treatment"

    mat <- SummarizedExperiment::assay(summarised_experiment, assay)
    mat <- mat[rowSums(mat) != 0, ] # remove rows with all zeros

    sample_sheet <- SummarizedExperiment::colData(summarised_experiment) |> as.data.frame()

    # colour the samples by `col_by`; If not given, colour by condition; if there is no condition,
    # colour by treatment; if there is no treatment, colour by sample
    if (is.null(col_by)) {
        condition <- dplyr::case_when(
            "condition" %in% colnames(sample_sheet) ~ "condition",
            "treatment" %in% colnames(sample_sheet) ~ "treatment",
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

    elbow <- PCAtools::findElbowPoint(p$variance)
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
    ) + labs(colour = condition)

    results$pairs_plot <- PCAtools::pairsplot(p,
        colby = condition,
        trianglelabSize = 8,
        plotaxes = F,
        margingaps = unit(c(0.0001, 0.0001, 0.0001, 0.0001), "cm")
    )

    results$correlation_plot <- PCAtools::eigencorplot(
        p,
        components = PCAtools::getComponents(p, 1:elbow),
        metavars = corr_metadata,
        col = c("white", "cornsilk1", "gold", "forestgreen", "darkgreen"),
        cexCorval = 0.7,
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
        signifSymbols = c("****", "***", "**", "*", ""),
        signifCutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1)
    )

    return(results)
}
