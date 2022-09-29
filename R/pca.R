# TODO https://github.com/pkimes/pkimes/blob/master/R/ggpcaplot.R

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

all_mice <- get_pin("mmu_mrna_2016_2022_GENCODEm28_HLT_qc.rds")
summarised_experiment <- subset(all_mice,
    select = grepl("2018", SummarizedExperiment::colData(all_mice)$batch))
assay <- "abundance"
col_by <- "treatment"

    mat <- SummarizedExperiment::assay(summarised_experiment, assay) |>
        as.matrix()

    sample_sheet <- SummarizedExperiment::colData(summarised_experiment) |>
        as.data.frame()

    # colour the samples by `col_by`; If not given, colour by condition;
    # if there is no condition, colour by treatment; if there is no treatment,
    # colour by batch colour by batch; if there is no treatment, colour by
    # sample
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

groups <- sample_sheet[[condition]]

col <- c("red", "black")

pca_fit <- mat |>
    genefilter::varFilter(var.cutoff = 0.10) |> 
    t() |>
    scale() |>
    prcomp()

pc1_var_exp <- broom::tidy(pca_fit, matrix = "d")[1, 3] * 100 |> round(1)
pc2_var_exp <- broom::tidy(pca_fit, matrix = "d")[2, 3] * 100 |> round(1)

biplot <- pca_fit$x |>
    as.data.frame() |>
    ggplot2::ggplot(ggplot2::aes(PC1, PC2, color = groups)) +
    ggplot2::geom_point(size = 1.5) +
    ggplot2::scale_color_manual(values = col) +
    ggplot2::ylab(glue::glue("PC1, {pc1_var_exp}% variation")) +
    ggplot2::xlab(glue::glue("PC2, {pc2_var_exp}% variation"))
biplot

pairs_plot <- pca_fit$x |>
    as.data.frame() |>
    GGally::ggpairs(
        columns = 1:5,
        mapping = ggplot2::aes(
            colour = sample_sheet$treatment,
            alpha = 0.5)) +
    ggplot2::scale_fill_manual(values = col) +
    ggplot2::scale_color_manual(values = col)
pairs_plot


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
        ggplot2::ggplot(ggplot2::aes(y = log2(value + 0.1), x = sample_weeks, group = mouse_id, col = eval(rlang::sym(group_by)))) +
        ggplot2::geom_line(alpha = 0.3) +
        ggplot2::geom_point(alpha = 0.3) +
        ggplot2::geom_smooth(ggplot2::aes(group = eval(rlang::sym(group_by))), span = 0.4) +
        ggplot2::facet_wrap("gene_id") +
        ggplot2::ggtitle("Abundance of human transgenes") +
        ggplot2::ylab("log2(TPM)") +
        ggplot2::xlab("Weeks") +
        ggsci::scale_color_d3("category20") +
        ggplot2::theme_minimal() +
        ggplot2::theme(
            axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
            legend.title = element_blank()
        )
}


    p <- PCAtools::pca(as.matrix(mat),
        metadata = sample_sheet,
        scale = TRUE,
        removeVar = 0.1
    )

    # Use the elbow point to determine the number of components to show, if fewer than 10
    num_components <- min(PCAtools::findElbowPoint(p$variance), 10)

    corr_metadata <- sample_sheet |>
        # dplyr::select(-c("sample", "fastq_1", "fastq_2")) |>
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
        signifSymbols = c("**", "*", ""),
        signifCutpoints = c(0, 0.01, 0.05, 1)
    )

    return(results)
}



#' @title PCA Scores Plot Using ggplot2
#'
#' This function is a simple wrapper for plotting two-way scatter plots
#' of principal component scores using \code{GGally::ggpairs}.
#'
#' @param x a matrix of data (n x p)
#' @param npc an integer number of principal components
#' @param labs labels for the different clusters
#' @param idiag a boolean whether to include
#'
#' @return
#' a \code{ggpairs} plot of the PC scores
#'
#' @import GGally
#' @author Patrick Kimes
#' @export

ggpcaplot <- function(x, npc, labs = NULL, idiag = TRUE) {
    pca <- prcomp(x)
    if (is.null(labs)) {
        ggp <- ggpairs(data.frame(pca$x[, 1:npc]),
            params = c(alpha = 1, color = "blue", fill = "red"),
            upper = "blank", lower = list(continuous = "points")
        ) # ,
        ## diag=list(continuous='blank', discrete='blank'))#ifelse(idiag, 'density', 'blank')))
    } else {
        ggp <- ggpairs(data.frame(pca$x[, 1:npc], labs = as.factor(labs)),
            params = c(alpha = 1), color = "labs", columns = 1:npc,
            upper = "blank", lower = list(continuous = "points")
        ) # ,
        ## diag=list(continuous='blank', discrete='blank'))#ifelse(idiag, 'density', 'blank')))
    }
    for (ip in 1:npc^2) {
        ggp$plots[[ip]] <- paste(ggp$plots[[ip]], "+ theme_bw()")
        ggp$plots[[ip]] <- paste(ggp$plots[[ip]], "+ scale_color_brewer(palette='Set1')")
    }

    ## remove diagonal if not wanted
    if (!idiag) {
        for (ip in 1:npc) {
            ggp$plots[[(ip - 1) * (npc + 1) + 1]] <- "ggally_blank()"
        }
    }
    ggp
}
