#' Plot human transgene expression
#'
#' Plots expression of genes matching the pattern `HSA_.*_gene` from a SummarisedExpression object.
#' Log2TPM is values are plotted from the `abundance` assay.
#'
#' @name plot_transgenes 
#' @param summarised_experiment a SummarisedExperiment object with an `abundance` assay
#' @param group_by a string indicating the column to group by
#' @return a ggplot2 plot
#' @author Denis O'Meally
#' @export
plot_transgenes <- function(summarised_experiment, assay = "abundance", group_by = "treatment") {

    mat <- SummarizedExperiment::assay(summarised_experiment, assay)

    col_data <- SummarizedExperiment::colData(summarised_experiment) |> as.data.frame() |> droplevels()

    mat |>
        tibble::rownames_to_column("gene_id") |>
        dplyr::filter(stringr::str_detect(gene_id, "HSA_.*_gene")) |>
        tidyr::pivot_longer(!gene_id, names_to = "sample") |>
        dplyr::left_join(col_data, by = "sample") |>
        ggplot(aes(y = log2(value + 0.1), x = weeks, group = mouse_id, col = !!sym(group_by))) +
        geom_line(alpha = 0.3) +
        geom_point(alpha = 0.3) +
        geom_smooth(aes(group = !!sym(group_by)), span = 0.4) +
        facet_wrap("gene_id") +
        ggtitle("Abundance of human transgenes") +
        ylab("log2(TPM)") +
        xlab("time point") +
        ggsci::scale_color_d3("category20") +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
}
