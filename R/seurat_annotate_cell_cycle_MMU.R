
#' Annotate cell cycle
#'
#' Get the mouse cell cycle markers from the package `scran`, and annotate the
#' Seutat object with the cell cycle labels for each cell.
#'
#' See <https://ucdavis-bioinformatics-training.github.io/2019-single-cell-RNA-sequencing-Workshop-UCD_UCSF/scrnaseq_analysis/scRNA_Workshop-PART1.html>
#' for more information.
#'
#' @name seurat_annotate_cell_cycle
#' @param seurat_object A Seurat object to annotate
#' @return A Surat object with cell cycle annotations (adds `S.Score`, `G2M.Score` and `Phase` to the cell metadata columns)
#' @author Denis O'Meally
#' @export
seurat_annotate_cell_cycle <- function(seurat_object) {
    # https://rdrr.io/cran/Seurat/man/CellCycleScoring.html
    # https://github.com/satijalab/seurat/issues/2493

    # Previously used this with success:
    # https://ucdavis-bioinformatics-training.github.io/2019-single-cell-RNA-sequencing-Workshop-UCD_UCSF/scrnaseq_analysis/scRNA_Workshop-PART1.html
    # But its a bit hacky and overly complex. Seurat has a built in function for this.

    # ############################################################################
    # ### This is the recommenced way to convert to mouse ids, but ensembl is down at time of testing
    # # Basic function to convert human to mouse gene names
    # convert_human_gene_list <- function(x) {
    #     # https://www.r-bloggers.com/2016/10/converting-mouse-to-human-gene-names-with-biomart-package/

    #     human <- biomaRt::useMart("ensembl", dataset = "hsapiens_gene_ensembl")
    #     mouse <- biomaRt::useMart("ensembl", dataset = "mmusculus_gene_ensembl")

    #     genesV2 <- biomaRt::getLDS(
    #         attributes = c("hgnc_symbol"),
    #         filters = "hgnc_symbol",
    #         values = x,
    #         mart = human,
    #         attributesL = c("mgi_symbol"),
    #         martL = mouse,
    #         uniqueRows = T
    #     )

    #     humanx <- unique(genesV2[, 2])

    #     # Print the first 6 genes found to the screen
    #     print(head(humanx))
    #     return(humanx)
    # }

    # g2m.genes <- convert_human_gene_list(Seurat::cc.genes$g2m.genes)
    # s.genes <- convert_human_gene_list(Seurat::cc.genes$s.genes)

    # In the interim, just make them all look like mice symbols
    g2m_genes <- gsub("(?<=\\b)([a-z])", "\\U\\1", tolower(Seurat::cc.genes$g2m.genes), perl = TRUE)
    s_genes <- gsub("(?<=\\b)([a-z])", "\\U\\1", tolower(Seurat::cc.genes$s.genes), perl = TRUE)

    seurat_object <- Seurat::CellCycleScoring(
        object = seurat_object,
        g2m.features = g2m_genes,
        s.features = s_genes
    )
    
    return(seurat_object)
}
