#' Annotate mouse bone marrow cell type
#'
#' Use the [`clustifyr::clustify()`] function to label cell types according
#' to expression profiles identified by Dolgalev & Tikhonova (2021) in a
#' meta-analysis of five mouse bone marrow scRNAseq datasets.
#'
#' @name seurat_annotate_type
#' @param seurat_object a list of seurat objects
#' @return a seurat object
#' @source https://www.frontiersin.org/articles/10.3389/fcell.2021.622519/full
#' @author Denis O'Meally
#' @export
seurat_annotate_type <- function(seurat_object) {
    dolgalev_2021 <- read.csv(row.names = 1, "data-raw/Dolgalev_2021_expression-mean-labels-harmonized.csv")
    
    Seurat::DefaultAssay(seurat_object) <- "integrated"

    seurat_object <- clustifyr::clustify(
        input = seurat_object,
        cluster_col = "seurat_clusters",
        ref_mat = dolgalev_2021,
        query_genes = Seurat::VariableFeatures(seurat_object)[1:750]
)
# copy type column to clustifyr_dolgalev_2021
Seurat::Idents(seurat_object) <- "type"
Seurat::DimPlot(
    seurat_object,
    reduction = "umap",
    label = TRUE,
    repel = TRUE) &
    theme(
    legend.position = c(0, 0),
    legend.justification = c("left", "bottom"),
    legend.box.just = "left")

    return(seurat_object)

}
