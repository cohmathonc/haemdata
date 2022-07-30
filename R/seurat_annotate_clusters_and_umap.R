#' Annotate UMAP and knn clusters
#'
#' Runs the workflow as described in the
#' [sctransform vignette](https://satijalab.org/seurat/articles/sctransform_vignette.html)
#' post `sctransform`. The number of clusters is tuned using the `cluster_res` parameter 
#' within the [`Seurat::FindClusters()`] function, and the chosen resolution positioned 
#' last in the dplyr pipe so as to set `seurat_identity` accordingly.
#'
#' https://satijalab.org/seurat/articles/sctransform_vignette.html
#' https://satijalab.org/seurat/articles/integration_introduction.html#perform-an-integrated-analysis-1
#'
#' @name seurat_annotate_umap_and_clusters
#' @param seurat_object a Seurat object
#' @return a Seurat object with UMAP and knn clusters annotated
#' @author Denis O'Meally
#' @export
seurat_annotate_clusters_and_umap <- function(seurat_object) {

    # TODO make a parameter to iterate over cluster resolution values, and one to set it explicitly

    # Run the standard workflow for visualization and clustering -------------------

    seurat_object <- seurat_object |>
        Seurat::RunPCA() |>
        Seurat::FindNeighbors() |>
        Seurat::FindClusters(resolution = 0.6) |>
        Seurat::FindClusters(resolution = 0.8) |> # 30 clusters
        Seurat::FindClusters(resolution = 1.0) |>
        Seurat::FindClusters(resolution = 1.2) |> # 45 clusters
        Seurat::FindClusters(resolution = 0.4) |> # 22 clusters for GENCODEm28 integrated
        Seurat::RunUMAP(
            reduction = "pca",
            dims = 1:30,
            seed.use = 1448145
        ) |>
        rotate_umap()

    # Explicitly set the default assay & identity -------------------------
    Seurat::DefaultAssay(seurat_object) <- "SCT" 
    Seurat::Idents(seurat_object) <- "seurat_clusters"

    return(seurat_object)
}
