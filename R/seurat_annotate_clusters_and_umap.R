#' Annotate UMAP and knn clusters
#'
#' Runs the 'standard workflow' using 2000 features, but differs in that the number of PCs used for clustering is
#' is first determined using using the [`find_elbow()`] function. The number of clusters is
#' tuned using the `cluster_res` parameter within the [`Seurat::FindClusters()`] function.
#'
#' https://satijalab.org/seurat/articles/integration_introduction.html#perform-an-integrated-analysis-1
#'
#' @name seurat_annotate_umap_and_clusters
#' @param seurat_object a Seurat object
#' @return a Seurat object with UMAP and knn clusters
#' @author Denis O'Meally
#' @export
seurat_annotate_clusters_and_umap <- function(seurat_object) {
    options(future.globals.maxSize = 10 * 1024^3) # 10GB

#TODO make a parameter to iterate over cluster resolution values, and one to set it explicitly

    # Run the standard workflow for visualization and clustering -------------------
    seed_use <- 24

    npcs <- find_elbow(seurat_object)

    seurat_object <- seurat_object |>
        Seurat::ScaleData() |>
        Seurat::RunPCA(npcs = npcs) |>
        Seurat::FindNeighbors() |>
        # Seurat::FindClusters(resolution = 0.4) |> # 22 clusters for GENCODEm28 integrated
        # Seurat::FindClusters(resolution = 0.6) |>
        Seurat::FindClusters(resolution = 0.8) |> # 30 clusters
        #Seurat::FindClusters(resolution = 1.0) |>
        #Seurat::FindClusters(resolution = 1.2) |> # 45 clusters
        Seurat::RunUMAP(
            reduction = "pca",
            dims = 1:npcs,
            seed.use = seed_use
        ) |>
        rotate_umap()

    # set clusters as the default identity -------------------------
    Seurat::DefaultAssay(seurat_object) <- "integrated" # setting explicitly
    Seurat::Idents(seurat_object) <- "seurat_clusters"

    return(seurat_object)

}
