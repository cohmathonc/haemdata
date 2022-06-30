#' Annotate UMAP and knn clusters
#'
#' Runs the 'standard workflow' using 2000 features, but differs in that the number of PCs used for clustering is
#' is first determined using using the [`find_elbow()`] function. The number of clusters is
#' tuned using the `cluster_res` parameter within the function.
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

    # seurat_object <- qs::qread(glue::glue("{nf_core_cache}/tmp/GENCODEm28_HLT_seurat_qc_rpca.qs"), nthreads = future::availableCores())

    # Run the standard workflow for visualization and clustering -------------------
    nfeatures <- 2000
    cluster_res <- 0.5
    seed_use <- 24

    npcs <- find_elbow(seurat_object)

    seurat_object <- seurat_object |>
        Seurat::ScaleData() |>
        Seurat::RunPCA(npcs = npcs) |>
        Seurat::FindNeighbors() |>
        Seurat::FindClusters(resolution = cluster_res) |>
        Seurat::RunUMAP(
            reduction = "pca",
            dims = 1:npcs,
            seed.use = seed_use
        ) |>
        rotate_umap()

    # set default identity -----------------------------------------------------
    Seurat::DefaultAssay(seurat_object) <- "integrated"
    Seurat::Idents(seurat_object) <- "seurat_clusters"

    return(seurat_object)

}
