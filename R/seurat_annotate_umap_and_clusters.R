#' Annotate UMAP and knn clusters
#'
#' Runs the 'standard workflow' using 2000 features and 50 PCs to cluster the data.
#'
#' @name seurat_annotate_umap_and_clusters
#' @param seurat_object a Seurat object
#' @return a Seurat object with UMAP and knn clusters
#' @author Denis O'Meally
#' @export
seurat_annotate_umap_and_clusters <- function(seurat_object) {
    # future::plan("multisession")
    # options(future.globals.maxSize = 10 * 1024^3) # 10GB
    # options(future.seed = TRUE)

    # Run the standard workflow for visualization and clustering -------------------
    nfeatures <- 2000
    dims <- 1:50
    cluster_res <- 0.5
    npcs <- 20
    seed_use <- 24
    seurat_object <- seurat_object %>%
        Seurat::FindVariableFeatures(selection.method = "vst", nfeatures = nfeatures) %>%
        Seurat::ScaleData() %>%
        Seurat::RunPCA(npcs = npcs) %>%
        Seurat::FindNeighbors() %>%
        Seurat::FindClusters(resolution = cluster_res) %>%
        Seurat::RunUMAP(
            reduction = "pca",
            dims = dims,
            seed.use = seed_use
        ) %>%
        rotate_umap()

    # seurat_object <- JackStraw(seurat_object, num.replicate = 100)
    # seurat_object <- ScoreJackStraw(seurat_object, dims = 1:15)
    # JackStrawPlot(seurat_object, dims = 1:15)
    # DimHeatmap(seurat_object, dims = 1:6, cells = 500, balanced = TRUE)
    # DimHeatmap(seurat_object, dims = 7:12, cells = 500, balanced = TRUE)
    # DimHeatmap(seurat_object, dims = 13:18, cells = 500, balanced = TRUE)


    npcs <- Seurat::find_elbow(seurat_object)

    seurat_object <- seurat_object %>%
        Seurat::RunPCA(npcs = npcs) %>%
        Seurat::FindNeighbors() %>%
        Seurat::FindClusters(resolution = cluster_res) %>%
        Seurat::RunUMAP(
            reduction = "pca",
            dims = 1:npcs,
            seed.use = seed_use
        ) %>%
        rotate_umap()

    # set default Assay & identity -----------------------------------------------------
    Seurat::DefaultAssay(integrated) <- "RNA"
    Seurat::Idents(integrated) <- "identity"

    # Normalise and scale  -----------------------------------------------------------
    integrated <- integrated %>%
        Seurat::NormalizeData() %>%
        Seurat::ScaleData()

    return(integrated)

}
