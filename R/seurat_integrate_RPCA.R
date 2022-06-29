#' RPCA integration for Seurat objects
#'
#' Integrate across orig.ident using RPCA.
#' https://satijalab.org/seurat/articles/integration_large_datasets.html
#'
#' @title integrate_RPCA
#' @param seurat_object_list a list of seurat objects
#' @return a seurat object
#' @author Denis O'Meally
#' @export
seurat_integrate_RPCA <- function(seurat_object_list) {
    requireNamespace("SeuratObject")
    # https://satijalab.org/seurat/articles/integration_large_datasets.html
    # Integrate across orig.ident
    # Reciprocal PCA
    options(future.globals.maxSize = 20 * 1024^3) # 20GB
    options(future.seed = TRUE)
    # seurat_object_list <- GENCODEm28_HLT_seurat_qc
    # seurat_object_list <- GRCm38_HLT_seurat_qc

    #Normalise and find variable features
    seurat_object_list <- future.apply::future_lapply(seurat_object_list, function(x) {
#    seurat_object_list <- lapply(seurat_object_list, function(x) {
        x <- Seurat::NormalizeData(x) |>
                Seurat::FindVariableFeatures()
    })
    # Get the most variable features across all samples in the list
    features <- Seurat::SelectIntegrationFeatures(seurat_object_list)
    # scale the normalised data and run PCA on variable features, using scaled data
    seurat_object_list <- future.apply::future_lapply(seurat_object_list, function(x) {
        #    seurat_object_list <- lapply(seurat_object_list, function(x) {
        x <- Seurat::ScaleData(x, features = features) |>
            Seurat::RunPCA(features = features)
    })
    #some mem management
    rm(features)
    # FindIntegrationAnchors
    anchors <- Seurat::FindIntegrationAnchors(
        seurat_object_list,
        reduction = "rpca",
        dims = 1:20
    )
    # Find a set of anchors between a list of Seurat objects.
    seurat_object_rpca <- Seurat::IntegrateData(anchorset = anchors, dims = 1:20) |>
        Seurat::ScaleData() |>
        Seurat::RunPCA() |>
        Seurat::RunUMAP(dims = 1:20)
    # some mem management
    rm(anchors)
#    qs::qsave(seurat_object_rpca, file = glue::glue("{nf_core_cache}/tmp/GENCODEm28_HLT_seurat_qc_rpca.qs", nthreads = 8L))
    qs::qsave(seurat_object_rpca, file = glue::glue("{nf_core_cache}/tmp/GRCm28_HLT_seurat_qc_rpca.qs", nthreads = 8L))
    
    return(seurat_object_rpca)

}
