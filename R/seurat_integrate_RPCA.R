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
    #seurat_object_list <- GENCODEm28_HLT_seurat_qc

    #Normalise and find variable features
    seurat_object_list <- future.apply::future_lapply(seurat_object_list, function(x) {
        x <- Seurat::NormalizeData(x, verbose = FALSE)
        x <- Seurat::FindVariableFeatures(x, verbose = FALSE)
    })
    # Get the most variable features across all samples in the list
    features <- Seurat::SelectIntegrationFeatures(seurat_object_list)
    # scale the normalised data and run PCA on variable features, using scaled data
    seurat_object_list <- future.apply::future_lapply(seurat_object_list,  function(x) {
        x <- Seurat::ScaleData(x, features = features, verbose = FALSE)
        x <- Seurat::RunPCA(x, features = features, verbose = FALSE)
    })
    # FindIntegrationAnchors
    anchors <- Seurat::FindIntegrationAnchors(
        seurat_object_list,
        reduction = "rpca",
        #reference = c(1, 8, 14),
        dims = 1:20,

    )
    # Find a set of anchors between a list of Seurat objects.
    seurat_object_rpca <- Seurat::IntegrateData(anchorset = anchors, dims = 1:20)
    message(">>>DefaultAssay", date())
    seurat_object_rpca <- Seurat::DefaultAssay(seurat_object_rpca, "integrated")
    message(">>>ScaleData", date())
    seurat_object_rpca <- Seurat::ScaleData(seurat_object_rpca, verbose = FALSE)
    message(">>>RunPCA", date())
    seurat_object_rpca <- Seurat::RunPCA(seurat_object_rpca, verbose = FALSE)
    message(">>>RunUMAP", date())
    seurat_object_rpca <- Seurat::RunUMAP(seurat_object_rpca, dims = 1:50)

    #<<<<<
    message("<<<", date())
    # Restore output to console
    sink(type = c("warning", "message"))
    sink()
    #<<<<<
    return(seurat_object_rpca)

}
