#' scTransform integration for Seurat objects
#'
#' Integrate a Seurat object using `sctransform`.
#' https://satijalab.org/seurat/articles/integration_introduction.html#performing-integration-on-datasets-normalized-with-sctransform-1
#'
#' @title Integrate a Seurat object using `sctransform`
#' @param seurat_object a Seurat object
#' @return Returns a Seurat object with a new assay (named SCT by default) with counts
#' being (corrected) counts, data being log1p(counts), scale.data being pearson residuals;
#' `sctransform::vst` intermediate results are saved in misc slot of the new assay.
#' @author Denis O'Meally
#' @export
seurat_sctransform <- function(seurat_object) {
    # scTransform
    Seurat::SCTransform(
        seurat_object,
        #vars.to.regress = "percent.mt",
        method = "glmGamPoi",
        seed.use = 1448145,
        verbose = TRUE)
}
