
#' scTransform integration for Seurat objects
#'
#' Integrate across orig.ident using scTransform.
#' https://satijalab.org/seurat/articles/integration_introduction.html#performing-integration-on-datasets-normalized-with-sctransform-1
#'
#' @title seurat_integrate_scTransform
#' @param seurat_object a Seurat object
#' @return a Seurat object
#' @author Denis O'Meally
#' @export
seurat_integrate_scTransform <- function(seurat_object_list) {
    # Integrate across orig.ident
    # scTransform
    options(future.globals.maxSize = 40 * 1024^3) # 40GB
    # seurat_object_list <- GENCODEm28_HLT_seurat_qc
    # seurat_object_list <- GRCm38_HLT_seurat_qc
    # seurat_object_list <- mmu_10x_2022_1_GENCODEm28_HLT_qc

    # Apply scTransform to each object in the list
    seurat_object_list <- future.apply::future_lapply(
        seurat_object_list,
        Seurat::SCTransform,
        future.seed = TRUE
    )
    # Select features for integration
    features <- Seurat::SelectIntegrationFeatures(
        seurat_object_list,
        nfeatures = 3000
    )
    # Prepare Objects for integration
    seurat_object_list <- Seurat::PrepSCTIntegration(
        seurat_object_list,
        anchor.features = features
    )
    # Find anchors
    anchors <- Seurat::FindIntegrationAnchors(
        seurat_object_list,
        normalization.method = "SCT",
        anchor.features = features
    )
    # Integrate & DimReduction
    seurat_object_scT <- IntegrateData(
        anchorset = anchors,
        normalization.method = "SCT") |>
        Seurat::RunPCA(seurat_object_scT, verbose = FALSE) |>
        Seurat::RunUMAP(seurat_object_scT,verbose = FALSE)



    # Find a set of anchors between a list of Seurat objects.
    seurat_object_rpca <- Seurat::IntegrateData(anchorset = anchors, dims = 1:50) |>
        Seurat::ScaleData() |>
        Seurat::RunPCA() |>
        Seurat::RunUMAP(dims = 1:50)
    # some mem management
    rm(anchors)

    # add Misc metadata to the seurat object
    name <- gsub(
        "(?<=GENCODEm28_HLT|GRCm38_HLT).*$",
        "",
        deparse(substitute(seurat_object_list)),
        perl = TRUE
    )
    description <- glue::glue(
        "A haemdata Seurat object created using the RPCA integration method. See {package_url}/reference/{name}.html for more information."
    )

    Seurat::Misc(seurat_object_rpca, slot = "description") <- description
    Seurat::Misc(seurat_object_rpca, slot = "name") <- name



    #! FIXME Work around....see _targets.R for details
    ref <- ifelse(stringr::str_detect(seurat_object_list[[1]]$ref_genome[1], "GENCODEm28"), "GENCODEm28", "GRCm38")
    qs::qsave(seurat_object_rpca, file = glue::glue("{nf_core_cache}/tmp/{ref}_HLT_seurat_qc_rpca.qs", nthreads = future::availableCores()))
    #!

    return(seurat_object_rpca)

}
