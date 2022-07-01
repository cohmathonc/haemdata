# Pin a Seurat object as an `h5ad` file, with name, description, and metadata
write_seurat_h5ad_pin <- function(seurat_object) {
    requireNamespace(SeuratDisk)
    requireNamespace(SeuratObject)

    tmp <- tempdir()

    name <- ifelse(
        stringr::str_detect(seurat_object$ref_genome[1], "GENCODEm28"),
            "mmu_aml_seurat_GENCODEm28_HLT",
            "mmu_aml_seurat_GRCm38_HLT"
    )

    description <- paste0(
        "A Seurat object exported to h5ad format. See http://cgt.coh.org/haemdata/articles/scRNAseq.html for more information."
    )
    metadata <- list(
        "cell_metadata" = list(seurat_object@meta.data |> names())
    )

    SeuratDisk::SaveH5Seurat(
        seurat_object,
        filename = glue::glue("{tmp}/{name}.h5seurat"),
        verbose = TRUE,
        overwrite = TRUE
    )

    SeuratDisk::Convert(glue::glue("{tmp}/{name}.h5seurat"), dest = "h5ad", overwrite = TRUE)

    h5ad_pin <- pin_board |>
            pins::pin_upload(
                paths = glue::glue("{tmp}/{name}.h5ad"),
                name = name,
                title = name,
                description = description,
                metadata = metadata
            )


    return(h5ad_pin) # return a pin URI
}
