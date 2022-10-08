#' Annotate mouse cell type
#'
#' Uses the [`SingleR::SingleR()`] function to label cell types according
#' to expression profiles from the [Immunological Genome Project](https://www.immgen.org/) 
#' made available via the the `{celldex}` package. The Seurat object per-cell metadata
#' is updated with `cell_type` and `cell_type_fine` columns.
#'
#' @name seurat_annotate_type
#' @param seurat_object a Seurat object
#' @return a seurat object
#' @source https://www.frontiersin.org/articles/10.3389/fcell.2021.622519/full
#' @author Denis O'Meally, Yu-Husan Fu
#' @export
seurat_annotate_type <- function(seurat_object) {
    #seurat_object <- get_pin("mmu_10x_2022_1_GENCODEm28_HLT.rds")

    # Setup parallel processing
    ncores <- parallelly::availableCores()
    BiocParallel::register(
        BiocParallel::MulticoreParam(ncores, ncores * 2, progressbar = TRUE)
    )

    Seurat::DefaultAssay(seurat_object) <- "RNA"

    mmu_ImmGenData <- celldex::ImmGenData()

    mmu_ImmGenData_label_fine <- mmu_ImmGenData$label.fine

    cell_labels <- SingleR::SingleR(
        Seurat::GetAssayData(seurat_object, assay = "RNA", slot = "data"),
        mmu_ImmGenData,
        mmu_ImmGenData_label_fine,
        BPPARAM = BiocParallel::bpparam()
    )

    cell_types <- cell_labels |>
        as.data.frame() |>
        dplyr::mutate(
            cell_type_fine = pruned.labels,
            cell_type = stringr::str_replace(pruned.labels, " \\s*\\([^\\)]+\\)", "")) |>
        dplyr::select(
            cell_type,
            cell_type_fine)

    seurat_object <- Seurat::AddMetaData(seurat_object, cell_types)

    return(seurat_object)

}
