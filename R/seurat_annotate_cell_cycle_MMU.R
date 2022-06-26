
#' Annotate cell cycle
#'
#' Get the mouse cell cycle markers from the package `scran`, and annotate the
#' Seutat object with the cell cycle labels for each cell.
#' 
#' See <https://ucdavis-bioinformatics-training.github.io/2019-single-cell-RNA-sequencing-Workshop-UCD_UCSF/scrnaseq_analysis/scRNA_Workshop-PART1.html>
#' for more information.
#'
#' @name seurat_annotate_cell_cycle
#' @param integrated A Seurat object to annotate
#' @return A Surat object with cell cycle annotations
#' @author Denis O'Meally
#' @export
seurat_annotate_cell_cycle <- function(integrated) {
    requireNamespace("org.Mm.eg.db")

    # https://ucdavis-bioinformatics-training.github.io/2019-single-cell-RNA-sequencing-Workshop-UCD_UCSF/scrnaseq_analysis/scRNA_Workshop-PART1.html
    # Calculate cell cycle, add to meta data
    # Using the package scran, get the mouse cell cycle markers and a mapping of m

    # set DefaultAssay
    Seurat::DefaultAssay(integrated) <- "RNA"

    # Convert to matrix for use by cycle
    mat <- as.matrix(Seurat::GetAssayData(integrated))

    # Get cell cycle markers
    mm_pairs <- readRDS(system.file("exdata", "mouse_cycle_markers.rds", package = "scran"))

    #match EnsemblIDs to MGI symbols
    anno_data <- clusterProfiler::bitr(rownames(mat), fromType = "SYMBOL", toType = "ENSEMBL", OrgDb = "org.Mm.eg.db")
    row_order <- match(rownames(mat), anno_data$SYMBOL)
    rownames(mat) <- anno_data$ENSEMBL[row_order]

    #drop rows with no Ensembl ID
    drop <- which(is.na(rownames(mat)))
    mat <- mat[-drop, ]

    #calculate cell cycle
    cycles <- scran::cyclone(mat, pairs = mm_pairs, BPPARAM = BiocParallel::MulticoreParam(future::availableCores()))

    #add to meta data
    tmp <- data.frame(cell_cycle = cycles$phases)
    rownames(tmp) <- colnames(mat)
    integrated <- Seurat::AddMetaData(integrated, tmp)

    #Make a DimPLot
    # set default Assay & identity -----------------------------------------------------
    Seurat::DefaultAssay(integrated) <- "RNA"
    Seurat::Idents(integrated) <- "cell_cycle"

    return(integrated)
}
