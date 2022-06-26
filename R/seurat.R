# Functions for dealing with scRNAseq data in the Surat package.

#' Rotate UMAP
#'
#' Makes the 1st cell in a Seurat object plot in the top left quadrant (ie -1,-1)
#' This is useful for plotting UMAPs in a consistent way across platforms/runs.
#' The "x" and "y" parameters rotate the axis by 90 degrees.
#'
#' @name rotate_umap
#' @param integrated a Seurat object
#' @param x flip the x-axis LOGICAL (TRUE or FALSE)
#' @param y flip the y-axis LOGICAL (TRUE or FALSE)
#' @return A Seurat object
#' @author Denis O'Meally
#' @export
rotate_umap <- function(integrated, x = FALSE, y = FALSE) {
    umap_coord <- (Seurat::Embeddings(integrated[["umap"]]))

    # check X coords
    if (umap_coord[1, 1] > 0) {
        umap_coord[, 1] <- umap_coord[, 1] * -1
    }
    if (x) {
        umap_coord[, 1] <- umap_coord[, 1] * -1
    }

    # check Y coords
    if (umap_coord[1, 2] > 0) {
        umap_coord[, 2] <- umap_coord[, 2] * -1
    }
    if (y) {
        umap_coord[, 2] <- umap_coord[, 2] * -1
    }

    integrated@reductions$umap <- Seurat::CreateDimReducObject(
        embeddings = umap_coord,
        assay = "RNA"
    )
    return(integrated)
}

#' @title Make a Seurat object from 10X Cellranger CellPlex data
#' @name make_seurat_objects
#' @description
#' Loads 10X Cellranger CellPlex h5 data, builds Seurat objects and adds the following metatdata columns:
#' `percent_mt`, `percent_ribo`, `percent_hb`, `percent_platelet`, `percent_xist`, `chrY_counts`, `percent_myh11`.
#'
#' Data are read in from the "cellranger" directory `/net/isi-dcnl/ifs/user_data/rrockne/MHO/AML.scRNA.2022/cellranger`
#' by searching for files named `sample_feature_bc_matrix.h5`. This can take some time, so after the 1st run, 
#' the paths are written to `data-raw/cellranger_h5_paths.txt` and read from there for subsequent runs.
#' Similarly, sex chromosome genes are parsed from the GTF and cached in `inst/extdata`.
#'
#' TODO: Include eg for accessing Y genes from package 
#'
#' @param path_regex string, optional, regex pattern to match in path. Use this to select different
#' genome builds, eg "GRCm38".
#' @return a list of Seurat objects
#' @author Denis O'Meally

# Make Seurat objects
seurat_import_objects <- function(path_regex, cellranger_folder = "/net/isi-dcnl/ifs/user_data/rrockne/MHO/AML.scRNA.2022/cellranger") {
    requireNamespace("hdf5r")
    # future::plan("multisession")
    # options(future.seed = TRUE)
    # options(future.globals.maxSize = 10 * 1024^3) # 10GB

    # get a list of ChrY genes from GTF file
    # Parse GTF: https://www.biostars.org/p/140471/

    if (file.exists("inst/extdata/mmu_chrY_genes.txt")) {
        chrY_genes <- scan("inst/extdata/mmu_chrY_genes.txt", character())
    } else {
        gtf <- glue::glue("{cellranger_folder}/ref/GENCODEm28_human_genes.filtered.gtf")
        gtf_genes <- rtracklayer::import(gtf) %>%
            as.data.frame() %>%
            unique()
        chrY_genes <- gtf_genes$gene_name[grepl("^chrY|^Y", gtf_genes$seqnames)] %>% unique()
        chrX_genes <- gtf_genes$gene_name[grepl("^chrX|^X", gtf_genes$seqnames)] %>% unique()
        # Remove PAR genes
        chrY_genes <- chrY_genes[!(chrY_genes %in% chrX_genes)]
        write.table(chrY_genes, "inst/extdata/mmu_chrY_genes.txt",
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE
        )
    }

    # find paths of sparse matrices

    if (file.exists("data-raw/cellranger_h5_paths.txt")) {
        h5_paths <- scan("data-raw/cellranger_h5_paths.txt", character())
    } else {
        # list.files(path = cellranger_folder,
        # pattern = "sample_feature_bc_matrix.h5", full.names = TRUE, recursive = TRUE)
        h5_paths <- system(
            glue::glue("find {cellranger_folder} -name sample_feature_bc_matrix.h5"),
            intern = TRUE
        )
        writeLines(h5_paths, "inst/extdata/cellranger_h5_paths.txt")
    }

    h5_paths <- grep(path_regex, h5_paths, value = TRUE)

    # Load all the samples
    seurat_object_list <- future.apply::future_lapply(X = h5_paths, FUN = function(x) {

        # extract metadata
        sample_name <- stringr::str_extract(x, "(?<=per_sample_outs\\/)(.*?)(?=\\/)")

        ref_genome <- gsub("_lib.*", "", stringr::str_extract(x, "(?<=cellranger\\/)(.*?)(?=\\/)"))
        tissue <- dplyr::case_when(
            stringr::str_detect(x, "PB") ~ "PBMC",
            stringr::str_detect(x, "BM") ~ "BM",
            TRUE ~ NA_character_
        )
        ckit <- dplyr::case_when(
            stringr::str_detect(x, "ckit") ~ TRUE,
            TRUE ~ FALSE
        )

        message(glue::glue("Loading sample: {sample_name}"))

        # assemble Seurat object
        x <- Seurat::CreateSeuratObject(
            counts = Seurat::Read10X_h5(x)$`Gene Expression`,
            project = sample_name,
            min.cells = 3,
            min.features = 200
        )

        # add metadata
        x <- Seurat::AddMetaData(x,
            metadata = c(ref_genome, tissue, ckit),
            col.name = c("ref_genome", "tissue", "ckit")
        )

        # add QC data
        x <- Seurat::AddMetaData(x,
            metadata = data.frame(
                Seurat::PercentageFeatureSet(x, pattern = "^mt-"),
                Seurat::PercentageFeatureSet(x, pattern = "^Rp[sl]"),
                Seurat::PercentageFeatureSet(x, pattern = "^Hb[^(p)]"),
                Seurat::PercentageFeatureSet(x, pattern = "^Pecam1|^Pf4"),
                Seurat::PercentageFeatureSet(x, pattern = "^Xist"),
                # some logic to manage the edge case where none or one Y gene(s) has counts resulting
                # in a zero length character or a vector rather than a matrix, as required by Matrix::colSums
                # colSums is required because PercentageFeatureSet() fails when supplied a full list of Y genes
                # Maybe https://github.com/satijalab/seurat/issues/1665
                ifelse(
                    sum(rownames(x@assays$RNA@counts) %in% chrY_genes) > 2,
                    (Matrix::colSums(x@assays$RNA@counts[rownames(x@assays$RNA@counts) %in% chrY_genes, ]) / Matrix::colSums(x@assays$RNA@counts)) %>% as.data.frame() * 100,
                    ifelse(
                        sum(rownames(x@assays$RNA@counts) %in% chrY_genes) == 1,
                        Seurat::PercentageFeatureSet(x, pattern = chrY_genes[chrY_genes %in% rownames(x@assays$RNA@counts)]),
                        Seurat::PercentageFeatureSet(x, pattern = "ZERO EXPRESSION")
                    )
                ),
                Seurat::PercentageFeatureSet(x, pattern = "^HSA-MYH11-gene")
            ),
            col.name = c(
                "percent_mt", "percent_ribo", "percent_hb", "percent_platelet",
                "percent_xist", "chrY_counts", "percent_myh11"
            )
        )
    })

    return(seurat_object_list)
}

#' Basic single cell QC
#'
#' Applies the following filters to the data:
#' 1. Removes genes expressed in fewer than 3 cells
#' 2. Removes cells with fewer than 200 genes expressed
#' 3. Removes cells where mitochondrial genes account for more than 20% of expressed transcripts
#' 4. Removes cells where ribosomal genes account for fewer than 5% of expressed transcripts
#'
#'
#' @name seurat_perform_cell_qc
#' @param raw_seurat_objects a list of Seurat objects, fresh in from Cellranger
#' @return a list of Seurat objects
#' @author Denis O'Meally
#' @export
seurat_perform_cell_qc <- function(raw_seurat_objects) {

    # loop through objects in list
    filtered_seurat_objects <- future.apply::future_lapply(
        X = raw_seurat_objects, FUN = function(x) {

            # Detection-based filtering
            # Genes expressed in at least 3 cells
            x <- subset(x, features = rownames(x)[Matrix::rowSums(x) > 3])

            # filter on mt and ribo reads, and cells with at least 200 genes ----------
            x <- subset(x, subset = nFeature_RNA > 200 & percent_mt < 20 & percent_ribo > 5)
        }, future.seed = TRUE
    )
}