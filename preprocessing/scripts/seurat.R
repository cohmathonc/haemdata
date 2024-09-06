#' @title Make a Seurat object from 10X Cell Ranger
#' @name seurat_import_cellranger
#' @description
#' Loads 10X Cellranger h5 data from Isilon, builds Seurat objects and adds the following metatdata columns:
#' `percent_mt`, `percent_ribo`, `percent_hb`, `percent_platelet`, `percent_xist`, `chrY_counts`, `percent_myh11`.
#'
#' Data are read in from the run folder identified by an nf-core/scrnaseq Cell Ranger multiqc report
#'
#' Sex chromosome genes are parsed from the GTF and cached in `inst/extdata`.
#'
#' TODO: Include eg for accessing Y genes from package
#'
#' @param multiqc_path string, required when cellplex = FALSE; full Isilon path to an nf-core/scrnaseq multiqc report
#' @return a list of Seurat objects
#' @author Denis O'Meally
#' @export
# Make Seurat objects
seurat_import_cellranger <- function(multiqc_path) {
  # get a list of ChrY genes from GTF file
  # Parse GTF: https://www.biostars.org/p/140471/

  if (file.exists("inst/extdata/mmu_chrY_genes.txt")) {
    chrY_genes <- scan("inst/extdata/mmu_chrY_genes.txt", character())
  } else {
    gtf <- ifelse(hprcc::get_cluster() == "apollo",
      "/ref_genome/igenomes/Mus_musculus/Cellranger/GENCODEm28_human_genes.filtered.gtf",
      stop("reference not available on gemini"))
    gtf_genes <- rtracklayer::import(gtf) |>
      as.data.frame() |>
      unique()
    chrY_genes <- gtf_genes$gene_name[grepl("^chrY|^Y", gtf_genes$seqnames)] |> unique()
    chrX_genes <- gtf_genes$gene_name[grepl("^chrX|^X", gtf_genes$seqnames)] |> unique()
    # Remove PAR genes
    chrY_genes <- chrY_genes[!(chrY_genes %in% chrX_genes)]
    write.table(chrY_genes, "inst/extdata/mmu_chrY_genes.txt",
      sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE
    )
  }

  run_path <- sub("multiqc/multiqc_report.html", "cellranger/count", multiqc_path)

  h5_paths <- list.files(run_path, pattern = "filtered_feature_bc_matrix.h5", full.names = TRUE, recursive = TRUE)

  # Load all the samples
  seurat_object_list <- future.apply::future_lapply(X = h5_paths, FUN = function(x) {

    # extract metadata
    sample_name <- stringr::str_extract(x, "(COHP_\\d{5})")

    ref_genome <- "GENCODEm28_HLT"

    message(glue::glue("Loading sample: {sample_name}"))

    # assemble Seurat object
    x <- Seurat::CreateSeuratObject(
      counts = Seurat::Read10X_h5(x),
      project = sample_name,
      min.cells = 3,
      min.features = 200
    )

    # add metadata
    x <- Seurat::AddMetaData(x,
      metadata = c(ref_genome),
      col.name = c("ref_genome")
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
          sum(rownames(x@assays$RNA$counts) %in% chrY_genes) > 2,
          (Matrix::colSums(x@assays$RNA$counts[rownames(x@assays$RNA$counts) %in% chrY_genes, ]) / Matrix::colSums(x@assays$RNA$counts)) %>% as.data.frame() * 100,
          ifelse(
            sum(rownames(x@assays$RNA$counts) %in% chrY_genes) == 1,
            Seurat::PercentageFeatureSet(x, pattern = chrY_genes[chrY_genes %in% rownames(x@assays$RNA$counts)]),
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


#' @title Make a Seurat object from 10X Cellranger CellPlex data
#' @name seurat_import_cellplex
#' @description
#' Loads 10X Cellranger CellPlex h5ad data, builds Seurat objects and adds the following metatdata columns:
#' `percent_mt`, `percent_ribo`, `percent_hb`, `percent_platelet`, `percent_xist`, `chrY_counts`, `percent_myh11`.
#'
#' Data are read in from file paths recorded in the metadata_mmu table, column `hdf5`.
#' the `cohort_regex` parameter is used to select the cohort from by filtering on the
#' `cohort` column of the metadata_mmu table.
#' Sex chromosome genes are parsed from the GTF and cached in `inst/extdata`.
#'
#' TODO: Include eg for accessing Y genes from package
#'
#' @param cohort_regex string, optional, regex pattern to match in `cohort`. Use this to select different
#' cohorts.
#' @param metadata_mmu_prepub data.frame, required; metadata table with columns `cohort` and `hdf5`
#' @return a list of Seurat objects
#' @author Denis O'Meally
#' @export
seurat_import_cellplex <- function(cohort_regex, metadata_mmu_prepub) {
  # get a list of ChrY genes from GTF file
  # Parse GTF: https://www.biostars.org/p/140471/

# TODO This is failing Aug 20204 - need to fix with absolute paths or use here::here()
  if (file.exists("inst/extdata/mmu_chrY_genes.txt")) {
    chrY_genes <- scan("inst/extdata/mmu_chrY_genes.txt", character())
  } else {
    gtf <- ifelse(hprcc::get_cluster() == "apollo",
      "/ref_genome/igenomes/Mus_musculus/Cellranger/GENCODEm28_human_genes.filtered.gtf",
      stop("reference not available on gemini"))
    gtf_genes <- rtracklayer::import(gtf) |>
      as.data.frame() |>
      unique()
    chrY_genes <- gtf_genes$gene_name[grepl("^chrY|^Y", gtf_genes$seqnames)] |> unique()
    chrX_genes <- gtf_genes$gene_name[grepl("^chrX|^X", gtf_genes$seqnames)] |> unique()
    # Remove PAR genes
    chrY_genes <- chrY_genes[!(chrY_genes %in% chrX_genes)]
    write.table(chrY_genes, "inst/extdata/mmu_chrY_genes.txt",
      sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE
    )
  }
  h5_paths <- metadata_mmu_prepub |>
    dplyr::filter(stringr::str_detect(cohort, {{ cohort_regex }})) |>
    dplyr::filter(stringr::str_detect(assay, "scRNA")) |>
    dplyr::pull(hdf5)
  message("loaded paths from pinboard")
  if (length(h5_paths) == 0) {
    stop("No matching samples have h5ad file paths in the metadata_mmu table; check the cohort_regex.")
  }

  # Load all the samples
  seurat_object_list <- future.apply::future_lapply(X = h5_paths, FUN = function(x) {

    # extract metadata
    sample_name <- stringr::str_extract(x, "(?<=per_sample_outs\\/)(.*?)(?=\\/)")

    ref_genome <- gsub("_lib.*", "", stringr::str_extract(x, "(?<=cellranger\\/)(.*?)(?=\\/)"))
    tissue <- dplyr::case_when(
      stringr::str_detect(x, "PB_ckit") ~ "PBMC_CKIT",
      stringr::str_detect(x, "BM_ckit") ~ "BM_CKIT",
      stringr::str_detect(x, "BM") ~ "BM",
      stringr::str_detect(x, "PB") ~ "PBMC",
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
          sum(rownames(x@assays$RNA$counts) %in% chrY_genes) > 2,
          (Matrix::colSums(x@assays$RNA$counts[rownames(x@assays$RNA$counts) %in% chrY_genes, ]) / Matrix::colSums(x@assays$RNA$counts)) %>% as.data.frame() * 100,
          ifelse(
            sum(rownames(x@assays$RNA$counts) %in% chrY_genes) == 1,
            Seurat::PercentageFeatureSet(x, pattern = chrY_genes[chrY_genes %in% rownames(x@assays$RNA$counts)]),
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
