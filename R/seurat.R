# Functions for dealing with scRNAseq data in the Surat package.

#' Find elbow in PCA for Seurat objects.
#'
#' Inspired by
#' https://hbctraining.github.io/scRNA-seq/lessons/elbow_plot_metric.html
#' Determine percent of variation associated with each PC,
#' Calculate cumulative percent variation for each PC,
#' Determine which PC exhibits cumulative variation greater than 90%,
#' and % variation associated with the PC is less than 5%, then
#' Choose the minimum of the two
#'
#' Optionally makes a plot of the cumulative percent variation.
#'
#' @title find_elbow
#' @name find_elbow
#' @param seurat_object a Seurat object
#' @param plot TRUE/FALSE to make a plot of the elbow, Default = FALSE
#' @return either an `integer` or a `ggplot2` object.
#' @author Denis O'Meally
#' @export
find_elbow <- function(seurat_object, plot = FALSE) {
    # Determine percent of variation associated with each PC
    pct <- seurat_object[["pca"]]@stdev / sum(seurat_object[["pca"]]@stdev) * 100

    # Calculate cumulative percent variation for each PC
    cumu <- cumsum(pct)

    # Determine which PC exhibits cumulative variation greater than 90%,
    # and % variation associated with the PC is less than 5%
    co1 <- which(cumu > 90 & pct < 5)[1]

    # Determine the difference between variation of PC and subsequent PC
    # and get the last point where change of % of variation is more than 0.1%.
    co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1

    # Choose the minimum of the two
    pcs <- min(co1, co2)

    if (plot) {
        # Create a data.frame with values to plot
        plot_df <- data.frame(
            pct = pct,
            cumu = cumu,
            rank = 1:length(pct)
        )
        # Elbow plot to visualize
        elbow_plot <- ggplot2::ggplot(plot_df, ggplot2::aes(cumu, pct, label = rank, color = rank > pcs)) +
            ggplot2::geom_point() +
            ggplot2::geom_text(hjust = 1, vjust = -0.5) +
            geom_vline(xintercept = 90, color = "grey") +
            ggplot2::geom_hline(yintercept = min(pct[pct > 5]), color = "grey") +
            ggplot2::theme_bw() +
            ggplot2::labs(
                title = paste0("Elbow Plot of Principal Components\nChoose the first ", pcs, " PCs to retain for clustering"),
                x = "Cumulative variation explained (%)",
                y = "Percent variation explained (%)"
            )
        return(elbow_plot)
    }
    return(pcs)
}

#' Rotate UMAP
#'
#' Makes the 1st cell in a Seurat object plot in the top left quadrant (ie -1,-1)
#' This is useful for plotting UMAPs in a consistent way across platforms/runs.
#' The "x" and "y" parameters rotate the axis by 90 degrees.
#'
#' @name rotate_umap
#' @param seurat_object a Seurat object
#' @param x flip the x-axis LOGICAL (TRUE or FALSE)
#' @param y flip the y-axis LOGICAL (TRUE or FALSE)
#' @return A Seurat object
#' @author Denis O'Meally
#' @export
rotate_umap <- function(seurat_object, x = FALSE, y = FALSE) {
    umap_coord <- (Seurat::Embeddings(seurat_object[["umap"]]))

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

    seurat_object@reductions$umap <- Seurat::CreateDimReducObject(
        embeddings = umap_coord,
        assay = "RNA"
    )
    return(seurat_object)
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
#' @param cellranger_folder string, path to the cellranger folder. Default
#' is `/net/isi-dcnl/ifs/user_data/rrockne/MHO/AML.scRNA.2022/cellranger`
#' @return a list of Seurat objects
#' @author Denis O'Meally
#' @export
# Make Seurat objects
seurat_import_objects <- function(path_regex,
                                  cellranger_folder = "/net/isi-dcnl/ifs/user_data/rrockne/MHO/AML.scRNA.2022/cellranger") {
    # get a list of ChrY genes from GTF file
    # Parse GTF: https://www.biostars.org/p/140471/

    if (file.exists("inst/extdata/mmu_chrY_genes.txt")) {
        chrY_genes <- scan("inst/extdata/mmu_chrY_genes.txt", character())
    } else {
        gtf <- glue::glue("{cellranger_folder}/ref/GENCODEm28_human_genes.filtered.gtf")
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
#' Applies the following filters to each Seurat object in the list:
#' 1. Removes genes expressed in fewer than 3 cells
#' 2. Removes cells with fewer than 200 genes expressed
#' 3. Removes cells where mitochondrial genes account for more than 20% of expressed transcripts
#' 4. Removes cells where ribosomal genes account for fewer than 5% of expressed transcripts
#'
#' The resulting list of Seurat objects is returned as a merged Seurat object, ready for downstream analysis.
#'
#' @name seurat_perform_cell_qc
#' @param raw_seurat_objects a list of Seurat objects, fresh in from Cellranger
#' @return a Seurat object
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

    seurat_object <- merge(
        x = filtered_seurat_objects[[1]],
        y = filtered_seurat_objects[2:length(filtered_seurat_objects)]
    )

    # add Misc metadata to the seurat object
    name <- gsub(
        "(?<=GENCODEm28_HLT|GRCm38_HLT).*$",
        "",
        deparse(substitute(raw_seurat_objects)),
        perl = TRUE
    )
    description <- glue::glue(
        "A haemdata Seurat object created using the SCTransform v2 integration method. See {haemdata_env$package_url}/reference/{name}.html for more information."
    )
    Seurat::Misc(seurat_object, slot = "description") <- description
    Seurat::Misc(seurat_object, slot = "name") <- name

    return(seurat_object)
}

#' Annotate mouse cell cycle
#'
#' Uses the function [`Seurat::CellCycleScoring()`] to annotate the cells with their cell cycle `Phase`,
#' `S.score` and `G2M.score`. The Seurat function uses an internal list of human gene symbols from
#' Tirosh et al 2016.
#'
#' Here we substitute equivalent mouse symbols from `biomaRt`,
#' but the service was down at time of writing. In a quick hack, we simply
#' convert gene symbols by making them lowercase, with the 1st letter in caps.
#'
#' https://rdrr.io/cran/Seurat/man/CellCycleScoring.html
#' https://github.com/satijalab/seurat/issues/2493
#'
#' @name seurat_annotate_cell_cycle
#' @source https://www.science.org/doi/abs/10.1126/science.aad0501
#' @param seurat_object A Seurat object to annotate
#' @return A Surat object with cell cycle annotations added to cell `meta.data` (`S.Score`, `G2M.Score` & `Phase`)
#' @author Denis O'Meally
#' @export
seurat_annotate_cell_cycle <- function(seurat_object) {

    # Previously used this with success:
    # https://ucdavis-bioinformatics-training.github.io/2019-single-cell-RNA-sequencing-Workshop-UCD_UCSF/scrnaseq_analysis/scRNA_Workshop-PART1.html
    # But its a bit hacky and overly complex.

    # ############################################################################
    # ### This is the recommended way to convert to mouse ids, but ensembl is down at time of testing
    # # Basic function to convert human to mouse gene names
    # convert_human_gene_list <- function(x) {
    #     # https://www.r-bloggers.com/2016/10/converting-mouse-to-human-gene-names-with-biomart-package/

    #     human <- biomaRt::useMart("ensembl", dataset = "hsapiens_gene_ensembl")
    #     mouse <- biomaRt::useMart("ensembl", dataset = "mmusculus_gene_ensembl")

    #     genesV2 <- biomaRt::getLDS(
    #         attributes = c("hgnc_symbol"),
    #         filters = "hgnc_symbol",
    #         values = x,
    #         mart = human,
    #         attributesL = c("mgi_symbol"),
    #         martL = mouse,
    #         uniqueRows = T
    #     )

    #     humanx <- unique(genesV2[, 2])

    #     # Print the first 6 genes found to the screen
    #     print(head(humanx))
    #     return(humanx)
    # }

    # g2m.genes <- convert_human_gene_list(Seurat::cc.genes.updated.2019$g2m.genes)
    # s.genes <- convert_human_gene_list(Seurat::cc.genes.updated.2019$s.genes)

    # In the interim, just make them all look like mouse symbols
    g2m_genes <- gsub("(?<=\\b)([a-z])", "\\U\\1", tolower(Seurat::cc.genes.updated.2019$g2m.genes), perl = TRUE)
    s_genes <- gsub("(?<=\\b)([a-z])", "\\U\\1", tolower(Seurat::cc.genes.updated.2019$s.genes), perl = TRUE)

    seurat_object <- Seurat::CellCycleScoring(
        object = seurat_object,
        g2m.features = g2m_genes,
        s.features = s_genes
    )

    return(seurat_object)
}