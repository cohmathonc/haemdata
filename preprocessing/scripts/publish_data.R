# Functions for publishing data to the Heamdata pin_board --- ####
#' @include haemdata.R

# --- Publishing versioned objects --- ####

#' @title Publish a SummarisedExperiment
#' @description Publishes a SummarisedExperiment to the current pinboard
#' (in `.rda` format) along with an expression matrix
#' (in `.csv` format). Files are named according to the SummarisedExperiment's
#' `metadata$project_name` variable
#' @param summarised_experiment a SummarisedExperiment ready to publish
#' @return the path to the published object
#' @author Denis O'Meally
#' @rdname publish_se
#' @export
publish_se <- function(summarised_experiment) {
    # Write the TPM expression matrix csv file
    csv_pin <- write_se2tpm_pin(summarised_experiment)
    csv_version <- get_latest_pin_version(csv_pin)

    # Write the SummarisedExperiment rds file
    rds_pin <- write_se_pin(summarised_experiment)
    rds_version <- get_latest_pin_version(rds_pin)

    return(data.frame("pin_name" = c(csv_pin, rds_pin), "version" = c(csv_version, rds_version)))
}

#' Publish a Seurat object
#'
#' Publishes a Seurat object to the current pinboard (in `.rda` format)
#' along with a ScanPy object (in `.h5ad` format).
#' Files are named according to the Seurat object's `object_name` slot. (see )
#'
#' @name publish_seurat
#' @param seurat_object a Seurat object ready to publish
#' @return the path to the published object
#' @author Denis O'Meally
#' @export
publish_seurat <- function(seurat_object) {
    # Write the Seurat rds file
    rds_pin <- write_seurat_pin(seurat_object)
    rds_version <- get_latest_pin_version(rds_pin)

    # Write the ScanPy h5ad file
    h5ad_pin <- write_seurat_h5ad_pin(seurat_object)
    h5ad_version <- get_latest_pin_version(rds_pin)

    return(data.frame("pin_name" = c(rds_pin, h5ad_pin), "version" = c(rds_version, h5ad_version)))
}
#' Publish sample metadata
#'
#' Publishes a metadata table to the  Haemdata Teams channel or to the `MHO/haemdata`
#' folder on devel storage in `.csv` format.
#'
#' @name publish_metadata
#' @param metadata a data.frame ready to publish
#' @return a data.frame of the pins and version numbers for the published object
#' @author Denis O'Meally
#' @export
publish_metadata <- function(metadata) {
    if (is.null(haemdata::haemdata_env$pin_board)) {
        stop(haemdata::haemdata_env$pin_board_msg)
    } else {
        # extract name
        name <- deparse(substitute(metadata))

        # make description
        description <- glue::glue(
            "A table describing sample metadata. See {haemdata::haemdata_env$package_url}/reference/{name}.html for more information."
        )

        csv_pin <- pins::pin_write(
            haemdata::haemdata_env$pin_board,
            metadata,
            name = glue::glue("{name}.csv"),
            type = "csv",
            title = name,
            description = description
        )

        csv_version <- get_latest_pin_version(csv_pin)

        return(data.frame("pin_name" = c(csv_pin), "version" = c(csv_version)))
    }
}
#' Publish mirtop counts
#'
#' Publishes mirtopcounts table to the  Haemdata Teams channel or to the `MHO/haemdata`
#' folder on devel storage in `.csv` format.
#'
#' @name publish_mirtop_counts
#' @param mirtop_counts a data.frame ready to publish
#' @return a data.frame of the pins and version numbers for the published object
#' @author Denis O'Meally
#' @export
publish_mirtop_counts <- function(mirtop_counts) {
    if (is.null(haemdata::haemdata_env$pin_board)) {
        stop(haemdata::haemdata_env$pin_board_msg)
    } else {
        # extract name
        name <- deparse(substitute(mirtop_counts))

        # make description
        description <- glue::glue(
            "A table of miR counts from the nfcore/smrnaseq pipeline [script](https://github.com/nf-core/smrnaseq/blob/master/bin/collapse_mirtop.r) that collapses mirtop output to a matrix"
        )

        csv_pin <- pins::pin_write(
            haemdata::haemdata_env$pin_board,
            mirtop_counts,
            name = glue::glue("{name}.csv"),
            type = "csv",
            title = name,
            description = description
        )

        csv_version <- get_latest_pin_version(csv_pin)

        return(data.frame("pin_name" = c(csv_pin), "version" = c(csv_version)))
    }
}
# --- Writing pins to the pin_board --- ####

# write a SummarisedExperiment pin, with name, description, and metadata
write_se_pin <- function(summarised_experiment) {
    if (is.null(haemdata::haemdata_env$pin_board)) {
        stop(haemdata::haemdata_env$pin_board_msg)
    } else {
        name <- summarised_experiment@metadata$object_name

        description <- glue::glue(
            "A SummarisedExperiment object. See {haemdata::haemdata_env$package_url}/reference/{name}.html for more information."
        )

        haemdata::haemdata_env$pin_board |>
            pins::pin_write(
                summarised_experiment,
                name = glue::glue("{name}.rds"),
                type = "rds",
                title = name,
                description = description,
                metadata = as.list(summarised_experiment@metadata)
            )
    }
}

# write an expression matrix pin, with name, description, and metadata from the summarised experment
# calls make_tpm_matrix() with defaults for the expression matrix
write_se2tpm_pin <- function(summarised_experiment) {
    if (is.null(haemdata::haemdata_env$pin_board)) {
        stop(haemdata::haemdata_env$pin_board_msg)
    } else {
        name <- summarised_experiment@metadata$object_name

        description <- glue::glue(
            "An expression matrix of genes expressed > 1 TPM in > 5 samples. See {haemdata::haemdata_env$package_url}/reference/{name}.html for more information."
        )

        expn_mat <- haemdata::make_tpm_matrix(summarised_experiment)

        haemdata::haemdata_env$pin_board |>
            pins::pin_write(
                expn_mat$tpm_matrix,
                name = glue::glue("{expn_mat$name}.csv"),
                type = "csv",
                title = name,
                description = description,
                metadata = as.list(summarised_experiment@metadata)
            )
    }
}

# Pin a Seurat object as an `h5ad` file, with name, description, and metadata
#'
write_seurat_h5ad_pin <- function(seurat_object) {
    if (is.null(haemdata::haemdata_env$pin_board)) {
        stop(haemdata::haemdata_env$pin_board_msg)
    } else {
        tmp <- tempdir()

        name <- Seurat::Misc(seurat_object)[["name"]]

        description <- glue::glue(
            "A Seurat object exported to h5ad format. See {haemdata::haemdata_env$package_url}/reference/{name}.html for more information."
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

        h5ad_pin <- haemdata::haemdata_env$pin_board |>
            pins::pin_upload(
                paths = glue::glue("{tmp}/{name}.h5ad"),
                name = glue::glue("{name}.h5ad"),
                title = name,
                description = description,
                metadata = metadata
            )

        system(glue::glue("rm {tmp}/{name}.*"))

        return(h5ad_pin)
    }
}

# Pin a Seurat object, with name, description, and metadata
#'
write_seurat_pin <- function(seurat_object) {
    if (is.null(haemdata::haemdata_env$pin_board)) {
        stop(haemdata::haemdata_env$pin_board_msg)
    } else {
        name <- Seurat::Misc(seurat_object)[["name"]]

        description <- glue::glue(
            "A Seurat object saved in rds format. See {haemdata::haemdata_env$package_url}/reference/{name}.html for more information."
        )

        metadata <- list(
            "cell_metadata" = list(seurat_object@meta.data |> names())
        )

        seurat_pin <- haemdata::haemdata_env$pin_board |>
            pins::pin_write(
                seurat_object,
                name = glue::glue("{name}.rds"),
                type = "rds",
                title = name,
                description = description,
                metadata = metadata
            )

        return(seurat_pin)
    }
}