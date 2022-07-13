# Functions for publishing data to the Heamdata pin_board

# Setup the pin_board to use
# Post data to OneDrive pin board (TRUE) or Isilon pin board (FALSE)
onedrive <- TRUE
# Location for posting data files, if not OneDrive
haemdata_folder <- "/net/isi-dcnl/ifs/user_data/rrockne/MHO/haemdata"

#' Publish a SummarisedExperiment
#'
#' Publishes a SummarisedExperiment to the Haemdata Teams channel and to the `MHO/haemdata` folder
#' on the COH Isilon storage (in `.rda` format) along with an expression matrix
#' (in `.csv` format). Files are named according to the SummarisedExperiment's
#' `metadata$oject_name` variable
#'
#' @name publish_se
#' @param summarised_experiment a SummarisedExperiment ready to publish
#' @return the path to the published object
#' @author Denis O'Meally
#' @export
publish_se <- function(summarised_experiment) {
    name <- summarised_experiment@metadata$object_name

    csv_pin <- write_se2tpm_pin(summarised_experiment)
    csv_version <- get_latest_pin_version(csv_pin)

    # Publish the SummarisedExperiment rds file
    rds_pin <- write_se_pin(summarised_experiment)
    rds_version <- get_latest_pin_version(rds_pin)

    # Publish a note to Teams
    if (teams == TRUE) {
        post_teams_chat(glue::glue(
            "A new <em>SummarisedExperiment</em> and expression matrix have been posted to Teams: <br>
        {name}.rds<br>
        {name}.csv"
        ))
    }

    return(data.frame("pin_name" = c(csv_pin, rds_pin), "version" = c(csv_version, rds_version)))
}

#' Publish sample metadata
#'
#' Publishes a metadata table to the `MHO/haemdata` folder on Isilon storage in `.csv` format.
#' If `teams` is `TRUE`, a note is posted to the Teams channel.
#'
#' @name publish_metadata
#' @param metadata a data.frame ready to publish
#' @return a data.frame of the pins and version numbers for the published object
#' @author Denis O'Meally
#' @export
publish_metadata <- function(metadata) {

    # extract name
    name <- deparse(substitute(metadata))

    # make description
    description <- paste0(
        "A table describing sample metadata. See http://cgt.coh.org/haemdata/reference/", name, ".html for more information."
    )

    csv_pin <- pins::pin_write(
        pin_board,
        metadata,
        name = glue::glue("{name}.csv"),
        type = "csv",
        title = name,
        description = description
    )

    csv_version <- get_latest_pin_version(csv_pin)

    # Publish a note to Teams
    if (teams == TRUE) {
        post_teams_chat(glue::glue("A new metadata table has been posted to MHO/haemdata on Isilon:<br>
        {name}.csv"))
    }

    return(data.frame("pin_name" = c(csv_pin), "version" = c(csv_version)))
}

# set up pin_board
if (onedrive == TRUE) {
    # # OneDrive pin_board
    pin_board <- pins::board_ms365(
        drive = Microsoft365R::get_team("PSON AML State-Transition")$get_drive(),
        path = "haemdata",
        versioned = TRUE
    )
} else {
    pin_board <- pins::board_folder(
        haemdata_folder,
        versioned = TRUE
    )
}

# send a message to the Haemdata Teams channel
post_teams_chat <- function(message) {
    team <- Microsoft365R::get_team("PSON AML State-Transition")
    channel <- team$get_channel("haemdata")
    channel$send_message(message, content_type = "html")
}

# get the latest version of a pin
get_latest_pin_version <- function(pin_name) {
    pin_board |>
        pins::pin_versions(pin_name) |>
        dplyr::pull(version) |>
        dplyr::first()
}

# get names of all heamdata pins
get_pin_list <- function() {
    pin_board |>
        pins::pin_list()
}

# write a SummarisedExperiment pin, with name, description, and metadata
write_se_pin <- function(summarised_experiment) {
    name <- summarised_experiment@metadata$object_name
    description <- paste0(
        "A SummarisedExperiment object. See http://cgt.coh.org/haemdata/reference/", name, ".html for more information."
    )
    pin_board |>
        pins::pin_write(
            summarised_experiment,
            name = glue::glue("{name}.rds"),
            type = "rds",
            title = name,
            description = description,
            metadata = as.list(summarised_experiment@metadata)
        )
}
# write an expression matrix pin, with name, description, and metadata from the summarised experment
# calls make_tpm_matrix() with defaults for the expression matrix
write_se2tpm_pin <- function(summarised_experiment) {
    name <- summarised_experiment@metadata$object_name
    description <- paste0(
        "An expression matrix of genes expressed > 1 TPM in > 5 samples. See http://cgt.coh.org/haemdata/reference/", name, ".html for more information."
    )
    expn_mat <- make_tpm_matrix(summarised_experiment)
    pin_board |>
        pins::pin_write(
            expn_mat$tpm_matrix,
            name = glue::glue("{expn_mat$name}.csv"),
            type = "csv",
            title = name,
            description = description,
            metadata = as.list(summarised_experiment@metadata)
        )
}

# Pin a Seurat object as an `h5ad` file, with name, description, and metadata
#' @import SeuratDisk SeuratObject
#'
write_seurat_h5ad_pin <- function(seurat_object) {

    tmp <- tempdir()

    name <- ifelse(
        stringr::str_detect(seurat_object$ref_genome[1], "GENCODEm28"),
            "mmu_10x_2022_1_GENCODEm28_HLT_seurat",
            "mmu_10x_2022_1_GRCm38_HLT_seurat"
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
                name = glue::glue("{name}.h5ad"),
                title = name,
                description = description,
                metadata = metadata
            )

    system(glue::glue("rm {tmp}/{name}.*"))
    return(h5ad_pin)
}
