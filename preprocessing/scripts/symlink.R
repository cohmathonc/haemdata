#' symlink_10x_fastqs
#'
#' Given a sample_sheet with library_id, fastq_1 and fastq_2, this function will symlink
#' files to the {nf_core_cache} in the `folder`, required when original file paths have spaces.
#'
#' @param folder
#' @param sample_sheet a data.frame with the columns library_id, fastq_1 and fastq_2
#' @return an nf-core sample sheet with paths of fastq symlinks, without spaces
symlink_10x_fastqs <- function(folder, sample_sheet) {
    # nf_core_cache <- "/labs/rrockne/MHO/haemdata-nf-core-cache"
    # sample_sheet <- cml_mir142_ko_sample_sheet
    # folder <- "mmu_10X_mir142_ko/fastqs"

    # add a lane number column
    sample_sheet <- sample_sheet |>
        arrange(fastq_1) |>
        group_by(library_id) |>
        mutate(lane = paste0("L00", row_number())) |>
        ungroup()

    output_folder <- glue::glue("{nf_core_cache}/{folder}")

    # remove the output folder if it exists
    unlink(output_folder, recursive = TRUE, force = TRUE)

    # Create the output folder
    dir.create(output_folder, recursive = TRUE, showWarnings = FALSE)

    # Use map_df to apply the symlinking logic to each row of the dataframe
    sample_sheet_out <- sample_sheet %>%
        purrr::pmap_dfr(function(library_id, fastq_1, fastq_2, lane) {

            # Create symlinks for fastq_1 and fastq_2
            new_fastq_1 <- file.path(output_folder, glue::glue("{library_id}_S1_{lane}_R1_001.fastq.gz"))
            new_fastq_2 <- file.path(output_folder, glue::glue("{library_id}_S1_{lane}_R2_001.fastq.gz"))

            if (!file.symlink(fastq_1, new_fastq_1)) stop("Error creating symlink for fastq_1 for", new_fastq_1)
            if (!file.symlink(fastq_2, new_fastq_2)) stop("Error creating symlink for fastq_2 for", new_fastq_2)

            # Return a data.frame with the new file paths
            data.frame(
                library_id = library_id,
                fastq_1 = new_fastq_1,
                fastq_2 = new_fastq_2
            )
        })
    return(sample_sheet_out)
}

#' symlink_mrna_fastqs
#'
#' Given a sample_sheet with library_id, fastq_1, fastq_2 and strandedness this function will symlink
#' files to the {nf_core_cache} in the `folder`, required when original file paths have spaces.
#'
#' @param folder
#' @param sample_sheet a data.frame with the columns library_id, fastq_1, fastq_2 and strandeness
#' @return an nf-core sample sheet with paths of fastq symlinks, without spaces
symlink_mrna_fastqs <- function(folder, sample_sheet) {
    # nf_core_cache <- "/labs/rrockne/MHO/haemdata-nf-core-cache"
    # sample_sheet <- sample_sheet_mrna_mir142ko
    # folder <- "mmu_mrna_mir142_ko/fastqs"

        output_folder <- glue::glue("{nf_core_cache}/{folder}")

    # remove the output folder if it exists
    unlink(output_folder, recursive = TRUE, force = TRUE)

    # Create the output folder
    dir.create(output_folder, recursive = TRUE, showWarnings = FALSE)

    # Use map_df to apply the symlinking logic to each row of the dataframe
    sample_sheet_out <- sample_sheet %>%
        purrr::pmap_dfr(function(library_id, fastq_1, fastq_2, strandedness) {
            # Create symlinks for fastq_1 and fastq_2
            new_fastq_1 <- file.path(output_folder, glue::glue("{library_id}_R1.fastq.gz"))
            new_fastq_2 <- file.path(output_folder, glue::glue("{library_id}_R2.fastq.gz"))

            if (!file.symlink(fastq_1, new_fastq_1)) stop("Error creating symlink for fastq_1 for", new_fastq_1)
            if (!file.symlink(fastq_2, new_fastq_2)) stop("Error creating symlink for fastq_2 for", new_fastq_2)

            # Return a data.frame with the new file paths
            data.frame(
                library_id = library_id,
                fastq_1 = new_fastq_1,
                fastq_2 = new_fastq_2,
                strandedness = strandedness
            )
        })
    return(sample_sheet_out)
}
