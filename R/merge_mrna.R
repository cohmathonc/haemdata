#' Merge a pair of `SummarisedExperiments`
#'
#' Takes a pair of `SummarisedExperiments` produced by \code{\link{get_rnaseq_se()}},
#' and merges them into a single `SummarisedExperiment`
#'
#' @title merge_mrna
#' @param summarised_experiment1, name of a `SummarisedExperiment` object to merge
#' @param summarised_experiment2, name of the other `SummarisedExperiment` object to merge
#' @param drop a vector of sample names to remove from the merged `SummarisedExperiment`
#' @return a `SummarisedExperiment` object
#' @author Denis O'Meally
#' @export
merge_mrna <- function(summarised_experiment1, summarised_experiment2, drop=NULL) {

    # metadata1 <- S4Vectors::metadata(cml_mrna_2021_m38)
    # metadata2 <- S4Vectors::metadata(cml_mrna_2022_m38)
    # drop <- c("F545.7", "X490.17")

    metadata1 <- S4Vectors::metadata(summarised_experiment1)
    metadata2 <- S4Vectors::metadata(summarised_experiment2)

    # Check the $reference_genome version matches
    identical(
        metadata1$reference_genome,
        metadata2$reference_genome
    )
    reference_genome <- metadata1$reference_genome

    # Check the pipeline version matches
#    identical(metadata1$pipeline, metadata2$pipeline)
    pipeline <- metadata1$pipeline

    # Check the project name matches
    identical(
        stringr::str_match(metadata1$project, "$*(.*_.*?)_")[, 2],
        stringr::str_match(metadata2$project, "$*(.*_.*?)_")[, 2]
    )
    project <- stringr::str_match(metadata1$project, "$*(.*_.*?)_")[, 2]

    # update the metadata
    metadata <- list(
        project = project,
        reference_genome = reference_genome,
        pipeline = pipeline,
        date = date()
    )

    # merge the datasets
    summarised_experiment <- cbind(summarised_experiment1, summarised_experiment2)
    S4Vectors::metadata(summarised_experiment) <- metadata

    # remove samples that failed QC
    summarised_experiment[, !(colnames(summarised_experiment) %in% drop)]

    # name the output
    out_name <- paste(tolower(project), reference_genome, sep = "_")

    # write to the package data folder
    write_data(out_name, summarised_experiment)

    # write a CSV to extdata
    tpm_matrix_csv <- write_tpm_matrix(summarised_experiment, tpm = 1, samples = 5)

    # publish to Teams
    team <- Microsoft365R::get_team("PSON AML State-Transition")
    channel <- team$get_channel("haemdata")
    channel$send_message(paste0("New dataset available for ", out_name, ": ", basename(tpm_matrix_csv)),
        attachments = tpm_matrix_csv)

    # PINS https://pins.rstudio.com/reference/board_ms365.html
    # # A board in a SharePoint Online document library
    # sp <- Microsoft365R::get_sharepoint_site("PSON AML State-Transition")
    # chan_folder <- chan$get_folder()
    # board <- pins::board_ms365(chan_folder, "haemdata")

    return(summarised_experiment)

}
