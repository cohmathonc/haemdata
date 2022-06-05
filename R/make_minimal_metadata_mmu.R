#' Make metadata for mouse AML samples
#'
#' This function prepares the `minimal_metadata_mmu` object for all RNAseq libraries from AML and CML mice.
#' Minimal metadata fields include sample, fastq_1, fastq_2, strandedness,
#' mouse_id, tissue, timepoint, batch, treatment, genotype, sex, dob, project.
#'
#' Raise an \href{https://github.com/drejom/haemdata/issues}{issue on GitHub}
#' to report erroneous or missing records.
#'
#' @name make_minimal_metadata_mmu
#' @param sample_sheet_all_mice a data.frame produced by row binding all mouse sample sheets
#' derived from [R/parse_metadata.R](https://github.com/drejom/haemdata/blob/main/R/parse_metadata.R).
#' @return a data.frame
#' @author Denis O'Meally

make_minimal_metadata_mmu <- function(sample_sheet_all_mice) {
    # consolidate sample metadata where possible
    sample_sheet <- sample_sheet_all_mice |>
        dplyr::group_by(mouse_id) |>
        tidyr::fill(c("genotype", "sex", "dob"), .direction = "updown") |>
        dplyr::ungroup() |>
        dplyr::mutate(timepoint_project = glue::glue("{timepoint}_{project}"))

    # Add weeks column
    xls <- "data-raw/timepoints.xlsx"
    weeks <- readxl::read_excel(xls, sheet = "weeks") |>
        dplyr::select(-c(timepoint, project))

    minimal_metadata <- dplyr::left_join(sample_sheet, weeks, by = "timepoint_project") |>
        dplyr::select(-c(timepoint_project)) |>
        dplyr::relocate(weeks, .before = "timepoint") |>
        dplyr::mutate(weeks = dplyr::case_when(
            mouse_id == "2684" & timepoint == "L" ~ 7,
            mouse_id == "2690" & timepoint == "L" ~ 4,
            mouse_id == "2708" & timepoint == "L" ~ 4,
            mouse_id == "2718" & timepoint == "L" ~ 6,
            mouse_id == "2719" & timepoint == "L" ~ 4,
            mouse_id == "2731" & timepoint == "L" ~ 3,
            mouse_id == "4309" & timepoint == "END" ~ 13,
            mouse_id == "4321" & timepoint == "END" ~ 7,
            mouse_id == "4324" & timepoint == "END" ~ 9,
            mouse_id == "4329" & timepoint == "END" ~ 13,
            mouse_id == "4319" & timepoint == "END" ~ 10,
            mouse_id == "4506" & timepoint == "END" ~ 10,
            stringr::str_detect(sample, "COHP_3849[123]") ~ 100,
            TRUE ~ weeks
        ))

    # write to the package data folder
    helpeRs::write_data("minimal_metadata_mmu", minimal_metadata)

    # write a CSV to extdata
    write.csv(minimal_metadata, file = "inst/extdata/minimal_metadata_samples.csv")

    # publish to Teams
    team <- Microsoft365R::get_team("PSON AML State-Transition")
    channel <- team$get_channel("haemdata")
    channel$send_message(paste0("Minimal metadata for all mice samples: minimal_metadata_samples.csv"),
        attachments = "inst/extdata/minimal_metadata_samples.csv"
    )

    return(minimal_metadata)
}