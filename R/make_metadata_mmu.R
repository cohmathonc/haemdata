#' Make metadata for mouse AML samples
#'
#' This function prepares the `metadata_mmu` object for all RNAseq libraries from AML and CML mice.
#' Minimal metadata fields include sample, fastq_1, fastq_2, strandedness,
#' mouse_id, tissue, week, timepoint, batch, treatment, genotype, sex, dob, project.
#'
#' Raise an \href{https://github.com/drejom/haemdata/issues}{issue on GitHub}
#' to report erroneous or missing records.
#' @details `weeks` are read in from `data-raw/timepoints.xlsx` and `left_joined` to the
#' `sample_sheet_all_mice` by `timepoint_project`. For samples with the time point `L` (leukemia) or `END`,
#' `weeks` is set to 1 + the penultimate sample for a mouse. Bone marrow samples are assigned the
#' `timepoint` `NA` and `weeks` is also set to 1 + the penultimate sample for that mouse.
#' @name make_metadata_mmu
#' @param sample_sheet_all_mice a data.frame produced by row binding all mouse sample sheets
#' derived from [R/parse_metadata.R](https://github.com/drejom/haemdata/blob/main/R/parse_metadata.R).
#' @return a data.frame
#' @author Denis O'Meally

make_metadata_mmu <- function(sample_sheet_all_mice) {
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

    metadata <- dplyr::left_join(sample_sheet, weeks, by = "timepoint_project") |>
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
            TRUE ~ weeks
        ))
    # TODO fox weeks for BM samples

    return(metadata)
}