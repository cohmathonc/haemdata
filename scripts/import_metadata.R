# Functions for importing metadata
# TODO Remove "dob" from the parse_metadata_* functions, as it's added
# TODO in make_metadata_mmu()
#### metadata consolidation -----

#' Make minimal metadata for HSA
#'
#' Collects minimal metadata for human samples from the "sample summary"
#' excel sheet for FLT3 samples from COH Biobank and for MDS patients from the
#' EGA ([EGAD00001003891](https://ega-archive.org/datasets/EGAD00001003891)).
#'
#' The resulting table is sorted by `batch`, `patient_id` and `date`.
#'
#' @name make_metadata_hsa
#' @param sample_sheet A data.frame with a sample sheet
#' @return a data.frame sorted by batch, patient_id, date
#' @author Denis O'Meally
#' @export
make_metadata_hsa <- function(sample_sheet) {
    metadata_hsa <- sample_sheet |>
        dplyr::arrange(batch, patient_id, date)

    return(metadata_hsa)
}

#' Make metadata for mouse AML samples
#'
#' This function prepares the `metadata_mmu` object for all RNAseq libraries from AML and CML mice.
#' Minimal metadata fields include sample, fastq_1, fastq_2, strandedness,
#' mouse_id, tissue, sample_date, week, timepoint, batch, treatment, genotype,
#' sex, dob, dod, project.
#'
#' Raise an \href{https://github.com/drejom/haemdata/issues}{issue on GitHub}
#' to report erroneous or missing records.
#' @details Sample dates are read in from Teams (General|AML.Seq.Samples_dates.xlsx). For samples with
#' the timepoint label `L` (leukemia) or `END`, `weeks` for the sample is
#' calculated from the latest date in the sample sheet for that mouse.
#' For bone marrow samples, `weeks` is set to the date of death (`dod`).
#' @name make_metadata_mmu
#' @param sample_sheet_all_mice a data.frame produced by row binding all mouse sample sheets
#' detailed in [R/import_metadata.R](https://github.com/drejom/haemdata/blob/main/R/import_metadata.R).
#' @return a data.frame
#' @author Denis O'Meally

make_metadata_mmu <- function(sample_sheet_all_mice) {
    all_mice <- sample_sheet_all_mice
    # Add in columns to update
    cols_to_update_character <- c(
        "mouse_id", "tissue", "timepoint", "project",
        "treatment", "genotype", "sex"
    )
    cols_to_update_date <- c(
        "sample_date", "dob", "dod"
    )
    all_mice[c(cols_to_update_character, cols_to_update_date)] <- NA
    all_mice[cols_to_update_character] <- lapply(all_mice[cols_to_update_character], as.character)
    all_mice[cols_to_update_date] <- lapply(all_mice[cols_to_update_date], as.Date)

    # AML sample and mouse metadata
    # download the file
    get_teams_file("General/Copy of matadata_mmu_pivoted_AMLmice.YK.xlsx")

    # read in the worksheets
    wide_dates <- readxl::read_excel(
        "data-raw/Copy of matadata_mmu_pivoted_AMLmice.YK.xlsx",
        col_names = TRUE,
        col_types = "text"
    ) |>
        dplyr::select(-c("Helper", "cohp"))

    missing_t6_t7_dates <- tidyr::pivot_longer(
        data = wide_dates,
        cols = !c("sample", "project", "mouse_id", "tissue"),
        names_to = "timepoint",
        values_to = "sample_date"
    ) |>
        tidyr::drop_na("sample_date") |>
        dplyr::mutate(sample_date = as.Date(sample_date)) |>
        dplyr::arrange(mouse_id, sample_date) |>
        dplyr::select("sample", "mouse_id", "tissue", "timepoint", "project", "sample_date")

    # T6 was left off the sheet sent to Ya-Huei, so we need to add it in separately
    t6_t7_dates <- tibble::tribble(
        ~sample, ~mouse_id, ~tissue, ~timepoint, ~project, ~sample_date,
        "COHP_11843", "2683", "PBMC", "T6", "AML.mRNA.2016", "2/23/2016",
        "COHP_11844", "2685", "PBMC", "T6", "AML.mRNA.2016", "2/23/2016",
        "COHP_11845", "2686", "PBMC", "T6", "AML.mRNA.2016", "2/23/2016",
        "COHP_11846", "2689", "PBMC", "T6", "AML.mRNA.2016", "2/23/2016",
        "COHP_11847", "2692", "PBMC", "T6", "AML.mRNA.2016", "2/23/2016",
        "COHP_11848", "2700", "PBMC", "T6", "AML.mRNA.2016", "2/23/2016",
        "COHP_11849", "2702", "PBMC", "T6", "AML.mRNA.2016", "2/23/2016",
        "COHP_11850", "2705", "PBMC", "T6", "AML.mRNA.2016", "2/23/2016",
        "COHP_11851", "2709", "PBMC", "T6", "AML.mRNA.2016", "2/23/2016",
        "COHP_11852", "2720", "PBMC", "T6", "AML.mRNA.2016", "2/23/2016",
        "COHP_20898", "3335", "PBMC", "T6", "AML.mRNA.2018.all_samples", "8/15/2017",
        "COHP_20978", "3336", "PBMC", "T6", "AML.mRNA.2018.all_samples", "8/15/2017",
        "COHP_20919", "3336", "PBMC", "T7", "AML.mRNA.2018.all_samples", "9/12/2017",
        "COHP_20922", "3341", "PBMC", "T6", "AML.mRNA.2018.all_samples", "8/15/2017",
        "COHP_20923", "3357", "PBMC", "T6", "AML.mRNA.2018.all_samples", "8/15/2017",
        "COHP_20882", "3357", "PBMC", "T7", "AML.mRNA.2018.all_samples", "9/12/2017",
        "COHP_20982", "3368", "PBMC", "T6", "AML.mRNA.2018.all_samples", "8/15/2017",
        "COHP_21001", "3368", "PBMC", "T7", "AML.mRNA.2018.all_samples", "9/12/2017",
        "COHP_20965", "3370", "PBMC", "T6", "AML.mRNA.2018.all_samples", "8/15/2017",
        "COHP_20985", "3370", "PBMC", "T7", "AML.mRNA.2018.all_samples", "9/12/2017",
        "COHP_21006", "3339", "PBMC", "T6", "AML.mRNA.2018.all_samples", "8/15/2017",
        "COHP_20986", "3339", "PBMC", "T7", "AML.mRNA.2018.all_samples", "9/12/2017",
        "COHP_20907", "3342", "PBMC", "T6", "AML.mRNA.2018.all_samples", "8/15/2017",
        "COHP_21008", "3342", "PBMC", "T7", "AML.mRNA.2018.all_samples", "9/12/2017",
        "COHP_20970", "3346", "PBMC", "T6", "AML.mRNA.2018.all_samples", "8/15/2017",
        "COHP_20910", "3346", "PBMC", "T7", "AML.mRNA.2018.all_samples", "9/12/2017",
        "COHP_20950", "3349", "PBMC", "T6", "AML.mRNA.2018.all_samples", "8/15/2017",
        "COHP_20971", "3349", "PBMC", "T7", "AML.mRNA.2018.all_samples", "9/12/2017",
        "COHP_20972", "3334", "PBMC", "T6", "AML.mRNA.2018.all_samples", "8/15/2017",
        "COHP_20993", "3335", "PBMC", "T7", "AML.mRNA.2018.all_samples", "9/12/2017",
        "COHP_20894", "3340", "PBMC", "T6", "AML.mRNA.2018.all_samples", "8/15/2017",
        "COHP_20915", "3340", "PBMC", "T7", "AML.mRNA.2018.all_samples", "9/12/2017",
        "COHP_20916", "3350", "PBMC", "T6", "AML.mRNA.2018.all_samples", "8/15/2017",
        "COHP_20956", "3350", "PBMC", "T7", "AML.mRNA.2018.all_samples", "9/12/2017",
        "COHP_41243", "3694", "PBMC", "T6", "AML.mRNA.2020", "5/18/2018",
        "COHP_41244", "3694", "PBMC", "T7", "AML.mRNA.2020", "5/25/2018",
        "COHP_41251", "3695", "PBMC", "T6", "AML.mRNA.2020", "5/18/2018",
        "COHP_41252", "3695", "PBMC", "T7", "AML.mRNA.2020", "5/25/2018",
        "COHP_41259", "3696", "PBMC", "T6", "AML.mRNA.2020", "5/18/2018",
        "COHP_41260", "3696", "PBMC", "T7", "AML.mRNA.2020", "5/25/2018",
        "COHP_41267", "3697", "PBMC", "T6", "AML.mRNA.2020", "5/18/2018",
        "COHP_41268", "3697", "PBMC", "T7", "AML.mRNA.2020", "5/25/2018",
        "COHP_41275", "3698", "PBMC", "T6", "AML.mRNA.2020", "5/18/2018",
        "COHP_41276", "3698", "PBMC", "T7", "AML.mRNA.2020", "5/25/2018",
        "COHP_41283", "3700", "PBMC", "T6", "AML.mRNA.2020", "5/18/2018",
        "COHP_41284", "3700", "PBMC", "T7", "AML.mRNA.2020", "5/25/2018",
        "COHP_41291", "3701", "PBMC", "T6", "AML.mRNA.2020", "5/18/2018",
        "COHP_41292", "3701", "PBMC", "T7", "AML.mRNA.2020", "5/25/2018",
        "COHP_41299", "3702", "PBMC", "T6", "AML.mRNA.2020", "5/18/2018",
        "COHP_41300", "3702", "PBMC", "T7", "AML.mRNA.2020", "5/25/2018",
        "COHP_41307", "3706", "PBMC", "T6", "AML.mRNA.2020", "5/18/2018",
        "COHP_41308", "3706", "PBMC", "T7", "AML.mRNA.2020", "5/25/2018",
        "COHP_41315", "3709", "PBMC", "T6", "AML.mRNA.2020", "5/18/2018",
        "COHP_41316", "3709", "PBMC", "T7", "AML.mRNA.2020", "5/25/2018"
    ) |>
        dplyr::mutate(sample_date = lubridate::mdy(sample_date)) |>
        dplyr::select("sample", "mouse_id", "tissue", "timepoint", "project", "sample_date")


    aml_sample_metadata <- rbind(t6_t7_dates, missing_t6_t7_dates)
    # Get the per-mouse metadata for aml mice, emailed from Ya-Huei Kuo, 2022-9-2

    aml_mouse_metadata <- tibble::tribble(
        ~mouse_id, ~treatment, ~genotype, ~sex, ~dob, ~dod,
        "2683", "Ctrl", "WT", "M", "6/2/2015", "2/23/2016",
        "2684", "CM", "CHW", "M", "6/2/2015", "1/19/2016",
        "2685", "CM", "CHW", "M", "6/2/2015", "2/23/2016",
        "2686", "CM", "CHW", "M", "6/2/2015", "2/23/2016",
        "2689", "Ctrl", "CW", "F", "6/2/2015", "2/23/2016",
        "2690", "CM", "CHW", "F", "6/2/2015", "12/2/2015",
        "2692", "Ctrl", "CW", "F", "6/2/2015", "2/23/2016",
        "2700", "Ctrl", "CW", "F", "6/3/2015", "2/23/2016",
        "2702", "Ctrl", "CW", "F", "6/3/2015", "2/23/2016",
        "2705", "Ctrl", "HW", "M", "6/5/2015", "2/23/2016",
        "2708", "CM", "CHW", "F", "6/5/2015", "12/2/2015",
        "2709", "CM", "CHW", "F", "6/5/2015", "2/23/2016",
        "2718", "CM", "CHW", "M", "6/5/2015", "2/24/2016",
        "2719", "CM", "CHW", "M", "6/5/2015", "12/2/2015",
        "2720", "Ctrl", "HW", "M", "6/5/2015", "2/23/2016",
        "2731", "CM", "CHW", "F", "6/12/2015", "10/24/2015",
        "3127", "Ctrl", "CW", "M", "4/1/2016", "2/8/2017",
        "3130", "CM", "CHW", "F", "4/1/2016", "11/11/2016",
        "3131", "CM", "CHW", "F", "4/1/2016", "12/1/2016",
        "3200", "CM", "CHW", "M", "6/7/2016", "1/4/2017",
        "3202", "Ctrl", "WT", "M", "6/7/2016", "2/8/2017",
        "3334", "CM", "CHW", "M", "10/31/2016", "9/4/2017",
        "3335", "Ctrl", "HW", "M", "10/31/2016", "11/1/2017",
        "3336", "CM", "CHW", "F", "11/1/2016", "10/6/2017",
        "3338", "CM", "CHW", "F", "11/2/2016", "7/3/2017",
        "3339", "Ctrl", "CW", "M", "10/29/2016", "11/1/2017",
        "3340", "Ctrl", "HW", "M", "10/29/2016", "11/1/2017",
        "3341", "CM", "CHW", "F", "10/29/2016", "8/17/2017",
        "3342", "Ctrl", "CW", "M", "11/10/2016", "11/1/2017",
        "3346", "Ctrl", "CW", "M", "11/12/2016", "11/1/2017",
        "3349", "Ctrl", "CW", "F", "11/12/2016", "11/1/2017",
        "3350", "Ctrl", "WT", "F", "11/12/2016", "11/1/2017",
        "3357", "CM", "CHW", "F", "11/15/2016", "10/6/2017",
        "3368", "CM", "CHW", "M", "12/7/2016", "11/1/2017",
        "3370", "CM", "CHW", "M", "12/8/2016", "11/1/2017",
        "3694", "Ctrl", "WT", "M", "2/1/2018", "8/16/2018",
        "3695", "Ctrl", "WT", "M", "2/1/2018", "8/16/2018",
        "3696", "Ctrl", "WT", "M", "2/1/2018", "8/16/2018",
        "3697", "3Gy", "WT", "F", "2/1/2018", "8/16/2018",
        "3698", "3Gy", "WT", "F", "2/1/2018", "8/16/2018",
        "3700", "3Gy", "WT", "F", "2/1/2018", "8/16/2018",
        "3701", "3Gy", "WT", "M", "2/1/2018", "8/16/2018",
        "3702", "3Gy", "WT", "M", "2/1/2018", "8/16/2018",
        "3706", "Ctrl", "WT", "F", "2/1/2018", "8/16/2018",
        "3709", "Ctrl", "WT", "F", "2/1/2018", "8/16/2018",
        "4309", "CM", "CHW", "F", "11/10/2019", "9/29/2020",
        "4311", "Ctrl", "CW", "F", "11/10/2019", "7/12/2021",
        "4321", "CM", "CHW", "F", "11/13/2019", "7/17/2020",
        "4324", "CM", "CHW", "M", "11/13/2019", "7/31/2020",
        "4329", "CM", "CHW", "F", "11/19/2019", "9/9/2020",
        "4419", "CM", "CHW", "F", "7/14/2020", "4/2/2021",
        "4433", "CM", "CHW", "M", "8/21/2020", "4/19/2021",
        "4436", "CM", "CHW", "F", "8/21/2020", "3/21/2021",
        "4443", "CM", "CHW", "F", "8/25/2020", "5/26/2021",
        "4498", "Ctrl", "CW", "M", "3/3/2021", "10/13/2021",
        "4501", "Ctrl", "CW", "F", "3/3/2021", "12/7/2021",
        "4502", "CM", "CHW", "M", "3/3/2021", "10/13/2021",
        "4506", "CM", "CHW", "F", "3/3/2021", "11/17/2021",
        "4510", "CM", "CHW", "F", "3/11/2021", "11/8/2021",
        "4512", "Ctrl", "WT", "M", "3/25/2021", "10/6/2021",
        "4520", "Ctrl", "CW", "F", "4/9/2021", "9/29/2021",
        "4521", "CM", "CHW", "M", "4/26/2021", "10/6/2021",
        "4522", "CM", "CHW", "F", "4/26/2021", "9/29/2021",
        "4534", "Ctrl", "WT", "F", "5/22/2021", "12/7/2021",
        "4535", "CM", "CHW", "F", "5/22/2021", "12/7/2021"
    ) |>
        dplyr::mutate(
            dob = lubridate::mdy(dob),
            dod = lubridate::mdy(dod)
        )

    aml_mice <- dplyr::left_join(aml_sample_metadata, aml_mouse_metadata, by = "mouse_id")

    # CML sample and mouse metadata
    use_pinboard("onedrive")
    cml_mice <- get_pin("metadata_mmu.csv", "20220904T234748Z-f5f2a") |>
        dplyr::filter(grepl("CML", project)) |>
        dplyr::select("sample", "fastq_1", "fastq_2", "strandedness", "mouse_id", "tissue", "timepoint", "project", "treatment", "genotype", "sex", "sample_date", "dob", "dod")

    cml_mice[cols_to_update_character] <- lapply(cml_mice[cols_to_update_character], as.character)
    cml_mice[cols_to_update_date] <- lapply(cml_mice[cols_to_update_date], as.Date)

    # consolidate sample metadata where possible
    # using dplyr::rows_update() & tidyr::fill() to fill in missing values
    sample_sheet <- dplyr::rows_update(all_mice, aml_mice, by = "sample", unmatched = "ignore") |>
        dplyr::rows_update(cml_mice, by = "sample", unmatched = "ignore") |>
        dplyr::group_by(mouse_id) |>
        tidyr::fill(c("treatment", "genotype", "sex", "dob", "dod"), .direction = "downup") |>
        dplyr::mutate(
            # add columns: assay, sample_weeks, age_at_end, age_at_start, age_at_sample
            assay = "mRNA",
            sample_weeks = difftime(sample_date, min(sample_date), units = "weeks"),
            age_at_end = difftime(max(sample_date), dob, units = "weeks"),
            age_at_start = difftime(min(sample_date), dob, units = "weeks"),
            age_at_sample = difftime(sample_date, dob, units = "weeks"),
            dplyr::across(dplyr::starts_with(c("age", "sample_weeks")), round, 1)
        ) |>
        dplyr::ungroup()

        # add batch
        use_pinboard("onedrive")
        batch <- get_pin("metadata_mmu.csv", "20220904T234748Z-f5f2a") |>
                dplyr::select("sample", "batch")
        sample_sheet <- dplyr::left_join(sample_sheet, batch, by = "sample") |>
            dplyr::select(
                sample, fastq_1, fastq_2, strandedness, assay, mouse_id,
                tissue, timepoint, project, batch, treatment, genotype, sex,
                sample_date, dob, dod, sample_weeks, age_at_end, age_at_start,
                age_at_sample
            ) |>
            dplyr::distinct()

        # add percent_ckit
        get_teams_file("General/cKit+ AllTimes.csv")
        percent_ckit <- read.csv("data-raw/cKit+ AllTimes.csv", colClasses = c("mouse_id" = "character")) |>
            tidyr::pivot_longer(
                cols = -mouse_id,
                names_to = "timepoint",
                values_to = "percent_ckit"
            ) |>
            na.omit()

        sample_sheet <- dplyr::left_join(sample_sheet, percent_ckit, by = c("mouse_id", "timepoint")) |>
            dplyr::select(
                sample, fastq_1, fastq_2, strandedness, assay, mouse_id,
                tissue, timepoint, project, batch, treatment, genotype, sex,
                sample_date, percent_ckit, dob, dod, sample_weeks, age_at_end, age_at_start,
                age_at_sample
            ) |>
            dplyr::distinct()

    return(sample_sheet)
}
#### mRNA -----

# ├ AML.mRNA.2016
parse_metadata_AML.mRNA.2016 <- function() {
    # TODO update to use file paths from new isilon - symlinks are broken
    sample_sheet <- read.csv("/net/nfs-irwrsrchnas01/labs/rrockne/MHO/AML.mRNA.2016/config/new_name_key.tsv", sep = "\t") |>
        dplyr::mutate(
            fastq_1 = paste0("/net/nfs-irwrsrchnas01/labs/rrockne/MHO/AML.mRNA.2016/data/", newname, ".gz"),
            sample = paste0("COHP_", ID),
            fastq_2 = "",
            strandedness = dplyr::case_when(
                stringr::str_detect(ID, "11548") ~ "reverse",
                stringr::str_detect(TimePoint, "L|T0|T1") ~ "unstranded",
                TRUE ~ "reverse"
            )
        ) |>
        dplyr::select(sample, fastq_1, fastq_2, strandedness)
    return(sample_sheet)
}

# AML.mRNA.2018.all_samples
parse_metadata_AML.mRNA.2018.all_samples <- function() {
    # TODO update to use file paths from new isilon - symlinks are broken
    sample_sheet <- read.csv("/net/nfs-irwrsrchnas01/labs/rrockne/MHO/AML.mRNA.2018.all_samples/config/new_name_key.tsv", sep = "\t") |>
        dplyr::mutate(
            sample = paste0("COHP_", gsub(".*_|\\.fq", "", Newname)),
            fastq_1 = paste0("/net/nfs-irwrsrchnas01/labs/rrockne/MHO/AML.mRNA.2018.all_samples/data/fastq/", Newname, ".gz"),
            fastq_2 = "", strandedness = "reverse"
        ) |>
        dplyr::select(sample, fastq_1, fastq_2, strandedness)
    return(sample_sheet)
}

# AML.mRNA.2020
parse_metadata_AML.mRNA.2020 <- function() {
    # Copied from "Seq Samples Dates.xlsx":  https://github.com/drejom/haemdata/issues/6
    sample_sheet <- data.frame(
        fastq_1 = system(
            paste0("find '/net/nfs-irwrsrchnas01/labs/ykuo/Seq/201124_TGen' -name '*.gz'"),
            intern = TRUE
        ) %>% grep("_R1_", ., value = TRUE),
        fastq_2 = system(
            paste0("find '/net/nfs-irwrsrchnas01/labs/ykuo/Seq/201124_TGen' -name '*.gz'"),
            intern = TRUE
        ) %>% grep("_R2_", ., value = TRUE)
    ) |>
        dplyr::mutate(sample = stringr::str_extract(fastq_1, "COHP_\\d{5}"), strandedness = "reverse") |>
        dplyr::select(sample, fastq_1, fastq_2, strandedness)
}
# AML.mRNA.2021.RxGroup1
parse_metadata_AML.mRNA.2021.RxGroup1 <- function() {
    sample_sheet <- data.frame(
        fastq_1 = system(
            paste0("find '/net/nfs-irwrsrchnas01/labs/rrockne/Seq/210412_Tgen' -name '*.gz'"),
            intern = TRUE
        ) %>% grep("_R1_", ., value = TRUE),
        fastq_2 = system(
            paste0("find '/net/nfs-irwrsrchnas01/labs/rrockne/Seq/210412_Tgen' -name '*.gz'"),
            intern = TRUE
        ) %>% grep("_R2_", ., value = TRUE)
    ) |>
        dplyr::mutate(sample = stringr::str_extract(fastq_1, "COHP_\\d{5}"), strandedness = "reverse") |>
        dplyr::select(sample, fastq_1, fastq_2, strandedness)
    return(sample_sheet)
}
# AML.mRNA.2021.RxGroup2
parse_metadata_AML.mRNA.2021.RxGroup2 <- function() {
    sample_sheet <- data.frame(
        fastq_1 = system(
            paste0("find '/net/nfs-irwrsrchnas01/labs/rrockne/Seq/210507_Tgen' -name '*.gz'"),
            intern = TRUE
        ) %>% grep("_R1_", ., value = TRUE),
        fastq_2 = system(
            paste0("find '/net/nfs-irwrsrchnas01/labs/rrockne/Seq/210507_Tgen' -name '*.gz'"),
            intern = TRUE
        ) %>% grep("_R2_", ., value = TRUE)
    ) |>
        dplyr::mutate(sample = stringr::str_extract(fastq_1, "COHP_\\d{5}"), strandedness = "reverse") |>
        dplyr::select(sample, fastq_1, fastq_2, strandedness)
    return(sample_sheet)
}

# AML.mRNA.2021.RxGroup2_pt2
parse_metadata_AML.mRNA.2021.RxGroup2_pt2 <- function() {
    sample_sheet <- data.frame(
        fastq_1 = system(
            paste0("find '/net/nfs-irwrsrchnas01/labs/ykuo/Seq/210910_TGen' -name '*.gz'"),
            intern = TRUE
        ) %>% grep("_R1_", ., value = TRUE),
        fastq_2 = system(
            paste0("find '/net/nfs-irwrsrchnas01/labs/ykuo/Seq/210910_TGen' -name '*.gz'"),
            intern = TRUE
        ) %>% grep("_R2_", ., value = TRUE)
    ) |>
        dplyr::mutate(sample = stringr::str_extract(fastq_1, "COHP_\\d{5}"), strandedness = "reverse") |>
        dplyr::select(sample, fastq_1, fastq_2, strandedness)
    return(sample_sheet)
}
# AML.mRNA.2022.RxGroup3
parse_metadata_AML.mRNA.2022.RxGroup3 <- function() {
    sample_sheet <- data.frame(
        fastq_1 = system(
            paste0("find '/net/nfs-irwrsrchnas01/labs/ykuo/Seq/220415_IGC-LZ-20205' -name '*.gz'"),
            intern = TRUE
        ) %>% grep("_R1_", ., value = TRUE),
        fastq_2 = system(
            paste0("find '/net/nfs-irwrsrchnas01/labs/ykuo/Seq/220415_IGC-LZ-20205' -name '*.gz'"),
            intern = TRUE
        ) %>% grep("_R2_", ., value = TRUE)
    ) |>
        dplyr::mutate(sample = stringr::str_extract(fastq_1, "COHP_\\d{5}"), strandedness = "reverse") |>
        dplyr::select(sample, fastq_1, fastq_2, strandedness)
    return(sample_sheet)
}
# AML.mRNA.novaseq_validation.2020
parse_metadata_AML.mRNA.novaseq_validation.2020 <- function() {
    sample_sheet <- data.frame(
        fastq_1 = system(
            paste0("find '/net/nfs-irwrsrchnas01/labs/ykuo/Seq/200928_TGen' -name '*.gz'"),
            intern = TRUE
        ) %>% grep("_R1_", ., value = TRUE),
        fastq_2 = system(
            paste0("find '/net/nfs-irwrsrchnas01/labs/ykuo/Seq/200928_TGen' -name '*.gz'"),
            intern = TRUE
        ) %>% grep("_R2_", ., value = TRUE)
    ) |>
        dplyr::mutate(sample = stringr::str_extract(fastq_1, "COHP_\\d{5}"), strandedness = "reverse") |>
        dplyr::select(sample, fastq_1, fastq_2, strandedness)
    return(sample_sheet)
}
# AML.validation.2017
parse_metadata_AML.validation.2017 <- function() {
    sample_sheet <- data.frame(
        fastq_1 = system(
            paste0("find '/net/nfs-irwrsrchnas01/labs/rrockne/MHO/AML.validation.2017/data/fastq' -name '*.gz'"),
            intern = TRUE
        ) %>% grep("_R1_", ., value = TRUE)
    ) |>
        dplyr::mutate(sample = paste0("COHP_", substr(basename(fastq_1), 1, 5)), fastq_2 = "", strandedness = "reverse") |>
        dplyr::select(sample, fastq_1, fastq_2, strandedness)
    return(sample_sheet)
}

# AML.scRNAseq.2022
parse_metadata_AML.scRNAseq.2022 <- function() {
    sample_sheet <- data.frame(
        fastq_1 = system(
            paste0("find '/net/nfs-irwrsrchnas01/labs/rrockne/Seq/211210_IGC-LZ-19773' -name '*.gz'"),
            intern = TRUE
        ) %>% grep("_R1_", ., value = TRUE),
        fastq_2 = system(
            paste0("find '/net/nfs-irwrsrchnas01/labs/rrockne/Seq/211210_IGC-LZ-19773' -name '*.gz'"),
            intern = TRUE
        ) %>% grep("_R2_", ., value = TRUE)
    ) |>
        dplyr::mutate(sample = stringr::str_extract(fastq_1, "COHP_\\d{5}"), strandedness = "reverse") |>
        dplyr::select(sample, fastq_1, fastq_2, strandedness)
    return(sample_sheet)
}

# CML.mRNA.2021
parse_metadata_CML.mRNA.2021 <- function() {
    sample_sheet <- data.frame(
        fastq_1 = system(
            paste0("find /net/nfs-irwrsrchnas01/labs/rrockne/Seq/210910_TGen /net/nfs-irwrsrchnas01/labs/rrockne/Seq/210907_TGen -name '*.gz'"),
            intern = TRUE
        ) %>% grep("_R1_", ., value = TRUE),
        fastq_2 = system(
            paste0("find /net/nfs-irwrsrchnas01/labs/rrockne/Seq/210910_TGen /net/nfs-irwrsrchnas01/labs/rrockne/Seq/210907_TGen -name '*.gz'"),
            intern = TRUE
        ) %>% grep("_R2_", ., value = TRUE)
    ) |>
        dplyr::mutate(sample = stringr::str_extract(fastq_1, "COHP_\\d{5}"), strandedness = "reverse") |>
        dplyr::select(sample, fastq_1, fastq_2, strandedness)
    return(sample_sheet)
}

# CML.mRNA.2022
parse_metadata_CML.mRNA.2022 <- function() {
    sample_sheet <- data.frame(
        fastq_1 = system(
            paste0("find /net/nfs-irwrsrchnas01/labs/gmarcucci/Seq/220210_IGC-LZ-19931\\ CML\\ PBMC\\ RNA-seq /net/nfs-irwrsrchnas01/labs/gmarcucci/Seq/220202_IGC-LZ-19931\\ CML\\ PBMC\\ RNA-seq -name '*.gz'"),
            intern = TRUE
        ) %>% grep("_R1_", ., value = TRUE),
        fastq_2 = system(
            paste0("find /net/nfs-irwrsrchnas01/labs/gmarcucci/Seq/220210_IGC-LZ-19931\\ CML\\ PBMC\\ RNA-seq /net/nfs-irwrsrchnas01/labs/gmarcucci/Seq/220202_IGC-LZ-19931\\ CML\\ PBMC\\ RNA-seq -name '*.gz'"),
            intern = TRUE
        ) %>% grep("_R2_", ., value = TRUE)
    ) |>
        dplyr::mutate(sample = stringr::str_extract(fastq_1, "COHP_\\d{5}"), strandedness = "reverse") |>
        dplyr::select(sample, fastq_1, fastq_2, strandedness)
    return(sample_sheet)
}

# CML.mRNA.2022_pt2
parse_metadata_CML.mRNA.2022_pt2 <- function() {
    sample_sheet <- data.frame(
        fastq_1 = system(
            paste0("find /net/nfs-irwrsrchnas01/labs/gmarcucci/Seq/220811_0471_IGC-LZ-20847 /net/nfs-irwrsrchnas01/labs/gmarcucci/Seq/220811_0472_IGC-LZ-20847 -name '*.gz'"),
            intern = TRUE
        ) %>% grep("_R1_", ., value = TRUE),
        fastq_2 = system(
            paste0("find /net/nfs-irwrsrchnas01/labs/gmarcucci/Seq/220811_0471_IGC-LZ-20847 /net/nfs-irwrsrchnas01/labs/gmarcucci/Seq/220811_0472_IGC-LZ-20847 -name '*.gz'"),
            intern = TRUE
        ) %>% grep("_R2_", ., value = TRUE)
    ) |>
        dplyr::mutate(sample = stringr::str_extract(fastq_1, "COHP_\\d{5}"), strandedness = "reverse") |>
        dplyr::select(sample, fastq_1, fastq_2, strandedness)
    return(sample_sheet)
}
# AML.mRNA.HSA_FLT3.2022
parse_metadata_AML.mRNA.HSA_FLT3.2022 <- function() {
    project <- "AML.mRNA.HSA_FLT3.2022"
    fastqs <- data.frame(
        fastq_1 = system(
            paste0("find /net/nfs-irwrsrchnas01/labs/ykuo/Seq/220511_IGC-LZ-20342 -name '*.gz'"),
            intern = TRUE
        ) %>% grep("_R1_", ., value = TRUE),
        fastq_2 = system(
            paste0("find /net/nfs-irwrsrchnas01/labs/ykuo/Seq/220511_IGC-LZ-20342 -name '*.gz'"),
            intern = TRUE
        ) %>% grep("_R2_", ., value = TRUE)
    ) |> dplyr::mutate(library_id = stringr::str_extract(fastq_1, "COHP_\\d{5}"))
    xls <- "data-raw/FLT3 AML samples information (IGC-LZ-20342).xlsx"
    xls_ss <- "data-raw/sample summary_IGC-LZ-20342.xlsx"
    sample_metadata <- readxl::read_excel(xls) |>
        janitor::clean_names("lower_camel") |>
        janitor::remove_constant(na.rm = TRUE) |>
        tidyr::fill(number, .direction = "down") |>
        tidyr::separate(htbIdPb, into = c(NA, NA, "sample_id"), sep = "-") |>
        dplyr::select(patient_id = number, sample_id, date, cell_amount = cellAmt) |>
        tidyr::drop_na(sample_id) |>
        dplyr::arrange(patient_id, date) |>
        dplyr::group_by(patient_id) |>
        dplyr::mutate(
            patient_id = gsub("#", "_", patient_id, fixed = TRUE),
            weeks = as.numeric(round((date - min(date)) / 604800, 2)),
            first_sample_date = min(date)
        ) |>
        dplyr::ungroup()
    sequencing_metadata <- readxl::read_excel(xls_ss) |>
        tidyr::separate(Sample_ID, into = c(NA, NA, "sample_id"), sep = "_") |>
        dplyr::mutate(
            project = project,
            library_id = TGen_Sample_Name,
            batch = "2022_C",
            strandedness = "reverse",
            sex = NA_character_,
            dob = NA,
            treatment = NA_character_,
            tissue = "PBMC"
        ) |>
        dplyr::select(sample_id, project, library_id, batch, strandedness, sex, dob, treatment, tissue) |>
        dplyr::left_join(fastqs)
    sample_sheet <- dplyr::left_join(sequencing_metadata, sample_metadata, by = "sample_id") |>
        dplyr::arrange(patient_id, weeks) |>
        dplyr::select(sample = library_id, fastq_1, fastq_2, strandedness, dplyr::everything())
    return(sample_sheet)
}
# MDS.rnaseq.EGAD00001003891
parse_metadata_MDS.rnaseq.EGAD00001003891 <- function() {
    project <- "MDS.rnaseq.EGAD00001003891"
    fastqs <- data.frame(
        fastq_1 = list.files("/net/nfs-irwrsrchnas01/labs/rrockne/MHO/MDS.rnaseq.EGAD00001003891/data/fastq2/reads",
            pattern = "1.fq.gz", full.names = TRUE
        ),
        fastq_2 = list.files("/net/nfs-irwrsrchnas01/labs/rrockne/MHO/MDS.rnaseq.EGAD00001003891/data/fastq2/reads",
            pattern = "2.fq.gz", full.names = TRUE
        )
    ) |> dplyr::mutate(
        ega_run_accession_id = gsub("_.*$", "", basename(fastq_1)),
        sample_name = gsub("^.*?_|.1.fq.gz", "", basename(fastq_1))
    )
    sex <- read.table("data-raw/EGAD00001003891/delimited_maps/Run_Sample_meta_info.map", sep = "\t", header = FALSE) |>
        dplyr::mutate(
            patient_id = stringr::str_extract(V1, "(?<=subject_id\\=)(.*?)(?=[;])"),
            sex = stringr::str_extract(V1, "(?<=gender\\=)(.*?)(?=[;])"),
            sex = dplyr::recode(sex, male = "M", female = "F", unknown = NA_character_)
        ) |>
        dplyr::select(patient_id, sex) |>
        dplyr::distinct()
    files <- read.table("data-raw/EGAD00001003891/delimited_maps/Sample_File.map",
        sep = "\t", header = FALSE,
        col.names = c("patient_id", "ega_sample_accession_id", "bam", "ega_file_unique_accession_id")
    ) |> dplyr::mutate(sample_name = gsub(".bam.cip", "", bam))

    experiment <- read.csv("data-raw/EGAD00001003891/delimited_maps/Study_Experiment_Run_sample.map",
        sep = "\t", header = FALSE,
        col.names = c(
            "study_accession", "study_description", "existing_study_type", "platform",
            "instrument_model", "library_layout", "library_name", "library_strategy",
            "library_source", "library_selection", "ega_experiment_id", "ega_run_accession_id",
            "submitter_id", "na", "ega_sample_accession_id"
        )
    )
    ega <- dplyr::left_join(
        dplyr::left_join(
            fastqs,
            dplyr::left_join(files, sex, by = "patient_id"),
            by = "sample_name"
        ),
        experiment,
        by = "ega_run_accession_id"
    )

    sample_sheet <- ega |>
        dplyr::mutate(
            batch = "2022_D",
            strandedness = "unstranded",
            treatment = NA_character_,
            tissue = dplyr::case_when(
                stringr::str_detect(sample_name, "Ctrl|Mut|WT") ~ "HEK293T",
                stringr::str_detect(sample_name, "DNA") ~ "BM",
                stringr::str_detect(sample_name, "BMMNC") ~ "BMMNC",
                stringr::str_detect(sample_name, "293T") ~ "HEK293T",
                stringr::str_detect(sample_name, "day7|day14") ~ "CD34",
                stringr::str_detect(sample_name, "CD34") ~ "CD34"
            ),
            project = project
        ) |>
        dplyr::select(sample = ega_run_accession_id, fastq_1, fastq_2, strandedness, dplyr::everything(), project) |>
        dplyr::filter(library_strategy == "RNA-Seq")
    return(sample_sheet)
}

# AML.PRJEB27973 Kim et al 2020 SciRep
parse_metadata_AML.PRJEB27973 <- function() {
    sample_sheet <- readr::read_csv(
        "/net/nfs-irwrsrchnas01/labs/rrockne/MHO/AML.PRJEB27973/samplesheet/samplesheet.csv",
        show_col_types = FALSE
    ) |>
        tidyr::separate(sample_alias, into = c("genotype", "patient_id", "timepoint"), sep = "-") |>
        dplyr::mutate(
            project = "AML.PRJEB27973",
            sample = sample,
            patient_id = patient_id,
            timepoint = timepoint,
            batch = "PRJEB27973",
            strandedness = "reverse",
            sex = NA_character_,
            dob = NA,
            treatment = "",
            genotype = genotype,
            tissue = "BM"
        ) |>
        dplyr::select(sample, fastq_1, fastq_2, strandedness, patient_id, tissue, timepoint, batch, treatment, genotype, sex, dob, project, everything())
    return(sample_sheet)
}

#### microRNA -----

# ├ AML.mRNA.2016
parse_metadata_AML.miRNA.2016 <- function() {
    project <- "AML.miRNA.2016"
    assay <- "microRNA"
    tissue <- "PBMC"
    sample_sheet <- read.csv("data-raw/2016miRNA_new_name_key.csv") |>
        dplyr::mutate(
            project = project,
            assay = assay,
            tissue = tissue,
            fastq_1 = oldname,
            library_id = paste0("COHP_", stringr::str_replace(oldname, ".*/(\\d*)_.*", "\\1")),
            mouse_id = stringr::str_replace(newname, "(\\d)-.*", "\\1"),
            timepoint = stringr::str_replace(newname, "^\\d*-(.)\\.fq", "\\1"),
            seq_batch = as.numeric(as.factor(stringr::str_replace(oldname, ".*(?>Seq/)(\\w*)/.*", "\\1"))),
            batch = paste0("2016_", seq_batch)
        ) |>
        dplyr::select(sample = library_id, fastq_1, mouse_id, timepoint, tissue, batch, project, assay)
    return(sample_sheet)
}
# ├ AML.mRNA.2018
parse_metadata_AML.miRNA.2018 <- function() {
    project <- "AML.miRNA.2018"
    assay <- "microRNA"
    tissue <- "PBMC"
    fastq_paths <- list.files(
        path = "/net/nfs-irwrsrchnas01/labs/ykuo/Seq/AML2018_new_samples/2018_Aug_miRNA/original_fastq",
        pattern = ".fastq", full.names = TRUE
    ) |>
        dplyr::as_tibble() |>
        dplyr::mutate(
            library_id = paste0("COHP_", stringr::str_replace(value, ".*/(\\d*)_.*", "\\1"))
        )
    sample_sheet <- read.csv("data-raw/2018miRNA_new_name_key.txt", sep = "\t") |>
        dplyr::mutate(
            project = project,
            assay = assay,
            tissue = tissue,
            library_id = paste0("COHP_", stringr::str_replace(OldName, "(\\d*)_.*", "\\1")),
            mouse_id = stringr::str_replace(NewName, "[:alnum:]+_(\\d+)_.*", "\\1"),
            timepoint = stringr::str_replace(NewName, "([:alnum:]+)_.*", "\\1"),
            batch = stringr::str_replace(NewName, "(?:[^_]+_){6}([^_]*).*", "\\1")
        ) |>
        dplyr::left_join(fastq_paths, by = "library_id") |>
        dplyr::select(sample = library_id, fastq_1 = value, mouse_id, timepoint, tissue, batch, project, assay)
    return(sample_sheet)
}
# ├ AML.mRNA.2020
# * NOTE: Unique to this miRNA project, the samples
# * were paired end, but only fastq_1 was used.
parse_metadata_AML.miRNA.2020 <- function() {
    project <- "AML.miRNA.2020"
    assay <- "microRNA"
    batch <- "2020_U"
    tissue <- "PBMC"
    sample_sheet <- list.files(
        path = "/net/nfs-irwrsrchnas01/labs/ykuo/Seq/210211_miRNA",
        pattern = ".fastq", full.names = TRUE
    ) |>
        dplyr::as_tibble() |>
        dplyr::mutate(
            library_id = paste0("COHP_", stringr::str_replace(value, ".*/(\\d*)_.*", "\\1")),
            project = project,
            assay = assay,
            fastq_1 = value,
            mouse_id = stringr::str_replace(value, ".*/(?:[^_]+_){2}([^_]*).*", "\\1"),
            timepoint = stringr::str_replace(value, ".*/(?:[^_]+_){3}([^_]*).*", "\\1"),
            batch = batch,
            tissue = dplyr::case_when(
                str_detect(timepoint, "S..") ~ "BM",
                TRUE ~ tissue
            ),
            timepoint = dplyr::case_when(
                tissue == "BM" ~ NA_character_,
                TRUE ~ timepoint
            )
        ) |>
        dplyr::select(sample = library_id, fastq_1, mouse_id, timepoint, tissue, batch, project, assay)
    return(sample_sheet)
}

# ├ AML.miRNA.2021.RxGroup1
parse_metadata_AML.miRNA.2021.RxGroup1 <- function() {
    project <- "AML.miRNA.2021.RxGroup1"
    assay <- "microRNA"
    tissue <- "PBMC"
    batch <- "2021_G"
    sample_sheet <- list.files(
        path = "/net/nfs-irwrsrchnas01/labs/rrockne/Seq/210429_B",
        pattern = ".fastq", full.names = TRUE
    ) |>
        dplyr::as_tibble() |>
        dplyr::mutate(
            library_id = paste0("COHP_", stringr::str_replace(value, ".*/(\\d*)_.*", "\\1")),
            project = project,
            assay = assay,
            fastq_1 = value,
            mouse_id = stringr::str_replace(value, ".*/(?:[^_]+_){2}(\\d{4}).*", "\\1"),
            timepoint = stringr::str_replace(value, ".*/(?:[^_]+_){2}\\d{4}([^_]*)_.*", "\\1"),
            batch = batch,
            tissue = tissue,
            tissue = dplyr::case_when(
                timepoint == "BM" ~ "BM",
                TRUE ~ tissue
            ),
            timepoint = sub("BM", NA_character_, timepoint)
        ) |>
        dplyr::select(sample = library_id, fastq_1, mouse_id, timepoint, tissue, batch, project, assay)
    return(sample_sheet)
}

# ├ AML.miRNA.2021.RxGroups1and2
parse_metadata_AML.miRNA.2021.RxGroups1and2 <- function() {
    project <- "AML.miRNA.2021.RxGroups1and2"
    assay <- "microRNA"
    tissue <- "PBMC"
    batch <- "2021_H"
    sample_sheet <- list.files(
        path = "/net/nfs-irwrsrchnas01/labs/rrockne/Seq/210528_A",
        pattern = ".fastq", full.names = TRUE
    ) |>
        dplyr::as_tibble() |>
        dplyr::mutate(
            library_id = paste0("COHP_", stringr::str_replace(value, ".*/(\\d*)_.*", "\\1")),
            project = project,
            assay = assay,
            fastq_1 = value,
            mouse_id = stringr::str_replace(value, ".*/(?:[^_]+_){2}(\\d{4}).*", "\\1"),
            timepoint = stringr::str_replace(value, ".*/(?:[^_]+_){2}\\d{4}([^_]*)_.*", "\\1"),
            batch = batch,
            tissue = tissue,
            tissue = dplyr::case_when(
                timepoint == "BM" ~ "BM",
                TRUE ~ tissue
            ),
            timepoint = sub("BM", NA_character_, timepoint),
            timepoint = dplyr::case_when(
                stringr::str_detect(fastq_1, "_T1") ~ "T1",
                stringr::str_detect(fastq_1, "_T2") ~ "T2",
                TRUE ~ timepoint
            )
        ) |>
        dplyr::select(sample = library_id, fastq_1, mouse_id, timepoint, tissue, batch, project, assay)
    return(sample_sheet)
}

# ├ AML.miRNA.2021.RxGroup2_pt2
parse_metadata_AML.miRNA.2021.RxGroup2_pt2 <- function() {
    project <- "AML.miRNA.2021.RxGroup2_pt2"
    assay <- "microRNA"
    tissue <- "PBMC"
    batch <- "2021_I"
    sample_sheet <- list.files(
        path = "/net/nfs-irwrsrchnas01/labs/ykuo/Seq/210907",
        pattern = ".fastq", full.names = TRUE
    ) |>
        dplyr::as_tibble() |>
        dplyr::mutate(
            library_id = paste0("COHP_", stringr::str_replace(value, ".*/(\\d*)_.*", "\\1")),
            project = project,
            assay = assay,
            fastq_1 = value,
            mouse_id = stringr::str_replace(value, ".*/(?:[^_]+_){2}(\\d{4}).*", "\\1"),
            timepoint = stringr::str_replace(value, ".*/(?:[^_]+_){2}\\d{4}-([^_]*)_.*", "\\1"),
            batch = batch,
            tissue = tissue,
            tissue = dplyr::case_when(
                timepoint == "BM" ~ "BM",
                TRUE ~ tissue
            ),
            timepoint = sub("BM", NA_character_, timepoint)
        ) |>
        dplyr::select(sample = library_id, fastq_1, mouse_id, timepoint, tissue, batch, project, assay)
    return(sample_sheet)
}


# ├ AML.miRNA.2021.RxGroup3
parse_metadata_AML.miRNA.2022.RxGroup3 <- function() {
    project <- "AML.miRNA.2022.RxGroup3"
    assay <- "microRNA"
    tissue <- "PBMC"
    batch <- "2022_D"
    fastq_paths <- list.files(
        path = "/net/nfs-irwrsrchnas01/labs/ykuo/Seq/220620_IGC-LZ-20205",
        pattern = ".fastq", full.names = TRUE
    ) |>
        dplyr::as_tibble() |>
        dplyr::mutate(library_id = paste0("COHP_", stringr::str_replace(value, ".*/(\\d*)_.*", "\\1")))


    # sample_sheet <-
    readxl::read_excel("data-raw/sample summary_IGC-LZ-20205_miRNA.xlsx") |>
        dplyr::mutate(
            library_id = paste0("COHP_", stringr::str_replace(Sample_ID, "(\\d*)_.*", "\\1")),
            project = project,
            assay = assay,
            mouse_id = stringr::str_replace(Sample_ID, "(?:[^_]+_){2}(\\d{4}).*", "\\1"),
            timepoint = stringr::str_replace(Sample_ID, "(?:[^_]+_){2}\\d{4}(.*)", "\\1"),
            batch = batch,
            tissue = tissue,
            tissue = dplyr::case_when(
                timepoint == "BM" ~ "BM",
                TRUE ~ tissue
            ),
            timepoint = sub("BM", NA_character_, timepoint)
        ) |>
        dplyr::left_join(fastq_paths, by = "library_id") |>
        dplyr::select(sample = library_id, fastq_1 = value, mouse_id, timepoint, tissue, batch, project, assay)
    return(sample_sheet)
}
