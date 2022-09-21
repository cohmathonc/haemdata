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

    all_mice <- sample_sheet_all_mice |>
        dplyr::mutate(
            # add in new columns to facilitate join() and update_rows():
            mouse_tp = glue::glue("{mouse_id}_{timepoint}_{tissue}"),
            dob = as.Date(NA),
            dod = as.Date(NA),
            sample_date = as.Date(NA)
        )


    # AML sample and mouse metadata

    # download the file 
    get_teams_file("General/Copy of matadata_mmu_pivoted_AMLmice.YK.xlsx")

    # read in the worksheets
    wide_dates <- readxl::read_excel(
        "data-raw/Copy of matadata_mmu_pivoted_AMLmice.YK.xlsx",
        col_names = TRUE,
        col_types = "text") |>
        dplyr::select(-c("Helper", "cohp"))

    missing_t6_t7_dates <- tidyr::pivot_longer(
        data = wide_dates,
        cols = !c("sample", "project", "mouse_id", "tissue"),
        names_to = "timepoint",
        values_to = "sample_date") |>
        tidyr::drop_na("sample_date") |> 
        dplyr::mutate(sample_date = as.Date(sample_date)) |> 
        dplyr::arrange(mouse_id, sample_date) 

    # T6 was left off the sheet sent to Ya-Huei, so we need to add it in separately
    t6_t7_dates <- tibble::tribble(
                                ~sample, ~mouse_id, ~tissue, ~timepoint,                    ~project, ~sample_date,
                            "COHP_11843",     "2683",  "PBMC",       "T6",             "AML.mRNA.2016",    "2/23/2016",
                            "COHP_11844",     "2685",  "PBMC",       "T6",             "AML.mRNA.2016",    "2/23/2016",
                            "COHP_11845",     "2686",  "PBMC",       "T6",             "AML.mRNA.2016",    "2/23/2016",
                            "COHP_11846",     "2689",  "PBMC",       "T6",             "AML.mRNA.2016",    "2/23/2016",
                            "COHP_11847",     "2692",  "PBMC",       "T6",             "AML.mRNA.2016",    "2/23/2016",
                            "COHP_11848",     "2700",  "PBMC",       "T6",             "AML.mRNA.2016",    "2/23/2016",
                            "COHP_11849",     "2702",  "PBMC",       "T6",             "AML.mRNA.2016",    "2/23/2016",
                            "COHP_11850",     "2705",  "PBMC",       "T6",             "AML.mRNA.2016",    "2/23/2016",
                            "COHP_11851",     "2709",  "PBMC",       "T6",             "AML.mRNA.2016",    "2/23/2016",
                            "COHP_11852",     "2720",  "PBMC",       "T6",             "AML.mRNA.2016",    "2/23/2016",
                            "COHP_20898",     "3335",  "PBMC",       "T6", "AML.mRNA.2018.all_samples",    "8/15/2017",
                            "COHP_20978",     "3336",  "PBMC",       "T6", "AML.mRNA.2018.all_samples",    "8/15/2017",
                            "COHP_20919",     "3336",  "PBMC",       "T7", "AML.mRNA.2018.all_samples",    "9/12/2017",
                            "COHP_20922",     "3341",  "PBMC",       "T6", "AML.mRNA.2018.all_samples",    "8/15/2017",
                            "COHP_20923",     "3357",  "PBMC",       "T6", "AML.mRNA.2018.all_samples",    "8/15/2017",
                            "COHP_20882",     "3357",  "PBMC",       "T7", "AML.mRNA.2018.all_samples",    "9/12/2017",
                            "COHP_20982",     "3368",  "PBMC",       "T6", "AML.mRNA.2018.all_samples",    "8/15/2017",
                            "COHP_21001",     "3368",  "PBMC",       "T7", "AML.mRNA.2018.all_samples",    "9/12/2017",
                            "COHP_20965",     "3370",  "PBMC",       "T6", "AML.mRNA.2018.all_samples",    "8/15/2017",
                            "COHP_20985",     "3370",  "PBMC",       "T7", "AML.mRNA.2018.all_samples",    "9/12/2017",
                            "COHP_21006",     "3339",  "PBMC",       "T6", "AML.mRNA.2018.all_samples",    "8/15/2017",
                            "COHP_20986",     "3339",  "PBMC",       "T7", "AML.mRNA.2018.all_samples",    "9/12/2017",
                            "COHP_20907",     "3342",  "PBMC",       "T6", "AML.mRNA.2018.all_samples",    "8/15/2017",
                            "COHP_21008",     "3342",  "PBMC",       "T7", "AML.mRNA.2018.all_samples",    "9/12/2017",
                            "COHP_20970",     "3346",  "PBMC",       "T6", "AML.mRNA.2018.all_samples",    "8/15/2017",
                            "COHP_20910",     "3346",  "PBMC",       "T7", "AML.mRNA.2018.all_samples",    "9/12/2017",
                            "COHP_20950",     "3349",  "PBMC",       "T6", "AML.mRNA.2018.all_samples",    "8/15/2017",
                            "COHP_20971",     "3349",  "PBMC",       "T7", "AML.mRNA.2018.all_samples",    "9/12/2017",
                            "COHP_20972",     "3334",  "PBMC",       "T6", "AML.mRNA.2018.all_samples",    "8/15/2017",
                            "COHP_20993",     "3335",  "PBMC",       "T7", "AML.mRNA.2018.all_samples",    "9/12/2017",
                            "COHP_20894",     "3340",  "PBMC",       "T6", "AML.mRNA.2018.all_samples",    "8/15/2017",
                            "COHP_20915",     "3340",  "PBMC",       "T7", "AML.mRNA.2018.all_samples",    "9/12/2017",
                            "COHP_20916",     "3350",  "PBMC",       "T6", "AML.mRNA.2018.all_samples",    "8/15/2017",
                            "COHP_20956",     "3350",  "PBMC",       "T7", "AML.mRNA.2018.all_samples",    "9/12/2017",
                            "COHP_41243",     "3694",  "PBMC",       "T6",             "AML.mRNA.2020",    "5/18/2018",
                            "COHP_41243",     "3694",  "PBMC",       "T6",             "AML.mRNA.2020",    "5/18/2018",
                            "COHP_41244",     "3694",  "PBMC",       "T7",             "AML.mRNA.2020",    "5/25/2018",
                            "COHP_41244",     "3694",  "PBMC",       "T7",             "AML.mRNA.2020",    "5/25/2018",
                            "COHP_41251",     "3695",  "PBMC",       "T6",             "AML.mRNA.2020",    "5/18/2018",
                            "COHP_41251",     "3695",  "PBMC",       "T6",             "AML.mRNA.2020",    "5/18/2018",
                            "COHP_41252",     "3695",  "PBMC",       "T7",             "AML.mRNA.2020",    "5/25/2018",
                            "COHP_41252",     "3695",  "PBMC",       "T7",             "AML.mRNA.2020",    "5/25/2018",
                            "COHP_41259",     "3696",  "PBMC",       "T6",             "AML.mRNA.2020",    "5/18/2018",
                            "COHP_41259",     "3696",  "PBMC",       "T6",             "AML.mRNA.2020",    "5/18/2018",
                            "COHP_41260",     "3696",  "PBMC",       "T7",             "AML.mRNA.2020",    "5/25/2018",
                            "COHP_41260",     "3696",  "PBMC",       "T7",             "AML.mRNA.2020",    "5/25/2018",
                            "COHP_41267",     "3697",  "PBMC",       "T6",             "AML.mRNA.2020",    "5/18/2018",
                            "COHP_41267",     "3697",  "PBMC",       "T6",             "AML.mRNA.2020",    "5/18/2018",
                            "COHP_41268",     "3697",  "PBMC",       "T7",             "AML.mRNA.2020",    "5/25/2018",
                            "COHP_41268",     "3697",  "PBMC",       "T7",             "AML.mRNA.2020",    "5/25/2018",
                            "COHP_41275",     "3698",  "PBMC",       "T6",             "AML.mRNA.2020",    "5/18/2018",
                            "COHP_41275",     "3698",  "PBMC",       "T6",             "AML.mRNA.2020",    "5/18/2018",
                            "COHP_41276",     "3698",  "PBMC",       "T7",             "AML.mRNA.2020",    "5/25/2018",
                            "COHP_41276",     "3698",  "PBMC",       "T7",             "AML.mRNA.2020",    "5/25/2018",
                            "COHP_41283",     "3700",  "PBMC",       "T6",             "AML.mRNA.2020",    "5/18/2018",
                            "COHP_41283",     "3700",  "PBMC",       "T6",             "AML.mRNA.2020",    "5/18/2018",
                            "COHP_41284",     "3700",  "PBMC",       "T7",             "AML.mRNA.2020",    "5/25/2018",
                            "COHP_41284",     "3700",  "PBMC",       "T7",             "AML.mRNA.2020",    "5/25/2018",
                            "COHP_41291",     "3701",  "PBMC",       "T6",             "AML.mRNA.2020",    "5/18/2018",
                            "COHP_41291",     "3701",  "PBMC",       "T6",             "AML.mRNA.2020",    "5/18/2018",
                            "COHP_41292",     "3701",  "PBMC",       "T7",             "AML.mRNA.2020",    "5/25/2018",
                            "COHP_41292",     "3701",  "PBMC",       "T7",             "AML.mRNA.2020",    "5/25/2018",
                            "COHP_41299",     "3702",  "PBMC",       "T6",             "AML.mRNA.2020",    "5/18/2018",
                            "COHP_41299",     "3702",  "PBMC",       "T6",             "AML.mRNA.2020",    "5/18/2018",
                            "COHP_41300",     "3702",  "PBMC",       "T7",             "AML.mRNA.2020",    "5/25/2018",
                            "COHP_41300",     "3702",  "PBMC",       "T7",             "AML.mRNA.2020",    "5/25/2018",
                            "COHP_41307",     "3706",  "PBMC",       "T6",             "AML.mRNA.2020",    "5/18/2018",
                            "COHP_41307",     "3706",  "PBMC",       "T6",             "AML.mRNA.2020",    "5/18/2018",
                            "COHP_41308",     "3706",  "PBMC",       "T7",             "AML.mRNA.2020",    "5/25/2018",
                            "COHP_41308",     "3706",  "PBMC",       "T7",             "AML.mRNA.2020",    "5/25/2018",
                            "COHP_41315",     "3709",  "PBMC",       "T6",             "AML.mRNA.2020",    "5/18/2018",
                            "COHP_41315",     "3709",  "PBMC",       "T6",             "AML.mRNA.2020",    "5/18/2018",
                            "COHP_41316",     "3709",  "PBMC",       "T7",             "AML.mRNA.2020",    "5/25/2018",
                            "COHP_41316",     "3709",  "PBMC",       "T7",             "AML.mRNA.2020",    "5/25/2018"
                            ) |>
                    dplyr::mutate(sample_date = lubridate::mdy(sample_date))

        aml_sample_metadata <- rbind(t6_t7_dates, missing_t6_t7_dates) |>
            dplyr::mutate(mouse_tp = glue::glue("{mouse_id}_{timepoint}_{tissue}"))
        # Get the per-mouse metadata for aml mice, emailed from Ya-Huei Kuo, 2022-9-2
        
        aml_mouse_metadata <- tibble::tribble(
                                    ~mouse_id, ~treatment, ~genotype, ~sex,       ~dob,       ~dod,
                                        "2683",     "Ctrl",      "WT",  "M",   "6/2/2015",  "2/23/2016",
                                        "2684",       "CM",     "CHW",  "M",   "6/2/2015",  "1/19/2016",
                                        "2685",       "CM",     "CHW",  "M",   "6/2/2015",  "2/23/2016",
                                        "2686",       "CM",     "CHW",  "M",   "6/2/2015",  "2/23/2016",
                                        "2689",     "Ctrl",      "CW",  "F",   "6/2/2015",  "2/23/2016",
                                        "2690",       "CM",     "CHW",  "F",   "6/2/2015",  "12/2/2015",
                                        "2692",     "Ctrl",      "CW",  "F",   "6/2/2015",  "2/23/2016",
                                        "2700",     "Ctrl",      "CW",  "F",   "6/3/2015",  "2/23/2016",
                                        "2702",     "Ctrl",      "CW",  "F",   "6/3/2015",  "2/23/2016",
                                        "2705",     "Ctrl",      "HW",  "M",   "6/5/2015",  "2/23/2016",
                                        "2708",       "CM",     "CHW",  "F",   "6/5/2015",  "12/2/2015",
                                        "2709",       "CM",     "CHW",  "F",   "6/5/2015",  "2/23/2016",
                                        "2718",       "CM",     "CHW",  "M",   "6/5/2015",  "2/24/2016",
                                        "2719",       "CM",     "CHW",  "M",   "6/5/2015",  "12/2/2015",
                                        "2720",     "Ctrl",      "HW",  "M",   "6/5/2015",  "2/23/2016",
                                        "2731",       "CM",     "CHW",  "F",  "6/12/2015", "10/24/2015",
                                        "3127",     "Ctrl",      "CW",  "M",   "4/1/2016",   "2/8/2017",
                                        "3130",       "CM",     "CHW",  "F",   "4/1/2016", "11/11/2016",
                                        "3131",       "CM",     "CHW",  "F",   "4/1/2016",  "12/1/2016",
                                        "3200",       "CM",     "CHW",  "M",   "6/7/2016",   "1/4/2017",
                                        "3202",     "Ctrl",      "WT",  "M",   "6/7/2016",   "2/8/2017",
                                        "3334",       "CM",     "CHW",  "M", "10/31/2016",   "9/4/2017",
                                        "3335",     "Ctrl",      "HW",  "M", "10/31/2016",  "11/1/2017",
                                        "3336",       "CM",     "CHW",  "F",  "11/1/2016",  "10/6/2017",
                                        "3338",       "CM",     "CHW",  "F",  "11/2/2016",   "7/3/2017",
                                        "3339",     "Ctrl",      "CW",  "M", "10/29/2016",  "11/1/2017",
                                        "3340",     "Ctrl",      "HW",  "M", "10/29/2016",  "11/1/2017",
                                        "3341",       "CM",     "CHW",  "F", "10/29/2016",  "8/17/2017",
                                        "3342",     "Ctrl",      "CW",  "M", "11/10/2016",  "11/1/2017",
                                        "3346",     "Ctrl",      "CW",  "M", "11/12/2016",  "11/1/2017",
                                        "3349",     "Ctrl",      "CW",  "F", "11/12/2016",  "11/1/2017",
                                        "3350",     "Ctrl",      "WT",  "F", "11/12/2016",  "11/1/2017",
                                        "3357",       "CM",     "CHW",  "F", "11/15/2016",  "10/6/2017",
                                        "3368",       "CM",     "CHW",  "M",  "12/7/2016",  "11/1/2017",
                                        "3370",       "CM",     "CHW",  "M",  "12/8/2016",  "11/1/2017",
                                        "3694",     "Ctrl",      "WT",  "M",   "2/1/2018",  "8/16/2018",
                                        "3695",     "Ctrl",      "WT",  "M",   "2/1/2018",  "8/16/2018",
                                        "3696",     "Ctrl",      "WT",  "M",   "2/1/2018",  "8/16/2018",
                                        "3697",      "3Gy",      "WT",  "F",   "2/1/2018",  "8/16/2018",
                                        "3698",      "3Gy",      "WT",  "F",   "2/1/2018",  "8/16/2018",
                                        "3700",      "3Gy",      "WT",  "F",   "2/1/2018",  "8/16/2018",
                                        "3701",      "3Gy",      "WT",  "M",   "2/1/2018",  "8/16/2018",
                                        "3702",      "3Gy",      "WT",  "M",   "2/1/2018",  "8/16/2018",
                                        "3706",     "Ctrl",      "WT",  "F",   "2/1/2018",  "8/16/2018",
                                        "3709",     "Ctrl",      "WT",  "F",   "2/1/2018",  "8/16/2018",
                                        "4309",       "CM",     "CHW",  "F", "11/10/2019",  "9/29/2020",
                                        "4311",     "Ctrl",      "CW",  "F", "11/10/2019",  "7/12/2021",
                                        "4321",       "CM",     "CHW",  "F", "11/13/2019",  "7/17/2020",
                                        "4324",       "CM",     "CHW",  "M", "11/13/2019",  "7/31/2020",
                                        "4329",       "CM",     "CHW",  "F", "11/19/2019",   "9/9/2020",
                                        "4419",       "CM",     "CHW",  "F",  "7/14/2020",   "4/2/2021",
                                        "4433",       "CM",     "CHW",  "M",  "8/21/2020",  "4/19/2021",
                                        "4436",       "CM",     "CHW",  "F",  "8/21/2020",  "3/21/2021",
                                        "4443",       "CM",     "CHW",  "F",  "8/25/2020",  "5/26/2021",
                                        "4498",     "Ctrl",      "CW",  "M",   "3/3/2021", "10/13/2021",
                                        "4501",     "Ctrl",      "CW",  "F",   "3/3/2021",  "12/7/2021",
                                        "4502",       "CM",     "CHW",  "M",   "3/3/2021", "10/13/2021",
                                        "4506",       "CM",     "CHW",  "F",   "3/3/2021", "11/17/2021",
                                        "4510",       "CM",     "CHW",  "F",  "3/11/2021",  "11/8/2021",
                                        "4512",     "Ctrl",      "WT",  "M",  "3/25/2021",  "10/6/2021",
                                        "4520",     "Ctrl",      "CW",  "F",   "4/9/2021",  "9/29/2021",
                                        "4521",       "CM",     "CHW",  "M",  "4/26/2021",  "10/6/2021",
                                        "4522",       "CM",     "CHW",  "F",  "4/26/2021",  "9/29/2021",
                                        "4534",     "Ctrl",      "WT",  "F",  "5/22/2021",  "12/7/2021",
                                        "4535",       "CM",     "CHW",  "F",  "5/22/2021",  "12/7/2021"
                                    ) |>
                    dplyr::mutate(
                        dob = lubridate::mdy(dob),
                        dod = lubridate::mdy(dod))

    aml_mice <- dplyr::left_join(aml_sample_metadata, aml_mouse_metadata, by = "mouse_id")

    # CML sample and mouse metadata
    # download the file "General/Copy of mmu_metadata_CML.xlsx"
    get_teams_file("General/Copy of mmu_metadata_CML.xlsx")

    # read in the worksheets
    cml_mice <- left_join(
        # Samples
        readxl::read_excel("data-raw/Copy of mmu_metadata_CML.xlsx",
            sheet = "CML samples", col_types = "text"
        ) |>
            dplyr::group_by(mouse_id) |>
            dplyr::mutate(start_date = as.Date(as.numeric(sample_date), origin = "1899-12-30")) |>
            tidyr::fill(start_date) |>
            dplyr::mutate(sample_date = ((7 * as.numeric(timepoint)) + start_date)),
        # Mice
        readxl::read_excel("data-raw/Copy of mmu_metadata_CML.xlsx",
            sheet = "CML mice", col_types = "text"
        ) |>
            dplyr::mutate(
                dob = as.Date(as.numeric(dob), origin = "1899-12-30"),
                dod = as.Date(as.numeric(dod), origin = "1899-12-30")
            )) |>
        dplyr::mutate(
            mouse_tp = glue::glue("{mouse_id}_{timepoint}_{tissue}")) |>
            dplyr::select(-c(start_date, `...1`))

    # consolidate sample metadata where possible
    # using dplyr::rows_update() & tidyr::fill() to fill in missing values
    sample_sheet <- dplyr::rows_update(all_mice, aml_mice, by = "mouse_tp", unmatched = "ignore") |>
        dplyr::rows_update(cml_mice, by = "mouse_tp", unmatched = "ignore") |>
        dplyr::group_by(mouse_id) |>
        tidyr::fill(c("treatment", "genotype", "sex", "dob", "dod"), .direction = "downup") |>
        dplyr::mutate(
            # # Use "dod" for "sample_date" for some samples
            # sample_date = dplyr::case_when(
            #     tissue == "BM" ~ dod,
            #     timepoint == "L" ~ dod,
            #     timepoint == "END" ~ dod,
            #     TRUE ~ sample_date
            # ),
            # add columns: sample_weeks, age_at_end, age_at_start, age_at_sample
            sample_weeks = difftime(sample_date, min(sample_date), units = "weeks"),
            age_at_end = difftime(max(sample_date), dob, units = "weeks"),
            age_at_start = difftime(min(sample_date), dob, units = "weeks"),
            age_at_sample = difftime(sample_date, dob, units = "weeks"),
            dplyr::across(dplyr::starts_with(c("age", "sample_weeks")), round, 1)
        ) |>
        dplyr::ungroup() |>
        dplyr::select(-c(mouse_tp))

    return(sample_sheet)
}

#### mRNA -----

# ├ AML.mRNA.2016
parse_metadata_AML.mRNA.2016 <- function() {
    project <- "AML.mRNA.2016"
    sample_sheet <- read.csv("/net/nfs-irwrsrchnas01/labs/rrockne/MHO/AML.mRNA.2016/config/new_name_key.tsv", sep = "\t") |>
        dplyr::mutate(
            project = project,
            fastq_1 = paste0("/net/nfs-irwrsrchnas01/labs/rrockne/MHO/AML.mRNA.2016/data/", newname, ".gz"),
            library_id = paste0("COHP_", ID),
            tissue = "PBMC",
            fastq_2 = "",
            dob = NA,
            Batch = paste0("2016_", Batch),
            strandedness = dplyr::case_when(
                stringr::str_detect(ID, "11548") ~ "reverse",
                stringr::str_detect(TimePoint, "L|T0|T1") ~ "unstranded",
                TRUE ~ "reverse"
            )
        ) |>
        dplyr::select(sample = library_id, fastq_1, fastq_2, strandedness, mouse_id = Mouse, tissue, timepoint = TimePoint, batch = Batch, treatment = Treatment, genotype = Genotype, sex = Sex, dob, project)
    return(sample_sheet)
}

# AML.mRNA.2018.all_samples
parse_metadata_AML.mRNA.2018.all_samples <- function() {
    project <- "AML.mRNA.2018.all_samples"
    sample_sheet <- read.csv("/net/nfs-irwrsrchnas01/labs/rrockne/MHO/AML.mRNA.2018.all_samples/config/new_name_key.tsv", sep = "\t") |>
        dplyr::mutate(
            project = project,
            fastq_1 = paste0("/net/nfs-irwrsrchnas01/labs/rrockne/MHO/AML.mRNA.2018.all_samples/data/fastq/", Newname, ".gz"),
            library_id = paste0("COHP_", gsub(".*_|\\.fq", "", Newname)),
            tissue = "PBMC",
            fastq_2 = "",
            strandedness = "reverse",
            dob = NA
        ) |>
        tidyr::separate(Newname, sep = "_", into = c("timepoint", "mouse_id", NA, "treatment", "genotype", "sex", "batch", NA)) |>
        dplyr::mutate(treatment = case_when(
            treatment == "CTL" ~ "Ctrl",
            TRUE ~ treatment
        )) |>
        dplyr::select(sample = library_id, fastq_1, fastq_2, strandedness, mouse_id, tissue, timepoint, batch, treatment, genotype, sex, dob, project) |>
        dplyr::mutate(batch = paste0("2018_", batch))
    return(sample_sheet)
}

# AML.mRNA.2020
parse_metadata_AML.mRNA.2020 <- function() {
    project <- "AML.mRNA.2020"
    # Copied from "Seq Samples Dates.xlsx":  https://github.com/drejom/haemdata/issues/6
    sample_data <- tibble::tribble(
        ~mouse_id, ~sex, ~genotype, ~treatment, ~dob,
        "3694", "M", "WT", "Ctrl", NA,
        "3695", "M", "WT", "Ctrl", NA,
        "3696", "M", "WT", "Ctrl", NA,
        "3697", "F", "WT", "3Gy", NA,
        "3698", "F", "WT", "3Gy", NA,
        "3700", "F", "WT", "3Gy", NA,
        "3701", "M", "WT", "3Gy", NA,
        "3702", "M", "WT", "3Gy", NA,
        "3706", "F", "WT", "Ctrl", NA,
        "3709", "F", "WT", "Ctrl", NA
    )
    fastqs <- data.frame(
        fastq_1 = system(
            paste0("find '/net/nfs-irwrsrchnas01/labs/ykuo/Seq/201124_TGen' -name '*.gz'"),
            intern = TRUE
        ) %>% grep("_R1_", ., value = TRUE),
        fastq_2 = system(
            paste0("find '/net/nfs-irwrsrchnas01/labs/ykuo/Seq/201124_TGen' -name '*.gz'"),
            intern = TRUE
        ) %>% grep("_R2_", ., value = TRUE)
    ) |> dplyr::mutate(library_id = stringr::str_extract(fastq_1, "COHP_\\d{5}"))
    sample_sheet <- read.csv("/net/nfs-irwrsrchnas01/labs/rrockne/MHO/AML.mRNA.2020/config/rename_fastqs.tsv", sep = "\t", header = FALSE) |>
        tidyr::separate(V2, sep = "_", c("group", "mouse_id", "timepoint")) |>
        dplyr::mutate(
            project = project,
            library_id = V1,
            tissue = "PBMC",
            strandedness = "reverse",
            batch = "2020_A",
            tissue = dplyr::case_when(
                group == "BM" ~ "BM",
                TRUE ~ tissue
            ),
            timepoint = dplyr::case_when(
                tissue == "BM" ~ NA_character_,
                TRUE ~ timepoint
            )
        ) |>
        dplyr::left_join(sample_data) |>
        dplyr::left_join(fastqs) |>
        dplyr::select(sample = library_id, fastq_1, fastq_2, strandedness, mouse_id, tissue, timepoint, batch, treatment, genotype, sex, dob, project)
    return(sample_sheet)
}

# AML.mRNA.2021.RxGroup1
parse_metadata_AML.mRNA.2021.RxGroup1 <- function() {
    project <- "AML.mRNA.2021.RxGroup1"
    fastqs <- data.frame(
        fastq_1 = system(
            paste0("find '/net/nfs-irwrsrchnas01/labs/rrockne/Seq/210412_Tgen' -name '*.gz'"),
            intern = TRUE
        ) %>% grep("_R1_", ., value = TRUE),
        fastq_2 = system(
            paste0("find '/net/nfs-irwrsrchnas01/labs/rrockne/Seq/210412_Tgen' -name '*.gz'"),
            intern = TRUE
        ) %>% grep("_R2_", ., value = TRUE)
    ) |> dplyr::mutate(library_id = stringr::str_extract(fastq_1, "COHP_\\d{5}"))
    xls <- "/net/nfs-irwrsrchnas01/labs/rrockne/MHO/AML.mRNA.2021.RxGroup1/config/run_summary_IGC-LZ-18757.xlsx"
    sample_sheet <- readxl::read_excel(xls) |>
        dplyr::mutate(
            project = project,
            library_id = TGen_Sample_Name,
            mouse_id = gsub("_.*", "", name),
            tissue = "PBMC",
            timepoint = gsub("^.*_", "", name),
            batch = "2021_A",
            strandedness = "reverse",
            genotype = NA_character_,
            tissue = dplyr::case_when(
                timepoint == "BM" ~ "BM",
                TRUE ~ tissue
            ),
            timepoint = sub("BM", NA_character_, timepoint),
            timepoint = dplyr::case_when(
                is.na(timepoint) ~ timepoint,
                !grepl("\\D", timepoint) ~ paste0("W", timepoint),
                TRUE ~ timepoint
            ),
        ) |>
        dplyr::left_join(fastqs)
    xls_cm <- "data-raw/CM mice chemo Rx survival.xlsx"
    sample_sheet <- dplyr::left_join(sample_sheet, readxl::read_excel(xls_cm) |>
        dplyr::select(
            mouse_id = ID,
            treatment = Mice,
            sex = Gender,
            dob = DOB
        ) |>
        dplyr::mutate(mouse_id = as.character(mouse_id))) |>
        dplyr::select(sample = library_id, fastq_1, fastq_2, strandedness, mouse_id, tissue, timepoint, batch, treatment, genotype, sex, dob, project)
    return(sample_sheet)
}

# AML.mRNA.2021.RxGroup2
parse_metadata_AML.mRNA.2021.RxGroup2 <- function() {
    project <- "AML.mRNA.2021.RxGroup2"
    fastqs <- data.frame(
        fastq_1 = system(
            paste0("find '/net/nfs-irwrsrchnas01/labs/rrockne/Seq/210507_Tgen' -name '*.gz'"),
            intern = TRUE
        ) %>% grep("_R1_", ., value = TRUE),
        fastq_2 = system(
            paste0("find '/net/nfs-irwrsrchnas01/labs/rrockne/Seq/210507_Tgen' -name '*.gz'"),
            intern = TRUE
        ) %>% grep("_R2_", ., value = TRUE)
    ) |> dplyr::mutate(library_id = stringr::str_extract(fastq_1, "COHP_\\d{5}"))
    xls <- "/net/nfs-irwrsrchnas01/labs/rrockne/MHO/AML.mRNA.2021.RxGroup2/config/sequencing summary_IGC-LZ-18862.xlsx"
    sample_sheet <- readxl::read_excel(xls) |>
        dplyr::mutate(
            project = project,
            library_id = TGen_Sample_Name,
            Sample_ID = gsub("-", "", Sample_ID),
            Sample_ID = gsub("PBend", "END", Sample_ID),
            mouse_id = gsub(".*_", "", Sample_ID) |> substr(1, 4),
            tissue = "PBMC",
            timepoint = gsub(".*_\\d{4}", "", Sample_ID),
            batch = "2021_B",
            strandedness = "reverse",
            genotype = NA_character_,
            tissue = dplyr::case_when(
                timepoint == "BM" ~ "BM",
                TRUE ~ tissue
            ),
            timepoint = sub("BM", NA_character_, timepoint)
        ) |>
        dplyr::left_join(fastqs)
    xls_cm <- "data-raw/CM mice chemo Rx survival.xlsx"
    sample_sheet <- dplyr::left_join(sample_sheet, readxl::read_excel(xls_cm) |>
        dplyr::select(
            mouse_id = ID,
            treatment = Mice,
            sex = Gender,
            dob = DOB
        ) |>
        dplyr::mutate(mouse_id = as.character(mouse_id))) |>
        dplyr::select(sample = library_id, fastq_1, fastq_2, strandedness, mouse_id, tissue, timepoint, batch, treatment, genotype, sex, dob, project)
    return(sample_sheet)
}
# AML.mRNA.2021.RxGroup2_pt2
parse_metadata_AML.mRNA.2021.RxGroup2_pt2 <- function() {
    project <- "AML.mRNA.2021.RxGroup2_pt2"
    fastqs <- data.frame(
        fastq_1 = system(
            paste0("find '/net/nfs-irwrsrchnas01/labs/ykuo/Seq/210910_TGen' -name '*.gz'"),
            intern = TRUE
        ) %>% grep("_R1_", ., value = TRUE),
        fastq_2 = system(
            paste0("find '/net/nfs-irwrsrchnas01/labs/ykuo/Seq/210910_TGen' -name '*.gz'"),
            intern = TRUE
        ) %>% grep("_R2_", ., value = TRUE)
    ) |> dplyr::mutate(library_id = stringr::str_extract(fastq_1, "COHP_\\d{5}"))
    sample_sheet <- read.csv("/net/nfs-irwrsrchnas01/labs/rrockne/MHO/AML.mRNA.2021.RxGroup2_pt2/config/rename_fastqs.tsv", sep = "\t", header = FALSE) |>
        dplyr::filter(stringr::str_detect(V2, "_R1\\.")) |>
        dplyr::mutate(
            project = project,
            library_id = stringr::str_extract(V1, "COHP_\\d{5}"),
            Sample_ID = gsub("_.*", "", V2),
            mouse_id = gsub("-.*", "", V2),
            tissue = "PBMC",
            timepoint = gsub(".*-", "", Sample_ID),
            batch = "2021_C",
            strandedness = "reverse",
            genotype = NA_character_,
            tissue = dplyr::case_when(
                timepoint == "BM" ~ "BM",
                TRUE ~ tissue
            ),
            timepoint = sub("BM", NA_character_, timepoint)
        ) |>
        dplyr::select(-V1) |>
        dplyr::left_join(fastqs)
    xls_cm <- "data-raw/CM mice chemo Rx survival.xlsx"
    sample_sheet <- dplyr::left_join(sample_sheet, readxl::read_excel(xls_cm) |>
        dplyr::select(
            mouse_id = ID,
            treatment = Mice,
            sex = Gender,
            dob = DOB
        ) |>
        dplyr::mutate(mouse_id = as.character(mouse_id))) |>
        dplyr::select(sample = library_id, fastq_1, fastq_2, strandedness, mouse_id, tissue, timepoint, batch, treatment, genotype, sex, dob, project)
    return(sample_sheet)
}
# AML.mRNA.2022.RxGroup3
parse_metadata_AML.mRNA.2022.RxGroup3 <- function() {
    project <- "AML.mRNA.2022.RxGroup3"
    fastqs <- data.frame(
        fastq_1 = system(
            paste0("find '/net/nfs-irwrsrchnas01/labs/ykuo/Seq/220415_IGC-LZ-20205' -name '*.gz'"),
            intern = TRUE
        ) %>% grep("_R1_", ., value = TRUE),
        fastq_2 = system(
            paste0("find '/net/nfs-irwrsrchnas01/labs/ykuo/Seq/220415_IGC-LZ-20205' -name '*.gz'"),
            intern = TRUE
        ) %>% grep("_R2_", ., value = TRUE)
    ) |> dplyr::mutate(library_id = stringr::str_extract(fastq_1, "COHP_\\d{5}"))

    xls <- "data-raw/sample summary_IGC-LZ-20205.xlsx"
    sample_sheet <- readxl::read_excel(xls) |>
        dplyr::transmute(
            sample_name = gsub(".*_", "", Sample_ID),
            library_id = paste0("COHP_", gsub("_.*", "", Sample_ID))
        ) |>
        dplyr::mutate(
            project = project,
            strandedness = "reverse",
            tissue = "PBMC",
            genotype = NA,
            batch = "2021_D",
            sample_name = gsub("-", "_", sample_name),
            sample_name = dplyr::case_when(
                sample_name == "T1" ~ "4534_T1",
                sample_name == "T2" ~ "4534_T2",
                stringr::str_detect(sample_name, "\\d{4}w\\d") ~ gsub("w", "_W", sample_name),
                stringr::str_detect(sample_name, "\\d{4}end") ~ gsub("end", "_END", sample_name),
                TRUE ~ sample_name
            ),
            sample_name = gsub("_2$", "", sample_name)
        ) |>
        tidyr::separate(sample_name, sep = "_", into = c("mouse_id", "timepoint")) |>
        dplyr::mutate(
            tissue = dplyr::case_when(
                timepoint == "BM" ~ "BM",
                TRUE ~ tissue
            ),
            timepoint = gsub("BM", NA, timepoint)
        ) |>
        dplyr::left_join(fastqs)
    xls_cm <- "data-raw/CM mice chemo Rx survival.xlsx"
    sample_sheet <- dplyr::left_join(sample_sheet, readxl::read_excel(xls_cm) |>
        dplyr::select(
            mouse_id = ID,
            treatment = Mice,
            sex = Gender,
            dob = DOB
        ) |>
        dplyr::mutate(mouse_id = as.character(mouse_id))) |>
        dplyr::select(sample = library_id, fastq_1, fastq_2, strandedness, mouse_id, tissue, timepoint, batch, treatment, genotype, sex, dob, project)
    return(sample_sheet)
}
# AML.mRNA.novaseq_validation.2020
parse_metadata_AML.mRNA.novaseq_validation.2020 <- function() {
    project <- "AML.mRNA.novaseq_validation.2020"
    ## fastqs generated at COH
    # fastqs <- data.frame(
    #     fastq_1 = system(
    #         paste0("find '/net/nfs-irwrsrchnas01/labs/ykuo/Seq/200829' -name '*.gz'"),
    #         intern = TRUE
    #     ) %>% grep("_R1_", ., value = TRUE),
    #     fastq_2 = system(
    #         paste0("find '/net/nfs-irwrsrchnas01/labs/ykuo/Seq/200829' -name '*.gz'"),
    #         intern = TRUE
    #     ) %>% grep("_R2_", ., value = TRUE)
    # ) |> dplyr::mutate(library_id = paste0("IGCP_", substr(basename(fastq_1), 1, 5)))
    ## fastqs generated at TGen
    fastqs <- data.frame(
        fastq_1 = system(
            paste0("find '/net/nfs-irwrsrchnas01/labs/ykuo/Seq/200928_TGen' -name '*.gz'"),
            intern = TRUE
        ) %>% grep("_R1_", ., value = TRUE),
        fastq_2 = system(
            paste0("find '/net/nfs-irwrsrchnas01/labs/ykuo/Seq/200928_TGen' -name '*.gz'"),
            intern = TRUE
        ) %>% grep("_R2_", ., value = TRUE)
    ) |> dplyr::mutate(library_id = stringr::str_extract(fastq_1, "COHP_\\d{5}"))
    # https://github.com/drejom/haemdata/issues/5
    mouse_data <- data.frame(sample = c("COHP_38491", "COHP_38493", "COHP_38492"), mouse_id = c("3341", "3341", "3341"), timepoint = c("T2", "T5", "T6"))

    sample_sheet <- fastqs |>
        dplyr::mutate(
            project = project,
            sample = library_id,
            tissue = "PBMC",
            batch = "2020_B",
            strandedness = "reverse",
            sex = NA_character_,
            dob = NA,
            treatment = NA_character_,
            genotype = NA_character_
        ) |>
        dplyr::left_join(mouse_data, by = "sample") |>
        dplyr::select(sample = library_id, fastq_1, fastq_2, strandedness, mouse_id, tissue, timepoint, batch, treatment, genotype, sex, dob, project)
    return(sample_sheet)
}
# AML.validation.2017
parse_metadata_AML.validation.2017 <- function() {
    project <- "AML.validation.2017"
    fastqs <- data.frame(
        fastq_1 = system(
            paste0("find '/net/nfs-irwrsrchnas01/labs/rrockne/MHO/AML.validation.2017/data/fastq' -name '*.gz'"),
            intern = TRUE
        ) %>% grep("_R1_", ., value = TRUE)
    ) |> dplyr::mutate(library_id = paste0("COHP_", substr(basename(fastq_1), 1, 5)))

    xls <- "/net/nfs-irwrsrchnas01/labs/rrockne/MHO/AML.validation.2017/data/validation_name_key.xlsx"
    sample_sheet <- readxl::read_excel(xls) |>
        dplyr::mutate(
            project = project,
            library_id = paste0("COHP_", substr(Filename, 1, 5)),
            mouse_id = Mouse,
            tissue = "PBMC",
            timepoint = TimePoint,
            batch = "2017_A",
            strandedness = "reverse",
            sex = Sex,
            dob = NA,
            fastq_2 = "",
            treatment = Treatment,
            genotype = Genotype
        ) |>
        dplyr::left_join(fastqs) |>
        dplyr::select(sample = library_id, fastq_1, fastq_2, strandedness, mouse_id, tissue, timepoint, batch, treatment, genotype, sex, dob, project)
    return(sample_sheet)
}
# AML.scRNAseq.2022
parse_metadata_AML.scRNAseq.2022 <- function() {
    project <- "AML.scRNAseq.2022"
    fastqs <- data.frame(
        fastq_1 = system(
            paste0("find '/net/nfs-irwrsrchnas01/labs/rrockne/Seq/211210_IGC-LZ-19773' -name '*.gz'"),
            intern = TRUE
        ) %>% grep("_R1_", ., value = TRUE),
        fastq_2 = system(
            paste0("find '/net/nfs-irwrsrchnas01/labs/rrockne/Seq/211210_IGC-LZ-19773' -name '*.gz'"),
            intern = TRUE
        ) %>% grep("_R2_", ., value = TRUE)
    ) |> dplyr::mutate(library_id = stringr::str_extract(fastq_1, "COHP_\\d{5}"))
    sample_data <- tibble::tribble(
        ~mouse_id, ~treatment, ~genotype, ~sex, ~dob,
        "4520", "Ctrl", "C", "F", "2021-04-09",
        "4522", "CM", "CHW", "F", "2021-04-26",
        "4512", "Ctrl", "WT", "M", "2021-03-25",
        "4521", "CM", "CHW", "M", "2021-04-26",
        "4498", "Ctrl", "C", "M", "2021-03-03",
        "4502", "CM", "CHW", "M", "2021-03-03"
    )
    xls <- "/net/nfs-irwrsrchnas01/labs/rrockne/MHO/AML.scRNA.2022/bulkseq_targets/config/sequencing summary_IGC-LZ-19773.xlsx"
    sample_sheet <- readxl::read_excel(xls) |>
        dplyr::mutate(
            project = project,
            library_id = paste0("COHP_", gsub("_.*", "", Sample_Name)),
            mouse_id = gsub(".*_", "", Sample_Name) |> substr(1, 4),
            timepoint = "T0",
            batch = "2022_B",
            strandedness = "reverse",
            tissue = gsub(".*_\\d{4}", "", Sample_Name),
            tissue = dplyr::case_when(
                tissue == "BM" ~ "BM",
                stringr::str_detect(Sample_Name, "CKIT") ~ "BM_CKIT",
                TRUE ~ "PBMC"
            )
        ) |>
        dplyr::left_join(fastqs) |>
        dplyr::left_join(sample_data) |>
        dplyr::select(sample = library_id, fastq_1, fastq_2, strandedness, mouse_id, tissue, timepoint, batch, treatment, genotype, sex, dob, project)
    return(sample_sheet)
}
# CML.mRNA.2021
parse_metadata_CML.mRNA.2021 <- function() {
    cohort <- "CML.mRNA.2021"
    fastqs <- data.frame(
        fastq_1 = system(
            paste0("find /net/nfs-irwrsrchnas01/labs/rrockne/Seq/210910_TGen /net/nfs-irwrsrchnas01/labs/rrockne/Seq/210907_TGen -name '*.gz'"),
            intern = TRUE
        ) %>% grep("_R1_", ., value = TRUE),
        fastq_2 = system(
            paste0("find /net/nfs-irwrsrchnas01/labs/rrockne/Seq/210910_TGen /net/nfs-irwrsrchnas01/labs/rrockne/Seq/210907_TGen -name '*.gz'"),
            intern = TRUE
        ) %>% grep("_R2_", ., value = TRUE)
    ) |> dplyr::mutate(library_id = stringr::str_extract(fastq_1, "COHP_\\d{5}"))

    xls <- "/net/nfs-irwrsrchnas01/labs/rrockne/MHO/CML.mRNA.2021/config/Lianjun_IGC-LZ-19411.xlsx"
    sample_sheet <- readxl::read_excel(xls) |>
        dplyr::mutate(
            project = cohort,
            library_id = TGen_Sample_Name,
            mouse_id = gsub("-.*", "", Sample_ID),
            timepoint = gsub(".*-", "", Sample_ID),
            batch = "2021_E",
            strandedness = "reverse",
            sex = NA_character_,
            dob = NA,
            treatment = dplyr::case_when(
                stringr::str_detect(mouse_id, regex("^473$|^474$|^475$|^476$|^485$|^482$|^483$|^545$")) ~ "TET_OFF_ON_A",
                stringr::str_detect(mouse_id, regex("^479$|^478$|^480$|^481$|^484$|^477$|^487$|^540$|^547$|^541$|^542$")) ~ "TET_OFF_B",
                stringr::str_detect(mouse_id, regex("^490$|^491$|^492$|^493$|^494$|^488$|^489$|^546$")) ~ "TET_ON_C",
                stringr::str_detect(mouse_id, regex("^502$|^511$|^512$|^513$|^507$|^508$|^505$|^510$|^514$")) ~ "TET_OFF_NIL_ON_D",
                TRUE ~ "Ctrl"
            ),
            genotype = NA_character_,
            tissue = "PBMC"
        ) |>
        dplyr::left_join(fastqs) |>
        dplyr::select(sample = library_id, fastq_1, fastq_2, strandedness, mouse_id, tissue, timepoint, batch, treatment, genotype, sex, dob, project)
    return(sample_sheet)
}
# CML.mRNA.2022
parse_metadata_CML.mRNA.2022 <- function() {
    project <- "CML.mRNA.2022"
    fastqs <- data.frame(
        fastq_1 = system(
            paste0("find /net/nfs-irwrsrchnas01/labs/gmarcucci/Seq/220210_IGC-LZ-19931 /net/nfs-irwrsrchnas01/labs/gmarcucci/Seq/220202_IGC-LZ-19931 -name '*.gz'"),
            intern = TRUE
        ) %>% grep("_R1_", ., value = TRUE),
        fastq_2 = system(
            paste0("find /net/nfs-irwrsrchnas01/labs/gmarcucci/Seq/220210_IGC-LZ-19931 /net/nfs-irwrsrchnas01/labs/gmarcucci/Seq/220202_IGC-LZ-19931 -name '*.gz'"),
            intern = TRUE
        ) %>% grep("_R2_", ., value = TRUE)
    ) |> dplyr::mutate(library_id = stringr::str_extract(fastq_1, "COHP_\\d{5}"))

    sample_sheet <- read.csv("/net/nfs-irwrsrchnas01/labs/rrockne/MHO/CML.mRNA.2022/config/sample_sheet.csv") |>
        dplyr::mutate(
            project = project,
            library_id = stringr::str_extract(fastq_1, "COHP_\\d{5}"),
            mouse_id = gsub("^.(\\d*)-.*", "\\1", sample),
            timepoint = gsub(".*-", "", sample),
            batch = "2022_A",
            strandedness = "reverse",
            sex = substr(sample, 1, 1),
            dob = NA,
            treatment = dplyr::case_when(
                stringr::str_detect(mouse_id, regex("^473$|^474$|^475$|^476$|^485$|^482$|^483$|^545$")) ~ "TET_OFF_ON_A",
                stringr::str_detect(mouse_id, regex("^479$|^478$|^480$|^481$|^484$|^477$|^487$|^540$|^547$|^541$|^542$")) ~ "TET_OFF_B",
                stringr::str_detect(mouse_id, regex("^490$|^491$|^492$|^493$|^494$|^488$|^489$|^546$")) ~ "TET_ON_C",
                stringr::str_detect(mouse_id, regex("^502$|^511$|^512$|^513$|^507$|^508$|^505$|^510$|^514$")) ~ "TET_OFF_NIL_ON_D",
                TRUE ~ "Ctrl"
            ),
            genotype = NA_character_,
            tissue = "PBMC"
        ) |>
        dplyr::left_join(fastqs) |>
        dplyr::select(sample = library_id, fastq_1, fastq_2, strandedness, mouse_id, tissue, timepoint, batch, treatment, genotype, sex, dob, project)
    return(sample_sheet)
}
# CML.mRNA.2022_pt2
parse_metadata_CML.mRNA.2022_pt2 <- function() {
    fastqs <- data.frame(
        fastq_1 = system(
            paste0("find /net/nfs-irwrsrchnas01/labs/gmarcucci/Seq/220811_0471_IGC-LZ-20847 /net/nfs-irwrsrchnas01/labs/gmarcucci/Seq/220811_0472_IGC-LZ-20847 -name '*.gz'"),
            intern = TRUE
        ) %>% grep("_R1_", ., value = TRUE),
        fastq_2 = system(
            paste0("find /net/nfs-irwrsrchnas01/labs/gmarcucci/Seq/220811_0471_IGC-LZ-20847 /net/nfs-irwrsrchnas01/labs/gmarcucci/Seq/220811_0472_IGC-LZ-20847 -name '*.gz'"),
            intern = TRUE
        ) %>% grep("_R2_", ., value = TRUE)
    ) |> dplyr::mutate(library_id = stringr::str_extract(fastq_1, "COHP_\\d{5}"))

    sample_sheet <- readxl::read_xlsx("/home/domeally/workspaces/haemdata/data-raw/sample summary_IGC--LZ-20847.xlsx") |>
        dplyr::mutate(
            project = "CML.mRNA.2022_pt2",
            library_id = paste0("COHP_", gsub("^(\\d{5})_(\\d\\d?)_([MF])(\\d{3}\\d?)-(\\d.*)", "\\1", Sample_ID)),
            mouse_id = gsub("^(\\d{5})_(\\d\\d?)_([MF])(\\d{3}\\d?)-(\\d.*)", "\\4", Sample_ID),
            timepoint = gsub("^(\\d{5})_(\\d\\d?)_([MF])(\\d{3}\\d?)-(\\d.*)", "\\5", Sample_ID),
            batch = "2022_B",
            strandedness = "reverse",
            sex = gsub("^(\\d{5})_(\\d\\d?)_([MF])(\\d{3}\\d?)-(\\d.*)", "\\3", Sample_ID),
            dob = NA,
            treatment = NA_character_,
            genotype = NA_character_,
            tissue = "PBMC"
        ) |>
        dplyr::left_join(fastqs) |>
        dplyr::select(sample = library_id, fastq_1, fastq_2, strandedness, mouse_id, tissue, timepoint, batch, treatment, genotype, sex, dob, project)
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
