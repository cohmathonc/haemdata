# Functions for importing metadata

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
        dplyr::arrange(batch, patient_id)

    return(metadata_hsa)
}

#' Update metadata for mouse AML samples
#'
#' This function prepares the `metadata_mmu` object for all RNAseq libraries from AML and CML mice.
#' Minimal metadata fields include sample_id, library_id, fastq_1, fastq_2, strandedness,
#' mouse_id, tissue, sample_date, week, timepoint, batch, treatment, genotype,
#' sex, dob, dod, cohort.
#'
#' Raise an \href{https://github.com/drejom/haemdata/issues}{issue on GitHub}
#' to report erroneous or missing records.
#'
#' @details This function updates the `metadata_mmu` object for all RNAseq libraries
#' by using a previously published version from a pin, excel file, csv, etc. The metadata
#' table was assembled by crawling a range of disparate datasources including emails,
#' CSVs, excel files, etc. and was originally put together by the make_metadata_mmu()
#' function (which has been retired since version 0.0.0.0.9008).
#'
#' This function is really just a place holder for code that manipulates the metadata table
#' for each release. See previous releases on GihUb to track down how raw metadata was prepared.
#' #'
#' @name update_metadata_mmu
#' @return a data.frame
#' @author Denis O'Meally

update_metadata_mmu <- function() {
    # Update sample and mouse metadata
    # Load in metadata_mmu from the devel pinboard and modify as needed

    # v0.0.0.9011 interim
    sample_sheet <- pins::pin_read(
        pins::board_folder(
            "/labs/rrockne/MHO/haemdata",
            versioned = FALSE
        ),
        "metadata_mmu.csv"
    ) |>
        purrr::modify_if(is.factor, as.character)

    # # # change sample_id prefix to MHO from PSON
    sample_sheet <- sample_sheet |>
        mutate(sample_id = str_replace(sample_id, "PSON", "MHO"))

    # # # update Cellplex file paths to Cellranger v7.0.0
    sample_sheet <- sample_sheet |>
        mutate(
            hdf5 = str_replace(hdf5, "(GENCODEm28_HLT_lib[1-9])/outs", "\\1_v7/outs"),
            hdf5 = str_replace(hdf5, "/net/nfs-irwrsrchnas01/labs", "/labs"),
            hdf5 = str_replace(hdf5, "sample_feature_bc_matrix.h5", "sample_filtered_feature_bc_matrix.h5"))

    # # # simplify fastq_1 & fastq_2 paths
    sample_sheet <- sample_sheet |>
        mutate(
            fastq_1 = str_replace(fastq_1, "/net/nfs-irwrsrchnas01/labs", "/labs"),
            fastq_2 = str_replace(fastq_2, "/net/nfs-irwrsrchnas01/labs", "/labs")
        )

    # # # Save a copy ----
    sample_sheet |>
        rio::export(here::here("data-raw/metadata_mmu.rds"))

    return(sample_sheet)
}
#### MMU mRNA -----
# ├ AML.mRNA.2016
# AML.mRNA.2018.all_samples
# TODO update to use file paths from new isilon - symlinks are broken

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
        dplyr::mutate(library_id = stringr::str_extract(fastq_1, "COHP_\\d{5}"), strandedness = "reverse") |>
        dplyr::select(library_id, fastq_1, fastq_2, strandedness)
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
        dplyr::mutate(library_id = paste0("COHP_", substr(basename(fastq_1), 1, 5)), fastq_2 = "", strandedness = "reverse") |>
        dplyr::select(library_id, fastq_1, fastq_2, strandedness)
    return(sample_sheet)
}

#### HSA mRNA -----

# AML.mRNA.HSA_FLT3.2022
parse_metadata_AML.mRNA.HSA_FLT3.2022 <- function() {
    cohort <- "AML.mRNA.HSA_FLT3.2022"

    fastqs <- data.frame(
        fastq_1 = list.files(c("/net/nfs-irwrsrchnas01/labs/ykuo/Seq/220511_IGC-LZ-20342", "/net/nfs-irwrsrchnas01/labs/ykuo/Seq/220715_IGC-LZ-20821"),
                pattern = "_R1_.*\\.gz", full.names = TRUE, recursive = TRUE
            ),
        fastq_2 = list.files(c("/net/nfs-irwrsrchnas01/labs/ykuo/Seq/220511_IGC-LZ-20342", "/net/nfs-irwrsrchnas01/labs/ykuo/Seq/220715_IGC-LZ-20821"),
            pattern = "_R2_.*\\.gz", full.names = TRUE, recursive = TRUE
        )
    ) |> dplyr::mutate(library_id = stringr::str_extract(fastq_1, "COHP_\\d{5}"))

    xls <- here::here("data-raw/FLT3 AML samples information (IGC-LZ-20342).xlsx")
    xls_ss <- here::here("data-raw/sample summary_IGC-LZ-20342.xlsx")

    sample_metadata <- readxl::read_excel(xls) |>
        janitor::clean_names("lower_camel") |>
        janitor::remove_constant(na.rm = TRUE) |>
        tidyr::fill(number, .direction = "down") |>
        tidyr::separate(htbIdPb, into = c(NA, NA, "sample_id"), sep = "-") |>
        dplyr::select(patient_id = number, sample_id, sample_date = date, cell_amount = cellAmt) |>
        tidyr::drop_na(sample_id) |>
        dplyr::arrange(patient_id, sample_date) |>
        dplyr::mutate(
            patient_id = gsub("#", "_", patient_id, fixed = TRUE)
        )

    sequencing_metadata <- readxl::read_excel(xls_ss) |>
        tidyr::separate(Sample_ID, into = c(NA, NA, "sample_id"), sep = "_") |>
        dplyr::mutate(
            cohort = cohort,
            library_id = TGen_Sample_Name,
            batch = "2022_C",
            strandedness = "reverse",
            sex = NA_character_,
            dob = NA,
            treatment = NA_character_,
            tissue = "PBMC"
        ) |>
        dplyr::select(sample_id, cohort, library_id, batch, strandedness, sex, dob, treatment, tissue) |>
        dplyr::left_join(fastqs)

    sample_sheet <- dplyr::left_join(sequencing_metadata, sample_metadata, by = "sample_id") |>
        dplyr::arrange(patient_id, sample_date) |>
        dplyr::select(library_id, fastq_1, fastq_2, strandedness, dplyr::everything())
    return(sample_sheet)
}
# MDS.rnaseq.EGAD00001003891
parse_metadata_MDS.rnaseq.EGAD00001003891 <- function() {
    cohort <- "MDS.rnaseq.EGAD00001003891"
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
    sex <- read.table(here::here("data-raw/EGAD00001003891/delimited_maps/Run_Sample_meta_info.map"), sep = "\t", header = FALSE) |>
        dplyr::mutate(
            patient_id = stringr::str_extract(V1, "(?<=subject_id\\=)(.*?)(?=[;])"),
            sex = stringr::str_extract(V1, "(?<=gender\\=)(.*?)(?=[;])"),
            sex = dplyr::recode(sex, male = "M", female = "F", unknown = NA_character_)
        ) |>
        dplyr::select(patient_id, sex) |>
        dplyr::distinct()
    files <- read.table(here::here("data-raw/EGAD00001003891/delimited_maps/Sample_File.map"),
        sep = "\t", header = FALSE,
        col.names = c("patient_id", "ega_sample_accession_id", "bam", "ega_file_unique_accession_id")
    ) |> dplyr::mutate(sample_name = gsub(".bam.cip", "", bam))

    experiment <- read.csv(here::here("data-raw/EGAD00001003891/delimited_maps/Study_Experiment_Run_sample.map"),
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
            cohort = cohort
        ) |>
        dplyr::select(library_id = ega_run_accession_id, fastq_1, fastq_2, strandedness, dplyr::everything(), cohort) |>
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
            cohort = "AML.PRJEB27973",
            library_id = sample,
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
        dplyr::select(library_id, fastq_1, fastq_2, strandedness, patient_id, tissue, timepoint, batch, treatment, genotype, sex, dob, cohort, everything())
    return(sample_sheet)
}
