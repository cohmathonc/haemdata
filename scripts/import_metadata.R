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

#' Update metadata for mouse AML samples
#'
#' This function prepares the `metadata_mmu` object for all RNAseq libraries from AML and CML mice.
#' Minimal metadata fields include sample, fastq_1, fastq_2, strandedness,
#' mouse_id, tissue, sample_date, week, timepoint, batch, treatment, genotype,
#' sex, dob, dod, project.
#'
#' Raise an \href{https://github.com/drejom/haemdata/issues}{issue on GitHub}
#' to report erroneous or missing records.
#' @details Updates the `metadata_mmu` object for all RNAseq libraries by using
#' a published pin, previously assembled by the make_metadata_mmu() function,
#' which was retired in version 0.0.0.0.9008.
#' @name update_metadata_mmu
#' @param sample_sheet_all_mice a data.frame produced by row binding all mouse sample sheets
#' detailed in [R/import_metadata.R](https://github.com/drejom/haemdata/blob/main/R/import_metadata.R).
#' @return a data.frame
#' @author Denis O'Meally

update_metadata_mmu <- function() {
    # Update sample and mouse metadata
    use_pinboard("onedrive")
    all_mice <- get_pin("metadata_mmu.csv", "20221003T041923Z-a93aa") |>
        purrr::modify_if(is.factor, as.character)

    # Fix bulk scRNAseq samples
    # download the file
    get_teams_file("General/sequencing summary_IGC-LZ-19773.xlsx")

    corrected_tissue <- readxl::read_excel(
        "data-raw/sequencing summary_IGC-LZ-19773.xlsx",
        col_names = TRUE,
        col_types = "text"
    ) |>
        dplyr::select(Sample_Name) |>
            dplyr::mutate(
                sample = paste0("COHP_",
                stringr::str_extract(Sample_Name, "\\d{5}")), tissue = "PBMC_CKIT") |>
            dplyr::filter(stringr::str_detect(Sample_Name, "PBCKIT")) |>
            dplyr::select(-Sample_Name)

    # Fix DOD for AML samples from @yufu1120
    # download the file
    get_teams_file("General/metadata_mmu_correct_YHF_20221028.xlsx")
    corrected_dod <- readxl::read_excel(
        "data-raw/metadata_mmu_correct_YHF_20221028.xlsx",
        col_names = TRUE,
        col_types = "text"
    ) |>
        dplyr::select(sample, dod) |>
        dplyr::mutate(dod = as.Date(as.numeric(dod), origin = "1899-12-30") |> as.character()) |>
        dplyr::distinct()

    # Update the sample metadata
    sample_sheet <- dplyr::rows_update(all_mice, corrected_tissue, by = "sample", unmatched = "ignore") |>
                    dplyr::rows_update(corrected_dod, by = "sample", unmatched = "ignore")

    # Add in cellranger h5 files (for Gencode reference)
    h5_paths <- scan("data-raw/cellranger_h5_paths.txt", character()) |>
        as_tibble() |>
        dplyr::filter(str_detect(value, "GENCODE")) |>    dplyr::mutate(
            mouse_id = stringr::str_extract(value, "(?<=per_sample_outs/)(.*?)(?=[_])") |> as.numeric()    ,
            tissue = stringr::str_extract(value, "(?<=per_sample_outs/\\d{4}_)(.*?)(?=[/])"),
            tissue = stringr::str_replace(tissue, "PB", "PBMC"),
            tissue = stringr::str_replace(tissue, "ckit", "CKIT"),
            hdf5 = value
        ) |>
        dplyr::select(-value)

    sample_sheet <- dplyr::left_join(sample_sheet, h5_paths, by = c("mouse_id", "tissue")) |>
        dplyr::reloacte(hdf5, .after = fastq_2)

    return(sample_sheet)
}
#### mRNA -----
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
