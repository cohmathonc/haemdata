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

    # Create a vector of colnames to remove
    aml_metadata <- metadata_hsa |>
    # Select cols of interest
        dplyr::select("library_id", "fastq_1", "fastq_2", "strandedness", "patient_id", "sex", 
            "study_accession", "batch", "treatment", "tissue", "cohort", "sample_id",
            "sample_date", "age_at_diagnosis", diagnosis_year = "diagnosis year", "clinical_treatment", 
            date_of_diagnosis = "date of diagnosis", "timepoint", "genotype", "dob", "study_title",
            "sample_title", "sample_description", "sample_name") |>
    # remove "/net/nfs-irwrsrchnas01" from fastq paths
        dplyr::mutate(
            fastq_1 = stringr::str_remove(fastq_1, "^/net/nfs-irwrsrchnas01"),
            fastq_2 = stringr::str_remove(fastq_2, "^/net/nfs-irwrsrchnas01"),
    # tidy up redundant sample identifiers
            sample_id = dplyr::coalesce(sample_name, sample_id, sample_title),
    # fix FLT3 date_of_diagnosis dates
            date_of_diagnosis = as.character(date_of_diagnosis),
            date_of_diagnosis = suppressWarnings(
                dplyr::case_when(
                str_detect(date_of_diagnosis,"/", negate = TRUE) ~ as.character(as.Date(as.numeric(date_of_diagnosis), origin = "1899-12-30")),
                TRUE ~ as.character(lubridate::parse_date_time(date_of_diagnosis, c("%m/%Y", "%m/%d/%Y"))))
            ),
    # extract patient id where obvious
            patient_id = dplyr::case_when(
                study_accession == "PRJEB27973" ~ str_extract(sample_id, "[a-zA-Z0-9]{6}"),
                study_accession == "EGAS00001002346" ~ str_extract(sample_id, "PV\\d+|normal\\d"),
                TRUE ~ patient_id 
            ),
    # extract data from MDS samples
            timepoint = dplyr::case_when(
                study_accession == "EGAS00001002346" ~ str_extract(sample_id, "day\\d+"),
                TRUE ~ timepoint),
            genotype = dplyr::case_when(
                study_accession == "EGAS00001002346" ~ str_extract(sample_id, "Ctrl|Mut\\d|WT\\d"),
                TRUE ~ genotype),
            treatment = dplyr::case_when(
                study_accession == "EGAS00001002346" ~ str_extract(sample_id, "CHX|untreated|None"),
                TRUE ~ treatment),
            treatment = dplyr::case_when(
                study_accession == "EGAS00001002346" ~ ifelse(stringr::str_detect(treatment, "None"), "untreated", treatment),
                TRUE ~ treatment),
            treatment = dplyr::case_when(
                study_accession == "EGAS00001002346" ~ ifelse(stringr::str_detect(sample_id, "normal"), "Ctrl", treatment),
                TRUE ~ treatment),
            study_title = dplyr::case_when(
                study_accession == "EGAS00001002346" ~ "Aberrant splicing and defective mRNA production induced by somatic spliceosome mutations in myelodysplasia; Shiozawa et al Nat Comms 2018",
                TRUE ~ study_title),
            cohort = dplyr::case_when(
                study_accession == "EGAS00001002346" ~ "MDS.mRNA.EGAD00001003891",
                TRUE ~ cohort),
    # Fix cohort for Kim et al
            cohort = dplyr::case_when(
                study_accession == "PRJEB27973" ~ "AML.mRNA.PRJEB27973",
                TRUE ~ cohort),
    # relabel CD+ cells
            tissue = ifelse(stringr::str_detect(tissue, "CD34"), "BM CD34+", tissue)
        ) |>
        dplyr::select(-c(sample_name, sample_title, sample_description)) 

    # healthy PBMCs & BM fastqs
    PE_reads <- fastqs <- data.frame(
        fastq_1 = list.files(c(
            "/labs/rrockne/MHO/AML.human.data/data/GSE58335/fastq",
            "/labs/rrockne/MHO/AML.human.data/data/fastq/PE-reads"),
                pattern = "_1\\.fastq\\.gz$", full.names = TRUE, recursive = TRUE),
        fastq_2 = list.files(c(
            "/labs/rrockne/MHO/AML.human.data/data/GSE58335/fastq",
            "/labs/rrockne/MHO/AML.human.data/data/fastq/PE-reads"),
            pattern = "_2\\.fastq\\.gz$", full.names = TRUE, recursive = TRUE
        )
        ) |> dplyr::mutate(library_id = stringr::str_extract(fastq_1, "SRR\\d{7}"),
        strandedness = case_when(
            grepl("marrow", fastq_1) ~ "unstranded",
            TRUE ~ "reverse"
        ))

    SE_reads <- fastqs <- data.frame(
        fastq_1 = list.files(c(
            "/labs/rrockne/MHO/AML.human.data/data/fastq/SE-reads"),
                pattern = "fastq\\.gz$", full.names = TRUE, recursive = TRUE),
        fastq_2 = ""
        ) |> dplyr::mutate(library_id = stringr::str_extract(fastq_1, "SRR\\d{7}"),
        strandedness = case_when(
            grepl("PBMC", fastq_1) ~ "unstranded",
            grepl("BM", fastq_1) ~ "forward",
            TRUE ~ "unstranded"
        ))
    
    reads <- rbind(PE_reads, SE_reads)

    # metadata for healthy PBMCs & BM
    non_aml_metadata <- readr::read_tsv("/labs/rrockne/MHO/AML.human.data/combined.tsv", show_col_types = FALSE) |>
    # drop cols
        dplyr::select(-matches(c("sample_accession", "single_end", "scientific", "experiment", "description", "fastq", "md5", "tax_id", "species", "id", "submission", "secondary", "alias", "library", "instrument", "count")))

    # merge metadata and reads
    non_aml_sample_sheet <- dplyr::left_join(
        non_aml_metadata,
        reads,
        by = c("run_accession" = "library_id")
    ) |>
        dplyr::group_by(study_accession) %>%
        dplyr::mutate(
            batch = paste0("non_aml.mRNA_", LETTERS[as.integer(cur_group_id())]),
            tissue = dplyr::case_when(
                stringr::str_detect(sample_title, "PBMC|Control|GBS") ~ "PBMC",
                stringr::str_detect(sample_title, "CD34") ~ "BM CD34+",
                TRUE ~ "BM"
            ),
            patient_id = dplyr::case_when(
                study_accession == "PRJNA294808" ~ if_else(stringr::str_detect(sample_title, "GBS"), "GBS_twin", "Ctrl_twin"),
                study_accession == "PRJNA487456" ~ str_remove(sample_title, "CD34\\+ BM donor "),
                study_accession == "PRJNA252189" ~ str_remove(sample_title, "PBMC\\.none\\.0h\\.X.\\.0"),
                study_accession == "PRJNA493081" ~ sample_title,
                TRUE ~ str_extract(sample_title, "\\d+") 
            ),
            sex = dplyr::case_when(
                stringr::str_detect(sample_title, "XX")  ~ "F",
                stringr::str_detect(sample_title, "XY")  ~ "M",
                TRUE ~ "" 
            ),
            diagnosis_year = "",
            clinical_treatment = "",
            date_of_diagnosis = "",
            age_at_diagnosis = "",
            sample_date = "",
            timepoint = "",
            genotype = "",
            treatment = "Ctrl",
            dob = "",
            cohort = "AML.mRNA.non_aml") |> ungroup() |>
            dplyr::select(library_id = run_accession, fastq_1, fastq_2, strandedness, sample_id = "sample_title", dplyr::everything())

        # put the two together
        combined_df <- rbind(non_aml_sample_sheet, aml_metadata) |>
            # replace all "" with NA_character_
            dplyr::select(-matches("fastq")) |>
            dplyr::mutate_all(list(~ dplyr::na_if(., "")))
            
    return(combined_df)
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

    # v0.0.0.9012 interim
    sample_sheet <- pins::pin_read(
        pins::board_folder(
            "/labs/rrockne/MHO/haemdata",
            versioned = FALSE
        ),
        "metadata_mmu.csv"
    ) |>
        purrr::modify_if(is.factor, as.character) |>
        dplyr::filter(cohort != "CML.blastcrisis.2023")

    # # # fix an issue with some CML mice have incorrect dod ----
    dod_fix <- tibble(
        mouse_id = as.character(c("487", "488", "489", "490", "541", "545")),
        dod = rep(NA, 6)
    )
    sample_sheet <- rows_update(sample_sheet, dod_fix, by = "mouse_id")

    # # # add survival column ----
    sample_sheet <- sample_sheet |>
        add_survival_columns() |>
        select(-(matches("age_at_end|age_at_sample")))

    # # # relocate some columns ----
    sample_sheet <- sample_sheet |>
        relocate(dead, .after = age_at_start)|>
        relocate(qc_pass_mapping, .after = everything())

    # # # add "blast crisis" mouse metadata
    blast_crisis_sample_sheet <- readxl::read_excel(here::here("data-raw/metadata_mmu_blastcrisis.xlsx")) |>
        mutate(
            dod = if_else(dod == "no die", NA_character_, dod),
            dod = as.Date(as.numeric(dod), origin = "1899-12-30"),
            mouse_id = as.character(mouse_id)
        ) |>
    # sequencing info
    left_join(readxl::read_excel(here::here("data-raw/sequencing summary_IGC-BZ-21029.xlsx")) |>
        select(library_id = "TGen_Sample_Name", prep_date = "sample  prep date", seq_date = "sequencing date", seq_length = "sequencing length...17") |>
        na.omit(), by = "library_id") |>
    # fastqs
    left_join(read.csv("data-raw/mmu_blastcrisis_10Xfastqs.csv") |>
        dplyr::filter(!grepl("counts|_WT_|_KO_", fastq)) |> ### Filter out analyses files from the IGC
        dplyr::mutate(
            library_id = stringr::str_replace(fastq, ".*(COHP_\\d*)_.*", "\\1"),
            fastq_1 = case_when(
                grepl("R1_001", fastq) ~ fastq, TRUE ~ NA_character_
            ),
            fastq_2 = case_when(
                grepl("R2_001", fastq) ~ fastq, TRUE ~ NA_character_
            ),
            basename = gsub("R[12]_001.fastq.gz$", "", fastq)
        ) |>
        group_by(basename) |>
        tidyr::fill(fastq_1, fastq_2, .direction = "downup") |>
        ungroup() |>
        select(-fastq, -basename) |>
        distinct(), by = "library_id") |>
        mutate(batch = paste(
            "2023",
            factor(prep_date, labels = paste0("L", seq_along(unique(prep_date)))),
            factor(seq_date, labels = paste0("S", seq_along(unique(seq_date)))),
            factor(seq_length, labels = paste0("R", seq_along(unique(seq_length)))),
            sep = "_")) |>
        add_survival_columns() |>
            select(-c(label, prep_date, seq_date, seq_length))

    # # # Add to sample_sheet, assign new sample_id
        sample_sheet <- rows_insert(sample_sheet, blast_crisis_sample_sheet, by = "library_id") |>
            arrange(sample_id) |>
            assign_new_sample_ids()

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

    fastqs <- data.frame(
        fastq_1 = list.files(c("/labs/ykuo/Seq/220511_IGC-LZ-20342", "/labs/ykuo/Seq/220715_IGC-LZ-20821"),
                pattern = "_R1_.*\\.gz$", full.names = TRUE, recursive = TRUE
            ),
        fastq_2 = list.files(c("/labs/ykuo/Seq/220511_IGC-LZ-20342", "/labs/ykuo/Seq/220715_IGC-LZ-20821"),
            pattern = "_R2_.*\\.gz$", full.names = TRUE, recursive = TRUE
        )
    ) |> dplyr::mutate(library_id = stringr::str_extract(fastq_1, "COHP_\\d{5}"))

    xls <- here::here("data-raw/flt3_sample_metadata.xlsx")

    sample_sheet <- readxl::read_excel(xls) |>
        dplyr::arrange(patient_id, sample_date) |>
        dplyr::mutate(
            sample_date = as.character(sample_date),
            strandedness = "reverse") |>
        dplyr::left_join(fastqs, by = "library_id") |>
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
