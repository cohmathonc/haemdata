z# Functions for importing metadata

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
#' detailed in [R/import_metadata.R](https://github.com/drejom/haemdata/blob/main/R/import_metadata.R).
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
    # TODO fix weeks for BM samples

    return(metadata)
}

# ├ AML.mRNA.2016 
parse_metadata_AML.mRNA.2016 <- function() {
    project <- "AML.mRNA.2016"
    sample_sheet <- read.csv("/net/isi-dcnl/ifs/user_data/rrockne/MHO/AML.mRNA.2016/config/new_name_key.tsv", , sep = "\t") |>
        dplyr::mutate(
            project = project,
            fastq_1 = paste0("/net/isi-dcnl/ifs/user_data/rrockne/MHO/AML.mRNA.2016/data/", newname, ".gz"),
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
    sample_sheet <- read.csv("/net/isi-dcnl/ifs/user_data/rrockne/MHO/AML.mRNA.2018.all_samples/config/new_name_key.tsv", sep = "\t") |>
        dplyr::mutate(
            project = project,
            fastq_1 = paste0("/net/isi-dcnl/ifs/user_data/rrockne/MHO/AML.mRNA.2018.all_samples/data/fastq/", Newname, ".gz"),
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
    fastqs <- data.frame(
        fastq_1 = system(
            paste0("find '/net/isi-dcnl/ifs/user_data/ykuo/Seq/201124_TGen' -name '*.gz'"),
            intern = TRUE
        ) %>% grep("_R1_", ., value = TRUE),
        fastq_2 = system(
            paste0("find '/net/isi-dcnl/ifs/user_data/ykuo/Seq/201124_TGen' -name '*.gz'"),
            intern = TRUE
        ) %>% grep("_R2_", ., value = TRUE)
    ) |> dplyr::mutate(library_id = stringr::str_extract(fastq_1, "COHP_\\d{5}"))
    sample_sheet <- read.csv("/net/isi-dcnl/ifs/user_data/rrockne/MHO/AML.mRNA.2020/config/rename_fastqs.tsv", sep = "\t", header = FALSE) |>
        tidyr::separate(V2, sep = "_", c("treatment", "mouse_id", "timepoint")) |>
        dplyr::mutate(
            project = project,
            library_id = V1,
            tissue = "PBMC",
            strandedness = "reverse",
            dob = NA,
            sex = NA_character_,
            genotype = NA_character_,
            batch = "2020_A",
            tissue = dplyr::case_when(
                treatment == "BM" ~ "BM",
                TRUE ~ tissue
            ),
            treatment = dplyr::case_when(
                treatment == "BM" ~ timepoint,
                TRUE ~ treatment
            ),
            timepoint = dplyr::case_when(
                tissue == "BM" ~ NA_character_,
                TRUE ~ timepoint
            )
        ) |>
        dplyr::left_join(fastqs) |>
        dplyr::select(sample = library_id, fastq_1, fastq_2, strandedness, mouse_id, tissue, timepoint, batch, treatment, genotype, sex, dob, project)
    return(sample_sheet)
}

# AML.mRNA.2021.RxGroup1 
parse_metadata_AML.mRNA.2021.RxGroup1 <- function() {
    project <- "AML.mRNA.2021.RxGroup1"
    fastqs <- data.frame(
        fastq_1 = system(
            paste0("find '/net/isi-dcnl/ifs/user_data/rrockne/Seq/210412_Tgen' -name '*.gz'"),
            intern = TRUE
        ) %>% grep("_R1_", ., value = TRUE),
        fastq_2 = system(
            paste0("find '/net/isi-dcnl/ifs/user_data/rrockne/Seq/210412_Tgen' -name '*.gz'"),
            intern = TRUE
        ) %>% grep("_R2_", ., value = TRUE)
    ) |> dplyr::mutate(library_id = stringr::str_extract(fastq_1, "COHP_\\d{5}"))
    xls <- "/net/isi-dcnl/ifs/user_data/rrockne/MHO/AML.mRNA.2021.RxGroup1/config/run_summary_IGC-LZ-18757.xlsx"
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
            timepoint = sub("BM", NA_character_, timepoint)
        ) |>
        dplyr::left_join(fastqs)
    xls_cm <- "/home/domeally/workspaces/haemdata/data-raw/CM mice chemo Rx survival.xlsx"
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
            paste0("find '/net/isi-dcnl/ifs/user_data/rrockne/Seq/210507_Tgen' -name '*.gz'"),
            intern = TRUE
        ) %>% grep("_R1_", ., value = TRUE),
        fastq_2 = system(
            paste0("find '/net/isi-dcnl/ifs/user_data/rrockne/Seq/210507_Tgen' -name '*.gz'"),
            intern = TRUE
        ) %>% grep("_R2_", ., value = TRUE)
    ) |> dplyr::mutate(library_id = stringr::str_extract(fastq_1, "COHP_\\d{5}"))
    xls <- "/net/isi-dcnl/ifs/user_data/rrockne/MHO/AML.mRNA.2021.RxGroup2/config/sequencing summary_IGC-LZ-18862.xlsx"
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
    xls_cm <- "/home/domeally/workspaces/haemdata/data-raw/CM mice chemo Rx survival.xlsx"
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
            paste0("find '/net/isi-dcnl/ifs/user_data/ykuo/Seq/210910_TGen' -name '*.gz'"),
            intern = TRUE
        ) %>% grep("_R1_", ., value = TRUE),
        fastq_2 = system(
            paste0("find '/net/isi-dcnl/ifs/user_data/ykuo/Seq/210910_TGen' -name '*.gz'"),
            intern = TRUE
        ) %>% grep("_R2_", ., value = TRUE)
    ) |> dplyr::mutate(library_id = stringr::str_extract(fastq_1, "COHP_\\d{5}"))
    sample_sheet <- read.csv("/net/isi-dcnl/ifs/user_data/rrockne/MHO/AML.mRNA.2021.RxGroup2_pt2/config/rename_fastqs.tsv", sep = "\t", header = FALSE) |>
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
    xls_cm <- "/home/domeally/workspaces/haemdata/data-raw/CM mice chemo Rx survival.xlsx"
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
            paste0("find '/net/isi-dcnl/ifs/user_data/ykuo/Seq/220415_IGC-LZ-20205' -name '*.gz'"),
            intern = TRUE
        ) %>% grep("_R1_", ., value = TRUE),
        fastq_2 = system(
            paste0("find '/net/isi-dcnl/ifs/user_data/ykuo/Seq/220415_IGC-LZ-20205' -name '*.gz'"),
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
    xls_cm <- "/home/domeally/workspaces/haemdata/data-raw/CM mice chemo Rx survival.xlsx"
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
    #         paste0("find '/net/isi-dcnl/ifs/user_data/ykuo/Seq/200829' -name '*.gz'"),
    #         intern = TRUE
    #     ) %>% grep("_R1_", ., value = TRUE),
    #     fastq_2 = system(
    #         paste0("find '/net/isi-dcnl/ifs/user_data/ykuo/Seq/200829' -name '*.gz'"),
    #         intern = TRUE
    #     ) %>% grep("_R2_", ., value = TRUE)
    # ) |> dplyr::mutate(library_id = paste0("IGCP_", substr(basename(fastq_1), 1, 5)))
    ## fastqs generated at TGen
    fastqs <- data.frame(
        fastq_1 = system(
            paste0("find '/net/isi-dcnl/ifs/user_data/ykuo/Seq/200928_TGen' -name '*.gz'"),
            intern = TRUE
        ) %>% grep("_R1_", ., value = TRUE),
        fastq_2 = system(
            paste0("find '/net/isi-dcnl/ifs/user_data/ykuo/Seq/200928_TGen' -name '*.gz'"),
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
            paste0("find '/net/isi-dcnl/ifs/user_data/rrockne/MHO/AML.validation.2017/data/fastq' -name '*.gz'"),
            intern = TRUE
        ) %>% grep("_R1_", ., value = TRUE)
    ) |> dplyr::mutate(library_id = paste0("COHP_", substr(basename(fastq_1), 1, 5)))

    xls <- "/net/isi-dcnl/ifs/user_data/rrockne/MHO/AML.validation.2017/data/validation_name_key.xlsx"
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
            paste0("find '/net/isi-dcnl/ifs/user_data/rrockne/Seq/211210_IGC-LZ-19773' -name '*.gz'"),
            intern = TRUE
        ) %>% grep("_R1_", ., value = TRUE),
        fastq_2 = system(
            paste0("find '/net/isi-dcnl/ifs/user_data/rrockne/Seq/211210_IGC-LZ-19773' -name '*.gz'"),
            intern = TRUE
        ) %>% grep("_R2_", ., value = TRUE)
    ) |> dplyr::mutate(library_id = stringr::str_extract(fastq_1, "COHP_\\d{5}"))

    xls <- "/net/isi-dcnl/ifs/user_data/rrockne/MHO/AML.scRNA.2022/bulkseq_targets/config/sequencing summary_IGC-LZ-19773.xlsx"
    sample_sheet <- readxl::read_excel(xls) |>
        dplyr::mutate(
            project = project,
            library_id = paste0("COHP_", gsub("_.*", "", Sample_Name)),
            mouse_id = gsub(".*_", "", Sample_Name) |> substr(1, 4),
            timepoint = NA_character_,
            batch = "2022_B",
            strandedness = "reverse",
            sex = NA_character_,
            dob = NA,
            treatment = NA_character_,
            genotype = NA_character_,
            tissue = gsub(".*_\\d{4}", "", Sample_Name),
            tissue = dplyr::case_when(
                tissue == "BM" ~ "BM",
                stringr::str_detect(Sample_Name, "CKIT") ~ "BM_CKIT",
                TRUE ~ "PBMC"
            )
        ) |>
        dplyr::left_join(fastqs) |>
        dplyr::select(sample = library_id, fastq_1, fastq_2, strandedness, mouse_id, tissue, timepoint, batch, treatment, genotype, sex, dob, project)
    return(sample_sheet)
}
# CML.mRNA.2021
parse_metadata_CML.mRNA.2021 <- function() {
    cohort <- "CML.mRNA.2021"
    fastqs <- data.frame(
        fastq_1 = system(
            paste0("find /net/isi-dcnl/ifs/user_data/rrockne/Seq/210910_TGen /net/isi-dcnl/ifs/user_data/rrockne/Seq/210907_TGen -name '*.gz'"),
            intern = TRUE
        ) %>% grep("_R1_", ., value = TRUE),
        fastq_2 = system(
            paste0("find /net/isi-dcnl/ifs/user_data/rrockne/Seq/210910_TGen /net/isi-dcnl/ifs/user_data/rrockne/Seq/210907_TGen -name '*.gz'"),
            intern = TRUE
        ) %>% grep("_R2_", ., value = TRUE)
    ) |> dplyr::mutate(library_id = stringr::str_extract(fastq_1, "COHP_\\d{5}"))

    xls <- "/net/isi-dcnl/ifs/user_data/rrockne/MHO/CML.mRNA.2021/config/Lianjun_IGC-LZ-19411.xlsx"
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
            paste0("find /net/isi-dcnl/ifs/user_data/gmarcucci/Seq/220210_IGC-LZ-19931 /net/isi-dcnl/ifs/user_data/gmarcucci/Seq/220202_IGC-LZ-19931 -name '*.gz'"),
            intern = TRUE
        ) %>% grep("_R1_", ., value = TRUE),
        fastq_2 = system(
            paste0("find /net/isi-dcnl/ifs/user_data/gmarcucci/Seq/220210_IGC-LZ-19931 /net/isi-dcnl/ifs/user_data/gmarcucci/Seq/220202_IGC-LZ-19931 -name '*.gz'"),
            intern = TRUE
        ) %>% grep("_R2_", ., value = TRUE)
    ) |> dplyr::mutate(library_id = stringr::str_extract(fastq_1, "COHP_\\d{5}"))

    sample_sheet <- read.csv("/net/isi-dcnl/ifs/user_data/rrockne/MHO/CML.mRNA.2022/config/sample_sheet.csv") |>
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
# AML.mRNA.HSA_FLT3.2022 
parse_metadata_AML.mRNA.HSA_FLT3.2022 <- function() {
    project <- "AML.mRNA.HSA_FLT3.2022"
    fastqs <- data.frame(
        fastq_1 = system(
            paste0("find /net/isi-dcnl/ifs/user_data/ykuo/Seq/220511_IGC-LZ-20342 -name '*.gz'"),
            intern = TRUE
        ) %>% grep("_R1_", ., value = TRUE),
        fastq_2 = system(
            paste0("find /net/isi-dcnl/ifs/user_data/ykuo/Seq/220511_IGC-LZ-20342 -name '*.gz'"),
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
        fastq_1 = list.files("/net/isi-dcnl/ifs/user_data/rrockne/MHO/MDS.rnaseq.EGAD00001003891/data/fastq2/reads",
            pattern = "1.fq.gz", full.names = TRUE
        ),
        fastq_2 = list.files("/net/isi-dcnl/ifs/user_data/rrockne/MHO/MDS.rnaseq.EGAD00001003891/data/fastq2/reads",
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
    by = "ega_run_accession_id")

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
                stringr::str_detect(sample_name, "CD34") ~ "CD34"),
            project = project
        ) |>
        dplyr::select(sample = ega_run_accession_id, fastq_1, fastq_2, strandedness, dplyr::everything(), project) |>
        dplyr::filter(library_strategy == "RNA-Seq")
    return(sample_sheet)
}

# AML.PRJEB27973 Kim et al 2020 SciRep
parse_metadata_AML.PRJEB27973 <- function() {
    sample_sheet <- readr::read_csv(
        "/net/isi-dcnl/ifs/user_data/rrockne/MHO/AML.PRJEB27973/samplesheet/samplesheet.csv",
        show_col_types = FALSE) |>
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
