#### miRNA -----

# ├ AML.mRNA.2016
parse_metadata_AML.miRNA.2016 <- function() {
    project <- "AML.miRNA.2016"
    assay <- "miRNA"
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
        dplyr::select(library_id, fastq_1, mouse_id, timepoint, tissue, batch, project, assay)
    return(sample_sheet)
}
# ├ AML.miRNA.2018
parse_metadata_AML.miRNA.2018 <- function() {
    project <- "AML.miRNA.2018"
    assay <- "miRNA"
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
        dplyr::select(library_id, fastq_1 = value, mouse_id, timepoint, tissue, batch, project, assay)
    return(sample_sheet)
}
# ├ AML.miRNA.2020
# * NOTE: Unique to this miRNA project, the samples
# * were paired end, but only fastq_1 was used.
parse_metadata_AML.miRNA.2020 <- function() {
    project <- "AML.miRNA.2020"
    assay <- "miRNA"
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
        dplyr::select(library_id, fastq_1, mouse_id, timepoint, tissue, batch, project, assay)
    return(sample_sheet)
}

# ├ AML.miRNA.2021.RxGroup1
parse_metadata_AML.miRNA.2021.RxGroup1 <- function() {
    project <- "AML.miRNA.2021.RxGroup1"
    assay <- "miRNA"
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
        dplyr::select(library_id, fastq_1, mouse_id, timepoint, tissue, batch, project, assay)
    return(sample_sheet)
}

# ├ AML.miRNA.2021.RxGroups1and2
parse_metadata_AML.miRNA.2021.RxGroups1and2 <- function() {
    project <- "AML.miRNA.2021.RxGroups1and2"
    assay <- "miRNA"
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
        dplyr::select(library_id, fastq_1, mouse_id, timepoint, tissue, batch, project, assay)
    return(sample_sheet)
}

# ├ AML.miRNA.2021.RxGroup2_pt2
parse_metadata_AML.miRNA.2021.RxGroup2_pt2 <- function() {
    project <- "AML.miRNA.2021.RxGroup2_pt2"
    assay <- "miRNA"
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
        dplyr::select(library_id, fastq_1, mouse_id, timepoint, tissue, batch, project, assay)
    return(sample_sheet)
}


# ├ AML.miRNA.2021.RxGroup3
parse_metadata_AML.miRNA.2022.RxGroup3 <- function() {
    project <- "AML.miRNA.2022.RxGroup3"
    assay <- "miRNA"
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
        dplyr::select(library_id, fastq_1 = value, mouse_id, timepoint, tissue, batch, project, assay)
    return(sample_sheet)
}
