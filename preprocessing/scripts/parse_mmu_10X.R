#### 10X -----

# ├ cml_mir142_ko
parse_10x_mir142_ko <- function() {

    cml_mir142_ko_fastq_paths <- c(
        "/labs/gmarcucci/Seq/221122_IGC-BZ-21028_Test_RUN",
        "/labs/gmarcucci/Seq/221216_IGC-BZ-21028_full_run",
        "/labs/gmarcucci/Seq/221019_IGC-BZ-20945_full_run PBMC and BM-T from miR-142 KO vs WT",
        "/labs/gmarcucci/Seq/221010_IGC-BZ-20945_test_run PBMC and BM-T from miR-142 KO vs WT",
        "/labs/gmarcucci/Seq/221007_IGC-BZ-20895 Spleen-T from miR-142 KO vs WT",
        "/labs/gmarcucci/Seq/220610_IGC-BZ-20758 PBMC and BM-T scRNA-seq",
        "/labs/gmarcucci/Seq/220610_IGC-BZ-20758-test_run KO vs WT CML BM-T scRNA-seq test run",
        "/labs/gmarcucci/Seq/220106_BinZhang Spleen-T scRNA-seq",
        "/labs/gmarcucci/Seq/211210_BinZhang Spleen-T scRNA-seq test run",
        "/labs/gmarcucci/Seq/221109_IGC-BZ-20933_full_run Spleen LMPP scRNA-seq"
    )

    if (file.exists("data-raw/mmu_mir142_ko_10Xfastqs.csv")) {
        fastq_paths <- read.csv("data-raw/mmu_mir142_ko_10Xfastqs.csv")
    } else {
        fastq_paths <- data.frame(fastq = list.files(
            path = cml_mir142_ko_fastq_paths,
            pattern = "R[12]_001.fastq.gz$", full.names = TRUE, recursive = TRUE
        ))
        rio::export(fastq_paths, "data-raw/mmu_mir142_ko_10Xfastqs.csv")
    }

    sample_sheet <- fastq_paths |>
        dplyr::filter(!grepl("counts|_WT_|_KO_|COHP_\\d{5}/Project", fastq)) |>    ### Filter out analyses files from the IGC
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
            distinct()


    return(sample_sheet)
}
