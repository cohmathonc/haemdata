#### 10X -----

# ├ cml_mir142_ko
parse_10x_mir142_ko <- function() {
    cml_mir142_ko_fastq_paths <- c(
        # Full runs
        "/labs/gmarcucci/Seq/221216_IGC-BZ-21028_full_run BM LMPP scRNA-seq",
        "/labs/gmarcucci/Seq/221109_IGC-BZ-20933_full_run Spleen LMPP scRNA-seq",
        "/labs/gmarcucci/Seq/221019_IGC-BZ-20945_full_run_PBMC_and_BM-T_from_miR-142_KO_vs_WT",
        "/labs/gmarcucci/Seq/220919_IGC-BZ-20895",
        "/labs/gmarcucci/Seq/220610_IGC-BZ-20758 PBMC and BM-T scRNA-seq",
        "/labs/gmarcucci/Seq/220106_BinZhang Spleen-T scRNA-seq",
        # Test runs
        "/labs/gmarcucci/Seq/221122_IGC-BZ-21028_Test_RUN",
        "/labs/gmarcucci/Seq/221010_IGC-BZ-20933_test_run Spleen LMPP scRNA-seq",
        "/labs/gmarcucci/Seq/221010_IGC-BZ-20945_test_run PBMC and BM-T from miR-142 KO vs WT",
        "/labs/gmarcucci/Seq/220520_IGC-BZ-20758 KO vs WT CML PBMC scRNA-seq test run",
        "/labs/gmarcucci/Seq/220610_IGC-BZ-20758-test_run KO vs WT CML BM-T scRNA-seq test run",
        "/labs/gmarcucci/Seq/211210_BinZhang Spleen-T scRNA-seq test run")

    if (file.exists(here::here(work_dir, "data-raw/mmu_mir142_ko_10Xfastqs.csv"))) {
        fastq_paths <- read.csv(here::here(work_dir, "data-raw/mmu_mir142_ko_10Xfastqs.csv"))
    } else {
        ncpus <- as.integer(Sys.getenv("SLURM_CPUS_ON_NODE", 1))

        tmp_paths <- tempfile()
        writeLines(cml_mir142_ko_fastq_paths, tmp_paths)

        cmd <- sprintf(
            "cat %s | tr '\\n' '\\0' | xargs -0 -P %d -I{} find {} -type f -name '*R[12]_001.fastq.gz'",
            tmp_paths, ncpus
        )

        fastq_paths <- data.frame(fastq = system(cmd, intern = TRUE))
        unlink(tmp_paths)

        rio::export(fastq_paths, here::here(work_dir, "data-raw/mmu_mir142_ko_10Xfastqs.csv"))
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

# ├ cml_blastcrisis
parse_10x_blastcrisis <- function() {
    cml_blastcrisis_fastq_paths <- c(
        "/labs/gmarcucci/Seq/230210_IGC-BZ-21029_FullRun",
        "/labs/gmarcucci/Seq/230119_IGC-BZ-21029_test_run"
    )

    if (file.exists(here::here(work_dir, "data-raw/mmu_blastcrisis_10Xfastqs.csv"))) {
        fastq_paths <- read.csv(here::here(work_dir, "data-raw/mmu_blastcrisis_10Xfastqs.csv"))
    } else {
        ncpus <- as.integer(Sys.getenv("SLURM_CPUS_ON_NODE", 1))

        tmp_paths <- tempfile()
        writeLines(cml_blastcrisis_fastq_paths, tmp_paths)

        cmd <- sprintf(
            "cat %s | tr '\\n' '\\0' | xargs -0 -P %d -I{} find {} -type f -name '*R[12]_001.fastq.gz'",
            tmp_paths, ncpus
        )

        fastq_paths <- data.frame(fastq = system(cmd, intern = TRUE))
        unlink(tmp_paths)

        rio::export(fastq_paths, here::here(work_dir, "data-raw/mmu_blastcrisis_10Xfastqs.csv"))
    }


    sample_sheet <- fastq_paths |>
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
        distinct()


    return(sample_sheet)
}
