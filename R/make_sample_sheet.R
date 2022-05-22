#' .. content for \description{} (no empty lines) ..
#'
#' .. content for \details{} ..
#'
#' @title Make an nf-core/rnaseq sample sheet
#' @param project the project ID for which a `sample_sheet` will be created
#' @return
#' @author Denis O'Meally
#' @export
make_sample_sheet <- function(project) {
    # sample = library_id, fastq1, fastq2, strandedness, mouse_id, tissue, timepoint, batch, treatment, genotype, sex, dob
    # AML.mRNA.2016 ├-----------------------------------------------------------
    if (project == "AML.mRNA.2016") {
        sample_sheet <- read.csv("/net/isi-dcnl/ifs/user_data/rrockne/MHO/AML.mRNA.2016/config/new_name_key.tsv", , sep = "\t") |>
            dplyr::mutate(
                fastq1 = paste0("/net/isi-dcnl/ifs/user_data/rrockne/MHO/AML.mRNA.2016/data/", newname, ".gz"),
                library_id = paste0("COHP_", ID),
                tissue = "PBMC",
                fastq2 = NA_character_,
                dob = NA_character_,
                Batch = paste0("2016_", Batch),
                strandedness = "reverse"
            ) |>
            dplyr::select(sample = library_id, fastq1, fastq2, strandedness, mouse_id = Mouse, tissue, timepoint = TimePoint, batch = Batch, treatment = Treatment, genotype = Genotype, sex = Sex, dob)
    } else if (project == "AML.mRNA.2018.all_samples") {
        # AML.mRNA.2018.all_samples ├-----------------------------------------------------------
        sample_sheet <- read.csv("/net/isi-dcnl/ifs/user_data/rrockne/MHO/AML.mRNA.2018.all_samples/config/new_name_key.tsv", sep = "\t") |>
            dplyr::mutate(
                fastq1 = paste0("/net/isi-dcnl/ifs/user_data/rrockne/MHO/AML.mRNA.2018.all_samples/data/fastq/", Newname, ".gz"),
                library_id = paste0("COHP_", gsub(".*_|\\.fq", "", Newname)),
                tissue = "PBMC",
                fastq2 = NA_character_,
                dob = NA_character_
            ) |>
            tidyr::separate(Newname, sep = "_", into = c("timepoint", "mouse_id", NA, "treatment", "genotype", "sex", "batch", NA)) |>
            dplyr::mutate(strandedness = dplyr::case_when(
                stringr::str_detect(timepoint, "L|T0|T1") ~ "unstranded",
                TRUE ~ "reverse"
            )) |>
            dplyr::select(sample = library_id, fastq1, fastq2, strandedness, mouse_id, tissue, timepoint, batch, treatment, genotype, sex, dob) |>
            dplyr::mutate(batch = paste0("2018_", batch))
    } else if (project == "AML.mRNA.2020") {
        # AML.mRNA.2020 ├-----------------------------------------------------------
        fastqs <- data.frame(
            fastq1 = system(
                paste0("find '/net/isi-dcnl/ifs/user_data/ykuo/Seq/201124_TGen' -name '*.gz'"),
                intern = TRUE
            ) %>% grep("_R1_", ., value = TRUE),
            fastq2 = system(
                paste0("find '/net/isi-dcnl/ifs/user_data/ykuo/Seq/201124_TGen' -name '*.gz'"),
                intern = TRUE
            ) %>% grep("_R2_", ., value = TRUE)
        ) |> dplyr::mutate(library_id = stringr::str_extract(fastq1, "COHP_\\d{5}"))
        sample_sheet <- read.csv("/net/isi-dcnl/ifs/user_data/rrockne/MHO/AML.mRNA.2020/config/rename_fastqs.tsv", sep = "\t", header = FALSE) |>
            tidyr::separate(V2, sep = "_", c("treatment", "mouse_id", "timepoint")) |>
            dplyr::mutate(
                library_id = V1,
                tissue = "PBMC",
                strandedness = "reverse",
                dob = NA_character_,
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
            dplyr::select(sample = library_id, fastq1, fastq2, strandedness, mouse_id, tissue, timepoint, batch, treatment, genotype, sex, dob)
    } else if (project == "AML.mRNA.2021.RxGroup1") {
        # AML.mRNA.2021.RxGroup1 ├-----------------------------------------------------------
        fastqs <- data.frame(
            fastq1 = system(
                paste0("find '/net/isi-dcnl/ifs/user_data/rrockne/Seq/210412_Tgen' -name '*.gz'"),
                intern = TRUE
            ) %>% grep("_R1_", ., value = TRUE),
            fastq2 = system(
                paste0("find '/net/isi-dcnl/ifs/user_data/rrockne/Seq/210412_Tgen' -name '*.gz'"),
                intern = TRUE
            ) %>% grep("_R2_", ., value = TRUE)
        ) |> dplyr::mutate(library_id = stringr::str_extract(fastq1, "COHP_\\d{5}"))

        xls <- "/net/isi-dcnl/ifs/user_data/rrockne/MHO/AML.mRNA.2021.RxGroup1/config/run_summary_IGC-LZ-18757.xlsx"
        sample_sheet <- readxl::read_excel(xls) |>
            dplyr::mutate(
                library_id = TGen_Sample_Name,
                mouse_id = gsub("_.*", "", name),
                tissue = "PBMC",
                timepoint = gsub("^.*_", "", name),
                batch = "2021_A",
                strandedness = "reverse",
                sex = NA_character_,
                dob = NA_character_,
                treatment = NA_character_,
                genotype = NA_character_,
                tissue = dplyr::case_when(
                    timepoint == "BM" ~ "BM",
                    TRUE ~ tissue
                ),
                timepoint = sub("BM", NA_character_, timepoint)
            ) |>
            dplyr::left_join(fastqs) |>
            dplyr::select(sample = library_id, fastq1, fastq2, strandedness, mouse_id, tissue, timepoint, batch, treatment, genotype, sex, dob)
    } else if (project == "AML.mRNA.2021.RxGroup2") {
        # AML.mRNA.2021.RxGroup2 ├-----------------------------------------------------------
        fastqs <- data.frame(
            fastq1 = system(
                paste0("find '/net/isi-dcnl/ifs/user_data/rrockne/Seq/210507_Tgen' -name '*.gz'"),
                intern = TRUE
            ) %>% grep("_R1_", ., value = TRUE),
            fastq2 = system(
                paste0("find '/net/isi-dcnl/ifs/user_data/rrockne/Seq/210507_Tgen' -name '*.gz'"),
                intern = TRUE
            ) %>% grep("_R2_", ., value = TRUE)
        ) |> dplyr::mutate(library_id = stringr::str_extract(fastq1, "COHP_\\d{5}"))

        xls <- "/net/isi-dcnl/ifs/user_data/rrockne/MHO/AML.mRNA.2021.RxGroup2/config/sequencing summary_IGC-LZ-18862.xlsx"
        sample_sheet <- readxl::read_excel(xls) |>
            dplyr::mutate(
                library_id = TGen_Sample_Name,
                Sample_ID = gsub("-", "", Sample_ID),
                Sample_ID = gsub("PBend", "END", Sample_ID),
                mouse_id = gsub(".*_", "", Sample_ID) |> substr(1, 4),
                tissue = "PBMC",
                timepoint = gsub(".*_\\d{4}", "", Sample_ID),
                batch = "2021_B",
                strandedness = "reverse",
                sex = NA_character_,
                dob = NA_character_,
                treatment = NA_character_,
                genotype = NA_character_,
                tissue = dplyr::case_when(
                    timepoint == "BM" ~ "BM",
                    TRUE ~ tissue
                ),
                timepoint = sub("BM", NA_character_, timepoint)
            ) |>
            dplyr::left_join(fastqs) |>
            dplyr::select(sample = library_id, fastq1, fastq2, strandedness, mouse_id, tissue, timepoint, batch, treatment, genotype, sex, dob)
    } else if (project == "AML.mRNA.2021.RxGroup2_pt2") {
        # AML.mRNA.2021.RxGroup2_pt2 ├-----------------------------------------------------------
        fastqs <- data.frame(
            fastq1 = system(
                paste0("find '/net/isi-dcnl/ifs/user_data/ykuo/Seq/210910_TGen' -name '*.gz'"),
                intern = TRUE
            ) %>% grep("_R1_", ., value = TRUE),
            fastq2 = system(
                paste0("find '/net/isi-dcnl/ifs/user_data/ykuo/Seq/210910_TGen' -name '*.gz'"),
                intern = TRUE
            ) %>% grep("_R2_", ., value = TRUE)
        ) |> dplyr::mutate(library_id = stringr::str_extract(fastq1, "COHP_\\d{5}"))

        sample_sheet <- read.csv("/net/isi-dcnl/ifs/user_data/rrockne/MHO/AML.mRNA.2021.RxGroup2_pt2/config/rename_fastqs.tsv", sep = "\t", header = FALSE) |>
            dplyr::mutate(
                library_id = stringr::str_extract(V1, "COHP_\\d{5}"),
                Sample_ID = gsub("_.*", "", V2),
                mouse_id = gsub("-.*", "", V2),
                tissue = "PBMC",
                timepoint = gsub(".*-", "", Sample_ID),
                batch = "2021_C",
                strandedness = "reverse",
                sex = NA_character_,
                dob = NA_character_,
                treatment = NA_character_,
                genotype = NA_character_,
                tissue = dplyr::case_when(
                    timepoint == "BM" ~ "BM",
                    TRUE ~ tissue
                ),
                timepoint = sub("BM", NA_character_, timepoint)
            ) |>
            dplyr::left_join(fastqs) |>
            dplyr::select(sample = library_id, fastq1, fastq2, strandedness, mouse_id, tissue, timepoint, batch, treatment, genotype, sex, dob)
    } else if (project == "AML.mRNA.2022.RxGroup3") {
        # AML.mRNA.2022.RxGroup3 ├-----------------------------------------------------------
        fastqs <- data.frame(
            fastq1 = system(
                paste0("find '/net/isi-dcnl/ifs/user_data/ykuo/Seq/220415_IGC-LZ-20205' -name '*.gz'"),
                intern = TRUE
            ) %>% grep("_R1_", ., value = TRUE),
            fastq2 = system(
                paste0("find '/net/isi-dcnl/ifs/user_data/ykuo/Seq/220415_IGC-LZ-20205' -name '*.gz'"),
                intern = TRUE
            ) %>% grep("_R2_", ., value = TRUE)
        ) |> dplyr::mutate(library_id = stringr::str_extract(fastq1, "COHP_\\d{5}"))

        xls <- "data-raw/sample summary_IGC-LZ-20205 1.xlsx"
        sample_sheet <- readxl::read_excel(xls) |>
            dplyr::transmute(
                sample_name = gsub(".*_", "", Sample_ID),
                library_id = paste0("COHP_", gsub("_.*", "", Sample_ID))
            ) |>
            dplyr::mutate(
                strandedness = "reverse",
                tissue = "PBMC",
                treatment = NA,
                genotype = NA,
                sex = NA,
                dob = NA,
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
            dplyr::left_join(fastqs) |>
            dplyr::select(sample = library_id, fastq1, fastq2, strandedness, mouse_id, tissue, timepoint, batch, treatment, genotype, sex, dob)
    } else if (project == "AML.mRNA.novaseq_validation.2020") {
        # AML.mRNA.novaseq_validation.2020 ├-----------------------------------------------------------
        ## fastqs generated at COH
        # fastqs <- data.frame(
        #     fastq1 = system(
        #         paste0("find '/net/isi-dcnl/ifs/user_data/ykuo/Seq/200829' -name '*.gz'"),
        #         intern = TRUE
        #     ) %>% grep("_R1_", ., value = TRUE),
        #     fastq2 = system(
        #         paste0("find '/net/isi-dcnl/ifs/user_data/ykuo/Seq/200829' -name '*.gz'"),
        #         intern = TRUE
        #     ) %>% grep("_R2_", ., value = TRUE)
        # ) |> dplyr::mutate(library_id = paste0("IGCP_", substr(basename(fastq1), 1, 5)))
        ## fastqs generated at TGen
        fastqs <- data.frame(
            fastq1 = system(
                paste0("find '/net/isi-dcnl/ifs/user_data/ykuo/Seq/200928_TGen' -name '*.gz'"),
                intern = TRUE
            ) %>% grep("_R1_", ., value = TRUE),
            fastq2 = system(
                paste0("find '/net/isi-dcnl/ifs/user_data/ykuo/Seq/200928_TGen' -name '*.gz'"),
                intern = TRUE
            ) %>% grep("_R2_", ., value = TRUE)
        ) |> dplyr::mutate(library_id = stringr::str_extract(fastq1, "COHP_\\d{5}"))

        sample_sheet <- fastqs |>
            dplyr::mutate(
                sample = library_id,
                mouse_id = NA_character_,
                tissue = "PBMC",
                timepoint = NA_character_,
                batch = "2020_B",
                strandedness = "reverse",
                sex = NA_character_,
                dob = NA_character_,
                treatment = NA_character_,
                genotype = NA_character_
            ) |>
            dplyr::select(sample = library_id, fastq1, fastq2, strandedness, mouse_id, tissue, timepoint, batch, treatment, genotype, sex, dob)
    } else if (project == "AML.validation.2017") {
        # AML.validation.2017 ├-----------------------------------------------------------
        fastqs <- data.frame(
            fastq1 = system(
                paste0("find '/net/isi-dcnl/ifs/user_data/rrockne/MHO/AML.validation.2017/data/fastq' -name '*.gz'"),
                intern = TRUE
            ) %>% grep("_R1_", ., value = TRUE)
        ) |> dplyr::mutate(library_id = paste0("COHP_", substr(basename(fastq1), 1, 5)))

        xls <- "/net/isi-dcnl/ifs/user_data/rrockne/MHO/AML.validation.2017/data/validation_name_key.xlsx"
        sample_sheet <- readxl::read_excel(xls) |>
            dplyr::mutate(
                library_id = paste0("COHP_", substr(Filename, 1, 5)),
                mouse_id = Mouse,
                tissue = "PBMC",
                timepoint = TimePoint,
                batch = "2017_A",
                strandedness = "reverse",
                sex = Sex,
                dob = NA_character_,
                fastq2 = NA_character_,
                treatment = Treatment,
                genotype = Genotype
            ) |>
            dplyr::left_join(fastqs) |>
            dplyr::select(sample = library_id, fastq1, fastq2, strandedness, mouse_id, tissue, timepoint, batch, treatment, genotype, sex, dob)
    } else if (project == "AML.scRNAseq.2022") {
        # AML.scRNA.2022 ├-----------------------------------------------------------
        fastqs <- data.frame(
            fastq1 = system(
                paste0("find '/net/isi-dcnl/ifs/user_data/rrockne/Seq/211210_IGC-LZ-19773' -name '*.gz'"),
                intern = TRUE
            ) %>% grep("_R1_", ., value = TRUE),
            fastq2 = system(
                paste0("find '/net/isi-dcnl/ifs/user_data/rrockne/Seq/211210_IGC-LZ-19773' -name '*.gz'"),
                intern = TRUE
            ) %>% grep("_R2_", ., value = TRUE)
        ) |> dplyr::mutate(library_id = stringr::str_extract(fastq1, "COHP_\\d{5}"))

        xls <- "/net/isi-dcnl/ifs/user_data/rrockne/MHO/AML.scRNA.2022/bulkseq_targets/config/sequencing summary_IGC-LZ-19773.xlsx"
        sample_sheet <- readxl::read_excel(xls) |>
            dplyr::mutate(
                library_id = paste0("COHP_", gsub("_.*", "", Sample_Name)),
                mouse_id = gsub(".*_", "", Sample_Name) |> substr(1, 4),
                timepoint = NA_character_,
                batch = "2022_SC",
                strandedness = "reverse",
                sex = NA_character_,
                dob = NA_character_,
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
            dplyr::select(sample = library_id, fastq1, fastq2, strandedness, mouse_id, tissue, timepoint, batch, treatment, genotype, sex, dob)
    } else if (project == "CML.mRNA.2021") {
        # CML.mRNA.2021 ├-----------------------------------------------------------
        fastqs <- data.frame(
            fastq1 = system(
                paste0("find /net/isi-dcnl/ifs/user_data/rrockne/Seq/210910_TGen /net/isi-dcnl/ifs/user_data/rrockne/Seq/210907_TGen -name '*.gz'"),
                intern = TRUE
            ) %>% grep("_R1_", ., value = TRUE),
            fastq2 = system(
                paste0("find /net/isi-dcnl/ifs/user_data/rrockne/Seq/210910_TGen /net/isi-dcnl/ifs/user_data/rrockne/Seq/210907_TGen -name '*.gz'"),
                intern = TRUE
            ) %>% grep("_R2_", ., value = TRUE)
        ) |> dplyr::mutate(library_id = stringr::str_extract(fastq1, "COHP_\\d{5}"))

        xls <- "/net/isi-dcnl/ifs/user_data/rrockne/MHO/CML.mRNA.2021/config/Lianjun_IGC-LZ-19411.xlsx"
        sample_sheet <- readxl::read_excel(xls) |>
            dplyr::mutate(
                library_id = TGen_Sample_Name,
                mouse_id = gsub("-.*", "", Sample_ID),
                timepoint = gsub(".*-", "", Sample_ID),
                batch = "2021_CML",
                strandedness = "reverse",
                sex = NA_character_,
                dob = NA_character_,
                treatment = dplyr::case_when(
                    stringr::str_detect(mouse_id, regex("473|474|475|476|485|482|483|545")) ~ "A",
                    stringr::str_detect(mouse_id, regex("479|478|480|481|484|477|487|540|547|541|542")) ~ "B",
                    stringr::str_detect(mouse_id, regex("490|491|492|493|494|488|489|546")) ~ "C",
                    stringr::str_detect(mouse_id, regex("502|511|512|513|507|508|505|510|514")) ~ "D",
                    TRUE ~ "Ctrl"
                ),
                genotype = NA_character_,
                tissue = "PBMC"
            ) |>
            dplyr::left_join(fastqs) |>
            dplyr::select(sample = library_id, fastq1, fastq2, strandedness, mouse_id, tissue, timepoint, batch, treatment, genotype, sex, dob)
    } else if (project == "CML.mRNA.2022") {
        # CML.mRNA.2022 ├-----------------------------------------------------------
        fastqs <- data.frame(
            fastq1 = system(
                paste0("find /net/isi-dcnl/ifs/user_data/gmarcucci/Seq/220210_IGC-LZ-19931 /net/isi-dcnl/ifs/user_data/gmarcucci/Seq/220202_IGC-LZ-19931 -name '*.gz'"),
                intern = TRUE
            ) %>% grep("_R1_", ., value = TRUE),
            fastq2 = system(
                paste0("find /net/isi-dcnl/ifs/user_data/gmarcucci/Seq/220210_IGC-LZ-19931 /net/isi-dcnl/ifs/user_data/gmarcucci/Seq/220202_IGC-LZ-19931 -name '*.gz'"),
                intern = TRUE
            ) %>% grep("_R2_", ., value = TRUE)
        ) |> dplyr::mutate(library_id = stringr::str_extract(fastq1, "COHP_\\d{5}"))

        sample_sheet <- read.csv("/net/isi-dcnl/ifs/user_data/rrockne/MHO/CML.mRNA.2022/config/sample_sheet.csv") |>
            dplyr::mutate(
                library_id = stringr::str_extract(fastq_1, "COHP_\\d{5}"),
                mouse_id = gsub("^.(\\d*)-.*", "\\1", sample),
                timepoint = gsub(".*-", "", sample),
                batch = "2022_CML",
                strandedness = "reverse",
                sex = substr(sample, 1, 1),
                dob = NA_character_,
                treatment = dplyr::case_when(
                    stringr::str_detect(mouse_id, regex("473|474|475|476|485|482|483|545")) ~ "TET_OFF_ON_A",
                    stringr::str_detect(mouse_id, regex("479|478|480|481|484|477|487|540|547|541|542")) ~ "TET_ON_B",
                    stringr::str_detect(mouse_id, regex("490|491|492|493|494|488|489|546")) ~ "TET_OFF_C",
                    stringr::str_detect(mouse_id, regex("502|511|512|513|507|508|505|510|514")) ~ "TET_OFF_NIL_ON_D",
                    TRUE ~ "Ctrl"
                ),
                genotype = NA_character_,
                tissue = "PBMC"
            ) |>
            dplyr::left_join(fastqs) |>
            dplyr::select(sample = library_id, fastq1, fastq2, strandedness, mouse_id, tissue, timepoint, batch, treatment, genotype, sex, dob)
        # AML.mRNA.HSA_FLT3.2022 ├-----------------------------------------------------------
    } else if (project == "AML.mRNA.HSA_FLT3.2022") {
        fastqs <- data.frame(
            fastq1 = system(
                paste0("find /net/isi-dcnl/ifs/user_data/ykuo/Seq/220511_IGC-LZ-20342 -name '*.gz'"),
                intern = TRUE
            ) %>% grep("_R1_", ., value = TRUE),
            fastq2 = system(
                paste0("find /net/isi-dcnl/ifs/user_data/ykuo/Seq/220511_IGC-LZ-20342 -name '*.gz'"),
                intern = TRUE
            ) %>% grep("_R2_", ., value = TRUE)
        ) |> dplyr::mutate(library_id = stringr::str_extract(fastq1, "COHP_\\d{5}"))

        xls <- "data-raw/sample summary_IGC-LZ-20342.xlsx"
        sample_sheet <- readxl::read_excel(xls) |>
            dplyr::mutate(
                library_id = TGen_Sample_Name,
                patient_id = gsub(".*_", "", Sample_ID),
                timepoint = NA_character_,
                batch = "2022_HSA_FLT3",
                strandedness = "reverse",
                sex = NA_character_,
                dob = NA_character_,
                treatment = NA_character_,
                genotype = NA_character_,
                tissue = NA_character_
            ) |>
            dplyr::left_join(fastqs) |>
            dplyr::select(sample = library_id, fastq1, fastq2, strandedness, patient_id, tissue, timepoint, batch, treatment, genotype, sex, dob)
    } else {
    stop("Please provide a valid project id:
        AML.mRNA.2016
        AML.mRNA.2018.all_samples
        AML.mRNA.2020
        AML.mRNA.2021.RxGroup1
        AML.mRNA.2021.RxGroup2
        AML.mRNA.2021.RxGroup2_pt2
        AML.mRNA.2022.RxGroup3
        AML.validation.2017
        AML.mRNA.novaseq_validation.2020
        CML.mRNA.2021
        CML.mRNA.2022
        AML.scRNAseq.2022
        AML.mRNA.HSA_FLT3.2022")}

return(sample_sheet)

}
