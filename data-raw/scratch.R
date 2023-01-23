####pivoted metadata table
library(haemdata)
library(tidyverse)
use_pinboard("haemdata")
#metadata<-
get_pin("metadata_mmu.csv") |> dplyr::filter(str_detect(project, "AML")) |>
    select(sample, project, mouse_id, tissue, timepoint, sample_date) |>
    mutate(across(where(is.factor), as.character)) |>
    distinct() |>
    arrange(mouse_id) |>
    pivot_wider(names_from = timepoint, values_from = sample_date, values_fill = "" ) |>
    write_csv("matadata_mmu_pivoted.csv") |>
    openxlsx::write.xlsx("matadata_mmu_pivoted.xlsx", keepNA = TRUE)

#####


library(gtsummary)

dat <- read.csv("inst/extdata/minimal_metadata_samples.csv")

dat |>
    select(tissue, timepoint, batch, treatment, genotype, sex, dob, project) |>
    tbl_summary()

libraries_with_no_mouse_id <-
dat |>
    dplyr::filter(is.na(mouse_id))

# get QC data
dat <- read.csv("/net/nfs-irwrsrchnas01/labs/rrockne/MHO/haemdata/nextflow/AML.mRNA.2016/nfcore-rnaseq-v3.7_GENCODEm28_HLT/multiqc/star_salmon/multiqc_data/multiqc_general_stats.txt", sep = "\t") |>
    dplyr::select(!Sample)

# remove invariant columns - ie, where min=max
dat <- dat[, !apply(dat, MARGIN = 2, function(x) max(x, na.rm = TRUE) == min(x, na.rm = TRUE))]

#run association tests
assoc<-mixed_assoc(dat)

#make a plot
corrplot::corrplot(assoc$assoc_matrix, is.corr = FALSE)

heatmap(assoc$assoc_matrix, na.rm = TRUE)

# Remove highly correlated variables

tmp <- assoc$assoc_matrix
tmp[!lower.tri(tmp)] <- 0

data.new <-
    assoc$assoc_matrix[, !apply(tmp, 2, function(x) any(abs(x) > 0.99, na.rm = TRUE))]

dim(data.new)
[1] 26 23

colnames(data.new) <- gsub(".*generalstats\\.", "", colnames(data.new))

 [1] "samtools.mapped_passed"             "biotype_counts.percent_rRNA"        "dupradar.dupRadar_intercept"        "picard.PERCENT_DUPLICATION"
 [5] "qualimap.5_3_bias"                  "salmon.percent_mapped"              "salmon.num_mapped"                  "samtools.error_rate"
 [9] "samtools.non_primary_alignments"    "samtools.reads_MQ0_percent"         "samtools.raw_total_sequences"       "star.uniquely_mapped_percent"
[13] "star.uniquely_mapped"               "fastqc_raw.percent_duplicates"      "fastqc_raw.percent_gc"              "fastqc_raw.percent_fails"
[17] "fastqc_raw.total_sequences"         "cutadapt.percent_trimmed"           "fastqc_trimmed.percent_duplicates"  "fastqc_trimmed.percent_gc"
[21] "fastqc_trimmed.avg_sequence_length" "fastqc_trimmed.percent_fails"       "fastqc_trimmed.total_sequences"

## try the same, but with human data - do we get the same columns? Nearly - 23 vs 20
# get QC data
dat <- read.csv("/net/nfs-irwrsrchnas01/labs/rrockne/MHO/haemdata/nextflow/AML.mRNA.HSA_FLT3.2022/nfcore-rnaseq-v3.7_GENCODEr40/multiqc/star_salmon/multiqc_data/multiqc_general_stats.txt", sep = "\t") |>
    dplyr::mutate(sample = gsub("_1$|_2$", "", Sample)) |>
    dplyr::group_by(sample) %>%
    dplyr::mutate(
            across(starts_with("FastQC"), ~ replace_na(.x, mean(.x, na.rm = TRUE)))
        )|>
    dplyr::filter(!grepl("_1$|_2$", Sample)) |>
    dplyr::ungroup() |>
    janitor::remove_empty("cols")|>
    dplyr::select(!c(Sample, sample))

# remove invariant columns - ie, where min=max
dat <- dat[, !apply(dat, MARGIN = 2, function(x) max(x, na.rm = TRUE) == min(x, na.rm = TRUE))]

# run association tests
assoc <- mixed_assoc(dat)

# make a plot
corrplot::corrplot(assoc$assoc_matrix, is.corr = FALSE)

heatmap(assoc$assoc_matrix, na.rm = TRUE)

# Remove highly correlated variables

tmp <- assoc$assoc_matrix
tmp[!lower.tri(tmp)] <- 0

data.new <-
    assoc$assoc_matrix[, !apply(tmp, 2, function(x) any(abs(x) > 0.99, na.rm = TRUE))]

dim(data.new)

colnames(data.new) <- gsub(".*generalstats\\.", "", colnames(data.new))
[1] "samtools.mapped_passed"                 "biotype_counts.percent_rRNA"            "dupradar.dupRadar_intercept"            "picard.PERCENT_DUPLICATION"
 [5] "qualimap.5_3_bias"                      "rseqc.proper_pairs_percent"             "salmon.percent_mapped"                  "salmon.num_mapped"
 [9] "samtools.error_rate"                    "samtools.non_primary_alignments"        "samtools.reads_properly_paired_percent" "samtools.reads_MQ0_percent"
[13] "star.uniquely_mapped_percent"           "star.uniquely_mapped"                   "fastqc_raw.percent_fails"               "fastqc_trimmed.percent_duplicates"
[17] "fastqc_trimmed.percent_gc"              "fastqc_trimmed.avg_sequence_length"     "fastqc_trimmed.percent_fails"           "fastqc_trimmed.total_sequences"

## SCONE
library(scone)

dat <- read.csv("/net/nfs-irwrsrchnas01/labs/rrockne/MHO/haemdata/nextflow/AML.mRNA.HSA_FLT3.2022/nfcore-rnaseq-v3.7_GENCODEr40/multiqc/star_salmon/multiqc_data/multiqc_general_stats.txt", sep = "\t") |>
    dplyr::mutate(sample = gsub("_1$|_2$", "", Sample)) |>
    dplyr::group_by(sample) %>%
    dplyr::mutate(
        across(starts_with("FastQC"), ~ replace_na(.x, mean(.x, na.rm = TRUE)))
    ) |>
    dplyr::filter(!grepl("_1$|_2$", Sample)) |>
    dplyr::ungroup() |>
    janitor::remove_constant(na.rm = TRUE) |>
    dplyr::select(!c(Sample))


AML.mRNA.HSA_FLT3.2022_metadata <- left_join(dat, metadata_mmu) |>
    dplyr::select(sample, everything())



#library("future.batchtools")
plan(batchtools_slurm, template = "future.tmpl", resources = list(ncpus = 2, memory = "80G", walltime = "24:00:00"))
demo("mandelbrot", package = "future", ask = FALSE)


### OUTRIDER

A<-AML.mRNA.2016_qc_se_outliers$aberrant_per_sample
B<-AML.all_mice.mRNA_qc_se_outliers$aberrant_per_sample
C <- all_mice.mRNA_qc_se_outliers$aberrant_per_sample


A <- A[names(C)]
B <- B[names(C)]
C <- C[names(C)]

cbind(A, B, C) |>
    `rownames<-`(names(C)) |>
    na.omit() |>
    heatmap()

D <- CML.all_mice.mRNA_qc_se_outliers$aberrant_per_sample
E <- all_mice.mRNA_qc_se_outliers$aberrant_per_sample

D <- D[names(E)]
E <- E[names(E)]

cbind(D,E) |>
    `rownames<-`(names(E)) #|>
    na.omit() |>
    heatmap()


    # devel board
    MHO_haemdata <- pins::board_folder(
        "/net/nfs-irwrsrchnas01/labs/rrockne/MHO/haemdata/",
        versioned = TRUE
    )
    data(published_pins)
    board <- MHO_haemdata


    name <- published_pins[1,1]

    pin_versions(board, name)

    version <- published_pins[1, 2]
    hash <- gsub(".*-", "", version)

    pins::pin_read(board, name, version = version, hash = hash)

pin_list(board)


look ahead

patient_id = stringr::str_extract(V1, "(?<=subject_id\\=)(.*?)(?=[;])"),

test <- "Gwet's AC1 Coefficient
======================
Percent agreement: 0.5666667 Percent chance agreement: 0.4965278
AC1 coefficient: 0.1393103 Standard error: 0.1617574
95 % Confidence Interval: ( -0.2431852 , 1 )
P-value:  0.4176327
"

pc_agreement <- stringr::str_extract(test, "(?<=Percent agreement: )(\\d+.?\\d+)")
pc_chance_agreement <- stringr::str_extract(test, "(?<=Percent chance agreement: )(\\d+.?\\d+)")
ac1_coefficient <-stringr::str_extract(test, "(?<=AC1 coefficient: )(\\d+.?\\d+)")
std_err <- stringr::str_extract(test, "(?<=Standard error: )(\\d+.?\\d+)")
ci_lower <- stringr::str_extract(test, "(?<=Confidence Interval: \\()( -?\\d+.?\\d+)")
ci_upper <- stringr::str_extract(test, "(?<= , )(\\d+.?\\d+)")
p_value <- stringr::str_extract(test, "(?<=P-value:  )(\\d+.?\\d+)")

data.frame(
    pc_agreement = pc_agreement,
    pc_chance_agreement = pc_chance_agreement,
    ac1_coefficient = ac1_coefficient,
    std_err = std_err,
    ci_lower = ci_lower,
    ci_upper = ci_upper,
    p_value = p_value
)

stringr::str_extract(test, "(?<= , )(\\d+.?(\\d+)?)")


5283175 domeally    GENCODEm28_HLT_seurat_rpca                                            PENDING 0:00      12:00:00  N/A                 N/A                 8    400G    fast                 jkaddis
5283144 domeally    GENCODEm28_HLT_seurat_rpca                                            RUNNING 26:20     11:33:40  2022-06-17T15:22:09 2022-06-18T03:22:09 8    400G    all     ppxhpcnode38 jkaddis
5283142 domeally    GRCm38_HLT_seurat_rpca                                         	  RUNNING 33:19     11:26:41  2022-06-17T15:15:10 2022-06-18T03:15:10 8    400G    all     ppxhpcnode37 jkaddis
5283171 domeally    GRCm38_HLT_seurat_rpca


library(rjson)

json <- "/net/nfs-irwrsrchnas01/labs/rrockne/MHO/haemdata-nf-core-cache/MDS.rnaseq.EGAD00001003891/nfcore-rnaseq-v3.7_GENCODEr40/multiqc/star_salmon/multiqc_data/multiqc_data.json"

myData <- fromJSON(file = json)
# Print the result.
head(myData)

myData<-jsonlite::read_json(json)

myData <- jsonlite::fromJSON(json)

mv  MDS.rnaseq.EGAD00001003891 hsa_mrna_mds


x <- x |>
    Seurat::ScaleData() |>
    Seurat::FindVariableFeatures(nfeatures = 1000) |>
    Seurat::RunPCA(npcs = 20) |>
    Seurat::FindNeighbors() |>
    Seurat::FindClusters(resolution = 0.4) |>
    Seurat::FindClusters(resolution = 0.6) |>
    Seurat::FindClusters(resolution = 0.8) |>
    Seurat::FindClusters(resolution = 1.0) |>
    Seurat::FindClusters(resolution = 1.2) |>
    Seurat::RunUMAP(
        reduction = "pca",
        dims = 1:2,
        seed.use = 3
    )




packages<-c("Seurat", "targets")
lapply(packages, require, character.only = TRUE)
targets::tar_load(mmu_10x_2022_1_GENCODEm28_HLT_sct_clust_cc)

DimPlot(mmu_10x_2022_1_GENCODEm28_HLT_sct_clust_cc, reduction = "umap")

DimPlot(mmu_10x_2022_1_GENCODEm28_HLT_sct_clust_cc, reduction = "umap", group.by = "tissue")

DimPlot(mmu_10x_2022_1_GENCODEm28_HLT_sct_clust_cc, reduction = "umap", group.by = "ckit")

DimPlot(mmu_10x_2022_1_GENCODEm28_HLT_sct_clust_cc, reduction = "umap", group.by = "Phase")

#Imports:
w <- renv::dependencies("./R/access_data.R")


#Imports:
x <- renv::dependencies("R")

x |>
    dplyr::pull("Package") |>
    unique() |>
    paste(collapse = ", ")
#Suggests:
y <- renv::dependencies("scripts")

y |>
    dplyr::pull("Package") |>
    unique() |>
    paste(collapse = ", ")


# Combined
z <- renv::dependencies(c("R", "scripts"))

z |>
    dplyr::pull("Package") |>
    unique() |>
    paste(collapse = ", ")





#####
Imports:
    digest,
    dplyr,
    future,
    future.apply,
    ggplot2,
    glmGamPoi,
    glue,
    helpeRs,
    janitor,
    magrittr,
    Microsoft365R,
    PCAtools,
    pins,
    readr,
    rlang (>= 0.4.11),
    S4Vectors,
    sctransform,
    Seurat,
    stats,
    stringr,
    SummarizedExperiment,
    tibble,
    tidyr,
    utils
Suggests:
dplyr, glue, haemdata, Microsoft365R, pins, rlang, readr, stringr, Seurat, clustifyr, future.apply, ggplot2, Matrix, rtracklayer, BiocParallel, ggsci, grid, helpeRs, janitor, OUTRIDER, parallelly, PCAtools, S4Vectors, SummarizedExperiment, tibble, tidyr, utils, future, future.batchtools, targets, readxl, SeuratDisk, hexSticker, magick, usethis

    BiocParallel,
    clustifyr,
    forcats,
    future.batchtools,
    ggsci,
    gtools,
    gtsummary,
    hdf5r,
    hexSticker,
    knitr,
    magick,
    Matrix,
    naniar,
    OUTRIDER,
    parallelly,
    plyr,
    qs,
    rcompanion,
    readxl,
    rmarkdown,
    rtracklayer,
    SeuratDisk,
    SeuratObject,
    targets,
    tidybulk,
    usethis
VignetteBuilder:
    knitr
Remotes:
    drejom/helpeRs,
    mojaveazure/seurat-disk,
    satijalab/sctransform@develop


    # if the summarised expt is from Kim et al, we need to trim down the rows to keep only those genes in
    # the targeted panel, Supplementary Table 2 (83 genes)
    if (!S4Vectors::metadata(summarised_experiment)["object_name"] == "hsa_mrna_kim_GENCODEr40_qc") {
        message("Trimming summarised_experiment to 83 genes")
        target_genes <- readxl::read_excel("data-raw/41598_2020_76933_MOESM8_ESM.xlsx", sheet = "TableS2") |> janitor::clean_names()

        keep_genes <- (rowData(filtered_se)$gene_name %in% target_genes$gene_name)

        filtered_se <- filtered_se[keep_genes, ]
        assay(filtered_se, "counts") |> head()
    }


#' Delete the latest pins on a board
#'
#' @param board The name of the board to delete pins from
#' @param days_to_delete The number of days to delete pins from
#' @export
delete_latest_pins <- function(board, days_to_delete = 7) {
    library(pins)
    # Get the pins on the specified board
    pins <- pin_list(board)

    # Create a list to store the pin IDs that will be deleted
    pin_ids_to_delete <- list()

    # Get the current time
    current_time <- Sys.time()

    # Loop through each pin and add the latest one to the list if it was created within the last week
    for (pin in pins) {
        # Get the metadata for the current pin
        pin_metadata <- pin_meta(board, pin)

        # Check if the pin was created within the last week
        pin_created_time <- as.POSIXct(pin_metadata$created)
        time_difference <- current_time - pin_created_time
        if (time_difference <= days_to_delete) {
            # Add the pin ID to the list
            pin_ids_to_delete <- c(pin_ids_to_delete, pin)
        }
    }

    # Print the list of pin IDs that will be deleted
    print(pin_ids_to_delete |> unlist())

    # Prompt the user to confirm the deletion
    confirm_delete <- readline(prompt = "Do you want to delete the most recent versions of these pins? (y/n) ")

    # If the user confirms, delete the pins
    if (confirm_delete == "y") {
        for (pin_id in pin_ids_to_delete) {
            pin_metadata <- pin_meta(board, pin_id)
            pin_version_delete(board, pin_id, pin_metadata$local$version)
        }
    }
}


prune_all_pins <- function(board, days_to_keep = 7) {
    library(pins)
    # Get the pins on the specified board
    pins <- pin_list(board)

    for (pin in pins) {
        pin_versions_prune(board, pin, days = days_to_delete)
    }
}
#### sample ID PSON

all_mice <- pins::pin_read(
    pins::board_folder(
        "/net/nfs-irwrsrchnas01/labs/rrockne/MHO/haemdata",
        versioned = FALSE
    ),
    "metadata_mmu.csv"
) |>
    purrr::modify_if(is.factor, as.character)


sample_ids <- all_mice |>
    select(mouse_id, tissue, sample_date) |>
    distinct() |>
    arrange(lubridate::date(sample_date)) |>
    group_by(sample_date, tissue, mouse_id) |>
    mutate(sample_id = sprintf("PSON_%04d", cur_group_id()))


all_mice |>
    select(-sample_id) |>
    full_join(sample_ids) |>
    relocate(sample_id) |>
    View()
    
    distinct() |>
    rio::export("data-raw/mouse_sample_metadata.xlsx")
 sample_sheet <- all_mice |>
     left_join(sample_ids, by = c("mouse_id", "tissue", "sample_date")) |>
     head()
    relocate(sample_id) |>
    head()

### FLT3
sample_sheet |>
    select(-c(fastq_1, fastq_2, strandedness)) |>
    relocate(library_id, patient_id, sample_id, sample_date, tissue, everything()) |>
    distinct() |>
    rio::export("data-raw/flt3_sample_metadata.xlsx")


# miRNA

miRNA_sample_sheet <- rbind(
    parse_metadata_AML.miRNA.2016(),
    parse_metadata_AML.miRNA.2018(),
    parse_metadata_AML.miRNA.2020(),
    parse_metadata_AML.miRNA.2021.RxGroup1(),
    parse_metadata_AML.miRNA.2021.RxGroups1and2(),
    parse_metadata_AML.miRNA.2021.RxGroup2_pt2(),
    parse_metadata_AML.miRNA.2022.RxGroup3()
) |> mutate(mouse_id = as.integer(mouse_id))

mRNA_sample_sheet <- pins::pin_read(
    pins::board_folder(
        "/net/nfs-irwrsrchnas01/labs/rrockne/MHO/haemdata",
        versioned = FALSE
    ),
    "metadata_mmu.csv"
) |>
    purrr::modify_if(is.factor, as.character)


combined <- full_join(mRNA_sample_sheet, miRNA_sample_sheet, by = c(mouse_id, tissue, timepoint)) |>
    arrange(mouse_id, tissue, timepoint) |>
    mutate(mouse_id = as.character(mouse_id)) |>
    fill(
        sample_id,
        treatment,
        genotype,
        sex,
        sample_date,
        percent_ckit,
        dob,
        dod,
        sample_weeks,
        age_at_end,
        age_at_start,
        age_at_sample
)

combined |>
    rio::export("data-raw/combined_sample_metadata.xlsx", overwrite = TRUE)


## scCustomise

scCustomize::Read_Metrics_10X(
    base_path,
    secondary_path = NULL,
    default_10X = TRUE,
    lib_list = NULL,
    lib_names = NULL
)



nextflow run \
    -c /net/nfs-irwrsrchnas01/labs/rrockne/MHO/haemdata-nf-core-cache/nextflow.apollo \
    nf-core/smrnaseq -r 2.1.0 -profile test,singularity \
    --email domeally@coh.org 

# Shared sample dates
#2015-11-06
12  PSON_0073 COHP_11208  T4       /labs/ykuo/Seq/160122/miRNA/11208_2684-5_CGATGT_L999_cutadapt.fastq.gz
13  PSON_0073 COHP_11536  T5       /labs/ykuo/Seq/160304/miRNA/orig/11536_2684-6_TAGCTT_L999_R1_001.fastq.gz

#2015-12-25
76  PSON_0062 COHP_10878 T3     /labs/ykuo/Seq/151219/Jing_Qi_miRNA/orig/10878_2709-4_CCAACA_L999_R1_001.fastq.gz
77  PSON_0062 COHP_11216 T4     /labs/ykuo/Seq/160122/miRNA/11216_2709-5_TAGCTT_L999_cutadapt.fastq.gz

#TODO MISSING METADATA MIRNA
1
/labs/ykuo/Seq/160122/miRNA/11216_2709-5_TAGCTT_L999_cutadapt.fastq.gz
COHP_11216
PSON_0852
2
/labs/ykuo/Seq/160304/miRNA/orig/11536_2684-6_TAGCTT_L999_R1_001.fastq.gz
COHP_11536
PSON_0853


Haemdata is composed of two parts: a targets pipeline for preprocessing raw sequencing reads and functions for accessing the processed data. The pipeline is designed to be efficient and reproducible, allowing users to focus on their data analysis rather than the tedious task of data preprocessing. With Haemdata, you can easily import and manipulate large sequencing datasets, streamlining your workflow and allowing you to quickly answer important biological questions. Give Haemdata a try and see how it can improve your research.
