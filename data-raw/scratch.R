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

