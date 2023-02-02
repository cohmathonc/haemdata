# Roxygen markdown for data/

# metadata_mmu -----------------------------------------------------------
#' Minimal metadata for mouse samples
#'
#' [Download](../metadata_mmu_template.xlsx) a template for new samples
#'  | Column | Description |
#'  |---|---|
#'  | `sample_id` | Unique identifier with the format `PSON_XXXX` that identifies a tissue sample with respect to the mouse, tissue type and sample date. Used to match samples across assays. |
#'  | `library_id` | Integrative Genomics Core *CLN* with the format `COHP_XXXXX`. Matches sequence reads, QC metrics, raw counts and metadata for each sequencing library. Not unique as some libraries are split across multiple sequencing runs.|
#'  | `fastq_1` | Full path of fastq1, unique in the table |
#'  | `fastq_2` | Full path of fastq2, unique in the table, or blank for single-end reads |
#'  | `hdf5` | Full path of h5 file, Cell Rangers's HDF5 Feature-Barcode Matrix Format |
#'  | `strandedness` | Sequencing library protocol (`reverse`, `forward` or `unstranded`) |
#'  | `assay` | one of `mRNA`, `miRNA`, `scRNA` |
#'  | `mouse_id` | A unique 3 or 4 digit identifier for each mouse; database held by [Kuo lab](mailto:YKuo@coh.org?subject=Question%20about%20AML%20mice%20from%20PSON) for AML mice and [Zhang lab](mailto:YKuo@coh.org?subject=Question%20about%20AML%20mice%20from%20PSON) for CML. |
#'  | `tissue` | Tissue type: `PBMC` peripheral blood mononuclear cells; `BM` bone marrow aspirant; `BM_CKIT` ckit+ flow sorted `BM`; `PBMC_CKIT` ckit+ flow sorted `PBMC` |
#'  | `timepoint` | Timepoint label; formatting is inconsistent across cohorts |
#'  | `treatment` | Experimental treatment group |
#'  | `genotype` | Mouse genotype |
#'  | `sex` | Sex (`M`/`F`) |
#'  | `dob` | Date of birth, YYYY-MM-DD |
#'  | `cohort` | Logical groups of samples, generally corresponding to experimental cohorts |
#'  | `batch` | Groups of samples from the sequencing run |
#'  | `dod` | Date of death, YYYY-MM-DD. `NA` if the mouse survived the experiment. |
#'  | `sample_date` | Date the sample was collected, YYYY-MM-DD |
#'  | `percent_ckit` | Percentage of c-KIT+ cells, measured by flow cytometry (CD117) |
#'  | `sample_weeks` | Timepoint in weeks (`sample_date - min(sample_date)`). Post treatment chemo samples begin at week 0, pretreatment samples < 0. |
#'  | `age_at_start` | Age at start of the experiment, in weeks (`min(sample_date) - dob`) |
#'  | `age_at_sample` | Age at sample collection, in weeks (`sample_date - dob`) |
#'  | `age_at_end` | Age at end of of the experiment, in weeks (`max(sample_date) - dob`) |
#'  | `qc_pass_mapping` | `mRNA` samples only. `TRUE` if STAR uniquely mapped reads >= `mapping_threshold` (5% by default), or `FALSE` if not. |
#'  | `ref_dim1` | reference dimension 1 - typically PC1; can be UMAP or any other dimension reduction coordinate. |
#'  | `ref_dim1` | reference dimension 2 - typically PC2 |
#'
#' The [`update_metadata_mmu()`](https://github.com/drejom/haemdata/blob/main/scripts/import_metadata.R#L24)
#' function assembles the metadata for all RNAseq libraries from AML and CML mice, by consolidating
#' data scraped from multiple sequencing run sheets, directly from sequencing folders, emails, and so forth.
#' The code is complex and ugly and undoubtedly some errors will have made it through.
#'
#' Raise an [issue on GitHub](https://github.com/drejom/haemdata/issues)
#' to report erroneous or missing records.
#'
#' @name metadata_mmu
#' @docType data
#' @source [`update_metadata_mmu()`](https://github.com/drejom/haemdata/blob/main/scripts/import_metadata.R#L24)
#' @author Denis O'Meally
NULL

# metadata_hsa  -----------------------------------------------------------
#' Minimal metadata for human samples
#'
#' The [`make_metadata_hsa()`](https://github.com/drejom/haemdata/blob/cf03cf0a3eb420a8ee6276c7ec0a9186a55c0e2b/scripts/import_metadata.R#L3)
#' function assembles the metadata for all RNAseq libraries
#' from patient samples. Minimal metadata fields include library_id, fastq_1, fastq_2, strandedness,
#' sample_id, tissue, weeks, timepoint, batch, treatment, genotype, sex, dob, project.
#'
#' For human samples, metadata are sourced from the EGA and supplied excel sheets.
#'
#' #TODO Describe the studies: MDS & COH Biobank FLT3
#'
#' Raise an [issue on GitHub](https://github.com/drejom/haemdata/issues)
#' to report erroneous or missing records.
#'
#' @name metadata_hsa
#' @docType data
#' @source [`make_metadata_hsa()`](https://github.com/drejom/haemdata/blob/cf03cf0a3eb420a8ee6276c7ec0a9186a55c0e2b/scripts/import_metadata.R#L3)
#' @author Denis O'Meally
NULL

# published_pins -----------------------------------------------------------├
#' Datasets published to the Haemdata pin board
#'
#' Lists the name and version of each dataset published to the Haemdata pinboard in the current release.
#' Pins can be retrieved with the [`get_pin()`] function.
#'
#' @format A `data.frame` with the following columns:
#' \describe{
#'   \item{pin_name}{Name of the dataset on the pin board, including the file extension}
#'   \item{version}{Version of the pin to which the dataset belongs}
#' }
#' @name published_pins
#' @docType data
#' @usage data(published_pins)
#' @examples print(published_pins, n = Inf)
NULL
# # hsa_mrna_mds_GENCODEr40_qc -----------------------------------------------------
# se <- get_pin("hsa_mrna_mds_GENCODEr40_qc.rds")
# #' SummarisedExperiment of `r length(se$sample_id)` samples
# #'
# #' A `SummarisedExperiment` generated by the nf-core/rnaseq `r rnaseq_release` pipeline using Salmon
# #' and tximport. See the pipeline [documentation](https://nf-co.re/rnaseq/3.5/output#pseudo-alignment-and-quantification)
# #' for details. The [reference](../articles/genomes.html)
# #' is GRCm38 using Ensembl annotation v81 (July 2015) with human leukemic transcripts (MYH11, BCR & ABL1).
# #'
# #' Only the Salmon *[transcript length bias corrected](https://nf-co.re/rnaseq/3.7/output#pseudo-alignment-and-quantification)* counts are available in this package;
# #' for other datasets, see the [full pipeline output](http://cgt.coh.org/MHO/CML.mRNA.2021/results/rnaseq3.5.GRCm38.HLT/)
# #'
# #' @format A `SummarisedExperiment` object containing `r length(se$sample_id)` samples.
# #' \describe{
# #'   \item{counts}{`mat` Salmon merged gene counts}
# #'   \item{abundance}{`mat` Salmon merged TPM}
# #' }
# #' @source \url{`r transmute(S4Vectors::metadata(se), url = paste0("http://cgt.coh.org/MHO/haemdata-nf-core-cache/", .data$run_folder, "/nfcore-rnaseq-", rnaseq_release, "_", .data$reference_genome))`}
# #' @name `r S4Vectors::metadata(se)$object_name`
# #' @docType data
# #' @keywords cml mrna 2021
# #' @usage get_pin(cml_mrna_2021_GRCm38_HLT)
# NULL
