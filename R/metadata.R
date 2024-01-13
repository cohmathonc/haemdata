# Roxygen markdown for data/

# metadata_mmu -----------------------------------------------------------
#' @title Minimal metadata for mouse samples
#'
#' @description 
#'  [Download](../metadata_mmu_template.xlsx) a template for new samples
#'  | Column | Description |
#'  |---|---|
#'  | `sample_id` | Unique identifier with the format `MHO_XXXX` that identifies a tissue sample with respect to the mouse, tissue type and sample date. Used to match samples across assays. |
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
#'  | `batch` | An arbitrary string to indicate samples in the same sequencing run; unique within cohorts |
#'  | `dod` | Date of death, YYYY-MM-DD. `NA` if the mouse survived the experiment. |
#'  | `sample_date` | Date the sample was collected, YYYY-MM-DD |
#'  | `percent_ckit` | Percentage of c-KIT+ cells, measured by flow cytometry (CD117) |
#'  | `sample_weeks` | Timepoint in weeks (`sample_date - min(sample_date)`). Post treatment chemo samples begin at week 0, pretreatment samples < 0. |
#'  | `age_at_start` | Age at start of the 1st sample, in weeks (`min(sample_date) - dob`) |
#'  | `dead` | Dead or alive at time of sampling; `0` = alive, `1` = dead |
#'  | `ref_dim1` | reference dimension 1 - typically PC1; can be UMAP or any other dimension reduction coordinate. |
#'  | `ref_dim2` | reference dimension 2 - typically PC2 |
#'  | `qc_pass_mapping` | `mRNA` samples only. `TRUE` if STAR uniquely mapped reads >= `mapping_threshold` (5% by default), or `FALSE` if not. |
#'
#' @details
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
#' @examples published_pins |> knitr::kable()
NULL

#' Plot cohort survival
#'
#' @param sample_sheet A data frame containing sample metadata.
#' @param cohort_regex A regular expression to filter the cohort of interest.
#' @param assay_regex A regular expression to filter the assay of interest.
#'
#' @details Uses the `survfit` function from the {survival} package to fit a
#' survival model to the data, and uses the `autoplot.survfit` function from
#' the {ggfortify} package to create a ggplot to visualize the fitted model.
#' The plot shows the cohort survival over time, with different colors indicating
#' different treatments.
#'
#' @return A ggplot object showing the effect of treatment on survival for a given cohort and assay.
#'
#' @examples
#' \dontrun{
#' use_pinboard("devel")
#' get_pin("metadata_mmu.csv") |>
#' plot_cohort_survival(cohort_regex = "AML")}
#' @export
plot_cohort_survival <- function(sample_sheet = metadata_mmu_prepub,
                                cohort_regex,
                                assay_regex = "mRNA") {
    cohort_sample_sheet <- sample_sheet |>
        dplyr::filter(grepl({{ cohort_regex }}, cohort)) |>
        dplyr::filter(grepl({{ assay_regex }}, assay)) |>
        dplyr::select(sample_id, sample_weeks, dead, treatment)

    if (nrow(cohort_sample_sheet) == 0) {
        print("No samples to plot...")
        return()
    }

    ggfortify:::autoplot.survfit(
        survival::survfit(
            survival::Surv(as.numeric(sample_weeks), dead) ~ treatment,
            data = cohort_sample_sheet
        ),
        xlab = "Weeks",
        ylab = "Cohort survival",
        main = "Effect of treatment on survival") +
        ggplot2::labs(
            fill = "Treatment",
            color = "Treatment",
        subtitle = glue::glue("Cohort: {cohort_regex}"))
}
