# Roxygen markdown for data/

# metadata_mmu -----------------------------------------------------------
#' Minimal metadata for mouse samples
#'
#'  |              | Description                                                                                                                                           |
#'  |--------------|-------------------------------------------------------------------------------------------------------------------------------------------------------|
#'  | `sample`       | Integrative Genomics Core `library_id` with the format `COHP_XXXXX`. We use this unique identifier to match sequence reads, QC metrics, raw counts and metadata, for each sample|
#'  | `fastq_1`      | Full path of fastq1                                                                                                                                   |
#'  | `fastq_2`      | Full path of fastq2, or blank for single-end reads                                                                                                     |
#'  | `strandedness` | Sequencing library protocol (`reverse` or `unstranded`)                                                                                               |
#'  | `mouse_id`     | A unique identifier for each mouse, a 3 or 4 digit number; database held by Kuo lab                                                                   |
#'  | `tissue`       | Tissue type of the sample being sequenced: `PBMC` peripheral blood mononuclear cells; `BM` bone marrow aspirant; `BM_CKIT` ckit+ flow sorted `BM`|
#'  | `weeks`        | Time in weeks a sample was taken                                                            |
#'  | `timepoint`    | Timepoint label; formatting is inconsistent and needs to be harmonised                                                             |
#'  | `batch`        | Sequencing groups                                                                                                      |
#'  | `treatment`    | Experimental treatment groups                                                                                                       |
#'  | `genotype`     | Mouse's genotype                                                                                                                |
#'  | `sex`          | Sex of the mouse                                                                                                                                    |
#'  | `dob`          | Date of birth of the mouse                                                                                                                          |
#'  | `project`      | Indicates which samples were processed together; a project may contain multiple `batches`                                                             |
#'
#' The [`make_metadata_mmu()`] function assembles the metadata for all RNAseq libraries from AML and CML mice.
#' Minimal metadata fields include sample, fastq_1, fastq_2, strandedness,
#' mouse_id, tissue, weeks, timepoint, batch, treatment, genotype, sex, dob, project.
#'
#' Raise an [issue on GitHub](https://github.com/drejom/haemdata/issues)
#' to report erroneous or missing records.
#'
#' @name metadata_mmu
#' @docType data
#' @source [`make_metadata_mmu()`]
#' @author Denis O'Meally
NULL

# metadata_hsa  -----------------------------------------------------------
#' Minimal metadata for human samples
#'
#' The [`make_metadata_hsa()`] function assembles the metadata for all RNAseq libraries
#' from patient samples. Minimal metadata fields include sample, fastq_1, fastq_2, strandedness,
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
#' @source [`make_metadata_hsa()`]
#' @author Denis O'Meally
NULL

# published_pins -----------------------------------------------------------├
#' Datasets published to the Haemdata pin board
#'
#' Lists the name and most recent version of each dataset published to the Haemdata pin board.
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
#' @examples published_pins
NULL