# Utility functions for the package

# working with pins -----------------------------------------------------------
#' Get an object from the pin_board
#'
#' @param pin_name The name of the pin to retrieve
#' @param version The version of the pin to retrieve; if NULL, the latest version is used. Default: NULL
#' @seealso
#'  \code{\link[pins]{pin_read}}
#' @return An object from the pin_board
#' @export

get_data <- function(pin_name, version = NULL) {

    available_pins <- get_pin_list()
    if (pin_name %in% available_pins) {

    if(is.null(version)) {
        version <- get_latest_pin_version(pin_name)
    }
    hash <- gsub(".*-", "", version)

    pins::pin_read(pin_board, pin_name, version = version, hash = hash)
} else {
    stop("Pin not found in pin_board; use get_pin_list() to see available pins")
}
}

#' @title get_pin_list
#' @description Get names of all available Heamdata pin_board pins
#' @return a vector of pin names
#' @details The `pin_board` is set at the top of in `_targets.R`
#' @examples
#' \dontrun{
#' if (interactive()) {
#'     get_pin_list()
#' }
#' }
#' @seealso
#'  \code{\link[pins]{pin_list}}
#' @rdname get_pin_list
#' @export
get_pin_list <- function() {
    pin_board |>
        pins::pin_list()
}

# run nf-core pipelines --------------------------------------------------------
#' Run nf-core/rnaseq pipeline on Apollo
#'
#' Makes an `sbatch` job script for the nf-core/rnaseq pipeline and submits it to the Apollo cluster.
#' The pipeline version is specified by the `rnaseq_release` parameter in the `_targets.R` file. The [default
#' parameters](https://nf-co.re/rnaseq/3.7/usage) are used for the pipeline, in addition to `--pseudo_aligner salmon` and `--deseq2_vst`.
#'
#' Execution of the pipeline is independent of the package build script, allowing asynchronous development.
#' If no `run_folder` exists, it will be created along with a `sample_sheet.csv` and `run_{ref_genome}.sh`,
#' which is submitted via `sbatch`. A successful run is detected by the presence of a a multicqc report.
#' An existing run script wont be overwritten by this function so it must be deleted manually.
#'
#' @name run_nf_core_rnaseq
#' @param run_folder subfolder in which to run nextflow, keeps the pipeline
#' from overwriting previous runs and enables the `-resume` feature of Nextflow.
#' @param sample_sheet nf-core sample sheet with the columns `sample, fastq_1, fastq_2, strandedness` (see [`R/import_metadata.R`](https://github.com/drejom/haemdata/blob/HEAD/R/import_metadata.R))
#' @param ref_genome nf-core reference genome, see \href{../articles/genomes.html}{Reference Genomes} for supported references
#' @param qc run the full QC pipeline with STAR_Salmon, or Salmon pseudoalignment only, skipping all QC steps (default: TRUE)
#' @return a path to the run script, or the multiqc report if it exists.
#' @author Denis O'Meally
#' @export

run_nf_core_rnaseq <- function(run_folder, sample_sheet, ref_genome, qc = TRUE) {
    run_path <- glue::glue("{nf_core_cache}/{run_folder}")
    out_folder <- glue::glue("nfcore-rnaseq-v{rnaseq_release}_{ref_genome}")

    if (!dir.exists(run_path)) {
        dir.create(run_path, recursive = TRUE)
    }
    if (qc == TRUE) {
        run_mode <- "qc"
        skip <- ""
    } else {
        run_mode <- "salmon"
        skip <- "--skip_alignment --skip_fastqc"
    }
    if (stringr::str_detect(ref_genome, "GENCODE")) {
        gencode <- "--gencode"
    } else {
        gencode <- ""
    }
    if (nrow(sample_sheet) > 20) {
        deseq_transform <- "--deseq2_vst"
    } else {
        deseq_transform <- ""
    }

    # if a multicqc report exists, use that
    if (file.exists(glue::glue("{run_path}/{out_folder}/multiqc/multiqc_report.html"))) {
        print(glue::glue("Found existing multiqc report, skipping nf-core/rnaseq run for {run_folder}..."))
        return(glue::glue("{run_path}/{out_folder}/multiqc/multiqc_report.html"))
    } else if (file.exists(glue::glue("{run_path}/{out_folder}/multiqc/star_salmon/multiqc_report.html"))) {
        print(glue::glue("Found existing multiqc report, skipping nf-core/rnaseq run for {run_folder}..."))
        return(glue::glue("{run_path}/{out_folder}/multiqc/star_salmon/multiqc_report.html"))
        # if theres no run script, make one and submit it to the cluster
    } else if (!file.exists(glue::glue("{run_path}/run_{run_mode}_{ref_genome}.sh"))) {
        sample_sheet |>
            dplyr::select(sample, fastq_1, fastq_2, strandedness) |>
            readr::write_csv(paste0(run_path, "/sample_sheet.csv"))

        # make a sbatch script
        cat(glue::glue("#!/bin/bash
#SBATCH --job-name={run_folder}_{run_mode}
#SBATCH --time=96:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=8G
#SBATCH --export=NONE
#SBATCH --get-user-env=L
#SBATCH -e {run_path}/slurm-%j.err
#SBATCH -o {run_path}/slurm-%j.out

cd {run_path}
module load singularity/3.7.0
module load Java/13.0.2
export NXF_ANSI_LOG=FALSE

nextflow run \\
    -c {nf_core_cache}/nextflow.apollo \\
    -c {nf_core_cache}/igenomes.apollo \\
    nf-core/rnaseq -r {rnaseq_release} -resume \\
    --input {run_path}/sample_sheet.csv \\
    --outdir {out_folder} {skip} {deseq_transform} \\
    --genome {ref_genome} --igenomes_base /ref_genome/igenomes \\
    --email domeally@coh.org --pseudo_aligner salmon {gencode}
", .trim = FALSE),
            file = glue::glue("{run_path}/run_{run_mode}_{ref_genome}.sh")
        )
        # submit to the cluster & return the path of the run script
        system(glue::glue("sbatch {run_path}/run_{run_mode}_{ref_genome}.sh"), wait = FALSE)
        return(glue::glue("{run_path}/run_{run_mode}_{ref_genome}.sh"))
    } else {
        # If there's no multiqc report but there's a run script, just return the run script path
        # with a message about what to do next
        print(glue::glue("The nextflow run for \"{run_folder}\" has already been submitted.
        Check the SLURM job queue or navigate to
        {run_path}
        and submit the `run_{run_mode}_{ref_genome}.sh` script to resume the run."))
        return(glue::glue("{run_path}/run_{run_mode}_{ref_genome}.sh"))
    }
}

# package logo -----------------------------------------------------------
#' Make the package logo
#'
#' makes a hex sticker given a url
#' @param url a URL to a picture
#' @param text the name on the sticker
#' @return a hex sticker
#' @export

# Hex sticker
# https://nelson-gon.github.io/12/06/2020/hex-sticker-creation-r/
make_logo <- function(url=NULL, text=NULL) {
    if (is.null(url)) {
        img <- magick::image_read("https://www.maxpixel.net/static/photo/2x/Blood-Group-0-Blood-Rh-factor-Positive-2781421.jpg")
    } else {
        img <- magick::image_read(url)
    }

    if (is.null(text)) {
        text <- "haemdata"
    }

    # make the logo 600dpi
    img |>
        magick::image_convert("png") |>
        magick::image_resize("900 x 900") |>
        hexSticker::sticker(
            package = text,
            p_size = 25, p_y = 1,
            s_x = 1, s_y = 0.8,
            s_width = 10, s_height = 10,
            white_around_sticker = TRUE,
            h_color = "#FE2B3F",
            filename = "data-raw/haemdata_icon.png",
            dpi = 600
        )

    # use the logo
    usethis::use_logo("data-raw/haemdata_icon.png", retina = TRUE)

    return(img)
}
