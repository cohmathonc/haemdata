#' .. content for \description{} (no empty lines) ..
#'
#' .. content for \details{} ..
#'
#' @title Run nf-core/rnaseq pipeline
#' @param run_folder subfolder in which to run nextflow, keeps the pipeline 
#' from overwriting previous runs nad enables the `-resume` feature of nextflow
#' @param sample_sheet nf-core sample sheet with the columns `sample, fastq_1, fastq_2, strandedness` (see \code{\link{parse_metadata_AML()}})
#' @param ref_genome nf-core reference genome, see Reference Genomes for supported references
#' @param qc run the full QC pipeline with STAR_Salmon or just the Salmon pseudoalignment (default: TRUE)
#' @return
#' @author Denis O'Meally
#' @export

run_nf_core_rnaseq <- function(run_folder, sample_sheet, ref_genome, qc = TRUE) {

    run_path <- glue::glue("{here::here()}/nextflow/{run_folder}")
    out_folder <- glue::glue("nfcore-rnaseq-v{release}_{ref_genome}")

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
        genocode <- "--gencode"
    } else {
        genocode <- ""
    }
    if (nrow(sample_sheet) > 20 ) {
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
module load singularity
nextflow run -c {here::here()}/nextflow/nextflow.apollo \\
    -c {here::here()}/nextflow/igenomes.apollo \\
    nf-core/rnaseq -r {release} -resume \\
    --input {run_path}/sample_sheet.csv \\
    --outdir {out_folder} {skip} {deseq_transform} \\
    --genome {ref_genome} --igenomes_base /ref_genome/igenomes \\
    --email domeally@coh.org --pseudo_aligner salmon {genocode}", .trim = FALSE),
            file = glue::glue("{run_path}/run_{run_mode}_{ref_genome}.sh")
        )
        # submit it to the cluster & return the path of the run script
        system(glue::glue("sbatch {run_path}/run_{run_mode}_{ref_genome}.sh"), wait = FALSE)
        return(glue::glue("{run_path}/run_{run_mode}_{ref_genome}.sh"))
    } else {
    # If there's no multiqc report, but there's a run script, just return the run script path
    # with a message about what to do next
        print(glue::glue("The nextflow run for \"{run_folder}\" has already been submitted.
        Check the job queue or navigate to
        {run_path}
        and submit the `run_{run_mode}_{ref_genome}.sh` script to resume the run."))
        return(glue::glue("{run_path}/run_{run_mode}_{ref_genome}.sh"))
    }
}
