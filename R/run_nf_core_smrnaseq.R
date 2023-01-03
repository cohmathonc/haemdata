#' Run nf-core/smrnaseq pipeline on Apollo
#'
#' Makes an `sbatch` job script for the nf-core/smrnaseq pipeline and submits it to the Apollo cluster.
#' The pipeline version is specified by the `smrnaseq_release` parameter in the `_targets.R` file. The [default
#' parameters](https://nf-co.re/smrnaseq/2.1.0/usage) are used for the pipeline.
#'
#' Execution of the pipeline is independent of the package build script, allowing asynchronous development.
#' If no `run_folder` exists, it will be created along with a `sample_sheet.csv` and `run_{ref_genome}.sh`,
#' which is submitted via `sbatch`. A successful run is detected by the presence of a a multiqc report.
#' An existing run script wont be overwritten by this function so it must be deleted manually.
#'
#' @name run_nf_core_rnaseq
#' @param run_folder subfolder in which to run nextflow, keeps the pipeline
#' from overwriting previous runs and enables the `-resume` feature of Nextflow.
#' @param sample_sheet nf-core sample sheet with the columns `library_id, fastq_1, fastq_2, strandedness` (see [`R/import_metadata.R`](https://github.com/drejom/haemdata/blob/HEAD/R/import_metadata.R))
#' @param species mmu or hsa. If *mmu*, additionally sets the 3' adapter to TCTGGAATTCTCGGGTGCCAAGGAACTCC. For *hsa*, teh 3' adapter is auto-discovered.
#' @return a path to the run script, or the multiqc report if it exists.
#' @author Denis O'Meally
#' @export

run_nf_core_smrnaseq <- function(run_folder, sample_sheet, species = "mmu") {
    run_path <- glue::glue("{nf_core_cache}/{run_folder}")
    out_folder <- glue::glue("nfcore-smrnaseq-v{smrnaseq_release}-{species}")

    if (!dir.exists(run_path)) {
        dir.create(run_path, recursive = TRUE)
    }
    if (species == "mmu") {
            ref_genome <- "GRCm38"
            three_prime_adapter <- "--protocol 'custom' --three_prime_adapter 'TCTGGAATTCTCGGGTGCCAAGGAACTCC'"
        } else if (species == "hsa") {
            ref_genome <- "GRCh38"
            three_prime_adapter <- "--protocol 'custom' --three_prime_adapter ''"
        } else {
        stop("Species not supported")
        }

    # if a multiqc report exists, use that
    if (file.exists(glue::glue("{run_path}/{out_folder}/multiqc/multiqc_report.html"))) {
        print(glue::glue("Found existing multiqc report, skipping nf-core/smrnaseq run for {run_folder}..."))
        return(glue::glue("{run_path}/{out_folder}/multiqc/multiqc_report.html"))
    } else if (file.exists(glue::glue("{run_path}/{out_folder}/multiqc/star_salmon/multiqc_report.html"))) {
        print(glue::glue("Found existing multiqc report, skipping nf-core/smrnaseq run for {run_folder}..."))
        return(glue::glue("{run_path}/{out_folder}/multiqc/star_salmon/multiqc_report.html"))
        # if theres no run script, make one and submit it to the cluster
    } else if (!file.exists(glue::glue("{run_path}/run_{ref_genome}.sh"))) {
        sample_sheet |>
            dplyr::select(sample = sample_id, fastq_1) |>
            readr::write_csv(paste0(run_path, "/sample_sheet.csv"))

        # make a sbatch script
        cat(glue::glue("#!/bin/bash
#SBATCH --job-name={run_folder}_{species}
#SBATCH --time=96:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=8G
#SBATCH --export=NONE
#SBATCH --get-user-env=L
#SBATCH -e {run_path}/slurm-%j.err
#SBATCH -o {run_path}/slurm-%j.out

source /opt/modulefiles/bashrc.local
PATH=~/bin:$PATH

cd {run_path}
module load singularity/3.7.0
module load Java/13.0.2
export NXF_ANSI_LOG=FALSE

nextflow run \\
    -c {nf_core_cache}/nextflow.apollo \\
    -c {nf_core_cache}/igenomes.apollo \\
    nf-core/smrnaseq -r {smrnaseq_release} -resume \\
    --input {run_path}/sample_sheet.csv \\
    --outdir {out_folder} \\
    --genome {ref_genome} --igenomes_base /ref_genome/igenomes \\
    --email domeally@coh.org {three_prime_adapter} 
", .trim = FALSE),
            file = glue::glue("{run_path}/run_{ref_genome}.sh")
        )
        # submit to the cluster & return the path of the run script
        system(glue::glue("sbatch {run_path}/run_{ref_genome}.sh"), wait = FALSE)
        return(glue::glue("{run_path}/run_{ref_genome}.sh"))
    } else {
        # If there's no multiqc report but there's a run script, just return the run script path
        # with a message about what to do next
        print(glue::glue("The nextflow run for \"{run_folder}\" has already been submitted.
        Check the SLURM job queue or navigate to
        {run_path}
        and submit the `run_{ref_genome}.sh` script to resume the run."))
        return(glue::glue("{run_path}/run_{ref_genome}.sh"))
    }
}