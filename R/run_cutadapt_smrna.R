#' Run cutadapt for COH miRNAseq "homebrew" on Apollo
#'
#' Makes an `sbatch` job script for cutadapt and submits it to the Apollo cluster, one job per library.
#'
#' Execution of the pipeline is independent of the package build script, allowing asynchronous development.
#' If no `run_folder` exists, it will be created along with a `sample_sheet.csv` and `run_{ref_genome}.sh`,
#' which is submitted via `sbatch`. A successful run is detected by the presence of a a multiqc report.
#' An existing run script wont be overwritten by this function so it must be deleted manually.
#'
#' @name run_nf_core_rnaseq
#' @param sample_sheet nf-core sample sheet with the columns `sample_id, library_id, fastq_1`
#' @return a data.frame with the columns `library_id, fastq_1`, ready to be used as input to the nf-core pipeline
#' @author Denis O'Meally
#' @export
run_cutadapt_smrna <- function(sample_sheet) {
    run_folder <- rlang::hash(sample_sheet)
    run_path <- glue::glue("{nf_core_cache}/cutadapt_fastq/{run_folder}")

    if (!dir.exists(run_path)) {
        dir.create(run_path, recursive = TRUE)
    }

    # Test that cutadapt has run already - note this only check for the presence of the run script
    if (file.exists(glue::glue("{run_path}/run_cutadapt.sh"))) {
    return(data.frame(
        fastq_1 = list.files(run_path, pattern = "*.fq.gz", full.names = TRUE, include.dirs = TRUE)
    ) |> mutate(library_id = paste0("COHP_", stringr::str_replace(fastq_1, "(?:.*/)(\\d+)_.*", "\\1")), .before = fastq_1))

    # # if theres no run script, make one and submit it to the cluster
    } else if (!file.exists(glue::glue("{run_path}/run_cutadapt.sh"))) {
        sample_sheet |>
            dplyr::select(fastq_1) |>
            readr::write_delim(paste0(run_path, "/untrimmed_file_paths.txt"), col_names = FALSE)

        # make a sbatch script
        cat(glue::glue("#!/bin/bash
#SBATCH --job-name=cutadapt-{run_folder}
#SBATCH --time=96:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=8G
#SBATCH --export=NONE
#SBATCH --get-user-env=L
#SBATCH -e {run_path}/slurm-%j.err
#SBATCH -o {run_path}/slurm-%j.out

source /opt/modulefiles/bashrc.local
PATH=~/bin:$PATH

cd {run_path}
module load cutadapt/1.18-foss-2018b-Python-3.6.6

# read the filepaths into an array
filepaths=()
while read -r line; do
    filepaths+=(\"$line\")
done < <(cat untrimmed_file_paths.txt| grep -v \"^$\")

# submit a job for each filepath
for filepath in ${{filepaths[@]}}; do
    basename=$(basename $filepath | cut -d. -f1)
    sbatch <<EOF
#!/bin/bash
#SBATCH --job-name=cutadapt_$basename
#SBATCH --output=cutadapt_$basename.out
#SBATCH --error=cutadapt_$basename.err
#SBATCH --cpus-per-task=6


cutadapt -j 6 -a TCTGGAATTCTCGGGTGCCAAGGAACTCC -m 16 -u 3 -o {run_path}/$basename.cutadapt.fq.gz \"$filepath\" > {run_path}/$basename.cutadapt.txt

EOF
done
", .trim = FALSE),
            file = glue::glue("{run_path}/run_cutadapt.sh")
        )
        # submit to the cluster & return the path of the run script
        system(glue::glue("sbatch {run_path}/run_cutadapt.sh"), wait = FALSE)
        return(glue::glue("{run_path}/run_cutadapt.sh"))
    } else {
        # If job done report but there's a run script, just return the run script path
        # with a message about what to do next
        print(glue::glue("The sbatch script for \"{run_folder}\" has already been submitted.
        Check the SLURM job queue or navigate to
        {run_path}
        and submit the `run_cutadapt.sh` script to resume the run."))
        return(glue::glue("{run_path}/run_cutadapt.sh"))
    }
}