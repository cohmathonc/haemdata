#' Get the counts table from mirtop output
#'
#' Given the path of a multiqc pipeline report file, this function will read in `mirna.tsv` from the miRTop folder and return a tibble of counts.
#'
#' @title nfcore_mirtop_counts
#' @param multiqc_path the full path to a multiqc pipeline report from a nfcore/smrnaseq run
#' on COH Isilon storage
#' @return a tibble of counts, with genes as rows, and samples as columns
#' @author Denis O'Meally
#' @export
nfcore_mirtop_counts <- function(multiqc_path) {

  # Test if the multiqc is from a smrnaseq run
  if (!stringr::str_detect(multiqc_path, "nfcore-smrnaseq")) stop("The multiqc report is not from a nf-core/smrnaseq run")

  mirtop_counts_tsv <- gsub("multiqc/multiqc_report.html", "mirtop/mirna.tsv", multiqc_path)

  mirtop_counts <- readr::read_delim(mirtop_counts_tsv, show_col_types = FALSE) |>
    dplyr::rename_all(~ gsub("_seqcluster", "", .))

  return(mirtop_counts)

}

#' Get the mirTrace and samtools QC tables
#'
#' Given the path of a multiqc pipeline report file, this function will read in mirTrace and samtools stats from the multiqc folder and return a tibble.
#'
#' @title nfcore_smrna_qc
#' @param multiqc_path the full path to a multiqc pipeline report from a nfcore/smrnaseq run
#' on COH Isilon storage
#' @return a tibble of counts, with samples as rows, and QC metrics as columns
#' @author Denis O'Meally
#' @export
nfcore_smrna_qc <- function(multiqc_path) {
  # Test if the multiqc is from a smrnaseq run
  if (!stringr::str_detect(multiqc_path, "nfcore-smrnaseq-v2.1.0-mmu")) stop("The multiqc report is not from a nf-core/smrnaseq v2.1.0 run")

  mirtrace_summary_txt <- gsub("multiqc_report.html", "multiqc_data/multiqc_mirtrace_summary.txt", multiqc_path)
  mirtrace_contamination_txt <- gsub("multiqc_report.html", "multiqc_data/multiqc_mirtrace_contamination.txt", multiqc_path)
  samtools_summary_txt <- gsub("multiqc_report.html", "multiqc_data/multiqc_samtools_stats.txt", multiqc_path)

  mirtrace_summary <- readr::read_delim(mirtrace_summary_txt, show_col_types = FALSE) |>
    mutate(library_id = paste0("COHP_", gsub("_.*", "", Sample)), .before = dplyr::everything()) |>
    select(-c(Sample, filename))

  mirtrace_contamination <- readr::read_delim(mirtrace_contamination_txt, show_col_types = FALSE) |>
    mutate(library_id = paste0("COHP_", gsub("_.*", "", Sample)), .before = dplyr::everything()) |>
    select(-c(Sample))

  mirtrace_qc <- left_join(mirtrace_summary, mirtrace_contamination, by = "library_id") %>%
    setNames(paste0("mirtrace_", names(.))) |>
    dplyr::rename(library_id = "mirtrace_library_id")

  samtools_mature_summary <- readr::read_delim(samtools_summary_txt, show_col_types = FALSE) |>
    dplyr::filter(!grepl("hairpin", Sample)) |>
    mutate(Sample = gsub("_mature", "", Sample)) %>%
    setNames(paste0("samtools_", names(.))) |>
    dplyr::rename(library_id = "samtools_Sample")

  smrnaseq_qc <- left_join(mirtrace_qc, samtools_mature_summary, by = "library_id")

  return(smrnaseq_qc)
}
#' Get the counts table from edgeR CPM file
#'
#' Given the path of a multiqc pipeline report file, this function will read in `hairpin_normalized_CPM.txt` from the edger folder and return a tibble of counts per million.
#'
#' @title nfcore_smrna_cpm
#' @param multiqc_path the full path to a multiqc pipeline report from a nfcore/smrnaseq run
#' on COH Isilon storage
#' @return a data.frame of counts per million, with genes as rows, and samples as columns
#' @author Denis O'Meally
#' @export
nfcore_smrna_cpm <- function(multiqc_path) {
  # Test if the multiqc is from a smrnaseq run
  if (!stringr::str_detect(multiqc_path, "nfcore-smrnaseq")) stop("The multiqc report is not from a nf-core/smrnaseq run")

  edger_norm_cpm_txt <- gsub("multiqc/multiqc_report.html", "edger/mature_normalized_CPM.txt", multiqc_path)

  edger_norm_cpm <- read.table(edger_norm_cpm_txt, header = TRUE, row.names = 1, sep = "\t")

  return(edger_norm_cpm)
}

#' Get the isomiRs table from mirtop_rawData.tsv file
#'
#' Given the path of a multiqc pipeline report file, this function will read in `mirtop.tsv`
#' from the mirtop folder and return a tibble mirtop isomiRs suitable for input to the
#' [`isomiRs`](https://www.bioconductor.org/packages/release/bioc/html/isomiRs.html) package.
#'
#' @title nfcore_mirtop_isomir
#' @param multiqc_path the full path to a multiqc pipeline report from a nfcore/smrnaseq run
#' on COH Isilon storage
#' @return a data.frame of counts per million, with genes as rows, and samples as columns
#' @author Denis O'Meally
#' @export
nfcore_mirtop_isomir <- function(multiqc_path) {
  # Test if the multiqc is from a smrnaseq run
  if (!stringr::str_detect(multiqc_path, "nfcore-smrnaseq")) stop("The multiqc report is not from a nf-core/smrnaseq run")

  mirtop_tsv <- gsub("multiqc/multiqc_report.html", "mirtop/mirtop_rawData.tsv", multiqc_path)

  mirtop <- readr::read_tsv(mirtop_tsv, show_col_types = FALSE) |>
    dplyr::rename_all(~ gsub("_seqcluster", "", .))

  return(mirtop)
}
#' Make an isomiRs object from the mirtop isomir table and a metadata table
#'
#' Given a mirtop isomir table and a metadata table, this will return an
#' [`isomiRs`](https://www.bioconductor.org/packages/release/bioc/html/isomiRs.html)
#' SummarizedExperiment object.
#'
#' @title make_isomir
#' @param mirtop_isomirs a data.frame containing the mirtop isomiR table
#' @param metadata a data.frame containing the colData to bind
#' @return an IsomirDataSeq object (a subclass of SummarizedExperiment)
#' @author Denis O'Meally
#' @export
make_isomir <- function(mirtop_isomirs, metadata) {

  keep <- grep("COHP", names(mirtop_isomirs), value = TRUE)

  metadata_df <- metadata |>
    dplyr::filter(library_id %in% keep) |>
    tibble::column_to_rownames("library_id") |>
    select(sample_id, mouse_id, tissue, cohort, batch, treatment, genotype, sex, percent_ckit, sample_weeks)

  isomir <- isomiRs::IsomirDataSeqFromMirtop(mirtop_isomirs, metadata_df)

  return(isomir)
}

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
#' @param sample_sheet nf-core sample sheet with the columns `library_id, fastq_1, fastq_2, strandedness` (see [`R/import_metadata.R`](https://github.com/cohmathonc/haemdata/blob/HEAD/R/import_metadata.R))
#' @param species mmu or hsa. If *mmu*, additionally sets the 3' adapter to TCTGGAATTCTCGGGTGCCAAGGAACTCC. For *hsa*, the 3' adapter is auto-discovered.
#' @param clip_r1_3p number of bases to trim from reads.
#' @return a path to the run script, or the multiqc report if it exists.
#' @author Denis O'Meally
#' @export

run_nfcore_smrnaseq <- function(run_folder, sample_sheet, clip_r1_3p = NULL, species = "mmu") {
  run_path <- glue::glue("{nf_core_cache}/{run_folder}")
  out_folder <- glue::glue("nfcore-smrnaseq-v{smrnaseq_release}-{species}")

  if (!dir.exists(run_path)) {
    dir.create(run_path, recursive = TRUE)
  }

  if (!is.null(clip_r1_3p)) {
    clip_r1 <- glue::glue("--three_prime_clip_r1 '{clip_r1_3p}'")
  } else {
    clip_r1 <- ""
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
      dplyr::select(sample = library_id, fastq_1) |>
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
    --publish_dir_mode link \\
    --input {run_path}/sample_sheet.csv \\
    --outdir {out_folder} --save_reference \\
    --genome {ref_genome} --igenomes_base /ref_genome/igenomes \\
    --email domeally@coh.org {three_prime_adapter} {clip_r1}
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
