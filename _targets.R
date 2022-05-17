# Created by use_targets().
# Follow the comments below to fill in this target script.
# Then follow the manual to check and run the pipeline:
#   https://books.ropensci.org/targets/walkthrough.html#inspect-the-pipeline # nolint

# Load packages required to define the pipeline:
library(targets)
library(tarchetypes) # Load other packages as needed. # nolint

# Set target options:
tar_option_set(
  packages = c("tidyverse"), # packages that your targets need to run
  format = "rds" # default storage format
  # Set other options as needed.
)

# tar_make_clustermq() configuration (okay to leave alone):
options(clustermq.scheduler = "slurm")
options(clustermq.template = "clustermq.tmpl")

# tar_make_future() configuration (okay to leave alone):
future::plan(future.batchtools::batchtools_slurm, template = "future.tmpl")

# Load the R scripts with your custom functions:
for (file in list.files("R", full.names = TRUE)) source(file)

# Global variables:

rnaseq_version <- "3.5"

# Replace the target list below with your own:
list(
    # make the package logo
    tar_target(
        logo,
        make_logo("https://www.maxpixel.net/static/photo/2x/Blood-Group-0-Blood-Rh-factor-Positive-2781421.jpg")
    ),
    # get raw data
    #CML
    tar_target(cml_mrna_2021_m38, get_rnaseq_se("/net/isi-dcnl/ifs/user_data/rrockne/MHO/CML.mRNA.2021/results/rnaseq3.5.GRCm38.HLT")),
    tar_target(cml_mrna_2022_m38, get_rnaseq_se("/net/isi-dcnl/ifs/user_data/rrockne/MHO/CML.mRNA.2022/results/rnaseq-v3.5_GRCm38.HLT")),
    tar_target(cml_mrna_m38, merge_mrna(cml_mrna_2021_m38, cml_mrna_2022_m38, drop = c("F545.7", "X490.17")))


    #AML
    #tar_target(aml_mrna_2016_m38, run_rnaseq_salmon(project, sample_sheet, reference_genome, pipeline_version, )),

)
