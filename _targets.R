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
pipeline_version <- "3.7"

list(
    # make the package logo
    tar_target(logo, make_logo()
    ),
    # get raw data
    #CML
    tar_target(cml_mrna_2021_m38, get_rnaseq_se("/net/isi-dcnl/ifs/user_data/rrockne/MHO/CML.mRNA.2021/results/rnaseq3.5.GRCm38.HLT")),
    tar_target(cml_mrna_2022_m38, get_rnaseq_se("/net/isi-dcnl/ifs/user_data/rrockne/MHO/CML.mRNA.2022/results/rnaseq-v3.5_GRCm38.HLT")),
    tar_target(cml_mrna_m38, merge_mrna(cml_mrna_2021_m38, cml_mrna_2022_m38, drop = c("F545.7", "X490.17"))),

    # Build sample sheets
    tar_target(sample_sheet_2016_1, make_sample_sheet("AML.mRNA.2016")),
    tar_target(sample_sheet_2018_1, make_sample_sheet("AML.mRNA.2018.all_samples")),
    tar_target(sample_sheet_2020_1, make_sample_sheet("AML.mRNA.2020")),
    tar_target(sample_sheet_2021_1, make_sample_sheet("AML.mRNA.2021.RxGroup1")),
    tar_target(sample_sheet_2021_2, make_sample_sheet("AML.mRNA.2021.RxGroup2")),
    tar_target(sample_sheet_2021_3, make_sample_sheet("AML.mRNA.2021.RxGroup2_pt2")),
    tar_target(sample_sheet_2021_4, make_sample_sheet("AML.mRNA.2022.RxGroup3")),
    tar_target(sample_sheet_2017_1, make_sample_sheet("AML.validation.2017")),
    tar_target(sample_sheet_2020_2, make_sample_sheet("AML.mRNA.novaseq_validation.2020")),
    tar_target(sample_sheet_CML_1, make_sample_sheet("CML.mRNA.2021")),
    tar_target(sample_sheet_CML_2, make_sample_sheet("CML.mRNA.2022")),
    tar_target(sample_sheet_2022_1, make_sample_sheet("AML.scRNAseq.2022")),
    tar_target(sample_sheet_2022_2, make_sample_sheet("AML.mRNA.HSA_FLT3.2022")),

    # combine sample sheets for running

    tar_target(sample_sheet_2016_2018, rbind(sample_sheet_2016_1, sample_sheet_2018_1)),
    tar_target(sample_sheet_CML, rbind(sample_sheet_CML_1, sample_sheet_CML_2)),
    tar_target(sample_sheet_all, rbind(sample_sheet_2016_1,
                                        sample_sheet_2018_1,
                                        sample_sheet_2020_1,
                                        sample_sheet_2021_1,
                                        sample_sheet_2021_2,
                                        sample_sheet_2021_3,
                                        sample_sheet_2021_4,
                                        sample_sheet_2017_1,
                                        sample_sheet_2020_2,
                                        sample_sheet_CML_1,
                                        sample_sheet_CML_2,
                                        sample_sheet_2022_1))

    
    tar_target(aml_mrna_2016_m38, run_rnaseq(run_folder = "all", sample_sheet = sample_sheet_2022_2, reference_genome = "GENCODEr33")),

)
