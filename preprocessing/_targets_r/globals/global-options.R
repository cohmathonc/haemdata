options(tidyverse.quiet = TRUE)
tar_option_set(
    packages = c("tidyverse", "SummarizedExperiment"), # packages that targets need to run
    error = "abridge", # continue or stop on error
    # format = "qs", # default storage format
    storage = "worker",
    retrieval = "worker"
    # garbage_collection = TRUE,
    # memory = "transient"
)
# Load the R scripts & functions:
for (file in list.files("scripts", full.names = TRUE)) source(file)
for (file in list.files("R", pattern = "*.R", full.names = TRUE)) source(file)
published_metadata_mmu <- readRDS("data-raw/metadata_mmu.rds")
