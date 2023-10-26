# Configurations for Apollo {crew} & {crew.cluster} controllers

# path for running nf-core pipelines, slurm templates, sbatch scripts, targets cache, logs, etc...
cache_path <- nf_core_cache
crew_cache_path <- here::here(glue::glue("{cache_path}/_crew"))
log_slurm <- TRUE

# Setup
library("crew")
library("crew.cluster")
options(future.seed = TRUE)
options(future.cache.path = glue::glue("{crew_cache_path}/_future"))
targets::tar_config_set("workers" = 350)

# #-----------------------------------------------------------------------------
# tar_option_set(
#     workspace_on_error = TRUE,
#     format = "qs",
#     storage = "worker",
#     retrieval = "worker",
#     controller = crew::crew_controller_group(local, small, medium, large, bigmem),
#     resources = tar_resources(
#     # Run jobs locally by default
#      crew = tar_resources_crew(controller = "local")
#     )
# )


#-----------------------------------------------------------------------------
# Crew controllers
local <- crew::crew_controller_local(
    name = "local",
    workers = 4L,
    seconds_interval = 0.25,
    seconds_timeout = 10,
    seconds_launch = 30,
    seconds_idle = 30,
    tasks_max = Inf,
    tasks_timers = 0L,
)

small <- crew.cluster::crew_controller_slurm(
    name = "small",
    slurm_memory_gigabytes_per_cpu = 20L,
    slurm_cpus_per_task = 2L,
    slurm_time_minutes = 360L,
    host = Sys.info()["nodename"],
    workers = 2L,
    seconds_idle = 30,
    script_directory = here::here(glue::glue("{crew_cache_path}")),
    slurm_log_output = if (isTRUE(log_slurm)) here::here(glue::glue("{crew_cache_path}/slurm-%j.out")) else NULL,
    slurm_log_error  = if (isTRUE(log_slurm)) here::here(glue::glue("{crew_cache_path}/slurm-%j.err")) else NULL,
    script_lines = glue::glue("#SBATCH --partition=fast,all \
cd {here::here()} \
/opt/singularity/3.7.0/bin/singularity exec \\
    --env R_LIBS_USER=~/R/bioc-3.17 \\
    --env R_LIBS_SITE=/opt/singularity-images/rbioc/rlibs/bioc-3.17 \\
    -B /labs,/opt/singularity,/opt/singularity-images \\
    /opt/singularity-images/rbioc/vscode-rbioc_3.17.sif \\"),
    verbose = TRUE
)

medium <- crew.cluster::crew_controller_slurm(
    name = "medium",
    slurm_memory_gigabytes_per_cpu = 80L,
    slurm_cpus_per_task = 8L,
    slurm_time_minutes = 720L,
    host = Sys.info()["nodename"],
    workers = 8L,
    seconds_idle = 30,
    script_directory = here::here(glue::glue("{crew_cache_path}")),
    slurm_log_output = if (isTRUE(log_slurm)) here::here(glue::glue("{crew_cache_path}/slurm-%j.out")) else NULL,
    slurm_log_error  = if (isTRUE(log_slurm)) here::here(glue::glue("{crew_cache_path}/slurm-%j.err")) else NULL,
    script_lines = glue::glue("#SBATCH --partition=fast,all \
cd {here::here()} \
/opt/singularity/3.7.0/bin/singularity exec \\
    --env R_LIBS_USER=~/R/bioc-3.17 \\
    --env R_LIBS_SITE=/opt/singularity-images/rbioc/rlibs/bioc-3.17 \\
    -B /labs,/opt/singularity,/opt/singularity-images \\
    /opt/singularity-images/rbioc/vscode-rbioc_3.17.sif \\"),
    verbose = TRUE
)

large <- crew.cluster::crew_controller_slurm(
    name = "large",
    slurm_memory_gigabytes_per_cpu = 200L,
    slurm_cpus_per_task = 20L,
    slurm_time_minutes = 720L,
    host = Sys.info()["nodename"],
    workers = 20L,
    seconds_idle = 30,
    script_directory = here::here(glue::glue("{crew_cache_path}")),
    slurm_log_output = if (isTRUE(log_slurm)) here::here(glue::glue("{crew_cache_path}/slurm-%j.out")) else NULL,
    slurm_log_error  = if (isTRUE(log_slurm)) here::here(glue::glue("{crew_cache_path}/slurm-%j.err")) else NULL,
    script_lines = glue::glue("#SBATCH --partition=fast,all \
cd {here::here()} \
/opt/singularity/3.7.0/bin/singularity exec \\
    --env R_LIBS_USER=~/R/bioc-3.17 \\
    --env R_LIBS_SITE=/opt/singularity-images/rbioc/rlibs/bioc-3.17 \\
    -B /labs,/opt/singularity,/opt/singularity-images \\
    /opt/singularity-images/rbioc/vscode-rbioc_3.17.sif \\"),
    verbose = TRUE
)

bigmem <- crew.cluster::crew_controller_slurm(
    name = "bigmem",
    slurm_memory_gigabytes_per_cpu = 200L,
    slurm_cpus_per_task = 6L,
    slurm_time_minutes = 360L,
    host = Sys.info()["nodename"],
    workers = 6L,
    seconds_idle = 30,
    script_directory = here::here(glue::glue("{crew_cache_path}")),
    slurm_log_output = if (isTRUE(log_slurm)) here::here(glue::glue("{crew_cache_path}/slurm-%j.out")) else NULL,
    slurm_log_error  = if (isTRUE(log_slurm)) here::here(glue::glue("{crew_cache_path}/slurm-%j.err")) else NULL,
    script_lines = glue::glue("#SBATCH --partition=all \
cd {here::here()} \
/opt/singularity/3.7.0/bin/singularity exec \\
    --env R_LIBS_USER=~/R/bioc-3.17 \\
    --env R_LIBS_SITE=/opt/singularity-images/rbioc/rlibs/bioc-3.17 \\
    -B /labs,/opt/singularity,/opt/singularity-images \\
    /opt/singularity-images/rbioc/vscode-rbioc_3.17.sif \\"),
    verbose = TRUE
)

## Some shortcuts

apollo_small <- tar_resources(
      crew = tar_resources_crew(controller = "small")
    )

apollo_medium <- tar_resources(
      crew = tar_resources_crew(controller = "medium")
    )

apollo_large <- tar_resources(
      crew = tar_resources_crew(controller = "large")
    )

apollo_bigmem <- tar_resources(
      crew = tar_resources_crew(controller = "bigmem")
    )
