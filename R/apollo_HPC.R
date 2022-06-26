
# path for running nf-core pipelines, slurm templates, sbatch scripts, targets cache, logs, etc...
nf_core_cache <- "/net/isi-dcnl/ifs/user_data/rrockne/MHO/haemdata-nf-core-cache"
# tar_config_set(store = glue::glue("{nf_core_cache}/_targets"))

# tar_make_clustermq() configuration:
options(clustermq.scheduler = "slurm")
options(clustermq.template = "clustermq.tmpl")

# tar_make_future() configuration:
options(future.globals.maxSize = 20 * 1024^3) # 20GB
options(future.seed = TRUE)
options(future.cache.path = glue::glue("{nf_core_cache}/_future"))

# ├ default plan -----------------------------------------------------------
# batchtools configuration
# uncomment one of:

# 1. Run individual targets on a single machine (testing)
# future::plan("multisession")

# 2. Run individual targets on cluster nodes, with default resources (recommended)
future::plan(
    list(
        future::tweak(
            future.batchtools::batchtools_slurm,
            template = "future.tmpl",
            resources = list(ncpus = 1L, memory = "8G", walltime = "1:00:00")
        ),
        future::tweak("multisession", workers = 1L)
    )
)

# ├ tweaked plans -----------------------------------------------------------
# "tweak" the plan for individual targets by choosing one of the following
# arguments to tar_target(resources = ...)
apollo_large <- targets::tar_resources(
    future = targets::tar_resources_future(plan = list(
        future::tweak(
            future.batchtools::batchtools_slurm,
            template = "future.tmpl",
            resources = list(ncpus = 20L, memory = "200G", walltime = "12:00:00")
        ),
        future::tweak("multisession", workers = 20L)
    ))
)

apollo_medium <- targets::tar_resources(
    future = targets::tar_resources_future(plan = list(
        future::tweak(
            future.batchtools::batchtools_slurm,
            template = "future.tmpl",
            resources = list(ncpus = 8L, memory = "80G", walltime = "12:00:00")
        ),
        future::tweak("multisession", workers = 8L)
    ))
)

apollo_small <- targets::tar_resources(
    future = targets::tar_resources_future(plan = list(
        future::tweak(
            future.batchtools::batchtools_slurm,
            template = "future.tmpl",
            resources = list(ncpus = 2L, memory = "20G", walltime = "2:00:00")
        ),
        future::tweak("multisession", workers = 2L)
    ))
)

apollo_bigmem <- targets::tar_resources(
    future = targets::tar_resources_future(plan = list(
        future::tweak(
            future.batchtools::batchtools_slurm,
            template = "future.tmpl",
            resources = list(ncpus = 8L, memory = "400G", walltime = "12:00:00", partition = "fast")
        ),
        future::tweak("multisession", workers = 8L)
    ))
)
