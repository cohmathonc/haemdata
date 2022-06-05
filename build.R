# This is a helper script to run the pipeline.
# Choose how to execute the pipeline below.
# See https://books.ropensci.org/targets/hpc.html
# to learn about your options.

targets::tar_make()
# targets::tar_make_clustermq(workers = 4) # nolint
#targets::tar_make_future(workers = 5) # nolint

# Check the package
devtools::check(error_on = "error")

# install locally
devtools::install()

# build the pkgdown site
pkgdown::build_site()
