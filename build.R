# This is a helper script to run the pipeline.
# Choose how to execute the pipeline below.
# See https://books.ropensci.org/targets/hpc.html
# to learn about your options.

#targets::tar_make()
#targets::tar_make_clustermq()
targets::tar_make_future() 

# Check the package
devtools::check(error_on = "error")

# build the package
#devtools::build()

# install locally
devtools::install()

# build the pkgdown site
pkgdown::clean_site()
devtools::document()

pkgdown::build_site()
