# This is a helper script to run the pipeline.
# Choose how to execute the pipeline below.
# See https://books.ropensci.org/targets/hpc.html
# to learn about your options.


# clear all pins and upload a new version to the pinboard
targets::tar_invalidate(ends_with("_pins"))

# run the pipeline
targets::tar_make_future()

# Check the package
devtools::check(error_on = "error", vignettes = FALSE)

# Make a release on GitHub
## put fields in standard order and alphabetises dependencies
use_tidy_description()
#use_tidy_eval()
#use_version()

# build the package
tgz <- devtools::build(vignettes = FALSE)
# Publish to cgt.coh.org
drat::insertPackage(tgz, "/net/isi-dcnl/ifs/user_data/rrockne/MHO")

# install locally
#devtools::install()
install.packages("haemdata", repo = "http://cgt.coh.org/MHO")

# build the pkgdown site
pkgdown::clean_site()
devtools::document()

pkgdown::build_site()