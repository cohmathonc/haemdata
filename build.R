# This is a helper script to run the pipeline.
# Choose how to execute the pipeline below.
# See https://books.ropensci.org/targets/hpc.html
# to learn about your options.


# invalidate all pins
targets::tar_invalidate(ends_with("_pins"))

# invalidate all metadata pins
targets::tar_invalidate(matches("metadata.*_pins"))

# invalidate all mouse 10X pins
targets::tar_invalidate(matches("mmu_10x.*_pins"))

# invalidate all mouse mRNA pins
targets::tar_invalidate(matches("mmu_mrna.*_pins"))

# invalidate the package build
targets::tar_invalidate(matches("built_package"))

# run the pipeline
targets::tar_make_future()

# Check the package
devtools::check(error_on = "error", vignettes = FALSE)

# Make a release on GitHub
## put fields in standard order and alphabetises dependencies
usethis::use_tidy_description()
# use_tidy_eval()
# use_version()

# build the package
tgz <- devtools::build(vignettes = FALSE)
# Publish to cgt.coh.org
drat::insertPackage(tgz, "/net/nfs-irwrsrchnas01/labs/rrockne/MHO")

# install locally
# devtools::install()
install.packages("haemdata", repo = "http://cgt.coh.org/MHO")

# build the pkgdown site
pkgdown::clean_site()
devtools::document()
pkgdown::build_site()

# Draft a release for GitHub
usethis::use_github_release()

fs.copy

devtools::document()
pkgdown::build_reference()
pkgdown::build_home()
pkgdown::build_article(name = "missing_metadata")
pkgdown::build_article(name = "preprocessing")

