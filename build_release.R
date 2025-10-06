# Build Release
# This script creates a new package release on GitHub
# Run this when you're ready to release a new version

# TODO: Add host check to ensure running on apollo with NFS access
# if (!check_host()) {
#   cli::cli_abort("This script must be run on apollo (needs NFS access for data/multiqc)")
# }

# Check the package
message("Running package checks...")
devtools::check(error_on = "error", vignettes = FALSE)

# Standardize DESCRIPTION
message("Tidying DESCRIPTION...")
usethis::use_tidy_description()

# Bump version (interactive - will prompt for major/minor/patch)
message("Bumping version...")
# Uncomment when ready to bump version:
# usethis::use_version()

# Build package tarball (for validation and GitHub release)
message("Building package...")
tgz <- devtools::build(vignettes = FALSE)
message("Built: ", tgz)

# Build documentation
message("Building documentation...")
source(here::here("build_docs.R"))

# Create GitHub release (interactive)
message("Creating GitHub release...")
message("Next steps:")
message("  1. Review all changes")
message("  2. Commit and push all changes")
message("  3. Run: usethis::use_github_release()")
message("  4. GitHub will create release with tarball")
message("  5. Users install via: remotes::install_github('cohmathonc/haemdata@vX.X.X')")

# Manual steps for now:
# 1. git add .
# 2. git commit -m "Release vX.X.X"
# 3. git push
# 4. usethis::use_github_release()
