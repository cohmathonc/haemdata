# Build Documentation
# This script rebuilds the pkgdown site and quarto notebooks
# Run this when you want to update documentation without a full release

# TODO: Add host check to ensure running on apollo with NFS access
# if (!check_host()) {
#   cli::cli_abort("This script must be run on apollo (needs NFS access for data/multiqc)")
# }

# Clean and rebuild pkgdown site
message("Building pkgdown site...")
pkgdown::clean_site()
devtools::document()
pkgdown::build_site()

# Render quarto notebooks
message("Rendering quarto notebooks...")
quarto::quarto_render(input = here::here("notebooks/"), as_job = FALSE)

# Copy rendered notebooks to docs/
message("Copying notebooks to docs/notebooks/...")
fs::dir_copy(
  path = here::here("notebooks/public"),
  new_path = here::here("docs/notebooks"),
  overwrite = TRUE
)

message("Documentation build complete!")
message("Next steps:")
message("  1. Review changes: git status")
message("  2. Commit: git add docs/ && git commit -m 'Update documentation'")
message("  3. Push: git push")
message("  4. GitHub Pages will auto-deploy")
