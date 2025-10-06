# Build Documentation
# This script rebuilds the pkgdown site and quarto notebooks
# Run this when you want to update documentation without a full release
#
# Usage:
#   Rscript build_docs.R              # Build and auto-commit/push
#   Rscript build_docs.R --no-push    # Build only, manual commit/push

# TODO: Add host check to ensure running on apollo with NFS access
# if (!check_host()) {
#   cli::cli_abort("This script must be run on apollo (needs NFS access for data/multiqc)")
# }

# Parse command line args
args <- commandArgs(trailingOnly = TRUE)
auto_push <- !("--no-push" %in% args)

# Clean and rebuild pkgdown site
message("Building pkgdown site...")
pkgdown::clean_site()
devtools::document()
pkgdown::build_site()

# Render quarto notebooks (uses freeze: auto for caching - only rebuilds changed notebooks)
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

# Auto-commit and push unless --no-push flag
if (auto_push) {
  message("Committing and pushing changes...")
  system("git add docs/")
  system('git commit -m "Update documentation

🤖 Generated with build_docs.R"')
  system("git push")
  message("Pushed to GitHub. Pages will auto-deploy shortly.")
} else {
  message("Skipping auto-commit/push (--no-push flag set)")
  message("Manual steps:")
  message("  1. Review changes: git status")
  message("  2. Commit: git add docs/ && git commit -m 'Update documentation'")
  message("  3. Push: git push")
  message("  4. GitHub Pages will auto-deploy")
}
