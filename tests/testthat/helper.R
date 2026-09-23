# helper.R - auto-sourced by testthat before each test file.
# testthat sets the working directory to tests/testthat, so repo-root-relative paths break.
# Anchor to the repo root by walking up to the directory holding DESCRIPTION-free marker README.md
# alongside an R/ directory, so tests can source R/* regardless of where they are run from.

.find_repo_root <- function(start = getwd()) {
  d <- normalizePath(start, mustWork = FALSE)
  while (!dir.exists(file.path(d, "R")) && dirname(d) != d) d <- dirname(d)
  d
}
REPO_ROOT <- .find_repo_root()

# source a repo-root-relative R file
fn_source <- function(rel) source(file.path(REPO_ROOT, rel))
