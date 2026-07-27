.avian_skeleton_artifact_path <- function(name, label) {
  root <- Sys.getenv("BIFROST_ARTIFACT_DIR", unset = "")
  testthat::skip_if(
    !nzchar(root) || !dir.exists(root),
    paste("repository empirical artifact directory unavailable for", label)
  )
  bifrost_example_file(name)
}

.avian_skeleton_read_compact <- function(name, label) {
  path <- .avian_skeleton_artifact_path(name, label)
  readRDS(path)
}

.avian_skeleton_restore_search_class <- function(search) {
  if (!inherits(search, "bifrost_search")) {
    class(search) <- c("bifrost_search", class(search))
  }
  search
}

.avian_skeleton_compact_search <- function() {
  search <- .avian_skeleton_read_compact(
    "passerine-search",
    "focal search"
  )
  .avian_skeleton_restore_search_class(search)
}

.avian_skeleton_compact_sensitivity <- function() {
  .avian_skeleton_read_compact(
    "passerine-sensitivity",
    "sensitivity bundle"
  )
}

.avian_skeleton_compact_posthoc <- function() {
  .avian_skeleton_read_compact(
    "passerine-posthoc",
    "post-hoc object"
  )
}
