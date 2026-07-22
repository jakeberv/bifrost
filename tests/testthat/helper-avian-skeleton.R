.avian_skeleton_extdata_path <- function(filename) {
  installed_path <- system.file(
    "extdata",
    "avian-skeleton",
    filename,
    package = "bifrost"
  )
  if (nzchar(installed_path)) {
    return(installed_path)
  }

  local_path <- file.path(
    getwd(),
    "inst",
    "extdata",
    "avian-skeleton",
    filename
  )
  if (file.exists(local_path)) {
    return(local_path)
  }

  ""
}

.avian_skeleton_read_compact <- function(filename, label) {
  path <- .avian_skeleton_extdata_path(filename)
  testthat::skip_if(!nzchar(path), paste("compact avian skeleton", label, "unavailable"))
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
    "passerine_bodyplan_search_compact.RDS",
    "focal search"
  )
  .avian_skeleton_restore_search_class(search)
}

.avian_skeleton_compact_sensitivity <- function() {
  .avian_skeleton_read_compact(
    "passerine_bodyplan_search_sensitivity_compact.RDS",
    "sensitivity bundle"
  )
}

.avian_skeleton_compact_posthoc <- function() {
  .avian_skeleton_read_compact(
    "passerine_bodyplan_posthoc_integration_compact.RDS",
    "post-hoc object"
  )
}
