collect_serialized_character_values <- function(x) {
  if (is.environment(x) || is.function(x) || typeof(x) == "externalptr") {
    return(character())
  }

  values <- if (is.character(x)) unname(x) else character()
  if (is.list(x) || is.pairlist(x)) {
    values <- c(
      values,
      unlist(lapply(x, collect_serialized_character_values), use.names = FALSE)
    )
  } else if (is.language(x)) {
    values <- c(values, paste(deparse(x), collapse = " "))
  }

  attrs <- attributes(x)
  if (!is.null(attrs)) {
    values <- c(
      values,
      unlist(lapply(attrs, collect_serialized_character_values), use.names = FALSE)
    )
  }
  values
}

test_that("shipped serialized artifacts do not contain absolute paths", {
  artifact_dir <- system.file("extdata", package = "bifrost")
  testthat::expect_true(nzchar(artifact_dir))
  artifact_paths <- list.files(
    artifact_dir,
    pattern = "\\.[Rr][Dd][Ss]$",
    recursive = TRUE,
    full.names = TRUE
  )
  testthat::expect_gt(length(artifact_paths), 0L)

  leaks <- unlist(lapply(artifact_paths, function(path) {
    values <- collect_serialized_character_values(readRDS(path))
    absolute <- grepl(
      "^(?:/|~[/\\\\]|[A-Za-z]:[/\\\\]|\\\\\\\\)",
      values,
      perl = TRUE
    )
    if (!any(absolute)) {
      return(character())
    }
    paste0(basename(path), ": ", unique(values[absolute]))
  }), use.names = FALSE)

  testthat::expect_equal(
    length(leaks),
    0L,
    info = paste(leaks, collapse = "\n")
  )
})
