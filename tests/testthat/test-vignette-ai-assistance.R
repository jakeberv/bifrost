test_that("public vignettes end with the canonical AI assistance disclosure", {
  vignette_dir <- testthat::test_path("../../vignettes")
  paths <- sort(list.files(
    vignette_dir,
    pattern = "\\.Rmd$",
    full.names = TRUE,
    recursive = FALSE
  ))
  testthat::skip_if(
    length(paths) == 0L,
    "Vignette sources are unavailable in installed tests"
  )

  heading <- "## AI Assistance"
  disclosure <- paste(
    "This vignette was developed with assistance from OpenAI tools for",
    "drafting, editing, and figure refinement; all scientific content,",
    "interpretation, and final decisions were reviewed by the authors."
  )

  invalid <- vapply(paths, function(path) {
    lines <- readLines(path, warn = FALSE)
    text <- paste(lines, collapse = "\n")
    nonempty <- lines[nzchar(trimws(lines))]

    heading_count <- sum(trimws(lines) == heading)
    disclosure_count <- lengths(regmatches(
      text,
      gregexpr(disclosure, text, fixed = TRUE)
    ))
    footer <- tail(nonempty, 2L)

    heading_count != 1L ||
      disclosure_count != 1L ||
      !identical(footer, c(heading, disclosure))
  }, logical(1L))

  testthat::expect_equal(basename(paths[invalid]), character())
})
