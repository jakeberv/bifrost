test_that("new vignette figure captions are format-aware", {
  expected_figures <- c(
    "avian-skeleton-part-1" = 3L,
    "avian-skeleton-part-2" = 2L,
    "avian-skeleton-part-3" = 2L,
    "avian-skeleton-part-4" = 2L,
    "avian-skeleton-part-5" = 3L,
    "simulation-study-part-1" = 1L,
    "simulation-study-part-2" = 0L
  )

  for (slug in names(expected_figures)) {
    vignette_path <- testthat::test_path(
      "../../vignettes",
      paste0(slug, ".Rmd")
    )
    testthat::skip_if_not(
      file.exists(vignette_path),
      paste("Vignette not present:", slug)
    )

    vignette_source <- readLines(vignette_path, warn = FALSE)
    caption_lines <- grep(
      "fig[.]cap[[:space:]]*=",
      vignette_source,
      value = TRUE
    )

    testthat::expect_length(
      caption_lines,
      expected_figures[[slug]]
    )

    if (length(caption_lines) > 0L) {
      testthat::expect_true(
        any(grepl(
          "vignette_figure_caption <- function",
          vignette_source,
          fixed = TRUE
        )),
        info = slug
      )
      testthat::expect_true(
        all(grepl(
          "fig[.]cap[[:space:]]*=[[:space:]]*vignette_figure_caption[(]",
          caption_lines
        )),
        info = slug
      )
    }
  }
})
