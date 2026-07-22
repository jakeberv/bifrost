test_that("lineage-rate HTML controls retain their authored DOM", {
  source <- testthat::test_path(
    "../../vignettes/avian-skeleton-part-2.Rmd"
  )
  testthat::skip_if_not(
    file.exists(source),
    "Lineage-rate vignette source is unavailable in installed tests"
  )
  testthat::skip_if_not_installed("rmarkdown")
  testthat::skip_if_not(
    rmarkdown::pandoc_available(),
    "Pandoc is required for the HTML widget regression test"
  )

  output_dir <- tempfile("lineage-rate-widget-")
  dir.create(output_dir)
  rendered <- rmarkdown::render(
    input = source,
    output_format = "rmarkdown::html_vignette",
    output_file = "lineage-rate-widget.html",
    output_dir = output_dir,
    envir = new.env(parent = globalenv()),
    quiet = TRUE,
    clean = TRUE
  )
  html <- paste(readLines(rendered, warn = FALSE), collapse = "\n")

  rate_match <- regexpr(
    '(?s)<div class="ldw-rate-controls"[^>]*>.*?</div>',
    html,
    perl = TRUE
  )
  testthat::expect_gt(rate_match[[1L]], 0L)
  rate_controls <- regmatches(html, rate_match)
  testthat::expect_match(
    rate_controls,
    '^<div class="ldw-rate-controls"[^>]*>\\s*<label',
    perl = TRUE
  )
  testthat::expect_false(grepl("<p", rate_controls, fixed = TRUE))
  testthat::expect_equal(
    length(strsplit(
      rate_controls,
      'class="ldw-rate-control"',
      fixed = TRUE
    )[[1L]]) - 1L,
    3L
  )

  axis_match <- regexpr(
    '(?s)<div class="ldw-axis-slider"[^>]*>.*?</div>',
    html,
    perl = TRUE
  )
  testthat::expect_gt(axis_match[[1L]], 0L)
  axis_slider <- regmatches(html, axis_match)
  testthat::expect_false(grepl("<p", axis_slider, fixed = TRUE))
  testthat::expect_match(
    axis_slider,
    '(?s)</label>\\s*<input[^>]+type="range"',
    perl = TRUE
  )

  readout_ids <- gregexpr(
    'id="ldw-readout-part2"',
    html,
    fixed = TRUE
  )[[1L]]
  testthat::expect_equal(sum(readout_ids > 0L), 1L)
  testthat::expect_false(grepl("<p>Myr", html, fixed = TRUE))
})
