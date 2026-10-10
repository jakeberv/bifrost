test_that("verbose output restores live progress rows in one console write", {
  withr::local_options(
    cli.dynamic = TRUE, cli.hide_cursor = FALSE,
    cli.num_colors = 1, cli.width = 200
  )
  withr::local_seed(75)
  seed <- .Random.seed
  renderer <- .bifrost_search_cli_renderer()
  first <- renderer$create("Candidates", 2L)
  second <- renderer$create("Searching", 2L)
  invisible(utils::capture.output({
    renderer$update(first, "complete", 2L, 2L, "Candidates done")
    renderer$update(second, "active", 1L, 2L, "Searching")
  }, type = "message"))

  # Count real console writes: a separate erase exposes a blank frame in RStudio.
  writes <- 0L
  assign("cat", function(...) {
    writes <<- writes + 1L
    base::cat(...)
  }, envir = environment(renderer$output))
  rendered <- utils::capture.output(
    renderer$output("Accepted shift"), type = "message"
  )
  text <- paste(rendered, collapse = "\n")
  plain <- cli::ansi_strip(text)
  expect_equal(writes, 1L)
  expect_match(text, "\033[1A", fixed = TRUE)
  expect_match(text, "\033[?7hAccepted shift\n\033[?7l", fixed = TRUE)
  expect_match(plain, "Accepted shift\n.*Candidates done\n.*Searching")
  expect_identical(.Random.seed, seed)

  writes <- 0L
  refreshed <- utils::capture.output(
    renderer$update(second, "active", 2L, 2L, "Search complete"),
    type = "message"
  )
  expect_equal(writes, 1L)
  expect_match(paste(refreshed, collapse = "\n"), "Search complete", fixed = TRUE)
  expect_identical(.Random.seed, seed)
  invisible(utils::capture.output(renderer$done(), type = "message"))
})
