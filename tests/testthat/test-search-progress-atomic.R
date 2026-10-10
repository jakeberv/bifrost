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

test_that("plain console output keeps verbose text before persistent progress rows", {
  withr::local_options(
    cli.dynamic = FALSE, cli.num_colors = 1, cli.width = 200
  )
  session <- .bifrost_search_progress_session(TRUE)
  invisible(utils::capture.output({
    session$skip("[1/2] Candidate scoring", "complete")
    session$skip("[2/2] Greedy search", "complete")
  }, type = "message"))

  # Exercise the real renderer: the recording renderer cannot catch dropped
  # messages or terminal control sequences leaking into non-interactive logs.
  rendered <- utils::capture.output(
    session$output("Accepted shift at node 14"), type = "message"
  )
  expect_length(rendered, 3L)
  expect_identical(rendered[[1L]], "Accepted shift at node 14")
  expect_match(rendered[[2L]], "[1/2] complete", fixed = TRUE)
  expect_match(rendered[[3L]], "[2/2] complete", fixed = TRUE)
  expect_false(any(grepl("[\033\r]", rendered)))
  session$finalize()
})
