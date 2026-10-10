test_that("bar updates do not change the spinner's elapsed-time frames or RNG", {
  withr::local_options(cli.dynamic = FALSE, cli.num_colors = 1, cli.width = 200)
  withr::local_seed(73)
  seed <- .Random.seed
  testthat::local_mocked_bindings(
    get_spinner = function(...) list(frames = LETTERS[1:8], interval = 80),
    .package = "cli"
  )
  # The final two observations follow a pause: animation must catch up to the
  # elapsed-time frame, not count redraws or restart its clock at a bar update.
  heartbeats <- c(seq(0, .5, by = .1), 1.30, 1.31)
  run_timeline <- function(bar_updates) {
    renderer <- .bifrost_search_cli_renderer()
    time <- 0
    assign("now", function() time, envir = environment(renderer$create))
    row <- renderer$create("Fitting", 10L)
    current <- 0L
    frames <- character()
    for (t in sort(c(heartbeats, bar_updates))) {
      time <- t
      is_bar_update <- t %in% bar_updates
      if (is_bar_update) current <- current + 1L
      line <- utils::capture.output(
        renderer$update(row, "active", current, 10L, "Fitting"), type = "message"
      )
      if (!is_bar_update) frames <- c(frames, substr(cli::ansi_strip(line), 1L, 1L))
    }
    frames
  }
  plain <- run_timeline(numeric())
  with_bar_updates <- run_timeline(c(.081, .162, .243, .324))
  expected <- c("A", "B", "C", "D", "F", "G", "A", "A")
  expect_identical(plain, expected)
  expect_identical(with_bar_updates, expected)
  expect_identical(with_bar_updates, plain)
  expect_identical(.Random.seed, seed)
})

test_that("default polling uses 100 ms without changing results or RNG", {
  withr::local_seed(74)
  seed <- .Random.seed
  time <- 0
  ticks <- numeric()
  testthat::local_mocked_bindings(
    resolved = function(...) time >= .29,
    value = function(...) 42L,
    .package = "future"
  )
  testthat::local_mocked_bindings(
    Sys.sleep = function(seconds) time <<- time + seconds,
    .package = "base"
  )
  result <- .bifrost_search_await_futures(
    list("pending fit"), heartbeat = function() ticks <<- c(ticks, time)
  )
  expect_identical(result, list(42L))
  expect_length(ticks, 3L)
  expect_equal(ticks, c(0, .1, .2))
  expect_identical(.Random.seed, seed)
})
