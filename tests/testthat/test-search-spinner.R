test_that("spinner advances one frame after delays and respects its frame interval", {
  withr::local_options(cli.dynamic = FALSE, cli.num_colors = 1, cli.width = 200)
  withr::local_seed(73)
  seed <- .Random.seed
  testthat::local_mocked_bindings(
    get_spinner = function(...) list(frames = c("A", "B", "C"), interval = 80),
    .package = "cli"
  )
  renderer <- .bifrost_search_cli_renderer()
  time <- 0
  assign("now", function() time, envir = environment(renderer$create))
  row <- renderer$create("Fitting", 2L)
  frames <- vapply(c(0, .04, .40, .80, 1.20, 1.21), function(t) {
    time <<- t
    line <- utils::capture.output(
      renderer$update(row, "active", 0L, 2L, "Fitting"), type = "message"
    )
    substr(cli::ansi_strip(line), 1L, 1L)
  }, character(1))
  expect_identical(frames, c("A", "A", "B", "C", "A", "A"))
  expect_identical(.Random.seed, seed)
})

test_that("default polling keeps progress responsive without changing results or RNG", {
  withr::local_seed(74)
  seed <- .Random.seed
  time <- 0
  ticks <- numeric()
  testthat::local_mocked_bindings(
    resolved = function(...) time >= .19,
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
  expect_length(ticks, 4L)
  expect_equal(ticks, c(0, .05, .10, .15))
  expect_identical(.Random.seed, seed)
})
