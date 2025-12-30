# tests/testthat/test-plot_ic_acceptance_matrix.R

testthat::skip_on_cran()

# ---- Helpers -----------------------------------------------------------------

make_ic_matrix <- function(n = 8L, seed = 1L, accept_every = 2L, start_ic = -1000) {
  set.seed(seed)
  # Create a roughly monotone IC series with some noise
  steps <- cumsum(sample(c(-30, -15, -10, -5, 0, 5), n - 1L, replace = TRUE))
  ic <- c(start_ic, start_ic + steps)

  # Acceptance flags: mark baseline as accepted (1) for plotting,
  # then mark every k-th step as accepted.
  acc <- integer(n)
  acc[1] <- 1L
  if (accept_every > 0L) {
    acc[seq(2L, n, by = accept_every)] <- 1L
  }

  cbind(ic, acc)
}

open_null_device <- function() {
  # Route plots to a null PDF device so CI doesn't try to open a GUI
  grDevices::pdf(NULL)
}

close_device_quietly <- function() {
  try(grDevices::dev.off(), silent = TRUE)
}

# Save/restore par() because the function tweaks margins/mgp and uses par(new=TRUE)
with_par_safely <- function(expr) {
  op <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(op), add = TRUE)
  force(expr)
}

# Group: plotting scenarios and acceptance patterns
# Test: plot_ic_acceptance_matrix runs with default settings (overlay on) (uses make_ic_matrix and null device plotting)
test_that("plot_ic_acceptance_matrix runs with default settings (overlay on)", {
  mat <- make_ic_matrix(n = 10, seed = 42, accept_every = 2)

  open_null_device()
  on.exit(close_device_quietly(), add = TRUE)

  expect_invisible(
    with_par_safely(
      plot_ic_acceptance_matrix(
        matrix_data = mat,
        plot_title = "IC Acceptance (default)",
        plot_rate_of_improvement = TRUE
      )
    )
  )
})

# Test: plot_ic_acceptance_matrix runs with overlay off without warnings (calls graphics::plot.new() to avoid par(new=TRUE) warnings)
test_that("plot_ic_acceptance_matrix runs with overlay off without warnings", {
  mat <- make_ic_matrix(n = 12, seed = 7, accept_every = 3)

  open_null_device()
  on.exit(close_device_quietly(), add = TRUE)

  # Ensure a plotting context exists so a stray par(new=TRUE) won't warn
  graphics::plot.new()

  expect_invisible(
    with_par_safely(
      plot_ic_acceptance_matrix(
        matrix_data = mat,
        plot_title = "IC Acceptance (no overlay)",
        plot_rate_of_improvement = FALSE
      )
    )
  )
})

# Test: plot_ic_acceptance_matrix tolerates many accepted steps (accept_every=1 so nearly all steps accepted)
test_that("plot_ic_acceptance_matrix tolerates many accepted steps", {
  # Nearly all steps accepted (besides baseline already marked 1)
  mat <- make_ic_matrix(n = 9, seed = 9, accept_every = 1)

  open_null_device()
  on.exit(close_device_quietly(), add = TRUE)

  expect_invisible(
    with_par_safely(
      plot_ic_acceptance_matrix(
        matrix_data = mat,
        plot_title = "IC Acceptance (many accepted)",
        plot_rate_of_improvement = TRUE
      )
    )
  )
})

# Test: plot_ic_acceptance_matrix tolerates few accepted steps (accept_every=100 with one manual accept flag)
test_that("plot_ic_acceptance_matrix tolerates few accepted steps", {
  # Only baseline and one later acceptance
  mat <- make_ic_matrix(n = 8, seed = 11, accept_every = 100)
  mat[6, 2] <- 1L  # ensure at least one accepted beyond baseline

  open_null_device()
  on.exit(close_device_quietly(), add = TRUE)

  expect_invisible(
    with_par_safely(
      plot_ic_acceptance_matrix(
        matrix_data = mat,
        plot_title = "IC Acceptance (few accepted)",
        plot_rate_of_improvement = TRUE
      )
    )
  )
})

# Group: input validation and type handling
# Test: plot_ic_acceptance_matrix errors with malformed input (non-numeric ic column triggers error)
test_that("plot_ic_acceptance_matrix errors with malformed input", {
  # Non-numeric IC column should produce an error during diff/pretty/plot usage
  bad <- cbind(ic = as.character(letters[1:5]), acc = c(1, 0, 1, 0, 1))

  open_null_device()
  on.exit(close_device_quietly(), add = TRUE)

  expect_error(
    with_par_safely(
      plot_ic_acceptance_matrix(
        matrix_data = bad,
        plot_title = "Malformed"
      )
    )
  )
})

# Test: plot_ic_acceptance_matrix accepts data.frame input as well as matrix (matrix coerced to data.frame with named columns)
test_that("plot_ic_acceptance_matrix accepts data.frame input as well as matrix", {
  mat <- make_ic_matrix(n = 10, seed = 21, accept_every = 2)
  df <- as.data.frame(mat)
  names(df) <- c("ic", "accepted")

  open_null_device()
  on.exit(close_device_quietly(), add = TRUE)

  expect_invisible(
    with_par_safely(
      plot_ic_acceptance_matrix(
        matrix_data = df,
        plot_title = "Data frame input"
      )
    )
  )
})

# Test: plot_ic_acceptance_matrix handles zero accepted steps beyond baseline (forces all accept flags to 0 after baseline)
test_that("plot_ic_acceptance_matrix handles zero accepted steps beyond baseline", {
  # Baseline is 1; force all others to 0
  mat <- make_ic_matrix(n = 8, seed = 123, accept_every = 100)
  mat[-1, 2] <- 0L  # no accepted after baseline

  open_null_device()
  on.exit(close_device_quietly(), add = TRUE)
  # Keep overlay on to exercise ROI plotting even when no accepted markers exist
  expect_invisible(
    with_par_safely(
      plot_ic_acceptance_matrix(
        matrix_data = mat,
        plot_title = "Only baseline accepted",
        plot_rate_of_improvement = TRUE
      )
    )
  )
})

# Test: plot_ic_acceptance_matrix handles logical acceptance flags (acceptance column supplied as logical)
test_that("plot_ic_acceptance_matrix handles logical acceptance flags", {
  mat <- cbind(ic = c(-1000, -990, -995), acc = c(TRUE, FALSE, TRUE))

  open_null_device()
  on.exit(close_device_quietly(), add = TRUE)
  expect_invisible(
    with_par_safely(
      plot_ic_acceptance_matrix(
        matrix_data = mat,
        plot_title = "Logical acceptance",
        plot_rate_of_improvement = TRUE
      )
    )
  )
})

# Test: plot_ic_acceptance_matrix works with minimal length (n = 2) (n=2 matrix exercises overlay on/off)
test_that("plot_ic_acceptance_matrix works with minimal length (n = 2)", {
  # Two points: baseline + one step; mark second as rejected
  mat <- cbind(ic = c(-1000, -995), acc = c(1L, 0L))

  open_null_device()
  on.exit(close_device_quietly(), add = TRUE)
  # Cover both overlay on (diff length 1) and off
  expect_invisible(
    with_par_safely(
      plot_ic_acceptance_matrix(
        matrix_data = mat,
        plot_title = "n=2 overlay on",
        plot_rate_of_improvement = TRUE
      )
    )
  )
  expect_invisible(
    with_par_safely(
      plot_ic_acceptance_matrix(
        matrix_data = mat,
        plot_title = "n=2 overlay off",
        plot_rate_of_improvement = FALSE
      )
    )
  )
})

# Test: plot_ic_acceptance_matrix validates rate_limits when overlay is on (invalid limits should error)
test_that("plot_ic_acceptance_matrix validates rate_limits when overlay is on", {
  mat <- make_ic_matrix(n = 6, seed = 99, accept_every = 2)

  open_null_device()
  on.exit(close_device_quietly(), add = TRUE)

  expect_error(
    with_par_safely(
      plot_ic_acceptance_matrix(
        matrix_data = mat,
        plot_title = "Bad rate_limits",
        plot_rate_of_improvement = TRUE,
        rate_limits = c(NA_real_, 1)
      )
    ),
    "`rate_limits` must be a numeric vector of length 2 with finite values"
  )
})
