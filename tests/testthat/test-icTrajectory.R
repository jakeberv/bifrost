# tests/testthat/test-icTrajectory.R

testthat::skip_on_cran()

make_search_for_trajectory <- function() {
  obj <- list(
    baseline_ic = -1000,
    IC_used = "GIC",
    model_fit_history = list(
      fits = list(
        list(
          step = 1L,
          candidate_node = 42L,
          regime_id = "1",
          ic = -1010,
          accepted = TRUE,
          delta_ic = 10
        ),
        list(
          step = 2L,
          candidate_node = 57L,
          regime_id = "2",
          ic = -1008,
          accepted = FALSE,
          delta_ic = -2
        ),
        list(
          step = 3L,
          candidate_node = 91L,
          regime_id = "3",
          ic = NA_real_,
          accepted = FALSE,
          status = "error",
          error = "fit failed"
        )
      )
    )
  )
  class(obj) <- c("bifrost_search", "list")
  obj
}

open_null_device_trajectory <- function() {
  grDevices::pdf(NULL)
}

close_device_trajectory <- function() {
  try(grDevices::dev.off(), silent = TRUE)
}

test_that("icTrajectory extracts baseline and proposal rows from bifrost_search", {
  traj <- icTrajectory(make_search_for_trajectory())

  testthat::expect_s3_class(traj, "icTrajectory")
  testthat::expect_equal(
    names(traj),
    c(
      "step", "ic", "accepted", "best_ic", "delta_ic",
      "status", "candidate_node", "regime_id", "error"
    )
  )
  testthat::expect_equal(traj$step, 0:3)
  testthat::expect_equal(traj$ic, c(-1000, -1010, -1008, NA))
  testthat::expect_equal(traj$accepted, c(NA, TRUE, FALSE, FALSE))
  testthat::expect_equal(traj$best_ic, c(-1000, -1010, -1010, -1010))
  testthat::expect_equal(traj$delta_ic, c(NA, 10, -2, NA))
  testthat::expect_equal(traj$status, c("baseline", "accepted", "rejected", "error"))
  testthat::expect_equal(traj$candidate_node, c(NA_integer_, 42L, 57L, 91L))
  testthat::expect_equal(traj$regime_id, c("0", "1", "2", "3"))
  testthat::expect_equal(traj$error, c(NA, NA, NA, "fit failed"))
  testthat::expect_equal(attr(traj, "IC_used"), "GIC")
})

test_that("icTrajectory computes missing delta_ic from the running best", {
  obj <- make_search_for_trajectory()
  obj$model_fit_history$fits[[1]]$delta_ic <- NULL
  obj$model_fit_history$fits[[2]]$delta_ic <- NULL

  traj <- icTrajectory(obj)

  testthat::expect_equal(traj$delta_ic[2:3], c(10, -2))
  testthat::expect_equal(traj$best_ic, c(-1000, -1010, -1010, -1010))
})

test_that("icTrajectory accepts legacy search-like lists with only an IC matrix", {
  obj <- list(
    baseline_ic = -1000,
    IC_used = "GIC",
    model_fit_history = list(
      ic_acceptance_matrix = cbind(
        ic = c(-1010, -1008, NA_real_),
        accepted = c(1L, 0L, 0L)
      )
    )
  )

  traj <- icTrajectory(obj)

  testthat::expect_s3_class(traj, "icTrajectory")
  testthat::expect_equal(traj$step, 0:3)
  testthat::expect_equal(traj$ic, c(-1000, -1010, -1008, NA))
  testthat::expect_equal(traj$accepted, c(NA, TRUE, FALSE, FALSE))
  testthat::expect_equal(traj$best_ic, c(-1000, -1010, -1010, -1010))
  testthat::expect_equal(traj$delta_ic, c(NA, 10, -2, NA))
  testthat::expect_equal(traj$status, c("baseline", "accepted", "rejected", "error"))
  testthat::expect_equal(traj$candidate_node, c(NA_integer_, NA_integer_, NA_integer_, NA_integer_))
  testthat::expect_equal(traj$regime_id, c("0", "1", "2", "3"))
  testthat::expect_equal(attr(traj, "IC_used"), "GIC")
})

test_that("icTrajectory infers regime ids for legacy fits without metadata", {
  obj <- list(
    baseline_ic = -1000,
    IC_used = "GIC",
    model_fit_history = list(
      fits = list(
        list(ic = -1010, accepted = TRUE, delta_ic = 10),
        list(ic = -1008, accepted = FALSE, delta_ic = -2)
      )
    )
  )

  traj <- icTrajectory(obj)

  testthat::expect_equal(traj$candidate_node, c(NA_integer_, NA_integer_, NA_integer_))
  testthat::expect_equal(traj$regime_id, c("0", "1", "2"))
})

test_that("icTrajectory accepts a baseline_ic override for legacy search-like lists", {
  obj <- list(
    IC_used = "GIC",
    model_fit_history = list(
      ic_acceptance_matrix = cbind(
        ic = c(-1010, -1008),
        accepted = c(1L, 0L)
      )
    )
  )

  traj <- icTrajectory(obj, baseline_ic = -995)

  testthat::expect_equal(traj$ic, c(-995, -1010, -1008))
  testthat::expect_equal(traj$best_ic, c(-995, -1010, -1010))
  testthat::expect_equal(traj$delta_ic, c(NA, 15, -2))
  testthat::expect_equal(traj$status, c("baseline", "accepted", "rejected"))
})

test_that("icTrajectory baseline_ic argument overrides stored baseline", {
  obj <- make_search_for_trajectory()

  traj <- icTrajectory(obj, baseline_ic = -995)

  testthat::expect_equal(traj$ic[1], -995)
  testthat::expect_equal(traj$best_ic[1:2], c(-995, -1010))
  testthat::expect_equal(traj$delta_ic[2], 15)
})

test_that("icTrajectory requires stored model-fit history", {
  obj <- list(baseline_ic = -1000, IC_used = "GIC")
  class(obj) <- c("bifrost_search", "list")

  testthat::expect_error(
    icTrajectory(obj),
    "`x\\$model_fit_history` is required"
  )
})

test_that("plot.icTrajectory runs with delta overlay, panel, and no delta", {
  traj <- icTrajectory(make_search_for_trajectory())

  open_null_device_trajectory()
  on.exit(close_device_trajectory(), add = TRUE)

  testthat::expect_invisible(
    plot(traj, main = "With default overlay", delta_limits = c(-5, 15))
  )
  testthat::expect_invisible(
    plot(traj, main = "With delta panel", show_delta = "panel", delta_limits = c(-5, 15))
  )
  testthat::expect_invisible(
    plot(traj, main = "With legacy overlay", show_delta = TRUE, delta_limits = c(-5, 15))
  )
  testthat::expect_invisible(
    plot(traj, main = "Without delta", show_delta = "none")
  )
  testthat::expect_invisible(
    plot(traj, main = "Without delta legacy", show_delta = FALSE)
  )
  testthat::expect_invisible(
    plot(
      traj,
      main = "Custom style",
      accepted_pch = 19,
      rejected_pch = 4,
      baseline_pch = 8,
      delta_pch = 20,
      point_scale = 1.2,
      line_scale = 1.1,
      text_scale = 0.9
    )
  )
  testthat::expect_invisible(
    plot(
      traj,
      main = "Custom legend",
      legend_position = "bottomleft",
      legend_inset = 0.02,
      legend_bty = "o",
      legend_labels = c(
        running_best = "Best",
        accepted = "Accepted",
        rejected = "Rejected",
        baseline = "Baseline",
        delta = "Delta"
      )
    )
  )
  testthat::expect_invisible(
    plot(
      traj,
      main = "Unnamed legend labels",
      legend_labels = c("Best", "Accepted", "Rejected", "Baseline", "Delta")
    )
  )
  testthat::expect_invisible(
    plot(
      traj,
      main = "Unnamed legend labels without delta",
      show_delta = "none",
      legend_labels = c("Best", "Accepted", "Rejected", "Baseline")
    )
  )
  testthat::expect_invisible(
    plot(traj, main = "No legend", legend = FALSE)
  )
})

test_that("plot.icTrajectory preserves existing layout outside panel mode", {
  traj <- icTrajectory(make_search_for_trajectory())

  open_null_device_trajectory()
  on.exit(close_device_trajectory(), add = TRUE)

  layout(matrix(1:2, nrow = 1L))
  plot(traj, show_delta = "none", legend = FALSE)

  testthat::expect_equal(par("mfg")[3:4], c(1L, 2L))
})

test_that("plot.icTrajectory validates plotting arguments", {
  traj <- icTrajectory(make_search_for_trajectory())

  open_null_device_trajectory()
  on.exit(close_device_trajectory(), add = TRUE)

  testthat::expect_error(
    plot(traj, show_delta = NA),
    "`show_delta` must be one of"
  )
  testthat::expect_error(
    plot(traj, delta_limits = c(NA_real_, 1)),
    "`delta_limits` must be a numeric vector of length 2"
  )
  testthat::expect_error(
    plot(traj, point_scale = -1),
    "`point_scale` must be a finite non-negative numeric scalar"
  )
  testthat::expect_error(
    plot(traj, accepted_pch = c(1, 2)),
    "`accepted_pch` must be a single numeric or character plotting symbol"
  )
  testthat::expect_error(
    plot(traj, accepted_pch = list(1)),
    "`accepted_pch` must be a single numeric or character plotting symbol"
  )
  testthat::expect_error(
    plot(traj, legend = NA),
    "`legend` must be TRUE or FALSE"
  )
  testthat::expect_error(
    plot(traj, legend_position = 1),
    "`legend_position` must be one of"
  )
  testthat::expect_error(
    plot(traj, legend_position = "notaposition"),
    "`legend_position` must be one of"
  )
  testthat::expect_error(
    plot(traj, legend_bty = "x"),
    "`legend_bty` must be"
  )
  testthat::expect_error(
    plot(traj, legend_labels = c("Best", accepted = "Accepted")),
    "`legend_labels` must be either fully named or unnamed"
  )
  testthat::expect_error(
    plot(traj, legend_labels = c(foo = "Foo")),
    "Named `legend_labels` entries"
  )
})

test_that("plot_ic_acceptance_matrix remains as a compatibility wrapper", {
  mat <- cbind(
    ic = c(-1000, -1010, -1008, NA_real_),
    accepted = c(1L, 1L, 0L, 0L)
  )

  open_null_device_trajectory()
  on.exit(close_device_trajectory(), add = TRUE)

  testthat::expect_invisible(
    plot_ic_acceptance_matrix(
      mat,
      plot_title = "Legacy default",
      rate_limits = c(-5, 15),
      legend = FALSE
    )
  )
  testthat::expect_invisible(
    plot_ic_acceptance_matrix(
      mat,
      plot_title = "Legacy no overlay",
      plot_rate_of_improvement = FALSE,
      baseline_ic = -995,
      accepted_pch = 19,
      point_scale = 1.2,
      legend = FALSE
    )
  )
})

test_that("plot_ic_acceptance_matrix validates legacy arguments", {
  mat <- cbind(
    ic = c(-1000, -1010, -1008),
    accepted = c(1L, 1L, 0L)
  )

  open_null_device_trajectory()
  on.exit(close_device_trajectory(), add = TRUE)

  testthat::expect_error(
    plot_ic_acceptance_matrix(mat[, 1, drop = FALSE]),
    "`matrix_data` must have at least one row and two columns"
  )
  testthat::expect_error(
    plot_ic_acceptance_matrix(mat, plot_rate_of_improvement = NA),
    "`plot_rate_of_improvement` must be TRUE or FALSE"
  )
  testthat::expect_error(
    plot_ic_acceptance_matrix(mat, rate_limits = c(NA_real_, 1)),
    "`rate_limits` must be a numeric vector of length 2"
  )
  testthat::expect_error(
    plot_ic_acceptance_matrix(cbind(ic = c(NA_real_, -1010), accepted = c(1L, 1L))),
    "first IC value in `matrix_data` must be finite"
  )
})
