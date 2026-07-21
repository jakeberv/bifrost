# tests/testthat/test-icTrajectory.R

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

trajectory_test_internals <- c(
  "icTrajectory.default",
  "plot.icTrajectory",
  ".bifrost_search_ic_value",
  ".icTrajectory_axis_limits",
  ".icTrajectory_delta_limits",
  ".icTrajectory_entry_accepted",
  ".icTrajectory_entry_candidate_node",
  ".icTrajectory_entry_ic",
  ".icTrajectory_entry_regime_id",
  ".icTrajectory_entry_step",
  ".icTrajectory_legend_group",
  ".icTrajectory_legend_inset",
  ".icTrajectory_legend_label_vector",
  ".icTrajectory_legend_labels",
  ".icTrajectory_legend_labels_arg",
  ".icTrajectory_logical_vector",
  ".icTrajectory_named_group",
  ".icTrajectory_numeric_vector",
  ".icTrajectory_plot_delta_panel",
  ".icTrajectory_plot_style",
  ".icTrajectory_scalar_numeric",
  ".icTrajectory_show_delta_mode",
  ".icTrajectory_y_label"
)
trajectory_test_ns <- asNamespace("bifrost")
for (name in trajectory_test_internals) {
  value <- if (exists(name, mode = "function", inherits = TRUE)) {
    get(name, mode = "function", inherits = TRUE)
  } else if (identical(name, ".bifrost_search_ic_value")) {
    function(model_result, IC) {
      if (identical(IC, "GIC")) {
        model_result$GIC$GIC
      } else {
        model_result$BIC$BIC
      }
    }
  } else {
    get(name, envir = trajectory_test_ns, mode = "function", inherits = FALSE)
  }
  assign(name, value, envir = environment())
}

testthat::test_that("icTrajectory extracts baseline and proposal rows from bifrost_search", {
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

testthat::test_that("icTrajectory computes missing delta_ic from the running best", {
  obj <- make_search_for_trajectory()
  obj$model_fit_history$fits[[1]]$delta_ic <- NULL
  obj$model_fit_history$fits[[2]]$delta_ic <- NULL

  traj <- icTrajectory(obj)

  testthat::expect_equal(traj$delta_ic[2:3], c(10, -2))
  testthat::expect_equal(traj$best_ic, c(-1000, -1010, -1010, -1010))
})

testthat::test_that("icTrajectory accepts legacy search-like lists with only an IC matrix", {
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

testthat::test_that("icTrajectory supports skeleton-style hybrid history", {
  obj <- list(
    baseline_ic = -1000,
    IC_used = "GIC",
    model_fit_history = list(
      fits = list(
        list(model = list(GIC = list(GIC = -9999)), accepted = TRUE, delta_ic = 10),
        list(model = list(GIC = list(GIC = -9998)), delta_ic = -2)
      ),
      ic_acceptance_matrix = cbind(
        ic = c(-1010, -1008),
        accepted = c(1L, 0L)
      )
    )
  )

  traj <- icTrajectory(obj)

  testthat::expect_equal(traj$ic, c(-1000, -1010, -1008))
  testthat::expect_equal(traj$accepted, c(NA, TRUE, FALSE))
  testthat::expect_equal(traj$candidate_node, c(NA_integer_, NA_integer_, NA_integer_))
  testthat::expect_equal(traj$regime_id, c("0", "1", "2"))
})

testthat::test_that("icTrajectory accepts a baseline_ic override for legacy search-like lists", {
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

testthat::test_that("icTrajectory baseline_ic argument overrides stored baseline", {
  obj <- make_search_for_trajectory()

  traj <- icTrajectory(obj, baseline_ic = -995)

  testthat::expect_equal(traj$ic[1], -995)
  testthat::expect_equal(traj$best_ic[1:2], c(-995, -1010))
  testthat::expect_equal(traj$delta_ic[2], 15)
})

testthat::test_that("icTrajectory requires stored model-fit history", {
  obj <- list(baseline_ic = -1000, IC_used = "GIC")
  class(obj) <- c("bifrost_search", "list")

  testthat::expect_error(
    icTrajectory(obj),
    "`x\\$model_fit_history` is required"
  )
})

testthat::test_that("icTrajectory handles extractor edge cases", {
  traj <- icTrajectory(make_search_for_trajectory())
  testthat::expect_identical(icTrajectory.default(traj), traj)
  testthat::expect_error(
    icTrajectory.default(list(
      baseline_ic = -1,
      model_fit_history = list(fits = list())
    )),
    "`x` must be a `bifrost_search` object"
  )
  testthat::expect_error(
    icTrajectory.default(1),
    "`x` must be a `bifrost_search` object"
  )
  testthat::expect_error(
    icTrajectory(list(
      baseline_ic = NA_real_,
      model_fit_history = list(fits = list())
    )),
    "`baseline_ic` must be a finite numeric scalar"
  )
  testthat::expect_error(
    icTrajectory(list(
      baseline_ic = -1,
      model_fit_history = list()
    )),
    "`x\\$model_fit_history\\$fits` or"
  )
  testthat::expect_error(
    icTrajectory(list(
      baseline_ic = -1,
      model_fit_history = list(ic_acceptance_matrix = matrix(numeric(), nrow = 0, ncol = 2))
    )),
    "`x\\$model_fit_history\\$ic_acceptance_matrix` must have at least two columns"
  )

  bad_entry <- icTrajectory(list(
    baseline_ic = -1,
    model_fit_history = list(fits = list("not an entry"))
  ))
  testthat::expect_equal(bad_entry$status, c("baseline", "error"))

  model_entry <- .icTrajectory_entry_ic(
    list(model = list(GIC = list(GIC = -123))),
    "GIC"
  )
  testthat::expect_equal(
    list(
      model_ic = model_entry,
      bad_step = .icTrajectory_entry_step(list(step = "nope"), 7),
      no_acceptance = .icTrajectory_entry_accepted(list()),
      numeric_acceptance = .icTrajectory_entry_accepted(list(accepted = 1)),
      bad_acceptance = .icTrajectory_entry_accepted(list(accepted = "maybe")),
      bad_candidate_node = .icTrajectory_entry_candidate_node(list(candidate_node = "nope")),
      blank_regime_id = .icTrajectory_entry_regime_id(list(regime_id = ""))
    ),
    list(
      model_ic = -123,
      bad_step = 7L,
      no_acceptance = NA,
      numeric_acceptance = TRUE,
      bad_acceptance = NA,
      bad_candidate_node = NA_integer_,
      blank_regime_id = NA_character_
    )
  )
})

testthat::test_that("plot.icTrajectory runs core display modes", {
  traj <- icTrajectory(make_search_for_trajectory())

  open_null_device_trajectory()
  on.exit(close_device_trajectory(), add = TRUE)

  testthat::expect_invisible(
    plot(traj, main = "With default overlay", delta_limits = c(-5, 15))
  )
  testthat::expect_invisible(
    plot(
      traj,
      main = "Grouped reversed-limit overlay",
      ic_limits = c(-1020, -990),
      delta_limits = c(15, -5),
      symbols = c(
        accepted = 19,
        rejected = 4,
        baseline = 8,
        delta = 20
      ),
      scales = c(point = 1.2, line = 1.1, text = 0.9),
      point_sizes = c(accepted = 0.8, rejected = 0.5),
      line_widths = c(running_best = 1.4, rejected = 1.2),
      text_sizes = c(
        x_axis = 0.8,
        y_axis = 0.7,
        delta_axis = 0.65,
        axis_label = 0.85
      ),
      annotation = c(baseline_label_offset = 0.8),
      legend = list(
        position = "bottomleft",
        inset = 0.02,
        bty = "o",
        labels = c(
          running_best = "Best",
          accepted = "Accepted",
          rejected = "Rejected",
          baseline = "Baseline",
          delta = "Delta"
        )
      )
    )
  )
  testthat::expect_invisible(
    plot(traj, main = "With delta panel", show_delta = "panel", delta_limits = c(-5, 15))
  )
  testthat::expect_invisible(
    plot(traj, main = "Without delta or legend", show_delta = "none", legend = FALSE)
  )
  testthat::expect_invisible(
    plot(
      icTrajectory(list(
        baseline_ic = -1000,
        model_fit_history = list(fits = list(
          list(ic = -1005, accepted = FALSE, status = "error", delta_ic = -5)
        ))
      )),
      main = "Finite error point",
      show_delta = "none",
      legend = FALSE
    )
  )
})

testthat::test_that("plot.icTrajectory parses grouped display controls", {
  testthat::expect_equal(
    c(
      true_delta = .icTrajectory_show_delta_mode(TRUE),
      false_delta = .icTrajectory_show_delta_mode(FALSE)
    ),
    c(
      true_delta = "overlay",
      false_delta = "none"
    )
  )

  keys <- c("running_best", "accepted", "rejected", "baseline", "delta")
  unnamed_labels <- .icTrajectory_legend_labels(
    c("Best", "Accepted", "Rejected", "Baseline", "Delta"),
    keys
  )
  named_labels <- .icTrajectory_legend_labels(
    c(accepted = "Accepted shift"),
    keys
  )
  testthat::expect_equal(
    unnamed_labels,
    c("Best", "Accepted", "Rejected", "Baseline", "Delta")
  )
  testthat::expect_true(is.expression(named_labels))
  testthat::expect_equal(
    as.list(named_labels)[1:4],
    list(
      "Running best",
      "Accepted shift",
      "Rejected proposal",
      "Baseline IC"
    )
  )
  testthat::expect_identical(as.list(named_labels)[[5L]], quote(Delta ~ IC))
  expression_labels <- .icTrajectory_legend_label_vector(
    list("Plain label", expression(Delta ~ IC))
  )
  testthat::expect_true(is.expression(expression_labels))
  testthat::expect_identical(as.list(expression_labels)[[2L]], quote(Delta ~ IC))

  style <- .icTrajectory_plot_style(
    symbols = c(
      accepted = 19,
      rejected = 4,
      error = 8,
      baseline = 20,
      delta = 2
    ),
    scales = c(point = 1.2, line = 1.1, text = 0.9),
    point_sizes = c(
      accepted = 0.7,
      rejected = 0.4,
      error = 0.5,
      baseline = 1,
      delta = 0.3
    ),
    line_widths = c(
      accepted = 0.2,
      rejected = 0.3,
      error = 0.4,
      running_best = 1.5,
      delta = 0.6,
      zero = 0.7
    ),
    text_sizes = c(axis = 0.6, y_axis = 0.7, legend = 0.8),
    annotation = c(baseline_label_offset = 0.75),
    legend = list(
      show = FALSE,
      position = c(0.2, 0.8),
      inset = c(0.01, 0.02),
      bty = "o",
      labels = c(delta = "Delta")
    )
  )
  axis_style <- .icTrajectory_plot_style(text_sizes = c(axis = 0.6))
  axis_label_style <- .icTrajectory_plot_style(
    text_sizes = c(axis_label = 1.4, y_axis = 0.9)
  )

  testthat::expect_equal(
    c(
      style[c(
        "accepted_pch", "error_pch", "accepted_cex", "running_best_lwd",
        "x_axis_cex", "y_axis_cex", "delta_axis_cex",
        "baseline_label_offset", "legend", "delta_legend_label"
      )],
      list(
        legend_position = style$legend_position,
        legend_inset = style$legend_inset,
        null_legend_show = .icTrajectory_legend_group(NULL)$show,
        cascaded_y_axis = axis_style$y_axis_cex,
        cascaded_delta_axis = axis_style$delta_axis_cex,
        axis_label_x_axis = axis_label_style$x_axis_cex,
        axis_label_delta_axis = axis_label_style$delta_axis_cex
      )
    ),
    list(
      accepted_pch = 19,
      error_pch = 8,
      accepted_cex = 0.7,
      running_best_lwd = 1.5,
      x_axis_cex = 0.6,
      y_axis_cex = 0.7,
      delta_axis_cex = 0.7,
      baseline_label_offset = 0.75,
      legend = FALSE,
      delta_legend_label = quote(Delta ~ IC),
      legend_position = c(0.2, 0.8),
      legend_inset = c(0.01, 0.02),
      null_legend_show = TRUE,
      cascaded_y_axis = 0.6,
      cascaded_delta_axis = 0.6,
      axis_label_x_axis = 0.8,
      axis_label_delta_axis = 0.9
    )
  )
})

testthat::test_that("icTrajectory scalar and vector validators coerce valid inputs", {
  testthat::expect_equal(.icTrajectory_scalar_numeric(NULL, "x"), NA_real_)
  testthat::expect_equal(.icTrajectory_numeric_vector(factor(c("1", "2")), "x"), c(1, 2))
  testthat::expect_equal(.icTrajectory_logical_vector(c(TRUE, FALSE, NA), "x"), c(TRUE, FALSE, NA))
})

testthat::test_that("icTrajectory scalar and vector validators reject invalid inputs", {
  testthat::expect_error(.icTrajectory_scalar_numeric("x", "x"), "`x` must be numeric")
  testthat::expect_error(.icTrajectory_numeric_vector("x", "x"), "`x` must be numeric")
  testthat::expect_error(
    .icTrajectory_logical_vector(c(0, 2), "x"),
    "`x` must contain logical values"
  )
  testthat::expect_error(
    .icTrajectory_logical_vector("maybe", "x"),
    "`x` must contain logical values"
  )
  testthat::expect_error(
    .icTrajectory_logical_vector(c("1", "0"), "x"),
    "`x` must contain logical values"
  )
})

testthat::test_that("icTrajectory axis and label helpers apply documented fallbacks", {
  testthat::expect_equal(.icTrajectory_delta_limits(c(NA_real_, NA_real_), NULL), c(-1, 1))
  testthat::expect_equal(.icTrajectory_delta_limits(c(0, 0), NULL), c(-1, 1))
  testthat::expect_equal(.icTrajectory_delta_limits(c(-2, 4), c(4, -2)), c(4, -2))
  testthat::expect_equal(.icTrajectory_y_label(NULL, NA_character_), "IC score")
  testthat::expect_equal(.icTrajectory_axis_limits(c(0, 0)), c(-1, 1))
})

testthat::test_that("icTrajectory legend helpers validate grouped controls", {
  testthat::expect_error(
    .icTrajectory_named_group(environment(), "a", "group"),
    "`group` must be a named vector or list"
  )
  testthat::expect_error(
    .icTrajectory_legend_group(1),
    "`legend` must be TRUE, FALSE, or a named list"
  )
  testthat::expect_error(
    .icTrajectory_legend_group(list(TRUE)),
    "`legend` entries must be named"
  )
  testthat::expect_error(
    .icTrajectory_legend_group(list(foo = TRUE)),
    "Named `legend` entries must match known keys"
  )
  testthat::expect_error(
    .icTrajectory_legend_inset(NA_real_),
    "`legend_inset` must be a finite numeric vector"
  )
  testthat::expect_error(
    .icTrajectory_legend_labels_arg(NA_character_),
    "`legend_labels` must be a character vector"
  )
  testthat::expect_error(
    .icTrajectory_legend_labels(c("Only one"), c("running_best", "accepted")),
    "Unnamed `legend_labels` must match"
  )
})

testthat::test_that("icTrajectory plotting helpers handle sparse delta data", {
  traj <- icTrajectory(make_search_for_trajectory())
  traj$delta_ic[] <- NA_real_
  style <- .icTrajectory_plot_style(legend = FALSE)

  open_null_device_trajectory()
  on.exit(close_device_trajectory(), add = TRUE)

  testthat::expect_invisible(
    .icTrajectory_plot_delta_panel(
      x = traj,
      x_limits = c(0, 3),
      x_ticks = 0:3,
      delta_limits = NULL,
      xlab = "Step",
      style = style
    )
  )
  testthat::expect_error(
    plot(traj, show_delta = "none", legend = FALSE, ic_limits = c(NA_real_, NA_real_)),
    "`ic_limits` must be a numeric vector of length 2"
  )

  no_finite <- traj
  no_finite$ic[] <- NA_real_
  no_finite$best_ic[] <- NA_real_
  testthat::expect_error(
    plot(no_finite),
    "`x` must contain at least one finite IC value"
  )
})

testthat::test_that("plot.icTrajectory preserves existing layout outside panel mode", {
  traj <- icTrajectory(make_search_for_trajectory())

  open_null_device_trajectory()
  on.exit(close_device_trajectory(), add = TRUE)

  layout(matrix(1:2, nrow = 1L))
  plot(traj, show_delta = "none", legend = FALSE)

  testthat::expect_equal(par("mfg")[3:4], c(1L, 2L))
})

testthat::test_that("plot.icTrajectory validates plotting arguments", {
  traj <- icTrajectory(make_search_for_trajectory())

  open_null_device_trajectory()
  on.exit(close_device_trajectory(), add = TRUE)

  testthat::expect_error(
    plot(traj, show_delta = NA),
    "`show_delta` must be one of"
  )
  testthat::expect_error(
    plot(traj, delta_orientation = "down"),
    "Unused plotting arguments"
  )
  testthat::expect_error(
    plot(traj, delta_limits = c(NA_real_, 1)),
    "`delta_limits` must be a numeric vector of length 2"
  )
  testthat::expect_error(
    plot(traj, ic_limits = c(NA_real_, 1)),
    "`ic_limits` must be a numeric vector of length 2"
  )
  testthat::expect_error(
    plot(traj, scales = c(point = -1)),
    "`scales\\$point` must be a finite non-negative numeric scalar"
  )
  testthat::expect_error(
    plot(traj, text_sizes = c(y_axis = -1)),
    "`text_sizes\\$y_axis` must be a finite non-negative numeric scalar"
  )
  testthat::expect_error(
    plot(traj, annotation = c(baseline_label_offset = -1)),
    "`annotation\\$baseline_label_offset` must be a finite non-negative numeric scalar"
  )
  testthat::expect_error(
    plot(traj, symbols = list(accepted = c(1, 2))),
    "`symbols\\$accepted` must be a single numeric or character plotting symbol"
  )
  testthat::expect_error(
    plot(traj, symbols = list(accepted = list(1))),
    "`symbols\\$accepted` must be a single numeric or character plotting symbol"
  )
  testthat::expect_error(
    plot(traj, legend = NA),
    "`legend` must be TRUE or FALSE"
  )
  testthat::expect_error(
    plot(traj, legend = list(position = 1)),
    "`legend_position` must be one of"
  )
  testthat::expect_error(
    plot(traj, legend = list(position = "notaposition")),
    "`legend_position` must be one of"
  )
  testthat::expect_error(
    plot(traj, legend = list(bty = "x")),
    "`legend_bty` must be"
  )
  testthat::expect_error(
    plot(traj, legend = list(labels = c("Best", accepted = "Accepted"))),
    "`legend_labels` must be either fully named or unnamed"
  )
  testthat::expect_error(
    plot(traj, legend = list(labels = c(foo = "Foo"))),
    "Named `legend_labels` entries"
  )
  testthat::expect_error(
    plot(traj, line_widths = c(foo = 1)),
    "Named `line_widths` entries"
  )
  testthat::expect_error(
    plot(traj, point_sizes = c(1)),
    "`point_sizes` entries must be named"
  )
  testthat::expect_error(
    plot(traj, accepted_cex = 1),
    "Unused plotting arguments"
  )
  testthat::expect_error(
    plot.icTrajectory(make_search_for_trajectory(), legend = FALSE),
    "`x` must be an `icTrajectory` object"
  )
})

testthat::test_that("plot_ic_acceptance_matrix remains as a compatibility wrapper", {
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
      symbols = c(accepted = 19),
      scales = c(point = 1.2),
      legend = FALSE
    )
  )
})

testthat::test_that("plot_ic_acceptance_matrix validates legacy arguments", {
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
