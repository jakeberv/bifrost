#' Create an Information-Criterion Trajectory
#'
#' @description
#' Extract a standardized information-criterion (IC) trajectory from a
#' completed `bifrost_search` object or a compatible search-result list.
#'
#' @details
#' `icTrajectory()` is a lightweight extractor. It does not refit models or
#' recompute the search. The baseline IC is taken from `x$baseline_ic` unless
#' `baseline_ic` is supplied, and is always included as `step = 0`; later rows
#' summarize the stored proposal history in `x$model_fit_history$fits`. Legacy
#' objects that only contain `x$model_fit_history$ic_acceptance_matrix`, or
#' legacy `fits` entries without explicit proposal metadata, are also
#' supported. In those legacy cases, `candidate_node` cannot be recovered, and
#' `regime_id` is inferred from proposal order because Bifrost assigns
#' candidate shift regimes sequentially during the search.
#'
#' The returned object is a data frame with class
#' `c("icTrajectory", "data.frame")` and the following columns:
#' \itemize{
#'   \item `step`: search step, with `0` reserved for the baseline.
#'   \item `ic`: IC for the baseline or evaluated proposal.
#'   \item `accepted`: proposal decision; `NA` for the baseline.
#'   \item `best_ic`: running best IC after each step.
#'   \item `delta_ic`: proposal improvement over the current best before the
#'     proposal; positive values indicate improvement.
#'   \item `status`: `"baseline"`, `"accepted"`, `"rejected"`, or `"error"`.
#'   \item `candidate_node`: phylogenetic node where the shift was evaluated.
#'   \item `regime_id`: mapped-state/regime label assigned to the evaluated
#'     candidate shift. The baseline row uses `"0"`.
#'   \item `error`: optional error message.
#' }
#'
#' @param x A `bifrost_search` object returned by [searchOptimalConfiguration()]
#'   with `store_model_fit_history = TRUE`, or a compatible list with
#'   `baseline_ic` and `model_fit_history`.
#' @param baseline_ic Optional finite numeric baseline IC. When supplied, this
#'   value is used for `step = 0` instead of `x$baseline_ic`. This is useful for
#'   legacy search-like lists that did not store `baseline_ic`.
#' @param ... Reserved for future extensions.
#'
#' @return A data frame with class `c("icTrajectory", "data.frame")`.
#'
#' @examples
#' search <- list(
#'   baseline_ic = -1000,
#'   IC_used = "GIC",
#'   model_fit_history = list(
#'     fits = list(
#'       list(step = 1, candidate_node = 42, regime_id = "1", ic = -1010, accepted = TRUE, delta_ic = 10),
#'       list(step = 2, candidate_node = 57, regime_id = "2", ic = -1008, accepted = FALSE, delta_ic = -2)
#'     )
#'   )
#' )
#' class(search) <- c("bifrost_search", "list")
#' traj <- icTrajectory(search)
#' legacy_traj <- icTrajectory(search, baseline_ic = -995)
#' plot(traj)
#'
#' @export
icTrajectory <- function(x, baseline_ic = NULL, ...) {
  UseMethod("icTrajectory")
}

#' @export
icTrajectory.bifrost_search <- function(x, baseline_ic = NULL, ...) {
  .icTrajectory_from_search_like(x, baseline_ic = baseline_ic)
}

#' @export
icTrajectory.list <- function(x, baseline_ic = NULL, ...) {
  .icTrajectory_from_search_like(x, baseline_ic = baseline_ic)
}

#' @export
icTrajectory.default <- function(x, baseline_ic = NULL, ...) {
  if (inherits(x, "icTrajectory")) {
    return(x)
  }
  if (is.list(x)) {
    return(.icTrajectory_from_search_like(x, baseline_ic = baseline_ic))
  }
  stop("`x` must be a `bifrost_search` object or compatible search-result list.")
}

.icTrajectory_from_search_like <- function(x, baseline_ic = NULL) {
  baseline_overridden <- !is.null(baseline_ic)
  baseline_ic <- if (is.null(baseline_ic)) {
    .icTrajectory_scalar_numeric(x$baseline_ic, "x$baseline_ic")
  } else {
    .icTrajectory_scalar_numeric(baseline_ic, "baseline_ic")
  }
  if (!is.finite(baseline_ic)) {
    stop("`baseline_ic` must be a finite numeric scalar or `x$baseline_ic` must be present.")
  }

  IC <- if (is.null(x$IC_used) || length(x$IC_used) == 0L) {
    NA_character_
  } else {
    as.character(x$IC_used[1L])
  }

  history <- x$model_fit_history
  if (!is.list(history)) {
    stop(
      "`x$model_fit_history` is required. ",
      "Run `searchOptimalConfiguration()` with `store_model_fit_history = TRUE`."
    )
  }

  if (!is.list(history$fits)) {
    if (!is.null(history$ic_acceptance_matrix)) {
      return(.icTrajectory_from_legacy_matrix(
        history$ic_acceptance_matrix,
        baseline_ic = baseline_ic,
        IC = IC
      ))
    }
    stop(
      "`x$model_fit_history$fits` or ",
      "`x$model_fit_history$ic_acceptance_matrix` is required."
    )
  }

  entries <- history$fits
  n_entries <- length(entries)

  out <- data.frame(
    step = integer(n_entries + 1L),
    ic = numeric(n_entries + 1L),
    accepted = rep(NA, n_entries + 1L),
    best_ic = numeric(n_entries + 1L),
    delta_ic = rep(NA_real_, n_entries + 1L),
    status = character(n_entries + 1L),
    candidate_node = rep(NA_integer_, n_entries + 1L),
    regime_id = rep(NA_character_, n_entries + 1L),
    error = rep(NA_character_, n_entries + 1L),
    stringsAsFactors = FALSE
  )

  out[1L, ] <- list(
    step = 0L,
    ic = baseline_ic,
    accepted = NA,
    best_ic = baseline_ic,
    delta_ic = NA_real_,
    status = "baseline",
    candidate_node = NA_integer_,
    regime_id = "0",
    error = NA_character_
  )

  current_best <- baseline_ic
  if (n_entries > 0L) {
    for (i in seq_along(entries)) {
      entry <- entries[[i]]
      if (!is.list(entry)) {
        entry <- list()
      }

      step <- .icTrajectory_entry_step(entry, i)
      ic <- .icTrajectory_entry_ic(entry, IC)
      accepted <- .icTrajectory_entry_accepted(entry)
      error <- .icTrajectory_entry_error(entry)
      status <- .icTrajectory_entry_status(entry, accepted, ic, error)
      delta_ic <- .icTrajectory_entry_delta(
        entry,
        current_best,
        ic,
        prefer_stored = !baseline_overridden
      )
      candidate_node <- .icTrajectory_entry_candidate_node(entry)
      regime_id <- .icTrajectory_entry_regime_id(entry, fallback = as.character(i))

      if (identical(status, "accepted") && is.finite(ic)) {
        current_best <- ic
      }

      out[i + 1L, ] <- list(
        step = step,
        ic = ic,
        accepted = accepted,
        best_ic = current_best,
        delta_ic = delta_ic,
        status = status,
        candidate_node = candidate_node,
        regime_id = regime_id,
        error = error
      )
    }
  }

  class(out) <- c("icTrajectory", "data.frame")
  attr(out, "IC_used") <- IC
  out
}

.icTrajectory_scalar_numeric <- function(x, name) {
  if (is.null(x) || length(x) == 0L) {
    return(NA_real_)
  }
  out <- suppressWarnings(as.numeric(x[1L]))
  if (is.na(out) && !is.na(x[1L])) {
    stop("`", name, "` must be numeric.")
  }
  out
}

.icTrajectory_numeric_vector <- function(x, name) {
  if (is.factor(x)) {
    x <- as.character(x)
  }
  out <- suppressWarnings(as.numeric(x))
  invalid <- !is.na(x) & is.na(out)
  if (any(invalid)) {
    stop("`", name, "` must be numeric.")
  }
  out
}

.icTrajectory_logical_vector <- function(x, name) {
  if (is.logical(x)) {
    return(as.logical(x))
  }
  if (is.factor(x)) {
    x <- as.character(x)
  }
  if (is.numeric(x) || is.integer(x)) {
    out <- rep(NA, length(x))
    out[x == 1] <- TRUE
    out[x == 0] <- FALSE
    invalid <- !is.na(x) & !(x %in% c(0, 1))
    if (any(invalid)) {
      stop("`", name, "` must contain logical values or 0/1 indicators.")
    }
    return(out)
  }

  x_chr <- tolower(trimws(as.character(x)))
  out <- rep(NA, length(x_chr))
  out[x_chr %in% c("1", "t", "true", "yes", "accepted")] <- TRUE
  out[x_chr %in% c("0", "f", "false", "no", "rejected")] <- FALSE
  invalid <- !is.na(x_chr) & nzchar(x_chr) & is.na(out)
  if (any(invalid)) {
    stop("`", name, "` must contain logical values or 0/1 indicators.")
  }
  out
}

.icTrajectory_from_legacy_matrix <- function(matrix_data, baseline_ic, IC) {
  data <- as.data.frame(matrix_data, stringsAsFactors = FALSE)
  if (nrow(data) == 0L || ncol(data) < 2L) {
    stop("`x$model_fit_history$ic_acceptance_matrix` must have at least two columns.")
  }

  ic <- .icTrajectory_numeric_vector(
    data[[1L]],
    "x$model_fit_history$ic_acceptance_matrix[, 1]"
  )
  accepted <- .icTrajectory_logical_vector(
    data[[2L]],
    "x$model_fit_history$ic_acceptance_matrix[, 2]"
  )

  n_entries <- length(ic)
  out <- data.frame(
    step = 0:n_entries,
    ic = c(baseline_ic, ic),
    accepted = c(NA, accepted),
    best_ic = rep(NA_real_, n_entries + 1L),
    delta_ic = rep(NA_real_, n_entries + 1L),
    status = rep(NA_character_, n_entries + 1L),
    candidate_node = rep(NA_integer_, n_entries + 1L),
    regime_id = c("0", as.character(seq_len(n_entries))),
    error = rep(NA_character_, n_entries + 1L),
    stringsAsFactors = FALSE
  )

  out$status[1L] <- "baseline"
  out$best_ic[1L] <- baseline_ic
  current_best <- baseline_ic

  if (n_entries > 0L) {
    for (i in seq_len(n_entries)) {
      row <- i + 1L
      if (!is.finite(out$ic[row])) {
        out$status[row] <- "error"
        out$accepted[row] <- FALSE
      } else if (isTRUE(out$accepted[row])) {
        out$status[row] <- "accepted"
      } else {
        out$status[row] <- "rejected"
      }

      if (is.finite(current_best) && is.finite(out$ic[row])) {
        out$delta_ic[row] <- current_best - out$ic[row]
      }
      if (identical(out$status[row], "accepted") && is.finite(out$ic[row])) {
        current_best <- out$ic[row]
      }
      out$best_ic[row] <- current_best
    }
  }

  class(out) <- c("icTrajectory", "data.frame")
  attr(out, "IC_used") <- IC
  out
}

.icTrajectory_entry_step <- function(entry, fallback) {
  step <- if (is.null(entry$step)) fallback else entry$step[1L]
  step <- suppressWarnings(as.integer(step))
  if (is.na(step)) {
    fallback
  } else {
    step
  }
}

.icTrajectory_entry_ic <- function(entry, IC) {
  if (!is.null(entry$ic)) {
    ic <- suppressWarnings(as.numeric(entry$ic[1L]))
    if (is.finite(ic)) {
      return(ic)
    }
  }

  if (!is.null(entry$model) && !is.na(IC)) {
    return(.bifrost_search_ic_value(entry$model, IC))
  }

  NA_real_
}

.icTrajectory_entry_accepted <- function(entry) {
  accepted <- entry$accepted
  if (is.null(accepted) || length(accepted) == 0L) {
    return(NA)
  }
  if (is.logical(accepted)) {
    return(as.logical(accepted[1L]))
  }
  accepted_num <- suppressWarnings(as.numeric(accepted[1L]))
  if (is.na(accepted_num)) {
    return(NA)
  }
  accepted_num == 1
}

.icTrajectory_entry_error <- function(entry) {
  if (is.null(entry$error) || length(entry$error) == 0L) {
    return(NA_character_)
  }
  error <- as.character(entry$error[1L])
  if (!nzchar(error)) {
    NA_character_
  } else {
    error
  }
}

.icTrajectory_entry_status <- function(entry, accepted, ic, error) {
  if (!is.null(entry$status) && length(entry$status) > 0L) {
    status <- tolower(as.character(entry$status[1L]))
    if (status %in% c("accepted", "rejected", "error")) {
      return(status)
    }
  }
  if (!is.na(error) || !is.finite(ic)) {
    return("error")
  }
  if (isTRUE(accepted)) {
    return("accepted")
  }
  "rejected"
}

.icTrajectory_entry_delta <- function(entry, current_best, ic, prefer_stored = TRUE) {
  if (isTRUE(prefer_stored) && !is.null(entry$delta_ic)) {
    delta_ic <- suppressWarnings(as.numeric(entry$delta_ic[1L]))
    if (is.finite(delta_ic)) {
      return(delta_ic)
    }
  }
  if (is.finite(current_best) && is.finite(ic)) {
    return(current_best - ic)
  }
  NA_real_
}

.icTrajectory_entry_candidate_node <- function(entry) {
  candidate_node <- NULL
  if (!is.null(entry$candidate_node)) {
    candidate_node <- entry$candidate_node[1L]
  } else if (!is.null(entry$proposal)) {
    candidate_node <- entry$proposal[1L]
  } else if (!is.null(entry$node)) {
    candidate_node <- entry$node[1L]
  }

  if (is.null(candidate_node) || length(candidate_node) == 0L || is.na(candidate_node)) {
    return(NA_integer_)
  }
  candidate_node <- suppressWarnings(as.integer(candidate_node))
  if (is.na(candidate_node)) {
    NA_integer_
  } else {
    candidate_node
  }
}

.icTrajectory_entry_regime_id <- function(entry, fallback = NA_character_) {
  regime_id <- NULL
  if (!is.null(entry$regime_id)) {
    regime_id <- entry$regime_id[1L]
  } else if (!is.null(entry$shift_id)) {
    regime_id <- entry$shift_id[1L]
  }

  if (is.null(regime_id) || length(regime_id) == 0L || is.na(regime_id)) {
    return(as.character(fallback[1L]))
  }
  regime_id <- as.character(regime_id)
  if (!nzchar(regime_id)) {
    NA_character_
  } else {
    regime_id
  }
}

.icTrajectory_delta_limits <- function(delta_ic, delta_limits) {
  if (!is.null(delta_limits)) {
    if (!is.numeric(delta_limits) ||
        length(delta_limits) != 2L ||
        any(!is.finite(delta_limits))) {
      stop("`delta_limits` must be a numeric vector of length 2 with finite values.")
    }
    return(sort(as.numeric(delta_limits)))
  }

  finite_delta <- delta_ic[is.finite(delta_ic)]
  if (length(finite_delta) == 0L) {
    return(c(-1, 1))
  }

  limits <- range(c(0, finite_delta))
  if (!is.finite(diff(limits)) || diff(limits) == 0) {
    pad <- if (limits[1L] == 0) 1 else abs(limits[1L]) * 0.1
    limits <- limits + c(-pad, pad)
  }
  range(pretty(limits))
}

#' Plot a Legacy IC Acceptance Matrix
#'
#' @description
#' Superseded compatibility wrapper for plotting the historical two-column IC
#' acceptance matrix. New code should usually call [icTrajectory()] and then
#' [plot.icTrajectory()] directly.
#'
#' @details
#' This wrapper preserves the old `plot_ic_acceptance_matrix()` argument names
#' while delegating the actual drawing to [plot.icTrajectory()]. `matrix_data`
#' must have IC values in the first column and accepted/rejected indicators in
#' the second column. To match the legacy plotting convention, the first row is
#' treated as the plotted baseline row; when `baseline_ic` is supplied, it
#' replaces the first IC value for that baseline row.
#'
#' The historical `plot_rate_of_improvement` argument maps to the new
#' `show_delta` plotting mode: `TRUE` draws the `delta_ic` overlay and `FALSE`
#' suppresses it. `rate_limits` is retained as the legacy name for the overlay
#' y-limits and is passed through as `delta_limits`.
#'
#' @param matrix_data A two-column matrix or data frame. Column 1 must contain
#'   numeric IC scores in evaluation order; column 2 must contain logical values
#'   or 0/1 indicators for accepted proposals.
#' @param plot_title Plot title.
#' @param plot_rate_of_improvement Logical; if `TRUE`, draw the proposal
#'   `delta_ic` overlay on a secondary y-axis.
#' @param rate_limits Numeric vector of length 2 giving limits for the overlay
#'   axis. Used only when `plot_rate_of_improvement = TRUE`.
#' @param baseline_ic Optional finite numeric baseline IC. When supplied, this
#'   replaces the first IC value for the plotted baseline row.
#' @param ... Additional plotting arguments passed to [plot.icTrajectory()].
#'
#' @return Invisibly returns `NULL`. Called for its plotting side effects.
#'
#' @seealso [icTrajectory()], [plot.icTrajectory()]
#'
#' @export
plot_ic_acceptance_matrix <- function(matrix_data,
                                      plot_title = "IC Acceptance Matrix Scatter Plot",
                                      plot_rate_of_improvement = TRUE,
                                      rate_limits = c(-400, 150),
                                      baseline_ic = NULL,
                                      ...) {
  show_delta <- if (.icTrajectory_logical_scalar(
    plot_rate_of_improvement,
    "plot_rate_of_improvement"
  )) {
    "overlay"
  } else {
    "none"
  }
  delta_limits <- if (identical(show_delta, "overlay")) {
    .icTrajectory_legacy_rate_limits(rate_limits)
  } else {
    NULL
  }

  traj <- .icTrajectory_from_legacy_plot_matrix(
    matrix_data = matrix_data,
    baseline_ic = baseline_ic
  )
  plot(
    traj,
    main = plot_title,
    show_delta = show_delta,
    delta_limits = delta_limits,
    xlab = "Sub-model evaluated",
    ylab = "IC Score",
    ...
  )

  invisible(NULL)
}

.icTrajectory_from_legacy_plot_matrix <- function(matrix_data, baseline_ic = NULL) {
  data <- as.data.frame(matrix_data, stringsAsFactors = FALSE)
  if (nrow(data) == 0L || ncol(data) < 2L) {
    stop("`matrix_data` must have at least one row and two columns.")
  }

  ic <- .icTrajectory_numeric_vector(data[[1L]], "matrix_data[, 1]")
  accepted <- .icTrajectory_logical_vector(data[[2L]], "matrix_data[, 2]")

  baseline_value <- if (is.null(baseline_ic)) {
    ic[1L]
  } else {
    .icTrajectory_scalar_numeric(baseline_ic, "baseline_ic")
  }
  if (!is.finite(baseline_value)) {
    stop("`baseline_ic` must be finite or the first IC value in `matrix_data` must be finite.")
  }

  fits <- list()
  if (length(ic) > 1L) {
    fits <- vector("list", length(ic) - 1L)
    for (i in seq_along(fits)) {
      row <- i + 1L
      fits[[i]] <- list(
        step = i,
        ic = ic[row],
        accepted = accepted[row],
        regime_id = as.character(i)
      )
    }
  }

  traj <- icTrajectory(list(
    baseline_ic = baseline_value,
    model_fit_history = list(fits = fits)
  ))
  traj$step <- seq_len(nrow(traj))
  traj
}

.icTrajectory_legacy_rate_limits <- function(rate_limits) {
  if (!is.numeric(rate_limits) ||
      length(rate_limits) != 2L ||
      any(!is.finite(rate_limits))) {
    stop("`rate_limits` must be a numeric vector of length 2 with finite values.")
  }
  sort(as.numeric(rate_limits))
}

#' Plot an Information-Criterion Trajectory
#'
#' @description
#' Plot IC scores and the running best IC across a stored Bifrost search
#' trajectory.
#'
#' @param x An `icTrajectory` object.
#' @param main Plot title.
#' @param show_delta One of `"overlay"`, `"panel"`, or `"none"`. The default
#'   `"overlay"` draws proposal `delta_ic` values on a secondary y-axis.
#'   `"panel"` draws them in a lower panel. Logical values are accepted for
#'   compatibility: `TRUE` maps to `"overlay"` and `FALSE` maps to `"none"`.
#' @param delta_limits Optional numeric vector of length 2 giving y-limits for
#'   the `delta_ic` panel or secondary axis.
#' @param xlab,ylab Axis labels.
#' @param accepted_pch,rejected_pch,error_pch,baseline_pch,delta_pch Point
#'   symbols for accepted proposals, rejected proposals, errored proposals, the
#'   baseline IC, and `delta_ic` points.
#' @param point_scale,line_scale,text_scale Multipliers for point sizes, line
#'   widths, and text sizes.
#' @param accepted_cex,rejected_cex,error_cex,baseline_cex,delta_cex Point-size
#'   controls for accepted proposals, rejected proposals, errored proposals, the
#'   baseline IC, and `delta_ic` points.
#' @param accepted_lwd,rejected_lwd,error_lwd,running_best_lwd,delta_lwd,zero_lwd
#'   Line-width controls for proposal points, the running-best line, the
#'   `delta_ic` trace or stems, and the zero-reference line.
#' @param main_cex,axis_cex,axis_label_cex,annotation_cex,legend_cex Text-size
#'   controls for the title, tick labels, axis labels, point annotations, and
#'   legend.
#' @param legend Logical; if `TRUE`, draw a legend.
#' @param legend_position Legend position passed to [graphics::legend()]. Use a
#'   base R legend keyword such as `"topright"` or a numeric vector of length 2.
#' @param legend_inset Numeric inset passed to [graphics::legend()].
#' @param legend_bty Legend box type passed to [graphics::legend()], typically
#'   `"n"` or `"o"`.
#' @param legend_labels Optional custom legend labels. Provide either a full
#'   unnamed character vector in display order (`running_best`, `accepted`,
#'   `rejected`, `baseline`, and `delta` when shown), or a fully named character
#'   vector replacing any of `running_best`, `accepted`, `rejected`, `baseline`,
#'   and `delta`.
#' @param ... Unused; included for S3 compatibility.
#'
#' @return Invisibly returns `x`.
#'
#' @importFrom graphics axis layout legend lines mtext par plot points segments text title
#' @importFrom grDevices rgb
#' @export
plot.icTrajectory <- function(x,
                              main = "IC Trajectory",
                              show_delta = "overlay",
                              delta_limits = NULL,
                              xlab = "Search step",
                              ylab = NULL,
                              accepted_pch = 21,
                              rejected_pch = 3,
                              error_pch = 4,
                              baseline_pch = 19,
                              delta_pch = 16,
                              point_scale = 1,
                              line_scale = 1,
                              text_scale = 1,
                              accepted_cex = 0.68 * point_scale,
                              rejected_cex = 0.36 * point_scale,
                              error_cex = 0.55 * point_scale,
                              baseline_cex = 1.0 * point_scale,
                              delta_cex = 0.30 * point_scale,
                              accepted_lwd = 0.30 * line_scale,
                              rejected_lwd = 0.45 * line_scale,
                              error_lwd = 0.60 * line_scale,
                              running_best_lwd = 1.0 * line_scale,
                              delta_lwd = 0.80 * line_scale,
                              zero_lwd = 0.70 * line_scale,
                              main_cex = 1.0 * text_scale,
                              axis_cex = 0.80 * text_scale,
                              axis_label_cex = 0.60 * text_scale,
                              annotation_cex = 0.60 * text_scale,
                              legend_cex = 0.68 * text_scale,
                              legend = TRUE,
                              legend_position = "topright",
                              legend_inset = 0,
                              legend_bty = "n",
                              legend_labels = NULL,
                              ...) {
  if (!inherits(x, "icTrajectory")) {
    x <- icTrajectory(x)
  }
  show_delta <- .icTrajectory_show_delta_mode(show_delta)
  style <- .icTrajectory_plot_style(
    accepted_pch = accepted_pch,
    rejected_pch = rejected_pch,
    error_pch = error_pch,
    baseline_pch = baseline_pch,
    delta_pch = delta_pch,
    point_scale = point_scale,
    line_scale = line_scale,
    text_scale = text_scale,
    accepted_cex = accepted_cex,
    rejected_cex = rejected_cex,
    error_cex = error_cex,
    baseline_cex = baseline_cex,
    delta_cex = delta_cex,
    accepted_lwd = accepted_lwd,
    rejected_lwd = rejected_lwd,
    error_lwd = error_lwd,
    running_best_lwd = running_best_lwd,
    delta_lwd = delta_lwd,
    zero_lwd = zero_lwd,
    main_cex = main_cex,
    axis_cex = axis_cex,
    axis_label_cex = axis_label_cex,
    annotation_cex = annotation_cex,
    legend_cex = legend_cex,
    legend = legend,
    legend_position = legend_position,
    legend_inset = legend_inset,
    legend_bty = legend_bty,
    legend_labels = legend_labels
  )

  oldpar <- par(no.readonly = TRUE)
  used_layout <- FALSE
  on.exit({
    if (isTRUE(used_layout)) {
      layout(1)
    }
    par(oldpar)
  }, add = TRUE)

  x_values <- x$step
  y_values <- x$ic
  finite_y <- c(y_values, x$best_ic)
  finite_y <- finite_y[is.finite(finite_y)]
  if (length(finite_y) == 0L) {
    stop("`x` must contain at least one finite IC value.")
  }

  x_ticks <- pretty(x_values)
  x_limits <- range(x_ticks)
  y_limits <- .icTrajectory_axis_limits(finite_y)
  y_ticks <- pretty(y_limits)
  IC <- attr(x, "IC_used")
  ylab <- .icTrajectory_y_label(ylab, IC)

  draw_delta <- !identical(show_delta, "none") && any(is.finite(x$delta_ic))

  if (identical(show_delta, "panel") && draw_delta) {
    used_layout <- TRUE
    layout(matrix(c(1, 2), ncol = 1L), heights = c(3, 1.25))
    par(mar = c(0.8, 5.5, 3.5, 1.5), mgp = c(3, 0.6, 0))
    .icTrajectory_plot_main_panel(
      x = x,
      x_limits = x_limits,
      y_limits = y_limits,
      x_ticks = x_ticks,
      y_ticks = y_ticks,
      main = main,
      xlab = "",
      ylab = ylab,
      draw_x_labels = FALSE,
      include_delta_legend = FALSE,
      style = style
    )
    par(mar = c(4.2, 5.5, 0.8, 1.5), mgp = c(3, 0.6, 0))
    .icTrajectory_plot_delta_panel(
      x = x,
      x_limits = x_limits,
      x_ticks = x_ticks,
      delta_limits = delta_limits,
      xlab = xlab,
      style = style
    )
    return(invisible(x))
  }

  par(
    mar = if (identical(show_delta, "overlay") && draw_delta) c(5, 5.5, 4, 6) else c(5, 5.5, 4, 1.5),
    mgp = c(3, 0.6, 0)
  )
  if (identical(show_delta, "overlay") && draw_delta) {
    .icTrajectory_plot_delta_overlay(
      x = x,
      x_limits = x_limits,
      x_ticks = x_ticks,
      delta_limits = delta_limits,
      style = style
    )
    par(new = TRUE)
  }

  .icTrajectory_plot_main_panel(
    x = x,
    x_limits = x_limits,
    y_limits = y_limits,
    x_ticks = x_ticks,
    y_ticks = y_ticks,
    main = main,
    xlab = xlab,
    ylab = ylab,
    draw_x_labels = TRUE,
    include_delta_legend = identical(show_delta, "overlay") && draw_delta,
    style = style
  )

  invisible(x)
}

.icTrajectory_show_delta_mode <- function(show_delta) {
  if (is.logical(show_delta) && length(show_delta) == 1L && !is.na(show_delta)) {
    return(if (isTRUE(show_delta)) "overlay" else "none")
  }
  if (!is.character(show_delta) || length(show_delta) != 1L || is.na(show_delta)) {
    stop("`show_delta` must be one of \"panel\", \"overlay\", \"none\", TRUE, or FALSE.")
  }
  show_delta <- tolower(show_delta)
  match.arg(show_delta, c("overlay", "panel", "none"))
}

.icTrajectory_y_label <- function(ylab, IC) {
  if (!is.null(ylab)) {
    return(ylab)
  }
  if (!is.null(IC) && length(IC) > 0L && !is.na(IC[1L]) && nzchar(IC[1L])) {
    return(paste0(IC[1L], " score"))
  }
  "IC score"
}

.icTrajectory_plot_style <- function(...) {
  style <- list(...)
  pch_names <- c(
    "accepted_pch", "rejected_pch", "error_pch",
    "baseline_pch", "delta_pch"
  )
  size_names <- c(
    "point_scale", "line_scale", "text_scale",
    "accepted_cex", "rejected_cex", "error_cex", "baseline_cex", "delta_cex",
    "accepted_lwd", "rejected_lwd", "error_lwd", "running_best_lwd",
    "delta_lwd", "zero_lwd",
    "main_cex", "axis_cex", "axis_label_cex", "annotation_cex", "legend_cex"
  )

  for (name in pch_names) {
    style[[name]] <- .icTrajectory_scalar_pch(style[[name]], name)
  }
  for (name in size_names) {
    style[[name]] <- .icTrajectory_nonnegative_scalar(style[[name]], name)
  }
  style$legend <- .icTrajectory_logical_scalar(style$legend, "legend")
  style$legend_position <- .icTrajectory_legend_position(style$legend_position)
  style$legend_inset <- .icTrajectory_legend_inset(style$legend_inset)
  style$legend_bty <- .icTrajectory_legend_bty(style$legend_bty)
  style$legend_labels <- .icTrajectory_legend_labels_arg(style$legend_labels)
  style
}

.icTrajectory_scalar_pch <- function(x, name) {
  valid_numeric <- (is.numeric(x) || is.integer(x)) &&
    length(x) == 1L &&
    is.finite(x)
  valid_character <- is.character(x) &&
    length(x) == 1L &&
    !is.na(x) &&
    nzchar(x)
  if (!isTRUE(valid_numeric) && !isTRUE(valid_character)) {
    stop("`", name, "` must be a single numeric or character plotting symbol.")
  }
  x
}

.icTrajectory_nonnegative_scalar <- function(x, name) {
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x) || x < 0) {
    stop("`", name, "` must be a finite non-negative numeric scalar.")
  }
  as.numeric(x)
}

.icTrajectory_logical_scalar <- function(x, name) {
  if (!is.logical(x) || length(x) != 1L || is.na(x)) {
    stop("`", name, "` must be TRUE or FALSE.")
  }
  x
}

.icTrajectory_legend_position <- function(x) {
  keywords <- c(
    "bottomright", "bottom", "bottomleft", "left", "topleft",
    "top", "topright", "right", "center"
  )
  if (is.character(x) && length(x) == 1L && !is.na(x) && nzchar(x)) {
    if (!(x %in% keywords)) {
      stop(
        "`legend_position` must be one of the base legend keywords ",
        "or a finite numeric vector of length 2."
      )
    }
    return(x)
  }
  if (is.numeric(x) && length(x) == 2L && all(is.finite(x))) {
    return(as.numeric(x))
  }
  stop(
    "`legend_position` must be one of the base legend keywords ",
    "or a finite numeric vector of length 2."
  )
}

.icTrajectory_legend_inset <- function(x) {
  if (!is.numeric(x) || !(length(x) %in% c(1L, 2L)) || any(!is.finite(x))) {
    stop("`legend_inset` must be a finite numeric vector of length 1 or 2.")
  }
  as.numeric(x)
}

.icTrajectory_legend_bty <- function(x) {
  if (!is.character(x) || length(x) != 1L || is.na(x) || !(x %in% c("n", "o"))) {
    stop("`legend_bty` must be \"n\" or \"o\".")
  }
  x
}

.icTrajectory_legend_labels_arg <- function(x) {
  if (is.null(x)) {
    return(NULL)
  }
  if (!is.character(x) || any(is.na(x))) {
    stop("`legend_labels` must be a character vector.")
  }
  label_names <- names(x)
  if (!is.null(label_names)) {
    named <- nzchar(label_names)
    if (any(named) && !all(named)) {
      stop("`legend_labels` must be either fully named or unnamed.")
    }
  }
  x
}

.icTrajectory_legend_labels <- function(labels, keys) {
  defaults <- c(
    running_best = "Running best",
    accepted = "Accepted proposal",
    rejected = "Rejected proposal",
    baseline = "Baseline IC",
    delta = "Delta IC"
  )
  out <- defaults[keys]
  if (is.null(labels)) {
    return(unname(out))
  }

  label_names <- names(labels)
  named_labels <- !is.null(label_names) && any(nzchar(label_names))
  if (isTRUE(named_labels)) {
    valid <- nzchar(label_names)
    invalid <- setdiff(label_names[valid], names(defaults))
    if (length(invalid) > 0L) {
      stop("Named `legend_labels` entries must match known legend keys.")
    }
    displayed <- intersect(label_names[valid], keys)
    out[displayed] <- labels[match(displayed, label_names)]
    return(unname(out))
  }

  if (length(labels) != length(keys)) {
    stop("Unnamed `legend_labels` must match the number of displayed legend entries.")
  }
  labels
}

.icTrajectory_axis_limits <- function(values) {
  limits <- range(values[is.finite(values)])
  pad <- diff(limits) * 0.06
  if (!is.finite(pad) || pad == 0) {
    pad <- if (limits[1L] == 0) 1 else abs(limits[1L]) * 0.06
  }
  range(pretty(limits + c(-pad, pad)))
}

.icTrajectory_status_cols <- function(status, alpha = 1) {
  out <- rep(rgb(0, 0, 0, alpha = alpha), length(status))
  out[status == "accepted"] <- rgb(0.05, 0.26, 0.78, alpha = alpha)
  out[status == "rejected"] <- rgb(0.85, 0.16, 0.08, alpha = alpha)
  out[status == "error"] <- rgb(0.90, 0.48, 0.00, alpha = alpha)
  out
}

.icTrajectory_fillable_pch <- function(pch) {
  is.numeric(pch) && length(pch) == 1L && is.finite(pch) && pch %in% 21:25
}

.icTrajectory_plot_main_panel <- function(x,
                                          x_limits,
                                          y_limits,
                                          x_ticks,
                                          y_ticks,
                                          main,
                                          xlab,
                                          ylab,
                                          draw_x_labels,
                                          include_delta_legend,
                                          style) {
  x_values <- x$step
  y_values <- x$ic
  finite_y <- c(y_values, x$best_ic)
  finite_y <- finite_y[is.finite(finite_y)]

  plot(
    x_values, y_values,
    type = "n",
    xlab = xlab,
    ylab = "",
    xlim = x_limits,
    ylim = y_limits,
    xaxt = "n",
    yaxt = "n",
    cex.lab = style$axis_label_cex,
    bty = "n"
  )
  title(main = main, line = 2, cex.main = style$main_cex)
  axis(
    1,
    at = x_ticks,
    labels = if (isTRUE(draw_x_labels)) x_ticks else FALSE,
    cex.axis = style$axis_cex,
    tck = -0.02
  )
  axis(2, at = y_ticks, labels = y_ticks, las = 1, cex.axis = style$axis_cex, tck = -0.02)
  mtext(ylab, side = 2, line = 3.5, cex = style$axis_label_cex)

  best_ok <- is.finite(x$best_ic)
  if (sum(best_ok) >= 2L) {
    lines(
      x_values[best_ok], x$best_ic[best_ok],
      col = rgb(0.05, 0.26, 0.78, alpha = 0.9),
      type = "s",
      lwd = style$running_best_lwd
    )
  }

  rejected <- x$status == "rejected" & is.finite(y_values)
  if (any(rejected)) {
    points(
      x_values[rejected], y_values[rejected],
      col = rgb(0.85, 0.16, 0.08, alpha = 0.55),
      pch = style$rejected_pch,
      cex = style$rejected_cex,
      lwd = style$rejected_lwd
    )
  }

  accepted <- x$status == "accepted" & is.finite(y_values)
  if (any(accepted)) {
    accepted_col <- if (.icTrajectory_fillable_pch(style$accepted_pch)) {
      "white"
    } else {
      rgb(0.05, 0.26, 0.78, alpha = 0.95)
    }
    points(
      x_values[accepted], y_values[accepted],
      col = accepted_col,
      bg = rgb(0.05, 0.26, 0.78, alpha = 0.95),
      pch = style$accepted_pch,
      cex = style$accepted_cex,
      lwd = style$accepted_lwd
    )
  }

  errored <- x$status == "error" & is.finite(y_values)
  if (any(errored)) {
    points(
      x_values[errored], y_values[errored],
      col = rgb(0.90, 0.48, 0.00, alpha = 0.85),
      pch = style$error_pch,
      cex = style$error_cex,
      lwd = style$error_lwd
    )
  }

  baseline <- x$status == "baseline" & is.finite(y_values)
  if (any(baseline)) {
    baseline_i <- which(baseline)[1L]
    points(
      x_values[baseline_i], y_values[baseline_i],
      col = "black",
      pch = style$baseline_pch,
      cex = style$baseline_cex
    )
    text(
      x_values[baseline_i], y_values[baseline_i],
      labels = paste0(round(y_values[baseline_i], 2)),
      pos = 4, col = "black", cex = style$annotation_cex
    )
  }

  if (any(accepted)) {
    accepted_idx <- which(accepted)
    min_idx <- accepted_idx[which.min(y_values[accepted_idx])]
    text(
      x_values[min_idx],
      y_values[min_idx] - diff(range(finite_y)) * 0.02,
      labels = paste0(round(y_values[min_idx], 2)),
      pos = 1, col = "black", cex = style$annotation_cex
    )
  }

  if (!isTRUE(style$legend)) {
    return(invisible(NULL))
  }

  legend_keys <- c("running_best", "accepted", "rejected", "baseline")
  legend_cols <- c(
    rgb(0.05, 0.26, 0.78, alpha = 0.9),
    if (.icTrajectory_fillable_pch(style$accepted_pch)) "white" else rgb(0.05, 0.26, 0.78, alpha = 0.95),
    rgb(0.85, 0.16, 0.08, alpha = 0.55),
    "black"
  )
  legend_lty <- c(1, NA, NA, NA)
  legend_pch <- c(NA, style$accepted_pch, style$rejected_pch, style$baseline_pch)
  legend_pt_bg <- c(NA, rgb(0.05, 0.26, 0.78, alpha = 0.95), NA, NA)
  legend_pt_lwd <- c(NA, style$accepted_lwd, style$rejected_lwd, 0)

  if (isTRUE(include_delta_legend)) {
    legend_keys <- c(legend_keys, "delta")
    legend_cols <- c(legend_cols, rgb(0, 0, 0, alpha = 0.5))
    legend_lty <- c(legend_lty, 1)
    legend_pch <- c(legend_pch, NA)
    legend_pt_bg <- c(legend_pt_bg, NA)
    legend_pt_lwd <- c(legend_pt_lwd, NA)
  }

  legend_labels <- .icTrajectory_legend_labels(style$legend_labels, legend_keys)

  legend(
    style$legend_position,
    legend = legend_labels,
    col = legend_cols,
    lty = legend_lty,
    pch = legend_pch,
    pt.bg = legend_pt_bg,
    pt.lwd = legend_pt_lwd,
    cex = style$legend_cex,
    bty = style$legend_bty,
    inset = style$legend_inset
  )
}

.icTrajectory_plot_delta_panel <- function(x, x_limits, x_ticks, delta_limits, xlab, style) {
  x_values <- x$step
  delta_limits <- .icTrajectory_delta_limits(x$delta_ic, delta_limits)
  delta_ticks <- pretty(delta_limits)

  plot(
    x_values, rep(NA_real_, length(x_values)),
    type = "n",
    xlab = xlab,
    ylab = "",
    xlim = x_limits,
    ylim = delta_limits,
    xaxt = "n",
    yaxt = "n",
    cex.lab = style$axis_label_cex,
    bty = "n"
  )
  axis(1, at = x_ticks, labels = x_ticks, cex.axis = style$axis_cex, tck = -0.03)
  axis(2, at = delta_ticks, labels = delta_ticks, las = 1, cex.axis = style$axis_cex, tck = -0.02)
  mtext("Delta IC", side = 2, line = 3.5, cex = style$axis_label_cex)
  lines(
    x = c(min(x_limits), max(x_limits)),
    y = c(0, 0),
    col = rgb(0, 0, 0, alpha = 0.45),
    lwd = style$zero_lwd
  )

  delta_ok <- is.finite(x$delta_ic)
  if (!any(delta_ok)) {
    return(invisible(NULL))
  }
  delta_cols <- .icTrajectory_status_cols(x$status[delta_ok], alpha = 0.45)
  segments(
    x0 = x_values[delta_ok],
    y0 = 0,
    x1 = x_values[delta_ok],
    y1 = x$delta_ic[delta_ok],
    col = delta_cols,
    lwd = style$delta_lwd
  )
  points(
    x_values[delta_ok],
    x$delta_ic[delta_ok],
    col = .icTrajectory_status_cols(x$status[delta_ok], alpha = 0.7),
    pch = style$delta_pch,
    cex = style$delta_cex
  )
  invisible(NULL)
}

.icTrajectory_plot_delta_overlay <- function(x, x_limits, x_ticks, delta_limits, style) {
  x_values <- x$step
  delta_limits <- .icTrajectory_delta_limits(x$delta_ic, delta_limits)
  delta_ticks <- pretty(delta_limits)

  plot(
    x_values, rep(NA_real_, length(x_values)),
    type = "n",
    xlab = "", ylab = "",
    xlim = x_limits, ylim = delta_limits,
    xaxt = "n", yaxt = "n", bty = "n"
  )
  lines(
    x = c(min(x_values), max(x_ticks)),
    y = c(0, 0),
    col = rgb(0, 0, 0, alpha = 0.5), lwd = style$zero_lwd
  )

  delta_ok <- is.finite(x$delta_ic)
  lines(
    x_values[delta_ok], x$delta_ic[delta_ok],
    col = rgb(0, 0, 0, alpha = 0.45), lwd = style$delta_lwd
  )
  points(
    x_values[delta_ok], x$delta_ic[delta_ok],
    col = rgb(0, 0, 0, alpha = 0.5),
    pch = style$delta_pch,
    cex = style$delta_cex
  )
  axis(
    4, at = delta_ticks, labels = delta_ticks,
    las = 1, cex.axis = style$axis_cex, tck = -0.02,
    col = rgb(0, 0, 0, alpha = 0.5),
    col.axis = rgb(0, 0, 0, alpha = 0.5)
  )
}
