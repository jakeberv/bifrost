# -------------------------------------------------------------------------
# Internal helpers for bifrost_search printing (not exported)
# -------------------------------------------------------------------------

.bifrost_has <- function(x, nm) !is.null(x[[nm]])

.bifrost_ui <- function(x, nm) {
  ui <- x[["user_input"]]
  if (is.null(ui) || is.null(ui[[nm]])) return(NULL)
  ui[[nm]]
}

.bifrost_fmt_num <- function(v, digits = 2) {
  if (length(v) == 0 || is.null(v) || all(is.na(v))) return("NA")
  if (!is.numeric(v)) return(as.character(v)[1])
  formatC(v[1], format = "f", digits = digits)
}

.bifrost_fmt_int <- function(v) {
  if (length(v) == 0 || is.null(v) || is.na(v[1])) return("NA")
  as.character(as.integer(v[1]))
}

.bifrost_fmt_lgl <- function(v) {
  if (is.null(v) || length(v) == 0) return("NA")
  if (is.logical(v)) return(ifelse(isTRUE(v[1]), "TRUE", "FALSE"))
  as.character(v)[1]
}

.bifrost_fmt_val <- function(v, digits = 2, cutoff = 80) {
  if (is.null(v) || length(v) == 0) return(NULL)
  if (is.character(v)) return(as.character(v)[1])         # no extra quotes
  if (is.logical(v))   return(.bifrost_fmt_lgl(v))
  if (is.numeric(v))   return(.bifrost_fmt_num(v, digits))
  s <- paste(deparse(v, width.cutoff = cutoff), collapse = "")
  if (nchar(s) > cutoff) paste0(substr(s, 1, cutoff - 3), "...")
  else s
}

.bifrost_safe_tip_count <- function(tree) {
  if (is.null(tree)) return(NA_integer_)
  #if (!requireNamespace("ape", quietly = TRUE)) return(NA_integer_)
  tryCatch(ape::Ntip(tree), error = function(e) NA_integer_)
}

.bifrost_safe_regime_count <- function(tree) {
  if (is.null(tree)) return(NA_integer_)
  #if (!requireNamespace("phytools", quietly = TRUE)) return(NA_integer_)
  tryCatch({
    st <- phytools::getStates(tree, type = "both")
    length(unique(st))
  }, error = function(e) NA_integer_)
}

.bifrost_safe_trait_count <- function(model_obj) {
  if (is.null(model_obj)) return(NA_integer_)
  tryCatch({
    if (!is.null(model_obj$Y)) return(as.integer(ncol(as.matrix(model_obj$Y))))
    if (!is.null(model_obj$residuals)) {
      r <- model_obj$residuals
      if (is.matrix(r)) return(as.integer(ncol(r)))
      if (is.numeric(r)) return(1L)
    }
    NA_integer_
  }, error = function(e) NA_integer_)
}

.bifrost_wrap_nodes <- function(nodes, indent = 2L) {
  width <- getOption("width", 80)
  if (is.null(nodes) || length(nodes) == 0) {
    cat(strrep(" ", indent), "none\n", sep = "")
    return(invisible(NULL))
  }
  nodes <- as.integer(nodes)
  per_line <- max(8L, floor((width - indent) / 6L))
  chunks <- split(nodes, ceiling(seq_along(nodes) / per_line))
  for (ch in chunks) {
    cat(strrep(" ", indent), paste(ch, collapse = "  "), "\n", sep = "")
  }
  invisible(NULL)
}

.bifrost_pack_line <- function(indent = 2L, ..., sep = "   ") {
  parts <- list(...)
  parts <- parts[!vapply(parts, is.null, logical(1))]
  parts <- parts[vapply(parts, function(z) nzchar(z), logical(1))]
  if (length(parts) == 0) return(invisible(NULL))

  width <- getOption("width", 80)
  prefix <- strrep(" ", indent)
  cur <- prefix

  for (p in parts) {
    add <- if (identical(cur, prefix)) p else paste0(sep, p)
    if (nchar(cur) + nchar(add) > width && !identical(cur, prefix)) {
      cat(cur, "\n", sep = "")
      cur <- paste0(prefix, p)
    } else {
      cur <- paste0(cur, add)
    }
  }
  cat(cur, "\n", sep = "")
  invisible(NULL)
}

.bifrost_kv <- function(k, v) if (!is.null(v) && nzchar(v)) paste0(k, ": ", v) else NULL

.bifrost_pick_tree <- function(x) {
  if (.bifrost_has(x, "tree_no_uncertainty_untransformed")) {
    x$tree_no_uncertainty_untransformed
  } else if (.bifrost_has(x, "tree_no_uncertainty_transformed")) {
    x$tree_no_uncertainty_transformed
  } else {
    NULL
  }
}

.bifrost_collect_print_data <- function(x) {
  tree <- .bifrost_pick_tree(x)
  model_obj <- if (.bifrost_has(x, "model_no_uncertainty")) x$model_no_uncertainty else NULL

  IC_used <- if (.bifrost_has(x, "IC_used")) as.character(x$IC_used) else "IC"
  baseline_ic <- if (.bifrost_has(x, "baseline_ic")) x$baseline_ic else NA_real_
  optimal_ic  <- if (.bifrost_has(x, "optimal_ic"))  x$optimal_ic  else NA_real_
  delta_ic <- if (is.numeric(baseline_ic) && is.numeric(optimal_ic)) baseline_ic - optimal_ic else NA_real_

  shifts <- x$shift_nodes_no_uncertainty
  n_shifts <- if (is.null(shifts)) 0L else length(shifts)

  candidates <- if (.bifrost_has(x, "num_candidates")) x$num_candidates else NA_integer_
  regimes <- .bifrost_safe_regime_count(tree)
  tips <- .bifrost_safe_tip_count(tree)

  min_desc   <- .bifrost_ui(x, "min_descendant_tips")
  num_cores  <- .bifrost_ui(x, "num_cores")
  thr        <- .bifrost_ui(x, "shift_acceptance_threshold")
  plot       <- .bifrost_ui(x, "plot")
  verbose    <- .bifrost_ui(x, "verbose")
  store_hist <- .bifrost_ui(x, "store_model_fit_history")
  uw_s       <- .bifrost_ui(x, "uncertaintyweights")
  uw_p       <- .bifrost_ui(x, "uncertaintyweights_par")
  wmode <- if (isTRUE(uw_s)) "Serial" else if (isTRUE(uw_p)) "Parallel" else "Off"

  method  <- .bifrost_fmt_val(.bifrost_ui(x, "method"))
  formula <- .bifrost_fmt_val(.bifrost_ui(x, "formula"))
  error   <- .bifrost_fmt_val(.bifrost_ui(x, "error"))
  penalty <- .bifrost_fmt_val(.bifrost_ui(x, "penalty"))
  target  <- .bifrost_fmt_val(.bifrost_ui(x, "target"))

  if (is.null(method) && !is.null(model_obj) && !is.null(model_obj$call$method)) {
    method <- as.character(model_obj$call$method)
  }
  if (is.null(formula) && !is.null(model_obj) && !is.null(model_obj$formula)) {
    formula <- paste(deparse(model_obj$formula), collapse = "")
  }

  model_code <- NULL
  if (!is.null(model_obj) && !is.null(model_obj$call) && !is.null(model_obj$call$model)) {
    model_code <- as.character(model_obj$call$model)
  }
  if (is.null(model_code)) model_code <- if (!is.na(regimes) && regimes > 1) "BMM" else "BM"

  model_label <- if (identical(model_code, "BMM")) "BMM (multi-rate BM)" else "BM (Brownian motion)"
  ntraits <- .bifrost_safe_trait_count(model_obj)

  list(
    tree = tree,
    model_obj = model_obj,
    IC_used = IC_used,
    baseline_ic = baseline_ic,
    optimal_ic = optimal_ic,
    delta_ic = delta_ic,
    shifts = shifts,
    n_shifts = n_shifts,
    candidates = candidates,
    regimes = regimes,
    tips = tips,
    min_desc = min_desc,
    num_cores = num_cores,
    thr = thr,
    plot = plot,
    verbose = verbose,
    store_hist = store_hist,
    wmode = wmode,
    method = method,
    formula = formula,
    error = error,
    penalty = penalty,
    target = target,
    model_label = model_label,
    ntraits = ntraits
  )
}

.bifrost_print_header <- function() {
  ver <- as.character(utils::packageVersion("bifrost"))
  hdr <- paste0("Bifrost Search Result (bifrost ", ver, ")")

  cat(hdr, "\n", sep = "")
  cat(strrep("=", nchar(hdr, type = "width")), "\n\n", sep = "")
}

.bifrost_print_ic_block <- function(d) {
  cat("IC (", d$IC_used, ")\n", sep = "")
  cat(sprintf("  Baseline:  %s   Optimal:  %s   \u0394IC:  %s\n\n",
              .bifrost_fmt_num(d$baseline_ic, 2),
              .bifrost_fmt_num(d$optimal_ic, 2),
              .bifrost_fmt_num(d$delta_ic, 2)))
}

.bifrost_print_search_block <- function(d) {
  cat("Search\n")
  .bifrost_pack_line(
    2,
    .bifrost_kv("Candidates", ifelse(is.na(d$candidates), "NA", .bifrost_fmt_int(d$candidates))),
    .bifrost_kv("Shifts", .bifrost_fmt_int(d$n_shifts)),
    .bifrost_kv("Regimes", ifelse(is.na(d$regimes), "NA", .bifrost_fmt_int(d$regimes))),
    .bifrost_kv("MinTips", if (!is.null(d$min_desc)) .bifrost_fmt_int(d$min_desc) else NULL),
    .bifrost_kv("Cores", if (!is.null(d$num_cores)) .bifrost_fmt_int(d$num_cores) else NULL),
    .bifrost_kv("Threshold", if (!is.null(d$thr)) .bifrost_fmt_num(as.numeric(d$thr), 2) else NULL)
  )
  .bifrost_pack_line(
    2,
    .bifrost_kv("Plot", if (!is.null(d$plot)) .bifrost_fmt_lgl(d$plot) else NULL),
    .bifrost_kv("Verbose", if (!is.null(d$verbose)) .bifrost_fmt_lgl(d$verbose) else NULL),
    .bifrost_kv("FitHistory", if (!is.null(d$store_hist)) .bifrost_fmt_lgl(d$store_hist) else NULL),
    .bifrost_kv("Weights", d$wmode)
  )
  cat("\n")
}

.bifrost_print_mvgls_block <- function(d) {
  cat("mvgls\n")
  .bifrost_pack_line(
    2,
    .bifrost_kv("Model", d$model_label),
    .bifrost_kv("Method", d$method),
    .bifrost_kv("Error", d$error),
    .bifrost_kv("Traits", ifelse(is.na(d$ntraits), "NA", .bifrost_fmt_int(d$ntraits))),
    .bifrost_kv("Tips", ifelse(is.na(d$tips), "NA", .bifrost_fmt_int(d$tips)))
  )
  if (!is.null(d$formula)) {
    cat("  Formula: ", d$formula, "\n", sep = "")
  }
  # NOTE: unconditional call is behavior-preserving (prints nothing if both are absent)
  .bifrost_pack_line(2, .bifrost_kv("Penalty", d$penalty), .bifrost_kv("Target", d$target), sep = "   ")
  cat("\n")
}

.bifrost_print_shift_nodes <- function(d) {
  cat("Shift Nodes\n")
  .bifrost_wrap_nodes(d$shifts, indent = 2L)
  cat("\n")
}

.bifrost_print_history_plot <- function(x, d) {
  if (isTRUE(d$store_hist) &&
      .bifrost_has(x, "model_fit_history") &&
      is.list(x$model_fit_history) &&
      !is.null(x$model_fit_history$ic_acceptance_matrix) &&
      requireNamespace("txtplot", quietly = TRUE)) {

    m <- x$model_fit_history$ic_acceptance_matrix
    if (is.matrix(m) && ncol(m) >= 2) {
      ic_raw  <- suppressWarnings(as.numeric(m[, 1]))
      acc_raw <- m[, 2]

      acc <- if (is.logical(acc_raw)) {
        as.integer(acc_raw)
      } else {
        suppressWarnings(as.integer(as.numeric(acc_raw) == 1))
      }

      ok <- is.finite(ic_raw) & !is.na(acc)
      ic  <- ic_raw[ok]
      acc <- acc[ok]

      if (length(ic) >= 2L && is.finite(d$baseline_ic)) {

        best <- numeric(length(ic))
        cur  <- as.numeric(d$baseline_ic)
        for (i in seq_along(ic)) {
          if (acc[i] == 1L && is.finite(ic[i])) cur <- ic[i]
          best[i] <- cur
        }

        iter  <- 0:length(best)
        best2 <- c(as.numeric(d$baseline_ic), best)

        cat("\nIC History (Best IC by Iteration)\n")

        ylab_txt <- trimws(gsub("\\(.*\\)", "", d$IC_used))
        if (!nzchar(ylab_txt)) ylab_txt <- "IC"

        txt_w <- getOption("bifrost.txtplot.width", 60L)
        txt_h <- getOption("bifrost.txtplot.height", 13L)

        xr <- range(iter)
        yr <- range(best2, finite = TRUE)
        xlim_pad <- c(xr[1] - 0.5, xr[2])

        ydiff <- yr[2] - yr[1]
        ypad  <- if (is.finite(ydiff) && ydiff > 0) 0.10 * ydiff else 1
        ylim_pad <- c(yr[1] - ypad, yr[2] + ypad)

        txtplot::txtplot(
          x = iter,
          y = best2,
          pch = "*",
          width = as.integer(txt_w),
          height = as.integer(txt_h),
          xlab = "Iteration",
          ylab = ylab_txt,
          xlim = xlim_pad,
          ylim = ylim_pad
        )
        cat("\n")
      }
    }
  }
}

.bifrost_print_weights <- function(x, d) {
  if (.bifrost_has(x, "ic_weights") && is.data.frame(x$ic_weights)) {
    title <- paste0(d$IC_used, " Weights (Support)")
    cat(title, "\n")
    cat("  Node   Weight\n")
    cat("  ----   ------\n")

    w <- x$ic_weights
    if (nrow(w) == 0) {
      cat("  (Requested, but no shifts detected)\n\n")
    } else if (!all(c("node", "ic_weight_withshift") %in% names(w))) {
      cat("  (Present, but expected columns missing: node, ic_weight_withshift)\n\n")
    } else {
      ww <- w[, c("node", "ic_weight_withshift"), drop = FALSE]
      ww <- ww[order(ww$ic_weight_withshift, decreasing = TRUE), , drop = FALSE]
      for (i in seq_len(nrow(ww))) {
        cat(sprintf("  %4d   %s\n",
                    as.integer(ww$node[i]),
                    formatC(as.numeric(ww$ic_weight_withshift[i]), format = "f", digits = 6)))
      }
      cat("\n")
    }
  }
}

.bifrost_print_warnings <- function(x) {
  n_warn <- if (.bifrost_has(x, "warnings")) length(x$warnings) else 0L
  cat("Warnings\n")
  cat("  Captured: ", n_warn, "\n", sep = "")
}

# -------------------------------------------------------------------------
#' Print method for bifrost search results
#'
#' Prints a compact summary of a completed Bifrost search, including the baseline and
#' optimal information criterion (IC) values, the inferred shift node set, key search
#' settings, and (when present) optional diagnostics such as IC-history and IC-weight
#' support.
#'
#' @param x A `bifrost_search` object returned by [searchOptimalConfiguration()].
#' @param ... Unused (S3 compatibility).
#'
#' @return Invisibly returns `x`. Called for its printing side effects.
#'
#' @export
print.bifrost_search <- function(x, ...) {
  d <- .bifrost_collect_print_data(x)

  .bifrost_print_header()
  .bifrost_print_ic_block(d)

  # Search summary first, then shift nodes (so "what happened" comes before fit details)
  .bifrost_print_search_block(d)
  .bifrost_print_shift_nodes(d)

  # Fit details next
  .bifrost_print_mvgls_block(d)

  # Optional diagnostics
  .bifrost_print_history_plot(x, d)
  .bifrost_print_weights(x, d)

  # Warnings footer
  .bifrost_print_warnings(x)

  # Citation hint (one-liner)
  cat("\nTo cite this package: citation(\"bifrost\")\n")

  invisible(x)
}
