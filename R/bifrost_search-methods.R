#' Print method for bifrost search results
#'
#' @param x A `bifrost_search` object returned by [searchOptimalConfiguration()].
#' @param ... Unused (S3 compatibility).
#'
#' @export
print.bifrost_search <- function(x, ...) {

  # ---- helpers --------------------------------------------------------------
  .has <- function(nm) !is.null(x[[nm]])

  .ui <- function(nm) {
    ui <- x[["user_input"]]
    if (is.null(ui) || is.null(ui[[nm]])) return(NULL)
    ui[[nm]]
  }

  .fmt_num <- function(v, digits = 2) {
    if (length(v) == 0 || is.null(v) || all(is.na(v))) return("NA")
    if (!is.numeric(v)) return(as.character(v)[1])
    formatC(v[1], format = "f", digits = digits)
  }

  .fmt_int <- function(v) {
    if (length(v) == 0 || is.null(v) || is.na(v[1])) return("NA")
    as.character(as.integer(v[1]))
  }

  .fmt_lgl <- function(v) {
    if (is.null(v) || length(v) == 0) return("NA")
    if (is.logical(v)) return(ifelse(isTRUE(v[1]), "TRUE", "FALSE"))
    as.character(v)[1]
  }

  .fmt_val <- function(v, digits = 2, cutoff = 80) {
    if (is.null(v) || length(v) == 0) return(NULL)
    if (is.character(v)) return(as.character(v)[1])         # no extra quotes
    if (is.logical(v))   return(.fmt_lgl(v))
    if (is.numeric(v))   return(.fmt_num(v, digits))
    s <- paste(deparse(v, width.cutoff = cutoff), collapse = "")
    if (nchar(s) > cutoff) paste0(substr(s, 1, cutoff - 3), "...")
    else s
  }

  .safe_tip_count <- function(tree) {
    if (is.null(tree)) return(NA_integer_)
    if (!requireNamespace("ape", quietly = TRUE)) return(NA_integer_)
    tryCatch(ape::Ntip(tree), error = function(e) NA_integer_)
  }

  .safe_regime_count <- function(tree) {
    if (is.null(tree)) return(NA_integer_)
    if (!requireNamespace("phytools", quietly = TRUE)) return(NA_integer_)
    tryCatch({
      st <- phytools::getStates(tree, type = "both")
      length(unique(st))
    }, error = function(e) NA_integer_)
  }

  .safe_trait_count <- function(model_obj) {
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

  .wrap_nodes <- function(nodes, indent = 2L) {
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

  .pack_line <- function(indent = 2L, ..., sep = "   ") {
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

  .kv <- function(k, v) if (!is.null(v) && nzchar(v)) paste0(k, ": ", v) else NULL

  # ---- pick tree/model ------------------------------------------------------
  tree <- NULL
  if (.has("tree_no_uncertainty_untransformed")) {
    tree <- x$tree_no_uncertainty_untransformed
  } else if (.has("tree_no_uncertainty_transformed")) {
    tree <- x$tree_no_uncertainty_transformed
  }

  model_obj <- if (.has("model_no_uncertainty")) x$model_no_uncertainty else NULL

  # ---- gather core values ---------------------------------------------------
  IC_used <- if (.has("IC_used")) as.character(x$IC_used) else "IC"
  baseline_ic <- if (.has("baseline_ic")) x$baseline_ic else NA_real_
  optimal_ic  <- if (.has("optimal_ic"))  x$optimal_ic  else NA_real_
  delta_ic <- if (is.numeric(baseline_ic) && is.numeric(optimal_ic)) baseline_ic - optimal_ic else NA_real_

  shifts <- x$shift_nodes_no_uncertainty
  n_shifts <- if (is.null(shifts)) 0L else length(shifts)

  candidates <- if (.has("num_candidates")) x$num_candidates else NA_integer_
  regimes <- .safe_regime_count(tree)
  tips <- .safe_tip_count(tree)

  min_desc   <- .ui("min_descendant_tips")
  num_cores  <- .ui("num_cores")
  thr        <- .ui("shift_acceptance_threshold")
  plot       <- .ui("plot")
  verbose    <- .ui("verbose")
  store_hist <- .ui("store_model_fit_history")
  uw_s       <- .ui("uncertaintyweights")
  uw_p       <- .ui("uncertaintyweights_par")
  wmode <- if (isTRUE(uw_s)) "Serial" else if (isTRUE(uw_p)) "Parallel" else "Off"

  # mvgls settings (prefer user_input; fall back to fitted object)
  method  <- .fmt_val(.ui("method"))
  formula <- .fmt_val(.ui("formula"))
  error   <- .fmt_val(.ui("error"))
  penalty <- .fmt_val(.ui("penalty"))
  target  <- .fmt_val(.ui("target"))

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
  ntraits <- .safe_trait_count(model_obj)

  # ---- print ---------------------------------------------------------------
  cat("Bifrost Search Result\n")
  cat(strrep("=", 20), "\n\n", sep = "")

  cat("IC (", IC_used, ")\n", sep = "")
  cat(sprintf("  Baseline:  %s   Optimal:  %s   \u0394IC:  %s\n\n",
              .fmt_num(baseline_ic, 2), .fmt_num(optimal_ic, 2), .fmt_num(delta_ic, 2)))

  cat("Search\n")
  .pack_line(
    2,
    .kv("Candidates", ifelse(is.na(candidates), "NA", .fmt_int(candidates))),
    .kv("Shifts", .fmt_int(n_shifts)),
    .kv("Regimes", ifelse(is.na(regimes), "NA", .fmt_int(regimes))),
    .kv("MinTips", if (!is.null(min_desc)) .fmt_int(min_desc) else NULL),
    .kv("Cores", if (!is.null(num_cores)) .fmt_int(num_cores) else NULL),
    .kv("Threshold", if (!is.null(thr)) .fmt_num(as.numeric(thr), 2) else NULL)
  )
  .pack_line(
    2,
    .kv("Plot", if (!is.null(plot)) .fmt_lgl(plot) else NULL),
    .kv("Verbose", if (!is.null(verbose)) .fmt_lgl(verbose) else NULL),
    .kv("FitHistory", if (!is.null(store_hist)) .fmt_lgl(store_hist) else NULL),
    .kv("Weights", wmode)
  )
  cat("\n")

  cat("mvgls\n")
  .pack_line(
    2,
    .kv("Model", model_label),
    .kv("Method", method),
    .kv("Error", error),
    .kv("Traits", ifelse(is.na(ntraits), "NA", .fmt_int(ntraits))),
    .kv("Tips", ifelse(is.na(tips), "NA", .fmt_int(tips)))
  )
  if (!is.null(formula)) {
    cat("  Formula: ", formula, "\n", sep = "")
  }
  if (!is.null(penalty) || !is.null(target)) {
    .pack_line(2, .kv("Penalty", penalty), .kv("Target", target), sep = "   ")
  }
  cat("\n")

  cat("Shift Nodes\n")
  .wrap_nodes(shifts, indent = 2L)
  cat("\n")

  # ---- IC history plot (txtplot) --------------------------------------------
  if (isTRUE(store_hist) &&
      .has("model_fit_history") &&
      is.list(x$model_fit_history) &&
      !is.null(x$model_fit_history$ic_acceptance_matrix) &&
      requireNamespace("txtplot", quietly = TRUE)) {

    m <- x$model_fit_history$ic_acceptance_matrix
    if (is.matrix(m) && ncol(m) >= 2) {
      ic_raw  <- suppressWarnings(as.numeric(m[, 1]))
      acc_raw <- m[, 2]

      # clearer accept parsing across logical/0-1
      acc <- if (is.logical(acc_raw)) {
        as.integer(acc_raw)
      } else {
        suppressWarnings(as.integer(as.numeric(acc_raw) == 1))
      }

      ok <- is.finite(ic_raw) & !is.na(acc)
      ic  <- ic_raw[ok]
      acc <- acc[ok]

      if (length(ic) >= 2L && is.finite(baseline_ic)) {

        # best-so-far IC across all iterations (updates only when accepted)
        best <- numeric(length(ic))
        cur  <- as.numeric(baseline_ic)
        for (i in seq_along(ic)) {
          if (acc[i] == 1L && is.finite(ic[i])) cur <- ic[i]
          best[i] <- cur
        }

        # Include baseline as iteration 0
        iter  <- 0:length(best)
        best2 <- c(as.numeric(baseline_ic), best)

        cat("\nIC History (Best IC by Iteration)\n")

        # y-axis label: just GIC/BIC (no parenthetical)
        ylab_txt <- trimws(gsub("\\(.*\\)", "", IC_used))
        if (!nzchar(ylab_txt)) ylab_txt <- "IC"

        txt_w <- getOption("bifrost.txtplot.width", 60L)
        txt_h <- getOption("bifrost.txtplot.height", 13L)

        # padding so baseline point isn't glued to top-left corner
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

  # ---- IC weights -----------------------------------------------------------
  if (.has("ic_weights") && is.data.frame(x$ic_weights)) {
    title <- paste0(IC_used, " Weights (Support)")
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

  # ---- warnings -------------------------------------------------------------
  n_warn <- if (.has("warnings")) length(x$warnings) else 0L
  cat("Warnings\n")
  cat("  Captured: ", n_warn, "\n", sep = "")

  invisible(x)
}
