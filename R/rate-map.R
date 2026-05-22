#' Resolve a Color Palette for `rateMap()`
#'
#' Internal helper used by [rateMap()] to normalize different palette
#' specifications into a character vector of colors with a fixed length.
#'
#' @param ncolors Integer number of colors to generate.
#' @param palette Palette specification. This can be a length-1 character
#'   string naming a [grDevices::hcl.colors()] palette, a character vector of
#'   explicit colors, or a function that accepts `ncolors`.
#' @param reverse_palette Logical; if `TRUE`, reverse the resolved palette.
#'
#' @return A character vector of length `ncolors`.
#' @keywords internal
#' @noRd
.rateMap_colors <- function(ncolors, palette, reverse_palette = FALSE) {
  cols <- if (is.function(palette)) {
    palette(ncolors)
  } else if (is.character(palette) && length(palette) == 1L) {
    grDevices::hcl.colors(ncolors, palette = palette)
  } else if (is.character(palette) && length(palette) >= 2L) {
    grDevices::colorRampPalette(palette)(ncolors)
  } else {
    stop("'palette' must be either a palette name, a color vector, or a function.")
  }

  if (!is.character(cols) || length(cols) < 2L) {
    stop("Palette resolution failed to produce at least two colors.")
  }

  if (length(cols) != ncolors) {
    cols <- grDevices::colorRampPalette(cols)(ncolors)
  }

  if (isTRUE(reverse_palette)) {
    cols <- rev(cols)
  }

  cols
}

.rateMap_is_bifrost_search <- function(x) {
  inherits(x, "bifrost_search")
}

.rateMap_extract_tree <- function(fit) {
  if (.rateMap_is_bifrost_search(fit)) {
    tree <- fit$tree_no_uncertainty_untransformed
    if (is.null(tree)) {
      tree <- fit$tree_no_uncertainty_transformed
    }
    return(tree)
  }

  if (inherits(fit, "mvgls")) {
    return(fit$corrSt$phy)
  }

  if (is.list(fit) && inherits(fit$model, "mvgls")) {
    return(fit$model$corrSt$phy)
  }

  if (is.list(fit) && !is.null(fit$variables$tree)) {
    return(fit$variables$tree)
  }

  NULL
}

.rateMap_extract_param <- function(fit) {
  if (.rateMap_is_bifrost_search(fit)) {
    return(fit$model_no_uncertainty$param)
  }

  if (inherits(fit, "mvgls")) {
    return(fit$param)
  }

  if (is.list(fit) && inherits(fit$model, "mvgls")) {
    return(fit$model$param)
  }

  if (is.list(fit) && !is.null(fit$param)) {
    return(fit$param)
  }

  NULL
}

.rateMap_extract_ic <- function(fit) {
  if (.rateMap_is_bifrost_search(fit) ||
      (is.list(fit) && !is.null(fit$optimal_ic))) {
    return(list(
      ic = fit$optimal_ic,
      IC_used = if (is.null(fit$IC_used)) NA_character_ else as.character(fit$IC_used)
    ))
  }

  list(ic = NA_real_, IC_used = NA_character_)
}

.rateMap_mapped_states <- function(tree) {
  unique(unlist(lapply(tree$maps, names), use.names = FALSE))
}

.rateMap_validate_fit <- function(tree, param, index, log) {
  if (!inherits(tree, "phylo") || is.null(tree$maps)) {
    return(paste0(
      "Fit ", index,
      " does not contain a mapped simmap/phylo tree with a '$maps' component."
    ))
  }

  if (!is.numeric(param) || is.null(names(param)) || any(!nzchar(names(param)))) {
    return(paste0("Fit ", index, " does not contain a named numeric parameter vector."))
  }

  if (!all(is.finite(param))) {
    return(paste0("Fit ", index, " contains non-finite rate parameters."))
  }

  if (isTRUE(log) && !all(param > 0)) {
    return(paste0(
      "Fit ", index,
      " contains non-positive rate parameters and cannot be log-transformed."
    ))
  }

  mapped_states <- .rateMap_mapped_states(tree)
  missing_states <- setdiff(mapped_states, names(param))
  if (length(missing_states) > 0L) {
    return(paste0(
      "Fit ", index, " is missing rate parameters for mapped states: ",
      paste(missing_states, collapse = ", ")
    ))
  }

  NULL
}

.rateMap_order_state_names <- function(states) {
  nums <- suppressWarnings(as.numeric(states))
  if (all(is.finite(nums))) {
    states[order(nums)]
  } else {
    sort(states)
  }
}

.rateMap_make_mapped_edge <- function(edge, maps) {
  states <- .rateMap_order_state_names(unique(unlist(lapply(maps, names), use.names = FALSE)))
  out <- matrix(
    0,
    nrow = nrow(edge),
    ncol = length(states),
    dimnames = list(
      edge = apply(edge, 1L, paste, collapse = ","),
      state = states
    )
  )

  for (i in seq_along(maps)) {
    map_i <- maps[[i]]
    sums_i <- rowsum(as.numeric(map_i), group = names(map_i), reorder = FALSE)
    out[i, rownames(sums_i)] <- as.numeric(sums_i[, 1L])
  }

  out
}

.rateMap_default_legend_length <- function(type, heights) {
  if (identical(type, "arc")) {
    max(heights)
  } else {
    0.5 * max(heights)
  }
}

.rateMap_show_legend <- function(legend) {
  !isFALSE(legend) &&
    is.numeric(legend) &&
    length(legend) == 1L &&
    is.finite(legend) &&
    legend > 0
}

.rateMap_getYmult <- function() {
  if (grDevices::dev.cur() == 1L) {
    warning("No graphics device open.")
    return(1)
  }

  xyasp <- graphics::par("pin")
  xycr <- diff(graphics::par("usr"))[c(1L, 3L)]
  xyasp[1L] / xyasp[2L] * xycr[2L] / xycr[1L]
}

.rateMap_with_plotrix_getYmult <- function(expr) {
  global <- .GlobalEnv
  had_getYmult <- exists("getYmult", envir = global, inherits = FALSE)
  old_getYmult <- if (had_getYmult) {
    get("getYmult", envir = global, inherits = FALSE)
  } else {
    NULL
  }

  assign("getYmult", .rateMap_getYmult, envir = global)
  on.exit({
    if (had_getYmult) {
      assign("getYmult", old_getYmult, envir = global)
    } else if (exists("getYmult", envir = global, inherits = FALSE)) {
      rm("getYmult", envir = global)
    }
  }, add = TRUE)

  force(expr)
}

.rateMap_normalize_check <- function(check) {
  if (is.logical(check) && length(check) == 1L && !is.na(check)) {
    return(if (isTRUE(check)) "full" else "none")
  }
  if (is.character(check) && length(check) == 1L) {
    choices <- c("full", "topology", "none")
    if (check %in% choices) {
      return(check)
    }
  }
  stop("'check' must be TRUE, FALSE, 'full', 'topology', or 'none'.")
}

.rateMap_edge_clade_keys <- function(tree) {
  n_tip <- length(tree$tip.label)
  child_map <- split(tree$edge[, 2L], tree$edge[, 1L])
  cache <- new.env(parent = emptyenv())
  separator <- "\r"

  descendants <- function(node) {
    cache_key <- as.character(node)
    if (exists(cache_key, envir = cache, inherits = FALSE)) {
      return(get(cache_key, envir = cache, inherits = FALSE))
    }

    out <- if (node <= n_tip) {
      tree$tip.label[node]
    } else {
      children <- child_map[[cache_key]]
      sort(unlist(lapply(children, descendants), use.names = FALSE))
    }

    assign(cache_key, out, envir = cache)
    out
  }

  vapply(
    tree$edge[, 2L],
    function(node) paste(descendants(node), collapse = separator),
    character(1)
  )
}

.rateMap_mcc_index <- function(trees) {
  clade_keys <- lapply(trees, .rateMap_edge_clade_keys)
  clade_counts <- table(unlist(clade_keys, use.names = FALSE))
  clade_credibility <- clade_counts / length(trees)
  scores <- vapply(
    clade_keys,
    function(keys) sum(log(unname(clade_credibility[keys]))),
    numeric(1)
  )
  which.max(scores)
}

.rateMap_target_tree <- function(trees, target, target_tree = NULL) {
  if (!is.null(target_tree)) {
    if (!inherits(target_tree, "phylo")) {
      stop("'target_tree' must inherit from class 'phylo'.")
    }
    return(target_tree)
  }

  target <- match.arg(target, c("first", "mcc"))
  if (identical(target, "mcc")) {
    trees[[.rateMap_mcc_index(trees)]]
  } else {
    trees[[1L]]
  }
}

.rateMap_match_edges <- function(target_tree, trees) {
  target_keys <- .rateMap_edge_clade_keys(target_tree)
  # nocov start
  if (anyDuplicated(target_keys)) {
    stop("Target tree contains duplicated descendant-tip branch keys.")
  }
  # nocov end

  lapply(seq_along(trees), function(i) {
    tree_keys <- .rateMap_edge_clade_keys(trees[[i]])
    idx <- match(target_keys, tree_keys)
    if (anyNA(idx)) {
      stop(paste0(
        "Tree ", i,
        " is missing one or more target-tree branches; use matching topologies."
      ))
    }
    idx
  })
}

.rateMap_map_value <- function(map, fit_rates, start, end, tol = 1e-10) {
  branch_length <- sum(map)
  if (!is.finite(branch_length) || branch_length <= tol || (end - start) <= tol) {
    return(unname(fit_rates[names(map)[1L]]))
  }

  start <- max(0, min(start, branch_length))
  end <- max(0, min(end, branch_length))
  if ((end - start) <= tol) {
    return(unname(fit_rates[names(map)[1L]]))
  }

  segment_end <- cumsum(as.numeric(map))
  segment_start <- c(0, utils::head(segment_end, -1L))
  overlap <- pmax(0, pmin(segment_end, end) - pmax(segment_start, start))
  total_overlap <- sum(overlap)

  # nocov start
  if (!is.finite(total_overlap) || total_overlap <= tol) {
    return(unname(fit_rates[names(map)[1L]]))
  }
  # nocov end

  rate_values <- unname(fit_rates[names(map)])
  sum(overlap * rate_values) / total_overlap
}

.rateMap_quantile_names <- function(probs) {
  paste0("q", sprintf("%03d", as.integer(round(1000 * probs))))
}

.rateMap_weighted_mean <- function(x, weights) {
  ok <- is.finite(x) & is.finite(weights) & weights > 0
  if (!any(ok)) {
    return(NA_real_)
  }
  x <- x[ok]
  weights <- weights[ok] / sum(weights[ok])
  sum(weights * x)
}

.rateMap_weighted_sd <- function(x, weights) {
  ok <- is.finite(x) & is.finite(weights) & weights > 0
  if (sum(ok) <= 1L) {
    return(NA_real_)
  }
  x <- x[ok]
  weights <- weights[ok] / sum(weights[ok])
  if (all(abs(weights - weights[1L]) < sqrt(.Machine$double.eps))) {
    return(stats::sd(x))
  }
  mu <- sum(weights * x)
  sqrt(sum(weights * (x - mu)^2))
}

.rateMap_weighted_quantile <- function(x, weights, probs) {
  ok <- is.finite(x) & is.finite(weights) & weights > 0
  if (!any(ok)) {
    return(rep(NA_real_, length(probs)))
  }

  x <- x[ok]
  weights <- weights[ok] / sum(weights[ok])

  if (length(x) == 1L) {
    return(rep(x, length(probs)))
  }

  if (all(abs(weights - weights[1L]) < sqrt(.Machine$double.eps))) {
    return(as.numeric(stats::quantile(x, probs = probs, names = FALSE, type = 7)))
  }

  ord <- order(x)
  x <- x[ord]
  weights <- weights[ord]
  cumulative <- cumsum(weights)

  vapply(
    probs,
    function(prob) {
      if (prob <= 0) {
        return(x[1L])
      }
      if (prob >= 1) {
        return(x[length(x)])
      }
      x[which(cumulative >= prob)[1L]]
    },
    numeric(1)
  )
}

.rateMap_hpd <- function(x, weights, prob) {
  ok <- is.finite(x) & is.finite(weights) & weights > 0
  if (!any(ok)) {
    return(c(NA_real_, NA_real_))
  }

  x <- x[ok]
  weights <- weights[ok] / sum(weights[ok])
  if (length(x) == 1L) {
    return(c(x, x))
  }

  ord <- order(x)
  x <- x[ord]
  weights <- weights[ord]

  if (all(abs(weights - weights[1L]) < sqrt(.Machine$double.eps))) {
    width <- max(1L, ceiling(prob * length(x)))
    starts <- seq_len(length(x) - width + 1L)
    ranges <- x[starts + width - 1L] - x[starts]
    best <- starts[which.min(ranges)]
    return(c(x[best], x[best + width - 1L]))
  }

  cumulative <- c(0, cumsum(weights))
  best <- c(x[1L], x[length(x)])
  best_width <- Inf
  tol <- sqrt(.Machine$double.eps)

  for (i in seq_along(x)) {
    target <- cumulative[i] + prob
    if (target > 1 + tol) {
      break
    }
    j <- which(cumulative >= target - tol)[1L] - 1L
    if (!is.na(j) && j >= i) {
      width <- x[j] - x[i]
      if (width < best_width) {
        best_width <- width
        best <- c(x[i], x[j])
      }
    }
  }

  best
}

.rateMap_row_means <- function(values, weights) {
  if (all(is.finite(values))) {
    return(as.numeric(values %*% weights))
  }
  vapply(
    seq_len(nrow(values)),
    function(i) .rateMap_weighted_mean(values[i, ], weights),
    numeric(1)
  )
}

.rateMap_summarize_run_values <- function(values,
                                          weights,
                                          quantile_probs,
                                          hpd_prob) {
  q_names <- .rateMap_quantile_names(quantile_probs)
  quantiles <- t(vapply(
    seq_len(nrow(values)),
    function(i) .rateMap_weighted_quantile(values[i, ], weights, quantile_probs),
    numeric(length(quantile_probs))
  ))
  colnames(quantiles) <- q_names

  hpds <- t(vapply(
    seq_len(nrow(values)),
    function(i) .rateMap_hpd(values[i, ], weights, hpd_prob),
    numeric(2)
  ))
  colnames(hpds) <- c("hpd_low", "hpd_high")

  data.frame(
    mean = .rateMap_row_means(values, weights),
    median = vapply(
      seq_len(nrow(values)),
      function(i) .rateMap_weighted_quantile(values[i, ], weights, 0.5),
      numeric(1)
    ),
    sd = vapply(
      seq_len(nrow(values)),
      function(i) .rateMap_weighted_sd(values[i, ], weights),
      numeric(1)
    ),
    quantiles,
    hpds,
    n = rowSums(is.finite(values)),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}

.rateMap_add_derived_uncertainty <- function(x) {
  if (!all(c("q025", "q975", "hpd_low", "hpd_high", "mean", "sd") %in% names(x))) {
    return(x)
  }

  x$ci_width <- x$q975 - x$q025
  x$hpd_width <- x$hpd_high - x$hpd_low
  x$cv <- ifelse(is.finite(x$mean) & x$mean != 0, x$sd / abs(x$mean), NA_real_)
  x
}

.rateMap_recolor <- function(x, value, palette = NULL, reverse_palette = NULL) {
  if (!is.character(value) || length(value) != 1L || is.na(value) || !nzchar(value)) {
    stop("'value' must be a single non-empty character string.")
  }
  if (!value %in% names(x$intervals)) {
    stop("'", value, "' is not a column in 'x$intervals'.")
  }
  if (!is.numeric(x$intervals[[value]])) {
    stop("'", value, "' must be a numeric column in 'x$intervals'.")
  }

  vals <- lapply(
    seq_along(x$tree$maps),
    function(i) x$intervals[[value]][x$intervals$edge == i]
  )
  lims <- range(unlist(vals, use.names = FALSE), finite = TRUE)
  # nocov start
  if (!all(is.finite(lims))) {
    stop("Could not compute finite plotting limits for '", value, "'.")
  }
  # nocov end
  if (diff(lims) == 0) {
    lims <- lims + c(-0.5, 0.5) * max(abs(lims[1L]), 1) * 1e-6
  }

  ncolors <- length(x$cols)
  if (!is.null(palette) || !is.null(reverse_palette)) {
    selected_palette <- if (is.null(palette)) x$palette else palette
    selected_reverse <- if (is.null(reverse_palette)) {
      isTRUE(x$reverse_palette)
    } else {
      reverse_palette
    }
    x$cols <- .rateMap_colors(ncolors, selected_palette, selected_reverse)
    names(x$cols) <- as.character(seq_len(ncolors))
    x$palette <- selected_palette
    x$reverse_palette <- selected_reverse
  }

  breaks <- seq(lims[1L], lims[2L], length.out = ncolors + 1L)
  bins_by_edge <- vector("list", length(x$tree$maps))
  for (i in seq_along(x$tree$maps)) {
    bins <- findInterval(vals[[i]], breaks, all.inside = TRUE)
    bins_by_edge[[i]] <- bins
    names(x$tree$maps[[i]]) <- as.character(bins)
  }

  x$tree$mapped.edge <- .rateMap_make_mapped_edge(x$tree$edge, x$tree$maps)
  x$values <- vals
  x$lims <- lims
  x$breaks <- breaks
  x$intervals$value <- x$intervals[[value]]
  x$intervals$color_bin <- unlist(bins_by_edge, use.names = FALSE)
  x$plot_value <- value
  x$title <- switch(
    value,
    value = x$title,
    mean = if (isTRUE(x$log)) "Mean log fitted rate" else "Mean fitted rate",
    median = if (isTRUE(x$log)) "Median log fitted rate" else "Median fitted rate",
    sd = if (isTRUE(x$log)) "SD log fitted rate" else "SD fitted rate",
    ci_width = if (isTRUE(x$log)) "95% quantile width (log rate)" else "95% quantile width",
    hpd_width = if (isTRUE(x$log)) "95% HPD width (log rate)" else "95% HPD width",
    cv = "Coefficient of variation",
    value
  )
  x
}

.rateMap_normalize_weights <- function(weights) {
  if (!is.numeric(weights) || length(weights) == 0L || any(!is.finite(weights))) {
    stop("Fit weights must be a non-empty finite numeric vector.")
  }
  if (any(weights < 0)) {
    stop("Fit weights must be non-negative.")
  }
  weight_sum <- sum(weights)
  if (!is.finite(weight_sum) || weight_sum <= 0) {
    stop("At least one fit weight must be positive.")
  }
  as.numeric(weights) / weight_sum
}

.rateMap_subset_user_weights <- function(weights, n_fits, keep) {
  if (length(weights) == n_fits) {
    return(weights[keep])
  }
  if (length(weights) == length(keep)) {
    return(weights)
  }
  stop("Fit weights must have length equal to the input fits or the retained fits.")
}

.rateMap_resolve_weights <- function(fits,
                                     keep,
                                     weights,
                                     fit_weights = NULL) {
  n_fits <- length(fits)
  n_used <- length(keep)

  if (!is.null(fit_weights)) {
    resolved <- .rateMap_normalize_weights(
      .rateMap_subset_user_weights(fit_weights, n_fits, keep)
    )
    return(list(
      weights = resolved,
      mode = "custom",
      ic = rep(NA_real_, n_used),
      IC_used = rep(NA_character_, n_used)
    ))
  }

  if (is.numeric(weights)) {
    resolved <- .rateMap_normalize_weights(
      .rateMap_subset_user_weights(weights, n_fits, keep)
    )
    return(list(
      weights = resolved,
      mode = "custom",
      ic = rep(NA_real_, n_used),
      IC_used = rep(NA_character_, n_used)
    ))
  }

  weights <- match.arg(weights, c("equal", "ic"))
  if (identical(weights, "equal")) {
    return(list(
      weights = rep(1 / n_used, n_used),
      mode = "equal",
      ic = rep(NA_real_, n_used),
      IC_used = rep(NA_character_, n_used)
    ))
  }

  ic_info <- lapply(fits[keep], .rateMap_extract_ic)
  ic_values <- vapply(ic_info, function(x) x$ic[1L], numeric(1))
  ic_family <- vapply(ic_info, function(x) x$IC_used[1L], character(1))

  if (!all(is.finite(ic_values))) {
    stop("weights = 'ic' requires every retained fit to contain a finite 'optimal_ic'.")
  }

  present_family <- unique(ic_family[!is.na(ic_family) & nzchar(ic_family)])
  if (length(present_family) != 1L || any(is.na(ic_family) | !nzchar(ic_family))) {
    stop("weights = 'ic' requires every retained fit to contain the same non-missing 'IC_used'.")
  }

  delta_ic <- ic_values - min(ic_values)
  raw_weights <- exp(-0.5 * delta_ic)
  resolved <- .rateMap_normalize_weights(raw_weights)

  list(
    weights = resolved,
    mode = "ic",
    ic = ic_values,
    IC_used = ic_family
  )
}

#' Summarize Branchwise Rate Variation Across Multiple Runs
#'
#' Summarize fitted regime-specific rates across a list of stochastic-map-aware
#' model fits or completed `bifrost_search` results. By default, `rateMap()`
#' preserves the original same-tree workflow: the first retained tree is used as
#' the plotting scaffold, all retained trees must match in topology and branch
#' lengths, branches are sliced on a global depth grid, and runs are averaged
#' with equal weight.
#'
#' @details
#' **Algorithmic provenance.** `rateMap()` is inspired by
#' [phytools::densityMap()], which summarizes a set of stochastic maps by
#' slicing mapped branches on a shared depth grid and coloring a SIMMAP-style
#' tree by an averaged branchwise quantity. Here the averaged quantity is not a
#' posterior probability of a mapped state; instead, each mapped state is
#' translated through the corresponding fitted regime-rate parameter from each
#' run. The plotting interface likewise follows the `phytools` density-map
#' family by returning a colored SIMMAP tree and drawing it with
#' [phytools::plotSimmap()] plus a continuous color-bar legend.
#'
#' `rateMap()` can also summarize same-topology posterior or sensitivity samples
#' where branch lengths differ. In that case, supply `check = "topology"` and,
#' usually, an explicit `target_tree`. This can be any target or summary tree
#' with the same topology and tip labels as the inputs; an MCC tree is only one
#' possible choice. The target tree supplies the plotted topology and branch
#' lengths, while rates are matched from each input tree by descendant-tip clade
#' keys rather than by edge order.
#'
#' **Summary modes.** With `summary = "interval"`, each target-tree branch is
#' subdivided by the global depth grid controlled by `res`. If source branch
#' lengths differ from the target branch length, source stochastic-map segments
#' are projected onto target intervals by relative position along the matched
#' branch. With `summary = "branch"`, each edge receives one length-weighted
#' average rate from each run before the across-run summary is computed.
#'
#' **Tree checks and targets.** `check = TRUE` is equivalent to `check = "full"`
#' and requires topology and branch lengths to match the target tree. Use
#' `check = "topology"` when all inputs have the same topology and tip labels but
#' may have different branch lengths. `check = FALSE` or `check = "none"` skips
#' the upfront `ape::all.equal.phylo()` check, but every target branch must still
#' be recoverable in every input tree by descendant-tip set. This function is not
#' a mixed-topology posterior summarizer; clades absent from an input tree are
#' treated as an error rather than being marginalized over topology.
#'
#' When `target_tree` is supplied, it is used directly as the plotting scaffold
#' and need not contain SIMMAP maps. When `target_tree = NULL`, `target = "first"`
#' uses the first retained input tree, and `target = "mcc"` chooses the retained
#' input tree with the highest sum of log clade credibilities. For truly
#' same-topology inputs, the MCC score is usually tied, so `"mcc"` commonly
#' resolves to the first retained tree unless topologies differ. If you already
#' have a preferred consensus, chronogram, MCC, maximum-likelihood, or otherwise
#' curated target tree, pass it with `target_tree` instead of using `target`.
#'
#' **Fit weights.** `weights = "equal"` assigns the same weight to each retained
#' fit. `weights = "ic"` computes standard information-criterion weights from
#' each retained fit's `optimal_ic`, requiring all retained fits to share the
#' same non-missing `IC_used`. Numeric `weights` or `fit_weights` can be supplied
#' for custom weighting. Weights are subset to retained fits after
#' `na_action = "omit"` and are normalized to sum to one.
#'
#' **Uncertainty summaries.** When `uncertainty = TRUE`, the returned
#' `intervals` data frame includes across-fit summaries for each branch or
#' interval: weighted mean, weighted median, weighted standard deviation,
#' quantiles, HPD interval, quantile/HPD widths, coefficient of variation, and
#' the number of finite run-level values. The run-level values are also returned
#' in `run_values` as one matrix per edge. For posterior tree samples, use
#' `weights = "equal"` so the quantiles and HPDs are empirical posterior
#' summaries. `value_summary` controls whether the plotted `value` column uses
#' the weighted mean or weighted median.
#'
#' @param fits A non-empty list of completed runs or fitted model objects.
#'   Supported shapes are `bifrost_search` objects, `mvgls` objects,
#'   `list(model = <mvgls>)`, and scratch-style lists with
#'   `variables$tree` plus `param`.
#' @param res Integer resolution of the global depth grid used to subdivide
#'   target-tree branches when `summary = "interval"`. Ignored by
#'   `summary = "branch"`.
#' @param fsize Optional plotting font sizes passed to [plotRateMap()] when
#'   `plot = TRUE`.
#' @param ftype Optional plotting font type(s) passed to [plotRateMap()] when
#'   `plot = TRUE`.
#' @param lwd Line width passed through to [plotRateMap()] when `plot = TRUE`.
#' @param check Logical or character check mode. `TRUE` or `"full"` verifies
#'   that extracted trees and the target tree match in topology and branch
#'   lengths. `"topology"` verifies matching topology/tip labels while allowing
#'   branch lengths to differ. `FALSE` or `"none"` skips this check, but target
#'   branches still must be matchable by descendant-tip sets.
#' @param legend Legend length passed through to [plotRateMap()] when
#'   `plot = TRUE`.
#' @param outline Logical; if `TRUE`, draw branch outlines in the plotted map
#'   when `plot = TRUE`.
#' @param type Plot type to use when `plot = TRUE`. Supported values are
#'   `"phylogram"`, `"fan"`, and `"arc"`.
#' @param direction Plotting direction for `type = "phylogram"` when
#'   `plot = TRUE`.
#' @param plot Logical; if `TRUE`, plot the result immediately using
#'   [plotRateMap()]. If `FALSE`, only return the computed `"rateMap"` object.
#' @param tree_fun Optional function used to extract a mapped tree from each
#'   element of `fits`. If `NULL`, common Bifrost and mvgls shapes are
#'   auto-detected.
#' @param param_fun Optional function used to extract a named numeric vector of
#'   state-specific fitted rates from each element of `fits`. If `NULL`, common
#'   Bifrost and mvgls shapes are auto-detected.
#' @param tip_fsize Optional override for phylogeny tip-label font size when
#'   `plot = TRUE`.
#' @param legend_fsize Optional override for legend font size when `plot = TRUE`.
#' @param show_tip_labels Logical; if `FALSE`, suppress tip labels in the plotted
#'   map when `plot = TRUE`.
#' @param workers Optional number of `future` workers to use. If `NULL`, the
#'   current `future::plan()` is used as-is.
#' @param future_strategy Future backend used only when `workers` is supplied.
#'   Must be one of `"multisession"` or `"multicore"`.
#' @param future_seed Seed control passed to [future.apply::future_lapply()].
#' @param future_scheduling Scheduling control passed to
#'   [future.apply::future_lapply()].
#' @param future_chunk_size Chunk size passed to
#'   [future.apply::future_lapply()].
#' @param progress Logical; if `TRUE`, display a text progress bar via
#'   [progressr::with_progress()].
#' @param log Logical; if `TRUE`, transform extracted rate parameters with
#'   `log()` before averaging.
#' @param ncolors Number of colors in the rendered palette.
#' @param palette Palette specification. This can be an `hcl.colors()` palette
#'   name, a vector of colors, or a palette function.
#' @param reverse_palette Logical; if `TRUE`, reverse the resolved palette.
#' @param legend_title Optional legend title. If `NULL`, defaults to
#'   `"Mean fitted rate"` or `"Mean log fitted rate"` depending on `log`.
#' @param na_action What to do when a run has invalid parameters. `"error"`
#'   stops immediately. `"omit"` drops invalid runs before aggregation.
#' @param summary Character; `"interval"` slices branches on a global depth grid,
#'   while `"branch"` computes one length-weighted value per edge.
#' @param target Character target-tree selection when `target_tree = NULL`.
#'   `"first"` uses the first retained input tree. `"mcc"` chooses the retained
#'   input tree with the highest sum of log clade credibilities. For same-topology
#'   inputs this is usually tied and therefore returns the first retained tree.
#' @param target_tree Optional explicit target tree used as the plotting scaffold.
#'   This may be any target or summary tree with the same topology and tip labels
#'   as the inputs. It does not need to contain stochastic maps because
#'   `rateMap()` replaces maps with color-bin maps in the returned object. Branch
#'   lengths in `target_tree` define the geometry of the returned and plotted
#'   tree.
#' @param weights Fit-level weighting mode. `"equal"` gives every retained fit
#'   equal weight. `"ic"` computes standard IC weights from `optimal_ic` and
#'   requires all retained fits to have the same `IC_used`. A numeric vector is
#'   also accepted and treated as custom fit weights.
#' @param fit_weights Optional numeric custom fit weights. If supplied, these
#'   override `weights`.
#' @param uncertainty Logical; if `TRUE`, compute and return across-fit
#'   uncertainty summaries for every plotted branch interval or whole branch.
#' @param value_summary Character; central estimate stored in `intervals$value`
#'   and mapped to colors. `"mean"` preserves the legacy weighted-mean behavior.
#'   `"median"` uses the weighted median.
#' @param quantile_probs Numeric length-2 vector of quantile probabilities to
#'   report when `uncertainty = TRUE`.
#' @param hpd_prob Numeric scalar giving the HPD mass to report when
#'   `uncertainty = TRUE`.
#' @param ... Additional plotting arguments passed to [plotRateMap()] when
#'   `plot = TRUE`, such as `mar`, `offset`, `xlim`, `ylim`, `hold`,
#'   `underscore`, or `arc_height`.
#'
#' @return An object of class `"rateMap"` with components:
#' \describe{
#'   \item{`tree`}{A SIMMAP-style tree whose mapped segments encode color-bin
#'   indices.}
#'   \item{`cols`}{The resolved color palette.}
#'   \item{`lims`}{Numeric length-2 vector giving the plotted value range.}
#'   \item{`breaks`}{Numeric vector of palette bin boundaries.}
#'   \item{`values`}{List of branch-interval central values before color
#'   binning.}
#'   \item{`intervals`}{Data frame with one row per plotted branch interval, or
#'   one row per branch when `summary = "branch"`.}
#'   \item{`run_values`}{When `uncertainty = TRUE`, list of numeric matrices
#'   containing run-level values for each edge. Rows are plotted intervals and
#'   columns are retained fits. Otherwise `NULL`.}
#'   \item{`clade_key`}{Character descendant-tip key for each target-tree edge.}
#'   \item{`edge_matches`}{Integer matrix mapping target-tree edge rows to
#'   matched source-tree edge rows for each retained fit.}
#'   \item{`summary`}{The summary mode used, `"interval"` or `"branch"`.}
#'   \item{`uncertainty`}{Logical indicating whether uncertainty summaries were
#'   computed.}
#'   \item{`value_summary`}{Central estimate used for `intervals$value`.}
#'   \item{`quantile_probs`}{Quantile probabilities used for uncertainty
#'   summaries.}
#'   \item{`hpd_prob`}{HPD mass used for uncertainty summaries.}
#'   \item{`plot_value`}{Current interval column mapped to branch colors.}
#'   \item{`target`}{Target-tree selection mode used.}
#'   \item{`check`}{Tree compatibility check mode used.}
#'   \item{`weights`}{Normalized fit weights used for aggregation.}
#'   \item{`weight_mode`}{Weighting mode used: `"equal"`, `"ic"`, or
#'   `"custom"`.}
#'   \item{`weight_table`}{Data frame linking retained input indices, weights,
#'   and IC values when available.}
#'   \item{`palette`}{Original palette specification.}
#'   \item{`reverse_palette`}{Logical indicating whether the palette was
#'   reversed.}
#'   \item{`title`}{Legend title used for plotting.}
#'   \item{`n_fits`}{Number of fits used after validation or omission.}
#'   \item{`omitted`}{Integer indices of omitted fits when `na_action = "omit"`.}
#' }
#'
#' @seealso [plotRateMap()], [phytools::densityMap()],
#'   [phytools::plotSimmap()]
#'
#' @examples
#' \dontrun{
#' # A list of completed Bifrost searches can be summarized directly:
#' rm_obj <- rateMap(list(search_a, search_b, search_c), plot = FALSE)
#' plotRateMap(rm_obj, type = "arc", show_tip_labels = FALSE)
#'
#' # Scratch-style lists remain supported:
#' scratch_fit <- list(
#'   variables = list(tree = mapped_tree),
#'   param = c("0" = 0.12, "1" = 0.45)
#' )
#' rateMap(list(scratch_fit), plot = FALSE)
#'
#' # Use Bifrost model-level IC weights across comparable sensitivity runs:
#' rm_ic <- rateMap(list(search_a, search_b), weights = "ic", plot = FALSE)
#'
#' # Whole-branch summaries avoid interval subdivision:
#' branch_rates <- rateMap(list(search_a, search_b), summary = "branch", plot = FALSE)
#'
#' # Custom fit weights are normalized internally:
#' custom_weighted <- rateMap(
#'   list(search_a, search_b, search_c),
#'   weights = c(2, 1, 1),
#'   summary = "branch",
#'   plot = FALSE
#' )
#'
#' # Posterior trees with the same topology but different branch lengths can be
#' # summarized on any explicit target or summary tree:
#' posterior_target_rates <- rateMap(
#'   posterior_fit_list,
#'   check = "topology",
#'   target_tree = summary_tree,
#'   summary = "branch",
#'   weights = "equal",
#'   uncertainty = TRUE,
#'   plot = FALSE
#' )
#' plotRateMap(posterior_target_rates, value = "sd", type = "arc")
#'
#' # Or choose the retained input tree with the highest summed log clade
#' # credibility as the plotting scaffold:
#' posterior_mcc_rates <- rateMap(
#'   posterior_fit_list,
#'   check = "topology",
#'   target = "mcc",
#'   summary = "branch",
#'   plot = FALSE
#' )
#'
#' # Large sensitivity sets can be explicitly subsampled before plotting:
#' set.seed(1)
#' idx <- sample(seq_along(fit_list), size = 1000)
#' rm_sub <- rateMap(
#'   fit_list[idx],
#'   res = 100,
#'   workers = 8,
#'   future_strategy = "multisession",
#'   log = TRUE,
#'   palette = c("lightblue", "blue", "pink", "red"),
#'   plot = FALSE
#' )
#' plotRateMap(rm_sub, type = "arc", show_tip_labels = FALSE, lwd = 1)
#' }
#'
#' @export
rateMap <- function(
  fits,
  res = 100,
  fsize = NULL,
  ftype = NULL,
  lwd = 3,
  check = TRUE,
  legend = NULL,
  outline = FALSE,
  type = "phylogram",
  direction = "rightwards",
  plot = TRUE,
  tree_fun = NULL,
  param_fun = NULL,
  tip_fsize = NULL,
  legend_fsize = NULL,
  show_tip_labels = TRUE,
  workers = NULL,
  future_strategy = c("multisession", "multicore"),
  future_seed = FALSE,
  future_scheduling = 1,
  future_chunk_size = NULL,
  progress = TRUE,
  log = FALSE,
  ncolors = 256,
  palette = "YlOrRd",
  reverse_palette = FALSE,
  legend_title = NULL,
  na_action = c("error", "omit"),
  summary = c("interval", "branch"),
  target = c("first", "mcc"),
  target_tree = NULL,
  weights = c("equal", "ic"),
  fit_weights = NULL,
  uncertainty = FALSE,
  value_summary = c("mean", "median"),
  quantile_probs = c(0.025, 0.975),
  hpd_prob = 0.95,
  ...
) {
  # nocov start
  if (!requireNamespace("phytools", quietly = TRUE)) {
    stop("Package 'phytools' is required.")
  }
  if (!requireNamespace("ape", quietly = TRUE)) {
    stop("Package 'ape' is required.")
  }
  if (!requireNamespace("future", quietly = TRUE)) {
    stop("Package 'future' is required.")
  }
  if (!requireNamespace("future.apply", quietly = TRUE)) {
    stop("Package 'future.apply' is required.")
  }
  if (!requireNamespace("progressr", quietly = TRUE)) {
    stop("Package 'progressr' is required.")
  }
  # nocov end
  if (!is.list(fits) || length(fits) == 0L) {
    stop("'fits' must be a non-empty list of fitted run objects.")
  }
  if (!is.numeric(res) || length(res) != 1L || !is.finite(res) || res < 1) {
    stop("'res' must be a finite number >= 1.")
  }
  if (!is.numeric(ncolors) || length(ncolors) != 1L || !is.finite(ncolors) || ncolors < 2) {
    stop("'ncolors' must be a finite number >= 2.")
  }
  if (!is.null(workers) && (!is.numeric(workers) || length(workers) != 1L || workers < 1)) {
    stop("'workers' must be NULL or a positive number.")
  }
  if (!is.logical(uncertainty) || length(uncertainty) != 1L || is.na(uncertainty)) {
    stop("'uncertainty' must be TRUE or FALSE.")
  }
  if (!is.numeric(quantile_probs) ||
      length(quantile_probs) != 2L ||
      any(!is.finite(quantile_probs)) ||
      any(quantile_probs < 0 | quantile_probs > 1) ||
      quantile_probs[1L] >= quantile_probs[2L]) {
    stop("'quantile_probs' must be an increasing finite numeric vector of length 2 between 0 and 1.")
  }
  if (!is.numeric(hpd_prob) ||
      length(hpd_prob) != 1L ||
      !is.finite(hpd_prob) ||
      hpd_prob <= 0 ||
      hpd_prob > 1) {
    stop("'hpd_prob' must be a finite number in (0, 1].")
  }

  res <- as.integer(res)
  ncolors <- as.integer(ncolors)
  if (!is.null(workers)) {
    workers <- as.integer(workers)
  }
  future_strategy <- match.arg(future_strategy)
  na_action <- match.arg(na_action)
  summary <- match.arg(summary)
  target <- match.arg(target)
  check_mode <- .rateMap_normalize_check(check)
  value_summary <- match.arg(value_summary)
  type <- match.arg(type, c("phylogram", "fan", "arc"))
  target_mode <- if (is.null(target_tree)) target else "user"

  if (is.null(legend_title)) {
    legend_title <- if (isTRUE(log)) "Mean log fitted rate" else "Mean fitted rate"
  }

  tree_fun <- if (is.null(tree_fun)) .rateMap_extract_tree else tree_fun
  param_fun <- if (is.null(param_fun)) .rateMap_extract_param else param_fun

  trees <- lapply(fits, tree_fun)
  params <- lapply(fits, param_fun)

  validation <- vapply(
    seq_along(fits),
    function(i) {
      msg <- .rateMap_validate_fit(trees[[i]], params[[i]], i, log = log)
      if (is.null(msg)) "" else msg
    },
    character(1)
  )
  invalid <- which(nzchar(validation))

  omitted <- integer()
  keep <- seq_along(fits)
  if (length(invalid) > 0L) {
    if (identical(na_action, "error")) {
      stop(validation[invalid[1L]])
    }
    omitted <- invalid
    keep <- setdiff(seq_along(fits), invalid)
    if (length(keep) == 0L) {
      stop("No fits remain after omitting invalid rate-map inputs.")
    }
    trees <- trees[keep]
    params <- params[keep]
  }

  resolved_weights <- .rateMap_resolve_weights(
    fits = fits,
    keep = keep,
    weights = weights,
    fit_weights = fit_weights
  )
  fit_weights_used <- resolved_weights$weights

  if (isTRUE(log)) {
    params <- lapply(params, base::log)
  }

  target_tree <- .rateMap_target_tree(
    trees = trees,
    target = target,
    target_tree = target_tree
  )

  if (identical(check_mode, "full")) {
    ref_tree <- ape::as.phylo(target_tree)
    same_tree <- vapply(
      trees,
      function(tr) isTRUE(ape::all.equal.phylo(ref_tree, ape::as.phylo(tr))),
      logical(1)
    )
    if (!all(same_tree)) {
      stop("Some trees do not match in topology or branch lengths.")
    }
  } else if (identical(check_mode, "topology")) {
    ref_tree <- ape::as.phylo(target_tree)
    same_tree <- vapply(
      trees,
      function(tr) {
        isTRUE(ape::all.equal.phylo(
          ref_tree,
          ape::as.phylo(tr),
          use.edge.length = FALSE
        ))
      },
      logical(1)
    )
    if (!all(same_tree)) {
      stop("Some trees do not match in topology.")
    }
  }

  trees <- lapply(trees, function(tr) {
    if (!inherits(tr, "simmap")) {
      class(tr) <- c("simmap", class(tr))
    }
    tr
  })
  class(trees) <- c("multiSimmap", "multiPhylo")
  target_tree <- ape::as.phylo(target_tree)
  edge_matches <- .rateMap_match_edges(target_tree, trees)
  target_clade_keys <- .rateMap_edge_clade_keys(target_tree)

  heights <- phytools::nodeHeights(target_tree)
  max_height <- max(heights)
  steps <- if (identical(summary, "interval")) {
    (0:res / res) * max_height
  } else {
    numeric()
  }
  tree <- target_tree
  H <- heights
  tol <- 1e-10

  edge_indices <- seq_len(nrow(tree$edge))

  compute_edge <- function(i) {
    mids <- if (identical(summary, "interval")) {
      intersect(which(steps > H[i, 1L]), which(steps < H[i, 2L]))
    } else {
      integer()
    }
    yy <- cbind(
      c(H[i, 1L], steps[mids]),
      c(steps[mids], H[i, 2L])
    ) - H[i, 1L]
    lengths_i <- yy[, 2L] - yy[, 1L]

    values_by_fit <- matrix(
      NA_real_,
      nrow = nrow(yy),
      ncol = length(trees),
      dimnames = list(NULL, paste0("fit_", keep))
    )
    target_length <- H[i, 2L] - H[i, 1L]

    for (j in seq_along(trees)) {
      source_edge <- edge_matches[[j]][i]
      map_i <- trees[[j]]$maps[[source_edge]]
      source_length <- sum(map_i)
      fit_rates <- params[[j]]

      for (k in seq_len(nrow(yy))) {
        if (target_length > tol) {
          source_start <- (yy[k, 1L] / target_length) * source_length
          source_end <- (yy[k, 2L] / target_length) * source_length
        } else {
          source_start <- 0
          source_end <- source_length
        }

        value_k <- .rateMap_map_value(
          map = map_i,
          fit_rates = fit_rates,
          start = source_start,
          end = source_end,
          tol = tol
        )
        values_by_fit[k, j] <- value_k
      }
    }

    uncertainty_i <- if (isTRUE(uncertainty) || identical(value_summary, "median")) {
      .rateMap_summarize_run_values(
        values = values_by_fit,
        weights = fit_weights_used,
        quantile_probs = quantile_probs,
        hpd_prob = hpd_prob
      )
    } else {
      NULL
    }
    values_i <- if (identical(value_summary, "median")) {
      uncertainty_i$median
    } else {
      .rateMap_row_means(values_by_fit, fit_weights_used)
    }

    if (isTRUE(progress)) {
      p(sprintf("edge %d", i))
    }

    list(
      values = values_i,
      run_values = if (isTRUE(uncertainty)) values_by_fit else NULL,
      uncertainty = if (isTRUE(uncertainty)) uncertainty_i else NULL,
      lengths = lengths_i,
      depth_start = H[i, 1L] + yy[, 1L],
      depth_end = H[i, 1L] + yy[, 2L]
    )
  }

  if (!exists("p", inherits = FALSE)) {
    p <- NULL
  }

  if (!is.null(workers)) {
    old_plan <- future::plan()
    on.exit(future::plan(old_plan), add = TRUE)
    strategy_fun <- get(future_strategy, envir = asNamespace("future"), mode = "function")
    future::plan(strategy_fun, workers = workers)
  }

  compute_all_edges <- function() {
    future.apply::future_lapply(
      edge_indices,
      compute_edge,
      future.seed = future_seed,
      future.scheduling = future_scheduling,
      future.chunk.size = future_chunk_size
    )
  }

  edge_results <- if (isTRUE(progress)) {
    progressr::with_progress(
      {
        p <- progressr::progressor(along = edge_indices)
        compute_all_edges()
      },
      handlers = progressr::handler_txtprogressbar(),
      enable = TRUE
    )
  } else {
    compute_all_edges()
  }

  interval_values <- lapply(edge_results, `[[`, "values")
  interval_lengths <- lapply(edge_results, `[[`, "lengths")
  run_values <- if (isTRUE(uncertainty)) {
    lapply(edge_results, `[[`, "run_values")
  } else {
    NULL
  }

  lims <- range(unlist(interval_values, use.names = FALSE), finite = TRUE)
  # nocov start
  if (!all(is.finite(lims))) {
    stop("Could not compute finite rate limits from the supplied fits.")
  }
  # nocov end

  if (diff(lims) == 0) {
    lims <- lims + c(-0.5, 0.5) * max(abs(lims[1L]), 1) * 1e-6
  }

  breaks <- seq(lims[1L], lims[2L], length.out = ncolors + 1L)
  cols <- .rateMap_colors(
    ncolors = ncolors,
    palette = palette,
    reverse_palette = reverse_palette
  )
  names(cols) <- as.character(seq_len(ncolors))

  tree$maps <- vector("list", nrow(tree$edge))
  bins_by_edge <- vector("list", nrow(tree$edge))
  for (i in seq_len(nrow(tree$edge))) {
    bins <- findInterval(interval_values[[i]], breaks, all.inside = TRUE)
    bins_by_edge[[i]] <- bins
    tree$maps[[i]] <- interval_lengths[[i]]
    names(tree$maps[[i]]) <- as.character(bins)
  }

  tree$mapped.edge <- .rateMap_make_mapped_edge(tree$edge, tree$maps)
  class(tree) <- c("simmap", setdiff(class(tree), "simmap"))
  attr(tree, "map.order") <- "right-to-left"

  intervals <- do.call(
    rbind,
    lapply(edge_indices, function(i) {
      base_row <- data.frame(
        edge = i,
        parent = tree$edge[i, 1L],
        child = tree$edge[i, 2L],
        depth_start = edge_results[[i]]$depth_start,
        depth_end = edge_results[[i]]$depth_end,
        interval_length = interval_lengths[[i]],
        value = interval_values[[i]],
        color_bin = bins_by_edge[[i]],
        stringsAsFactors = FALSE
      )
      if (isTRUE(uncertainty)) {
        base_row <- cbind(
          base_row,
          .rateMap_add_derived_uncertainty(edge_results[[i]]$uncertainty)
        )
      }
      base_row
    })
  )
  rownames(intervals) <- NULL
  intervals$clade_key <- rep(
    target_clade_keys,
    times = vapply(interval_lengths, length, integer(1))
  )

  weight_table <- data.frame(
    input_index = keep,
    weight = fit_weights_used,
    ic = resolved_weights$ic,
    IC_used = resolved_weights$IC_used,
    stringsAsFactors = FALSE
  )

  out <- list(
    tree = tree,
    cols = cols,
    lims = lims,
    breaks = breaks,
    values = interval_values,
    intervals = intervals,
    run_values = run_values,
    clade_key = target_clade_keys,
    edge_matches = do.call(cbind, edge_matches),
    summary = summary,
    uncertainty = uncertainty,
    value_summary = value_summary,
    quantile_probs = quantile_probs,
    hpd_prob = hpd_prob,
    plot_value = "value",
    target = target_mode,
    check = check_mode,
    weights = fit_weights_used,
    weight_mode = resolved_weights$mode,
    weight_table = weight_table,
    palette = palette,
    reverse_palette = reverse_palette,
    log = isTRUE(log),
    title = legend_title,
    n_fits = length(trees),
    omitted = omitted
  )
  colnames(out$edge_matches) <- paste0("fit_", keep)
  class(out) <- "rateMap"

  if (isTRUE(plot)) {
    plotRateMap(
      out,
      fsize = fsize,
      tip_fsize = tip_fsize,
      legend_fsize = legend_fsize,
      ftype = ftype,
      show_tip_labels = show_tip_labels,
      lwd = lwd,
      legend = legend,
      outline = outline,
      type = type,
      direction = direction,
      ...
    )
  }

  invisible(out)
}

#' Plot a `rateMap` Object
#'
#' Render a rate-variation map produced by [rateMap()] using
#' [phytools::plotSimmap()] and a continuous color-bar legend.
#' The plotting controls intentionally mirror the `phytools` density-map
#' plotting style for phylogram, fan, and arc layouts.
#'
#' @param x An object of class `"rateMap"` returned by [rateMap()].
#' @param value Character interval column to map to branch colors. The default
#'   `"value"` plots the central estimate chosen by [rateMap()]. When
#'   `uncertainty = TRUE`, useful alternatives include `"mean"`, `"median"`,
#'   `"sd"`, `"ci_width"`, `"hpd_width"`, and `"cv"`.
#' @param palette Optional palette override used for this plot. This accepts the
#'   same values as [rateMap()]'s `palette` argument.
#' @param reverse_palette Optional logical override for palette reversal used
#'   for this plot.
#' @param ... Additional plotting controls. Supported options include
#'   `legend`, `fsize`, `tip_fsize`, `legend_fsize`, `ftype`,
#'   `show_tip_labels`, `outline`, `lwd`, `type`, `mar`, `direction`,
#'   `offset`, `xlim`, `ylim`, `hold`, `underscore`, and `arc_height`.
#'
#' @return Invisibly returns `x`.
#'
#' @seealso [rateMap()], [phytools::densityMap()], [phytools::plotSimmap()]
#'
#' @examples
#' \dontrun{
#' rm_obj <- rateMap(fits, plot = FALSE, progress = FALSE)
#' plotRateMap(
#'   rm_obj,
#'   type = "arc",
#'   show_tip_labels = FALSE,
#'   legend_fsize = 0.8
#' )
#'
#' # If rm_obj was built with uncertainty = TRUE, plot uncertainty directly:
#' plotRateMap(rm_obj, value = "sd", palette = "Inferno")
#' }
#'
#' @export
plotRateMap <- function(x,
                        value = "value",
                        palette = NULL,
                        reverse_palette = NULL,
                        ...) {
  if (!inherits(x, "rateMap")) {
    stop("'x' must be an object of class 'rateMap'.")
  }
  if (!identical(value, "value") || !is.null(palette) || !is.null(reverse_palette)) {
    x <- .rateMap_recolor(
      x,
      value = value,
      palette = palette,
      reverse_palette = reverse_palette
    )
  }

  tree <- x$tree
  cols <- x$cols
  H <- phytools::nodeHeights(tree)

  dots <- list(...)
  legend <- if ("legend" %in% names(dots)) dots$legend else NULL
  fsize <- if ("fsize" %in% names(dots)) dots$fsize else NULL
  tip_fsize <- if ("tip_fsize" %in% names(dots)) dots$tip_fsize else NULL
  legend_fsize <- if ("legend_fsize" %in% names(dots)) dots$legend_fsize else NULL
  ftype <- if ("ftype" %in% names(dots)) dots$ftype else NULL
  show_tip_labels <- if ("show_tip_labels" %in% names(dots)) dots$show_tip_labels else TRUE
  outline <- if ("outline" %in% names(dots)) dots$outline else FALSE
  lwd <- if ("lwd" %in% names(dots)) dots$lwd else 3
  type <- if ("type" %in% names(dots)) dots$type else "phylogram"
  mar <- if ("mar" %in% names(dots)) dots$mar else rep(0.3, 4)
  direction <- if ("direction" %in% names(dots)) dots$direction else "rightwards"
  offset <- if ("offset" %in% names(dots)) dots$offset else NULL
  xlim <- if ("xlim" %in% names(dots)) dots$xlim else NULL
  ylim <- if ("ylim" %in% names(dots)) dots$ylim else NULL
  hold <- if ("hold" %in% names(dots)) dots$hold else TRUE
  underscore <- if ("underscore" %in% names(dots)) dots$underscore else FALSE
  arc_height <- if ("arc_height" %in% names(dots)) dots$arc_height else 2

  type <- match.arg(type, c("phylogram", "fan", "arc"))

  if (length(lwd) == 1L) {
    lwd <- rep(lwd, 2L)
  } else if (length(lwd) > 2L) {
    lwd <- lwd[1:2]
  }

  if (is.null(legend) || isTRUE(legend)) {
    legend <- .rateMap_default_legend_length(type, H)
  }
  show_legend <- .rateMap_show_legend(legend)

  if (is.null(fsize)) {
    fsize <- c(1, 1)
  }
  if (length(fsize) == 1L) {
    fsize <- rep(fsize, 2L)
  }
  if (!is.null(tip_fsize)) {
    fsize[1L] <- tip_fsize
  }
  if (!is.null(legend_fsize)) {
    fsize[2L] <- legend_fsize
  }

  if (is.null(ftype)) {
    ftype <- c("i", "reg")
  }
  if (length(ftype) == 1L) {
    ftype <- c(ftype, "reg")
  }
  if (!isTRUE(show_tip_labels)) {
    ftype[1L] <- "off"
  }

  leg_txt <- c(
    formatC(x$lims[1L], digits = 3, format = "fg"),
    x$title,
    formatC(x$lims[2L], digits = 3, format = "fg")
  )

  if (show_legend && !identical(type, "arc") && legend > max(H)) {
    message("legend scale cannot be longer than total tree length; resetting")
    legend <- 0.5 * max(H)
  }

  if (isTRUE(hold)) {
    grDevices::dev.hold()
  }
  on.exit({
    if (isTRUE(hold)) {
      grDevices::dev.flush()
    }
  }, add = TRUE)

  if (identical(type, "phylogram")) {
    n_tips <- length(tree$tip.label)

    if (show_legend && is.null(ylim) && direction %in% c("rightwards", "leftwards")) {
      ylim <- c(1 - 0.12 * (n_tips - 1L), n_tips)
    }

    if (isTRUE(outline)) {
      old_col <- graphics::par()$col
      graphics::par(col = "transparent")
      outline_offset <- if (is.null(offset)) {
        0.2 * lwd[1L] / 3 + 0.2 / 3
      } else {
        offset + 0.2 * lwd[1L] / 3 + 0.2 / 3
      }
      phytools::plotTree(
        tree,
        fsize = fsize[1L],
        lwd = lwd[1L] + 2,
        offset = outline_offset,
        color = graphics::par()$fg,
        ftype = ftype[1L],
        xlim = xlim,
        ylim = ylim,
        mar = mar,
        direction = direction,
        hold = FALSE,
        add = direction %in% c("upwards", "downwards") && show_legend,
        underscore = underscore
      )
      graphics::par(col = old_col)
    }

    phytools::plotSimmap(
      tree,
      colors = cols,
      pts = FALSE,
      lwd = lwd[1L],
      fsize = fsize[1L],
      mar = mar,
      ftype = ftype[1L],
      add = isTRUE(outline),
      xlim = xlim,
      ylim = ylim,
      direction = direction,
      offset = offset,
      hold = FALSE,
      underscore = underscore
    )

    if (show_legend) {
      phytools::add.color.bar(
        legend,
        cols,
        title = leg_txt[2L],
        lims = as.numeric(leg_txt[c(1L, 3L)]),
        digits = 3,
        prompt = FALSE,
        x = if (identical(direction, "leftwards")) max(H) - legend else 0,
        y = 1 - 0.08 * (n_tips - 1L),
        lwd = lwd[2L],
        fsize = fsize[2L],
        outline = outline,
        direction = if (!is.null(xlim) && xlim[2L] < xlim[1L]) "leftwards" else "rightwards"
      )
    }
  } else {
    .rateMap_with_plotrix_getYmult({
      if (isTRUE(outline)) {
        old_col <- graphics::par()$col
        graphics::par(col = "white")
        invisible(utils::capture.output(
          phytools::plotTree(
            tree,
            type = type,
            lwd = lwd[1L] + 2,
            mar = mar,
            fsize = fsize[1L],
            color = graphics::par()$fg,
            ftype = ftype[1L],
            xlim = xlim,
            ylim = ylim,
            hold = FALSE,
            offset = offset,
            underscore = underscore,
            arc_height = arc_height
          )
        ))
        graphics::par(col = old_col)
      }

      invisible(utils::capture.output(
        phytools::plotSimmap(
          tree,
          colors = cols,
          lwd = lwd[1L],
          mar = mar,
          fsize = fsize[1L],
          add = isTRUE(outline),
          ftype = ftype[1L],
          type = type,
          xlim = xlim,
          ylim = ylim,
          hold = FALSE,
          offset = offset,
          underscore = underscore,
          arc_height = arc_height
        )
      ))

      if (show_legend) {
        if (identical(type, "arc")) {
          phytools::add.color.bar(
            legend,
            cols,
            title = leg_txt[2L],
            lims = as.numeric(leg_txt[c(1L, 3L)]),
            digits = 3,
            outline = outline,
            prompt = FALSE,
            x = mean(graphics::par()$usr[1:2]) - 0.5 * legend,
            y = graphics::par()$usr[3L] + 0.1 * diff(graphics::par()$usr[3:4]),
            lwd = lwd[2L],
            fsize = fsize[2L]
          )
        } else {
          phytools::add.color.bar(
            legend,
            cols,
            title = leg_txt[2L],
            lims = as.numeric(leg_txt[c(1L, 3L)]),
            digits = 3,
            outline = outline,
            prompt = FALSE,
            x = 0.9 * graphics::par()$usr[1L],
            y = 0.9 * graphics::par()$usr[3L],
            lwd = lwd[2L],
            fsize = fsize[2L]
          )
        }
      }
    })
  }

  invisible(x)
}

#' @rdname plotRateMap
#' @method plot rateMap
#' @export
plot.rateMap <- function(x, ...) {
  plotRateMap(x, ...)
}

#' Print a `rateMap` Object
#'
#' @param x An object of class `"rateMap"` returned by [rateMap()].
#' @param ... Ignored.
#'
#' @return Invisibly returns `x`.
#'
#' @export
print.rateMap <- function(x, ...) {
  if (!inherits(x, "rateMap")) {
    stop("'x' must be an object of class 'rateMap'.")
  }

  plot_value <- if (is.null(x$plot_value)) "value" else x$plot_value
  uncertainty <- if (isTRUE(x$uncertainty)) "TRUE" else "FALSE"
  cat("rateMap object\n")
  cat("- fits: ", x$n_fits, "\n", sep = "")
  cat("- summary: ", x$summary, "\n", sep = "")
  cat("- target: ", x$target, "\n", sep = "")
  cat("- check: ", x$check, "\n", sep = "")
  cat("- weights: ", x$weight_mode, "\n", sep = "")
  cat("- uncertainty: ", uncertainty, "\n", sep = "")
  cat("- plotted value: ", plot_value, "\n", sep = "")

  if (!is.null(x$lims) && length(x$lims) == 2L && all(is.finite(x$lims))) {
    cat(
      "- value range: ",
      formatC(x$lims[1L], digits = 4, format = "fg"),
      " to ",
      formatC(x$lims[2L], digits = 4, format = "fg"),
      "\n",
      sep = ""
    )
  }

  if (isTRUE(x$uncertainty) && "sd" %in% names(x$intervals)) {
    sd_range <- range(x$intervals$sd, finite = TRUE)
    if (all(is.finite(sd_range))) {
      cat(
        "- sd range: ",
        formatC(sd_range[1L], digits = 4, format = "fg"),
        " to ",
        formatC(sd_range[2L], digits = 4, format = "fg"),
        "\n",
        sep = ""
      )
    }
  }

  invisible(x)
}
