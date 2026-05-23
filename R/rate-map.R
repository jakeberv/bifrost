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

.rateMap_validate_dots <- function(dots, allowed) {
  if (length(dots) == 0L) {
    return(invisible(NULL))
  }

  dot_names <- names(dots)
  if (is.null(dot_names) || any(!nzchar(dot_names))) {
    stop("Additional arguments in '...' must be named.")
  }

  unknown_dots <- setdiff(dot_names, allowed)
  if (length(unknown_dots) > 0L) {
    stop(
      "Unsupported argument(s) in '...': ",
      paste(unknown_dots, collapse = ", "),
      "."
    )
  }

  invisible(NULL)
}

.rateMap_rate_dot_names <- function() {
  c(
    "value",
    "mar",
    "offset",
    "xlim",
    "ylim",
    "hold",
    "underscore",
    "arc_height",
    "legend_digits"
  )
}

.rateMap_plot_dot_names <- function() {
  c(
    "legend",
    "fsize",
    "tip_fsize",
    "legend_fsize",
    "ftype",
    "show_tip_labels",
    "outline",
    "lwd",
    "type",
    "mar",
    "direction",
    "offset",
    "xlim",
    "ylim",
    "hold",
    "underscore",
    "arc_height",
    "legend_digits"
  )
}

.rateMap_rate_arg_names <- function() {
  setdiff(names(formals(rateMap)), c("fits", "plot", "..."))
}

.rateMap_is_bifrost_search <- function(x) {
  inherits(x, "bifrost_search")
}

.rateMap_is_single_fit <- function(x) {
  .rateMap_is_bifrost_search(x) ||
    inherits(x, "mvgls") ||
    (is.list(x) && inherits(x$model, "mvgls")) ||
    (is.list(x) && !is.null(x$variables$tree) && !is.null(x$param))
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
    ic <- fit$optimal_ic
    if (is.null(ic) || length(ic) < 1L || !is.numeric(ic)) {
      ic <- NA_real_
    } else {
      ic <- as.numeric(ic[1L])
    }

    return(list(
      ic = ic,
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

.rateMap_category_lims <- function(x) {
  categories <- x$rate_categories
  if (is.null(categories) || nrow(categories) == 0L) {
    return(x$lims)
  }

  if (all(c("lower", "upper") %in% names(categories))) {
    bounds <- c(categories$lower, categories$upper)
    bounds <- bounds[is.finite(bounds)]
    if (length(bounds) > 0L) {
      lims <- range(bounds)
      if (diff(lims) > 0) {
        return(lims)
      }
    }
  }

  if ("value" %in% names(categories)) {
    vals <- categories$value[is.finite(categories$value)]
    if (length(vals) > 0L) {
      lims <- range(vals)
      if (diff(lims) > 0) {
        return(lims)
      }
    }
  }

  x$lims
}

.rateMap_add_category_legend <- function(x,
                                         legend,
                                         x_pos,
                                         y_pos,
                                         lwd = 3,
                                         fsize = 1,
                                         outline = FALSE,
                                         digits = .rateMap_legend_digits(x$lims),
                                         direction = "rightwards") {
  categories <- x$rate_categories
  if (is.null(categories) || nrow(categories) == 0L) {
    return(invisible(NULL))
  }

  bar_cols <- x$cols
  if (length(bar_cols) == 1L) {
    bar_cols <- rep(bar_cols, 2L)
  }

  lims <- .rateMap_category_lims(x)
  if (identical(direction, "leftwards")) {
    bar_cols <- rev(bar_cols)
    lims <- rev(lims)
  }

  op <- graphics::par(xpd = NA)
  on.exit(graphics::par(op), add = TRUE)

  usr <- graphics::par("usr")
  yr <- diff(usr[3:4])
  bar_height <- 0.012 * yr * max(lwd / 3, 0.75)
  label_gap <- 0.018 * yr * max(fsize, 0.75)
  title_gap <- 0.04 * yr * max(fsize, 0.75)

  x0 <- x_pos
  x1 <- x_pos + legend
  y0 <- y_pos
  y1 <- y_pos + bar_height
  xs <- seq(x0, x1, length.out = length(bar_cols) + 1L)

  graphics::rect(
    xs[-length(xs)],
    y0,
    xs[-1L],
    y1,
    col = unname(bar_cols),
    border = NA
  )
  if (isTRUE(outline)) {
    graphics::rect(x0, y0, x1, y1, border = "grey30")
  }
  graphics::text(
    c(x0, x1),
    rep(y1 + label_gap, 2L),
    labels = formatC(lims, digits = digits, format = "fg"),
    cex = fsize,
    adj = c(0.5, 0.5)
  )
  graphics::text(
    mean(c(x0, x1)),
    y1 + title_gap,
    labels = x$title,
    cex = fsize
  )
  invisible(NULL)
}

.rateMap_legend_digits <- function(lims) {
  if (!is.numeric(lims) || length(lims) != 2L || any(!is.finite(lims))) {
    return(3L)
  }

  nonzero <- abs(lims[lims != 0])
  if (length(nonzero) == 0L) {
    return(3L)
  }

  scale <- min(nonzero)
  # nocov start
  if (!is.finite(scale)) {
    return(3L)
  }
  # nocov end

  if (scale < 1) {
    return(as.integer(max(3L, min(12L, ceiling(-log10(scale)) + 2L))))
  }

  3L
}

.rateMap_format_rate <- function(x, lims = range(x, finite = TRUE)) {
  digits <- .rateMap_legend_digits(lims)
  trimws(formatC(x, digits = digits, format = "fg"))
}

.rateMap_validate_n_categories <- function(n_categories) {
  if (!is.numeric(n_categories) ||
      length(n_categories) != 1L ||
      !is.finite(n_categories) ||
      n_categories < 1) {
    stop("'n_categories' must be a finite number >= 1.")
  }

  as.integer(n_categories)
}

.rateMap_category_labels <- function(labels, n, user_supplied) {
  if (is.null(labels)) {
    return(NULL)
  }
  if (!is.character(labels) || length(labels) != n || anyNA(labels) || any(!nzchar(labels))) {
    stop("'category_labels' must be a non-missing character vector with one label per rate category.")
  }
  if (anyDuplicated(labels)) {
    if (isTRUE(user_supplied)) {
      stop("'category_labels' must contain unique labels.")
    }
    labels <- make.unique(labels, sep = "_")
  }

  labels
}

.rateMap_resolve_categories <- function(values,
                                        n_categories,
                                        category_breaks = NULL,
                                        category_labels = NULL,
                                        category_bin_method = c("pretty", "equal")) {
  values <- as.numeric(values)
  finite_values <- values[is.finite(values)]
  if (length(finite_values) == 0L) {
    stop("Could not compute finite rate categories from the supplied values.")
  }

  n_categories <- .rateMap_validate_n_categories(n_categories)
  category_bin_method <- match.arg(category_bin_method)
  user_labels <- !is.null(category_labels)

  if (!is.null(category_breaks)) {
    if (!is.numeric(category_breaks) ||
        length(category_breaks) < 2L ||
        any(!is.finite(category_breaks)) ||
        is.unsorted(category_breaks, strictly = TRUE)) {
      stop("'category_breaks' must be a strictly increasing finite numeric vector.")
    }
    if (min(finite_values) < min(category_breaks) || max(finite_values) > max(category_breaks)) {
      stop("'category_breaks' must cover the finite values being mapped.")
    }

    lower <- utils::head(category_breaks, -1L)
    upper <- utils::tail(category_breaks, -1L)
    labels <- .rateMap_category_labels(category_labels, length(lower), user_labels)
    if (is.null(labels)) {
      labels <- paste(
        .rateMap_format_rate(lower, range(category_breaks)),
        .rateMap_format_rate(upper, range(category_breaks)),
        sep = " to "
      )
      labels <- .rateMap_category_labels(labels, length(labels), user_supplied = FALSE)
    }

    bin <- findInterval(values, category_breaks, all.inside = TRUE)
    bin[values == max(category_breaks)] <- length(lower)

    return(list(
      index = bin,
      labels = labels,
      breaks = category_breaks,
      table = data.frame(
        color_bin = seq_along(labels),
        rate_category = labels,
        lower = lower,
        upper = upper,
        value = NA_real_,
        stringsAsFactors = FALSE
      )
    ))
  }

  unique_values <- sort(unique(finite_values))
  if (length(unique_values) <= n_categories) {
    labels <- .rateMap_category_labels(category_labels, length(unique_values), user_labels)
    if (is.null(labels)) {
      labels <- .rateMap_format_rate(unique_values, range(unique_values))
      labels <- .rateMap_category_labels(labels, length(labels), user_supplied = FALSE)
    }

    return(list(
      index = match(values, unique_values),
      labels = labels,
      breaks = unique_values,
      table = data.frame(
        color_bin = seq_along(labels),
        rate_category = labels,
        lower = unique_values,
        upper = unique_values,
        value = unique_values,
        stringsAsFactors = FALSE
      )
    ))
  }

  lims <- range(finite_values)
  breaks <- if (identical(category_bin_method, "equal")) {
    seq(lims[1L], lims[2L], length.out = n_categories + 1L)
  } else {
    pretty(lims, n = n_categories)
  }
  # nocov start
  if (length(breaks) < 2L ||
      min(breaks) > lims[1L] ||
      max(breaks) < lims[2L]) {
    breaks <- seq(lims[1L], lims[2L], length.out = n_categories + 1L)
  }
  # nocov end
  breaks[1L] <- min(breaks[1L], lims[1L])
  breaks[length(breaks)] <- max(breaks[length(breaks)], lims[2L])
  breaks <- unique(breaks)

  lower <- utils::head(breaks, -1L)
  upper <- utils::tail(breaks, -1L)
  labels <- .rateMap_category_labels(category_labels, length(lower), user_labels)
  if (is.null(labels)) {
    labels <- paste(
      .rateMap_format_rate(lower, range(breaks)),
      .rateMap_format_rate(upper, range(breaks)),
      sep = " to "
    )
    labels <- .rateMap_category_labels(labels, length(labels), user_supplied = FALSE)
  }

  bin <- findInterval(values, breaks, all.inside = TRUE)
  bin[values == max(breaks)] <- length(lower)

  list(
    index = bin,
    labels = labels,
    breaks = breaks,
    table = data.frame(
      color_bin = seq_along(labels),
      rate_category = labels,
      lower = lower,
      upper = upper,
      value = NA_real_,
      stringsAsFactors = FALSE
    )
  )
}

.rateMap_build_color_map <- function(values_by_edge,
                                     ncolors,
                                     palette,
                                     reverse_palette,
                                     color_mode,
                                     n_categories,
                                     category_breaks = NULL,
                                     category_labels = NULL,
                                     category_bin_method = "pretty") {
  flat_values <- unlist(values_by_edge, use.names = FALSE)
  lims <- range(flat_values, finite = TRUE)
  # nocov start
  if (!all(is.finite(lims))) {
    stop("Could not compute finite rate limits from the supplied fits.")
  }
  # nocov end

  if (diff(lims) == 0) {
    lims <- lims + c(-0.5, 0.5) * max(abs(lims[1L]), 1) * 1e-6
  }

  color_mode <- match.arg(color_mode, c("continuous", "category"))

  if (identical(color_mode, "continuous")) {
    breaks <- seq(lims[1L], lims[2L], length.out = ncolors + 1L)
    cols <- .rateMap_colors(
      ncolors = ncolors,
      palette = palette,
      reverse_palette = reverse_palette
    )
    names(cols) <- as.character(seq_len(ncolors))
    bins_by_edge <- lapply(values_by_edge, findInterval, vec = breaks, all.inside = TRUE)
    states_by_edge <- lapply(bins_by_edge, as.character)

    return(list(
      color_mode = color_mode,
      cols = cols,
      lims = lims,
      breaks = breaks,
      bins_by_edge = bins_by_edge,
      states_by_edge = states_by_edge,
      rate_categories = NULL,
      category_breaks = NULL,
      category_labels = NULL,
      category_bin_method = NULL
    ))
  }

  categories <- .rateMap_resolve_categories(
    values = flat_values,
    n_categories = n_categories,
    category_breaks = category_breaks,
    category_labels = category_labels,
    category_bin_method = category_bin_method
  )
  category_count <- length(categories$labels)
  cols <- .rateMap_colors(
    ncolors = max(category_count, 2L),
    palette = palette,
    reverse_palette = reverse_palette
  )[seq_len(category_count)]
  names(cols) <- categories$labels
  categories$table$color <- unname(cols[categories$table$rate_category])

  bins_by_edge <- vector("list", length(values_by_edge))
  states_by_edge <- vector("list", length(values_by_edge))
  start <- 1L
  for (i in seq_along(values_by_edge)) {
    n_i <- length(values_by_edge[[i]])
    idx <- categories$index[start:(start + n_i - 1L)]
    bins_by_edge[[i]] <- idx
    states_by_edge[[i]] <- categories$labels[idx]
    start <- start + n_i
  }

  list(
    color_mode = color_mode,
    cols = cols,
    lims = lims,
    breaks = categories$breaks,
    bins_by_edge = bins_by_edge,
    states_by_edge = states_by_edge,
    rate_categories = categories$table,
    category_breaks = categories$breaks,
    category_labels = categories$labels,
    category_bin_method = category_bin_method
  )
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

.rateMap_highest_density_interval <- function(x, weights, prob) {
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
                                          highest_density_interval_prob) {
  q_names <- .rateMap_quantile_names(quantile_probs)
  quantiles <- t(vapply(
    seq_len(nrow(values)),
    function(i) .rateMap_weighted_quantile(values[i, ], weights, quantile_probs),
    numeric(length(quantile_probs))
  ))
  colnames(quantiles) <- q_names

  highest_density_intervals <- t(vapply(
    seq_len(nrow(values)),
    function(i) .rateMap_highest_density_interval(values[i, ], weights, highest_density_interval_prob),
    numeric(2)
  ))
  colnames(highest_density_intervals) <- c("highest_density_interval_low", "highest_density_interval_high")

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
    highest_density_intervals,
    n = rowSums(is.finite(values)),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}

.rateMap_add_derived_uncertainty <- function(x) {
  if (!all(c("q025", "q975", "highest_density_interval_low", "highest_density_interval_high", "mean", "sd") %in% names(x))) {
    return(x)
  }

  x$ci_width <- x$q975 - x$q025
  x$highest_density_interval_width <- x$highest_density_interval_high - x$highest_density_interval_low
  x$cv <- ifelse(is.finite(x$mean) & x$mean != 0, x$sd / abs(x$mean), NA_real_)
  x
}

.rateMap_recolor <- function(x,
                             value,
                             palette = NULL,
                             reverse_palette = NULL,
                             color_mode = NULL,
                             n_categories = NULL,
                             category_breaks = NULL,
                             category_labels = NULL,
                             category_bin_method = NULL) {
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

  selected_palette <- if (is.null(palette)) x$palette else palette
  selected_reverse <- if (is.null(reverse_palette)) {
    isTRUE(x$reverse_palette)
  } else {
    reverse_palette
  }
  selected_color_mode <- if (is.null(color_mode)) {
    if (is.null(x$color_mode)) "continuous" else x$color_mode
  } else {
    match.arg(color_mode, c("continuous", "category"))
  }
  selected_n_categories <- if (is.null(n_categories)) {
    if (is.null(x$n_categories)) 6L else x$n_categories
  } else {
    .rateMap_validate_n_categories(n_categories)
  }
  selected_category_bin_method <- if (is.null(category_bin_method)) {
    if (is.null(x$category_bin_method)) "pretty" else x$category_bin_method
  } else {
    match.arg(category_bin_method, c("pretty", "equal"))
  }
  current_plot_value <- if (is.null(x$plot_value)) "value" else x$plot_value
  reuse_categories <- is.null(color_mode) &&
    identical(selected_color_mode, "category") &&
    identical(value, current_plot_value)
  reuse_exact_categories <- isTRUE(reuse_categories) &&
    !is.null(x$rate_categories) &&
    "value" %in% names(x$rate_categories) &&
    all(is.finite(x$rate_categories$value))
  selected_category_breaks <- if (is.null(category_breaks)) {
    if (isTRUE(reuse_categories) && !isTRUE(reuse_exact_categories)) {
      x$category_breaks
    } else {
      NULL
    }
  } else {
    category_breaks
  }
  selected_category_labels <- if (is.null(category_labels)) {
    if (isTRUE(reuse_categories)) {
      x$category_labels
    } else {
      NULL
    }
  } else {
    category_labels
  }

  ncolors <- if (identical(selected_color_mode, "continuous")) {
    length(x$cols)
  } else {
    max(length(x$cols), 2L)
  }
  color_map <- .rateMap_build_color_map(
    values_by_edge = vals,
    ncolors = ncolors,
    palette = selected_palette,
    reverse_palette = selected_reverse,
    color_mode = selected_color_mode,
    n_categories = selected_n_categories,
    category_breaks = selected_category_breaks,
    category_labels = selected_category_labels,
    category_bin_method = selected_category_bin_method
  )

  for (i in seq_along(x$tree$maps)) {
    names(x$tree$maps[[i]]) <- color_map$states_by_edge[[i]]
  }

  x$tree$mapped.edge <- .rateMap_make_mapped_edge(x$tree$edge, x$tree$maps)
  x$values <- vals
  x$cols <- color_map$cols
  x$lims <- color_map$lims
  x$breaks <- color_map$breaks
  x$intervals$value <- x$intervals[[value]]
  x$intervals$color_bin <- unlist(color_map$bins_by_edge, use.names = FALSE)
  if (identical(selected_color_mode, "category")) {
    x$intervals$rate_category <- unlist(color_map$states_by_edge, use.names = FALSE)
  } else if ("rate_category" %in% names(x$intervals)) {
    x$intervals$rate_category <- NULL
  }
  x$rate_categories <- color_map$rate_categories
  x$palette <- selected_palette
  x$reverse_palette <- selected_reverse
  x$color_mode <- selected_color_mode
  x$n_categories <- selected_n_categories
  x$category_breaks <- color_map$category_breaks
  x$category_labels <- color_map$category_labels
  x$category_bin_method <- color_map$category_bin_method
  x$plot_value <- value
  x$title <- switch(
    value,
    value = x$title,
    mean = if (isTRUE(x$log)) "Mean log fitted rate" else "Mean fitted rate",
    median = if (isTRUE(x$log)) "Median log fitted rate" else "Median fitted rate",
    sd = if (isTRUE(x$log)) "SD log fitted rate" else "SD fitted rate",
    ci_width = if (isTRUE(x$log)) "95% quantile width (log rate)" else "95% quantile width",
    highest_density_interval_width = if (isTRUE(x$log)) {
      "95% highest-density interval width (log rate)"
    } else {
      "95% highest-density interval width"
    },
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
#' follows `bifrost`'s branch-level framing: the first retained tree is used as
#' the plotting scaffold, all retained trees must match in topology and branch
#' lengths, each branch receives one summarized log-rate, and runs are averaged
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
#' [phytools::plotSimmap()] plus either a segmented rate-category color bar or
#' a continuous color-bar legend.
#'
#' Although `rateMap()` is designed for summarizing multiple fitted maps, it also
#' accepts a single completed `bifrost` search or supported fit. Single-fit input
#' is a convenience for applying the same plotting controls to one fitted model
#' before scaling up to multi-run summaries.
#'
#' `rateMap()` can also summarize same-topology posterior or sensitivity samples
#' where branch lengths differ. In that case, supply `check = "topology"` and,
#' usually, an explicit `target_tree`. This can be any target or summary tree
#' with the same topology and tip labels as the inputs; an MCC tree is only one
#' possible choice. The target tree supplies the plotted topology and branch
#' lengths, while rates are matched from each input tree by descendant-tip clade
#' keys rather than by edge order.
#'
#' **Summary modes.** With the default `summary = "branch"` and `log = TRUE`,
#' each edge receives one length-weighted average log-rate from each run before
#' the across-run summary is computed. This is the natural default for `bifrost`
#' searches because shifts are placed at nodes, so a `bifrost` branch is not
#' expected to change regimes internally. With `summary = "interval"`, each
#' target-tree branch is subdivided by the global depth grid controlled by
#' `res`. If source branch lengths differ from the target branch length, source
#' stochastic-map segments are projected onto target intervals by relative
#' position along the matched branch. Interval mode is useful for general
#' stochastic maps that can genuinely change state along a branch.
#'
#' **Log-rate averaging.** When `log = TRUE`, rate parameters are transformed
#' before branch-level, interval-level, and across-fit summaries are computed.
#' With weights `w_i`, the default plotted mean is therefore
#' `sum(w_i * log(rate_i))`, which is the log of a weighted geometric mean. It
#' is not `log(sum(w_i * rate_i))`, the log of a weighted arithmetic mean. Use
#' `log = FALSE` when downstream interpretation requires arithmetic summaries on
#' the original rate scale.
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
#' resolves to the first retained tree. MCC target selection does not make
#' `rateMap()` a mixed-topology summarizer: every target branch must still be
#' present in every retained input. If you already have a preferred consensus,
#' chronogram, MCC, maximum-likelihood, or otherwise curated target tree, pass it
#' with `target_tree` instead of using `target`.
#'
#' **Fit weights.** `weights = "equal"` assigns the same weight to each retained
#' fit. `weights = "ic"` computes standard information-criterion weights from
#' each retained fit's `optimal_ic`, requiring all retained fits to share the
#' same non-missing `IC_used`. Numeric `weights` or `fit_weights` can be supplied
#' for custom weighting. Weights are subset to retained fits after
#' `na_action = "omit"` and are normalized to sum to one.
#'
#' **Uncertainty summaries.** The returned `intervals` data frame is the plotted
#' summary table. With `summary = "branch"`, it has one row per branch. With
#' `summary = "interval"`, it has one row per plotted depth-grid interval. When
#' `uncertainty = TRUE`, this table includes across-fit summaries for each
#' plotted row: weighted mean, weighted median, weighted standard deviation,
#' quantiles, highest-density interval bounds, quantile and highest-density
#' interval widths, coefficient of variation, and the number of finite run-level
#' values. The run-level values are also returned in `run_values` as one matrix
#' per edge.
#' These are weighted empirical summaries of run-level values. For posterior
#' tree samples, use `weights = "equal"` if each retained run represents one
#' posterior draw. `value_summary` controls whether the plotted `value` column
#' uses the weighted mean or weighted median. These summaries are computed on
#' the log-rate scale when `log = TRUE`. Highest-density intervals use the same
#' shortest empirical interval calculation for both equal and unequal fit
#' weights. For `weights = "equal"`, `sd` is the ordinary sample standard
#' deviation of retained run-level values. For unequal weights, `sd` is the
#' square root of the normalized weighted variance, `sum(w * (x - mu)^2)`, with
#' weights normalized to sum to one.
#'
#' **Color modes.** `color_mode = "category"` is the default and assigns plotted
#' rows to ordered rate categories. For a single `bifrost` search with
#' `log = FALSE`, this is the closest formal analogue to the illustrative
#' `generateViridisColorScale()` plot: exact unique fitted rates are preserved
#' as categories when their count is at most `n_categories`, and each displayed
#' category receives one color. For multi-run summaries, the plotted values are
#' usually branch summaries rather than named regimes. If there are more unique
#' finite plotted values than `n_categories`, the values are binned into
#' automatic intervals using `category_bin_method = "pretty"` by default. Use
#' `category_bin_method = "equal"` for equal-width numeric intervals, or
#' `category_breaks` to supply custom category boundaries. These categories are
#' display bins, not additional inferred regimes. `color_mode = "continuous"`
#' uses the density-map-style behavior: summarized numeric values are binned
#' onto a many-color ramp controlled by `ncolors` and shown with a continuous
#' color-bar legend. `pretty()` uses `n_categories` as a target, so the final
#' number of bins can differ and the displayed break range can extend slightly
#' beyond the finite plotted values.
#'
#' @param fits A completed run, fitted model object, or non-empty list of these
#'   objects. Supported shapes are `bifrost_search` objects, `mvgls` objects,
#'   `list(model = <mvgls>)`, and scratch-style lists with `variables$tree` plus
#'   `param`. Single supported objects are wrapped automatically as a convenience
#'   for one-fit plotting and inspection.
#' @param res Integer resolution of the global depth grid used to subdivide
#'   target-tree branches when `summary = "interval"`. Ignored by
#'   `summary = "branch"`.
#' @param fsize Optional plotting font sizes passed to the `"rateMap"` plot
#'   method when `plot = TRUE`.
#' @param ftype Optional plotting font type(s) passed to the `"rateMap"` plot
#'   method when `plot = TRUE`.
#' @param lwd Line width passed through to the `"rateMap"` plot method when
#'   `plot = TRUE`.
#' @param check Logical or character check mode. `TRUE` or `"full"` verifies
#'   that extracted trees and the target tree match in topology and branch
#'   lengths. `"topology"` verifies matching topology/tip labels while allowing
#'   branch lengths to differ. `FALSE` or `"none"` skips this check, but target
#'   branches still must be matchable by descendant-tip sets.
#' @param legend Legend length passed through to the `"rateMap"` plot method when
#'   `plot = TRUE`.
#' @param outline Logical; if `TRUE`, draw branch outlines in the plotted map
#'   when `plot = TRUE`.
#' @param type Plot type to use when `plot = TRUE`. Supported values are
#'   `"phylogram"`, `"fan"`, and `"arc"`.
#' @param direction Plotting direction for `type = "phylogram"` when
#'   `plot = TRUE`.
#' @param plot Logical; if `TRUE`, plot the result immediately using the
#'   `"rateMap"` plot method. If `FALSE`, only return the computed `"rateMap"`
#'   object.
#' @param tree_fun Optional function used to extract a mapped tree from each
#'   element of `fits`. If `NULL`, common `bifrost` and mvgls shapes are
#'   auto-detected.
#' @param param_fun Optional function used to extract a named numeric vector of
#'   state-specific fitted rates from each element of `fits`. If `NULL`, common
#'   `bifrost` and mvgls shapes are auto-detected.
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
#' @param log Logical; if `TRUE` (the default), transform extracted rate
#'   parameters with `log()` before branch, interval, and across-fit averaging
#'   and plotting. Weighted means are then mean log-rates, equivalent to the log
#'   of weighted geometric mean rates. Set `log = FALSE` to summarize rates in
#'   their original units.
#' @param ncolors Number of colors in the rendered palette for
#'   `color_mode = "continuous"`. In `color_mode = "category"`, the rendered
#'   tree uses one color per exact value or category bin; `ncolors` does not
#'   increase the number of displayed category colors.
#' @param palette Palette specification. This can be an `hcl.colors()` palette
#'   name, a vector of colors, or a palette function.
#' @param reverse_palette Logical; if `TRUE`, reverse the resolved palette.
#' @param color_mode Character; `"category"` maps numeric summaries onto
#'   discrete ordered rate categories, while `"continuous"` maps them onto a
#'   continuous color ramp.
#' @param n_categories Maximum number of exact unique values to preserve as
#'   discrete rate categories, and target number of automatic bins when more
#'   unique values are present. Ignored by `color_mode = "continuous"` and by
#'   user-supplied `category_breaks`.
#' @param category_bin_method Character; automatic binning method for
#'   `color_mode = "category"` when there are more unique plotted values than
#'   `n_categories`. `"pretty"` (the default) uses [pretty()] breaks, treating
#'   `n_categories` as a target rather than an exact bin count. `"equal"` uses
#'   equal-width numeric intervals across the finite plotted range. Ignored by
#'   `color_mode = "continuous"` and by user-supplied `category_breaks`.
#' @param category_breaks Optional strictly increasing numeric vector of
#'   category boundaries for `color_mode = "category"`. When supplied, these
#'   breaks must cover all finite plotted values and override automatic
#'   binning.
#' @param category_labels Optional character labels for the displayed categories.
#'   Labels can be used with custom `category_breaks`, exact unique-value
#'   categories, or automatic bins, and must match the final category count.
#' @param legend_title Optional legend title. If `NULL`, defaults to
#'   `"Mean log fitted rate"` or `"Mean fitted rate"` depending on `log`.
#' @param na_action What to do when a run has invalid parameters. `"error"`
#'   stops immediately. `"omit"` drops invalid runs before aggregation.
#' @param summary Character; `"branch"` computes one length-weighted value per
#'   edge, while `"interval"` slices branches on a global depth grid.
#' @param target Character target-tree selection when `target_tree = NULL`.
#'   `"first"` uses the first retained input tree. `"mcc"` chooses the retained
#'   input tree with the highest sum of log clade credibilities. This only
#'   selects the plotting scaffold; every target branch must still be present in
#'   every retained input.
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
#'   uncertainty summaries for every plotted row: whole branches when
#'   `summary = "branch"` and depth-grid intervals when `summary = "interval"`.
#' @param value_summary Character; central estimate stored in the plotted
#'   summary table column `intervals$value` and mapped to colors. `"mean"` uses
#'   the weighted mean. `"median"` uses the weighted median. When `log = TRUE`,
#'   both summaries are computed on the log-rate scale.
#' @param quantile_probs Numeric length-2 vector of quantile probabilities to
#'   report when `uncertainty = TRUE`.
#' @param highest_density_interval_prob Numeric scalar giving the highest-density interval mass to
#'   report when `uncertainty = TRUE`.
#' @param ... Additional plotting arguments passed to the `"rateMap"` plot
#'   method when `plot = TRUE`, such as `value`, `mar`, `offset`, `xlim`,
#'   `ylim`, `hold`, `underscore`, `arc_height`, or `legend_digits`.
#'   Non-empty `...` is rejected when `plot = FALSE` because these controls do
#'   not change the computed `"rateMap"` object. Unsupported arguments are
#'   rejected.
#'
#' @return An object of class `"rateMap"` with components:
#' \describe{
#'   \item{`tree`}{A SIMMAP-style tree whose mapped segments encode color-bin
#'   indices in continuous mode, or named rate categories in category mode.}
#'   \item{`cols`}{The resolved color palette. In continuous mode this has
#'   length `ncolors`; in category mode it has one color per exact value or
#'   category bin.}
#'   \item{`lims`}{Numeric length-2 vector giving the plotted value range.}
#'   \item{`breaks`}{Numeric vector of palette bin boundaries in continuous
#'   mode, or category boundaries/values in category mode.}
#'   \item{`values`}{List of plotted-row central values by edge before color
#'   binning.}
#'   \item{`intervals`}{Plotted summary table. With `summary = "branch"`, this
#'   has one row per branch. With `summary = "interval"`, this has one row per
#'   plotted depth-grid interval.}
#'   \item{`rate_categories`}{Data frame describing discrete rate categories
#'   when `color_mode = "category"`; otherwise `NULL`.}
#'   \item{`run_values`}{When `uncertainty = TRUE`, list of numeric matrices
#'   containing run-level values for each edge. Matrix rows match the plotted
#'   rows for that edge and columns are retained fits. Otherwise `NULL`.}
#'   \item{`clade_key`}{Character descendant-tip key for each target-tree edge.}
#'   \item{`edge_matches`}{Integer matrix mapping target-tree edge rows to
#'   matched source-tree edge rows for each retained fit.}
#'   \item{`summary`}{The summary mode used, `"interval"` or `"branch"`.}
#'   \item{`uncertainty`}{Logical indicating whether uncertainty summaries were
#'   computed.}
#'   \item{`value_summary`}{Central estimate used for the plotted summary table
#'   column `intervals$value`.}
#'   \item{`quantile_probs`}{Quantile probabilities used for uncertainty
#'   summaries.}
#'   \item{`highest_density_interval_prob`}{Highest-density interval mass used for uncertainty summaries.}
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
#'   \item{`color_mode`}{Coloring mode used for the current tree.}
#'   \item{`n_categories`}{Category count target used when
#'   `color_mode = "category"`.}
#'   \item{`category_breaks`}{Category breaks or exact category values used
#'   when `color_mode = "category"`.}
#'   \item{`category_labels`}{Category labels used when
#'   `color_mode = "category"`.}
#'   \item{`category_bin_method`}{Automatic category-binning method used when
#'   `color_mode = "category"` and `category_breaks = NULL`.}
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
#' # A list of completed bifrost searches can be summarized directly:
#' rm_obj <- rateMap(list(search_a, search_b, search_c), plot = FALSE)
#' plot(rm_obj, type = "arc", show_tip_labels = FALSE)
#'
#' # A single completed bifrost search can be inspected with the same API:
#' one_search_map <- rateMap(search_a, plot = FALSE)
#' plot(one_search_map, type = "arc", show_tip_labels = FALSE)
#'
#' # plotRateMap() is the explicit helper and can plot one fitted object
#' # directly by first calling rateMap(..., plot = FALSE):
#' plotRateMap(search_a, type = "arc", show_tip_labels = FALSE)
#'
#' # Scratch-style lists remain supported:
#' scratch_fit <- list(
#'   variables = list(tree = mapped_tree),
#'   param = c("0" = 0.12, "1" = 0.45)
#' )
#' rateMap(list(scratch_fit), plot = FALSE)
#'
#' # Use bifrost model-level IC weights across comparable sensitivity runs:
#' rm_ic <- rateMap(list(search_a, search_b), weights = "ic", plot = FALSE)
#'
#' # Branch-level, discrete-category maps are the bifrost default:
#' branch_rates <- rateMap(list(search_a, search_b), plot = FALSE)
#' branch_rates$rate_categories
#'
#' # To mimic the old jaw-shape preview style for a single search, use raw
#' # fitted rates and category colors:
#' jaw_style <- rateMap(
#'   search_a,
#'   log = FALSE,
#'   color_mode = "category",
#'   palette = viridis::viridis,
#'   plot = FALSE
#' )
#'
#' # Category bins use pretty breaks by default. Use equal-width bins or custom
#' # boundaries when those are easier to compare across figures:
#' equal_bin_rates <- rateMap(
#'   list(search_a, search_b),
#'   n_categories = 5,
#'   category_bin_method = "equal",
#'   plot = FALSE
#' )
#' custom_bin_rates <- rateMap(
#'   list(search_a, search_b),
#'   category_breaks = c(-4, -2, 0, 2),
#'   category_labels = c("slow", "middle", "fast"),
#'   plot = FALSE
#' )
#'
#' # Continuous interval maps remain available for stochastic maps with
#' # along-branch changes:
#' interval_rates <- rateMap(
#'   list(search_a, search_b),
#'   summary = "interval",
#'   color_mode = "continuous",
#'   plot = FALSE
#' )
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
#' plot(posterior_target_rates, value = "sd", type = "arc")
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
#' # Large sensitivity sets can be explicitly subsampled before plotting.
#' # By default, rates are mapped on the log scale:
#' set.seed(1)
#' idx <- sample(seq_along(fit_list), size = 1000)
#' rm_sub <- rateMap(
#'   fit_list[idx],
#'   res = 100,
#'   workers = 8,
#'   future_strategy = "multisession",
#'   palette = c("lightblue", "blue", "pink", "red"),
#'   plot = FALSE
#' )
#'
#' # Use log = FALSE only when a raw-rate scale is preferred:
#' rm_raw <- rateMap(fit_list[idx], log = FALSE, plot = FALSE)
#' plot(rm_sub, type = "arc", show_tip_labels = FALSE, lwd = 1)
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
  log = TRUE,
  ncolors = 256,
  palette = "YlOrRd",
  reverse_palette = FALSE,
  color_mode = c("category", "continuous"),
  n_categories = 6,
  category_bin_method = c("pretty", "equal"),
  category_breaks = NULL,
  category_labels = NULL,
  legend_title = NULL,
  na_action = c("error", "omit"),
  summary = c("branch", "interval"),
  target = c("first", "mcc"),
  target_tree = NULL,
  weights = c("equal", "ic"),
  fit_weights = NULL,
  uncertainty = FALSE,
  value_summary = c("mean", "median"),
  quantile_probs = c(0.025, 0.975),
  highest_density_interval_prob = 0.95,
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
  dots <- list(...)
  .rateMap_validate_dots(dots, .rateMap_rate_dot_names())
  if (!is.logical(plot) || length(plot) != 1L || is.na(plot)) {
    stop("'plot' must be TRUE or FALSE.")
  }
  if (!isTRUE(plot) && length(dots) > 0L) {
    stop("Additional plotting arguments in '...' are only used when 'plot = TRUE'.")
  }

  if (.rateMap_is_single_fit(fits)) {
    fits <- list(fits)
  }
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
  if (!is.numeric(highest_density_interval_prob) ||
      length(highest_density_interval_prob) != 1L ||
      !is.finite(highest_density_interval_prob) ||
      highest_density_interval_prob <= 0 ||
      highest_density_interval_prob > 1) {
    stop("'highest_density_interval_prob' must be a finite number in (0, 1].")
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
  color_mode <- match.arg(color_mode)
  n_categories <- .rateMap_validate_n_categories(n_categories)
  category_bin_method <- match.arg(category_bin_method)
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
        highest_density_interval_prob = highest_density_interval_prob
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

  color_map <- .rateMap_build_color_map(
    values_by_edge = interval_values,
    ncolors = ncolors,
    palette = palette,
    reverse_palette = reverse_palette,
    color_mode = color_mode,
    n_categories = n_categories,
    category_breaks = category_breaks,
    category_labels = category_labels,
    category_bin_method = category_bin_method
  )

  tree$maps <- vector("list", nrow(tree$edge))
  for (i in seq_len(nrow(tree$edge))) {
    tree$maps[[i]] <- interval_lengths[[i]]
    names(tree$maps[[i]]) <- color_map$states_by_edge[[i]]
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
        color_bin = color_map$bins_by_edge[[i]],
        stringsAsFactors = FALSE
      )
      if (identical(color_mode, "category")) {
        base_row$rate_category <- color_map$states_by_edge[[i]]
      }
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
    cols = color_map$cols,
    lims = color_map$lims,
    breaks = color_map$breaks,
    values = interval_values,
    intervals = intervals,
    run_values = run_values,
    rate_categories = color_map$rate_categories,
    clade_key = target_clade_keys,
    edge_matches = do.call(cbind, edge_matches),
    summary = summary,
    uncertainty = uncertainty,
    value_summary = value_summary,
    quantile_probs = quantile_probs,
    highest_density_interval_prob = highest_density_interval_prob,
    plot_value = "value",
    target = target_mode,
    check = check_mode,
    weights = fit_weights_used,
    weight_mode = resolved_weights$mode,
    weight_table = weight_table,
    palette = palette,
    reverse_palette = reverse_palette,
    color_mode = color_mode,
    n_categories = n_categories,
    category_breaks = color_map$category_breaks,
    category_labels = color_map$category_labels,
    category_bin_method = color_map$category_bin_method,
    log = isTRUE(log),
    title = legend_title,
    n_fits = length(trees),
    omitted = omitted
  )
  colnames(out$edge_matches) <- paste0("fit_", keep)
  class(out) <- "rateMap"

  if (isTRUE(plot)) {
    plot_args <- c(list(
      x = out,
      fsize = fsize,
      tip_fsize = tip_fsize,
      legend_fsize = legend_fsize,
      ftype = ftype,
      show_tip_labels = show_tip_labels,
      lwd = lwd,
      legend = legend,
      outline = outline,
      type = type,
      direction = direction
    ), dots)
    do.call(graphics::plot, plot_args)
  }

  invisible(out)
}

#' Plot `rateMap` Objects
#'
#' Render a rate-variation map produced by [rateMap()]. Use `plot(x, ...)` for
#' ordinary S3 dispatch on a `"rateMap"` object. `plotRateMap(x, ...)` is the
#' explicit helper with the same plotting controls, and can also accept one
#' supported fitted object as a convenience.
#'
#' Both interfaces use [phytools::plotSimmap()] and draw either a continuous
#' color-bar legend or a segmented rate-category color bar.
#' The plotting controls intentionally mirror the `phytools` density-map
#' plotting style for phylogram, fan, and arc layouts.
#'
#' @param x An object of class `"rateMap"` returned by [rateMap()], or a single
#'   supported fitted object such as a completed `bifrost_search`. Single fitted
#'   objects are first passed to [rateMap()] with `plot = FALSE`, which is useful
#'   for applying the same plotting controls used for multi-run summaries to one
#'   fitted model.
#' @param value Character column in the plotted summary table (`x$intervals`)
#'   to map to branch colors. The default `"value"` plots the central estimate
#'   chosen by [rateMap()]. When `uncertainty = TRUE`, useful alternatives
#'   include `"mean"`, `"median"`, `"sd"`, `"ci_width"`,
#'   `"highest_density_interval_width"`, and `"cv"`.
#' @param palette Optional palette override used for this plot. This accepts the
#'   same values as [rateMap()]'s `palette` argument.
#' @param reverse_palette Optional logical override for palette reversal used
#'   for this plot.
#' @param color_mode Optional color-mode override. Use `"continuous"` for a
#'   numeric color ramp or `"category"` for ordered discrete rate categories.
#'   If `NULL`, the stored mode in `x` is used. Category mode draws one color
#'   per exact value or category bin; continuous mode draws a many-color ramp.
#' @param n_categories Optional category-count override for
#'   `color_mode = "category"`. This changes the target number of category
#'   bins, not the number of colors in continuous mode.
#' @param category_bin_method Optional category-binning override for
#'   `color_mode = "category"`. Use `"pretty"` for pretty breaks or `"equal"`
#'   for equal-width numeric intervals. Ignored when `category_breaks` is
#'   supplied.
#' @param category_breaks Optional category-break override for
#'   `color_mode = "category"`. Custom breaks override automatic bins.
#' @param category_labels Optional category-label override for
#'   `color_mode = "category"`.
#' @param ... Additional plotting controls. Supported options include
#'   `legend`, `fsize`, `tip_fsize`, `legend_fsize`, `ftype`,
#'   `show_tip_labels`, `outline`, `lwd`, `type`, `mar`, `direction`,
#'   `offset`, `xlim`, `ylim`, `hold`, `underscore`, `arc_height`, and
#'   `legend_digits`. If `legend_digits` is omitted, small-magnitude values use
#'   enough digits to avoid zero-valued legend endpoints. When `x` is a single
#'   fitted object rather than a `"rateMap"`, `...` may also include [rateMap()]
#'   arguments such as `summary`, `uncertainty`, `log`, or `target_tree`.
#'   Unsupported arguments are rejected.
#'
#' @details
#' For `"rateMap"` objects, prefer `plot(x, ...)` for ordinary S3 dispatch. It
#' calls `plot.rateMap()`, which delegates to `plotRateMap()`. The explicit
#' helper is kept so users can find the plotting workflow directly, and for
#' direct plotting of one completed `bifrost_search` or other supported fit.
#' Both interfaces use [phytools::plotSimmap()] for `"phylogram"`, `"fan"`, and
#' `"arc"` layouts. `legend` controls the length of the color bar; set
#' `legend = FALSE` to suppress it. `legend_digits` controls numeric endpoint
#' labels and defaults to enough precision to avoid rounding small values to
#' zero. If `x` is a single supported fitted object rather than a `"rateMap"`,
#' `plotRateMap()` first calls [rateMap()] with `plot = FALSE`, forwarding
#' compatible `rateMap()` arguments through `...`.
#'
#' @return Invisibly returns the plotted `"rateMap"` object.
#'
#' @seealso [rateMap()], [phytools::densityMap()], [phytools::plotSimmap()]
#'
#' @examples
#' \dontrun{
#' rm_obj <- rateMap(fits, plot = FALSE, progress = FALSE)
#' plot(
#'   rm_obj,
#'   type = "arc",
#'   show_tip_labels = FALSE,
#'   legend_fsize = 0.8
#' )
#'
#' # If rm_obj was built with uncertainty = TRUE, plot uncertainty directly:
#' plot(rm_obj, value = "sd", palette = "Inferno")
#'
#' # Direct plotting of one completed search is available as a convenience:
#' plotRateMap(search_a, type = "arc", show_tip_labels = FALSE)
#'
#' # Use a continuous ramp instead of the default ordered rate categories:
#' plot(rm_obj, color_mode = "continuous")
#'
#' # Or keep category colors but change the binning:
#' plot(rm_obj, n_categories = 5, category_bin_method = "equal")
#' plot(rm_obj, category_breaks = c(-4, -2, 0, 2))
#' }
#'
#' @export
plotRateMap <- function(x,
                        value = "value",
                        palette = NULL,
                        reverse_palette = NULL,
                        color_mode = NULL,
                        n_categories = NULL,
                        category_bin_method = NULL,
                        category_breaks = NULL,
                        category_labels = NULL,
                        ...) {
  dots <- list(...)
  if (!inherits(x, "rateMap")) {
    if (.rateMap_is_single_fit(x)) {
      allowed_dots <- union(.rateMap_plot_dot_names(), .rateMap_rate_arg_names())
      .rateMap_validate_dots(dots, allowed_dots)

      rate_arg_names <- setdiff(.rateMap_rate_arg_names(), .rateMap_plot_dot_names())
      plot_arg_names <- .rateMap_plot_dot_names()
      rate_args <- dots[names(dots) %in% rate_arg_names]
      plot_args <- dots[names(dots) %in% plot_arg_names]
      rate_args$fits <- x
      rate_args$plot <- FALSE
      if (!"progress" %in% names(rate_args)) {
        rate_args$progress <- FALSE
      }
      if (!is.null(palette)) {
        rate_args$palette <- palette
      }
      if (!is.null(reverse_palette)) {
        rate_args$reverse_palette <- reverse_palette
      }
      if (!is.null(color_mode)) {
        rate_args$color_mode <- color_mode
      }
      if (!is.null(n_categories)) {
        rate_args$n_categories <- n_categories
      }
      if (!is.null(category_bin_method)) {
        rate_args$category_bin_method <- category_bin_method
      }
      if (!is.null(category_breaks)) {
        rate_args$category_breaks <- category_breaks
      }
      if (!is.null(category_labels)) {
        rate_args$category_labels <- category_labels
      }

      rate_map <- do.call(rateMap, rate_args)
      return(do.call(plotRateMap, c(list(
        x = rate_map,
        value = value,
        palette = palette,
        reverse_palette = reverse_palette,
        color_mode = color_mode,
        n_categories = n_categories,
        category_bin_method = category_bin_method,
        category_breaks = category_breaks,
        category_labels = category_labels
      ), plot_args)))
    }
    stop("'x' must be an object of class 'rateMap' or a single supported fitted object.")
  }

  .rateMap_validate_dots(dots, .rateMap_plot_dot_names())

  if (!identical(value, "value") ||
      !is.null(palette) ||
      !is.null(reverse_palette) ||
      !is.null(color_mode) ||
      !is.null(n_categories) ||
      !is.null(category_bin_method) ||
      !is.null(category_breaks) ||
      !is.null(category_labels)) {
    x <- .rateMap_recolor(
      x,
      value = value,
      palette = palette,
      reverse_palette = reverse_palette,
      color_mode = color_mode,
      n_categories = n_categories,
      category_bin_method = category_bin_method,
      category_breaks = category_breaks,
      category_labels = category_labels
    )
  }

  tree <- x$tree
  cols <- x$cols
  H <- phytools::nodeHeights(tree)

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
  legend_digits <- if ("legend_digits" %in% names(dots)) {
    dots$legend_digits
  } else {
    .rateMap_legend_digits(x$lims)
  }

  if (!is.numeric(legend_digits) ||
      length(legend_digits) != 1L ||
      !is.finite(legend_digits) ||
      legend_digits < 0) {
    stop("'legend_digits' must be a non-negative finite number.")
  }
  legend_digits <- as.integer(legend_digits)

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
  categorical_legend <- show_legend && identical(x$color_mode, "category")

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

    if (show_legend && !categorical_legend) {
      phytools::add.color.bar(
        legend,
        cols,
        title = leg_txt[2L],
        lims = as.numeric(leg_txt[c(1L, 3L)]),
        digits = legend_digits,
        prompt = FALSE,
        x = if (identical(direction, "leftwards")) max(H) - legend else 0,
        y = 1 - 0.08 * (n_tips - 1L),
        lwd = lwd[2L],
        fsize = fsize[2L],
        outline = outline,
        direction = if (!is.null(xlim) && xlim[2L] < xlim[1L]) "leftwards" else "rightwards"
      )
    } else if (categorical_legend) {
      .rateMap_add_category_legend(
        x,
        legend = legend,
        x_pos = if (identical(direction, "leftwards")) max(H) - legend else 0,
        y_pos = 1 - 0.08 * (n_tips - 1L),
        lwd = lwd[2L],
        fsize = fsize[2L],
        outline = outline,
        digits = legend_digits,
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
        if (categorical_legend && identical(type, "arc")) {
          .rateMap_add_category_legend(
            x,
            legend = legend,
            x_pos = mean(graphics::par()$usr[1:2]) - 0.5 * legend,
            y_pos = graphics::par()$usr[3L] + 0.1 * diff(graphics::par()$usr[3:4]),
            lwd = lwd[2L],
            fsize = fsize[2L],
            outline = outline,
            digits = legend_digits
          )
        } else if (categorical_legend) {
          .rateMap_add_category_legend(
            x,
            legend = legend,
            x_pos = 0.9 * graphics::par()$usr[1L],
            y_pos = 0.9 * graphics::par()$usr[3L],
            lwd = lwd[2L],
            fsize = fsize[2L],
            outline = outline,
            digits = legend_digits
          )
        } else if (identical(type, "arc")) {
          phytools::add.color.bar(
            legend,
            cols,
            title = leg_txt[2L],
            lims = as.numeric(leg_txt[c(1L, 3L)]),
            digits = legend_digits,
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
            digits = legend_digits,
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
#' Print a concise summary of a `"rateMap"` object, including the number of fits,
#' summary mode, target/check mode, weighting mode, uncertainty status, color
#' mode, plotted value, value range, and uncertainty ranges when available.
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
  cat("- color mode: ", if (is.null(x$color_mode)) "continuous" else x$color_mode, "\n", sep = "")
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
