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
    return(fit$tree_no_uncertainty_untransformed)
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

  if (!all(param > 0)) {
    return(paste0(
      "Fit ", index,
      " contains non-positive rate parameters; fitted rates must be strictly positive."
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

.rateMap_validate_n_categories <- function(n_categories) {
  if (!is.numeric(n_categories) ||
      length(n_categories) != 1L ||
      !is.finite(n_categories) ||
      n_categories < 1) {
    stop("'n_categories' must be a finite number >= 1.")
  }

  as.integer(n_categories)
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

.rateMap_category_sd <- function(values) {
  if (length(values) <= 1L) {
    return(NA_real_)
  }

  span <- diff(range(values))
  tolerance <- sqrt(.Machine$double.eps) * max(1, max(abs(values)))
  if (span <= tolerance) {
    return(0)
  }

  stats::sd(values)
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
                                     category_bin_method = "pretty",
                                     display_flags_by_edge = NULL,
                                     special_categories = NULL) {
  flat_values <- unlist(values_by_edge, use.names = FALSE)
  lims <- range(flat_values, finite = TRUE)
  # nocov start
  if (!all(is.finite(lims))) {
    stop("Could not compute finite rate limits from the supplied fits.")
  }
  # nocov end
  if (any(!is.finite(flat_values))) {
    stop(
      "Values supplied for rate-map coloring must all be finite; ",
      "cannot assign colors to NA, NaN, or infinite values."
    )
  }

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

  flat_flags <- if (is.null(display_flags_by_edge)) {
    rep("regular", length(flat_values))
  } else {
    unlist(display_flags_by_edge, use.names = FALSE)
  }
  flat_flags[is.na(flat_flags) | !nzchar(flat_flags)] <- "regular"
  use_special <- !is.null(special_categories) && any(flat_flags != "regular")
  regular <- flat_flags == "regular"
  regular_values <- flat_values[regular]

  if (length(regular_values) > 0L) {
    categories <- .rateMap_resolve_categories(
      values = regular_values,
      n_categories = n_categories,
      category_breaks = category_breaks,
      category_labels = category_labels,
      category_bin_method = category_bin_method
    )
    regular_table <- categories$table
    regular_breaks <- categories$breaks
    regular_labels <- categories$labels
    regular_index <- categories$index
  } else {
    regular_table <- data.frame(
      color_bin = integer(),
      rate_category = character(),
      lower = numeric(),
      upper = numeric(),
      value = numeric(),
      stringsAsFactors = FALSE
    )
    regular_breaks <- numeric()
    regular_labels <- character()
    regular_index <- integer()
  }

  zero_present <- isTRUE(use_special) && any(flat_flags == special_categories$zero_label)
  high_present <- isTRUE(use_special) && any(flat_flags == special_categories$high_label)
  special_rows <- list()
  if (zero_present) {
    special_rows <- c(special_rows, list(data.frame(
      color_bin = NA_integer_,
      rate_category = special_categories$zero_label,
      lower = NA_real_,
      upper = NA_real_,
      value = NA_real_,
      stringsAsFactors = FALSE
    )))
  }
  if (high_present) {
    special_rows <- c(special_rows, list(data.frame(
      color_bin = NA_integer_,
      rate_category = special_categories$high_label,
      lower = NA_real_,
      upper = NA_real_,
      value = NA_real_,
      stringsAsFactors = FALSE
    )))
  }

  if (nrow(regular_table) > 0L) {
    regular_cols <- .rateMap_colors(
      ncolors = max(nrow(regular_table), 2L),
      palette = palette,
      reverse_palette = reverse_palette
    )[seq_len(nrow(regular_table))]
    names(regular_cols) <- regular_table$rate_category
    regular_table$color <- unname(regular_cols[regular_table$rate_category])
  } else {
    regular_cols <- character()
    regular_table$color <- character()
  }

  category_table <- regular_table
  if (zero_present) {
    zero_row <- special_rows[[1L]]
    zero_row$color <- special_categories$zero_color
    category_table <- rbind(zero_row, category_table)
  }
  if (high_present) {
    high_row <- special_rows[[length(special_rows)]]
    high_row$color <- special_categories$high_color
    category_table <- rbind(category_table, high_row)
  }
  category_table$color_bin <- seq_len(nrow(category_table))
  cols <- category_table$color
  names(cols) <- category_table$rate_category
  category_lookup <- stats::setNames(category_table$color_bin, category_table$rate_category)

  regular_bins_flat <- rep(NA_integer_, length(flat_values))
  regular_bins_flat[regular] <- regular_index
  if (zero_present || high_present) {
    regular_bins_flat[regular] <- category_lookup[regular_labels[regular_index]]
  }

  bins_by_edge <- vector("list", length(values_by_edge))
  states_by_edge <- vector("list", length(values_by_edge))
  start <- 1L
  for (i in seq_along(values_by_edge)) {
    n_i <- length(values_by_edge[[i]])
    row_idx <- start:(start + n_i - 1L)
    idx <- regular_bins_flat[row_idx]
    states <- flat_flags[row_idx]
    regular_i <- states == "regular"
    if (any(regular_i)) {
      states[regular_i] <- category_table$rate_category[idx[regular_i]]
    }
    if (any(!regular_i)) {
      idx[!regular_i] <- category_lookup[states[!regular_i]]
    }
    bins_by_edge[[i]] <- idx
    states_by_edge[[i]] <- states
    start <- start + n_i
  }

  list(
    color_mode = color_mode,
    cols = cols,
    lims = lims,
    breaks = regular_breaks,
    bins_by_edge = bins_by_edge,
    states_by_edge = states_by_edge,
    rate_categories = category_table,
    category_breaks = regular_breaks,
    category_labels = regular_labels,
    category_bin_method = category_bin_method
  )
}

.rateMap_summarize_category_values <- function(values, lengths) {
  finite_values <- values[is.finite(values)]
  finite_lengths <- lengths[is.finite(lengths)]

  data.frame(
    n_branches = length(values),
    value_mean = if (length(finite_values) > 0L) mean(finite_values) else NA_real_,
    value_median = if (length(finite_values) > 0L) stats::median(finite_values) else NA_real_,
    value_min = if (length(finite_values) > 0L) min(finite_values) else NA_real_,
    value_max = if (length(finite_values) > 0L) max(finite_values) else NA_real_,
    value_sd = .rateMap_category_sd(finite_values),
    total_branch_length = if (length(finite_lengths) > 0L) sum(finite_lengths) else 0,
    stringsAsFactors = FALSE
  )
}

.rateMap_add_category_summaries <- function(rate_categories,
                                            intervals,
                                            summary) {
  if (is.null(rate_categories) || !identical(summary, "branch")) {
    return(rate_categories)
  }

  summaries <- do.call(
    rbind,
    lapply(rate_categories$color_bin, function(bin) {
      rows <- intervals$color_bin == bin
      rows[is.na(rows)] <- FALSE
      .rateMap_summarize_category_values(
        values = intervals$value[rows],
        lengths = intervals$interval_length[rows]
      )
    })
  )
  rownames(summaries) <- NULL

  cbind(rate_categories, summaries)
}

.rateMap_fold_range <- function(rates) {
  rates <- rates[is.finite(rates)]
  if (length(rates) == 0L) {
    return(NA_real_)
  }
  rate_range <- range(rates)
  rate_range[2L] / rate_range[1L]
}

.rateMap_tail_cluster_flags <- function(rates,
                                        tail = c("lower", "upper"),
                                        min_log_gap,
                                        min_fold_reduction,
                                        max_tail_fraction,
                                        min_flagged,
                                        min_regular) {
  tail <- match.arg(tail)
  out <- rep(FALSE, length(rates))
  ok <- is.finite(rates) & rates > 0
  n_ok <- sum(ok)
  empty <- list(
    flag = out,
    cutoff = NA_real_,
    log_cutoff = NA_real_,
    gap = NA_real_,
    fold_reduction = NA_real_,
    score = NA_real_,
    tail_fraction = NA_real_,
    n_tail = 0L,
    n_regular = n_ok
  )
  if (n_ok < min_flagged + min_regular) {
    return(empty)
  }

  log_rates <- log(rates[ok])
  ord <- order(log_rates)
  sorted <- log_rates[ord]
  n <- length(sorted)
  candidates <- seq_len(n - 1L)
  gaps <- diff(sorted)
  total <- sum(sorted)
  cumulative <- cumsum(sorted)
  overall_mean <- mean(sorted)
  lower_n <- candidates
  upper_n <- n - candidates
  lower_mean <- cumulative[candidates] / lower_n
  upper_mean <- (total - cumulative[candidates]) / upper_n
  score <- lower_n * (lower_mean - overall_mean)^2 +
    upper_n * (upper_mean - overall_mean)^2
  full_log_range <- sorted[n] - sorted[1L]

  if (identical(tail, "lower")) {
    n_tail <- candidates
    n_regular <- n - candidates
    tail_fraction <- n_tail / n
    regular_log_range <- sorted[n] - sorted[candidates + 1L]
  } else {
    n_tail <- n - candidates
    n_regular <- candidates
    tail_fraction <- n_tail / n
    regular_log_range <- sorted[candidates] - sorted[1L]
  }
  log_reduction <- full_log_range - regular_log_range

  valid <- n_tail >= min_flagged &
    n_regular >= min_regular &
    tail_fraction <= max_tail_fraction &
    gaps >= min_log_gap &
    log_reduction >= log(min_fold_reduction) &
    is.finite(score)
  if (!any(valid)) {
    return(empty)
  }

  valid_candidates <- candidates[valid]
  best <- valid_candidates[which.max(score[valid])]
  ok_indices <- which(ok)[ord]
  if (identical(tail, "lower")) {
    selected <- ok_indices[seq_len(best)]
  } else {
    selected <- ok_indices[(best + 1L):n]
  }
  out[selected] <- TRUE

  log_cutoff <- mean(sorted[c(best, best + 1L)])
  list(
    flag = out,
    cutoff = exp(log_cutoff),
    log_cutoff = log_cutoff,
    gap = gaps[best],
    fold_reduction = exp(log_reduction[best]),
    score = score[best],
    tail_fraction = tail_fraction[best],
    n_tail = n_tail[best],
    n_regular = n_regular[best]
  )
}

.rateMap_compute_rate_flags <- function(values_by_edge,
                                        summary,
                                        log,
                                        rate_flags) {
  flat_values <- unlist(values_by_edge, use.names = FALSE)
  regular <- rep("regular", length(flat_values))
  flags_by_edge <- rep(list(character()), length(values_by_edge))

  if (isTRUE(rate_flags$disabled)) {
    return(list(
      rate_for_flagging = rep(NA_real_, length(flat_values)),
      rate_flag = regular,
      flags_by_edge = flags_by_edge,
      diagnostics = list(
        enabled = FALSE,
        reason = "rate flagging is disabled",
        method = "disabled"
      )
    ))
  }

  diagnostics <- list(
    enabled = FALSE,
    reason = NA_character_,
    method = rate_flags$method,
    zero_floor = if (is.null(rate_flags$zero_floor)) NA_real_ else rate_flags$zero_floor,
    cluster_min_log_gap = rate_flags$cluster_min_log_gap,
    cluster_min_fold_reduction = rate_flags$cluster_min_fold_reduction,
    cluster_max_tail_fraction = rate_flags$cluster_max_tail_fraction,
    cluster_min_flagged = rate_flags$cluster_min_flagged,
    cluster_min_regular = rate_flags$cluster_min_regular,
    zero_label = rate_flags$zero_label,
    high_label = rate_flags$high_label
  )

  if (!identical(summary, "branch")) {
    diagnostics$reason <- "rate flags are only computed for branch summaries"
    return(list(
      rate_for_flagging = rep(NA_real_, length(flat_values)),
      rate_flag = regular,
      flags_by_edge = flags_by_edge,
      diagnostics = diagnostics
    ))
  }

  rates <- if (isTRUE(log)) base::exp(flat_values) else flat_values
  near_zero <- rep(FALSE, length(rates))
  high_outlier <- rep(FALSE, length(rates))
  lower_cluster <- list(
    flag = near_zero, cutoff = NA_real_, log_cutoff = NA_real_,
    gap = NA_real_, fold_reduction = NA_real_, score = NA_real_,
    tail_fraction = NA_real_, n_tail = 0L, n_regular = length(rates)
  )
  upper_cluster <- list(
    flag = high_outlier, cutoff = NA_real_, log_cutoff = NA_real_,
    gap = NA_real_, fold_reduction = NA_real_, score = NA_real_,
    tail_fraction = NA_real_, n_tail = 0L, n_regular = length(rates)
  )

  if (isTRUE(rate_flags$near_zero) &&
      identical(rate_flags$method, "floor") &&
      !is.null(rate_flags$zero_floor)) {
    near_zero <- near_zero | (is.finite(rates) & rates <= rate_flags$zero_floor)
  }
  if (isTRUE(rate_flags$near_zero) && identical(rate_flags$method, "tail_cluster")) {
    lower_cluster <- .rateMap_tail_cluster_flags(
      rates,
      tail = "lower",
      min_log_gap = rate_flags$cluster_min_log_gap,
      min_fold_reduction = rate_flags$cluster_min_fold_reduction,
      max_tail_fraction = rate_flags$cluster_max_tail_fraction,
      min_flagged = rate_flags$cluster_min_flagged,
      min_regular = rate_flags$cluster_min_regular
    )
    near_zero <- near_zero | lower_cluster$flag
  }
  if (isTRUE(rate_flags$high_outlier) && identical(rate_flags$method, "tail_cluster")) {
    upper_cluster <- .rateMap_tail_cluster_flags(
      rates,
      tail = "upper",
      min_log_gap = rate_flags$cluster_min_log_gap,
      min_fold_reduction = rate_flags$cluster_min_fold_reduction,
      max_tail_fraction = rate_flags$cluster_max_tail_fraction,
      min_flagged = rate_flags$cluster_min_flagged,
      min_regular = rate_flags$cluster_min_regular
    )
    high_outlier <- high_outlier | upper_cluster$flag
  }

  high_outlier <- high_outlier & !near_zero
  flag <- regular
  flag[near_zero] <- rate_flags$zero_label
  flag[high_outlier] <- rate_flags$high_label
  groups <- rep(seq_along(values_by_edge), lengths(values_by_edge))
  regular_rates <- rates[flag == "regular"]

  diagnostics$enabled <- TRUE
  diagnostics$reason <- NA_character_
  diagnostics$n_total <- length(rates)
  diagnostics$n_regular <- sum(flag == "regular")
  diagnostics$n_near_zero <- sum(near_zero)
  diagnostics$n_high_outlier <- sum(high_outlier)
  diagnostics$min_rate <- if (any(is.finite(rates))) min(rates, na.rm = TRUE) else NA_real_
  diagnostics$max_rate <- if (any(is.finite(rates))) max(rates, na.rm = TRUE) else NA_real_
  diagnostics$min_regular_rate <- if (any(is.finite(regular_rates))) min(regular_rates, na.rm = TRUE) else NA_real_
  diagnostics$max_regular_rate <- if (any(is.finite(regular_rates))) max(regular_rates, na.rm = TRUE) else NA_real_
  diagnostics$full_fold_range <- .rateMap_fold_range(rates)
  diagnostics$regular_fold_range <- .rateMap_fold_range(regular_rates)
  diagnostics$zero_cluster_cutoff_rate <- lower_cluster$cutoff
  diagnostics$zero_cluster_log_cutoff <- lower_cluster$log_cutoff
  diagnostics$zero_cluster_log_gap <- lower_cluster$gap
  diagnostics$zero_cluster_fold_reduction <- lower_cluster$fold_reduction
  diagnostics$zero_cluster_score <- lower_cluster$score
  diagnostics$zero_cluster_tail_fraction <- lower_cluster$tail_fraction
  diagnostics$zero_cluster_n_tail <- lower_cluster$n_tail
  diagnostics$zero_cluster_n_regular <- lower_cluster$n_regular
  diagnostics$high_cluster_cutoff_rate <- upper_cluster$cutoff
  diagnostics$high_cluster_log_cutoff <- upper_cluster$log_cutoff
  diagnostics$high_cluster_log_gap <- upper_cluster$gap
  diagnostics$high_cluster_fold_reduction <- upper_cluster$fold_reduction
  diagnostics$high_cluster_score <- upper_cluster$score
  diagnostics$high_cluster_tail_fraction <- upper_cluster$tail_fraction
  diagnostics$high_cluster_n_tail <- upper_cluster$n_tail
  diagnostics$high_cluster_n_regular <- upper_cluster$n_regular
  diagnostics$value_scale <- if (isTRUE(log)) "rate (exp plotted values)" else "rate"

  list(
    rate_for_flagging = rates,
    rate_flag = flag,
    flags_by_edge = split(flag, groups),
    diagnostics = diagnostics
  )
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
                                     weights) {
  n_fits <- length(fits)
  n_used <- length(keep)

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

.rateMap_require_packages <- function() {
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
  invisible(NULL)
}

#' Normalize Summary Options for `rateMap()`
#'
#' @keywords internal
#' @noRd
.rateMap_normalize_summary_options <- function(summary = c("branch", "interval"),
                                               res = 100,
                                               log = TRUE,
                                               na_action = c("error", "omit")) {
  if (!is.numeric(res) || length(res) != 1L || !is.finite(res) || res < 1) {
    stop("'res' must be a finite number >= 1.")
  }

  list(
    mode = match.arg(summary),
    res = as.integer(res),
    log = isTRUE(log),
    na_action = match.arg(na_action)
  )
}

#' Normalize Target and Tree-Check Options for `rateMap()`
#'
#' @keywords internal
#' @noRd
.rateMap_normalize_target_options <- function(target = c("first", "mcc"),
                                              target_tree = NULL,
                                              check = TRUE) {
  list(
    mode = if (is.null(target_tree)) match.arg(target) else "user",
    selection = match.arg(target),
    tree = target_tree,
    check = .rateMap_normalize_check(check)
  )
}

#' Normalize Uncertainty Options for `rateMap()`
#'
#' @keywords internal
#' @noRd
.rateMap_normalize_uncertainty_options <- function(
  uncertainty = FALSE,
  value_summary = c("mean", "median"),
  quantile_probs = c(0.025, 0.975),
  highest_density_interval_prob = 0.95
) {
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

  value_summary <- match.arg(value_summary)
  list(
    enabled = isTRUE(uncertainty),
    value_summary = value_summary,
    quantile_probs = quantile_probs,
    highest_density_interval_prob = highest_density_interval_prob,
    summarize_runs = isTRUE(uncertainty) || identical(value_summary, "median")
  )
}

#' Normalize Color Options for `rateMap()`
#'
#' @keywords internal
#' @noRd
.rateMap_normalize_color_options <- function(ncolors = 256,
                                             palette = "YlOrRd",
                                             reverse_palette = FALSE,
                                             color_mode = c("category", "continuous"),
                                             n_categories = 6,
                                             category_bin_method = c("pretty", "equal"),
                                             category_breaks = NULL,
                                             category_labels = NULL,
                                             legend_title = NULL,
                                             log = TRUE) {
  if (!is.numeric(ncolors) || length(ncolors) != 1L || !is.finite(ncolors) || ncolors < 2) {
    stop("'ncolors' must be a finite number >= 2.")
  }

  color_mode <- match.arg(color_mode)
  n_categories <- .rateMap_validate_n_categories(n_categories)
  category_bin_method <- match.arg(category_bin_method)
  if (is.null(legend_title)) {
    legend_title <- if (isTRUE(log)) "Mean log fitted rate" else "Mean fitted rate"
  }

  list(
    mode = color_mode,
    ncolors = as.integer(ncolors),
    palette = palette,
    reverse_palette = reverse_palette,
    n_categories = n_categories,
    category_bin_method = category_bin_method,
    category_breaks = category_breaks,
    category_labels = category_labels,
    title = legend_title
  )
}

#' Normalize Compute Options for `rateMap()`
#'
#' @keywords internal
#' @noRd
.rateMap_normalize_compute_options <- function(workers = NULL,
                                               future_strategy = c("multisession", "multicore"),
                                               future_seed = FALSE,
                                               future_scheduling = 1,
                                               future_chunk_size = NULL,
                                               progress = TRUE) {
  if (!is.null(workers) &&
      (!is.numeric(workers) || length(workers) != 1L || !is.finite(workers) || workers < 1)) {
    stop("'workers' must be NULL or a positive number.")
  }

  list(
    workers = if (is.null(workers)) NULL else as.integer(workers),
    future_strategy = match.arg(future_strategy),
    future_seed = future_seed,
    future_scheduling = future_scheduling,
    future_chunk_size = future_chunk_size,
    progress = isTRUE(progress)
  )
}

.rateMap_default_control_values <- function() {
  list(
    res = 100,
    check = TRUE,
    target = "first",
    tree_fun = NULL,
    param_fun = NULL,
    na_action = "error",
    rate_flags = rateMapRateFlags(),
    quantile_probs = c(0.025, 0.975),
    highest_density_interval_prob = 0.95,
    future_strategy = "multisession",
    future_seed = FALSE,
    future_scheduling = 1,
    future_chunk_size = NULL
  )
}

#' Control Near-Zero and Tail Rate Diagnostics for `rateMap()`
#'
#' Build a rate-flagging control object for [rateMapControl()]. These options
#' identify branch-rate summaries that are effectively zero, or optionally in a
#' separated high-rate tail, without deleting or changing the fitted rates. Flags
#' are reported in the returned `intervals` table and summarized in
#' `rate_diagnostics`.
#'
#' @param near_zero Logical or `NULL`. If `TRUE`, flag lower-tail or floor-level
#'   rates. The default `NULL` resolves to `FALSE` for `method = "none"` and
#'   `TRUE` for `method = "floor"` or `method = "tail_cluster"`.
#' @param high_outlier Logical; if `TRUE`, flag isolated upper-tail rates.
#' @param method Optional detection method. `NULL` chooses `"floor"` when
#'   `zero_floor` is supplied and `"none"` otherwise. `"floor"` uses
#'   `zero_floor` for near-zero rates. `"tail_cluster"` uses an Otsu-style
#'   guarded two-class split of log rates to identify a separated lower-tail
#'   cluster, and, when `high_outlier = TRUE`, a separated upper-tail cluster.
#'   `"none"` computes diagnostics without adding special rate flags. Inactive
#'   switches are normalized to `FALSE`: `method = "none"` sets `near_zero` and
#'   `high_outlier` to `FALSE`, and `method = "floor"` sets `high_outlier` to
#'   `FALSE`.
#' @param zero_floor Optional non-negative rate floor. When supplied with
#'   `method = NULL`, the method resolves to `"floor"`. Finite rates less than
#'   or equal to this value are flagged as near-zero.
#' @param cluster_min_log_gap Minimum log-rate gap required between a
#'   cluster-defined tail and the remaining regular rates.
#' @param cluster_min_fold_reduction Minimum fold-range reduction required after
#'   separating a cluster-defined tail from the regular rates.
#' @param cluster_max_tail_fraction Maximum fraction of positive finite rates
#'   that can be assigned to a cluster-defined tail.
#' @param cluster_min_flagged Minimum number of branches required for a
#'   cluster-defined tail.
#' @param cluster_min_regular Minimum number of branches that must remain in the
#'   regular rate set after separating a cluster-defined tail.
#' @param zero_label,high_label Labels used for flagged display categories.
#' @param zero_color,high_color Colors used for flagged display categories in
#'   category mode.
#'
#' @details
#' The `"tail_cluster"` method is an Otsu-style diagnostic adapted to sorted
#' branch log rates. It evaluates two-class splits and keeps a tail only when
#' guardrails for log-rate gap size, fold-range reduction, tail fraction, and
#' minimum counts are all satisfied. It is display metadata, not data deletion,
#' model correction, or formal threshold selection.
#'
#' @return A normalized `"rateMap_rate_flags"` object for
#'   `rateMapControl(rate_flags = )`.
#'
#' @examples
#' \dontrun{
#' ctrl <- rateMapControl(rate_flags = rateMapRateFlags(zero_floor = 1e-8))
#' rm_obj <- rateMap(fits, control = ctrl)
#' rm_obj$rate_diagnostics
#'
#' cluster_ctrl <- rateMapControl(
#'   rate_flags = rateMapRateFlags(method = "tail_cluster")
#' )
#' rm_obj <- rateMap(fits, control = cluster_ctrl)
#' rm_obj$rate_diagnostics
#' }
#'
#' @references
#' Otsu, N. (1979). A threshold selection method from gray-level histograms.
#' *IEEE Transactions on Systems, Man, and Cybernetics*, 9(1), 62-66.
#' \doi{10.1109/TSMC.1979.4310076}
#'
#' @export
rateMapRateFlags <- function(near_zero = NULL,
                             high_outlier = FALSE,
                             method = NULL,
                             zero_floor = NULL,
                             cluster_min_log_gap = log(10),
                             cluster_min_fold_reduction = 10,
                             cluster_max_tail_fraction = 0.4,
                             cluster_min_flagged = 3,
                             cluster_min_regular = 10,
                             zero_label = "near-zero",
                             high_label = "high-outlier",
                             zero_color = "grey70",
                             high_color = "black") {
  if (!is.null(near_zero) &&
      (!is.logical(near_zero) || length(near_zero) != 1L || is.na(near_zero))) {
    stop("'near_zero' must be NULL, TRUE, or FALSE.")
  }
  if (!is.logical(high_outlier) || length(high_outlier) != 1L || is.na(high_outlier)) {
    stop("'high_outlier' must be TRUE or FALSE.")
  }
  if (!is.null(zero_floor) &&
      (!is.numeric(zero_floor) || length(zero_floor) != 1L ||
       !is.finite(zero_floor) || zero_floor < 0)) {
    stop("'zero_floor' must be NULL or a finite non-negative number.")
  }
  if (is.null(method)) {
    method <- if (is.null(zero_floor)) "none" else "floor"
  } else {
    if (!is.character(method) || length(method) != 1L || is.na(method)) {
      stop("'method' must be one of 'none', 'floor', or 'tail_cluster'.")
    }
    method <- match.arg(method, c("none", "floor", "tail_cluster"))
  }
  if (identical(method, "floor") && is.null(zero_floor)) {
    stop("method = 'floor' requires 'zero_floor'.")
  }
  if (!identical(method, "floor") && !is.null(zero_floor)) {
    stop("'zero_floor' can only be used with method = 'floor'.")
  }
  if (is.null(near_zero)) {
    near_zero <- !identical(method, "none")
  }
  if (identical(method, "none")) {
    near_zero <- FALSE
    high_outlier <- FALSE
  } else if (identical(method, "floor")) {
    high_outlier <- FALSE
  }
  if (!is.numeric(cluster_min_log_gap) || length(cluster_min_log_gap) != 1L ||
      !is.finite(cluster_min_log_gap) || cluster_min_log_gap <= 0) {
    stop("'cluster_min_log_gap' must be a finite positive number.")
  }
  if (!is.numeric(cluster_min_fold_reduction) ||
      length(cluster_min_fold_reduction) != 1L ||
      !is.finite(cluster_min_fold_reduction) ||
      cluster_min_fold_reduction <= 1) {
    stop("'cluster_min_fold_reduction' must be a finite number greater than 1.")
  }
  if (!is.numeric(cluster_max_tail_fraction) ||
      length(cluster_max_tail_fraction) != 1L ||
      !is.finite(cluster_max_tail_fraction) ||
      cluster_max_tail_fraction <= 0 ||
      cluster_max_tail_fraction > 0.5) {
    stop("'cluster_max_tail_fraction' must be a finite number in (0, 0.5].")
  }
  if (!is.numeric(cluster_min_flagged) || length(cluster_min_flagged) != 1L ||
      !is.finite(cluster_min_flagged) || cluster_min_flagged < 1) {
    stop("'cluster_min_flagged' must be a finite number >= 1.")
  }
  if (!is.numeric(cluster_min_regular) || length(cluster_min_regular) != 1L ||
      !is.finite(cluster_min_regular) || cluster_min_regular < 1) {
    stop("'cluster_min_regular' must be a finite number >= 1.")
  }
  labels <- c(zero_label, high_label)
  if (!is.character(labels) || anyNA(labels) || any(!nzchar(labels)) || anyDuplicated(labels)) {
    stop("'zero_label' and 'high_label' must be non-empty unique strings.")
  }
  colors <- c(zero_color, high_color)
  if (!is.character(colors) || anyNA(colors) || any(!nzchar(colors))) {
    stop("'zero_color' and 'high_color' must be non-empty strings.")
  }

  structure(
    list(
      near_zero = near_zero,
      high_outlier = high_outlier,
      method = method,
      zero_floor = zero_floor,
      cluster_min_log_gap = cluster_min_log_gap,
      cluster_min_fold_reduction = cluster_min_fold_reduction,
      cluster_max_tail_fraction = cluster_max_tail_fraction,
      cluster_min_flagged = as.integer(cluster_min_flagged),
      cluster_min_regular = as.integer(cluster_min_regular),
      zero_label = zero_label,
      high_label = high_label,
      zero_color = zero_color,
      high_color = high_color
    ),
    class = c("rateMap_rate_flags", "list")
  )
}

.rateMap_disabled_rate_flags <- function() {
  defaults <- unclass(rateMapRateFlags())
  defaults$disabled <- TRUE
  defaults$method <- "disabled"
  class(defaults) <- c("rateMap_rate_flags", "list")
  defaults
}

.rateMap_normalize_rate_flags <- function(rate_flags) {
  if (is.null(rate_flags)) {
    return(.rateMap_disabled_rate_flags())
  }
  if (!is.list(rate_flags)) {
    stop("'rate_flags' must be created by rateMapRateFlags() or be a named list.")
  }
  if (isTRUE(rate_flags$disabled)) {
    return(.rateMap_disabled_rate_flags())
  }
  rate_flag_names <- names(rate_flags)
  if (length(rate_flags) > 0L &&
      (is.null(rate_flag_names) || any(!nzchar(rate_flag_names)))) {
    stop("'rate_flags' must be a named list.")
  }

  defaults <- unclass(rateMapRateFlags())
  unknown <- setdiff(rate_flag_names, names(defaults))
  if (length(unknown) > 0L) {
    stop("Unsupported rateMapRateFlags option(s): ", paste(unknown, collapse = ", "), ".")
  }
  defaults[rate_flag_names] <- unclass(rate_flags)[rate_flag_names]
  if (!("near_zero" %in% rate_flag_names)) {
    defaults$near_zero <- NULL
  }
  if (!("method" %in% rate_flag_names) &&
      ("zero_floor" %in% rate_flag_names) &&
      !is.null(rate_flags$zero_floor)) {
    defaults$method <- NULL
  }
  do.call(rateMapRateFlags, defaults)
}

#' Control Advanced `rateMap()` Options
#'
#' Build a control object for less common `rateMap()` settings. These options are
#' useful for interval maps, same-topology tree samples, custom fit object
#' shapes, invalid-fit handling, and future-based parallel scheduling.
#'
#' @param res Integer resolution of the global depth grid used to subdivide
#'   target-tree branches when `summary = "interval"`. Ignored by
#'   `summary = "branch"`.
#' @param check Logical or character check mode. `TRUE` or `"full"` verifies
#'   that extracted trees and the target tree match in topology and branch
#'   lengths. `"topology"` verifies matching topology/tip labels while allowing
#'   branch lengths to differ. `FALSE` or `"none"` skips the upfront
#'   [ape::all.equal.phylo()] check, but target branches still must be
#'   matchable by descendant-tip sets.
#' @param target Character target-tree selection when `target_tree = NULL` in
#'   [rateMap()]. `"first"` uses the first retained input tree. `"mcc"` chooses
#'   the retained input tree with the highest sum of log clade credibilities.
#' @param tree_fun Optional function used to extract a mapped tree from each
#'   element of `fits`. If `NULL`, common `bifrost` and mvgls shapes are
#'   auto-detected. For `"bifrost_search"` objects, the default uses
#'   `tree_no_uncertainty_untransformed` to preserve the original branch-length
#'   scale; supply `tree_fun` explicitly to map a different tree field.
#' @param param_fun Optional function used to extract a named numeric vector of
#'   state-specific fitted rates from each element of `fits`. If `NULL`, common
#'   `bifrost` and mvgls shapes are auto-detected.
#' @param na_action What to do when a run has invalid parameters. `"error"`
#'   stops immediately. `"omit"` drops invalid runs before aggregation.
#' @param rate_flags A `"rateMap_rate_flags"` object from [rateMapRateFlags()],
#'   or a named list of rate-flagging options. Use `NULL` to disable rate
#'   diagnostics.
#' @param quantile_probs Numeric length-2 vector of quantile probabilities to
#'   report when `uncertainty = TRUE`.
#' @param highest_density_interval_prob Numeric scalar giving the
#'   highest-density interval mass to report when `uncertainty = TRUE`.
#' @param future_strategy Future backend used only when `workers` is supplied to
#'   [rateMap()]. Must be one of `"multisession"` or `"multicore"`.
#' @param future_seed Seed control passed to [future.apply::future_lapply()].
#' @param future_scheduling Scheduling control passed to
#'   [future.apply::future_lapply()].
#' @param future_chunk_size Chunk size passed to
#'   [future.apply::future_lapply()].
#'
#' @return A `"rateMap_control"` object for the `control` argument of
#'   [rateMap()].
#'
#' @examples
#' \dontrun{
#' # Same-topology tree samples with differing branch lengths:
#' ctrl <- rateMapControl(check = "topology")
#' posterior_rates <- rateMap(
#'   posterior_fit_list,
#'   target_tree = summary_tree,
#'   uncertainty = TRUE,
#'   control = ctrl
#' )
#'
#' # Interval maps with a finer depth grid:
#' interval_rates <- rateMap(
#'   fits,
#'   summary = "interval",
#'   control = rateMapControl(res = 200)
#' )
#' }
#'
#' @export
rateMapControl <- function(res = 100,
                           check = TRUE,
                           target = c("first", "mcc"),
                           tree_fun = NULL,
                           param_fun = NULL,
                           na_action = c("error", "omit"),
                           rate_flags = rateMapRateFlags(),
                           quantile_probs = c(0.025, 0.975),
                           highest_density_interval_prob = 0.95,
                           future_strategy = c("multisession", "multicore"),
                           future_seed = FALSE,
                           future_scheduling = 1,
                           future_chunk_size = NULL) {
  summary_options <- .rateMap_normalize_summary_options(
    summary = "branch",
    res = res,
    log = TRUE,
    na_action = na_action
  )
  uncertainty_options <- .rateMap_normalize_uncertainty_options(
    uncertainty = FALSE,
    value_summary = "mean",
    quantile_probs = quantile_probs,
    highest_density_interval_prob = highest_density_interval_prob
  )
  target <- match.arg(target)
  check <- .rateMap_normalize_check(check)
  future_strategy <- match.arg(future_strategy)

  if (!is.null(tree_fun) && !is.function(tree_fun)) {
    stop("'tree_fun' must be NULL or a function.")
  }
  if (!is.null(param_fun) && !is.function(param_fun)) {
    stop("'param_fun' must be NULL or a function.")
  }
  rate_flags <- .rateMap_normalize_rate_flags(rate_flags)

  structure(
    list(
      res = summary_options$res,
      check = check,
      target = target,
      tree_fun = tree_fun,
      param_fun = param_fun,
      na_action = summary_options$na_action,
      rate_flags = rate_flags,
      quantile_probs = uncertainty_options$quantile_probs,
      highest_density_interval_prob = uncertainty_options$highest_density_interval_prob,
      future_strategy = future_strategy,
      future_seed = future_seed,
      future_scheduling = future_scheduling,
      future_chunk_size = future_chunk_size
    ),
    class = c("rateMap_control", "list")
  )
}

.rateMap_normalize_control <- function(control) {
  if (is.null(control) || !is.list(control)) {
    stop("'control' must be created by rateMapControl() or be a named list.")
  }

  control_names <- names(control)
  if (length(control) > 0L &&
      (is.null(control_names) || any(!nzchar(control_names)))) {
    stop("'control' must be a named list.")
  }

  defaults <- .rateMap_default_control_values()
  unknown <- setdiff(control_names, names(defaults))
  if (length(unknown) > 0L) {
    stop(
      "Unsupported rateMapControl option(s): ",
      paste(unknown, collapse = ", "),
      "."
    )
  }

  merged <- defaults
  merged[control_names] <- unclass(control)[control_names]
  do.call(rateMapControl, merged)
}

#' Normalize Flat `rateMap()` Options into Internal Groups
#'
#' @keywords internal
#' @noRd
.rateMap_normalize_grouped_options <- function(
  res = 100,
  check = TRUE,
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
  rate_flags = rateMapRateFlags(),
  summary = c("branch", "interval"),
  target = c("first", "mcc"),
  target_tree = NULL,
  uncertainty = FALSE,
  value_summary = c("mean", "median"),
  quantile_probs = c(0.025, 0.975),
  highest_density_interval_prob = 0.95
) {
  summary_options <- .rateMap_normalize_summary_options(
    summary = summary,
    res = res,
    log = log,
    na_action = na_action
  )
  target_options <- .rateMap_normalize_target_options(
    target = target,
    target_tree = target_tree,
    check = check
  )
  uncertainty_options <- .rateMap_normalize_uncertainty_options(
    uncertainty = uncertainty,
    value_summary = value_summary,
    quantile_probs = quantile_probs,
    highest_density_interval_prob = highest_density_interval_prob
  )
  color_options <- .rateMap_normalize_color_options(
    ncolors = ncolors,
    palette = palette,
    reverse_palette = reverse_palette,
    color_mode = color_mode,
    n_categories = n_categories,
    category_bin_method = category_bin_method,
    category_breaks = category_breaks,
    category_labels = category_labels,
    legend_title = legend_title,
    log = summary_options$log
  )
  compute_options <- .rateMap_normalize_compute_options(
    workers = workers,
    future_strategy = future_strategy,
    future_seed = future_seed,
    future_scheduling = future_scheduling,
    future_chunk_size = future_chunk_size,
    progress = progress
  )
  rate_flags <- .rateMap_normalize_rate_flags(rate_flags)

  structure(
    list(
      summary = summary_options,
      target = target_options,
      uncertainty = uncertainty_options,
      color = color_options,
      compute = compute_options,
      rate_flags = rate_flags,
      flat = list(
        res = summary_options$res,
        check = target_options$check,
        workers = compute_options$workers,
        future_strategy = compute_options$future_strategy,
        future_seed = compute_options$future_seed,
        future_scheduling = compute_options$future_scheduling,
        future_chunk_size = compute_options$future_chunk_size,
        progress = compute_options$progress,
        log = summary_options$log,
        ncolors = color_options$ncolors,
        palette = color_options$palette,
        reverse_palette = color_options$reverse_palette,
        color_mode = color_options$mode,
        n_categories = color_options$n_categories,
        category_bin_method = color_options$category_bin_method,
        category_breaks = color_options$category_breaks,
        category_labels = color_options$category_labels,
        legend_title = color_options$title,
        na_action = summary_options$na_action,
        rate_flags = rate_flags,
        summary = summary_options$mode,
        target = target_options$selection,
        target_tree = target_options$tree,
        target_mode = target_options$mode,
        uncertainty = uncertainty_options$enabled,
        value_summary = uncertainty_options$value_summary,
        quantile_probs = uncertainty_options$quantile_probs,
        highest_density_interval_prob = uncertainty_options$highest_density_interval_prob
      )
    ),
    class = c("rateMap_options", "list")
  )
}

.rateMap_normalize_options <- function(res,
                                       workers,
                                       future_strategy,
                                       future_seed,
                                       future_scheduling,
                                       future_chunk_size,
                                       progress,
                                       na_action,
                                       summary,
                                       target,
                                       check,
                                       value_summary,
                                       uncertainty,
                                       quantile_probs,
                                       highest_density_interval_prob,
                                       log,
                                       target_tree,
                                       rate_flags) {
  grouped <- .rateMap_normalize_grouped_options(
    res = res,
    check = check,
    workers = workers,
    future_strategy = future_strategy,
    future_seed = future_seed,
    future_scheduling = future_scheduling,
    future_chunk_size = future_chunk_size,
    progress = progress,
    log = log,
    na_action = na_action,
    rate_flags = rate_flags,
    summary = summary,
    target = target,
    target_tree = target_tree,
    uncertainty = uncertainty,
    value_summary = value_summary,
    quantile_probs = quantile_probs,
    highest_density_interval_prob = highest_density_interval_prob
  )
  flat <- grouped$flat

  list(
    res = flat$res,
    ncolors = flat$ncolors,
    palette = flat$palette,
    reverse_palette = flat$reverse_palette,
    workers = flat$workers,
    future_strategy = flat$future_strategy,
    future_seed = flat$future_seed,
    future_scheduling = flat$future_scheduling,
    future_chunk_size = flat$future_chunk_size,
    progress = flat$progress,
    na_action = flat$na_action,
    summary = flat$summary,
    target = flat$target,
    target_mode = flat$target_mode,
    check_mode = flat$check,
    rate_flags = flat$rate_flags,
    value_summary = flat$value_summary,
    color_mode = flat$color_mode,
    n_categories = flat$n_categories,
    category_bin_method = flat$category_bin_method,
    category_breaks = flat$category_breaks,
    category_labels = flat$category_labels,
    uncertainty = flat$uncertainty,
    quantile_probs = flat$quantile_probs,
    highest_density_interval_prob = flat$highest_density_interval_prob,
    legend_title = flat$legend_title,
    log = flat$log,
    groups = grouped
  )
}

.rateMap_prepare_inputs <- function(fits,
                                    tree_fun,
                                    param_fun,
                                    log,
                                    na_action,
                                    weights,
                                    target,
                                    target_tree,
                                    check_mode) {
  if (.rateMap_is_single_fit(fits)) {
    fits <- list(fits)
  }
  if (!is.list(fits) || length(fits) == 0L) {
    stop("'fits' must be a non-empty list of fitted run objects.")
  }

  default_tree_fun <- is.null(tree_fun)
  tree_fun <- if (isTRUE(default_tree_fun)) .rateMap_extract_tree else tree_fun
  param_fun <- if (is.null(param_fun)) .rateMap_extract_param else param_fun

  trees <- lapply(fits, tree_fun)
  params <- lapply(fits, param_fun)

  validation <- vapply(
    seq_along(fits),
    function(i) {
      if (isTRUE(default_tree_fun) &&
          .rateMap_is_bifrost_search(fits[[i]]) &&
          is.null(fits[[i]]$tree_no_uncertainty_untransformed)) {
        return(paste0(
          "Fit ", i,
          " is a bifrost_search object but is missing ",
          "'tree_no_uncertainty_untransformed'. rateMap() uses the ",
          "original-length tree by default; use rateMapControl(tree_fun = ...) ",
          "to supply another mapped tree."
        ))
      }
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
    weights = weights
  )

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

  list(
    fits = fits,
    trees = trees,
    params = params,
    target_tree = target_tree,
    edge_matches = .rateMap_match_edges(target_tree, trees),
    target_clade_keys = .rateMap_edge_clade_keys(target_tree),
    keep = keep,
    omitted = omitted,
    resolved_weights = resolved_weights
  )
}

.rateMap_compute_edge <- function(i,
                                  heights,
                                  steps,
                                  summary,
                                  trees,
                                  params,
                                  edge_matches,
                                  keep,
                                  fit_weights,
                                  uncertainty,
                                  value_summary,
                                  quantile_probs,
                                  highest_density_interval_prob,
                                  progressor = NULL,
                                  tol = 1e-10) {
  mids <- if (identical(summary, "interval")) {
    intersect(which(steps > heights[i, 1L]), which(steps < heights[i, 2L]))
  } else {
    integer()
  }
  yy <- cbind(
    c(heights[i, 1L], steps[mids]),
    c(steps[mids], heights[i, 2L])
  ) - heights[i, 1L]
  lengths_i <- yy[, 2L] - yy[, 1L]

  values_by_fit <- matrix(
    NA_real_,
    nrow = nrow(yy),
    ncol = length(trees),
    dimnames = list(NULL, paste0("fit_", keep))
  )
  target_length <- heights[i, 2L] - heights[i, 1L]

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

      values_by_fit[k, j] <- .rateMap_map_value(
        map = map_i,
        fit_rates = fit_rates,
        start = source_start,
        end = source_end,
        tol = tol
      )
    }
  }

  uncertainty_i <- if (isTRUE(uncertainty) || identical(value_summary, "median")) {
    .rateMap_summarize_run_values(
      values = values_by_fit,
      weights = fit_weights,
      quantile_probs = quantile_probs,
      highest_density_interval_prob = highest_density_interval_prob
    )
  } else {
    NULL
  }
  values_i <- if (identical(value_summary, "median")) {
    uncertainty_i$median
  } else {
    .rateMap_row_means(values_by_fit, fit_weights)
  }

  if (!is.null(progressor)) {
    progressor(sprintf("edge %d", i))
  }

  list(
    values = values_i,
    run_values = if (isTRUE(uncertainty)) values_by_fit else NULL,
    uncertainty = if (isTRUE(uncertainty)) uncertainty_i else NULL,
    lengths = lengths_i,
    depth_start = heights[i, 1L] + yy[, 1L],
    depth_end = heights[i, 1L] + yy[, 2L]
  )
}

.rateMap_compute_edges <- function(target_tree,
                                   trees,
                                   params,
                                   edge_matches,
                                   keep,
                                   weights,
                                   summary,
                                   res,
                                   uncertainty,
                                   value_summary,
                                   quantile_probs,
                                   highest_density_interval_prob,
                                   progress,
                                   workers,
                                   future_strategy,
                                   future_seed,
                                   future_scheduling,
                                   future_chunk_size) {
  heights <- phytools::nodeHeights(target_tree)
  max_height <- max(heights)
  steps <- if (identical(summary, "interval")) {
    (0:res / res) * max_height
  } else {
    numeric()
  }
  edge_indices <- seq_len(nrow(target_tree$edge))

  if (!is.null(workers)) {
    old_plan <- future::plan()
    on.exit(future::plan(old_plan), add = TRUE)
    strategy_fun <- get(future_strategy, envir = asNamespace("future"), mode = "function")
    future::plan(strategy_fun, workers = workers)
  }

  compute_all_edges <- function(progressor = NULL) {
    future.apply::future_lapply(
      edge_indices,
      function(i) {
        .rateMap_compute_edge(
          i = i,
          heights = heights,
          steps = steps,
          summary = summary,
          trees = trees,
          params = params,
          edge_matches = edge_matches,
          keep = keep,
          fit_weights = weights,
          uncertainty = uncertainty,
          value_summary = value_summary,
          quantile_probs = quantile_probs,
          highest_density_interval_prob = highest_density_interval_prob,
          progressor = progressor
        )
      },
      future.seed = future_seed,
      future.scheduling = future_scheduling,
      future.chunk.size = future_chunk_size
    )
  }

  edge_results <- if (isTRUE(progress)) {
    progressr::with_progress(
      {
        p <- progressr::progressor(along = edge_indices)
        compute_all_edges(p)
      },
      handlers = progressr::handler_txtprogressbar(),
      enable = TRUE
    )
  } else {
    compute_all_edges()
  }

  list(
    edge_results = edge_results,
    heights = heights
  )
}

.rateMap_build_result <- function(target_tree,
                                  edge_results,
                                  edge_matches,
                                  target_clade_keys,
                                  keep,
                                  resolved_weights,
                                  omitted,
                                  n_fits,
                                  summary,
                                  uncertainty,
                                  value_summary,
                                  quantile_probs,
                                  highest_density_interval_prob,
                                  target_mode,
                                  check_mode,
                                  palette,
                                  reverse_palette,
                                  color_mode,
                                  ncolors,
                                  n_categories,
                                  category_breaks,
                                  category_labels,
                                  category_bin_method,
                                  log,
                                  legend_title,
                                  rate_flags) {
  interval_values <- lapply(edge_results, `[[`, "values")
  interval_lengths <- lapply(edge_results, `[[`, "lengths")
  run_values <- if (isTRUE(uncertainty)) {
    lapply(edge_results, `[[`, "run_values")
  } else {
    NULL
  }
  rate_flag_info <- .rateMap_compute_rate_flags(
    values_by_edge = interval_values,
    summary = summary,
    log = log,
    rate_flags = rate_flags
  )
  display_flags <- if (identical(color_mode, "category") &&
                       isTRUE(rate_flag_info$diagnostics$enabled)) {
    rate_flag_info$flags_by_edge
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
    category_bin_method = category_bin_method,
    display_flags_by_edge = display_flags,
    special_categories = list(
      zero_label = rate_flags$zero_label,
      high_label = rate_flags$high_label,
      zero_color = rate_flags$zero_color,
      high_color = rate_flags$high_color
    )
  )

  tree <- target_tree
  tree$maps <- vector("list", nrow(tree$edge))
  for (i in seq_len(nrow(tree$edge))) {
    tree$maps[[i]] <- interval_lengths[[i]]
    names(tree$maps[[i]]) <- color_map$states_by_edge[[i]]
  }

  tree$mapped.edge <- .rateMap_make_mapped_edge(tree$edge, tree$maps)
  class(tree) <- c("simmap", setdiff(class(tree), "simmap"))
  attr(tree, "map.order") <- "right-to-left"

  edge_indices <- seq_len(nrow(tree$edge))
  interval_counts <- vapply(interval_lengths, length, integer(1))
  interval_offsets <- c(0L, cumsum(interval_counts))
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
      if (isTRUE(rate_flag_info$diagnostics$enabled)) {
        idx <- seq_len(nrow(base_row)) + interval_offsets[i]
        base_row$rate_for_flagging <- rate_flag_info$rate_for_flagging[idx]
        base_row$rate_flag <- rate_flag_info$rate_flag[idx]
        base_row$rate_flag_source <- "value"
        base_row$is_near_zero <- base_row$rate_flag == rate_flags$zero_label
        base_row$is_high_outlier <- base_row$rate_flag == rate_flags$high_label
      }
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
  rate_categories <- .rateMap_add_category_summaries(
    rate_categories = color_map$rate_categories,
    intervals = intervals,
    summary = summary
  )

  weight_table <- data.frame(
    input_index = keep,
    weight = resolved_weights$weights,
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
    rate_categories = rate_categories,
    clade_key = target_clade_keys,
    edge_matches = do.call(cbind, edge_matches),
    summary = summary,
    uncertainty = uncertainty,
    value_summary = value_summary,
    quantile_probs = quantile_probs,
    highest_density_interval_prob = highest_density_interval_prob,
    rate_diagnostics = rate_flag_info$diagnostics,
    rate_flags = rate_flags,
    rate_flag_source = if (isTRUE(rate_flag_info$diagnostics$enabled)) "value" else NA_character_,
    plot_value = "value",
    target = target_mode,
    check = check_mode,
    weights = resolved_weights$weights,
    weight_mode = resolved_weights$mode,
    weight_table = weight_table,
    palette = palette,
    reverse_palette = reverse_palette,
    ncolors = as.integer(ncolors),
    color_mode = color_mode,
    n_categories = n_categories,
    category_breaks = color_map$category_breaks,
    category_labels = color_map$category_labels,
    category_bin_method = color_map$category_bin_method,
    log = isTRUE(log),
    title = legend_title,
    n_fits = n_fits,
    omitted = omitted
  )
  colnames(out$edge_matches) <- paste0("fit_", keep)
  class(out) <- "rateMap"
  out
}

#' Compute Branchwise Rate Maps Across Runs
#'
#' Summarize fitted regime-specific rates across a list of stochastic-map-aware
#' model fits or completed `bifrost_search` results. By default, `rateMap()`
#' follows `bifrost`'s branch-level framing: the first retained tree is used as
#' the plotting scaffold, all retained trees must match in topology and branch
#' lengths, each branch receives one summarized log-rate, and runs are averaged
#' with equal weight.
#'
#' @details
#' **Everyday workflow.** `rateMap()` is a compute function: call it to build a
#' reusable `"rateMap"` object, then call `plot(x, ...)` to choose what to
#' display. Common compute controls are `fits`, `weights`, `uncertainty`,
#' `summary`, `log`, and `value_summary`. Common display controls belong in
#' `plot()` or [rateMapView()], including `value`, `type`, `palette`,
#' `color_mode`, `n_categories`, `show_tip_labels`, and legend sizing. Controls
#' such as `target`, `check`, `res`, `tree_fun`, `param_fun`, and the `future_*`
#' arguments live in [rateMapControl()] because they are mainly advanced tools
#' for same-topology tree samples, interval summaries, custom fit objects, and
#' larger run sets.
#'
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
#' is a convenience for inspecting one fitted model before scaling up to
#' multi-run summaries.
#'
#' `rateMap()` can also summarize same-topology posterior or sensitivity samples
#' where branch lengths differ. In that case, usually supply an explicit
#' `target_tree` and use `control = list(check = "topology")`. The target can
#' be any tree with the same topology and tip labels as the inputs; an MCC tree
#' is only one possible choice. The target tree supplies the plotted topology
#' and branch lengths, while rates are matched from each input tree by
#' descendant-tip clade keys rather than by edge order.
#'
#' **Summary modes.** With the default `summary = "branch"` and `log = TRUE`,
#' each edge receives one length-weighted average log-rate from each run before
#' the across-run summary is computed. This is the natural default for `bifrost`
#' searches because shifts are placed at nodes, so a `bifrost` branch is not
#' expected to change regimes internally. With `summary = "interval"`, each
#' target-tree branch is subdivided by the global depth grid controlled by
#' `control = list(res = ...)`. If source branch lengths differ from the target
#' branch length, source stochastic-map segments are projected onto target
#' intervals by relative position along the matched branch. Interval mode is
#' useful for general stochastic maps that can genuinely change state along a
#' branch.
#'
#' **Log-rate averaging.** When `log = TRUE`, rate parameters are transformed
#' before branch-level, interval-level, and across-fit summaries are computed.
#' With weights `w_i`, the default plotted mean is therefore
#' `sum(w_i * log(rate_i))`, which is the log of a weighted geometric mean. It
#' is not `log(sum(w_i * rate_i))`, the log of a weighted arithmetic mean. Use
#' `log = FALSE` when downstream interpretation requires arithmetic summaries on
#' the original rate scale. Raw fitted rate parameters must be strictly positive
#' in either mode. Negative plotted values are therefore valid in log-rate maps
#' whenever the positive raw fitted rates are less than one.
#'
#' **Tree checks and targets.** In [rateMapControl()], `check = TRUE` is
#' equivalent to `check = "full"` and requires topology and branch lengths to
#' match the target tree. Use `check = "topology"` when all inputs have the same
#' topology and tip labels but may have different branch lengths. `check = FALSE`
#' or `check = "none"` skips the upfront `ape::all.equal.phylo()` check, but
#' every target branch must still be recoverable in every input tree by
#' descendant-tip set. This function is not a mixed-topology posterior summarizer;
#' clades absent from an input tree are treated as an error rather than being
#' marginalized over topology.
#'
#' When `target_tree` is supplied, it is used directly as the plotting scaffold
#' and need not contain SIMMAP maps. When `target_tree = NULL`,
#' `rateMapControl(target = "first")` uses the first retained input tree, and
#' `rateMapControl(target = "mcc")` chooses the retained input tree with the
#' highest sum of log clade credibilities. For truly same-topology inputs, the
#' MCC score is usually tied, so `"mcc"` commonly resolves to the first retained
#' tree. MCC target selection does not make `rateMap()` a mixed-topology
#' summarizer: every target branch must still be present in every retained input.
#' If you already have a preferred consensus, chronogram, MCC, maximum-likelihood,
#' or otherwise curated target tree, pass it with `target_tree`.
#'
#' **Fit weights.** `weights = "equal"` assigns the same weight to each retained
#' fit. `weights = "ic"` computes standard information-criterion weights from
#' each retained fit's `optimal_ic`, requiring all retained fits to share the
#' same non-missing `IC_used`. A numeric `weights` vector can be supplied for
#' custom weighting. Weights are subset to retained fits after
#' `rateMapControl(na_action = "omit")` and are normalized to sum to one. IC
#' weights are descriptive fit-level weights for comparable retained searches;
#' they do not choose a formal threshold or make incomparable searches
#' comparable.
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
#' These are weighted empirical summaries of run-level values, not posterior
#' uncertainty or model-internal uncertainty unless the retained inputs
#' themselves have that interpretation. For posterior tree samples, use
#' `weights = "equal"` if each retained run represents one posterior draw.
#' `value_summary` controls whether the plotted `value` column uses the weighted
#' mean or weighted median. These summaries are computed on the log-rate scale
#' when `log = TRUE`. Highest-density intervals use the same shortest empirical
#' interval calculation for both equal and unequal fit weights. For
#' `weights = "equal"`, `sd` is the ordinary sample standard deviation of
#' retained run-level values. For unequal weights, `sd` is the square root of the
#' normalized weighted variance, `sum(w * (x - mu)^2)`, with weights normalized
#' to sum to one.
#'
#' **Display views.** Returned objects include a default category-style display
#' mapping so they can be plotted immediately. Palette, category-bin, legend
#' title, and alternative-value choices are intentionally controlled by
#' `plot(x, ...)` or [rateMapView()], not by `rateMap()`. For a single `bifrost`
#' search with `log = FALSE`, category display is the closest formal analogue to
#' the illustrative `generateViridisColorScale()` plot. For multi-run summaries,
#' displayed categories are bins for summarized branch values; they should not be
#' read as newly inferred `bifrost` regimes. When `summary = "branch"` and
#' `color_mode = "category"`, the `rate_categories` table also reports
#' branch-level summaries for the plotted values assigned to each bin. These
#' summaries are recomputed by [rateMapView()] or `plot()` whenever the plotted
#' value or category breaks change. They are not computed for `summary =
#' "interval"` because interval rows depend on the plotting grid rather than on
#' discrete biological branch units.
#'
#' **Rate diagnostics.** Branch-summary maps compute rate diagnostics through
#' `rateMapControl(rate_flags = rateMapRateFlags())`. The default
#' `rateMapRateFlags()` setting records the fitted-rate range and fold range but
#' applies no special near-zero or high-outlier rule. Supplying `zero_floor`, equivalently
#' `rateMapRateFlags(method = "floor", zero_floor = ...)`, flags finite rates at
#' or below an explicit manual floor. Use `rateMapRateFlags(method =
#' "tail_cluster")` to apply a deterministic, Otsu-style guarded two-class split
#' on sorted log rates when the near-zero values form a broader separated
#' lower-tail cluster rather than one extreme adjacent gap. The split is kept
#' only when gap, fold-reduction, tail-fraction, and minimum-count guardrails are
#' satisfied. These flags never remove branches and never alter
#' `intervals$value`; they add `rate_flag` metadata, `rate_flag_source`
#' provenance, object-level `rate_diagnostics`, and, in category display mode,
#' special display categories such as `"near-zero"`. Special categories are
#' colored outside the ordered palette, so the remaining regular categories
#' still use the full low-to-high palette range. In category legends, diagnostic
#' cutoffs mark where special categories end and regular bins begin. Set
#' `rate_flags = NULL` to turn off rate-flag metadata and diagnostics entirely.
#'
#' @param fits A completed run, fitted model object, or non-empty list of these
#'   objects. Supported shapes are `bifrost_search` objects, `mvgls` objects,
#'   `list(model = <mvgls>)`, and scratch-style lists with `variables$tree` plus
#'   `param`. Single supported objects are wrapped automatically as a convenience
#'   for one-fit plotting and inspection.
#' @param weights Fit-level weighting mode. `"equal"` gives every retained fit
#'   equal weight. `"ic"` computes standard IC weights from `optimal_ic` and
#'   requires all retained fits to have the same `IC_used`. A numeric vector is
#'   also accepted and treated as custom fit weights.
#' @param uncertainty Logical; if `TRUE`, compute and return across-fit
#'   uncertainty summaries for every summary row: whole branches when
#'   `summary = "branch"` and depth-grid intervals when `summary = "interval"`.
#' @param summary Character; `"branch"` computes one length-weighted value per
#'   edge, while `"interval"` slices branches on a global depth grid.
#' @param log Logical; if `TRUE` (the default), transform extracted rate
#'   parameters with `log()` before branch, interval, and across-fit averaging.
#'   Weighted means are then mean log-rates, equivalent to the log of weighted
#'   geometric mean rates. Set `log = FALSE` to summarize rates in their original
#'   units.
#' @param value_summary Character; central estimate stored in the summary table
#'   column `intervals$value`. `"mean"` uses the weighted mean. `"median"` uses
#'   the weighted median. When `log = TRUE`, both summaries are computed on the
#'   log-rate scale.
#' @param target_tree Optional explicit target tree used as the summary and
#'   plotting scaffold. This may be any target or summary tree with the same
#'   topology and tip labels as the inputs. It does not need to contain
#'   stochastic maps because `rateMap()` replaces maps with display maps in the
#'   returned object. Branch lengths in `target_tree` define the geometry of the
#'   returned and plotted tree.
#' @param workers Optional number of `future` workers to use. If `NULL`, the
#'   current `future::plan()` is used as-is.
#' @param progress Logical; if `TRUE`, display a text progress bar via
#'   [progressr::with_progress()].
#' @param control A `"rateMap_control"` object from [rateMapControl()], or a
#'   named list of advanced control options.
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
#'   plotted depth-grid interval. When branch-level rate diagnostics are
#'   enabled, this table also includes `rate_for_flagging`, `rate_flag`,
#'   `rate_flag_source`, `is_near_zero`, and `is_high_outlier`.}
#'   \item{`rate_categories`}{Data frame describing discrete rate categories
#'   when `color_mode = "category"`; otherwise `NULL`. With `summary =
#'   "branch"`, this table also includes bin-level summaries of the plotted
#'   branch values, including `n_branches`, `value_mean`, `value_median`,
#'   `value_min`, `value_max`, `value_sd`, and `total_branch_length`.}
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
#'   \item{`rate_diagnostics`}{List summarizing rate-flag settings, counts,
#'   detected tail cutoffs, and fold-rate ranges with and without flagged
#'   branches.}
#'   \item{`rate_flags`}{The normalized `"rateMap_rate_flags"` control object
#'   used for rate diagnostics.}
#'   \item{`rate_flag_source`}{Character name of the rate-valued column used to
#'   compute or preserve `rate_flag` metadata, or `NA` when diagnostics are
#'   disabled.}
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
#'   \item{`ncolors`}{Stored continuous-ramp resolution used when recoloring
#'   with `color_mode = "continuous"` and no explicit `ncolors`.}
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
#' @seealso [plot.rateMap()], [rateMapView()], [rateMapControl()],
#'   [phytools::densityMap()], [phytools::plotSimmap()]
#'
#' @references
#' Paradis, E., and Schliep, K. (2019). ape 5.0: an environment for modern
#' phylogenetics and evolutionary analyses in R. *Bioinformatics*, 35, 526-528.
#' \doi{10.1093/bioinformatics/bty633}
#'
#' Clavel, J., Aristide, L., and Morlon, H. (2019). A penalized likelihood
#' framework for high-dimensional phylogenetic comparative methods and an
#' application to new-world monkeys brain evolution. *Systematic Biology*, 68,
#' 93-116. \doi{10.1093/sysbio/syy045}
#'
#' Revell, L. J. (2013). Two new graphical methods for mapping trait evolution
#' on phylogenies. *Methods in Ecology and Evolution*, 4, 754-759.
#'
#' Revell, L. J. (2024). phytools 2.0: an updated R ecosystem for phylogenetic
#' comparative methods (and other things). *PeerJ*, 12, e16505.
#' \doi{10.7717/peerj.16505}
#'
#' @examples
#' \dontrun{
#' # A list of completed bifrost searches can be summarized directly:
#' rm_obj <- rateMap(list(search_a, search_b, search_c))
#' plot(rm_obj, type = "arc", show_tip_labels = FALSE)
#'
#' # A single completed bifrost search can be inspected with the same API:
#' one_search_map <- rateMap(search_a)
#' plot(one_search_map, type = "arc", show_tip_labels = FALSE)
#'
#' # Scratch-style lists remain supported:
#' scratch_fit <- list(
#'   variables = list(tree = mapped_tree),
#'   param = c("0" = 0.12, "1" = 0.45)
#' )
#' rateMap(list(scratch_fit))
#'
#' # Use bifrost model-level IC weights across comparable sensitivity runs:
#' rm_ic <- rateMap(list(search_a, search_b), weights = "ic")
#'
#' # Branch-level, discrete-category maps are the bifrost default:
#' branch_rates <- rateMap(list(search_a, search_b))
#' branch_rates$rate_categories
#'
#' # To mimic the old jaw-shape preview style for a single search, use raw
#' # fitted rates and category colors:
#' jaw_view <- rateMapView(rateMap(search_a, log = FALSE), palette = viridis::viridis)
#'
#' # Category bins use pretty breaks by default. Use equal-width bins or custom
#' # boundaries when those are easier to compare across figures:
#' equal_bin_rates <- rateMapView(
#'   branch_rates,
#'   n_categories = 5,
#'   category_bin_method = "equal"
#' )
#' custom_bin_rates <- rateMapView(
#'   branch_rates,
#'   category_breaks = c(-4, -2, 0, 2),
#'   category_labels = c("slow", "middle", "fast")
#' )
#'
#' # Continuous interval maps remain available for stochastic maps with
#' # along-branch changes:
#' interval_rates <- rateMap(
#'   list(search_a, search_b),
#'   summary = "interval",
#'   control = list(res = 200)
#' )
#' plot(interval_rates, color_mode = "continuous")
#'
#' # Custom fit weights are normalized internally:
#' custom_weighted <- rateMap(
#'   list(search_a, search_b, search_c),
#'   weights = c(2, 1, 1),
#'   summary = "branch"
#' )
#'
#' # Posterior trees with the same topology but different branch lengths can be
#' # summarized on any explicit target or summary tree:
#' posterior_target_rates <- rateMap(
#'   posterior_fit_list,
#'   target_tree = summary_tree,
#'   summary = "branch",
#'   weights = "equal",
#'   uncertainty = TRUE,
#'   control = list(check = "topology")
#' )
#' plot(posterior_target_rates, value = "sd", type = "arc")
#'
#' # Or choose the retained input tree with the highest summed log clade
#' # credibility as the plotting scaffold:
#' posterior_mcc_rates <- rateMap(
#'   posterior_fit_list,
#'   summary = "branch",
#'   control = list(check = "topology", target = "mcc")
#' )
#'
#' # Large sensitivity sets can be explicitly subsampled before plotting.
#' # By default, rates are mapped on the log scale:
#' set.seed(1)
#' idx <- sample(seq_along(fit_list), size = 1000)
#' rm_sub <- rateMap(
#'   fit_list[idx],
#'   workers = 8,
#'   control = list(res = 100, future_strategy = "multisession")
#' )
#'
#' # Use log = FALSE only when a raw-rate scale is preferred:
#' rm_raw <- rateMap(fit_list[idx], log = FALSE)
#' plot(
#'   rm_sub,
#'   type = "arc",
#'   show_tip_labels = FALSE,
#'   lwd = 1,
#'   palette = c("lightblue", "blue", "pink", "red")
#' )
#' }
#'
#' @export
rateMap <- function(
  fits,
  weights = c("equal", "ic"),
  uncertainty = FALSE,
  summary = c("branch", "interval"),
  log = TRUE,
  value_summary = c("mean", "median"),
  target_tree = NULL,
  workers = NULL,
  progress = TRUE,
  control = rateMapControl()
) {
  .rateMap_require_packages()
  control <- .rateMap_normalize_control(control)
  opts <- .rateMap_normalize_options(
    res = control$res,
    workers = workers,
    future_strategy = control$future_strategy,
    future_seed = control$future_seed,
    future_scheduling = control$future_scheduling,
    future_chunk_size = control$future_chunk_size,
    progress = progress,
    na_action = control$na_action,
    summary = summary,
    target = control$target,
    check = control$check,
    value_summary = value_summary,
    uncertainty = uncertainty,
    quantile_probs = control$quantile_probs,
    highest_density_interval_prob = control$highest_density_interval_prob,
    log = log,
    target_tree = target_tree,
    rate_flags = control$rate_flags
  )

  input <- .rateMap_prepare_inputs(
    fits = fits,
    tree_fun = control$tree_fun,
    param_fun = control$param_fun,
    log = opts$log,
    na_action = opts$na_action,
    weights = weights,
    target = opts$target,
    target_tree = target_tree,
    check_mode = opts$check_mode
  )

  edge_computation <- .rateMap_compute_edges(
    target_tree = input$target_tree,
    trees = input$trees,
    params = input$params,
    edge_matches = input$edge_matches,
    keep = input$keep,
    weights = input$resolved_weights$weights,
    summary = opts$summary,
    res = opts$res,
    uncertainty = opts$uncertainty,
    value_summary = opts$value_summary,
    quantile_probs = opts$quantile_probs,
    highest_density_interval_prob = opts$highest_density_interval_prob,
    progress = opts$progress,
    workers = opts$workers,
    future_strategy = opts$future_strategy,
    future_seed = opts$future_seed,
    future_scheduling = opts$future_scheduling,
    future_chunk_size = opts$future_chunk_size
  )

  out <- .rateMap_build_result(
    target_tree = input$target_tree,
    edge_results = edge_computation$edge_results,
    edge_matches = input$edge_matches,
    target_clade_keys = input$target_clade_keys,
    keep = input$keep,
    resolved_weights = input$resolved_weights,
    omitted = input$omitted,
    n_fits = length(input$trees),
    summary = opts$summary,
    uncertainty = opts$uncertainty,
    value_summary = opts$value_summary,
    quantile_probs = opts$quantile_probs,
    highest_density_interval_prob = opts$highest_density_interval_prob,
    target_mode = opts$target_mode,
    check_mode = opts$check_mode,
    palette = opts$palette,
    reverse_palette = opts$reverse_palette,
    color_mode = opts$color_mode,
    ncolors = opts$ncolors,
    n_categories = opts$n_categories,
    category_breaks = opts$category_breaks,
    category_labels = opts$category_labels,
    category_bin_method = opts$category_bin_method,
    log = opts$log,
    legend_title = opts$legend_title,
    rate_flags = opts$rate_flags
  )

  out
}
