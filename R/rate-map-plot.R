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
  if (isTRUE(.rateMap_has_special_categories(x))) {
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

.rateMap_has_special_categories <- function(x) {
  categories <- x$rate_categories
  flags <- x$rate_flags
  labels <- if (is.null(categories)) character() else categories$rate_category
  any(labels %in% c(flags$zero_label, flags$high_label))
}

.rateMap_rate_marker_to_value <- function(rate, log) {
  if (!is.numeric(rate) || length(rate) != 1L || !is.finite(rate)) {
    return(NA_real_)
  }
  if (isTRUE(log)) {
    if (rate <= 0) {
      return(NA_real_)
    }
    return(base::log(rate))
  }
  rate
}

.rateMap_category_marker_values <- function(x) {
  categories <- x$rate_categories
  flags <- x$rate_flags
  diagnostics <- x$rate_diagnostics
  if (is.null(categories) || is.null(flags) || is.null(diagnostics)) {
    return(numeric())
  }

  markers <- numeric()
  labels <- categories$rate_category
  if (flags$zero_label %in% labels) {
    marker <- .rateMap_rate_marker_to_value(diagnostics$zero_floor, x$log)
    if (!is.finite(marker)) {
      marker <- .rateMap_rate_marker_to_value(diagnostics$zero_cluster_cutoff_rate, x$log)
    }
    if (is.finite(marker)) {
      markers[flags$zero_label] <- marker
    }
  }
  if (flags$high_label %in% labels) {
    marker <- .rateMap_rate_marker_to_value(diagnostics$high_cluster_cutoff_rate, x$log)
    if (is.finite(marker)) {
      markers[flags$high_label] <- marker
    }
  }

  markers
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

  lims <- .rateMap_category_lims(x)

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
  markers <- .rateMap_category_marker_values(x)
  draw_special_scale <- isTRUE(.rateMap_has_special_categories(x)) &&
    length(markers) > 0L &&
    is.finite(diff(x$lims)) &&
    diff(x$lims) > 0

  if (isTRUE(draw_special_scale)) {
    full_lims <- x$lims
    markers <- pmin(pmax(markers, full_lims[1L]), full_lims[2L])
    value_to_x <- function(value) {
      frac <- (value - full_lims[1L]) / diff(full_lims)
      if (identical(direction, "leftwards")) {
        x1 - frac * legend
      } else {
        x0 + frac * legend
      }
    }
    rect_by_value <- function(left, right, col) {
      xs <- sort(c(value_to_x(left), value_to_x(right)))
      graphics::rect(xs[1L], y0, xs[2L], y1, col = col, border = NA)
    }

    flags <- x$rate_flags
    labels <- categories$rate_category
    regular <- !labels %in% c(flags$zero_label, flags$high_label)
    regular_cols <- categories$color[regular]
    zero_marker <- markers[flags$zero_label]
    high_marker <- markers[flags$high_label]
    regular_min <- if (is.finite(zero_marker)) zero_marker else full_lims[1L]
    regular_max <- if (is.finite(high_marker)) high_marker else full_lims[2L]

    if (is.finite(zero_marker)) {
      rect_by_value(full_lims[1L], zero_marker, flags$zero_color)
    }
    if (length(regular_cols) > 0L && regular_max > regular_min) {
      regular_x <- c(value_to_x(regular_min), value_to_x(regular_max))
      draw_cols <- if (regular_x[1L] <= regular_x[2L]) regular_cols else rev(regular_cols)
      xs <- seq(min(regular_x), max(regular_x), length.out = length(draw_cols) + 1L)
      graphics::rect(
        xs[-length(xs)],
        y0,
        xs[-1L],
        y1,
        col = unname(draw_cols),
        border = NA
      )
    }
    if (is.finite(high_marker)) {
      rect_by_value(high_marker, full_lims[2L], flags$high_color)
    }

    label_values <- unique(c(full_lims[1L], unname(markers), full_lims[2L]))
    label_values <- label_values[is.finite(label_values)]
    tick_x <- value_to_x(label_values)
    graphics::segments(tick_x, y0, tick_x, y1, col = "grey20", lwd = 0.5)
    graphics::text(
      tick_x,
      rep(y1 + label_gap, length(tick_x)),
      labels = formatC(label_values, digits = digits, format = "fg"),
      cex = fsize,
      adj = c(0.5, 0.5)
    )
  } else {
    bar_cols <- x$cols
    if (length(bar_cols) == 1L) {
      bar_cols <- rep(bar_cols, 2L)
    }
    if (identical(direction, "leftwards")) {
      bar_cols <- rev(bar_cols)
      lims <- rev(lims)
    }
    xs <- seq(x0, x1, length.out = length(bar_cols) + 1L)
    graphics::rect(
      xs[-length(xs)],
      y0,
      xs[-1L],
      y1,
      col = unname(bar_cols),
      border = NA
    )
    graphics::text(
      c(x0, x1),
      rep(y1 + label_gap, 2L),
      labels = formatC(lims, digits = digits, format = "fg"),
      cex = fsize,
      adj = c(0.5, 0.5)
    )
  }
  if (isTRUE(outline)) {
    graphics::rect(x0, y0, x1, y1, border = "grey30")
  }
  graphics::text(
    mean(c(x0, x1)),
    y1 + title_gap,
    labels = x$title,
    cex = fsize
  )
  invisible(NULL)
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

.rateMap_stored_ncolors <- function(x) {
  ncolors <- x$ncolors
  if (is.numeric(ncolors) && length(ncolors) == 1L &&
      is.finite(ncolors) && ncolors >= 2) {
    return(as.integer(ncolors))
  }
  256L
}

.rateMap_validate_ncolors <- function(ncolors) {
  if (!is.numeric(ncolors) || length(ncolors) != 1L ||
      !is.finite(ncolors) || ncolors < 2) {
    stop("'ncolors' must be a finite number >= 2.")
  }
  as.integer(ncolors)
}

.rateMap_recolor <- function(x,
                             value,
                             palette = NULL,
                             reverse_palette = NULL,
                             color_mode = NULL,
                             n_categories = NULL,
                             ncolors = NULL,
                             category_breaks = NULL,
                             category_labels = NULL,
                             category_bin_method = NULL,
                             legend_title = NULL) {
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
  flat_vals <- unlist(vals, use.names = FALSE)
  if (!any(is.finite(flat_vals))) {
    stop(
      "Column '", value, "' contains no finite values to map. ",
      "Choose a 'rateMap' interval column with at least one finite value."
    )
  }
  rate_value <- value %in% c("value", "mean", "median")

  selected_palette <- if (is.null(palette)) x$palette else palette
  selected_reverse <- if (is.null(reverse_palette)) {
    isTRUE(x$reverse_palette)
  } else {
    reverse_palette
  }
  selected_color_mode <- if (is.null(color_mode)) {
    x$color_mode
  } else {
    match.arg(color_mode, c("continuous", "category"))
  }
  selected_n_categories <- if (is.null(n_categories)) {
    x$n_categories
  } else {
    .rateMap_validate_n_categories(n_categories)
  }
  stored_ncolors <- .rateMap_stored_ncolors(x)
  selected_ncolors <- if (is.null(ncolors)) {
    stored_ncolors
  } else {
    .rateMap_validate_ncolors(ncolors)
  }
  selected_category_bin_method <- if (is.null(category_bin_method)) {
    x$category_bin_method
  } else {
    match.arg(category_bin_method, c("pretty", "equal"))
  }
  current_plot_value <- x$plot_value
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
    selected_ncolors
  } else {
    max(length(x$cols), 2L)
  }
  rate_flag_info <- if (isTRUE(rate_value)) {
    .rateMap_compute_rate_flags(
      values_by_edge = vals,
      summary = x$summary,
      log = x$log,
      rate_flags = x$rate_flags
    )
  } else {
    NULL
  }
  display_flags <- if (identical(selected_color_mode, "category") &&
                       !is.null(rate_flag_info) &&
                       isTRUE(rate_flag_info$diagnostics$enabled)) {
    rate_flag_info$flags_by_edge
  } else {
    NULL
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
    category_bin_method = selected_category_bin_method,
    display_flags_by_edge = display_flags,
    special_categories = list(
      zero_label = x$rate_flags$zero_label,
      high_label = x$rate_flags$high_label,
      zero_color = x$rate_flags$zero_color,
      high_color = x$rate_flags$high_color
    )
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
  if (!is.null(rate_flag_info) && isTRUE(rate_flag_info$diagnostics$enabled)) {
    x$intervals$rate_for_flagging <- rate_flag_info$rate_for_flagging
    x$intervals$rate_flag <- rate_flag_info$rate_flag
    x$intervals$rate_flag_source <- value
    x$intervals$is_near_zero <- x$intervals$rate_flag == x$rate_flags$zero_label
    x$intervals$is_high_outlier <- x$intervals$rate_flag == x$rate_flags$high_label
    x$rate_diagnostics <- rate_flag_info$diagnostics
    x$rate_flag_source <- value
  } else if (!is.null(rate_flag_info)) {
    x$rate_diagnostics <- rate_flag_info$diagnostics
  } else if ("rate_flag" %in% names(x$intervals)) {
    x$intervals$rate_flag_source <- x$rate_flag_source
  }
  if (identical(selected_color_mode, "category")) {
    x$intervals$rate_category <- unlist(color_map$states_by_edge, use.names = FALSE)
  } else if ("rate_category" %in% names(x$intervals)) {
    x$intervals$rate_category <- NULL
  }
  x$rate_categories <- .rateMap_add_category_summaries(
    rate_categories = color_map$rate_categories,
    intervals = x$intervals,
    summary = x$summary
  )
  x$palette <- selected_palette
  x$reverse_palette <- selected_reverse
  x$color_mode <- selected_color_mode
  x$ncolors <- if (identical(selected_color_mode, "continuous")) {
    selected_ncolors
  } else {
    stored_ncolors
  }
  x$n_categories <- selected_n_categories
  x$category_breaks <- color_map$category_breaks
  x$category_labels <- color_map$category_labels
  x$category_bin_method <- color_map$category_bin_method
  x$plot_value <- value
  default_title <- switch(
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
  x$title <- if (is.null(legend_title)) default_title else legend_title
  x
}

#' Create a Display View of a `rateMap` Object
#'
#' Recompute the display mapping for a computed `"rateMap"` object without
#' drawing it. This is optional in the usual `rateMap()` then `plot()` workflow;
#' use it when you want to save or inspect category tables, color bins, or an
#' uncertainty-valued map object and then reuse the same display mapping in one
#' or more plots. The numeric branch or interval summaries computed by
#' [rateMap()] are not recomputed.
#'
#' @param x An object of class `"rateMap"` returned by [rateMap()].
#' @param value Name of the numeric column in `x$intervals` to map to colors.
#'   The default, `"value"`, uses the central summary computed by [rateMap()].
#' @param palette Optional palette override. This can be an `hcl.colors()`
#'   palette name, a vector of colors, or a palette function.
#' @param reverse_palette Optional logical override for palette reversal.
#' @param color_mode Optional color-mode override. Use `"category"` for
#'   discrete ordered bins or `"continuous"` for a continuous ramp.
#' @param n_categories Optional target category count for
#'   `color_mode = "category"`.
#' @param category_bin_method Optional automatic category-binning method:
#'   `"pretty"` or `"equal"`.
#' @param category_breaks Optional strictly increasing category boundaries.
#' @param category_labels Optional labels for displayed categories.
#' @param ncolors Optional number of colors for continuous ramps. If omitted,
#'   the stored `x$ncolors` value is used, falling back to `256` for older
#'   objects.
#' @param legend_title Optional legend title stored on the returned object.
#'
#' @return A `"rateMap"` object with updated tree maps, color palette,
#'   `intervals$value`, `rate_categories`, and legend title. For branch-summary
#'   category views, `rate_categories` includes bin-level summaries of the
#'   plotted branch values and is recomputed whenever the view changes. For
#'   branch-summary maps with active rate diagnostics, diagnostic flags are
#'   recomputed for rate-valued views (`"value"`, `"mean"`, or `"median"`) and
#'   preserved as metadata for non-rate views such as `"sd"`. When preserved for
#'   a non-rate view, `rate_flag_source` identifies the
#'   rate-valued column that the flags classify; the flags do not classify the
#'   displayed uncertainty value. The selected `value` column must contain at
#'   least one finite value; uncertainty columns such as `"sd"` may be all `NA`
#'   for single-fit objects.
#'
#' @examples
#' \dontrun{
#' rm_obj <- rateMap(fits, uncertainty = TRUE)
#' display_obj <- rateMapView(
#'   rm_obj,
#'   palette = "Viridis",
#'   n_categories = 5,
#'   legend_title = "Mean log fitted rate"
#' )
#' plot(display_obj, type = "arc", show_tip_labels = FALSE)
#'
#' sd_obj <- rateMapView(
#'   rm_obj,
#'   value = "sd",
#'   palette = "Inferno",
#'   color_mode = "continuous"
#' )
#' }
#'
#' @export
rateMapView <- function(x,
                        value = "value",
                        palette = NULL,
                        reverse_palette = NULL,
                        color_mode = NULL,
                        n_categories = NULL,
                        category_bin_method = NULL,
                        category_breaks = NULL,
                        category_labels = NULL,
                        ncolors = NULL,
                        legend_title = NULL) {
  if (!inherits(x, "rateMap")) {
    stop("'x' must be an object of class 'rateMap'.")
  }

  .rateMap_recolor(
    x,
    value = value,
    palette = palette,
    reverse_palette = reverse_palette,
    color_mode = color_mode,
    n_categories = n_categories,
    ncolors = ncolors,
    category_bin_method = category_bin_method,
    category_breaks = category_breaks,
    category_labels = category_labels,
    legend_title = legend_title
  )
}

.plot_rate_map <- function(x,
                           value = "value",
                           palette = NULL,
                           reverse_palette = NULL,
                           color_mode = NULL,
                           n_categories = NULL,
                           category_bin_method = NULL,
                           category_breaks = NULL,
                           category_labels = NULL,
                           ncolors = NULL,
                           legend_title = NULL,
                           legend = NULL,
                           fsize = NULL,
                           tip_fsize = NULL,
                           legend_fsize = NULL,
                           ftype = NULL,
                           show_tip_labels = TRUE,
                           outline = FALSE,
                           lwd = 3,
                           type = c("phylogram", "fan", "arc"),
                           mar = rep(0.3, 4),
                           direction = "rightwards",
                           offset = NULL,
                           xlim = NULL,
                           ylim = NULL,
                           hold = TRUE,
                           underscore = FALSE,
                           arc_height = 2,
                           legend_digits = NULL,
                           ...) {
  dots <- list(...)
  if (!inherits(x, "rateMap")) {
    stop("'x' must be an object of class 'rateMap'.")
  }

  .rateMap_validate_dots(dots, character(0))

  if (!identical(value, "value") ||
      !is.null(palette) ||
      !is.null(reverse_palette) ||
      !is.null(color_mode) ||
      !is.null(n_categories) ||
      !is.null(ncolors) ||
      !is.null(category_bin_method) ||
      !is.null(category_breaks) ||
      !is.null(category_labels) ||
      !is.null(legend_title)) {
    x <- .rateMap_recolor(
      x,
      value = value,
      palette = palette,
      reverse_palette = reverse_palette,
      color_mode = color_mode,
      n_categories = n_categories,
      ncolors = ncolors,
      category_bin_method = category_bin_method,
      category_breaks = category_breaks,
      category_labels = category_labels,
      legend_title = legend_title
    )
  }

  tree <- x$tree
  cols <- x$cols
  H <- phytools::nodeHeights(tree)
  legend_digits <- if (is.null(legend_digits)) .rateMap_legend_digits(x$lims) else legend_digits

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

#' Plot `rateMap` Objects
#'
#' Render a computed `"rateMap"` object. The usual workflow is
#' `x <- rateMap(...)` followed by `plot(x, ...)`.
#'
#' The plot method uses [phytools::plotSimmap()] and draws either a continuous
#' color-bar legend or a segmented rate-category color bar.
#' The plotting controls intentionally mirror the `phytools` density-map
#' plotting style for phylogram, fan, and arc layouts.
#' In branch-summary category mode, near-zero or high-outlier rate diagnostics
#' are drawn as special categories when rate-valued columns are plotted; these
#' special categories do not consume positions in the ordered palette used for
#' regular rate bins. When special categories are present, the category legend
#' spans the full plotted value range and marks the diagnostic cutoff separating
#' special and regular bins. When plotting non-rate columns such as `"sd"`,
#' diagnostic columns are preserved only as metadata with `rate_flag_source`
#' provenance; special rate categories are not drawn for those non-rate values.
#' The selected `value` column must contain at least one finite value;
#' uncertainty columns such as `"sd"` may be all `NA` for single-fit objects.
#'
#' @param x An object of class `"rateMap"` returned by [rateMap()] or
#'   [rateMapView()].
#' @param value Character column in the plotted summary table (`x$intervals`)
#'   to map to branch colors. The default `"value"` plots the central estimate
#'   chosen by [rateMap()]. When `uncertainty = TRUE`, useful alternatives
#'   include `"mean"`, `"median"`, `"sd"`, `"ci_width"`,
#'   `"highest_density_interval_width"`, and `"cv"`.
#' @param palette Optional palette override used for this plot. This can be an
#'   `hcl.colors()` palette name, a vector of colors, or a palette function.
#' @param reverse_palette Optional logical override for palette reversal used
#'   for this plot.
#' @param ncolors Optional number of colors to use when recoloring this plot
#'   with `color_mode = "continuous"`. If omitted, the stored `x$ncolors` value
#'   is used, falling back to `256` for older objects.
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
#' @param legend_title Optional legend title override for this plot.
#' @param legend Legend length. If `NULL` or `TRUE`, a layout-specific default is
#'   used. Set `legend = FALSE` to suppress the legend. Continuous legends are
#'   drawn with [phytools::add.color.bar()]; category legends are drawn by
#'   `rateMap` so bin labels can reflect rate categories and diagnostics.
#' @param fsize Numeric font-size vector. The first element is passed as the tip
#'   label `fsize` to [phytools::plotSimmap()] and, when `outline = TRUE`, to
#'   [phytools::plotTree()]. The second element is used by the legend.
#' @param tip_fsize Optional override for the first `fsize` element passed to
#'   the underlying tree plot.
#' @param legend_fsize Optional override for the legend font size.
#' @param ftype Font type. The first element is passed as `ftype` to
#'   [phytools::plotSimmap()] and, when `outline = TRUE`, to
#'   [phytools::plotTree()]. The second element is reserved for legends.
#' @param show_tip_labels Logical; if `FALSE`, `rateMap` sets the tree-plot
#'   `ftype` to `"off"` before calling the underlying `phytools` plotter.
#' @param outline Logical; if `TRUE`, draw a branch outline beneath the rate map.
#'   The outline pass is drawn with [phytools::plotTree()] before the colored
#'   [phytools::plotSimmap()] pass.
#' @param lwd Branch and legend line widths. The first element is passed to
#'   [phytools::plotSimmap()] and, with `+ 2`, to [phytools::plotTree()] for
#'   outlines. The second element is used for the legend.
#' @param type Plot type: `"phylogram"`, `"fan"`, or `"arc"`. This is passed to
#'   [phytools::plotSimmap()] for fan and arc layouts and to
#'   [phytools::plotTree()] for outline passes.
#' @param mar Plot margins passed to [phytools::plotSimmap()] and, when
#'   `outline = TRUE`, to [phytools::plotTree()].
#' @param direction Plotting direction for `type = "phylogram"`; passed to
#'   [phytools::plotSimmap()] and [phytools::plotTree()].
#' @param offset Tip-label offset passed to [phytools::plotSimmap()] and, when
#'   `outline = TRUE`, to [phytools::plotTree()].
#' @param xlim,ylim Optional plot limits passed to [phytools::plotSimmap()] and,
#'   when `outline = TRUE`, to [phytools::plotTree()].
#' @param hold Logical controlling `rateMap`'s device hold/flush guard.
#'   `rateMap` passes `hold = FALSE` to the internal `phytools` calls so the
#'   two-pass outline/color drawing is controlled in one place.
#' @param underscore Logical; if `FALSE`, underscores in tip labels may be shown
#'   as spaces by the underlying [phytools::plotSimmap()] and
#'   [phytools::plotTree()] calls.
#' @param arc_height Arc height passed through to [phytools::plotSimmap()] and,
#'   for outlines, [phytools::plotTree()] when `type = "arc"`.
#' @param legend_digits Optional number of digits for legend endpoint labels. If
#'   omitted, small-magnitude values use enough digits to avoid zero-valued
#'   legend endpoints.
#' @param ... Additional arguments are rejected. Include all display choices as
#'   named `plot()` arguments.
#'
#' @details
#' The `plot()` method handles `"phylogram"`, `"fan"`, and `"arc"` layouts.
#' `legend` controls the length of the color bar; set `legend = FALSE` to
#' suppress it. `legend_digits` controls numeric endpoint labels and defaults
#' to enough precision to avoid rounding small values to zero. Single fitted
#' objects should be converted explicitly with [rateMap()] before plotting, for
#' example `plot(rateMap(search_a), ...)`.
#'
#' **Relationship to `phytools` plotting arguments.** `plot.rateMap()` keeps the
#' tree-layout argument names close to `phytools`: `type`, `fsize`, `ftype`,
#' `lwd`, `mar`, `direction`, `offset`, `xlim`, `ylim`, `underscore`, and
#' `arc_height` are forwarded to [phytools::plotSimmap()] or
#' [phytools::plotTree()] as described above. `rateMap` fixes
#' `phytools::plotSimmap()` options such as `colors`, `pts`, `node.numbers`,
#' `add`, and `hold` internally, because those are determined by the computed
#' `"rateMap"` object and by the optional outline pass. `palette`,
#' `color_mode`, `n_categories`, `category_breaks`, `category_labels`,
#' `legend_title`, `legend_digits`, and the rate-flag display behavior are
#' `rateMap` controls, not `phytools` arguments.
#'
#' @return Invisibly returns the plotted `"rateMap"` object.
#'
#' @seealso [rateMap()], [phytools::densityMap()], [phytools::plotSimmap()],
#'   [phytools::plotTree()]
#'
#' @references
#' Revell, L. J. (2013). Two new graphical methods for mapping trait evolution
#' on phylogenies. *Methods in Ecology and Evolution*, 4, 754-759.
#'
#' Revell, L. J. (2024). phytools 2.0: an updated R ecosystem for phylogenetic
#' comparative methods (and other things). *PeerJ*, 12, e16505.
#' \doi{10.7717/peerj.16505}
#'
#' @examples
#' \dontrun{
#' rm_obj <- rateMap(fits, progress = FALSE)
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
#' # Use a continuous ramp instead of the default ordered rate categories:
#' plot(rm_obj, color_mode = "continuous")
#'
#' # Or keep category colors but change the binning:
#' plot(rm_obj, n_categories = 5, category_bin_method = "equal")
#' plot(rm_obj, category_breaks = c(-4, -2, 0, 2))
#' }
#'
#' @method plot rateMap
#' @export
plot.rateMap <- function(x,
                         value = "value",
                         palette = NULL,
                         reverse_palette = NULL,
                         color_mode = NULL,
                         n_categories = NULL,
                         category_bin_method = NULL,
                         category_breaks = NULL,
                         category_labels = NULL,
                         ncolors = NULL,
                         legend_title = NULL,
                         legend = NULL,
                         fsize = NULL,
                         tip_fsize = NULL,
                         legend_fsize = NULL,
                         ftype = NULL,
                         show_tip_labels = TRUE,
                         outline = FALSE,
                         lwd = 3,
                         type = c("phylogram", "fan", "arc"),
                         mar = rep(0.3, 4),
                         direction = "rightwards",
                         offset = NULL,
                         xlim = NULL,
                         ylim = NULL,
                         hold = TRUE,
                         underscore = FALSE,
                         arc_height = 2,
                         legend_digits = NULL,
                         ...) {
  .plot_rate_map(
    x,
    value = value,
    palette = palette,
    reverse_palette = reverse_palette,
    color_mode = color_mode,
    n_categories = n_categories,
    category_bin_method = category_bin_method,
    category_breaks = category_breaks,
    category_labels = category_labels,
    ncolors = ncolors,
    legend_title = legend_title,
    legend = legend,
    fsize = fsize,
    tip_fsize = tip_fsize,
    legend_fsize = legend_fsize,
    ftype = ftype,
    show_tip_labels = show_tip_labels,
    outline = outline,
    lwd = lwd,
    type = type,
    mar = mar,
    direction = direction,
    offset = offset,
    xlim = xlim,
    ylim = ylim,
    hold = hold,
    underscore = underscore,
    arc_height = arc_height,
    legend_digits = legend_digits,
    ...
  )
}

#' Print a `rateMap` Object
#'
#' Print a concise summary of a `"rateMap"` object, including the number of fits,
#' summary mode, target/check mode, weighting mode, uncertainty status, color
#' mode, plotted value, value range, rate-flag diagnostics when active, and
#' uncertainty ranges when available.
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

  if (!is.null(x$rate_diagnostics) &&
      isTRUE(x$rate_diagnostics$enabled)) {
    cat(
      "- rate flags: ",
      x$rate_diagnostics$n_near_zero,
      " near-zero, ",
      x$rate_diagnostics$n_high_outlier,
      " high-outlier (method: ",
      x$rate_diagnostics$method,
      ")\n",
      sep = ""
    )
    if (!is.na(x$rate_diagnostics$zero_floor)) {
      cat(
        "- near-zero floor: ",
        formatC(x$rate_diagnostics$zero_floor, digits = 4, format = "fg"),
        "\n",
        sep = ""
      )
    }
    if (is.finite(x$rate_diagnostics$full_fold_range) &&
        is.finite(x$rate_diagnostics$regular_fold_range)) {
      cat(
        "- fold range: ",
        trimws(formatC(x$rate_diagnostics$full_fold_range, digits = 4, format = "fg")),
        " full, ",
        trimws(formatC(x$rate_diagnostics$regular_fold_range, digits = 4, format = "fg")),
        " regular\n",
        sep = ""
      )
    }
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
