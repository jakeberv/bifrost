rate_map_artifact <- function(file) {
  file.path("rate-map-jaw-shape", file)
}

rate_map_print_file <- function(file) {
  cat(readLines(rate_map_artifact(file)), sep = "\n")
}

rate_map_cache_path <- function() {
  cache_rel <- file.path(
    "local-cache",
    "rate-map-jaw-shape",
    "jaw_threshold_runs_full.rds"
  )
  candidates <- c(
    file.path("..", cache_rel),
    cache_rel,
    "jaw_threshold_runs.rds"
  )
  existing <- candidates[file.exists(candidates)]
  if (length(existing) > 0L) {
    return(existing[1L])
  }
  if (file.exists(file.path("..", "DESCRIPTION"))) {
    return(file.path("..", cache_rel))
  }
  cache_rel
}

rate_map_caption <- function(lead, ...) {
  paste(c(paste0("<strong>", lead, "</strong>"), ...), collapse = " ")
}

rate_map_split_clade_key <- function(clade_key) {
  if (!is.character(clade_key) || length(clade_key) != 1L || is.na(clade_key)) {
    return(character())
  }

  taxa <- unlist(
    strsplit(clade_key, "\r\n|\n|\r|\\s*\\|\\s*", perl = TRUE),
    use.names = FALSE
  )
  taxa <- trimws(taxa)
  taxa[nzchar(taxa)]
}

rate_map_format_clade_key <- function(clade_key) {
  taxa <- rate_map_split_clade_key(clade_key)
  if (length(taxa) == 0L) {
    return(NA_character_)
  }

  paste(taxa, collapse = " | ")
}

rate_map_format_taxon_examples <- function(clade_key, max_examples = 2) {
  taxa <- rate_map_split_clade_key(clade_key)
  if (length(taxa) == 0L) {
    return(NA_character_)
  }

  taxa <- gsub("_", " ", taxa, fixed = TRUE)
  shown <- utils::head(taxa, max_examples)
  if (length(taxa) > length(shown)) {
    shown <- c(shown, "...")
  }

  paste(shown, collapse = "; ")
}

rate_map_apply_labels <- function(tab, labels = NULL) {
  if (is.null(labels)) {
    return(tab)
  }

  labels <- labels[names(labels) %in% names(tab)]
  if (length(labels) == 0L) {
    return(tab)
  }

  names(tab)[match(names(labels), names(tab))] <- unname(labels)
  tab
}

rate_map_table_data <- function(file, digits = 4, keep = NULL, drop = NULL,
                                n = NULL, arrange = NULL,
                                decreasing = FALSE, labels = NULL) {
  tab <- utils::read.csv(rate_map_artifact(file), check.names = FALSE)
  if ("clade_key" %in% names(tab)) {
    tab$descendant_examples <- vapply(
      tab$clade_key,
      rate_map_format_taxon_examples,
      character(1)
    )
  }
  if (!is.null(arrange) && arrange %in% names(tab)) {
    tab <- tab[order(tab[[arrange]], decreasing = decreasing), , drop = FALSE]
  }
  if (!is.null(keep)) {
    missing <- setdiff(keep, names(tab))
    if (length(missing) > 0L) {
      stop(
        "Static table artifact '", file, "' is missing column(s): ",
        paste(missing, collapse = ", ")
      )
    }
    tab <- tab[, intersect(keep, names(tab)), drop = FALSE]
  }
  if (!is.null(drop)) {
    tab <- tab[, setdiff(names(tab), drop), drop = FALSE]
  }
  if (!is.null(n)) {
    tab <- utils::head(tab, n)
  }
  numeric_cols <- vapply(tab, is.numeric, logical(1))
  tab[numeric_cols] <- lapply(tab[numeric_cols], function(x) {
    if (all(is.na(x) | x == round(x))) {
      return(x)
    }
    vapply(x, function(z) {
      if (is.na(z)) {
        return(NA_character_)
      }
      if (abs(z) >= 1000) {
        return(formatC(z, digits = 1, format = "f"))
      }
      if (z == 0 || abs(z) >= 1e-3) {
        return(as.character(signif(z, digits)))
      }
      formatC(z, digits = 3, format = "e")
    }, character(1))
  })
  rate_map_apply_labels(tab, labels)
}

rate_map_table <- function(file, digits = 4, keep = NULL, drop = NULL, n = NULL,
                           arrange = NULL, decreasing = FALSE,
                           caption = NULL, full_width_caption = FALSE,
                           labels = NULL) {
  tab <- rate_map_table_data(
    file,
    digits = digits,
    keep = keep,
    drop = drop,
    n = n,
    arrange = arrange,
    decreasing = decreasing,
    labels = labels
  )
  table_caption <- if (isTRUE(full_width_caption) && knitr::is_html_output()) {
    NULL
  } else {
    caption
  }
  table <- knitr::kable(tab, caption = table_caption, row.names = FALSE)
  if (knitr::is_html_output()) {
    table <- c(
      '<div class="rate-map-table-spacer" aria-hidden="true"></div>',
      '<div class="rate-map-table-scroller">',
      table,
      '</div>'
    )
    if (isTRUE(full_width_caption) && !is.null(caption)) {
      table <- c(
        table,
        paste0('<div class="rate-map-table-caption-full">', caption, '</div>')
      )
    }
  }
  knitr::asis_output(paste(c("", "", table), collapse = "\n"))
}

rate_map_table_html <- function(file, digits = 4, keep = NULL, drop = NULL,
                                n = NULL, arrange = NULL,
                                decreasing = FALSE,
                                caption = NULL, labels = NULL) {
  tab <- rate_map_table_data(
    file,
    digits = digits,
    keep = keep,
    drop = drop,
    n = n,
    arrange = arrange,
    decreasing = decreasing,
    labels = labels
  )
  table <- knitr::kable(
    tab,
    format = "html",
    caption = caption,
    escape = FALSE,
    row.names = FALSE
  )
  table <- c(
    '<div class="rate-map-table-spacer" aria-hidden="true"></div>',
    '<div class="rate-map-table-scroller">',
    table,
    '</div>'
  )
  paste(c("", "", table), collapse = "\n")
}

rate_map_legend_limits <- function(x) {
  categories <- x$rate_categories
  if (!is.null(categories) && nrow(categories) > 0L) {
    flags <- x$rate_flags
    if (!is.null(flags) &&
        any(categories$rate_category %in% c(flags$zero_label, flags$high_label))) {
      return(x$lims)
    }
    return(range(c(categories$lower, categories$upper), finite = TRUE))
  }
  x$lims
}

rate_map_format_legend <- function(x) {
  formatC(x, digits = 3, format = "fg")
}

rate_map_rate_marker_to_value <- function(rate, log) {
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

rate_map_legend_markers <- function(x) {
  categories <- x$rate_categories
  flags <- x$rate_flags
  diagnostics <- x$rate_diagnostics
  intervals <- x$intervals
  if (is.null(categories) || is.null(flags) || is.null(diagnostics)) {
    return(numeric())
  }

  markers <- numeric()
  if (flags$zero_label %in% categories$rate_category) {
    marker <- rate_map_rate_marker_to_value(diagnostics$zero_floor, x$log)
    if (!is.finite(marker)) {
      marker <- rate_map_rate_marker_to_value(diagnostics$zero_cluster_cutoff_rate, x$log)
    }
    if (!is.finite(marker) && "is_near_zero" %in% names(intervals)) {
      vals <- intervals$value[intervals$is_near_zero %in% TRUE]
      vals <- vals[is.finite(vals)]
      if (length(vals) > 0L) {
        marker <- max(vals)
      }
    }
    if (is.finite(marker)) {
      markers[flags$zero_label] <- marker
    }
  }
  if (flags$high_label %in% categories$rate_category) {
    marker <- rate_map_rate_marker_to_value(diagnostics$high_cluster_cutoff_rate, x$log)
    if (!is.finite(marker) && "is_high_outlier" %in% names(intervals)) {
      vals <- intervals$value[intervals$is_high_outlier %in% TRUE]
      vals <- vals[is.finite(vals)]
      if (length(vals) > 0L) {
        marker <- min(vals)
      }
    }
    if (is.finite(marker)) {
      markers[flags$high_label] <- marker
    }
  }

  markers
}

rate_map_add_color_bar <- function(x, legend, y_frac = 0.19, fsize = 0.8,
                                   lwd = 3) {
  usr <- graphics::par("usr")
  yr <- diff(usr[3:4])
  x0 <- mean(usr[1:2]) - 0.5 * legend
  x1 <- x0 + legend
  y0 <- usr[3] + y_frac * yr

  cols <- unname(x$cols)
  if (length(cols) == 1L) {
    cols <- rep(cols, 2L)
  }

  bar_height <- 0.012 * yr * max(lwd / 3, 0.75)
  label_gap <- 0.018 * yr * max(fsize, 0.75)
  title_gap <- 0.06 * yr * max(fsize, 0.75)
  lims <- rate_map_legend_limits(x)
  markers <- rate_map_legend_markers(x)

  op <- graphics::par(xpd = NA)
  on.exit(graphics::par(op), add = TRUE)

  if (length(markers) > 0L && is.finite(diff(x$lims)) && diff(x$lims) > 0) {
    flags <- x$rate_flags
    categories <- x$rate_categories
    full_lims <- x$lims
    markers <- pmin(pmax(markers, full_lims[1L]), full_lims[2L])
    value_to_x <- function(value) {
      x0 + (value - full_lims[1L]) / diff(full_lims) * legend
    }
    rect_by_value <- function(left, right, col) {
      xs <- sort(c(value_to_x(left), value_to_x(right)))
      graphics::rect(xs[1L], y0, xs[2L], y0 + bar_height, col = col, border = NA)
    }

    regular <- !categories$rate_category %in% c(flags$zero_label, flags$high_label)
    regular_cols <- categories$color[regular]
    zero_marker <- markers[flags$zero_label]
    high_marker <- markers[flags$high_label]
    regular_min <- if (is.finite(zero_marker)) zero_marker else full_lims[1L]
    regular_max <- if (is.finite(high_marker)) high_marker else full_lims[2L]

    if (is.finite(zero_marker)) {
      rect_by_value(full_lims[1L], zero_marker, flags$zero_color)
    }
    if (length(regular_cols) > 0L && regular_max > regular_min) {
      xs <- seq(
        value_to_x(regular_min),
        value_to_x(regular_max),
        length.out = length(regular_cols) + 1L
      )
      graphics::rect(
        xs[-length(xs)],
        y0,
        xs[-1L],
        y0 + bar_height,
        col = unname(regular_cols),
        border = NA
      )
    }
    if (is.finite(high_marker)) {
      rect_by_value(high_marker, full_lims[2L], flags$high_color)
    }

    label_values <- unique(c(full_lims[1L], unname(markers), full_lims[2L]))
    tick_x <- value_to_x(label_values)
    graphics::segments(tick_x, y0, tick_x, y0 + bar_height, col = "grey20", lwd = 0.5)
    graphics::text(
      tick_x,
      rep(y0 + bar_height + label_gap, length(tick_x)),
      labels = rate_map_format_legend(label_values),
      cex = fsize,
      adj = c(0.5, 0.5)
    )
  } else {
    xs <- seq(x0, x1, length.out = length(cols) + 1L)
    graphics::rect(
      xs[-length(xs)],
      y0,
      xs[-1L],
      y0 + bar_height,
      col = cols,
      border = NA
    )
    graphics::text(
      c(x0, x1),
      rep(y0 + bar_height + label_gap, 2L),
      labels = rate_map_format_legend(lims),
      cex = fsize,
      adj = c(0.5, 0.5)
    )
  }
  graphics::text(
    mean(c(x0, x1)),
    y0 + bar_height + title_gap,
    labels = x$title,
    cex = fsize,
    font = 2
  )
}

rate_map_trim_bottom <- function(path, pixels) {
  if (!requireNamespace("png", quietly = TRUE)) {
    stop("The png package is needed to rebuild static vignette figures.")
  }
  if (!is.numeric(pixels) || length(pixels) != 1L || pixels <= 0) {
    return(invisible(path))
  }
  img <- png::readPNG(path)
  keep <- seq_len(max(1L, dim(img)[1L] - as.integer(pixels)))
  png::writePNG(img[keep, , , drop = FALSE], path)
  invisible(path)
}

rate_map_save_arc_figure <- function(file, x, ..., width = 8, height = 5.2,
                                     dpi = 225) {
  path <- rate_map_artifact(file)
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  grDevices::png(path, width = width, height = height, units = "in", res = dpi)
  device_open <- TRUE
  on.exit(if (isTRUE(device_open)) grDevices::dev.off(), add = TRUE)
  legend <- 0.7 * max(phytools::nodeHeights(x$tree))
  plotted <- plot(
    x,
    type = "arc",
    show_tip_labels = FALSE,
    legend = FALSE,
    legend_fsize = 0.8,
    arc_height = 0.5,
    mar = rep(1.2, 4),
    ylim = c(-25, 80),
    ...
  )
  rate_map_add_color_bar(plotted, legend = legend)
  grDevices::dev.off()
  device_open <- FALSE
  rate_map_trim_bottom(path, pixels = round(0.7 * dpi))
  invisible(path)
}

rate_map_include_style <- function(file = "rate-map-style.css") {
  if (knitr::is_html_output()) {
    css <- readLines(rate_map_artifact(file), warn = FALSE)
    cat("<style>\n", paste(css, collapse = "\n"), "\n</style>\n", sep = "")
  }
}
