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

#' Summarize Branchwise Rate Variation Across Multiple Runs
#'
#' Computes a branch-interval summary of fitted regime-specific rates across a
#' list of stochastic-map-aware model fits or completed `bifrost_search`
#' results. Each branch is sliced on a global depth grid, the rate associated
#' with each mapped state is looked up for each run, and interval values are
#' averaged across runs.
#'
#' @param fits A non-empty list of completed runs or fitted model objects.
#'   Supported shapes are `bifrost_search` objects, `mvgls` objects,
#'   `list(model = <mvgls>)`, and scratch-style lists with
#'   `variables$tree` plus `param`.
#' @param res Integer resolution of the global depth grid used to subdivide
#'   branches.
#' @param fsize Optional plotting font sizes passed to [plotRateMap()] when
#'   `plot = TRUE`.
#' @param ftype Optional plotting font type(s) passed to [plotRateMap()] when
#'   `plot = TRUE`.
#' @param lwd Line width passed through to [plotRateMap()] when `plot = TRUE`.
#' @param check Logical; if `TRUE`, verify that extracted trees match in
#'   topology and branch lengths before aggregation.
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
#'   \item{`values`}{List of branch-interval mean values before color binning.}
#'   \item{`intervals`}{Data frame with one row per plotted branch interval.}
#'   \item{`title`}{Legend title used for plotting.}
#'   \item{`n_fits`}{Number of fits used after validation or omission.}
#'   \item{`omitted`}{Integer indices of omitted fits when `na_action = "omit"`.}
#' }
#'
#' @seealso [plotRateMap()], [phytools::plotSimmap()]
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
  ...
) {
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

  res <- as.integer(res)
  ncolors <- as.integer(ncolors)
  if (!is.null(workers)) {
    workers <- as.integer(workers)
  }
  future_strategy <- match.arg(future_strategy)
  na_action <- match.arg(na_action)
  type <- match.arg(type, c("phylogram", "fan", "arc"))

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

  if (isTRUE(log)) {
    params <- lapply(params, base::log)
  }

  if (isTRUE(check)) {
    ref_tree <- ape::as.phylo(trees[[1L]])
    same_tree <- vapply(
      trees,
      function(tr) isTRUE(ape::all.equal.phylo(ref_tree, ape::as.phylo(tr))),
      logical(1)
    )
    if (!all(same_tree)) {
      stop("Some trees do not match in topology or branch lengths.")
    }
  }

  if (all(vapply(trees, inherits, logical(1), what = "simmap"))) {
    class(trees) <- c("multiSimmap", "multiPhylo")
  } else {
    class(trees) <- "multiPhylo"
  }

  heights <- vapply(trees, function(x) max(phytools::nodeHeights(x)), numeric(1))
  steps <- (0:res / res) * max(heights)
  trees <- phytools::rescaleSimmap(trees, totalDepth = max(heights))
  tree <- trees[[1L]]
  H <- phytools::nodeHeights(tree)
  tol <- 1e-10

  edge_indices <- seq_len(nrow(tree$edge))

  compute_edge <- function(i) {
    mids <- intersect(which(steps > H[i, 1L]), which(steps < H[i, 2L]))
    yy <- cbind(
      c(H[i, 1L], steps[mids]),
      c(steps[mids], H[i, 2L])
    ) - H[i, 1L]
    lengths_i <- yy[, 2L] - yy[, 1L]

    values_i <- numeric(nrow(yy))

    for (j in seq_along(trees)) {
      map_i <- trees[[j]]$maps[[i]]
      state_names <- names(map_i)
      fit_rates <- params[[j]]

      xx <- matrix(
        0,
        nrow = length(map_i),
        ncol = 2L,
        dimnames = list(state_names, c("start", "end"))
      )
      xx[1L, 2L] <- map_i[1L]

      if (length(map_i) > 1L) {
        for (k in 2:length(map_i)) {
          xx[k, 1L] <- xx[k - 1L, 2L]
          xx[k, 2L] <- xx[k, 1L] + map_i[k]
        }
      }

      for (k in seq_len(nrow(yy))) {
        lower <- which(xx[, 1L] <= yy[k, 1L])
        lower <- lower[length(lower)]
        upper <- which(xx[, 2L] >= (yy[k, 2L] - tol))[1L]
        names(lower) <- names(upper) <- NULL

        if (!all(xx == 0)) {
          value_k <- 0
          for (m in lower:upper) {
            overlap <- min(xx[m, 2L], yy[k, 2L]) - max(xx[m, 1L], yy[k, 1L])
            frac <- overlap / (yy[k, 2L] - yy[k, 1L])
            value_k <- value_k + frac * unname(fit_rates[rownames(xx)[m]])
          }
        } else {
          value_k <- unname(fit_rates[rownames(xx)[1L]])
        }

        values_i[k] <- values_i[k] + value_k / length(trees)
      }
    }

    if (isTRUE(progress)) {
      p(sprintf("edge %d", i))
    }

    list(
      values = values_i,
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

  lims <- range(unlist(interval_values, use.names = FALSE), finite = TRUE)
  if (!all(is.finite(lims))) {
    stop("Could not compute finite rate limits from the supplied fits.")
  }

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
      data.frame(
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
    })
  )
  rownames(intervals) <- NULL

  out <- list(
    tree = tree,
    cols = cols,
    lims = lims,
    breaks = breaks,
    values = interval_values,
    intervals = intervals,
    title = legend_title,
    n_fits = length(trees),
    omitted = omitted
  )
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
#'
#' @param x An object of class `"rateMap"` returned by [rateMap()].
#' @param ... Additional plotting controls. Supported options include
#'   `legend`, `fsize`, `tip_fsize`, `legend_fsize`, `ftype`,
#'   `show_tip_labels`, `outline`, `lwd`, `type`, `mar`, `direction`,
#'   `offset`, `xlim`, `ylim`, `hold`, `underscore`, and `arc_height`.
#'
#' @return Invisibly returns `x`.
#'
#' @seealso [rateMap()]
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
#' }
#'
#' @export
plotRateMap <- function(x, ...) {
  if (!inherits(x, "rateMap")) {
    stop("'x' must be an object of class 'rateMap'.")
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
  }

  invisible(x)
}

#' @rdname plotRateMap
#' @method plot rateMap
#' @export
plot.rateMap <- function(x, ...) {
  plotRateMap(x, ...)
}
