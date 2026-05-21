# tests/testthat/test-rate-map.R

.rate_map_skip_if_missing_deps <- function() {
  testthat::skip_if_not_installed("ape")
  testthat::skip_if_not_installed("phytools")
  testthat::skip_if_not_installed("future")
  testthat::skip_if_not_installed("future.apply")
  testthat::skip_if_not_installed("progressr")
}

.rate_map_fixture <- function() {
  set.seed(20260521)
  tr <- ape::rtree(8)
  base <- phytools::paintSubTree(
    tree = ape::as.phylo(tr),
    node = ape::Ntip(tr) + 1L,
    state = "0"
  )
  shifted_a <- phytools::paintSubTree(
    tree = base,
    node = ape::Ntip(base) + 2L,
    state = "1"
  )
  shifted_b <- phytools::paintSubTree(
    tree = base,
    node = ape::Ntip(base) + 3L,
    state = "1"
  )
  list(base = base, shifted_a = shifted_a, shifted_b = shifted_b)
}

.rate_map_fit <- function(tree, param = c("0" = 0.1, "1" = 0.5)) {
  list(
    variables = list(tree = tree),
    param = param
  )
}

.rate_map_two_tip <- function(edge_length = c(1, 1), maps = NULL) {
  tr <- ape::read.tree(text = "(a:1,b:1);")
  tr$edge.length <- edge_length

  if (is.null(maps)) {
    maps <- lapply(edge_length, function(len) c("0" = len))
  }

  tr$maps <- maps
  tr$mapped.edge <- .rateMap_make_mapped_edge(tr$edge, tr$maps)
  class(tr) <- c("simmap", "phylo")
  tr
}

test_that("rateMap summarizes small synthetic stochastic-map fits", {
  .rate_map_skip_if_missing_deps()
  fx <- .rate_map_fixture()

  fits <- list(
    .rate_map_fit(fx$shifted_a, c("0" = 0.10, "1" = 0.50)),
    .rate_map_fit(fx$shifted_b, c("0" = 0.20, "1" = 0.80))
  )

  out <- rateMap(
    fits,
    res = 8,
    plot = FALSE,
    progress = FALSE,
    ncolors = 16,
    palette = "Viridis"
  )

  testthat::expect_s3_class(out, "rateMap")
  testthat::expect_s3_class(out$tree, "simmap")
  testthat::expect_length(out$cols, 16)
  testthat::expect_length(out$lims, 2)
  testthat::expect_true(is.data.frame(out$intervals))
  testthat::expect_true(all(c(
    "edge",
    "parent",
    "child",
    "depth_start",
    "depth_end",
    "interval_length",
    "value",
    "color_bin"
  ) %in% names(out$intervals)))
  testthat::expect_equal(out$n_fits, 2L)
})

test_that("rateMap auto-detects bifrost_search-style objects", {
  .rate_map_skip_if_missing_deps()
  fx <- .rate_map_fixture()

  fit <- structure(
    list(
      tree_no_uncertainty_untransformed = fx$shifted_a,
      model_no_uncertainty = list(param = c("0" = 0.25, "1" = 0.75))
    ),
    class = c("bifrost_search", "list")
  )

  out <- rateMap(
    list(fit, fit),
    res = 4,
    plot = FALSE,
    progress = FALSE,
    ncolors = 8
  )

  testthat::expect_s3_class(out, "rateMap")
  testthat::expect_equal(out$n_fits, 2L)
  testthat::expect_equal(range(out$intervals$value), out$lims)
})

test_that("rateMap validates mapped state and parameter compatibility", {
  .rate_map_skip_if_missing_deps()
  fx <- .rate_map_fixture()

  testthat::expect_error(
    rateMap(
      list(list(param = c("0" = 0.1))),
      plot = FALSE,
      progress = FALSE
    ),
    "mapped simmap"
  )

  testthat::expect_error(
    rateMap(
      list(list(variables = list(tree = fx$base))),
      plot = FALSE,
      progress = FALSE
    ),
    "named numeric"
  )

  testthat::expect_error(
    rateMap(
      list(.rate_map_fit(fx$shifted_a, c("0" = 0.1))),
      plot = FALSE,
      progress = FALSE
    ),
    "missing rate parameters"
  )

  testthat::expect_error(
    rateMap(
      list(.rate_map_fit(fx$shifted_a, c("0" = 0.1, "1" = 0))),
      log = TRUE,
      plot = FALSE,
      progress = FALSE
    ),
    "non-positive"
  )

  testthat::expect_error(
    rateMap(
      list(.rate_map_fit(fx$shifted_a, c("0" = 0.1, "1" = Inf))),
      plot = FALSE,
      progress = FALSE
    ),
    "non-finite"
  )
})

test_that("rateMap errors on no-shift NA params by default and can omit them", {
  .rate_map_skip_if_missing_deps()
  fx <- .rate_map_fixture()

  invalid <- .rate_map_fit(fx$base, NA_real_)
  valid <- .rate_map_fit(fx$shifted_a, c("0" = 0.1, "1" = 0.5))

  testthat::expect_error(
    rateMap(list(invalid), plot = FALSE, progress = FALSE),
    "named numeric parameter"
  )
  testthat::expect_error(
    rateMap(list(invalid), plot = FALSE, progress = FALSE, na_action = "omit"),
    "No fits remain"
  )

  out <- rateMap(
    list(invalid, valid),
    res = 4,
    plot = FALSE,
    progress = FALSE,
    na_action = "omit"
  )

  testthat::expect_equal(out$n_fits, 1L)
  testthat::expect_equal(out$omitted, 1L)
})

test_that("rateMap supports all default extraction shapes", {
  .rate_map_skip_if_missing_deps()
  fx <- .rate_map_fixture()

  mvgls_stub <- structure(
    list(corrSt = list(phy = fx$base), param = c("0" = 0.2)),
    class = "mvgls"
  )
  wrapped_mvgls <- list(model = mvgls_stub)
  transformed_only <- structure(
    list(
      tree_no_uncertainty_transformed = fx$base,
      model_no_uncertainty = list(param = c("0" = 0.4))
    ),
    class = c("bifrost_search", "list")
  )

  out <- rateMap(
    list(mvgls_stub, wrapped_mvgls, transformed_only),
    summary = "branch",
    plot = FALSE,
    progress = FALSE
  )

  testthat::expect_equal(out$n_fits, 3L)
  testthat::expect_equal(out$weight_mode, "equal")
  testthat::expect_true(all(abs(out$intervals$value - mean(c(0.2, 0.2, 0.4))) < 1e-12))
})

test_that("rateMap validates top-level arguments and tree checks", {
  .rate_map_skip_if_missing_deps()
  fx <- .rate_map_fixture()
  fit <- .rate_map_fit(fx$base, c("0" = 1))

  testthat::expect_error(rateMap(NULL), "non-empty list")
  testthat::expect_error(rateMap(list(fit), res = 0), "res")
  testthat::expect_error(rateMap(list(fit), ncolors = 1), "ncolors")
  testthat::expect_error(rateMap(list(fit), workers = 0), "workers")
  testthat::expect_error(rateMap(list(fit), summary = "node"), "'arg' should be")
  testthat::expect_error(rateMap(list(fit), check = "branch"), "check")
  testthat::expect_error(
    rateMap(list(fit), target_tree = list(), plot = FALSE, progress = FALSE),
    "target_tree"
  )

  changed <- fx$base
  changed$edge.length[1L] <- changed$edge.length[1L] * 2
  testthat::expect_error(
    rateMap(
      list(fit, .rate_map_fit(changed, c("0" = 1))),
      plot = FALSE,
      progress = FALSE
    ),
    "do not match"
  )

  phylo_with_maps <- fx$base
  class(phylo_with_maps) <- "phylo"
  out <- rateMap(
    list(.rate_map_fit(phylo_with_maps, c("0" = 1))),
    summary = "branch",
    plot = FALSE,
    progress = FALSE
  )
  testthat::expect_s3_class(out, "rateMap")
})

test_that("rateMap projects same-topology branch-length variation onto a target tree", {
  .rate_map_skip_if_missing_deps()

  target <- .rate_map_two_tip(edge_length = c(1, 1))
  edge_a <- which(target$edge[, 2L] == match("a", target$tip.label))
  edge_b <- which(target$edge[, 2L] == match("b", target$tip.label))

  sampled <- target
  sampled$edge.length[edge_a] <- 2
  sampled$edge.length[edge_b] <- 0.5
  sampled$maps[[edge_a]] <- c("0" = 1, "1" = 1)
  sampled$maps[[edge_b]] <- c("1" = 0.5)
  sampled$mapped.edge <- .rateMap_make_mapped_edge(sampled$edge, sampled$maps)

  testthat::expect_error(
    rateMap(
      list(.rate_map_fit(target, c("0" = 1, "1" = 3)), .rate_map_fit(sampled, c("0" = 1, "1" = 3))),
      summary = "branch",
      plot = FALSE,
      progress = FALSE
    ),
    "do not match"
  )

  branch_out <- rateMap(
    list(.rate_map_fit(sampled, c("0" = 1, "1" = 3))),
    target_tree = ape::as.phylo(target),
    check = "topology",
    summary = "branch",
    plot = FALSE,
    progress = FALSE
  )

  testthat::expect_equal(branch_out$target, "user")
  testthat::expect_equal(branch_out$check, "topology")
  testthat::expect_equal(branch_out$tree$edge.length, target$edge.length)
  testthat::expect_equal(branch_out$intervals$value[branch_out$intervals$edge == edge_a], 2)
  testthat::expect_equal(branch_out$intervals$value[branch_out$intervals$edge == edge_b], 3)

  interval_out <- rateMap(
    list(.rate_map_fit(sampled, c("0" = 1, "1" = 3))),
    target_tree = ape::as.phylo(target),
    check = "topology",
    summary = "interval",
    res = 2,
    plot = FALSE,
    progress = FALSE
  )
  edge_a_rows <- interval_out$intervals[interval_out$intervals$edge == edge_a, ]
  testthat::expect_equal(edge_a_rows$interval_length, c(0.5, 0.5))
  testthat::expect_equal(edge_a_rows$value, c(1, 3))
  testthat::expect_equal(
    as.numeric(rowSums(interval_out$tree$mapped.edge)),
    target$edge.length,
    tolerance = 1e-10
  )

  permuted <- sampled
  ord <- rev(seq_len(nrow(permuted$edge)))
  permuted$edge <- permuted$edge[ord, , drop = FALSE]
  permuted$edge.length <- permuted$edge.length[ord]
  permuted$maps <- permuted$maps[ord]
  permuted$mapped.edge <- .rateMap_make_mapped_edge(permuted$edge, permuted$maps)

  permuted_out <- rateMap(
    list(.rate_map_fit(permuted, c("0" = 1, "1" = 3))),
    target_tree = ape::as.phylo(target),
    check = FALSE,
    summary = "branch",
    plot = FALSE,
    progress = FALSE
  )
  testthat::expect_equal(permuted_out$check, "none")
  testthat::expect_equal(permuted_out$intervals$value[permuted_out$intervals$edge == edge_a], 2)

  mcc_out <- rateMap(
    list(.rate_map_fit(sampled, c("0" = 1, "1" = 3)), .rate_map_fit(permuted, c("0" = 1, "1" = 3))),
    target = "mcc",
    check = "topology",
    summary = "branch",
    plot = FALSE,
    progress = FALSE
  )
  testthat::expect_equal(mcc_out$target, "mcc")
  testthat::expect_equal(mcc_out$tree$edge.length, sampled$edge.length)

  bad_target <- ape::read.tree(text = "(a:1,c:1);")
  testthat::expect_error(
    rateMap(
      list(.rate_map_fit(sampled, c("0" = 1, "1" = 3))),
      target_tree = bad_target,
      check = "topology",
      summary = "branch",
      plot = FALSE,
      progress = FALSE
    ),
    "topology"
  )
  testthat::expect_error(
    rateMap(
      list(.rate_map_fit(sampled, c("0" = 1, "1" = 3))),
      target_tree = bad_target,
      check = FALSE,
      summary = "branch",
      plot = FALSE,
      progress = FALSE
    ),
    "missing one or more target-tree branches"
  )
})

test_that("rateMap mapped.edge row sums match edge lengths", {
  .rate_map_skip_if_missing_deps()
  fx <- .rate_map_fixture()
  out <- rateMap(
    list(.rate_map_fit(fx$shifted_a), .rate_map_fit(fx$shifted_a)),
    res = 10,
    plot = FALSE,
    progress = FALSE
  )

  testthat::expect_equal(
    as.numeric(rowSums(out$tree$mapped.edge)),
    out$tree$edge.length,
    tolerance = 1e-10
  )
})

test_that("rateMap computes whole-branch summaries without splitting branches", {
  .rate_map_skip_if_missing_deps()
  fx <- .rate_map_fixture()

  out <- rateMap(
    list(
      .rate_map_fit(fx$base, c("0" = 0.1)),
      .rate_map_fit(fx$base, c("0" = 0.5))
    ),
    summary = "branch",
    weights = c(1, 3),
    plot = FALSE,
    progress = FALSE
  )

  testthat::expect_equal(out$summary, "branch")
  testthat::expect_equal(nrow(out$intervals), nrow(fx$base$edge))
  testthat::expect_equal(out$intervals$interval_length, out$tree$edge.length)
  testthat::expect_equal(out$weights, c(0.25, 0.75))
  testthat::expect_equal(out$weight_mode, "custom")
  testthat::expect_true(all(abs(out$intervals$value - 0.4) < 1e-12))
  testthat::expect_true(all(vapply(out$tree$maps, length, integer(1)) == 1L))
})

test_that("rateMap handles valid log transforms and multi-segment branch maps", {
  .rate_map_skip_if_missing_deps()
  fx <- .rate_map_fixture()

  segmented <- fx$base
  edge_len <- segmented$edge.length[1L]
  segmented$maps[[1L]] <- c("0" = edge_len / 2, "1" = edge_len / 2)
  segmented$mapped.edge <- .rateMap_make_mapped_edge(segmented$edge, segmented$maps)

  out <- rateMap(
    list(.rate_map_fit(segmented, c("0" = exp(1), "1" = exp(3)))),
    summary = "branch",
    log = TRUE,
    plot = FALSE,
    progress = FALSE
  )

  testthat::expect_equal(out$summary, "branch")
  testthat::expect_true(any(abs(out$intervals$value - 2) < 1e-12))
})

test_that("rateMap handles zero-length mapped branches", {
  .rate_map_skip_if_missing_deps()
  fx <- .rate_map_fixture()

  zero_branch <- fx$base
  zero_branch$edge.length[1L] <- 0
  zero_branch$maps[[1L]] <- c("0" = 0)
  zero_branch$mapped.edge <- .rateMap_make_mapped_edge(zero_branch$edge, zero_branch$maps)

  out <- rateMap(
    list(.rate_map_fit(zero_branch, c("0" = 1))),
    summary = "branch",
    plot = FALSE,
    progress = FALSE
  )

  testthat::expect_equal(out$intervals$interval_length[1L], 0)
  testthat::expect_equal(out$intervals$value[1L], 1)
})

test_that("rateMap supports IC weights for comparable bifrost_search objects", {
  .rate_map_skip_if_missing_deps()
  fx <- .rate_map_fixture()

  make_search <- function(rate, ic, family = "GIC") {
    structure(
      list(
        tree_no_uncertainty_untransformed = fx$base,
        model_no_uncertainty = list(param = c("0" = rate)),
        optimal_ic = ic,
        IC_used = family
      ),
      class = c("bifrost_search", "list")
    )
  }

  fit_a <- make_search(10, 100)
  fit_b <- make_search(20, 102)
  expected_weights <- exp(-0.5 * c(0, 2))
  expected_weights <- expected_weights / sum(expected_weights)

  out <- rateMap(
    list(fit_a, fit_b),
    summary = "branch",
    weights = "ic",
    plot = FALSE,
    progress = FALSE
  )

  testthat::expect_equal(out$weight_mode, "ic")
  testthat::expect_equal(out$weights, expected_weights, tolerance = 1e-12)
  testthat::expect_equal(out$weight_table$ic, c(100, 102))
  testthat::expect_equal(out$weight_table$IC_used, c("GIC", "GIC"))
  testthat::expect_true(all(abs(out$intervals$value - sum(expected_weights * c(10, 20))) < 1e-12))

  testthat::expect_error(
    rateMap(
      list(fit_a, make_search(20, 102, "BIC")),
      weights = "ic",
      plot = FALSE,
      progress = FALSE
    ),
    "same non-missing 'IC_used'"
  )

  no_family <- fit_a
  no_family$IC_used <- NULL
  testthat::expect_error(
    rateMap(
      list(fit_a, no_family),
      weights = "ic",
      plot = FALSE,
      progress = FALSE
    ),
    "same non-missing 'IC_used'"
  )

  missing_ic <- fit_a
  missing_ic$optimal_ic <- NA_real_
  testthat::expect_error(
    rateMap(
      list(fit_a, missing_ic),
      weights = "ic",
      plot = FALSE,
      progress = FALSE
    ),
    "finite 'optimal_ic'"
  )

  testthat::expect_error(
    rateMap(
      list(.rate_map_fit(fx$base, c("0" = 1))),
      weights = "ic",
      plot = FALSE,
      progress = FALSE
    ),
    "finite 'optimal_ic'"
  )
})

test_that("rateMap validates custom fit weights", {
  .rate_map_skip_if_missing_deps()
  fx <- .rate_map_fixture()
  fits <- list(
    .rate_map_fit(fx$base, c("0" = 1)),
    .rate_map_fit(fx$base, c("0" = 3))
  )

  out <- rateMap(
    fits,
    summary = "branch",
    weights = "ic",
    fit_weights = c(2, 2),
    plot = FALSE,
    progress = FALSE
  )
  testthat::expect_equal(out$weight_mode, "custom")
  testthat::expect_equal(out$weights, c(0.5, 0.5))
  testthat::expect_true(all(abs(out$intervals$value - 2) < 1e-12))

  testthat::expect_error(
    rateMap(fits, weights = c(1, -1), plot = FALSE, progress = FALSE),
    "non-negative"
  )
  testthat::expect_error(
    rateMap(fits, weights = c(0, 0), plot = FALSE, progress = FALSE),
    "positive"
  )
  testthat::expect_error(
    rateMap(fits, weights = c(1, 2, 3), plot = FALSE, progress = FALSE),
    "length"
  )
  testthat::expect_error(
    rateMap(fits, weights = c(1, NA), plot = FALSE, progress = FALSE),
    "finite numeric"
  )

  invalid <- .rate_map_fit(fx$base, NA_real_)
  out_omit <- rateMap(
    list(invalid, fits[[1]], fits[[2]]),
    summary = "branch",
    weights = c(100, 1, 3),
    na_action = "omit",
    plot = FALSE,
    progress = FALSE
  )
  testthat::expect_equal(out_omit$weights, c(0.25, 0.75))

  out_omit_used_weights <- rateMap(
    list(invalid, fits[[1]], fits[[2]]),
    summary = "branch",
    weights = c(1, 3),
    na_action = "omit",
    plot = FALSE,
    progress = FALSE
  )
  testthat::expect_equal(out_omit_used_weights$weights, c(0.25, 0.75))
})

test_that("rateMap palette resolution handles names, vectors, and reversal", {
  .rate_map_skip_if_missing_deps()

  vir <- .rateMap_colors(5, "Viridis")
  vir_rev <- .rateMap_colors(5, "Viridis", reverse_palette = TRUE)
  testthat::expect_equal(vir_rev, rev(vir))

  ramp <- .rateMap_colors(6, c("black", "white"))
  ramp_rev <- .rateMap_colors(6, c("black", "white"), reverse_palette = TRUE)
  testthat::expect_equal(ramp_rev, rev(ramp))

  fun_cols <- .rateMap_colors(4, function(n) grDevices::hcl.colors(n, "Inferno"))
  testthat::expect_length(fun_cols, 4)
  testthat::expect_type(fun_cols, "character")

  resized <- .rateMap_colors(4, function(n) c("black", "white"))
  testthat::expect_length(resized, 4)

  testthat::expect_error(.rateMap_colors(3, 1:3), "palette")
  testthat::expect_error(.rateMap_colors(3, function(n) "black"), "at least two")

  testthat::expect_equal(.rateMap_order_state_names(c("b", "a")), c("a", "b"))
  testthat::expect_true(.rateMap_show_legend(0.1))
  testthat::expect_false(.rateMap_show_legend(FALSE))
  testthat::expect_equal(
    .rateMap_map_value(c("0" = 1), c("0" = 7), start = 2, end = 3),
    7
  )
})

test_that("rateMap can plot directly during computation and use progress/workers", {
  .rate_map_skip_if_missing_deps()
  fx <- .rate_map_fixture()

  path <- tempfile(fileext = ".pdf")
  grDevices::pdf(path)
  on.exit(grDevices::dev.off(), add = TRUE)

  testthat::expect_error(
    rateMap(
      list(.rate_map_fit(fx$base, c("0" = 1))),
      summary = "branch",
      plot = TRUE,
      legend = FALSE,
      progress = TRUE,
      workers = 1,
      future_strategy = "multisession"
    ),
    NA
  )
})

test_that("plotRateMap draws rate maps on a graphics device", {
  .rate_map_skip_if_missing_deps()
  fx <- .rate_map_fixture()
  out <- rateMap(
    list(.rate_map_fit(fx$shifted_a), .rate_map_fit(fx$shifted_b)),
    res = 4,
    plot = FALSE,
    progress = FALSE
  )

  path <- tempfile(fileext = ".pdf")
  grDevices::pdf(path)
  on.exit(grDevices::dev.off(), add = TRUE)

  testthat::expect_error(
    plotRateMap(out, type = "arc", show_tip_labels = FALSE, legend = FALSE, hold = FALSE),
    NA
  )

  testthat::expect_error(plotRateMap(list()), "rateMap")
})

test_that("plotRateMap covers phylogram, fan, legend, outline, and S3 paths", {
  .rate_map_skip_if_missing_deps()
  fx <- .rate_map_fixture()
  out <- rateMap(
    list(.rate_map_fit(fx$shifted_a), .rate_map_fit(fx$shifted_b)),
    res = 4,
    plot = FALSE,
    progress = FALSE
  )

  path <- tempfile(fileext = ".pdf")
  grDevices::pdf(path)
  on.exit(grDevices::dev.off(), add = TRUE)

  testthat::expect_error(
    plotRateMap(out, type = "arc", legend = TRUE, hold = FALSE),
    NA
  )

  testthat::expect_error(
    plotRateMap(out, type = "fan", legend = TRUE, hold = FALSE),
    NA
  )

  testthat::expect_message(
    plotRateMap(
      out,
      type = "phylogram",
      outline = TRUE,
      legend = 1e9,
      fsize = 0.5,
      ftype = "i",
      lwd = c(1, 2, 3),
      hold = TRUE
    ),
    "legend scale"
  )

  testthat::expect_error(
    plotRateMap(
      out,
      type = "phylogram",
      outline = TRUE,
      offset = 0.1,
      direction = "leftwards",
      xlim = c(1, 0),
      legend = 0.1,
      tip_fsize = 0.4,
      legend_fsize = 0.6,
      hold = FALSE
    ),
    NA
  )

  testthat::expect_error(
    plotRateMap(
      out,
      type = "fan",
      outline = TRUE,
      legend = 0.2,
      fsize = c(0.4, 0.5),
      hold = FALSE
    ),
    NA
  )

  testthat::expect_error(
    plot(
      out,
      type = "arc",
      outline = TRUE,
      legend = 0.2,
      show_tip_labels = FALSE,
      hold = FALSE
    ),
    NA
  )
})
