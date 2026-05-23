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
  testthat::expect_equal(out$summary, "branch")
  testthat::expect_equal(out$color_mode, "category")
  testthat::expect_true(out$log)
  testthat::expect_s3_class(out$rate_categories, "data.frame")
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
    "color_bin",
    "rate_category"
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
  testthat::expect_true(out$log)
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
    log = FALSE,
    plot = FALSE,
    progress = FALSE
  )

  testthat::expect_equal(out$n_fits, 3L)
  testthat::expect_equal(out$weight_mode, "equal")
  testthat::expect_true(all(abs(out$intervals$value - mean(c(0.2, 0.2, 0.4))) < 1e-12))
})

test_that("rateMap and plotRateMap accept a single fitted object", {
  .rate_map_skip_if_missing_deps()
  fx <- .rate_map_fixture()

  scratch_fit <- .rate_map_fit(fx$base, c("0" = 0.2))
  scratch_out <- rateMap(
    scratch_fit,
    plot = FALSE,
    progress = FALSE
  )
  testthat::expect_s3_class(scratch_out, "rateMap")
  testthat::expect_equal(scratch_out$n_fits, 1L)
  testthat::expect_equal(scratch_out$summary, "branch")
  testthat::expect_true(scratch_out$log)

  bifrost_fit <- structure(
    list(
      tree_no_uncertainty_untransformed = fx$base,
      model_no_uncertainty = list(param = c("0" = 0.3))
    ),
    class = c("bifrost_search", "list")
  )
  bifrost_out <- rateMap(
    bifrost_fit,
    plot = FALSE,
    progress = FALSE
  )
  testthat::expect_equal(bifrost_out$n_fits, 1L)
  testthat::expect_true(bifrost_out$log)
  testthat::expect_equal(bifrost_out$intervals$value, rep(log(0.3), nrow(fx$base$edge)))

  path <- tempfile(fileext = ".pdf")
  grDevices::pdf(path)
  on.exit(grDevices::dev.off(), add = TRUE)
  drawn <- NULL
  testthat::expect_error(
    drawn <- plotRateMap(
      bifrost_fit,
      legend = FALSE,
      hold = FALSE,
      palette = "Viridis",
      reverse_palette = TRUE,
      color_mode = "category",
      n_categories = 2,
      category_bin_method = "equal",
      category_breaks = c(-2, 0),
      category_labels = "single",
      uncertainty = TRUE,
      summary = "branch"
    ),
    NA
  )
  testthat::expect_s3_class(drawn, "rateMap")
  testthat::expect_equal(drawn$n_fits, 1L)
  testthat::expect_equal(drawn$category_labels, "single")
  testthat::expect_true(drawn$uncertainty)
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
  testthat::expect_error(rateMap(list(fit), uncertainty = NA), "uncertainty")
  testthat::expect_error(rateMap(list(fit), plot = NA, progress = FALSE), "plot")
  testthat::expect_error(rateMap(list(fit), value_summary = "mode"), "'arg' should be")
  testthat::expect_error(rateMap(list(fit), color_mode = "state"), "'arg' should be")
  testthat::expect_error(rateMap(list(fit), n_categories = 0), "n_categories")
  testthat::expect_error(rateMap(list(fit), quantile_probs = c(0.9, 0.1)), "quantile_probs")
  testthat::expect_error(rateMap(list(fit), highest_density_interval_prob = 0), "highest_density_interval_prob")
  testthat::expect_error(
    rateMap(list(fit), unsupported_arg = TRUE, plot = FALSE, progress = FALSE),
    "Unsupported argument"
  )
  testthat::expect_error(
    rateMap(list(fit), mar = c(0, 0, 0, 0), plot = FALSE, progress = FALSE),
    "only used when 'plot = TRUE'"
  )
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

  generic_target <- ape::as.phylo(target)
  generic_target$edge.length <- c(0.75, 1.25)

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
    target_tree = generic_target,
    check = "topology",
    summary = "branch",
    log = FALSE,
    plot = FALSE,
    progress = FALSE
  )

  testthat::expect_equal(branch_out$target, "user")
  testthat::expect_equal(branch_out$check, "topology")
  testthat::expect_equal(branch_out$tree$edge.length, generic_target$edge.length)
  testthat::expect_equal(branch_out$intervals$value[branch_out$intervals$edge == edge_a], 2)
  testthat::expect_equal(branch_out$intervals$value[branch_out$intervals$edge == edge_b], 3)

  interval_out <- rateMap(
    list(.rate_map_fit(sampled, c("0" = 1, "1" = 3))),
    target_tree = ape::as.phylo(target),
    check = "topology",
    summary = "interval",
    log = FALSE,
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
    log = FALSE,
    plot = FALSE,
    progress = FALSE
  )
  testthat::expect_equal(permuted_out$check, "none")
  testthat::expect_equal(permuted_out$intervals$value[permuted_out$intervals$edge == edge_a], 2)
  testthat::expect_length(permuted_out$clade_key, nrow(target$edge))
  testthat::expect_true("clade_key" %in% names(permuted_out$intervals))
  testthat::expect_equal(dim(permuted_out$edge_matches), c(nrow(target$edge), 1L))
  testthat::expect_equal(colnames(permuted_out$edge_matches), "fit_1")
  testthat::expect_equal(
    permuted_out$clade_key[permuted_out$intervals$edge],
    permuted_out$intervals$clade_key
  )

  mcc_out <- rateMap(
    list(.rate_map_fit(sampled, c("0" = 1, "1" = 3)), .rate_map_fit(permuted, c("0" = 1, "1" = 3))),
    target = "mcc",
    check = "topology",
    summary = "branch",
    log = FALSE,
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
    log = FALSE,
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

test_that("rateMap can color by discrete rate categories", {
  .rate_map_skip_if_missing_deps()

  tr <- .rate_map_two_tip(
    edge_length = c(1, 1),
    maps = list(c("0" = 1), c("1" = 1))
  )
  fit <- .rate_map_fit(tr, c("0" = 1, "1" = 4))

  out <- rateMap(
    list(fit),
    summary = "branch",
    color_mode = "category",
    n_categories = 4,
    palette = c("blue", "red"),
    log = FALSE,
    plot = FALSE,
    progress = FALSE
  )

  testthat::expect_equal(out$color_mode, "category")
  testthat::expect_equal(out$category_labels, c("1", "4"))
  testthat::expect_equal(out$rate_categories$value, c(1, 4))
  testthat::expect_equal(names(out$cols), c("1", "4"))
  testthat::expect_equal(length(out$cols), nrow(out$rate_categories))
  testthat::expect_equal(out$intervals$rate_category, c("1", "4"))
  testthat::expect_equal(
    unname(unlist(lapply(out$tree$maps, names), use.names = FALSE)),
    c("1", "4")
  )

  custom <- rateMap(
    list(fit),
    summary = "branch",
    color_mode = "category",
    category_breaks = c(0, 2, 5),
    category_labels = c("slow", "fast"),
    log = FALSE,
    plot = FALSE,
    progress = FALSE
  )

  testthat::expect_equal(custom$category_labels, c("slow", "fast"))
  testthat::expect_equal(custom$intervals$rate_category, c("slow", "fast"))
  testthat::expect_equal(custom$rate_categories$lower, c(0, 2))
  testthat::expect_equal(custom$rate_categories$upper, c(2, 5))
  testthat::expect_equal(length(custom$cols), nrow(custom$rate_categories))

  equal_binned <- rateMap(
    list(fit),
    summary = "branch",
    color_mode = "category",
    n_categories = 1,
    category_bin_method = "equal",
    log = FALSE,
    plot = FALSE,
    progress = FALSE
  )
  testthat::expect_equal(equal_binned$category_bin_method, "equal")
  testthat::expect_equal(equal_binned$category_breaks, c(1, 4))
  testthat::expect_equal(equal_binned$category_labels, "1 to 4")

  many <- .rateMap_resolve_categories(1:10, n_categories = 3)
  testthat::expect_equal(many$breaks, c(0, 5, 10))
  testthat::expect_equal(many$labels, c("0 to 5", "5 to 10"))
  testthat::expect_equal(many$index, c(1L, 1L, 1L, 1L, 2L, 2L, 2L, 2L, 2L, 2L))

  equal_width <- .rateMap_resolve_categories(
    1:10,
    n_categories = 3,
    category_bin_method = "equal"
  )
  testthat::expect_equal(equal_width$breaks, c(1, 4, 7, 10))
  testthat::expect_equal(equal_width$labels, c("1 to 4", "4 to 7", "7 to 10"))
  testthat::expect_equal(
    equal_width$index,
    c(1L, 1L, 1L, 2L, 2L, 2L, 3L, 3L, 3L, 3L)
  )

  auto_breaks <- .rateMap_resolve_categories(
    c(1, 4),
    n_categories = 2,
    category_breaks = c(0, 2, 5)
  )
  testthat::expect_equal(auto_breaks$labels, c("0 to 2", "2 to 5"))
  testthat::expect_equal(auto_breaks$index, c(1L, 2L))

  testthat::expect_error(
    .rateMap_resolve_categories(c(NA_real_, Inf), n_categories = 2),
    "finite rate categories"
  )

  testthat::expect_error(
    rateMap(
      list(fit),
      summary = "branch",
      color_mode = "category",
      category_bin_method = "quantile",
      log = FALSE,
      plot = FALSE,
      progress = FALSE
    ),
    "'arg' should be"
  )
  testthat::expect_error(
    rateMap(
      list(fit),
      summary = "branch",
      color_mode = "category",
      category_breaks = c(2, 1),
      log = FALSE,
      plot = FALSE,
      progress = FALSE
    ),
    "category_breaks"
  )
  testthat::expect_error(
    rateMap(
      list(fit),
      summary = "branch",
      color_mode = "category",
      category_breaks = c(2, 3),
      log = FALSE,
      plot = FALSE,
      progress = FALSE
    ),
    "cover"
  )
  testthat::expect_error(
    rateMap(
      list(fit),
      summary = "branch",
      color_mode = "category",
      category_breaks = c(0, 2, 5),
      category_labels = c("same", "same"),
      log = FALSE,
      plot = FALSE,
      progress = FALSE
    ),
    "unique"
  )
})

test_that("rateMap computes optional uncertainty summaries", {
  .rate_map_skip_if_missing_deps()
  fx <- .rate_map_fixture()

  fits <- list(
    .rate_map_fit(fx$base, c("0" = 1)),
    .rate_map_fit(fx$base, c("0" = 3))
  )

  out <- rateMap(
    fits,
    summary = "branch",
    uncertainty = TRUE,
    log = FALSE,
    plot = FALSE,
    progress = FALSE
  )

  testthat::expect_true(out$uncertainty)
  testthat::expect_equal(out$value_summary, "mean")
  testthat::expect_length(out$run_values, nrow(fx$base$edge))
  testthat::expect_equal(dim(out$run_values[[1L]]), c(1L, 2L))
  testthat::expect_true(all(c(
    "mean", "median", "sd", "q025", "q975", "highest_density_interval_low", "highest_density_interval_high",
    "ci_width", "highest_density_interval_width", "cv", "n"
  ) %in% names(out$intervals)))
  testthat::expect_equal(out$intervals$value, out$intervals$mean)
  testthat::expect_true(all(abs(out$intervals$mean - 2) < 1e-12))
  testthat::expect_true(all(abs(out$intervals$median - 2) < 1e-12))
  testthat::expect_true(all(abs(out$intervals$sd - sqrt(2)) < 1e-12))
  testthat::expect_true(all(abs(out$intervals$q025 - 1.05) < 1e-12))
  testthat::expect_true(all(abs(out$intervals$q975 - 2.95) < 1e-12))
  testthat::expect_true(all(abs(out$intervals$highest_density_interval_low - 1) < 1e-12))
  testthat::expect_true(all(abs(out$intervals$highest_density_interval_high - 3) < 1e-12))
  testthat::expect_true(all(abs(out$intervals$ci_width - 1.9) < 1e-12))
  testthat::expect_true(all(abs(out$intervals$highest_density_interval_width - 2) < 1e-12))
  testthat::expect_true(all(abs(out$intervals$cv - sqrt(2) / 2) < 1e-12))
  testthat::expect_true(all(out$intervals$n == 2L))
  testthat::expect_equal(
    .rateMap_add_derived_uncertainty(data.frame(mean = 1)),
    data.frame(mean = 1)
  )
  sd_map <- .rateMap_recolor(out, "sd")
  testthat::expect_equal(sd_map$plot_value, "sd")
  testthat::expect_true(diff(sd_map$lims) > 0)

  median_out <- rateMap(
    list(
      .rate_map_fit(fx$base, c("0" = 1)),
      .rate_map_fit(fx$base, c("0" = 5)),
      .rate_map_fit(fx$base, c("0" = 9))
    ),
    summary = "branch",
    value_summary = "median",
    log = FALSE,
    plot = FALSE,
    progress = FALSE
  )

  testthat::expect_false(median_out$uncertainty)
  testthat::expect_null(median_out$run_values)
  testthat::expect_equal(median_out$value_summary, "median")
  testthat::expect_true(all(abs(median_out$intervals$value - 5) < 1e-12))

  weighted <- rateMap(
    fits,
    summary = "branch",
    weights = c(1, 3),
    uncertainty = TRUE,
    value_summary = "median",
    log = FALSE,
    plot = FALSE,
    progress = FALSE
  )

  testthat::expect_equal(weighted$weights, c(0.25, 0.75))
  testthat::expect_true(all(abs(weighted$intervals$mean - 2.5) < 1e-12))
  testthat::expect_true(all(abs(weighted$intervals$value - 3) < 1e-12))
  testthat::expect_true(all(abs(weighted$intervals$median - 3) < 1e-12))
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
    plot = FALSE,
    progress = FALSE
  )

  testthat::expect_equal(out$summary, "branch")
  testthat::expect_true(out$log)
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
    log = FALSE,
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
    log = FALSE,
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

  missing_ic_field <- fit_a
  missing_ic_field$optimal_ic <- NULL
  testthat::expect_error(
    rateMap(
      list(fit_a, missing_ic_field),
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
    log = FALSE,
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
  testthat::expect_error(
    .rateMap_validate_dots(list(TRUE), allowed = "x"),
    "must be named"
  )

  testthat::expect_equal(.rateMap_order_state_names(c("b", "a")), c("a", "b"))
  testthat::expect_true(.rateMap_show_legend(0.1))
  testthat::expect_false(.rateMap_show_legend(FALSE))
  testthat::expect_equal(.rateMap_legend_digits(c(1e-9, 1e-5)), 11L)
  testthat::expect_equal(.rateMap_legend_digits(c(-10, 10)), 3L)
  testthat::expect_equal(.rateMap_legend_digits(c(NA_real_, 1)), 3L)
  testthat::expect_equal(.rateMap_legend_digits(c(0, 0)), 3L)
  testthat::expect_equal(.rateMap_format_rate(1e-9, c(1e-9, 1e-5)), "0.000000001")
  testthat::expect_error(.rateMap_validate_n_categories(NA_real_), "n_categories")
  testthat::expect_error(.rateMap_category_labels(c("a", NA), 2, TRUE), "non-missing")
  testthat::expect_equal(
    .rateMap_category_labels(c("a", "a"), 2, FALSE),
    c("a", "a_1")
  )
  testthat::expect_equal(
    .rateMap_category_lims(list(rate_categories = NULL, lims = c(0, 1))),
    c(0, 1)
  )
  testthat::expect_equal(
    .rateMap_category_lims(list(
      rate_categories = data.frame(value = c(1, 4)),
      lims = c(0, 5)
    )),
    c(1, 4)
  )
  testthat::expect_equal(
    .rateMap_category_lims(list(
      rate_categories = data.frame(lower = 1, upper = 1, value = 1),
      lims = c(0.99, 1.01)
    )),
    c(0.99, 1.01)
  )
  testthat::expect_null(.rateMap_add_category_legend(list(rate_categories = NULL)))
  testthat::expect_equal(
    .rateMap_map_value(c("0" = 1), c("0" = 7), start = 2, end = 3),
    7
  )
  testthat::expect_equal(
    .rateMap_weighted_mean(c(NA, 4), c(1, 3)),
    4
  )
  testthat::expect_true(is.na(.rateMap_weighted_mean(c(NA, Inf), c(1, 1))))
  testthat::expect_true(is.na(.rateMap_weighted_sd(c(1, NA), c(1, 1))))
  testthat::expect_equal(
    .rateMap_weighted_quantile(c(NA, Inf), c(1, 1), c(0.25, 0.75)),
    c(NA_real_, NA_real_)
  )
  testthat::expect_equal(
    .rateMap_weighted_quantile(c(4, NA), c(1, 1), c(0.25, 0.75)),
    c(4, 4)
  )
  testthat::expect_equal(
    .rateMap_weighted_quantile(c(1, 3), c(1, 3), c(0, 1)),
    c(1, 3)
  )
  testthat::expect_equal(
    .rateMap_highest_density_interval(c(NA, Inf), c(1, 1), 0.95),
    c(NA_real_, NA_real_)
  )
  testthat::expect_equal(
    .rateMap_highest_density_interval(c(4, NA), c(1, 1), 0.95),
    c(4, 4)
  )
  testthat::expect_equal(
    .rateMap_highest_density_interval(c(3, 1, 2, NA), c(1, 1, 1, 1), 0.5),
    c(1, 2)
  )
  testthat::expect_equal(
    .rateMap_highest_density_interval(c(1, 2, 100), c(0.49, 0.49, 0.02), 0.9),
    c(1, 2)
  )
  testthat::expect_equal(
    .rateMap_highest_density_interval(c(1, 2, 3), c(1, 1, 1), 1),
    c(1, 3)
  )
  testthat::expect_equal(
    .rateMap_row_means(matrix(c(1, NA, 3, 5), nrow = 2, byrow = TRUE), c(0.25, 0.75)),
    c(1, 4.5)
  )
})

test_that("rateMap can plot directly during computation and use progress/workers", {
  .rate_map_skip_if_missing_deps()
  fx <- .rate_map_fixture()

  path <- tempfile(fileext = ".pdf")
  grDevices::pdf(path)
  on.exit(grDevices::dev.off(), add = TRUE)
  graphics::plot.new()
  graphics::plot.window(xlim = c(0, 1), ylim = c(0, 1))
  testthat::expect_null(.rateMap_add_category_legend(
    list(
      rate_categories = data.frame(lower = 1, upper = 1, value = 1),
      cols = c("1" = "red"),
      title = "one category",
      lims = c(0.99, 1.01)
    ),
    legend = 0.5,
    x_pos = 0.2,
    y_pos = 0.2
  ))

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

test_that("plotRateMap getYmult compatibility shim restores global state", {
  global <- .GlobalEnv
  had_getYmult <- exists("getYmult", envir = global, inherits = FALSE)
  old_getYmult <- if (had_getYmult) {
    get("getYmult", envir = global, inherits = FALSE)
  } else {
    NULL
  }
  on.exit({
    if (exists("getYmult", envir = global, inherits = FALSE)) {
      rm("getYmult", envir = global)
    }
    if (had_getYmult) {
      assign("getYmult", old_getYmult, envir = global)
    }
  }, add = TRUE)

  grDevices::graphics.off()
  no_device <- NULL
  testthat::expect_warning(
    no_device <- .rateMap_getYmult(),
    "No graphics device open"
  )
  testthat::expect_equal(no_device, 1)

  path <- tempfile(fileext = ".pdf")
  grDevices::pdf(path)
  on.exit(if (grDevices::dev.cur() > 1L) grDevices::dev.off(), add = TRUE)
  graphics::plot.new()
  graphics::plot.window(xlim = c(0, 2), ylim = c(0, 4))
  ymult <- .rateMap_getYmult()
  testthat::expect_true(is.finite(ymult))
  testthat::expect_gt(ymult, 0)
  grDevices::dev.off()

  if (exists("getYmult", envir = global, inherits = FALSE)) {
    rm("getYmult", envir = global)
  }

  testthat::expect_true(.rateMap_with_plotrix_getYmult(
    exists("getYmult", envir = global, inherits = FALSE)
  ))
  testthat::expect_false(exists("getYmult", envir = global, inherits = FALSE))

  assign("getYmult", "sentinel", envir = global)
  shimmed <- .rateMap_with_plotrix_getYmult(
    get("getYmult", envir = global, inherits = FALSE)
  )
  testthat::expect_true(is.function(shimmed))
  testthat::expect_identical(get("getYmult", envir = global, inherits = FALSE), "sentinel")
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

  legacy_like <- out
  legacy_like$color_mode <- NULL
  legacy_like$n_categories <- NULL
  legacy_like$category_bin_method <- NULL
  legacy_like$category_breaks <- NULL
  legacy_like$category_labels <- NULL
  legacy_recolored <- .rateMap_recolor(legacy_like, "value")
  testthat::expect_equal(legacy_recolored$color_mode, "continuous")

  binned_like <- out
  binned_like$color_mode <- "category"
  binned_like$n_categories <- 2L
  binned_like$category_breaks <- c(-3, -1, 0)
  binned_like$category_labels <- c("low", "high")
  binned_like$rate_categories <- data.frame(
    color_bin = 1:2,
    rate_category = c("low", "high"),
    lower = c(-3, -1),
    upper = c(-1, 0),
    value = c(NA_real_, NA_real_),
    color = c("black", "white"),
    stringsAsFactors = FALSE
  )
  binned_recolored <- .rateMap_recolor(binned_like, "value")
  testthat::expect_equal(binned_recolored$category_labels, c("low", "high"))

  path <- tempfile(fileext = ".pdf")
  grDevices::pdf(path)
  on.exit(grDevices::dev.off(), add = TRUE)

  drawn <- NULL
  testthat::expect_error(
    drawn <- plotRateMap(out, type = "arc", show_tip_labels = FALSE, legend = FALSE, hold = FALSE),
    NA
  )
  testthat::expect_s3_class(drawn, "rateMap")

  testthat::expect_error(plotRateMap(list()), "rateMap")
})

test_that("plotRateMap covers phylogram, fan, legend, outline, and S3 paths", {
  .rate_map_skip_if_missing_deps()
  fx <- .rate_map_fixture()
  out <- rateMap(
    list(.rate_map_fit(fx$shifted_a), .rate_map_fit(fx$shifted_b)),
    res = 4,
    uncertainty = TRUE,
    plot = FALSE,
    progress = FALSE
  )

  path <- tempfile(fileext = ".pdf")
  grDevices::pdf(path)
  on.exit(grDevices::dev.off(), add = TRUE)

  sd_drawn <- NULL
  testthat::expect_error(
    sd_drawn <- plotRateMap(
      out,
      value = "sd",
      palette = c("white", "black"),
      reverse_palette = TRUE,
      color_mode = "continuous",
      type = "arc",
      legend = FALSE,
      hold = FALSE
    ),
    NA
  )
  testthat::expect_equal(sd_drawn$plot_value, "sd")
  testthat::expect_equal(sd_drawn$intervals$value, sd_drawn$intervals$sd)
  testthat::expect_equal(
    sd_drawn$cols,
    stats::setNames(
      .rateMap_colors(length(out$cols), c("white", "black"), reverse_palette = TRUE),
      as.character(seq_along(out$cols))
    )
  )
  testthat::expect_equal(
    as.numeric(rowSums(sd_drawn$tree$mapped.edge)),
    sd_drawn$tree$edge.length,
    tolerance = 1e-10
  )

  default_recolored <- NULL
  testthat::expect_error(
    default_recolored <- plotRateMap(
      out,
      palette = c("navy", "gold"),
      legend = FALSE,
      hold = FALSE
    ),
    NA
  )
  testthat::expect_equal(default_recolored$plot_value, "value")
  testthat::expect_equal(default_recolored$intervals$value, out$intervals$value)

  category_drawn <- NULL
  testthat::expect_error(
    category_drawn <- plotRateMap(
      out,
      color_mode = "category",
      n_categories = 3,
      category_bin_method = "equal",
      type = "phylogram",
      legend = TRUE,
      hold = FALSE
    ),
    NA
  )
  testthat::expect_equal(category_drawn$color_mode, "category")
  testthat::expect_equal(category_drawn$category_bin_method, "equal")
  testthat::expect_true("rate_category" %in% names(category_drawn$intervals))
  testthat::expect_s3_class(category_drawn$rate_categories, "data.frame")

  continuous_again <- plotRateMap(
    category_drawn,
    color_mode = "continuous",
    legend = FALSE,
    hold = FALSE
  )
  testthat::expect_equal(continuous_again$color_mode, "continuous")
  testthat::expect_false("rate_category" %in% names(continuous_again$intervals))

  continuous_out <- rateMap(
    list(.rate_map_fit(fx$shifted_a), .rate_map_fit(fx$shifted_b)),
    res = 4,
    color_mode = "continuous",
    plot = FALSE,
    progress = FALSE
  )
  highest_density_recolored <- .rateMap_recolor(out, "highest_density_interval_width")
  testthat::expect_equal(
    highest_density_recolored$title,
    "95% highest-density interval width (log rate)"
  )

  raw_out <- rateMap(
    list(.rate_map_fit(fx$shifted_a), .rate_map_fit(fx$shifted_b)),
    res = 4,
    uncertainty = TRUE,
    log = FALSE,
    plot = FALSE,
    progress = FALSE
  )
  raw_highest_density_recolored <- .rateMap_recolor(raw_out, "highest_density_interval_width")
  testthat::expect_equal(
    raw_highest_density_recolored$title,
    "95% highest-density interval width"
  )

  testthat::expect_error(
    plotRateMap(
      continuous_out,
      type = "phylogram",
      legend = 0.2,
      hold = FALSE
    ),
    NA
  )
  testthat::expect_error(
    plotRateMap(
      continuous_out,
      type = "arc",
      legend = 0.2,
      hold = FALSE
    ),
    NA
  )
  testthat::expect_error(
    plotRateMap(
      continuous_out,
      type = "fan",
      outline = TRUE,
      legend = 0.2,
      hold = FALSE
    ),
    NA
  )

  testthat::expect_error(plotRateMap(out, value = "missing", legend = FALSE), "not a column")
  testthat::expect_error(plotRateMap(out, value = "clade_key", legend = FALSE), "numeric column")
  testthat::expect_error(plotRateMap(out, value = NA_character_, legend = FALSE), "single non-empty")
  testthat::expect_error(plotRateMap(out, pallete = "Viridis"), "Unsupported argument")
  testthat::expect_error(plotRateMap(out, legend = FALSE, legend_digits = -1), "legend_digits")
  testthat::expect_error(plotRateMap(out, legend = FALSE, n_categories = 0), "n_categories")
  testthat::expect_error(
    plotRateMap(
      out,
      color_mode = "category",
      category_breaks = c(0, 0.01),
      legend = FALSE
    ),
    "cover"
  )

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

  printed <- utils::capture.output(print(out))
  testthat::expect_true(any(grepl("rateMap object", printed)))
  testthat::expect_true(any(grepl("uncertainty: TRUE", printed)))
  testthat::expect_error(print.rateMap(list()), "rateMap")

  testthat::expect_error(
    plot(
      out,
      value = "sd",
      type = "arc",
      outline = TRUE,
      legend = 0.2,
      show_tip_labels = FALSE,
      hold = FALSE
    ),
    NA
  )
})
