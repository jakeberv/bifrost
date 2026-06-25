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

.rate_map_eval_no_error <- function(expr) {
  expr <- substitute(expr)
  env <- parent.frame()
  value <- NULL
  testthat::expect_error(value <- eval(expr, env), NA)
  value
}

.rate_map_with_pdf <- function(expr) {
  expr <- substitute(expr)
  env <- parent.frame()
  path <- tempfile(fileext = ".pdf")
  grDevices::pdf(path)
  on.exit(grDevices::dev.off(), add = TRUE)
  eval(expr, env)
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

.rate_map_state_fit <- function(rates, edge_length = 1) {
  tr <- ape::stree(length(rates), type = "star")
  tr$edge.length <- rep(edge_length, nrow(tr$edge))
  states <- paste0("s", seq_len(nrow(tr$edge)))
  tr$maps <- lapply(states, function(state) stats::setNames(edge_length, state))
  tr$mapped.edge <- .rateMap_make_mapped_edge(tr$edge, tr$maps)
  class(tr) <- c("simmap", "phylo")
  .rate_map_fit(tr, stats::setNames(rates, states))
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
    progress = FALSE,
    control = rateMapControl(res = 8)
  )

  testthat::expect_s3_class(out, "rateMap")
  testthat::expect_s3_class(out$tree, "simmap")
  testthat::expect_equal(out$summary, "branch")
  testthat::expect_equal(out$color_mode, "category")
  testthat::expect_true(out$log)
  testthat::expect_s3_class(out$rate_categories, "data.frame")
  testthat::expect_true(all(c(
    "n_branches",
    "value_mean",
    "value_median",
    "value_min",
    "value_max",
    "value_sd",
    "total_branch_length"
  ) %in% names(out$rate_categories)))
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
    "rate_flag_source",
    "rate_category"
  ) %in% names(out$intervals)))
  testthat::expect_equal(out$ncolors, 256L)
  testthat::expect_equal(out$rate_flag_source, "value")
  testthat::expect_true(all(out$intervals$rate_flag_source == "value"))
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
    progress = FALSE,
    control = rateMapControl(res = 4)
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
      progress = FALSE
    ),
    "mapped simmap"
  )

  testthat::expect_error(
    rateMap(
      list(list(variables = list(tree = fx$base))),
      progress = FALSE
    ),
    "named numeric"
  )

  testthat::expect_error(
    rateMap(
      list(.rate_map_fit(fx$shifted_a, c("0" = 0.1))),
      progress = FALSE
    ),
    "missing rate parameters"
  )

  testthat::expect_error(
    rateMap(
      list(.rate_map_fit(fx$shifted_a, c("0" = 0.1, "1" = 0))),
      progress = FALSE
    ),
    "non-positive"
  )
  testthat::expect_error(
    rateMap(
      list(.rate_map_fit(fx$shifted_a, c("0" = 0.1, "1" = -0.5))),
      progress = FALSE
    ),
    "strictly positive"
  )
  testthat::expect_error(
    rateMap(
      list(.rate_map_fit(fx$shifted_a, c("0" = 0.1, "1" = 0))),
      log = FALSE,
      progress = FALSE
    ),
    "strictly positive"
  )
  testthat::expect_error(
    rateMap(
      list(.rate_map_fit(fx$shifted_a, c("0" = 0.1, "1" = -0.5))),
      log = FALSE,
      progress = FALSE
    ),
    "strictly positive"
  )
  below_one <- rateMap(
    list(.rate_map_state_fit(c(0.1, 0.2))),
    log = TRUE,
    progress = FALSE
  )
  testthat::expect_true(all(below_one$intervals$value < 0))

  testthat::expect_error(
    rateMap(
      list(.rate_map_fit(fx$shifted_a, c("0" = 0.1, "1" = Inf))),
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
    rateMap(list(invalid), progress = FALSE),
    "named numeric parameter"
  )
  testthat::expect_error(
    rateMap(
      list(invalid),
      progress = FALSE,
      control = rateMapControl(na_action = "omit")
    ),
    "No fits remain"
  )

  out <- rateMap(
    list(invalid, valid),
    progress = FALSE,
    control = rateMapControl(res = 4, na_action = "omit")
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
  bifrost_fit <- structure(
    list(
      tree_no_uncertainty_untransformed = fx$base,
      model_no_uncertainty = list(param = c("0" = 0.4))
    ),
    class = c("bifrost_search", "list")
  )

  out <- rateMap(
    list(mvgls_stub, wrapped_mvgls, bifrost_fit),
    summary = "branch",
    log = FALSE,
    progress = FALSE
  )

  testthat::expect_equal(out$n_fits, 3L)
  testthat::expect_equal(out$weight_mode, "equal")
  testthat::expect_true(all(abs(out$intervals$value - mean(c(0.2, 0.2, 0.4))) < 1e-12))

  transformed_only <- structure(
    list(
      tree_no_uncertainty_transformed = fx$base,
      model_no_uncertainty = list(param = c("0" = 0.4))
    ),
    class = c("bifrost_search", "list")
  )
  testthat::expect_error(
    rateMap(list(transformed_only), progress = FALSE),
    "tree_no_uncertainty_untransformed"
  )
  transformed_override <- rateMap(
    list(transformed_only),
    summary = "branch",
    log = FALSE,
    progress = FALSE,
    control = rateMapControl(
      tree_fun = function(x) x$tree_no_uncertainty_transformed
    )
  )
  testthat::expect_s3_class(transformed_override, "rateMap")
  testthat::expect_equal(transformed_override$n_fits, 1L)
  testthat::expect_true(all(abs(transformed_override$intervals$value - 0.4) < 1e-12))
  omitted_tree <- rateMap(
    list(transformed_only, bifrost_fit),
    summary = "branch",
    log = FALSE,
    progress = FALSE,
    control = rateMapControl(na_action = "omit")
  )
  testthat::expect_equal(omitted_tree$n_fits, 1L)
  testthat::expect_equal(omitted_tree$omitted, 1L)
  testthat::expect_true(all(abs(omitted_tree$intervals$value - 0.4) < 1e-12))
})

test_that("rateMap accepts a single fitted object and plot() draws rateMap objects", {
  .rate_map_skip_if_missing_deps()
  fx <- .rate_map_fixture()

  scratch_fit <- .rate_map_fit(fx$base, c("0" = 0.2))
  scratch_out <- rateMap(
    scratch_fit,
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
    progress = FALSE
  )
  testthat::expect_equal(bifrost_out$n_fits, 1L)
  testthat::expect_true(bifrost_out$log)
  testthat::expect_equal(bifrost_out$intervals$value, rep(log(0.3), nrow(fx$base$edge)))

  bifrost_view <- rateMap(
    bifrost_fit,
    uncertainty = TRUE,
    summary = "branch",
    progress = FALSE
  )
  drawn <- .rate_map_with_pdf({
    .rate_map_eval_no_error(plot(
      bifrost_view,
      legend = FALSE,
      hold = FALSE,
      palette = "Viridis",
      reverse_palette = TRUE,
      color_mode = "category",
      n_categories = 2,
      category_bin_method = "equal",
      category_breaks = c(-2, 0),
      category_labels = "single"
    ))
  })
  testthat::expect_s3_class(drawn, "rateMap")
  testthat::expect_equal(drawn$n_fits, 1L)
  testthat::expect_equal(drawn$category_labels, "single")
  testthat::expect_true(drawn$uncertainty)
  testthat::expect_error(.plot_rate_map(bifrost_fit, legend = FALSE), "rateMap")
})

test_that("rateMap validates top-level arguments and tree checks", {
  .rate_map_skip_if_missing_deps()
  fx <- .rate_map_fixture()
  fit <- .rate_map_fit(fx$base, c("0" = 1))

  testthat::expect_error(rateMap(NULL), "non-empty list")
  testthat::expect_error(rateMapControl(res = 0), "res")
  testthat::expect_error(rateMap(list(fit), workers = 0), "workers")
  testthat::expect_error(rateMap(list(fit), summary = "node"), "'arg' should be")
  testthat::expect_error(rateMap(list(fit), uncertainty = NA), "uncertainty")
  testthat::expect_error(rateMap(list(fit), value_summary = "mode"), "'arg' should be")
  testthat::expect_error(rateMapControl(quantile_probs = c(0.9, 0.1)), "quantile_probs")
  testthat::expect_error(rateMapControl(highest_density_interval_prob = 0), "highest_density_interval_prob")
  testthat::expect_error(
    rateMap(list(fit), unsupported_arg = TRUE, progress = FALSE),
    "unused argument"
  )
  testthat::expect_error(rateMapControl(check = "branch"), "check")
  testthat::expect_error(rateMapControl(tree_fun = 1), "tree_fun")
  testthat::expect_error(rateMapControl(param_fun = 1), "param_fun")
  testthat::expect_error(rateMap(list(fit), control = list(bogus = TRUE)), "Unsupported rateMapControl")
  testthat::expect_error(rateMap(list(fit), control = list(TRUE)), "named list")
  testthat::expect_error(rateMap(list(fit), control = 1), "control")
  testthat::expect_error(rateMap(list(fit), control = NULL, progress = FALSE), "control")
  testthat::expect_error(
    rateMap(list(fit), target_tree = list(), progress = FALSE),
    "target_tree"
  )

  changed <- fx$base
  changed$edge.length[1L] <- changed$edge.length[1L] * 2
  testthat::expect_error(
    rateMap(
      list(fit, .rate_map_fit(changed, c("0" = 1))),
      progress = FALSE
    ),
    "do not match"
  )

  phylo_with_maps <- fx$base
  class(phylo_with_maps) <- "phylo"
  out <- rateMap(
    list(.rate_map_fit(phylo_with_maps, c("0" = 1))),
    summary = "branch",
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
      progress = FALSE
    ),
    "do not match"
  )

  branch_out <- rateMap(
    list(.rate_map_fit(sampled, c("0" = 1, "1" = 3))),
    target_tree = generic_target,
    summary = "branch",
    log = FALSE,
    progress = FALSE,
    control = list(check = "topology")
  )

  testthat::expect_equal(branch_out$target, "user")
  testthat::expect_equal(branch_out$check, "topology")
  testthat::expect_equal(branch_out$tree$edge.length, generic_target$edge.length)
  testthat::expect_equal(branch_out$intervals$value[branch_out$intervals$edge == edge_a], 2)
  testthat::expect_equal(branch_out$intervals$value[branch_out$intervals$edge == edge_b], 3)

  interval_out <- rateMap(
    list(.rate_map_fit(sampled, c("0" = 1, "1" = 3))),
    target_tree = ape::as.phylo(target),
    summary = "interval",
    log = FALSE,
    progress = FALSE,
    control = list(check = "topology", res = 2)
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
    summary = "branch",
    log = FALSE,
    progress = FALSE,
    control = list(check = FALSE)
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
    summary = "branch",
    log = FALSE,
    progress = FALSE,
    control = list(target = "mcc", check = "topology")
  )
  testthat::expect_equal(mcc_out$target, "mcc")
  testthat::expect_equal(mcc_out$tree$edge.length, sampled$edge.length)

  bad_target <- ape::read.tree(text = "(a:1,c:1);")
  testthat::expect_error(
    rateMap(
      list(.rate_map_fit(sampled, c("0" = 1, "1" = 3))),
      target_tree = bad_target,
      summary = "branch",
      progress = FALSE,
      control = list(check = "topology")
    ),
    "topology"
  )
  testthat::expect_error(
    rateMap(
      list(.rate_map_fit(sampled, c("0" = 1, "1" = 3))),
      target_tree = bad_target,
      summary = "branch",
      progress = FALSE,
      control = list(check = FALSE)
    ),
    "missing one or more target-tree branches"
  )
})

test_that("rateMap mapped.edge row sums match edge lengths", {
  .rate_map_skip_if_missing_deps()
  fx <- .rate_map_fixture()
  out <- rateMap(
    list(.rate_map_fit(fx$shifted_a), .rate_map_fit(fx$shifted_a)),
    progress = FALSE,
    control = rateMapControl(res = 10)
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

test_that("rateMap floor flags near-zero branch rates without changing values", {
  .rate_map_skip_if_missing_deps()

  default_no_rule_out <- rateMap(
    list(.rate_map_state_fit(c(0.001, 1, 2))),
    summary = "branch",
    log = FALSE,
    progress = FALSE,
    control = rateMapControl(
      rate_flags = rateMapRateFlags()
    )
  )
  testthat::expect_true(default_no_rule_out$rate_diagnostics$enabled)
  testthat::expect_equal(default_no_rule_out$rate_diagnostics$method, "none")
  testthat::expect_equal(default_no_rule_out$rate_diagnostics$n_near_zero, 0L)
  testthat::expect_equal(default_no_rule_out$rate_diagnostics$full_fold_range, 2000)
  testthat::expect_false("near-zero" %in% default_no_rule_out$rate_categories$rate_category)
  default_no_rule_print <- utils::capture.output(print(default_no_rule_out))
  testthat::expect_true(any(grepl(
    "rate flags: 0 near-zero, 0 high-outlier \\(method: none\\)",
    default_no_rule_print
  )))

  floor_out <- rateMap(
    list(.rate_map_state_fit(c(0.001, 0.005, 1, 2))),
    summary = "branch",
    log = FALSE,
    progress = FALSE,
    control = rateMapControl(
      rate_flags = rateMapRateFlags(
        zero_floor = 0.01,
        zero_color = "grey60"
      )
    )
  )

  testthat::expect_equal(floor_out$intervals$value, c(0.001, 0.005, 1, 2))
  testthat::expect_equal(
    floor_out$intervals$rate_flag,
    c("near-zero", "near-zero", "regular", "regular")
  )
  testthat::expect_equal(
    floor_out$intervals$rate_category,
    c("near-zero", "near-zero", "1", "2")
  )
  testthat::expect_equal(floor_out$rate_categories$rate_category[1L], "near-zero")
  testthat::expect_equal(floor_out$rate_categories$color[1L], "grey60")
  testthat::expect_equal(floor_out$rate_diagnostics$n_near_zero, 2L)
  testthat::expect_equal(floor_out$rate_diagnostics$n_high_outlier, 0L)
  testthat::expect_equal(floor_out$rate_diagnostics$full_fold_range, 2000)
  testthat::expect_equal(floor_out$rate_diagnostics$regular_fold_range, 2)
  testthat::expect_equal(.rateMap_category_lims(floor_out), c(0.001, 2))
  testthat::expect_equal(
    .rateMap_category_marker_values(floor_out)["near-zero"],
    c("near-zero" = 0.01)
  )
  testthat::expect_true(is.na(.rateMap_rate_marker_to_value(NA_real_, log = FALSE)))
  testthat::expect_true(is.na(.rateMap_rate_marker_to_value(0, log = TRUE)))
  floor_print <- utils::capture.output(print(floor_out))
  testthat::expect_true(any(grepl("rate flags: 2 near-zero", floor_print)))
  testthat::expect_true(any(grepl("near-zero floor:", floor_print)))

  palette_out <- rateMapView(
    floor_out,
    palette = c("blue", "red"),
    color_mode = "category"
  )
  regular_category_rows <- palette_out$rate_categories$rate_category != "near-zero"
  testthat::expect_equal(
    unname(palette_out$rate_categories$color[regular_category_rows]),
    grDevices::colorRampPalette(c("blue", "red"))(2)
  )

  log_floor_out <- rateMap(
    list(.rate_map_state_fit(c(1e-10, 1e-5, 1e-4))),
    summary = "branch",
    log = TRUE,
    progress = FALSE,
    control = rateMapControl(
      rate_flags = rateMapRateFlags(method = "floor", zero_floor = 1e-8)
    )
  )
  testthat::expect_equal(
    .rateMap_category_marker_values(log_floor_out)["near-zero"],
    c("near-zero" = log(1e-8))
  )
  plot_file <- tempfile(fileext = ".png")
  grDevices::png(plot_file)
  on.exit(if (grDevices::dev.cur() > 1L) grDevices::dev.off(), add = TRUE)
  plot(floor_out, legend = 1, show_tip_labels = FALSE)
  grDevices::dev.off()

  finite_floor_out <- rateMap(
    list(.rate_map_state_fit(c(0.001, 1, 2))),
    summary = "branch",
    log = FALSE,
    progress = FALSE,
    control = rateMapControl(
      rate_flags = rateMapRateFlags(method = "floor", zero_floor = 0.01)
    )
  )
  finite_floor_print <- utils::capture.output(print(finite_floor_out))
  testthat::expect_true(any(grepl("fold range: 2000 full, 2 regular", finite_floor_print)))

  all_floor_out <- rateMap(
    list(.rate_map_state_fit(c(0.001, 0.005))),
    summary = "branch",
    log = FALSE,
    progress = FALSE,
    control = rateMapControl(
      rate_flags = rateMapRateFlags(method = "floor", zero_floor = 0.01)
    )
  )
  testthat::expect_equal(unique(all_floor_out$intervals$rate_category), "near-zero")
  testthat::expect_true(is.na(all_floor_out$rate_diagnostics$regular_fold_range))
  testthat::expect_equal(nrow(all_floor_out$rate_categories), 1L)

})

test_that("rateMap rate flags support tail clustering, disabled, and interval-mode behavior", {
  .rate_map_skip_if_missing_deps()

  cluster_rates <- exp(c(
    -24, -22, -20,
    -17, -16.5, -16, -15.5, -15, -14.5, -14, -13.5, -13, -12.5, -12
  ))
  cluster_out <- rateMap(
    list(.rate_map_state_fit(cluster_rates)),
    summary = "branch",
    log = FALSE,
    progress = FALSE,
    control = rateMapControl(
      rate_flags = rateMapRateFlags(method = "tail_cluster")
    )
  )
  testthat::expect_equal(cluster_out$rate_diagnostics$n_near_zero, 3L)
  testthat::expect_equal(cluster_out$rate_diagnostics$zero_cluster_n_tail, 3L)
  testthat::expect_equal(cluster_out$rate_diagnostics$zero_cluster_n_regular, 11L)
  testthat::expect_true(is.finite(cluster_out$rate_diagnostics$zero_cluster_log_gap))
  testthat::expect_gt(cluster_out$rate_diagnostics$zero_cluster_log_gap, 0)
  testthat::expect_gt(
    cluster_out$rate_diagnostics$zero_cluster_cutoff_rate,
    max(cluster_rates[1:3])
  )
  testthat::expect_lt(
    cluster_out$rate_diagnostics$zero_cluster_cutoff_rate,
    min(cluster_rates[-(1:3)])
  )
  testthat::expect_true(all(cluster_out$intervals$is_near_zero[1:3]))
  testthat::expect_false(any(cluster_out$intervals$is_near_zero[-(1:3)]))
  testthat::expect_equal(
    .rateMap_category_marker_values(cluster_out)["near-zero"],
    c("near-zero" = cluster_out$rate_diagnostics$zero_cluster_cutoff_rate)
  )

  no_cluster_out <- rateMap(
    list(.rate_map_state_fit(exp(seq(-10, -4, length.out = 14)))),
    summary = "branch",
    log = FALSE,
    progress = FALSE,
    control = rateMapControl(
      rate_flags = rateMapRateFlags(method = "tail_cluster")
    )
  )
  testthat::expect_equal(no_cluster_out$rate_diagnostics$n_near_zero, 0L)
  testthat::expect_true(is.na(no_cluster_out$rate_diagnostics$zero_cluster_cutoff_rate))

  too_small_cluster_out <- rateMap(
    list(.rate_map_state_fit(exp(c(-20, -19, -10, -9)))),
    summary = "branch",
    log = FALSE,
    progress = FALSE,
    control = rateMapControl(
      rate_flags = rateMapRateFlags(method = "tail_cluster")
    )
  )
  testthat::expect_equal(too_small_cluster_out$rate_diagnostics$n_near_zero, 0L)
  testthat::expect_equal(too_small_cluster_out$rate_diagnostics$zero_cluster_n_tail, 0L)
  testthat::expect_equal(too_small_cluster_out$rate_diagnostics$zero_cluster_n_regular, 4L)

  high_cluster_out <- rateMap(
    list(.rate_map_state_fit(exp(c(
      -16, -15.5, -15, -14.5, -14, -13.5, -13, -12.5, -12, -11.5, -11,
      -8, -6, -4
    )))),
    summary = "branch",
    log = FALSE,
    progress = FALSE,
    control = rateMapControl(
      rate_flags = rateMapRateFlags(
        near_zero = FALSE,
        high_outlier = TRUE,
        method = "tail_cluster"
      )
    )
  )
  testthat::expect_equal(high_cluster_out$rate_diagnostics$n_high_outlier, 3L)
  testthat::expect_equal(high_cluster_out$rate_diagnostics$high_cluster_n_tail, 3L)
  testthat::expect_true(all(utils::tail(high_cluster_out$intervals$is_high_outlier, 3)))
  testthat::expect_equal(
    .rateMap_category_marker_values(high_cluster_out)["high-outlier"],
    c("high-outlier" = exp(-9.5))
  )
  plot_file <- tempfile(fileext = ".png")
  grDevices::png(plot_file)
  on.exit(if (grDevices::dev.cur() > 1L) grDevices::dev.off(), add = TRUE)
  plot(
    high_cluster_out,
    legend = 1,
    show_tip_labels = FALSE,
    xlim = rev(range(phytools::nodeHeights(high_cluster_out$tree)))
  )
  grDevices::dev.off()

  disabled <- rateMap(
    list(.rate_map_state_fit(c(1e-9, 1, 2))),
    summary = "branch",
    log = FALSE,
    progress = FALSE,
    control = rateMapControl(rate_flags = NULL)
  )
  testthat::expect_false("rate_flag" %in% names(disabled$intervals))
  testthat::expect_false(disabled$rate_diagnostics$enabled)
  testthat::expect_equal(disabled$rate_diagnostics$reason, "rate flagging is disabled")
  disabled_view <- rateMapView(disabled, color_mode = "category")
  testthat::expect_false(disabled_view$rate_diagnostics$enabled)
  disabled_print <- utils::capture.output(print(disabled))
  testthat::expect_false(any(grepl("rate flags:", disabled_print)))

  interval_out <- rateMap(
    list(.rate_map_state_fit(c(1e-9, 1, 2))),
    summary = "interval",
    log = FALSE,
    progress = FALSE,
    control = rateMapControl(
      res = 2,
      rate_flags = rateMapRateFlags(method = "tail_cluster")
    )
  )
  testthat::expect_false("rate_flag" %in% names(interval_out$intervals))
  testthat::expect_false(interval_out$rate_diagnostics$enabled)
  testthat::expect_equal(
    interval_out$rate_diagnostics$reason,
    "rate flags are only computed for branch summaries"
  )
})

test_that("rateMapView can color by discrete rate categories", {
  .rate_map_skip_if_missing_deps()

  tr <- .rate_map_two_tip(
    edge_length = c(1, 1),
    maps = list(c("0" = 1), c("1" = 1))
  )
  fit <- .rate_map_fit(tr, c("0" = 1, "1" = 4))

  base <- rateMap(
    list(fit),
    summary = "branch",
    log = FALSE,
    progress = FALSE
  )
  out <- rateMapView(
    base,
    color_mode = "category",
    n_categories = 4,
    palette = c("blue", "red")
  )

  testthat::expect_equal(out$color_mode, "category")
  testthat::expect_equal(out$category_labels, c("1", "4"))
  testthat::expect_equal(out$rate_categories$value, c(1, 4))
  testthat::expect_equal(out$rate_categories$n_branches, c(1L, 1L))
  testthat::expect_equal(out$rate_categories$value_mean, c(1, 4))
  testthat::expect_equal(out$rate_categories$value_median, c(1, 4))
  testthat::expect_equal(out$rate_categories$value_min, c(1, 4))
  testthat::expect_equal(out$rate_categories$value_max, c(1, 4))
  testthat::expect_true(all(is.na(out$rate_categories$value_sd)))
  testthat::expect_equal(out$rate_categories$total_branch_length, c(1, 1))
  testthat::expect_equal(names(out$cols), c("1", "4"))
  testthat::expect_equal(length(out$cols), nrow(out$rate_categories))
  testthat::expect_equal(out$intervals$rate_category, c("1", "4"))
  testthat::expect_equal(
    unname(unlist(lapply(out$tree$maps, names), use.names = FALSE)),
    c("1", "4")
  )

  custom <- rateMapView(
    base,
    color_mode = "category",
    category_breaks = c(0, 2, 5),
    category_labels = c("slow", "fast")
  )

  testthat::expect_equal(custom$category_labels, c("slow", "fast"))
  testthat::expect_equal(custom$intervals$rate_category, c("slow", "fast"))
  testthat::expect_equal(custom$rate_categories$lower, c(0, 2))
  testthat::expect_equal(custom$rate_categories$upper, c(2, 5))
  testthat::expect_equal(length(custom$cols), nrow(custom$rate_categories))

  empty_bin <- rateMapView(
    base,
    color_mode = "category",
    category_breaks = c(0, 2, 3, 5),
    category_labels = c("slow", "middle", "fast")
  )
  testthat::expect_equal(empty_bin$rate_categories$n_branches, c(1L, 0L, 1L))
  testthat::expect_equal(empty_bin$rate_categories$value_mean, c(1, NA, 4))
  testthat::expect_equal(empty_bin$rate_categories$value_median, c(1, NA, 4))
  testthat::expect_equal(empty_bin$rate_categories$value_min, c(1, NA, 4))
  testthat::expect_equal(empty_bin$rate_categories$value_max, c(1, NA, 4))
  testthat::expect_true(all(is.na(empty_bin$rate_categories$value_sd)))
  testthat::expect_equal(empty_bin$rate_categories$total_branch_length, c(1, 0, 1))

  equal_binned <- rateMapView(
    base,
    color_mode = "category",
    n_categories = 1,
    category_bin_method = "equal"
  )
  testthat::expect_equal(equal_binned$category_bin_method, "equal")
  testthat::expect_equal(equal_binned$category_breaks, c(1, 4))
  testthat::expect_equal(equal_binned$category_labels, "1 to 4")
  testthat::expect_equal(equal_binned$rate_categories$n_branches, 2L)
  testthat::expect_equal(equal_binned$rate_categories$value_mean, 2.5)
  testthat::expect_equal(equal_binned$rate_categories$value_median, 2.5)
  testthat::expect_equal(equal_binned$rate_categories$value_min, 1)
  testthat::expect_equal(equal_binned$rate_categories$value_max, 4)
  testthat::expect_equal(equal_binned$rate_categories$value_sd, stats::sd(c(1, 4)))
  testthat::expect_equal(equal_binned$rate_categories$total_branch_length, 2)

  near_constant <- rateMap(
    list(.rate_map_state_fit(c(1, 1 + 1e-12))),
    summary = "branch",
    log = FALSE,
    progress = FALSE
  )
  near_constant <- rateMapView(
    near_constant,
    color_mode = "category",
    category_breaks = c(0, 2),
    category_labels = "one"
  )
  testthat::expect_equal(near_constant$rate_categories$n_branches, 2L)
  testthat::expect_equal(near_constant$rate_categories$value_sd, 0)

  continuous_three <- rateMapView(
    base,
    color_mode = "continuous",
    ncolors = 3
  )
  testthat::expect_equal(continuous_three$color_mode, "continuous")
  testthat::expect_length(continuous_three$cols, 3)
  testthat::expect_equal(continuous_three$ncolors, 3L)
  testthat::expect_null(continuous_three$rate_categories)
  continuous_default <- rateMapView(base, color_mode = "continuous")
  testthat::expect_equal(continuous_default$ncolors, 256L)
  testthat::expect_length(continuous_default$cols, 256)
  old_base <- base
  old_base$ncolors <- NULL
  old_continuous <- rateMapView(old_base, color_mode = "continuous")
  testthat::expect_equal(old_continuous$ncolors, 256L)
  testthat::expect_length(old_continuous$cols, 256)
  category_again <- rateMapView(continuous_three, color_mode = "category")
  testthat::expect_equal(category_again$ncolors, 3L)

  interval_base <- rateMap(
    list(fit),
    summary = "interval",
    log = FALSE,
    progress = FALSE,
    control = rateMapControl(res = 2)
  )
  testthat::expect_false("n_branches" %in% names(interval_base$rate_categories))

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
    .rateMap_build_color_map(
      values_by_edge = list(c(1, NA_real_), 2),
      ncolors = 4,
      palette = "YlOrRd",
      reverse_palette = FALSE,
      color_mode = "category",
      n_categories = 2
    ),
    "must all be finite"
  )

  testthat::expect_error(rateMapView(list()), "rateMap")
  testthat::expect_error(
    rateMapView(base, color_mode = "continuous", ncolors = 1),
    "ncolors"
  )
  testthat::expect_error(
    rateMapView(
      base,
      color_mode = "category",
      category_bin_method = "quantile"
    ),
    "'arg' should be"
  )
  testthat::expect_error(
    rateMapView(
      base,
      color_mode = "category",
      category_breaks = c(2, 1)
    ),
    "category_breaks"
  )
  testthat::expect_error(
    rateMapView(
      base,
      color_mode = "category",
      category_breaks = c(2, 3)
    ),
    "cover"
  )
  testthat::expect_error(
    rateMapView(
      base,
      color_mode = "category",
      category_breaks = c(0, 2, 5),
      category_labels = c("same", "same")
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
  single_fit <- rateMap(
    list(.rate_map_fit(fx$base, c("0" = 1))),
    summary = "branch",
    uncertainty = TRUE,
    log = FALSE,
    progress = FALSE
  )
  testthat::expect_error(
    rateMapView(single_fit, value = "sd"),
    "contains no finite values"
  )
  mixed_cv <- rateMap(
    list(
      .rate_map_state_fit(c(1, 2)),
      .rate_map_state_fit(c(1, 4))
    ),
    summary = "branch",
    uncertainty = TRUE,
    log = TRUE,
    progress = FALSE
  )
  testthat::expect_true(any(is.na(mixed_cv$intervals$cv)))
  testthat::expect_true(any(is.finite(mixed_cv$intervals$cv)))
  testthat::expect_error(
    rateMapView(mixed_cv, value = "cv", color_mode = "category"),
    "contains non-finite values"
  )
  testthat::expect_error(
    rateMapView(mixed_cv, value = "cv", color_mode = "continuous"),
    "contains non-finite values"
  )

  flagged_uncertainty <- rateMap(
    list(
      .rate_map_state_fit(c(1e-9, 1, 2)),
      .rate_map_state_fit(c(1e-9, 1.2, 3))
    ),
    summary = "branch",
    uncertainty = TRUE,
    log = FALSE,
    progress = FALSE,
    control = rateMapControl(
      rate_flags = rateMapRateFlags(method = "floor", zero_floor = 1e-8)
    )
  )
  flagged_sd <- rateMapView(flagged_uncertainty, value = "sd", color_mode = "category")
  testthat::expect_equal(flagged_uncertainty$rate_diagnostics$n_near_zero, 1L)
  testthat::expect_equal(flagged_sd$rate_diagnostics$n_near_zero, 1L)
  testthat::expect_equal(flagged_uncertainty$rate_flag_source, "value")
  testthat::expect_equal(flagged_sd$rate_flag_source, "value")
  testthat::expect_true(all(flagged_sd$intervals$rate_flag_source == "value"))
  testthat::expect_true(flagged_sd$intervals$is_near_zero[1L])
  testthat::expect_false("near-zero" %in% flagged_sd$intervals$rate_category)

  median_out <- rateMap(
    list(
      .rate_map_fit(fx$base, c("0" = 1)),
      .rate_map_fit(fx$base, c("0" = 5)),
      .rate_map_fit(fx$base, c("0" = 9))
    ),
    summary = "branch",
    value_summary = "median",
    log = FALSE,
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
      progress = FALSE
    ),
    "finite 'optimal_ic'"
  )

  testthat::expect_error(
    rateMap(
      list(.rate_map_fit(fx$base, c("0" = 1))),
      weights = "ic",
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
    weights = c(2, 2),
    log = FALSE,
    progress = FALSE
  )
  testthat::expect_equal(out$weight_mode, "custom")
  testthat::expect_equal(out$weights, c(0.5, 0.5))
  testthat::expect_true(all(abs(out$intervals$value - 2) < 1e-12))

  testthat::expect_error(
    rateMap(fits, weights = c(1, -1), progress = FALSE),
    "non-negative"
  )
  testthat::expect_error(
    rateMap(fits, weights = c(0, 0), progress = FALSE),
    "positive"
  )
  testthat::expect_error(
    rateMap(fits, weights = c(1, 2, 3), progress = FALSE),
    "length"
  )
  testthat::expect_error(
    rateMap(fits, weights = c(1, NA), progress = FALSE),
    "finite numeric"
  )

  invalid <- .rate_map_fit(fx$base, NA_real_)
  out_omit <- rateMap(
    list(invalid, fits[[1]], fits[[2]]),
    summary = "branch",
    weights = c(100, 1, 3),
    progress = FALSE,
    control = rateMapControl(na_action = "omit")
  )
  testthat::expect_equal(out_omit$weights, c(0.25, 0.75))

  out_omit_used_weights <- rateMap(
    list(invalid, fits[[1]], fits[[2]]),
    summary = "branch",
    weights = c(1, 3),
    progress = FALSE,
    control = rateMapControl(na_action = "omit")
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
  testthat::expect_null(.rateMap_validate_dots(list(x = TRUE), allowed = "x"))

  testthat::expect_equal(.rateMap_order_state_names(c("b", "a")), c("a", "b"))
  testthat::expect_error(.rateMap_normalize_color_options(ncolors = 1), "ncolors")
  testthat::expect_true(.rateMap_show_legend(0.1))
  testthat::expect_false(.rateMap_show_legend(FALSE))
  testthat::expect_equal(.rateMap_legend_digits(c(1e-9, 1e-5)), 11L)
  testthat::expect_equal(.rateMap_legend_digits(c(-10, 10)), 3L)
  testthat::expect_equal(.rateMap_legend_digits(c(NA_real_, 1)), 3L)
  testthat::expect_equal(.rateMap_legend_digits(c(0, 0)), 3L)
  testthat::expect_equal(.rateMap_format_rate(1e-9, c(1e-9, 1e-5)), "0.000000001")
  testthat::expect_error(.rateMap_validate_n_categories(NA_real_), "n_categories")
  testthat::expect_error(.rateMap_category_labels(c("a", NA), 2, TRUE), "non-missing")
  testthat::expect_equal(rateMapRateFlags()$method, "none")
  testthat::expect_false(rateMapRateFlags()$near_zero)
  testthat::expect_false(rateMapRateFlags()$high_outlier)
  flag_control <- rateMapRateFlags(zero_floor = 1e-8)
  testthat::expect_s3_class(flag_control, "rateMap_rate_flags")
  testthat::expect_equal(flag_control$method, "floor")
  testthat::expect_true(flag_control$near_zero)
  testthat::expect_false(flag_control$high_outlier)
  testthat::expect_false(rateMapRateFlags(
    method = "floor",
    zero_floor = 1e-8,
    high_outlier = TRUE
  )$high_outlier)
  tail_control <- rateMapRateFlags(method = "tail_cluster", high_outlier = TRUE)
  testthat::expect_true(tail_control$near_zero)
  testthat::expect_true(tail_control$high_outlier)
  testthat::expect_equal(
    rateMapControl(rate_flags = list(zero_floor = 1e-8))$rate_flags$method,
    "floor"
  )
  testthat::expect_true(rateMapControl(rate_flags = list(zero_floor = 1e-8))$rate_flags$near_zero)
  testthat::expect_equal(
    rateMapControl(rate_flags = list(method = "none"))$rate_flags$method,
    "none"
  )
  testthat::expect_false(rateMapControl(rate_flags = list(method = "none"))$rate_flags$near_zero)
  testthat::expect_error(rateMapRateFlags(near_zero = NA), "near_zero")
  testthat::expect_error(rateMapRateFlags(high_outlier = NA), "high_outlier")
  testthat::expect_error(rateMapRateFlags(method = TRUE), "method")
  testthat::expect_error(rateMapRateFlags(method = "floor"), "zero_floor")
  testthat::expect_error(
    rateMapRateFlags(method = "tail_cluster", zero_floor = 1e-8),
    "zero_floor"
  )
  testthat::expect_error(rateMapRateFlags(method = "none", zero_floor = 1e-8), "zero_floor")
  testthat::expect_error(rateMapRateFlags(zero_floor = -1), "zero_floor")
  testthat::expect_error(rateMapRateFlags(cluster_min_log_gap = 0), "cluster_min_log_gap")
  testthat::expect_error(
    rateMapRateFlags(cluster_min_fold_reduction = 1),
    "cluster_min_fold_reduction"
  )
  testthat::expect_error(
    rateMapRateFlags(cluster_max_tail_fraction = 0.75),
    "cluster_max_tail_fraction"
  )
  testthat::expect_error(rateMapRateFlags(cluster_min_flagged = 0), "cluster_min_flagged")
  testthat::expect_error(rateMapRateFlags(cluster_min_regular = 0), "cluster_min_regular")
  testthat::expect_error(rateMapRateFlags(zero_label = "same", high_label = "same"), "unique")
  testthat::expect_error(rateMapRateFlags(zero_color = ""), "color")
  testthat::expect_error(rateMapControl(rate_flags = TRUE), "rate_flags")
  testthat::expect_error(rateMapControl(rate_flags = list(TRUE)), "rate_flags")
  testthat::expect_error(rateMapControl(rate_flags = list(nope = TRUE)), "Unsupported")
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
  fisher_breaks <- c(-14.22, -13.69, -13.21, -12.91, -12.57, -11.40, -10.57)
  testthat::expect_equal(
    .rateMap_category_legend_breaks(fisher_breaks, 6L),
    fisher_breaks
  )
  xs <- .rateMap_category_legend_x(fisher_breaks, 0, 1)
  testthat::expect_equal(diff(xs), diff(fisher_breaks) / diff(range(fisher_breaks)))
  leftward_xs <- .rateMap_category_legend_x(fisher_breaks, 0, 1, "leftwards")
  testthat::expect_equal(abs(diff(leftward_xs)), diff(xs))
  testthat::expect_equal(leftward_xs[c(1L, length(leftward_xs))], c(1, 0))
  testthat::expect_null(.rateMap_category_legend_breaks(c(1, 4), 2L))
  testthat::expect_null(.rateMap_category_legend_breaks(c(0, 2, 2), 2L))
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

test_that("rateMap uses progress/workers and plot() draws separately", {
  .rate_map_skip_if_missing_deps()
  fx <- .rate_map_fixture()

  .rate_map_with_pdf({
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

    testthat::expect_null(.rateMap_add_category_legend(
      list(
        rate_categories = data.frame(
          lower = c(0, 2),
          upper = c(2, 5),
          value = c(NA_real_, NA_real_)
        ),
        cols = c(slow = "red", fast = "blue"),
        category_breaks = c(0, 2, 5),
        title = "custom categories",
        lims = c(0, 5)
      ),
      legend = 0.5,
      x_pos = 0.2,
      y_pos = 0.3
    ))

    out <- .rate_map_eval_no_error(rateMap(
      list(.rate_map_fit(fx$base, c("0" = 1))),
      summary = "branch",
      progress = TRUE,
      workers = 1,
      control = rateMapControl(future_strategy = "multisession")
    ))
    testthat::expect_s3_class(out, "rateMap")
    .rate_map_eval_no_error(plot(out, legend = FALSE, hold = FALSE))
  })
})

test_that("plot.rateMap getYmult compatibility shim restores global state", {
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

  .rate_map_with_pdf({
    graphics::plot.new()
    graphics::plot.window(xlim = c(0, 2), ylim = c(0, 4))
    ymult <- .rateMap_getYmult()
    testthat::expect_true(is.finite(ymult))
    testthat::expect_gt(ymult, 0)
  })

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

test_that("plot.rateMap draws rate maps on a graphics device", {
  .rate_map_skip_if_missing_deps()
  fx <- .rate_map_fixture()
  out <- rateMap(
    list(.rate_map_fit(fx$shifted_a), .rate_map_fit(fx$shifted_b)),
    progress = FALSE,
    control = rateMapControl(res = 4)
  )

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

  drawn <- .rate_map_with_pdf({
    .rate_map_eval_no_error(plot(
      out,
      type = "arc",
      show_tip_labels = FALSE,
      legend = FALSE,
      hold = FALSE
    ))
  })
  testthat::expect_s3_class(drawn, "rateMap")

  testthat::expect_error(.plot_rate_map(list()), "rateMap")
})

test_that("plot.rateMap covers phylogram, fan, legend, outline, and S3 paths", {
  .rate_map_skip_if_missing_deps()
  fx <- .rate_map_fixture()
  out <- rateMap(
    list(.rate_map_fit(fx$shifted_a), .rate_map_fit(fx$shifted_b)),
    uncertainty = TRUE,
    progress = FALSE,
    control = rateMapControl(res = 4)
  )

  path <- tempfile(fileext = ".pdf")
  grDevices::pdf(path)
  on.exit(grDevices::dev.off(), add = TRUE)

  sd_drawn <- .rate_map_eval_no_error(plot(
      out,
      value = "sd",
      palette = c("white", "black"),
      reverse_palette = TRUE,
      color_mode = "continuous",
      type = "arc",
      legend = FALSE,
      hold = FALSE
  ))
  testthat::expect_equal(sd_drawn$plot_value, "sd")
  testthat::expect_equal(sd_drawn$intervals$value, sd_drawn$intervals$sd)
  testthat::expect_equal(
    sd_drawn$cols,
    stats::setNames(
      .rateMap_colors(out$ncolors, c("white", "black"), reverse_palette = TRUE),
      as.character(seq_len(out$ncolors))
    )
  )
  testthat::expect_equal(sd_drawn$ncolors, out$ncolors)
  testthat::expect_equal(
    as.numeric(rowSums(sd_drawn$tree$mapped.edge)),
    sd_drawn$tree$edge.length,
    tolerance = 1e-10
  )

  default_recolored <- .rate_map_eval_no_error(plot(
      out,
      palette = c("navy", "gold"),
      legend = FALSE,
      hold = FALSE
  ))
  testthat::expect_equal(default_recolored$plot_value, "value")
  testthat::expect_equal(default_recolored$intervals$value, out$intervals$value)

  category_drawn <- .rate_map_eval_no_error(plot(
      out,
      color_mode = "category",
      n_categories = 3,
      category_bin_method = "equal",
      type = "phylogram",
      legend = TRUE,
      hold = FALSE
  ))
  testthat::expect_equal(category_drawn$color_mode, "category")
  testthat::expect_equal(category_drawn$category_bin_method, "equal")
  testthat::expect_true("rate_category" %in% names(category_drawn$intervals))
  testthat::expect_s3_class(category_drawn$rate_categories, "data.frame")

  continuous_again <- plot(
    category_drawn,
    color_mode = "continuous",
    legend = FALSE,
    hold = FALSE
  )
  testthat::expect_equal(continuous_again$color_mode, "continuous")
  testthat::expect_false("rate_category" %in% names(continuous_again$intervals))
  testthat::expect_equal(continuous_again$ncolors, out$ncolors)
  testthat::expect_length(continuous_again$cols, out$ncolors)

  continuous_out <- rateMapView(
    rateMap(
      list(.rate_map_fit(fx$shifted_a), .rate_map_fit(fx$shifted_b)),
      progress = FALSE,
      control = rateMapControl(res = 4)
    ),
    color_mode = "continuous"
  )
  testthat::expect_equal(continuous_out$ncolors, 256L)
  testthat::expect_length(continuous_out$cols, 256)
  highest_density_recolored <- .rateMap_recolor(out, "highest_density_interval_width")
  testthat::expect_equal(
    highest_density_recolored$title,
    "95% highest-density interval width (log rate)"
  )

  raw_out <- rateMap(
    list(.rate_map_fit(fx$shifted_a), .rate_map_fit(fx$shifted_b)),
    uncertainty = TRUE,
    log = FALSE,
    progress = FALSE,
    control = rateMapControl(res = 4)
  )
  raw_highest_density_recolored <- .rateMap_recolor(raw_out, "highest_density_interval_width")
  testthat::expect_equal(
    raw_highest_density_recolored$title,
    "95% highest-density interval width"
  )

  .rate_map_eval_no_error(plot(
      continuous_out,
      type = "phylogram",
      legend = 0.2,
      hold = FALSE
  ))
  .rate_map_eval_no_error(plot(
      continuous_out,
      type = "arc",
      legend = 0.2,
      hold = FALSE
  ))
  .rate_map_eval_no_error(plot(
      continuous_out,
      type = "fan",
      outline = TRUE,
      legend = 0.2,
      hold = FALSE
  ))

  testthat::expect_error(plot(out, value = "missing", legend = FALSE), "not a column")
  testthat::expect_error(plot(out, value = "clade_key", legend = FALSE), "numeric column")
  testthat::expect_error(plot(out, value = NA_character_, legend = FALSE), "single non-empty")
  testthat::expect_error(plot(out, pallete = "Viridis"), "Unsupported argument")
  testthat::expect_error(plot(out, legend = FALSE, legend_digits = -1), "legend_digits")
  testthat::expect_error(plot(out, legend = FALSE, n_categories = 0), "n_categories")
  testthat::expect_error(
    plot(
      out,
      color_mode = "category",
      category_breaks = c(0, 0.01),
      legend = FALSE
    ),
    "cover"
  )

  .rate_map_eval_no_error(plot(out, type = "arc", legend = TRUE, hold = FALSE))
  .rate_map_eval_no_error(plot(out, type = "fan", legend = TRUE, hold = FALSE))

  testthat::expect_message(
    plot(
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

  .rate_map_eval_no_error(plot(
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
  ))

  .rate_map_eval_no_error(plot(
      out,
      type = "fan",
      outline = TRUE,
      legend = 0.2,
      fsize = c(0.4, 0.5),
      hold = FALSE
  ))

  printed <- utils::capture.output(print(out))
  testthat::expect_true(any(grepl("rateMap object", printed)))
  testthat::expect_true(any(grepl("uncertainty: TRUE", printed)))
  testthat::expect_error(print.rateMap(list()), "rateMap")

  .rate_map_eval_no_error(plot(
      out,
      value = "sd",
      type = "arc",
      outline = TRUE,
      legend = 0.2,
      show_tip_labels = FALSE,
      hold = FALSE
  ))
})
