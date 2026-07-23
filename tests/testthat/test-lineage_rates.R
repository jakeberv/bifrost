# tests/testthat/test-lineage_rates.R

.lineage_rates_skip_deps <- function() {
  testthat::skip_if_not_installed("ape")
  testthat::skip_if_not_installed("phytools")
  testthat::skip_if_not_installed("future")
  testthat::skip_if_not_installed("future.apply")
  testthat::skip_if_not_installed("progressr")
}

.lineage_rates_fixture <- function() {
  tr <- ape::read.tree(text = "((a:1,b:1):1,c:2);")
  tr <- phytools::paintSubTree(
    tree = tr,
    node = ape::Ntip(tr) + 1L,
    state = "0"
  )
  phytools::paintSubTree(
    tree = tr,
    node = ape::Ntip(tr) + 2L,
    state = "1"
  )
}

.lineage_rates_fit <- function(tree, rates = c("0" = 1, "1" = 4)) {
  structure(
    list(
      tree_no_uncertainty_untransformed = tree,
      model_no_uncertainty = list(
        model = "BMM",
        call = list(model = "BMM"),
        param = rates
      )
    ),
    class = c("bifrost_search", "list")
  )
}

.lineage_rates_single_fixture <- function() {
  tr <- ape::read.tree(text = "((a:1,b:1):1,c:2);")
  phytools::paintSubTree(
    tree = tr,
    node = ape::Ntip(tr) + 1L,
    state = "0"
  )
}

test_that("lineage_rates keeps numeric tip labels separate from node numbers", {
  .lineage_rates_skip_deps()
  tree <- ape::read.tree(text = "((4:1,b:1):1,c:2);")
  tree <- phytools::paintSubTree(
    tree = tree,
    node = ape::Ntip(tree) + 1L,
    state = "0"
  )
  tree <- phytools::paintSubTree(
    tree = tree,
    node = ape::Ntip(tree) + 2L,
    state = "1"
  )

  out <- lineage_rates(
    tree = tree,
    state_values = c("0" = 1, "1" = 4),
    log = FALSE,
    progress = FALSE
  )

  testthat::expect_equal(out$tip_state, c("1", "1", "0"))
  testthat::expect_equal(out$tip_rate, c(4, 4, 1))
  testthat::expect_equal(out$lineage_rate[[1L]], out$lineage_rate[[2L]])
})

.lineage_rates_hetero_fixture <- function() {
  tr <- ape::read.tree(text = "(((a:1,b:1):1,c:1):1,d:1);")
  tr <- phytools::paintSubTree(
    tree = tr,
    node = ape::Ntip(tr) + 1L,
    state = "0"
  )
  phytools::paintSubTree(
    tree = tr,
    node = ape::Ntip(tr) + 2L,
    state = "1"
  )
}

test_that("lineage_rates accepts bifrost_search output", {
  .lineage_rates_skip_deps()
  fit <- .lineage_rates_fit(.lineage_rates_fixture())

  out <- lineage_rates(
    fit,
    progress = FALSE
  )

  expected_weighted_rate <- {
    raw_weights <- 1 / (2 ^ (c(1, 2) / 5))
    sum(c(4, 1) * raw_weights / sum(raw_weights))
  }

  testthat::expect_s3_class(out, "data.frame")
  testthat::expect_identical(class(out), "data.frame")
  testthat::expect_named(out, c(
    "tip.label",
    "shift_count",
    "log_lineage_rate",
    "tip_state",
    "tip_rate",
    "log_tip_rate"
  ))
  testthat::expect_equal(out$tip.label, c("a", "b", "c"))
  testthat::expect_equal(out$shift_count, c(1, 1, 0))
  testthat::expect_equal(out$tip_state, c("1", "1", "0"))
  testthat::expect_equal(out$tip_rate, c(4, 4, 1))
  testthat::expect_equal(
    out$log_lineage_rate,
    log(c(expected_weighted_rate, expected_weighted_rate, 1))
  )
  testthat::expect_s3_class(attr(out, "branch_metrics"), "data.frame")
  testthat::expect_s3_class(attr(out, "shift_metrics"), "data.frame")
  testthat::expect_named(attr(out, "shift_metrics"), c(
    "Tip",
    "Shift_Count",
    "Total_Time",
    "Shift_Rate_Per_Time",
    "Shift_Rate_Per_Speciation",
    "Tip_State",
    "Tip_State_Value"
  ))
  testthat::expect_named(attr(out, "branch_metrics"), c(
    "Tip",
    "Shift_Count",
    "Total_Time",
    "Shift_Rate_Per_Time",
    "Shift_Rate_Per_Speciation",
    "Weighted_Lineage_Value",
    "Tip_State",
    "Tip_State_Value"
  ))
  testthat::expect_equal(attr(out, "settings")$model_family, "BMM")
  testthat::expect_equal(attr(out, "settings")$age_reference, "tree")

  raw_out <- lineage_rates(
    fit,
    log = FALSE,
    progress = FALSE
  )
  testthat::expect_named(raw_out, c(
    "tip.label",
    "shift_count",
    "lineage_rate",
    "tip_state",
    "tip_rate"
  ))
  testthat::expect_false(any(startsWith(names(raw_out), "log_")))
  testthat::expect_equal(
    raw_out$lineage_rate,
    c(expected_weighted_rate, expected_weighted_rate, 1)
  )
  testthat::expect_false(attr(raw_out, "settings")$log)

  named_out <- lineage_rates(
    bifrost_search = fit,
    progress = FALSE
  )
  testthat::expect_equal(named_out$log_lineage_rate, out$log_lineage_rate)
})

test_that("lineage_rates accepts plain lists with required components", {
  .lineage_rates_skip_deps()
  fit <- unclass(.lineage_rates_fit(.lineage_rates_fixture()))

  out <- lineage_rates(
    fit,
    progress = FALSE
  )

  testthat::expect_equal(names(out)[1L], "tip.label")
  testthat::expect_equal(out$tip.label, c("a", "b", "c"))
  testthat::expect_error(
    lineage_rates("not a fit", progress = FALSE),
    "`bifrost_search` must be a `bifrost_search` object or plain list"
  )
})

test_that("lineage_rates accepts generic input state values through tree", {
  .lineage_rates_skip_deps()
  tree <- .lineage_rates_fixture()

  out <- lineage_rates(
    tree = tree,
    state_values = c("0" = 10, "1" = 40),
    progress = FALSE
  )

  expected_weighted_value <- {
    raw_weights <- 1 / (2 ^ (c(1, 2) / 5))
    sum(c(40, 10) * raw_weights / sum(raw_weights))
  }

  testthat::expect_equal(out$tip.label, c("a", "b", "c"))
  testthat::expect_equal(out$tip_rate, c(40, 40, 10))
  testthat::expect_equal(
    out$log_lineage_rate,
    log(c(expected_weighted_value, expected_weighted_value, 10))
  )
  testthat::expect_equal(attr(out, "settings")$input_mode, "generic_input")
  testthat::expect_true(is.na(attr(out, "settings")$model_family))

  testthat::expect_error(
    lineage_rates(
      x = NULL,
      tree = tree,
      state_values = c("0" = 10, "1" = 40),
      progress = FALSE
    ),
    "unused argument"
  )

  testthat::expect_error(
    lineage_rates(
      tree,
      state_values = c("0" = 10, "1" = 40),
      progress = FALSE
    ),
    "Generic tree inputs must be supplied with `tree`"
  )

  explicit_null_out <- lineage_rates(
    bifrost_search = NULL,
    tree = tree,
    state_values = c("0" = 10, "1" = 40),
    progress = FALSE
  )
  testthat::expect_equal(explicit_null_out$tip.label, c("a", "b", "c"))
  testthat::expect_equal(explicit_null_out$log_lineage_rate, out$log_lineage_rate)
  testthat::expect_equal(attr(explicit_null_out, "settings")$tree_name, "tree")

  raw_out <- lineage_rates(
    tree = tree,
    state_values = c("0" = 0, "1" = -4),
    log = FALSE,
    progress = FALSE
  )
  testthat::expect_false(any(startsWith(names(raw_out), "log_")))
  testthat::expect_equal(raw_out$tip_rate, c(-4, -4, 0))
  testthat::expect_false(attr(raw_out, "settings")$log)

  testthat::expect_error(
    lineage_rates(tree = tree, progress = FALSE),
    "state_values"
  )
  testthat::expect_error(
    lineage_rates(
      tree = tree,
      state_values = c("0" = 10),
      progress = FALSE
    ),
    "missing values"
  )
  testthat::expect_error(
    lineage_rates(
      tree = tree,
      state_values = c("0" = 10, "1" = 0),
      progress = FALSE
    ),
    "strictly positive"
  )
})

test_that("lineage_rates rejects noncanonical root numbering clearly", {
  .lineage_rates_skip_deps()
  tree <- .lineage_rates_fixture()
  n_tip <- ape::Ntip(tree)
  internal_ids <- n_tip + seq_len(ape::Nnode(tree))
  replacement <- stats::setNames(rev(internal_ids), internal_ids)
  is_internal <- tree$edge > n_tip
  tree$edge[is_internal] <- unname(
    replacement[as.character(tree$edge[is_internal])]
  )
  storage.mode(tree$edge) <- "integer"

  testthat::expect_error(
    lineage_rates(
      tree = tree,
      state_values = c("0" = 1, "1" = 4),
      progress = FALSE
    ),
    "root node.*Ntip\\(tree\\) \\+ 1.*renumber"
  )
})

test_that("lineage_rates validates SIMMAP and fitted parameter inputs", {
  .lineage_rates_skip_deps()
  fit <- .lineage_rates_fit(.lineage_rates_fixture())

  testthat::expect_error(
    lineage_rates(
      fit,
      tree = .lineage_rates_fixture(),
      state_values = c("0" = 1, "1" = 4),
      progress = FALSE
    ),
    "only used when `bifrost_search = NULL`"
  )
  testthat::expect_error(
    lineage_rates(fit, log = NA, progress = FALSE),
    "`log` must be TRUE or FALSE"
  )
  testthat::expect_error(
    lineage_rates(fit, decay_base = 0.5, progress = FALSE),
    "`decay_base` must be a finite numeric scalar greater than or equal to 1"
  )
  testthat::expect_error(
    lineage_rates(fit, half_life = 0, progress = FALSE),
    "`half_life` must be NULL or a positive finite numeric scalar"
  )
  testthat::expect_error(
    lineage_rates(fit, normalize_weights = NA, progress = FALSE),
    "`normalize_weights` must be TRUE or FALSE"
  )
  testthat::expect_error(
    lineage_rates(fit, age_reference = "node", progress = FALSE),
    "`age_reference` must be \"tree\" or \"tip\""
  )
  testthat::expect_error(
    lineage_rates(fit, progress = NA),
    "`progress` must be TRUE or FALSE"
  )
  testthat::expect_error(
    lineage_rates(fit, cores = 0, progress = FALSE),
    "`cores` must be a positive integer"
  )
  testthat::expect_error(
    lineage_rates(fit, cores = 1.5, progress = FALSE),
    "`cores` must be a positive integer"
  )
  testthat::expect_error(
    lineage_rates(fit, cores = c(1, 2), progress = FALSE),
    "`cores` must be a positive integer"
  )
  outside_int <- as.double(.Machine$integer.max) + 1
  testthat::expect_error(
    lineage_rates(fit, cores = outside_int, progress = FALSE),
    "`cores` must be.*R integer range"
  )

  missing_tree_fit <- fit
  missing_tree_fit$tree_no_uncertainty_untransformed <- NULL
  testthat::expect_error(
    lineage_rates(missing_tree_fit, progress = FALSE),
    "must include `tree_no_uncertainty_untransformed`"
  )

  malformed_model_fit <- fit
  malformed_model_fit$model_no_uncertainty <- 1
  testthat::expect_error(
    lineage_rates(malformed_model_fit, progress = FALSE),
    "multi-regime BMM"
  )

  empty_model_family_fit <- fit
  empty_model_family_fit$model_no_uncertainty <- list(
    model = "",
    call = list(model = ""),
    param = c("0" = 1, "1" = 4)
  )
  testthat::expect_error(
    lineage_rates(empty_model_family_fit, progress = FALSE),
    "must identify a multi-regime BMM"
  )

  non_brownian_fit <- fit
  non_brownian_fit$model_no_uncertainty$model <- "OU"
  non_brownian_fit$model_no_uncertainty$call$model <- "OU"
  testthat::expect_error(
    lineage_rates(non_brownian_fit, progress = FALSE),
    "multi-regime BMM"
  )

  missing_rates_fit <- fit
  missing_rates_fit$model_no_uncertainty$param <- c("0" = 1)
  testthat::expect_error(
    lineage_rates(
      missing_rates_fit,
      progress = FALSE
    ),
    "missing rate values"
  )

  unnamed_rates_fit <- fit
  unnamed_rates_fit$model_no_uncertainty$param <- c(1, 4)
  testthat::expect_error(
    lineage_rates(
      unnamed_rates_fit,
      progress = FALSE
    ),
    "BMM regime rates"
  )

  duplicated_rates_fit <- fit
  duplicated_rates_fit$model_no_uncertainty$param <- c("0" = 1, "0" = 2, "1" = 4)
  testthat::expect_error(
    lineage_rates(
      duplicated_rates_fit,
      progress = FALSE
    ),
    "names must be unique"
  )

  zero_rate_fit <- fit
  zero_rate_fit$model_no_uncertainty$param <- c("0" = 1, "1" = 0)
  testthat::expect_error(
    lineage_rates(
      zero_rate_fit,
      progress = FALSE
    ),
    "strictly positive"
  )

  plain_tree_fit <- fit
  plain_tree_fit$tree_no_uncertainty_untransformed <- ape::as.phylo(
    plain_tree_fit$tree_no_uncertainty_untransformed
  )
  plain_tree_fit$tree_no_uncertainty_untransformed$maps <- NULL
  testthat::expect_error(
    lineage_rates(plain_tree_fit, progress = FALSE),
    "SIMMAP `maps`"
  )
  bm_fit <- structure(
    list(
      tree_no_uncertainty_untransformed = .lineage_rates_single_fixture(),
      model_no_uncertainty = list(
        model = "BM",
        call = list(model = "BM"),
        param = NA_real_,
        sigma = list(Pinv = diag(c(2, 4)))
      )
    ),
    class = c("bifrost_search", "list")
  )
  testthat::expect_error(
    lineage_rates(bm_fit, progress = FALSE),
    "multi-regime BMM"
  )

  single_state_bmm_fit <- .lineage_rates_fit(
    .lineage_rates_single_fixture(),
    rates = c("0" = 3)
  )
  testthat::expect_error(
    lineage_rates(single_state_bmm_fit, progress = FALSE),
    "at least two mapped states"
  )
})

test_that("lineage_rates respects weighting, progress, and parallel execution settings", {
  .lineage_rates_skip_deps()
  fit <- .lineage_rates_fit(.lineage_rates_fixture())

  unscaled <- lineage_rates(
    fit,
    half_life = NULL,
    normalize_weights = FALSE,
    log = FALSE,
    progress = FALSE
  )
  testthat::expect_s3_class(unscaled, "data.frame")
  testthat::expect_true(all(is.finite(unscaled$lineage_rate)))

  scaled <- lineage_rates(
    fit,
    half_life = 5,
    normalize_weights = FALSE,
    log = FALSE,
    progress = FALSE
  )
  testthat::expect_equal(
    scaled$lineage_rate[1L],
    4 * 2^(-1 / 5) + 2^(-2 / 5)
  )
  testthat::expect_equal(
    unscaled$lineage_rate[1L],
    4 * 2^-1 + 2^-2
  )
  testthat::expect_equal(
    attr(scaled, "shift_metrics")$Shift_Count,
    attr(unscaled, "shift_metrics")$Shift_Count
  )
  testthat::expect_equal(
    attr(scaled, "shift_metrics")$Shift_Rate_Per_Time,
    c(0.5, 0.5, 0)
  )
  testthat::expect_equal(
    attr(scaled, "shift_metrics")$Shift_Rate_Per_Speciation,
    c(1, 1, NA_real_)
  )

  equal_age_weights <- lineage_rates(
    fit,
    decay_base = 1,
    log = FALSE,
    progress = FALSE
  )
  testthat::expect_s3_class(equal_age_weights, "data.frame")
  testthat::expect_equal(attr(equal_age_weights, "settings")$decay_base, 1)

  tree_reference <- lineage_rates(
    fit,
    age_reference = "tree",
    progress = FALSE
  )
  tip_reference <- lineage_rates(
    fit,
    age_reference = "tip",
    progress = FALSE
  )
  testthat::expect_equal(tree_reference$log_lineage_rate, tip_reference$log_lineage_rate)
  testthat::expect_equal(attr(tip_reference, "settings")$age_reference, "tip")

  progress_out <- lineage_rates(
    fit,
    progress = TRUE
  )
  testthat::expect_s3_class(progress_out, "data.frame")

  parallel_out <- lineage_rates(
    fit,
    cores = 2,
    progress = FALSE
  )
  testthat::expect_s3_class(parallel_out, "data.frame")
  testthat::expect_equal(parallel_out, tree_reference, ignore_attr = TRUE)
  testthat::expect_equal(attr(parallel_out, "settings")$cores, 2L)
})

test_that("lineage_rates rejects age-decay weight underflow", {
  .lineage_rates_skip_deps()
  fit <- .lineage_rates_fit(.lineage_rates_fixture())

  testthat::expect_error(
    lineage_rates(
      fit,
      half_life = 1e-300,
      log = FALSE,
      progress = FALSE
    ),
    paste0(
      "Age-decay weights experienced numerical underflow.*",
      "Increase `half_life` or reduce `decay_base`"
    )
  )
})

test_that("lineage_rates can use tip-relative ages for non-ultrametric trees", {
  .lineage_rates_skip_deps()
  fit <- .lineage_rates_fit(.lineage_rates_hetero_fixture())

  normalized_tree <- lineage_rates(
    fit,
    age_reference = "tree",
    progress = FALSE
  )
  normalized_tip <- lineage_rates(
    fit,
    age_reference = "tip",
    progress = FALSE
  )
  testthat::expect_equal(normalized_tree$log_lineage_rate, normalized_tip$log_lineage_rate)

  unnormalized_tree <- lineage_rates(
    fit,
    age_reference = "tree",
    normalize_weights = FALSE,
    log = FALSE,
    progress = FALSE
  )
  unnormalized_tip <- lineage_rates(
    fit,
    age_reference = "tip",
    normalize_weights = FALSE,
    log = FALSE,
    progress = FALSE
  )

  expected_tree_c <- 4 * 2^(-2 / 5) + 1 * 2^(-3 / 5)
  expected_tip_c <- 4 * 2^(-1 / 5) + 1 * 2^(-2 / 5)
  testthat::expect_equal(
    unnormalized_tree$lineage_rate[unnormalized_tree$tip.label == "c"],
    expected_tree_c
  )
  testthat::expect_equal(
    unnormalized_tip$lineage_rate[unnormalized_tip$tip.label == "c"],
    expected_tip_c
  )
  testthat::expect_gt(
    unnormalized_tip$lineage_rate[unnormalized_tip$tip.label == "c"],
    unnormalized_tree$lineage_rate[unnormalized_tree$tip.label == "c"]
  )
})

test_that("lineage_rates reuses linear-time node depths", {
  .lineage_rates_skip_deps()
  fit <- .lineage_rates_fit(.lineage_rates_hetero_fixture())
  tree <- fit$tree_no_uncertainty_untransformed

  expected_tip_depths <- unname(diag(ape::vcv.phylo(tree)))
  expected <- lineage_rates(
    fit,
    log = FALSE,
    progress = FALSE
  )

  testthat::local_mocked_bindings(
    vcv.phylo = function(...) stop("quadratic covariance path used"),
    .package = "ape"
  )
  testthat::local_mocked_bindings(
    nodeheight = function(...) stop("repeated node-height path used"),
    .package = "phytools"
  )

  out <- lineage_rates(
    fit,
    log = FALSE,
    progress = FALSE
  )

  testthat::expect_equal(out, expected)
  testthat::expect_equal(
    attr(out, "branch_metrics")$Total_Time,
    expected_tip_depths
  )
})

test_that("lineage_rates validates generic input mode inputs", {
  .lineage_rates_skip_deps()
  tree <- .lineage_rates_fixture()

  testthat::expect_error(
    lineage_rates(progress = FALSE),
    "Provide either `bifrost_search` or both `tree` and `state_values`"
  )
  testthat::expect_error(
    lineage_rates(tree = list(), state_values = c("0" = 1)),
    "`tree` must be a SIMMAP/phylo tree"
  )
  testthat::expect_error(
    lineage_rates(tree, tree = tree, state_values = c("0" = 1, "1" = 4)),
    "Generic tree inputs must be supplied with `tree`"
  )

  plain_tree <- ape::as.phylo(tree)
  plain_tree$maps <- NULL
  testthat::expect_error(
    lineage_rates(
      tree = plain_tree,
      state_values = c("0" = 1, "1" = 4)
    ),
    "SIMMAP `maps`"
  )

  testthat::expect_error(
    lineage_rates(
      tree = tree,
      state_values = c("0" = Inf, "1" = 4),
      log = FALSE
    ),
    "`state_values` must be finite"
  )
  testthat::expect_error(
    lineage_rates(
      tree = tree,
      state_values = c("0" = 1, "0" = 2, "1" = 4)
    ),
    "`state_values` names must be unique"
  )
})

test_that("compact avian skeleton focal lineage rates are stable", {
  .lineage_rates_skip_deps()

  search <- .avian_skeleton_compact_search()
  provenance <- attr(search, "provenance")
  testthat::expect_type(provenance, "list")
  testthat::expect_match(provenance$source_file, "min10.ic20.gic.RDS", fixed = TRUE)

  out <- lineage_rates(search, progress = FALSE)
  settings <- attr(out, "settings")
  top <- out[
    order(out$log_lineage_rate, decreasing = TRUE),
    c("tip.label", "shift_count", "log_lineage_rate", "tip_state", "tip_rate")
  ]

  testthat::expect_equal(nrow(out), 2057L)
  testthat::expect_equal(settings$input_mode, "bifrost")
  testthat::expect_equal(settings$model_family, "BMM")
  testthat::expect_equal(settings$tree_name, "tree_no_uncertainty_untransformed")
  testthat::expect_false(settings$progress)
  testthat::expect_equal(
    top$tip.label[seq_len(3L)],
    c(
      "Menura_novaehollandiae",
      "Conopophaga_ardesiaca",
      "Conopophaga_lineata"
    )
  )
  testthat::expect_equal(round(top$log_lineage_rate[seq_len(3L)], 3), c(-4.350, -5.750, -5.750))
  testthat::expect_equal(top$shift_count[seq_len(3L)], c(1, 2, 2))
})

test_that("lineage root numbering rejects trees without a unique root", {
  .lineage_rates_skip_deps()

  tree <- ape::read.tree(text = "(a:1,b:1);")
  tree$edge[, 1L] <- c(3L, 4L)

  testthat::expect_error(
    .lineage_rates_validate_root_numbering(tree),
    "Cannot identify a unique root node"
  )
})

test_that("lineage dependency checks identify the missing package", {
  check_packages <- .lineage_rates_check_packages
  environment(check_packages) <- list2env(
    list(requireNamespace = function(package, ...) !identical(package, "ape")),
    parent = environment(.lineage_rates_check_packages)
  )

  testthat::expect_error(check_packages(), "Missing required package\\(s\\): ape")
})

test_that("lineage metric rows resolve an isolated tip-numbered ancestor path", {
  .lineage_rates_skip_deps()

  malformed_tree <- structure(
    list(
      edge = matrix(c(3L, 2L, 2L, 1L), ncol = 2L, byrow = TRUE),
      tip.label = c("a", "b"),
      Nnode = 1L,
      edge.length = c(1, 1),
      maps = list(c("0" = 1), c("1" = 1))
    ),
    class = "phylo"
  )

  testthat::local_mocked_bindings(
    getStates = function(tree, type) {
      if (identical(type, "tips")) c(a = "0", b = "1") else c("3" = "0")
    },
    nodeHeights = function(...) matrix(c(0, 1, 0, 1), ncol = 2L, byrow = TRUE),
    nodeheight = function(...) 0,
    .package = "phytools"
  )
  testthat::local_mocked_bindings(
    vcv.phylo = function(...) diag(c(a = 1, b = 1)),
    .package = "ape"
  )

  out <- .lineage_rates_metric_table(
    malformed_tree,
    mode = "branches",
    state_values = c("0" = 1, "1" = 2),
    cores = 1,
    progress = FALSE
  )

  testthat::expect_equal(out$Shift_Count, c(1, 0))
  testthat::expect_equal(out$Tip_State, c("0", "1"))
})

test_that("lineage metric rows fail cleanly when normalized weights cannot be summed", {
  .lineage_rates_skip_deps()

  metric_table <- .lineage_rates_metric_table
  environment(metric_table) <- list2env(
    list(sum = function(...) Inf),
    parent = environment(.lineage_rates_metric_table)
  )

  testthat::expect_error(
    metric_table(
      .lineage_rates_fixture(),
      mode = "branches",
      state_values = c("0" = 1, "1" = 4),
      cores = 1,
      progress = FALSE
    ),
    "Age-decay weights experienced numerical underflow"
  )
})
