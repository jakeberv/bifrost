.regime_integration_tree <- function() {
  tr <- ape::read.tree(text = "(((a:1,b:1):1,(c:1,d:1):1):1,e:3);")
  tr <- phytools::paintSubTree(
    tree = tr,
    node = ape::Ntip(tr) + 1L,
    state = "root"
  )
  tr <- phytools::paintSubTree(
    tree = tr,
    node = ape::Ntip(tr) + 2L,
    state = "wing"
  )
  phytools::paintSubTree(
    tree = tr,
    node = ape::Ntip(tr) + 3L,
    state = "hind"
  )
}

.regime_integration_data <- function(tree) {
  dat <- data.frame(
    trait1 = seq_along(tree$tip.label),
    trait2 = rev(seq_along(tree$tip.label)),
    row.names = tree$tip.label
  )
  as.matrix(dat)
}

.regime_integration_success_tree <- function() {
  tr <- ape::read.tree(text = paste0(
    "((((((a:1,b:1):1,c:2):1,d:3):1,e:4):1,",
    "((((f:1,g:1):1,h:2):1,i:3):1,j:4):1):1,k:7);"
  ))
  tr <- phytools::paintSubTree(
    tree = tr,
    node = ape::Ntip(tr) + 1L,
    state = "root"
  )
  tr <- phytools::paintSubTree(
    tree = tr,
    node = ape::getMRCA(tr, c("a", "e")),
    state = "wing"
  )
  phytools::paintSubTree(
    tree = tr,
    node = ape::getMRCA(tr, c("f", "j")),
    state = "hind"
  )
}

.regime_integration_success_data <- function(tree) {
  matrix(
    c(
      -0.62, 0.18, -0.84, 1.6, 0.33, -0.82, 0.49, 0.74, 0.58, -0.31, 1.51,
      0.39, -0.62, -2.21, 1.12, -0.04, -0.02, 0.94, 0.82, 0.59, 0.92, 0.78
    ),
    ncol = 2,
    dimnames = list(tree$tip.label, c("trait1", "trait2"))
  )
}

.regime_integration_cor <- function(r12, r13, r23) {
  matrix(
    c(1, r12, r13,
      r12, 1, r23,
      r13, r23, 1),
    nrow = 3,
    byrow = TRUE,
    dimnames = list(c("a", "b", "c"), c("a", "b", "c"))
  )
}

.regime_integration_pca_mats <- function() {
  list(
    r1 = .regime_integration_cor(0.1, 0.2, 0.3),
    r2 = .regime_integration_cor(0.2, 0.1, 0.4),
    r3 = .regime_integration_cor(-0.1, 0.3, 0.2),
    r4 = .regime_integration_cor(0.4, -0.2, 0.1)
  )
}

.regime_integration_relationship_summary <- function() {
  data.frame(
    regime = paste0("r", 1:6),
    rate = c(1.0, 1.4, 1.7, 2.5, 3.1, 4.2),
    mean_variance = c(0.7, 0.9, 1.4, 1.8, 2.6, 3.4),
    mean_abs_correlation = c(0.12, 0.18, 0.23, 0.31, 0.37, 0.44)
  )
}

test_that("fit_regime_covariances records skipped and failed regimes", {
  testthat::skip_if_not_installed("ape")
  testthat::skip_if_not_installed("phytools")
  testthat::skip_if_not_installed("mvMORPH")

  tree <- .regime_integration_tree()
  dat <- .regime_integration_data(tree)

  skipped <- fit_regime_covariances(
    tree = tree,
    trait_data = dat,
    min_tips = 3
  )

  testthat::expect_s3_class(skipped, "regime_covariances")
  testthat::expect_true(all(c("regime", "tip_count", "status", "message") %in% names(skipped$status)))
  testthat::expect_true("skipped" %in% skipped$status$status)
  testthat::expect_equal(skipped$status$tip_count[skipped$status$regime == "wing"], 2)

  failed <- fit_regime_covariances(
    tree = tree,
    trait_data = dat,
    formula = trait_data[, "missing"] ~ 1,
    min_tips = 2
  )

  testthat::expect_true("failed" %in% failed$status$status)
  failed_messages <- failed$status$message[failed$status$status == "failed"]
  testthat::expect_true(all(startsWith(
    failed_messages,
    "Formula/model fit failed for trait_data[, \"missing\"] ~ 1 : "
  )))
})

test_that("fit_regime_covariances defaults to manuscript two-tip refit attempts", {
  testthat::skip_if_not_installed("ape")
  testthat::skip_if_not_installed("phytools")
  testthat::skip_if_not_installed("mvMORPH")

  tree <- .regime_integration_tree()
  dat <- .regime_integration_data(tree)

  out <- fit_regime_covariances(
    tree = tree,
    trait_data = dat,
    formula = trait_data[, "missing"] ~ 1
  )
  status_by_regime <- stats::setNames(out$status$status, out$status$regime)

  testthat::expect_equal(out$min_tips, 2L)
  testthat::expect_equal(unname(status_by_regime["root"]), "skipped")
  testthat::expect_equal(unname(status_by_regime[c("hind", "wing")]), c("failed", "failed"))
})

test_that("fit_regime_covariances can return successful independent fits", {
  testthat::skip_if_not_installed("ape")
  testthat::skip_if_not_installed("phytools")
  testthat::skip_if_not_installed("mvMORPH")

  tree <- .regime_integration_success_tree()
  dat <- .regime_integration_success_data(tree)

  out <- fit_regime_covariances(
    tree = tree,
    trait_data = dat,
    min_tips = 5,
    error = FALSE,
    method = "LL"
  )

  status_by_regime <- stats::setNames(out$status$status, out$status$regime)
  testthat::expect_equal(unname(status_by_regime[c("hind", "wing", "root")]), c("ok", "ok", "skipped"))
  testthat::expect_true(all(vapply(out$covariances[c("hind", "wing")], is.matrix, logical(1))))

  summary <- summarize_regime_covariances(
    out,
    rates = c(hind = 2, root = 0.5, wing = 1)
  )

  rate_by_regime <- stats::setNames(summary$rate, summary$regime)
  testthat::expect_equal(unname(rate_by_regime[c("hind", "root", "wing")]), c(2, 0.5, 1))
  testthat::expect_true(all(is.finite(summary$mean_variance[summary$status == "ok"])))
  testthat::expect_true(all(is.finite(summary$mean_abs_correlation[summary$status == "ok"])))
})

test_that("fit_regime_covariances records unusable fitted matrices as failed", {
  testthat::skip_if_not_installed("ape")
  testthat::skip_if_not_installed("phytools")
  testthat::skip_if_not_installed("mvMORPH")

  tree <- .regime_integration_success_tree()
  dat <- .regime_integration_success_data(tree)
  testthat::local_mocked_bindings(
    .regime_extract_covariance = function(fit) {
      matrix(c(0, 0, 0, 1), nrow = 2)
    },
    .package = "bifrost"
  )

  out <- fit_regime_covariances(
    tree = tree,
    trait_data = dat,
    min_tips = 5,
    error = FALSE,
    method = "LL"
  )

  eligible <- out$status$tip_count >= 5
  testthat::expect_true(all(out$status$status[eligible] == "failed"))
  testthat::expect_true(all(vapply(out$covariances[eligible], is.null, logical(1))))
  testthat::expect_true(all(grepl(
    "positive diagonal variances",
    out$status$message[eligible],
    fixed = TRUE
  )))
})

test_that("fit_regime_covariances records the future parallel strategy", {
  testthat::skip_if_not_installed("ape")
  testthat::skip_if_not_installed("phytools")
  testthat::skip_if_not_installed("mvMORPH")
  testthat::skip_if_not_installed("future.apply")

  tree <- .regime_integration_success_tree()
  dat <- .regime_integration_success_data(tree)
  testthat::local_mocked_bindings(
    plan = function(...) list(),
    .package = "future"
  )
  testthat::local_mocked_bindings(
    future_lapply = function(X, FUN, ...) lapply(X, FUN),
    .package = "future.apply"
  )

  out <- fit_regime_covariances(
    tree = tree,
    trait_data = dat,
    min_tips = 5,
    cores = 2,
    error = FALSE,
    method = "LL"
  )

  testthat::expect_equal(out$cores, 2L)
  testthat::expect_equal(out$parallel_strategy, "future.apply::future_lapply")
  testthat::expect_true(all(vapply(out$covariances[c("hind", "wing")], is.matrix, logical(1))))
})

test_that("fit_regime_covariance_runs fits named runs and forwards parallel settings", {
  testthat::skip_if_not_installed("ape")
  testthat::skip_if_not_installed("phytools")
  testthat::skip_if_not_installed("mvMORPH")
  testthat::skip_if_not_installed("future.apply")

  tree <- .regime_integration_success_tree()
  dat <- .regime_integration_success_data(tree)
  runs <- list(alpha = tree, beta = tree)
  testthat::local_mocked_bindings(
    plan = function(...) list(),
    .package = "future"
  )
  testthat::local_mocked_bindings(
    future_lapply = function(X, FUN, ...) lapply(X, FUN),
    .package = "future.apply"
  )

  out <- fit_regime_covariance_runs(
    runs,
    trait_data = dat,
    min_tips = 5,
    cores = 2,
    error = FALSE,
    method = "LL"
  )

  testthat::expect_s3_class(out, "regime_covariance_runs")
  testthat::expect_named(out, c("alpha", "beta"))
  testthat::expect_true(all(vapply(out, inherits, logical(1), "regime_covariances")))
  testthat::expect_equal(vapply(out, `[[`, integer(1), "cores"), c(alpha = 2L, beta = 2L))
  testthat::expect_true(all(vapply(out, `[[`, character(1), "parallel_strategy") == "future.apply::future_lapply"))
  testthat::expect_true(all(vapply(out, function(x) {
    all(vapply(x$covariances[c("hind", "wing")], is.matrix, logical(1)))
  }, logical(1))))
})

test_that("summarize_regime_covariance_runs summarizes run fits by matching run names", {
  testthat::skip_if_not_installed("ape")
  testthat::skip_if_not_installed("phytools")
  testthat::skip_if_not_installed("mvMORPH")

  tree <- .regime_integration_success_tree()
  dat <- .regime_integration_success_data(tree)
  runs <- list(alpha = tree, beta = tree)
  fits <- fit_regime_covariance_runs(
    runs,
    trait_data = dat,
    min_tips = 5,
    error = FALSE,
    method = "LL"
  )
  rates <- list(
    alpha = c(hind = 2, root = 0.5, wing = 1),
    beta = c(hind = 3, root = 0.75, wing = 1.5)
  )

  out <- summarize_regime_covariance_runs(
    fits,
    rates = rates
  )
  manual_alpha <- summarize_regime_covariances(
    fits$alpha,
    rates = rates$alpha
  )

  testthat::expect_s3_class(out, "regime_covariance_run_summaries")
  testthat::expect_named(out, c("alpha", "beta"))
  testthat::expect_s3_class(out$alpha, "data.frame")
  testthat::expect_equal(out$alpha, manual_alpha)
  testthat::expect_equal(
    stats::setNames(out$beta$rate, out$beta$regime)[c("hind", "root", "wing")],
    rates$beta[c("hind", "root", "wing")]
  )
})

test_that("regime run wrappers reject duplicate normalized names before iteration", {
  fit_calls <- 0L
  summary_calls <- 0L
  testthat::local_mocked_bindings(
    fit_regime_covariances = function(...) {
      fit_calls <<- fit_calls + 1L
      stop("fit loop entered", call. = FALSE)
    },
    summarize_regime_covariances = function(...) {
      summary_calls <<- summary_calls + 1L
      stop("summary loop entered", call. = FALSE)
    },
    .package = "bifrost"
  )

  run <- list(tree_no_uncertainty_untransformed = NULL)
  normalized_collision <- list(run, run)
  names(normalized_collision) <- c("", "run1")
  testthat::expect_error(
    fit_regime_covariance_runs(
      normalized_collision,
      trait_data = matrix(0, nrow = 1L, ncol = 1L)
    ),
    "`x` contains duplicated run name\\(s\\): run1"
  )
  testthat::expect_identical(fit_calls, 0L)

  duplicate_fits <- stats::setNames(list(list(), list()), c("alpha", "alpha"))
  testthat::expect_error(
    summarize_regime_covariance_runs(duplicate_fits),
    "`x` contains duplicated run name\\(s\\): alpha"
  )
  testthat::expect_identical(summary_calls, 0L)
})

test_that("regime covariance run wrappers validate run-shaped inputs", {
  testthat::skip_if_not_installed("ape")
  testthat::skip_if_not_installed("phytools")

  tree <- .regime_integration_success_tree()
  dat <- .regime_integration_success_data(tree)
  search_like <- list(tree_no_uncertainty_untransformed = tree)
  runs <- list(search_like)

  testthat::expect_error(
    fit_regime_covariance_runs(tree, trait_data = dat),
    "list of runs"
  )
  testthat::expect_error(
    fit_regime_covariance_runs(list(), trait_data = dat),
    "non-empty list"
  )
  testthat::expect_error(
    fit_regime_covariance_runs(runs, trait_data = dat, verbose = NA),
    "verbose"
  )
  testthat::expect_message(
    unnamed_fit <- fit_regime_covariance_runs(
      runs,
      trait_data = dat,
      min_tips = 99,
      verbose = TRUE
    ),
    "run1"
  )
  testthat::expect_named(unnamed_fit, "run1")

  mixed_names <- list(search_like, named = search_like)
  names(mixed_names)[[1L]] <- ""
  mixed_fit <- fit_regime_covariance_runs(
    mixed_names,
    trait_data = dat,
    min_tips = 99
  )
  testthat::expect_named(mixed_fit, c("run1", "named"))

  whitespace_names <- list(search_like, named = search_like)
  names(whitespace_names)[[1L]] <- "   "
  whitespace_fit <- fit_regime_covariance_runs(
    whitespace_names,
    trait_data = dat,
    min_tips = 99
  )
  testthat::expect_named(whitespace_fit, c("run1", "named"))

  missing_names <- list(named = search_like, search_like)
  names(missing_names)[[2L]] <- NA_character_
  missing_name_fit <- fit_regime_covariance_runs(
    missing_names,
    trait_data = dat,
    min_tips = 99
  )
  testthat::expect_named(missing_name_fit, c("named", "run2"))

  skipped <- fit_regime_covariance_runs(
    list(first = search_like, second = search_like),
    trait_data = dat,
    min_tips = 99
  )
  testthat::expect_error(
    summarize_regime_covariance_runs(
      skipped,
      rates = list(first = c(hind = 1))
    ),
    "same run names"
  )
  recycled <- summarize_regime_covariance_runs(
    skipped,
    rates = c(hind = 1, root = 2, wing = 3)
  )
  testthat::expect_named(recycled, c("first", "second"))
  positional <- summarize_regime_covariance_runs(
    skipped,
    rates = list(c(hind = 1, root = 2, wing = 3), c(hind = 4, root = 5, wing = 6))
  )
  testthat::expect_named(positional, c("first", "second"))
})

test_that("regime covariance inputs validate trees, traits, and scalar settings", {
  testthat::skip_if_not_installed("ape")
  testthat::skip_if_not_installed("phytools")

  tree <- .regime_integration_tree()
  dat <- .regime_integration_data(tree)
  plain_tree <- tree
  plain_tree$maps <- NULL

  testthat::expect_error(
    fit_regime_covariances(trait_data = dat),
    "Supply a SIMMAP-style tree"
  )
  testthat::expect_error(
    fit_regime_covariances(tree = plain_tree, trait_data = dat),
    "SIMMAP-style"
  )
  testthat::expect_error(
    fit_regime_covariances(tree = tree, trait_data = dat[, 1]),
    "matrix or data frame"
  )
  no_names <- dat
  rownames(no_names) <- NULL
  testthat::expect_error(
    fit_regime_covariances(tree = tree, trait_data = no_names),
    "row names"
  )
  testthat::expect_error(
    fit_regime_covariances(tree = tree, trait_data = dat[-1, , drop = FALSE]),
    "missing rows"
  )
  testthat::expect_error(
    fit_regime_covariances(tree = tree, trait_data = dat, min_tips = 0),
    "min_tips"
  )
  testthat::expect_error(
    fit_regime_covariances(tree = tree, trait_data = dat, cores = 0),
    "cores"
  )
  outside_int <- as.double(.Machine$integer.max) + 1
  testthat::expect_error(
    fit_regime_covariances(tree = tree, trait_data = dat, min_tips = 1.5),
    "`min_tips` must be.*integer"
  )
  testthat::expect_error(
    fit_regime_covariances(tree = tree, trait_data = dat, min_tips = outside_int),
    "`min_tips` must be.*R integer range"
  )
  testthat::expect_error(
    fit_regime_covariances(tree = tree, trait_data = dat, cores = 1.5),
    "`cores` must be.*integer"
  )
  testthat::expect_error(
    fit_regime_covariances(tree = tree, trait_data = dat, cores = outside_int),
    "`cores` must be.*R integer range"
  )
})

test_that("regime covariance fitting rejects ambiguous tree and trait identifiers", {
  testthat::skip_if_not_installed("ape")
  testthat::skip_if_not_installed("phytools")

  tree <- .regime_integration_tree()
  dat <- .regime_integration_data(tree)

  duplicate_tree <- tree
  duplicate_tree$tip.label[[2L]] <- duplicate_tree$tip.label[[1L]]
  testthat::expect_error(
    fit_regime_covariances(tree = duplicate_tree, trait_data = dat),
    "tree tip labels.*duplicated identifier.*a"
  )

  empty_tree <- tree
  empty_tree$tip.label[[1L]] <- ""
  testthat::expect_error(
    fit_regime_covariances(tree = empty_tree, trait_data = dat),
    "tree tip labels.*non-empty"
  )

  blank_state_tree <- tree
  names(blank_state_tree$maps[[1L]])[[1L]] <- ""
  testthat::expect_error(
    fit_regime_covariances(tree = blank_state_tree, trait_data = dat),
    "mapped regime-state identifiers.*non-empty"
  )

  missing_state_tree <- tree
  names(missing_state_tree$maps[[1L]])[[1L]] <- NA_character_
  testthat::expect_error(
    fit_regime_covariances(tree = missing_state_tree, trait_data = dat),
    "mapped regime-state identifiers.*non-empty"
  )

  duplicate_traits <- rbind(dat, dat[1L, , drop = FALSE])
  rownames(duplicate_traits) <- c(rownames(dat), rownames(dat)[[1L]])
  testthat::expect_error(
    fit_regime_covariances(tree = tree, trait_data = duplicate_traits),
    "trait-data row names.*duplicated identifier.*a"
  )

  empty_traits <- dat
  rownames(empty_traits)[[1L]] <- ""
  testthat::expect_error(
    fit_regime_covariances(tree = tree, trait_data = empty_traits),
    "trait-data row names.*non-empty"
  )

  duplicate_columns <- dat
  colnames(duplicate_columns) <- c("trait", "trait")
  testthat::expect_error(
    fit_regime_covariances(tree = tree, trait_data = duplicate_columns),
    "trait-data column names.*duplicated identifier.*trait"
  )

  empty_columns <- dat
  colnames(empty_columns)[[1L]] <- ""
  testthat::expect_error(
    fit_regime_covariances(tree = tree, trait_data = empty_columns),
    "trait-data column names.*non-empty"
  )

  unnamed_columns <- as.matrix(dat)
  colnames(unnamed_columns) <- NULL
  testthat::expect_no_error(
    .regime_validate_trait_data(unnamed_columns, tree)
  )
})

test_that("summarize_regime_covariances computes hand variance and correlation summaries", {
  mats <- list(
    alpha = matrix(c(4, 3, 3, 9), nrow = 2),
    beta = matrix(c(1, 0, 0, 3), nrow = 2)
  )

  out <- summarize_regime_covariances(
    mats,
    rates = c(alpha = 2, beta = 5),
    tip_counts = c(alpha = 4, beta = 6),
    regime_ages = c(alpha = 1.5, beta = 3)
  )

  testthat::expect_s3_class(out, "data.frame")
  testthat::expect_equal(out$regime, c("alpha", "beta"))
  testthat::expect_equal(out$rate, c(2, 5))
  testthat::expect_equal(out$mean_variance, c(6.5, 2))
  testthat::expect_equal(out$mean_abs_correlation, c(0.5, 0))
  testthat::expect_equal(out$fisher_z_mean_abs_correlation[1], atanh(0.5))
  testthat::expect_equal(out$tip_count, c(4, 6))
  testthat::expect_equal(out$regime_age, c(1.5, 3))
  testthat::expect_false(any(is.na(out$status)))
})

test_that("summarize_regime_covariances marks single-trait correlations undefined", {
  single_trait <- list(
    alpha = matrix(
      4,
      nrow = 1L,
      dimnames = list("trait", "trait")
    )
  )

  out <- summarize_regime_covariances(
    single_trait,
    rates = c(alpha = 2)
  )

  testthat::expect_equal(out$mean_variance, 4)
  testthat::expect_identical(out$mean_abs_correlation, NA_real_)
  testthat::expect_identical(out$fisher_z_mean_abs_correlation, NA_real_)
  testthat::expect_equal(out$status, "ok")
})

test_that("regime covariance consumers reject invalid raw matrices by regime", {
  invalid <- matrix(c(0, 0, 0, 1), nrow = 2)

  testthat::expect_error(
    summarize_regime_covariances(list(zero_variance = invalid)),
    "Regime `zero_variance`.*positive diagonal variances"
  )
  testthat::expect_error(
    regime_correlation_pca(list(
      zero_variance = invalid,
      valid_a = diag(2),
      valid_b = diag(c(2, 3))
    )),
    "Regime `zero_variance`.*positive diagonal variances"
  )
  nonfinite <- diag(2)
  nonfinite[1, 2] <- nonfinite[2, 1] <- Inf
  testthat::expect_error(
    summarize_regime_covariances(list(nonfinite = nonfinite)),
    "Regime `nonfinite`.*finite"
  )

  asymmetric <- matrix(c(1, 0.2, 0.8, 1), nrow = 2)
  testthat::expect_error(
    summarize_regime_covariances(list(asymmetric = asymmetric)),
    "Regime `asymmetric`.*symmetric"
  )

  indefinite <- matrix(
    c(1, 2, 2, 1),
    nrow = 2,
    dimnames = list(c("wing", "leg"), c("wing", "leg"))
  )
  testthat::expect_error(
    summarize_regime_covariances(list(indefinite = indefinite)),
    "Regime `indefinite`.*positive semidefinite"
  )

  singular <- matrix(
    c(1, 1, 1, 1),
    nrow = 2,
    dimnames = list(c("wing", "leg"), c("wing", "leg"))
  )
  testthat::expect_no_error(
    summarize_regime_covariances(list(singular = singular))
  )

  nearly_symmetric <- matrix(c(1, 0.2 + 1e-10, 0.2, 1), nrow = 2)
  testthat::expect_no_error(
    summarize_regime_covariances(list(nearly_symmetric = nearly_symmetric))
  )

  duplicate_traits <- diag(2)
  dimnames(duplicate_traits) <- list(c("wing", "wing"), c("wing", "wing"))
  testthat::expect_error(
    summarize_regime_covariances(list(duplicate_traits = duplicate_traits)),
    "Regime `duplicate_traits`.*trait names.*duplicated identifier.*wing"
  )

  mismatched_traits <- diag(2)
  dimnames(mismatched_traits) <- list(c("wing", "leg"), c("wing", "tail"))
  testthat::expect_error(
    summarize_regime_covariances(list(mismatched_traits = mismatched_traits)),
    "Regime `mismatched_traits`.*matching row and column trait names"
  )
})

test_that("regime covariance matching rejects ambiguous regime identifiers", {
  mats <- list(one = diag(2), two = diag(c(2, 3)))

  testthat::expect_error(
    summarize_regime_covariances(unname(mats)),
    "covariance-list regime names.*non-empty"
  )
  names_with_empty <- mats
  names(names_with_empty)[[2L]] <- ""
  testthat::expect_error(
    summarize_regime_covariances(names_with_empty),
    "covariance-list regime names.*non-empty"
  )
  duplicated <- mats
  names(duplicated) <- c("one", "one")
  testthat::expect_error(
    summarize_regime_covariances(duplicated),
    "covariance-list regime names.*duplicated identifier.*one"
  )

  fitted <- structure(
    list(
      covariances = mats,
      status = data.frame(
        regime = c("one", "one"),
        status = c("ok", "ok"),
        message = c(NA_character_, NA_character_)
      ),
      rates = NULL
    ),
    class = c("regime_covariances", "list")
  )
  testthat::expect_error(
    summarize_regime_covariances(fitted),
    "status-table regime IDs.*duplicated identifier.*one"
  )
  fitted$status$regime <- c("one", "")
  testthat::expect_error(
    summarize_regime_covariances(fitted),
    "status-table regime IDs.*non-empty"
  )

  testthat::expect_error(
    summarize_regime_covariances(mats, rates = c(one = 1, one = 2)),
    "rate names.*duplicated identifier.*one"
  )
  empty_rates <- c(1, 2)
  names(empty_rates) <- c("one", "")
  testthat::expect_error(
    summarize_regime_covariances(mats, rates = empty_rates),
    "rate names.*non-empty"
  )
  testthat::expect_error(
    summarize_regime_covariances(mats, tip_counts = c(one = 1, one = 2)),
    "optional-vector regime names.*duplicated identifier.*one"
  )
  empty_optional <- c(1, 2)
  names(empty_optional) <- c("one", "")
  testthat::expect_error(
    summarize_regime_covariances(mats, regime_ages = empty_optional),
    "optional-vector regime names.*non-empty"
  )
})

test_that("invalid matrices in fitted covariance objects become failed rows", {
  fitted <- structure(
    list(
      covariances = list(
        valid = diag(2),
        zero_variance = matrix(c(0, 0, 0, 1), nrow = 2)
      ),
      status = data.frame(
        regime = c("valid", "zero_variance"),
        tip_count = c(5L, 5L),
        regime_age = c(2, 2),
        status = c("ok", "ok"),
        message = c(NA_character_, NA_character_)
      ),
      rates = NULL
    ),
    class = c("regime_covariances", "list")
  )

  out <- summarize_regime_covariances(fitted)

  testthat::expect_equal(out$status, c("ok", "failed"))
  testthat::expect_true(is.na(out$mean_variance[out$regime == "zero_variance"]))
  testthat::expect_match(
    out$message[out$regime == "zero_variance"],
    "positive diagonal variances"
  )
})

test_that("summarize_regime_covariances warns for proportional search VCVs", {
  search <- structure(
    list(
      VCVs = list(
        alpha = matrix(c(4, 1, 1, 2), nrow = 2),
        beta = matrix(c(8, 2, 2, 4), nrow = 2)
      ),
      model_no_uncertainty = list(param = c(alpha = 2, beta = 5))
    ),
    class = c("bifrost_search", "list")
  )

  testthat::expect_warning(
    out <- summarize_regime_covariances(
      search$VCVs,
      search = search
    ),
    "proportional `search\\$VCVs`"
  )
  testthat::expect_equal(out$regime, c("alpha", "beta"))
})

test_that("summarize_regime_covariances can apply manuscript high-correlation filter", {
  mats <- list(
    low = matrix(c(1, 0.5, 0.5, 1), nrow = 2),
    high = matrix(c(1, 0.99, 0.99, 1), nrow = 2)
  )

  out <- summarize_regime_covariances(
    mats,
    rates = c(low = 1, high = 2),
    remove_high_corr = TRUE,
    corr_threshold = 0.95
  )

  testthat::expect_equal(out$regime, "low")
  testthat::expect_true(all(abs(out$mean_abs_correlation) <= 0.95))

  unfiltered <- summarize_regime_covariances(
    mats,
    rates = c(low = 1, high = 2)
  )
  testthat::expect_equal(unfiltered$regime, c("low", "high"))
})

test_that("summarize_regime_covariances matches rates by regime names", {
  mats <- list(
    one = diag(c(1, 2)),
    two = diag(c(3, 4))
  )

  out <- summarize_regime_covariances(
    mats,
    rates = c(two = 20, one = 10)
  )

  testthat::expect_equal(out$regime, c("one", "two"))
  testthat::expect_equal(out$rate, c(10, 20))
})

test_that("regime integration summary helpers reject malformed inputs", {
  mats <- list(one = diag(2), two = diag(2))
  relationship_summary <- .regime_integration_relationship_summary()

  testthat::expect_error(
    summarize_regime_covariances(mats, rates = c(1, 2)),
    "named numeric"
  )
  testthat::expect_error(
    summarize_regime_covariances(mats, rates = c(one = 1, two = 2), corr_threshold = 1),
    "between 0 and 1"
  )
  testthat::expect_error(
    regime_integration_relationships(list()),
    "non-empty list"
  )
  testthat::expect_error(
    regime_integration_relationships(data.frame(rate = 1, vars = 1)),
    "missing required columns"
  )
  testthat::expect_error(
    regime_integration_relationships("not data"),
    "data frame"
  )
  testthat::expect_error(
    regime_integration_relationships(relationship_summary, n_boot = 0),
    "n_boot"
  )
  testthat::expect_error(
    regime_integration_relationships(relationship_summary, ci_level = 1),
    "ci_level"
  )
  testthat::expect_error(
    regime_integration_relationships(
      relationship_summary,
      resid_sd_threshold_vars = c(1, 2),
      n_boot = 1
    ),
    "resid_sd_threshold_vars.*non-negative finite number"
  )
  testthat::expect_error(
    regime_integration_relationships(
      relationship_summary,
      resid_sd_threshold_corrs = NA_real_,
      n_boot = 1
    ),
    "resid_sd_threshold_corrs.*non-negative finite number"
  )
  testthat::expect_error(
    regime_integration_relationships(
      relationship_summary,
      resid_sd_threshold_vars = -1,
      n_boot = 1
    ),
    "resid_sd_threshold_vars.*non-negative finite number"
  )
  testthat::expect_error(
    regime_integration_relationships(
      relationship_summary,
      resid_sd_threshold_corrs = Inf,
      n_boot = 1
    ),
    "resid_sd_threshold_corrs.*non-negative finite number"
  )
  testthat::expect_error(
    regime_integration_relationships(
      relationship_summary,
      resid_sd_threshold_vars = "2",
      n_boot = 1
    ),
    "resid_sd_threshold_vars.*non-negative finite number"
  )
  testthat::expect_identical(
    .regime_validate_resid_sd_threshold(0, "resid_sd_threshold_vars"),
    0
  )
  testthat::expect_error(fisher_z_transform(1), "undefined")

  duplicated_summary <- relationship_summary
  duplicated_summary$regime[[2L]] <- duplicated_summary$regime[[1L]]
  testthat::expect_error(
    regime_integration_relationships(duplicated_summary, n_boot = 1),
    "summary-data regime IDs.*duplicated identifier.*r1"
  )

  empty_summary <- relationship_summary
  empty_summary$regime[[1L]] <- ""
  testthat::expect_error(
    regime_integration_relationships(empty_summary, n_boot = 1),
    "summary-data regime IDs.*non-empty"
  )

  manuscript_summary <- data.frame(
    State = c("r1", "r1", paste0("r", 3:6)),
    rate = relationship_summary$rate,
    vars = relationship_summary$mean_variance,
    corrs = relationship_summary$mean_abs_correlation
  )
  testthat::expect_error(
    regime_integration_relationships(manuscript_summary, n_boot = 1),
    "summary-data regime IDs.*duplicated identifier.*r1"
  )
})

test_that("relationship models preserve rows with incomplete panel data", {
  summary <- .regime_integration_relationship_summary()

  missing_rate <- summary
  missing_rate$rate[[1L]] <- NA_real_
  rate_relationships <- regime_integration_relationships(
    missing_rate,
    resid_sd_threshold_vars = 1000,
    resid_sd_threshold_corrs = 1000,
    n_boot = 1
  )
  testthat::expect_equal(nrow(rate_relationships$combined), nrow(summary))
  testthat::expect_true(is.na(rate_relationships$combined$vars_resid[[1L]]))
  testthat::expect_true(is.na(rate_relationships$combined$corrs_resid[[1L]]))
  testthat::expect_false("r1" %in% rate_relationships$variance_points$regime)
  testthat::expect_false("r1" %in% rate_relationships$correlation_points$regime)

  missing_variance <- summary
  missing_variance$mean_variance[[1L]] <- NA_real_
  variance_relationships <- regime_integration_relationships(
    missing_variance,
    resid_sd_threshold_vars = 1000,
    resid_sd_threshold_corrs = 1000,
    n_boot = 1
  )
  testthat::expect_true(is.na(variance_relationships$combined$vars_resid[[1L]]))
  testthat::expect_false("r1" %in% variance_relationships$variance_points$regime)
  testthat::expect_true("r1" %in% variance_relationships$correlation_points$regime)

  missing_correlation <- summary
  missing_correlation$mean_abs_correlation[[1L]] <- NA_real_
  correlation_relationships <- regime_integration_relationships(
    missing_correlation,
    resid_sd_threshold_vars = 1000,
    resid_sd_threshold_corrs = 1000,
    n_boot = 1
  )
  testthat::expect_true(is.na(correlation_relationships$combined$corrs_resid[[1L]]))
  testthat::expect_true("r1" %in% correlation_relationships$variance_points$regime)
  testthat::expect_false("r1" %in% correlation_relationships$correlation_points$regime)
})

test_that("regime relationship bootstrap controls require usable R integers", {
  summary <- .regime_integration_relationship_summary()
  call_relationships <- function(n_boot = 1, seed = NULL) {
    regime_integration_relationships(
      summary,
      resid_sd_threshold_vars = 1000,
      resid_sd_threshold_corrs = 1000,
      n_boot = n_boot,
      seed = seed
    )
  }
  outside_int <- as.double(.Machine$integer.max) + 1

  testthat::expect_error(call_relationships(n_boot = 1.5), "n_boot.*whole number")
  testthat::expect_error(call_relationships(n_boot = outside_int), "n_boot.*R integer range")
  testthat::expect_error(call_relationships(seed = 1.5), "seed.*whole number")
  testthat::expect_error(call_relationships(seed = outside_int), "seed.*R integer range")

  set.seed(2026)
  seed_before <- .Random.seed
  testthat::expect_warning(
    endpoint_relationships <- call_relationships(seed = .Machine$integer.max),
    NA
  )
  testthat::expect_s3_class(endpoint_relationships, "regime_integration_relationships")
  testthat::expect_equal(.Random.seed, seed_before)

  testthat::expect_warning(
    safe_relationships <- call_relationships(
      seed = as.double(.Machine$integer.max) - 1
    ),
    NA
  )
  testthat::expect_s3_class(safe_relationships, "regime_integration_relationships")
  testthat::expect_equal(.Random.seed, seed_before)
})

test_that("relationship bootstrap streams remain distinct at the seed boundary", {
  bootstrap_seeds <- list()
  testthat::local_mocked_bindings(
    .regime_bootstrap_curve = function(..., seed) {
      bootstrap_seeds[[length(bootstrap_seeds) + 1L]] <<- seed
      data.frame(x = 1, fit = 1, lower = 1, upper = 1)
    },
    .package = "bifrost"
  )

  regime_integration_relationships(
    .regime_integration_relationship_summary(),
    resid_sd_threshold_vars = 1000,
    resid_sd_threshold_corrs = 1000,
    n_boot = 1,
    seed = .Machine$integer.max
  )

  testthat::expect_identical(
    bootstrap_seeds,
    list(.Machine$integer.max, 0L)
  )
})

test_that("regime integration summaries cover matrix-list and legacy status shapes", {
  testthat::skip_if_not_installed("ape")

  single <- summarize_regime_covariances(
    list(one = diag(2)),
    tip_counts = 4,
    regime_ages = 2
  )
  testthat::expect_equal(single$regime, "one")
  testthat::expect_true(is.na(single$rate))
  testthat::expect_equal(single$tip_count, 4)

  mixed <- summarize_regime_covariances(
    list(ok = diag(2), failed = NULL)
  )
  testthat::expect_equal(mixed$status, c("ok", "failed"))
  testthat::expect_true(is.na(mixed$mean_variance[mixed$regime == "failed"]))

  legacy <- structure(
    list(
      covariances = list(one = diag(2)),
      status = data.frame(regime = character()),
      rates = NULL
    ),
    class = c("regime_covariances", "list")
  )
  legacy_summary <- summarize_regime_covariances(legacy)
  testthat::expect_equal(legacy_summary$status, "ok")
  testthat::expect_true(is.na(legacy_summary$message))

  testthat::expect_error(
    summarize_regime_covariances(list(one = matrix(1:6, nrow = 2))),
    "square numeric"
  )
  testthat::expect_error(
    summarize_regime_covariances("bad input"),
    "regime_covariances"
  )
  testthat::expect_error(
    summarize_regime_covariances(list(one = diag(2)), tip_counts = c(1, 2)),
    "one value per regime"
  )
  testthat::expect_s3_class(
    .regime_formula_with_data("trait_data ~ 1", matrix(1, nrow = 1)),
    "formula"
  )
  null_env_formula <- stats::as.formula("trait_data ~ 1")
  environment(null_env_formula) <- NULL
  testthat::expect_s3_class(
    .regime_formula_with_data(null_env_formula, matrix(1, nrow = 1)),
    "formula"
  )
  testthat::expect_equal(
    .regime_extract_covariance(list(sigma = diag(2))),
    diag(2)
  )
  testthat::expect_null(
    .regime_extract_covariance(list(sigma = list()))
  )
  testthat::expect_equal(
    .regime_extract_rates(list(param = c(one = 1))),
    c(one = 1)
  )
  testthat::expect_null(
    .regime_extract_rates(list(param = c(1)))
  )
  testthat::expect_true(is.na(.regime_fisher_z(NA_real_, "NA")))
  testthat::expect_error(
    .regime_standardize_summary_data("bad input"),
    "data frame"
  )
  tree_without_lengths <- ape::read.tree(text = "((a,b),c);")
  tree_without_lengths$edge.length <- NULL
  testthat::expect_true(is.na(
    .regime_subtree_age(tree_without_lengths, c("a", "b"))
  ))
})

test_that("Fisher Z boundary behavior is explicit", {
  boundary <- list(perfect = matrix(c(1, 1, 1, 1), nrow = 2))

  safe <- summarize_regime_covariances(
    boundary,
    rates = c(perfect = 1)
  )
  testthat::expect_true(is.na(safe$fisher_z_mean_abs_correlation))

  testthat::expect_error(
    summarize_regime_covariances(
      boundary,
      rates = c(perfect = 1),
      fisher_boundary = "error"
    ),
    "undefined"
  )
})

test_that("regime_correlation_pca returns variance, scores, and loadings", {
  mats <- .regime_integration_pca_mats()

  pca <- regime_correlation_pca(
    mats,
    tip_counts = c(r1 = 5, r2 = 6, r3 = 7, r4 = 8),
    regime_ages = c(r1 = 1, r2 = 2, r3 = 3, r4 = 4)
  )

  testthat::expect_s3_class(pca, "regime_correlation_pca")
  testthat::expect_equal(dim(pca$vectorized_matrix), c(4, 3))
  testthat::expect_equal(dim(pca$scores), c(4, 3))
  testthat::expect_equal(dim(pca$loadings), c(3, 3))
  testthat::expect_equal(sum(pca$variance_explained), 1)
  testthat::expect_equal(pca$settings$use_correlation, TRUE)
  testthat::expect_null(pca$settings$tip_filter)

  scores <- as.data.frame(pca)
  loadings <- as.data.frame(pca, component = "loadings")
  variance <- as.data.frame(pca, component = "variance")
  diagnostics <- as.data.frame(pca, component = "diagnostics")

  testthat::expect_true(all(c("regime", "PC1", "PC2", "PC3") %in% names(scores)))
  testthat::expect_true(all(c("loading", "PC1") %in% names(loadings)))
  testthat::expect_true(all(c("PC", "variance_explained") %in% names(variance)))
  testthat::expect_equal(diagnostics$tip_count, c(5, 6, 7, 8))
})

test_that("regime_correlation_pca warns for scalar-proportional covariance lists", {
  proportional_vcvs <- list(
    low = matrix(c(4, 1, 1, 2), nrow = 2),
    medium = matrix(c(8, 2, 2, 4), nrow = 2),
    high = matrix(c(12, 3, 3, 6), nrow = 2)
  )

  testthat::expect_warning(
    pca <- regime_correlation_pca(proportional_vcvs),
    "scalar-proportional covariance matrices"
  )
  testthat::expect_s3_class(pca, "regime_correlation_pca")
  testthat::expect_true(pca$settings$proportional_covariances)
})

test_that("regime_correlation_pca handles raw covariance PCA without proportional warnings", {
  diagonal <- diag(3)
  dimnames(diagonal) <- list(letters[1:3], letters[1:3])
  first <- matrix(
    c(1, 0.1, 0.2, 0.1, 2, 0.3, 0.2, 0.3, 3),
    nrow = 3,
    dimnames = list(letters[1:3], letters[1:3])
  )
  second <- matrix(
    c(2, 0.4, -0.1, 0.4, 3, 0.5, -0.1, 0.5, 4),
    nrow = 3,
    dimnames = list(letters[1:3], letters[1:3])
  )

  testthat::expect_silent(
    pca <- regime_correlation_pca(
      list(diagonal = diagonal, first = first, second = second),
      use_correlation = FALSE,
      scale. = FALSE
    )
  )
  testthat::expect_false(pca$settings$proportional_covariances)
  testthat::expect_equal(dim(pca$vectorized_matrix), c(3, 3))

  modules <- list(anterior = c("a", "b"), posterior = "c")
  testthat::expect_error(
    regime_module_diagnostics(pca, modules = modules),
    "requires a correlation-mode PCA.*use_correlation = TRUE"
  )
})

test_that("regime_correlation_pca uses manuscript-style strict min_tips filtering", {
  mats <- list(
    r10 = .regime_integration_cor(0.1, 0.2, 0.3),
    r11 = .regime_integration_cor(0.2, 0.1, 0.4),
    r12 = .regime_integration_cor(-0.1, 0.3, 0.2)
  )

  pca <- regime_correlation_pca(
    mats,
    tip_counts = c(r10 = 10, r11 = 11, r12 = 12),
    min_tips = 10
  )

  testthat::expect_equal(pca$regime_ids, c("r11", "r12"))
  testthat::expect_equal(pca$settings$min_tips, 10)
  testthat::expect_equal(pca$settings$tip_filter, "tip_count > min_tips")
  testthat::expect_error(
    regime_correlation_pca(
      mats,
      tip_counts = c(r10 = 10, r11 = 11, r12 = 12),
      min_tips = c(5, 100)
    ),
    "min_tips.*whole number"
  )
  testthat::expect_error(
    regime_correlation_pca(mats, min_tips = -1),
    "min_tips.*positive integer"
  )
  testthat::expect_error(
    regime_correlation_pca(mats, min_tips = 1.5),
    "min_tips.*whole number"
  )
})

test_that("regime_correlation_pca aligns matrices by trait names", {
  base <- .regime_integration_cor(0.1, 0.2, 0.3)
  permuted <- base[c("c", "a", "b"), c("c", "a", "b")]
  other <- .regime_integration_cor(0.4, -0.2, 0.1)

  pca <- regime_correlation_pca(list(base = base, permuted = permuted, other = other))

  testthat::expect_equal(
    unname(pca$vectorized_matrix["base", ]),
    unname(pca$vectorized_matrix["permuted", ])
  )

  mismatched <- base
  dimnames(mismatched) <- list(c("a", "b", "d"), c("a", "b", "d"))
  testthat::expect_error(
    regime_correlation_pca(list(base = base, mismatched = mismatched)),
    "trait names"
  )
})

test_that("regime_correlation_pca validates inputs and supports plotting modes", {
  mats <- .regime_integration_pca_mats()
  base <- mats$r1
  pca <- regime_correlation_pca(mats)

  testthat::expect_error(regime_correlation_pca(mats, unused = TRUE), "Unused")
  testthat::expect_error(
    regime_correlation_pca(
      mats,
      TRUE,
      TRUE,
      TRUE,
      NULL,
      NULL,
      NULL,
      NULL,
      TRUE
    ),
    "Unused argument\\(s\\): \\.\\.1"
  )
  testthat::expect_error(regime_correlation_pca(mats[1], min_tips = 10), "At least two")
  testthat::expect_error(
    regime_correlation_pca(list(r1 = base, r2 = base[1:2, 1:2], r3 = base)),
    "same dimensions"
  )
  bad_traits <- base
  dimnames(bad_traits) <- list(c("a", "b", "c"), c("a", "b", "d"))
  testthat::expect_error(
    regime_correlation_pca(list(r1 = base, r2 = bad_traits, r3 = base * 0.9)),
    "matching row and column"
  )
  unnamed <- unname(base)
  testthat::expect_warning(
    unnamed_pca <- regime_correlation_pca(list(a = unnamed, b = unnamed * 0.9, c = unnamed * 1.1)),
    "scalar-proportional covariance matrices"
  )
  testthat::expect_equal(unnamed_pca$trait_labels, paste0("trait", 1:3))
  testthat::expect_equal(
    regime_correlation_pca(
      mats,
      trait_labels = c("A", "B", "C")
    )$trait_labels,
    c("A", "B", "C")
  )
  testthat::expect_error(
    regime_correlation_pca(mats, trait_labels = c(a = "A", b = "B")),
    "missing labels"
  )
  testthat::expect_error(
    regime_correlation_pca(mats, trait_labels = c("A", "B")),
    "match the matrix"
  )
  testthat::expect_error(
    regime_correlation_pca(
      mats,
      trait_labels = c(a = "X", b = "X", c = "Y")
    ),
    "trait labels.*duplicated identifier.*X"
  )
  testthat::expect_error(
    regime_correlation_pca(
      mats,
      trait_labels = c(a = "X", b = "   ", c = "Y")
    ),
    "trait labels.*non-empty"
  )
  testthat::expect_error(
    regime_correlation_pca(
      mats,
      trait_labels = c(a = "X", b = NA_character_, c = "Y")
    ),
    "trait labels.*non-empty"
  )

  pdf_file <- tempfile(fileext = ".pdf")
  grDevices::pdf(pdf_file)
  on.exit({
    grDevices::dev.off()
    unlink(pdf_file)
  }, add = TRUE)
  testthat::expect_s3_class(plot(pca, type = "scores"), "regime_correlation_pca")
  testthat::expect_s3_class(plot(pca, type = "variance", components = 1:2), "regime_correlation_pca")
  testthat::expect_s3_class(
    plot(pca, type = "loadings", components = 1:2, cluster = "global", heatmap_engine = "base"),
    "regime_correlation_pca"
  )
  testthat::expect_error(
    plot(pca, type = "loadings", components = "not_a_pc"),
    "Unknown PCA component"
  )
  testthat::expect_error(
    plot(pca, type = "variance", components = 99),
    "valid principal components"
  )
  testthat::expect_error(
    plot(pca, type = "loadings", components = 1:2, main = "one"),
    NA
  )
  testthat::expect_error(
    as.data.frame(pca, unused = TRUE),
    "Unused"
  )
})

test_that("regime PCA plotting supports clustered and zero-loading heatmap paths", {
  mats <- .regime_integration_pca_mats()
  pca <- regime_correlation_pca(mats)

  pdf_file <- tempfile(fileext = ".pdf")
  grDevices::pdf(pdf_file)
  on.exit({
    grDevices::dev.off()
    unlink(pdf_file)
  }, add = TRUE)

  testthat::skip_if_not_installed("ComplexHeatmap")
  testthat::skip_if_not_installed("circlize")
  testthat::expect_s3_class(
    plot(
      pca,
      type = "loadings",
      components = 1:2,
      cluster = "global",
      heatmap_engine = "ComplexHeatmap",
      show_dendrogram = FALSE,
      show_legend = FALSE
    ),
    "regime_correlation_pca"
  )
  testthat::expect_s3_class(
    plot(
      pca,
      type = "loadings",
      components = 1,
      cluster = "local",
      heatmap_engine = "ComplexHeatmap",
      show_dendrogram = FALSE,
      show_legend = FALSE
    ),
    "regime_correlation_pca"
  )
  zero_pca <- pca
  zero_pca$loadings[] <- 0
  testthat::expect_true(
    .regime_pca_plot_loadings_complex_heatmap(
      x = zero_pca,
      component_idx = 1L,
      main = NULL,
      palette = NULL,
      cluster = "none",
      show_dendrogram = FALSE,
      show_legend = FALSE,
      cex.axis = 0.62
    )
  )

  zero <- matrix(
    0,
    nrow = 2,
    ncol = 2,
    dimnames = list(c("a", "b"), c("a", "b"))
  )
  testthat::expect_error(
    .regime_pca_plot_loading_matrix(
      zero,
      main = "Zero loadings",
      palette = c("#3B4CC0", "#F7F7F7", "#B40426"),
      orders = NULL,
      cex.axis = 0.6
    ),
    NA
  )
  single_orders <- .regime_heatmap_orders(
    matrix(1, nrow = 1, ncol = 1)
  )
  testthat::expect_false(single_orders$row_dend)
  testthat::expect_false(single_orders$column_dend)
})

test_that("regime_module_diagnostics summarizes named module comparisons", {
  mats <- .regime_integration_pca_mats()
  pca <- regime_correlation_pca(mats)

  diagnostics <- regime_module_diagnostics(
    pca,
    modules = list(
      anterior = c("a", "b"),
      posterior = "c"
    ),
    comparisons = list(
      within_anterior = c("anterior", "anterior"),
      anterior_posterior = c("anterior", "posterior")
    )
  )

  testthat::expect_s3_class(diagnostics, "regime_module_diagnostics")
  scores <- as.data.frame(diagnostics, component = "scores")
  correlations <- as.data.frame(diagnostics, component = "correlations")
  plot_data <- as.data.frame(diagnostics, component = "plot_data")

  testthat::expect_equal(scores$within_anterior[1], 0.1)
  testthat::expect_equal(scores$anterior_posterior[1], mean(c(0.2, 0.3)))
  testthat::expect_true(all(c("PC", "module_comparison", "correlation") %in% names(correlations)))
  testthat::expect_true(all(c("regime", "PC1", "within_anterior") %in% names(plot_data)))

  pdf_file <- tempfile(fileext = ".pdf")
  grDevices::pdf(pdf_file)
  on.exit({
    grDevices::dev.off()
    unlink(pdf_file)
  }, add = TRUE)
  testthat::expect_s3_class(
    plot(
      diagnostics,
      pc = c("PC1", "PC2"),
      comparison = c("within_anterior", "anterior_posterior"),
      n_boot = 5,
      ci_level = 0.9,
      seed = 2
    ),
    "regime_module_diagnostics"
  )
  set.seed(123)
  seed_before <- .Random.seed
  testthat::expect_s3_class(
    plot(
      diagnostics,
      pc = "PC1",
      comparison = "within_anterior",
      n_boot = 5,
      seed = 3
    ),
    "regime_module_diagnostics"
  )
  testthat::expect_equal(.Random.seed, seed_before)
})

test_that("regime_module_diagnostics handles one and singleton modules", {
  pca <- regime_correlation_pca(.regime_integration_pca_mats())

  one_module <- regime_module_diagnostics(
    pca,
    modules = list(all = c("a", "b")),
    pcs = "PC1"
  )
  testthat::expect_named(one_module$comparisons, "within_all")
  testthat::expect_true(all(is.finite(one_module$module_scores$within_all)))

  singleton_modules <- regime_module_diagnostics(
    pca,
    modules = list(first = "a", second = "b"),
    pcs = "PC1"
  )
  testthat::expect_named(
    singleton_modules$comparisons,
    c("within_first", "within_second", "first_vs_second")
  )
  testthat::expect_true(all(is.na(singleton_modules$module_scores$within_first)))
  testthat::expect_true(all(is.na(singleton_modules$module_scores$within_second)))
  within_correlations <- singleton_modules$correlations$module_comparison %in%
    c("within_first", "within_second")
  testthat::expect_true(all(is.na(
    singleton_modules$correlations$correlation[within_correlations]
  )))

  pdf_file <- tempfile(fileext = ".pdf")
  grDevices::pdf(pdf_file)
  on.exit({
    grDevices::dev.off()
    unlink(pdf_file)
  }, add = TRUE)
  testthat::expect_warning(
    singleton_plot <- plot(
      singleton_modules,
      pc = "PC1",
      comparison = "within_first",
      n_boot = 0
    ),
    NA
  )
  testthat::expect_s3_class(singleton_plot, "regime_module_diagnostics")
})

test_that("regime_module_diagnostics validates modules, comparisons, and plot selections", {
  mats <- .regime_integration_pca_mats()
  pca <- regime_correlation_pca(mats)
  modules <- list(anterior = c("a", "b"), posterior = c("b", "c"))
  diagnostics <- regime_module_diagnostics(pca, modules = modules, comparisons = NULL, pcs = "PC1")

  testthat::expect_s3_class(diagnostics, "regime_module_diagnostics")
  testthat::expect_true(all(c("within_anterior", "anterior_vs_posterior") %in% names(diagnostics$module_scores)))
  testthat::expect_equal(
    as.data.frame(diagnostics, component = "comparisons")$module_comparison,
    names(diagnostics$comparisons)
  )
  testthat::expect_error(regime_module_diagnostics(list(), modules), "regime_correlation_pca")
  testthat::expect_error(regime_module_diagnostics(pca, modules = list()), "non-empty named list")
  testthat::expect_error(regime_module_diagnostics(pca, modules = list(c("a", "b"))), "must be named")
  whitespace_module <- list(c("a", "b"))
  names(whitespace_module) <- "   "
  testthat::expect_error(
    regime_module_diagnostics(pca, modules = whitespace_module),
    "module names.*non-empty"
  )
  missing_module_name <- list(c("a", "b"))
  names(missing_module_name) <- NA_character_
  testthat::expect_error(
    regime_module_diagnostics(pca, modules = missing_module_name),
    "module names.*non-empty"
  )
  duplicated_module_names <- list(c("a", "b"), c("b", "c"))
  names(duplicated_module_names) <- c("same", "same")
  testthat::expect_error(
    regime_module_diagnostics(pca, modules = duplicated_module_names),
    "module names.*duplicated identifier.*same"
  )
  testthat::expect_error(regime_module_diagnostics(pca, modules = list(a = character())), "Empty module")
  testthat::expect_error(
    regime_module_diagnostics(
      pca,
      modules = list(anterior = c("a", "a", "b"), posterior = "c")
    ),
    "Module `anterior`.*duplicated trait label.*a"
  )
  testthat::expect_error(regime_module_diagnostics(pca, modules = list(a = "missing")), "not present")
  testthat::expect_error(
    regime_module_diagnostics(pca, modules = modules, comparisons = list()),
    "non-empty"
  )
  mixed_comparison_names <- list(
    custom = c("anterior", "posterior"),
    c("posterior", "anterior")
  )
  names(mixed_comparison_names)[[2L]] <- "   "
  mixed_diagnostics <- regime_module_diagnostics(
    pca,
    modules = modules,
    comparisons = mixed_comparison_names,
    pcs = "PC1"
  )
  testthat::expect_named(
    mixed_diagnostics$comparisons,
    c("custom", "posterior_vs_anterior")
  )
  duplicated_comparison_names <- list(
    c("anterior", "posterior"),
    c("posterior", "anterior")
  )
  names(duplicated_comparison_names) <- c("same", "same")
  testthat::expect_error(
    regime_module_diagnostics(
      pca,
      modules = modules,
      comparisons = duplicated_comparison_names
    ),
    "comparison names.*duplicated identifier.*same"
  )
  generated_name_collision <- list(
    c("anterior", "posterior"),
    c("anterior", "posterior")
  )
  names(generated_name_collision) <- c("", NA_character_)
  testthat::expect_error(
    regime_module_diagnostics(
      pca,
      modules = modules,
      comparisons = generated_name_collision
    ),
    "comparison names.*duplicated identifier.*anterior_vs_posterior"
  )
  testthat::expect_error(
    regime_module_diagnostics(pca, modules = modules, comparisons = list(bad = "anterior")),
    "exactly two"
  )
  testthat::expect_error(
    regime_module_diagnostics(pca, modules = modules, comparisons = list(bad = c("anterior", "missing"))),
    "unknown module"
  )
  testthat::expect_error(
    as.data.frame(diagnostics, component = "scores", unused = TRUE),
    "Unused"
  )

  pdf_file <- tempfile(fileext = ".pdf")
  grDevices::pdf(pdf_file)
  on.exit({
    grDevices::dev.off()
    unlink(pdf_file)
  }, add = TRUE)
  testthat::expect_error(
    plot.regime_module_diagnostics(
      list(),
      pc = "PC1",
      comparison = "within_anterior"
    ),
    "regime_module_diagnostics"
  )
  testthat::expect_error(plot(diagnostics), "Supply both")
  testthat::expect_error(plot(diagnostics, pc = "PC9", comparison = "within_anterior"), "Unknown PC")
  testthat::expect_error(plot(diagnostics, pc = "PC1", comparison = "missing"), "Unknown module")
  testthat::expect_error(
    plot(diagnostics, pc = "PC1", comparison = "within_anterior", n_boot = -1),
    "n_boot"
  )
  outside_int <- as.double(.Machine$integer.max) + 1
  testthat::expect_error(
    plot(diagnostics, pc = "PC1", comparison = "within_anterior", n_boot = 1.5),
    "n_boot.*whole number"
  )
  testthat::expect_error(
    plot(diagnostics, pc = "PC1", comparison = "within_anterior", n_boot = outside_int),
    "n_boot.*R integer range"
  )
  testthat::expect_error(
    plot(diagnostics, pc = "PC1", comparison = "within_anterior", n_boot = 1, seed = 1.5),
    "seed.*whole number"
  )
  testthat::expect_error(
    plot(diagnostics, pc = "PC1", comparison = "within_anterior", n_boot = 1, seed = outside_int),
    "seed.*R integer range"
  )
  set.seed(2026)
  seed_before <- .Random.seed
  graphics_before <- graphics::par(c("mfrow", "usr"))
  testthat::expect_warning(
    testthat::expect_error(
      plot(
        diagnostics,
        pc = c("PC1", "PC1"),
        comparison = c("within_anterior", "within_posterior"),
        n_boot = 1,
        seed = .Machine$integer.max
      ),
      "seed.*R integer range"
    ),
    NA
  )
  testthat::expect_equal(.Random.seed, seed_before)
  testthat::expect_equal(graphics::par(c("mfrow", "usr")), graphics_before)
  testthat::expect_error(
    plot(diagnostics, pc = "PC1", comparison = "within_anterior", ci_level = 1),
    "ci_level"
  )
  testthat::expect_error(
    plot(diagnostics, pc = "PC1", comparison = "within_anterior", seed = NA_real_),
    "seed"
  )
  testthat::expect_error(
    plot(diagnostics, pc = c("PC1", "PC1"), comparison = c("within_anterior", "anterior_vs_posterior", "within_posterior")),
    "same length"
  )
  auto_named <- regime_module_diagnostics(
    pca,
    modules = modules,
    comparisons = list(c("anterior", "posterior")),
    pcs = "PC1"
  )
  testthat::expect_named(auto_named$comparisons, "anterior_vs_posterior")
  testthat::expect_true(is.na(
    .regime_module_comparison_score(
      matrix(1, nrow = 1, ncol = 1, dimnames = list("c", "c")),
      "c",
      "c"
    )
  ))
  testthat::expect_s3_class(
    plot(
      diagnostics,
      pc = "PC1",
      comparison = c("within_anterior", "within_posterior"),
      n_boot = 0
    ),
    "regime_module_diagnostics"
  )
  testthat::expect_s3_class(
    plot(
      diagnostics,
      pc = c("PC1", "PC1"),
      comparison = "within_anterior",
      xlab = c("PC1 custom", "PC1 repeat"),
      ylab = "Module custom",
      n_boot = 0
    ),
    "regime_module_diagnostics"
  )
  testthat::expect_s3_class(
    plot(
      diagnostics,
      pc = "PC1",
      comparison = "within_anterior",
      n_boot = 1,
      seed = NULL
    ),
    "regime_module_diagnostics"
  )
  set.seed(2027)
  endpoint_seed_before <- .Random.seed
  testthat::expect_warning(
    endpoint_plot <- plot(
      diagnostics,
      pc = "PC1",
      comparison = "within_anterior",
      n_boot = 1,
      seed = .Machine$integer.max
    ),
    NA
  )
  testthat::expect_s3_class(endpoint_plot, "regime_module_diagnostics")
  testthat::expect_equal(.Random.seed, endpoint_seed_before)
  testthat::expect_s3_class(
    plot(
      diagnostics,
      pc = c("PC1", "PC1"),
      comparison = c("within_anterior", "within_posterior"),
      n_boot = 0,
      seed = .Machine$integer.max
    ),
    "regime_module_diagnostics"
  )
})

test_that("summarize_regime_covariances can derive tip counts and ages from search", {
  compact <- .avian_skeleton_compact_posthoc()
  search <- .avian_skeleton_compact_search()

  out <- summarize_regime_covariances(
    compact$correlation_matrices,
    search = search
  )

  testthat::expect_false(all(is.na(out$tip_count)))
  testthat::expect_false(all(is.na(out$regime_age)))
  testthat::expect_true(all(out$tip_count[!is.na(out$tip_count)] > 10))
})

test_that("regime phylogeny collapse rejects duplicated output tip labels", {
  testthat::skip_if_not_installed("ape")
  testthat::skip_if_not_installed("phytools")

  tree <- ape::read.tree(text = "(((a:1,b:1):1,(c:1,d:1):1):1,e:3);")
  tree <- phytools::paintSubTree(
    tree,
    node = ape::Ntip(tree) + 1L,
    state = "root"
  )
  tree <- phytools::paintSubTree(
    tree,
    node = ape::Ntip(tree) + 2L,
    state = "wing"
  )
  tree <- phytools::paintSubTree(
    tree,
    node = ape::Ntip(tree) + 3L,
    state = "c"
  )

  testthat::expect_error(
    .collapse_regime_phylogeny(tree),
    "collapsed tree tip labels.*duplicated identifier.*c"
  )
})

test_that("regime_integration_pgls reproduces representative manuscript model", {
  testthat::skip_if_not_installed("phylolm")
  compact <- .avian_skeleton_compact_posthoc()
  sensitivity <- .avian_skeleton_compact_sensitivity()

  collapse_warnings <- character()
  fit <- withCallingHandlers(
    regime_integration_pgls(
      compact$vars_cors$min20.ic20.gic,
      search = sensitivity$min20.ic20.gic
    ),
    warning = function(w) {
      collapse_warnings <<- c(collapse_warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )

  pgls_coef <- stats::setNames(
    compact$pgls_summary$coefficients$estimate,
    compact$pgls_summary$coefficients$term
  )
  testthat::expect_s3_class(fit, "phylolm")
  testthat::expect_equal(fit$r.squared, compact$pgls_summary$r_squared, tolerance = 1e-6)
  testthat::expect_equal(unname(fit$coefficients["scale(log(vars))"]), unname(pgls_coef["scale(log(vars))"]), tolerance = 1e-6)
  testthat::expect_equal(
    unname(fit$coefficients["scale(fisher_z_transform(corrs))"]),
    unname(pgls_coef["scale(fisher_z_transform(corrs))"]),
    tolerance = 1e-6
  )
  testthat::expect_length(collapse_warnings, 1L)
  testthat::expect_match(collapse_warnings, "nonmonophyletic regime")
  testthat::expect_match(collapse_warnings, "186")
  testthat::expect_match(collapse_warnings, "52")
})

test_that("regime_integration_relationships mirrors manuscript residual filtering and bootstrap curves", {
  compact <- .avian_skeleton_compact_posthoc()

  relationships <- regime_integration_relationships(
    compact$vars_cors,
    n_boot = 25,
    seed = 1
  )

  testthat::expect_s3_class(relationships, "regime_integration_relationships")
  testthat::expect_named(relationships, c(
    "combined",
    "variance_points",
    "correlation_points",
    "variance_removed",
    "correlation_removed",
    "variance_curve",
    "correlation_curve",
    "variance_lm",
    "correlation_lm",
    "settings"
  ))
  testthat::expect_equal(nrow(relationships$combined), sum(vapply(compact$vars_cors, nrow, integer(1))))
  testthat::expect_true(nrow(relationships$variance_points) < nrow(relationships$combined))
  testthat::expect_true(nrow(relationships$correlation_points) < nrow(relationships$combined))
  testthat::expect_equal(nrow(relationships$variance_curve), 100)
  testthat::expect_equal(nrow(relationships$correlation_curve), 100)
  testthat::expect_s3_class(relationships$variance_lm, "lm")
  testthat::expect_s3_class(relationships$correlation_lm, "lm")

  combined <- as.data.frame(relationships, component = "combined")
  models <- as.data.frame(relationships, component = "models")
  testthat::expect_equal(nrow(combined), nrow(relationships$combined))
  testthat::expect_true(all(c("panel", "slope", "r_squared") %in% names(models)))
})

test_that("regime_integration_relationships can filter by user-supplied minimum tips", {
  summary <- .regime_integration_relationship_summary()
  summary$tip_count <- c(2, 3, 4, 5, 6, 7)

  relationships <- regime_integration_relationships(
    summary,
    min_tips = 4,
    resid_sd_threshold_vars = 1000,
    resid_sd_threshold_corrs = 1000,
    n_boot = 5,
    seed = 1
  )

  testthat::expect_equal(relationships$combined$regime, paste0("r", 3:6))
  testthat::expect_true(all(relationships$combined$tip_count >= 4))

  testthat::expect_error(
    regime_integration_relationships(
      summary[, setdiff(names(summary), "tip_count")],
      min_tips = 4,
      n_boot = 5
    ),
    "tip_count"
  )
})

test_that("plot methods draw regime integration panels and PCA loadings", {
  mats <- .regime_integration_pca_mats()
  pca <- regime_correlation_pca(mats)

  relationships <- regime_integration_relationships(
    .regime_integration_relationship_summary(),
    resid_sd_threshold_vars = 1000,
    resid_sd_threshold_corrs = 1000,
    n_boot = 5,
    seed = 1
  )

  pdf_file <- tempfile(fileext = ".pdf")
  grDevices::pdf(pdf_file)
  on.exit({
    grDevices::dev.off()
    unlink(pdf_file)
  }, add = TRUE)

  testthat::expect_s3_class(plot(relationships, panel = "variance"), "regime_integration_relationships")
  testthat::expect_s3_class(plot(pca, type = "loadings", components = 1:2), "regime_correlation_pca")
  testthat::expect_s3_class(
    plot(
      pca,
      type = "loadings",
      components = 1:2,
      cluster = "local",
      heatmap_engine = "base"
    ),
    "regime_correlation_pca"
  )
})

test_that("clustered PCA loading heatmaps request square bodies", {
  mat <- matrix(0, nrow = 12, ncol = 12)
  size <- .regime_heatmap_body_size(mat)

  testthat::expect_equal(
    grid::convertWidth(size$width, "mm", valueOnly = TRUE),
    grid::convertHeight(size$height, "mm", valueOnly = TRUE)
  )
})

test_that("regime_integration_pgls requires tip counts when downstream min_tips is requested", {
  testthat::skip_if_not_installed("ape")
  testthat::skip_if_not_installed("phytools")
  testthat::skip_if_not_installed("phylolm")

  summary <- data.frame(
    regime = c("hind", "root", "wing"),
    rate = c(1.2, 0.8, 1.7),
    mean_variance = c(0.4, 0.6, 1.1),
    mean_abs_correlation = c(0.2, 0.3, 0.4)
  )

  testthat::expect_error(
    regime_integration_pgls(
      summary,
      tree = .regime_integration_success_tree(),
      min_tips = 3
    ),
    "tip_count"
  )
})

test_that("regime_integration_pgls requires at least three collapsed regimes", {
  testthat::skip_if_not_installed("ape")
  testthat::skip_if_not_installed("phytools")
  testthat::skip_if_not_installed("phylolm")

  two_regimes <- data.frame(
    regime = c("hind", "wing"),
    rate = c(1.2, 1.7),
    mean_variance = c(0.4, 1.1),
    mean_abs_correlation = c(0.2, 0.4)
  )

  testthat::expect_error(
    regime_integration_pgls(
      two_regimes,
      tree = .regime_integration_success_tree()
    ),
    "At least three"
  )
})

test_that("regime integration plotting and data-frame methods cover alternate panels", {
  summary <- .regime_integration_relationship_summary()
  relationships <- regime_integration_relationships(
    summary,
    resid_sd_threshold_vars = 1000,
    resid_sd_threshold_corrs = 1000,
    n_boot = 1
  )

  testthat::expect_equal(as.data.frame(relationships, component = "variance_curve"), relationships$variance_curve)
  testthat::expect_equal(as.data.frame(relationships, component = "correlation_curve"), relationships$correlation_curve)
  testthat::expect_error(
    as.data.frame(relationships, component = "combined", unused = TRUE),
    "Unused"
  )

  pdf_file <- tempfile(fileext = ".pdf")
  grDevices::pdf(pdf_file)
  on.exit({
    grDevices::dev.off()
    unlink(pdf_file)
  }, add = TRUE)
  testthat::expect_s3_class(plot(relationships, panel = "both"), "regime_integration_relationships")
  testthat::expect_s3_class(plot(relationships, panel = "correlation"), "regime_integration_relationships")
  testthat::expect_error(plot(relationships, panel = "both", main = character()), "one value per requested panel")
})

test_that("relationship summaries accept manuscript-shaped tables and preserve bootstrap seed state", {
  vars_cors <- data.frame(
    State = paste0("r", 1:6),
    rate = c(1.0, 1.4, 1.7, 2.5, 3.1, 4.2),
    vars = c(0.7, 0.9, 1.4, 1.8, 2.6, 3.4),
    corrs = c(0.12, 0.18, 0.23, 0.31, 0.37, 0.44),
    run = rep("manual", 6),
    tip_count = c(2, 3, 4, 5, 6, 7)
  )
  relationships <- regime_integration_relationships(
    vars_cors,
    min_tips = 3,
    resid_sd_threshold_vars = 1000,
    resid_sd_threshold_corrs = 1000,
    n_boot = 1
  )
  testthat::expect_equal(unique(relationships$combined$run), "manual")
  testthat::expect_true(all(relationships$combined$tip_count >= 3))

  unnamed_runs <- regime_integration_relationships(
    list(vars_cors, vars_cors),
    resid_sd_threshold_vars = 1000,
    resid_sd_threshold_corrs = 1000,
    n_boot = 1
  )
  testthat::expect_equal(unique(unnamed_runs$combined$run), c("run1", "run2"))

  mixed_runs <- list(empirical = vars_cors, vars_cors)
  mixed_runs_relationships <- regime_integration_relationships(
    mixed_runs,
    resid_sd_threshold_vars = 1000,
    resid_sd_threshold_corrs = 1000,
    n_boot = 1
  )
  testthat::expect_equal(
    unique(mixed_runs_relationships$combined$run),
    c("empirical", "run2")
  )

  missing_name_runs <- list(empirical = vars_cors, vars_cors)
  names(missing_name_runs)[[2L]] <- NA_character_
  missing_name_relationships <- regime_integration_relationships(
    missing_name_runs,
    resid_sd_threshold_vars = 1000,
    resid_sd_threshold_corrs = 1000,
    n_boot = 1
  )
  testthat::expect_equal(
    unique(missing_name_relationships$combined$run),
    c("empirical", "run2")
  )

  whitespace_name_runs <- list(vars_cors, empirical = vars_cors)
  names(whitespace_name_runs)[[1L]] <- "   "
  whitespace_name_relationships <- regime_integration_relationships(
    whitespace_name_runs,
    resid_sd_threshold_vars = 1000,
    resid_sd_threshold_corrs = 1000,
    n_boot = 1
  )
  testthat::expect_equal(
    unique(whitespace_name_relationships$combined$run),
    c("run1", "empirical")
  )

  generated_collision_runs <- list(run2 = vars_cors, vars_cors)
  names(generated_collision_runs)[[2L]] <- ""
  testthat::expect_error(
    regime_integration_relationships(
      generated_collision_runs,
      resid_sd_threshold_vars = 1000,
      resid_sd_threshold_corrs = 1000,
      n_boot = 1
    ),
    "summaries.*duplicated run name.*run2"
  )

  old_seed <- if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
    get(".Random.seed", envir = .GlobalEnv)
  } else {
    NULL
  }
  if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
    rm(".Random.seed", envir = .GlobalEnv)
  }
  on.exit({
    if (is.null(old_seed)) {
      if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
        rm(".Random.seed", envir = .GlobalEnv)
      }
    } else {
      assign(".Random.seed", old_seed, envir = .GlobalEnv)
    }
  }, add = TRUE)
  curve <- .regime_bootstrap_curve(
    data = data.frame(x = 1:4, y = c(1, 1.5, 2.5, 4)),
    x_col = "x",
    y_col = "y",
    formula = y ~ x,
    n_boot = 1,
    ci_level = 0.95,
    seed = 1
  )
  testthat::expect_equal(nrow(curve), 100)
  testthat::expect_false(exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
  testthat::expect_error(
    .regime_bootstrap_curve(
      data = data.frame(x = 1:2, y = 1:2),
      x_col = "x",
      y_col = "y",
      formula = y ~ x,
      n_boot = 1,
      ci_level = 0.95,
      seed = NULL
    ),
    "At least three rows"
  )
  testthat::expect_equal(
    .regime_integration_plot_main(c("A", "B", "C"), c("x", "y")),
    as.list(c("A", "B"))
  )
})

test_that("compact avian object reproduces pGLS and PCA diagnostics", {
  compact <- .avian_skeleton_compact_posthoc()
  testthat::expect_identical(attr(compact, "provenance"), compact$provenance)
  testthat::expect_identical(
    compact$provenance$source_repo_path,
    "jakeberv/passerine-bodyplan-evolution"
  )
  testthat::expect_identical(
    compact$provenance$manuscript_script,
    "TemporalAnalyses.R"
  )

  pgls_coef <- stats::setNames(
    compact$pgls_summary$coefficients$estimate,
    compact$pgls_summary$coefficients$term
  )
  testthat::expect_equal(compact$pgls_summary$r_squared, 0.8411565, tolerance = 1e-6)
  testthat::expect_equal(unname(pgls_coef["scale(log(vars))"]), 0.5384966, tolerance = 1e-6)
  testthat::expect_equal(
    unname(pgls_coef["scale(fisher_z_transform(corrs))"]),
    -0.4910791,
    tolerance = 1e-6
  )

  pca <- regime_correlation_pca(
    compact$correlation_matrices,
    tip_counts = compact$pca_inputs$tip_counts,
    regime_ages = compact$pca_inputs$regime_ages,
    trait_labels = compact$trait_labels
  )
  testthat::expect_equal(round(100 * pca$variance_explained[1:4], 2), c(39.89, 11.96, 7.13, 6.52))

  diag_vals <- stats::setNames(
    compact$module_diagnostics$correlations$value,
    compact$module_diagnostics$correlations$diagnostic
  )
  testthat::expect_equal(unname(diag_vals["cor(PC1, within-wing mean)"]), 0.878513, tolerance = 1e-6)
  testthat::expect_equal(
    unname(diag_vals["cor(PC1, hindlimb+cranial vs wing mean)"]),
    -0.9341625,
    tolerance = 1e-6
  )
  testthat::expect_equal(
    unname(diag_vals["cor(PC2, cranial vs hindlimb mean)"]),
    0.6833727,
    tolerance = 1e-6
  )
  testthat::expect_equal(
    unname(diag_vals["cor(PC2, cranial vs wing mean)"]),
    -0.566403,
    tolerance = 1e-6
  )

  testthat::expect_equal(compact$settings$pca_min_tips, 10)
  testthat::expect_equal(compact$settings$pca_tip_filter, "tip_count > pca_min_tips")
  testthat::expect_equal(compact$settings$fit_min_tips, 2)
  testthat::expect_true(compact$settings$remove_high_corr)
  testthat::expect_equal(compact$settings$corr_threshold, 0.95)
})

test_that("Part 5 guardrail uses post-hoc matrices instead of proportional search VCVs", {
  compact <- .avian_skeleton_compact_posthoc()

  testthat::expect_true(any(grepl("not search\\$VCVs", compact$provenance$notes)))

  vignette_path <- testthat::test_path("../../vignettes/avian-skeleton-part-5.Rmd")
  testthat::skip_if_not(file.exists(vignette_path), "Part 5 vignette not present")
  vignette_text <- paste(readLines(vignette_path, warn = FALSE), collapse = "\n")

  testthat::expect_match(vignette_text, "correlation_matrices", fixed = TRUE)
  testthat::expect_match(vignette_text, "fit_regime_covariance_runs", fixed = TRUE)
  testthat::expect_match(vignette_text, "summarize_regime_covariance_runs", fixed = TRUE)
  testthat::expect_false(grepl("posthoc_fits_by_run <- lapply", vignette_text, fixed = TRUE))
  testthat::expect_false(grepl("posthoc_summaries <- Map", vignette_text, fixed = TRUE))
  testthat::expect_false(grepl("regime_correlation_pca\\(bodyplan_search\\$VCVs", vignette_text))
  testthat::expect_false(grepl("regime_correlation_pca\\([^\\n]*\\$VCVs", vignette_text))
})

test_that("regime tree resolution rejects unnamed mapped state segments", {
  testthat::skip_if_not_installed("ape")
  testthat::skip_if_not_installed("phytools")

  tree <- .regime_integration_tree()
  names(tree$maps[[1L]]) <- NULL

  testthat::expect_error(
    .regime_resolve_tree(tree = tree),
    "mapped regime-state identifiers.*non-empty"
  )
})

test_that("regime covariance containers validate required structure", {
  covariance <- diag(2)
  dimnames(covariance) <- list(c("a", "b"), c("a", "b"))

  empty <- structure(
    list(covariances = list(), status = data.frame(regime = character())),
    class = c("regime_covariances", "list")
  )
  testthat::expect_error(
    .regime_extract_covariance_list(empty),
    "x\\$covariances.*non-empty named list"
  )

  missing_status_id <- structure(
    list(covariances = list(r1 = covariance), status = data.frame(status = "ok")),
    class = c("regime_covariances", "list")
  )
  testthat::expect_error(
    .regime_extract_covariance_list(missing_status_id),
    "x\\$status.*regime"
  )
})

test_that("regime covariance helpers reject partial trait names and degenerate templates", {
  partial_names <- diag(2)
  rownames(partial_names) <- c("a", "b")
  testthat::expect_error(
    .regime_validate_covariance_matrix(partial_names, "r1"),
    "matching row and column trait names"
  )

  mismatched_names <- diag(2)
  rownames(mismatched_names) <- c("a", "b")
  colnames(mismatched_names) <- c("b", "a")
  testthat::expect_error(
    .regime_align_matrix_traits(mismatched_names, c("a", "b"), "r2"),
    "matching row and column trait names"
  )

  zero <- matrix(0, nrow = 2, ncol = 2)
  testthat::expect_false(
    .regime_covariance_list_is_scalar_proportional(list(zero, zero))
  )
})

test_that("regime optional-package guards report missing dependencies and retain base fallbacks", {
  pgls <- regime_integration_pgls
  environment(pgls) <- list2env(
    list(requireNamespace = function(...) FALSE),
    parent = environment(regime_integration_pgls)
  )
  testthat::expect_error(
    pgls(data.frame()),
    "Package `phylolm` is required"
  )

  palette <- .regime_pca_loading_palette
  environment(palette) <- list2env(
    list(requireNamespace = function(...) FALSE),
    parent = environment(.regime_pca_loading_palette)
  )
  testthat::expect_length(palette(7L), 7L)

  complex_draw <- .regime_pca_plot_loadings_complex_heatmap
  environment(complex_draw) <- list2env(
    list(requireNamespace = function(...) FALSE),
    parent = environment(.regime_pca_plot_loadings_complex_heatmap)
  )
  testthat::expect_false(complex_draw(
    x = NULL,
    component_idx = 1L,
    main = NULL,
    palette = NULL,
    cluster = "none",
    show_dendrogram = FALSE,
    show_legend = FALSE,
    cex.axis = 1
  ))

  pca <- regime_correlation_pca(.regime_integration_pca_mats())
  plot_pca <- plot.regime_correlation_pca
  environment(plot_pca) <- list2env(
    list(.regime_pca_plot_loadings_complex_heatmap = function(...) FALSE),
    parent = environment(plot.regime_correlation_pca)
  )
  testthat::expect_error(
    plot_pca(
      pca,
      type = "loadings",
      components = 1L,
      heatmap_engine = "ComplexHeatmap"
    ),
    "requires packages.*ComplexHeatmap.*circlize"
  )
})
