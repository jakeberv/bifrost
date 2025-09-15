# tests/testthat/test-searchOptimalConfiguration.R

testthat::skip_on_cran()

# ---- dependency guard --------------------------------------------------------
skip_if_missing_deps <- function() {
  testthat::skip_if_not_installed("ape")
  testthat::skip_if_not_installed("phytools")
  testthat::skip_if_not_installed("mvMORPH")
  testthat::skip_if_not_installed("future")
}

# ---- locate and load fixture -------------------------------------------------
load_simdata_fixture <- function() {
  # Expect the file at tests/testthat/fixtures/simdata.RDS
  rds_path <- testthat::test_path("fixtures", "simdata.RDS")
  if (!file.exists(rds_path)) {
    testthat::skip(paste("Fixture not found:", rds_path))
  }
  obj <- readRDS(rds_path)

  # Some fixtures are a list-of-lists; unwrap first element if needed
  bundle <- if (is.list(obj) && !all(c("paintedTree", "simulatedData") %in% names(obj))) obj[[1]] else obj

  # Basic sanity
  if (!all(c("paintedTree", "simulatedData") %in% names(bundle))) {
    testthat::skip("Fixture does not contain expected names: paintedTree, simulatedData")
  }
  bundle
}

# ---- helper to build 100-tip baseline + aligned data -------------------------
build_baseline_and_data <- function(simdata) {
  # Coerce paintedTree to phylo if needed
  tree0 <-
    if (inherits(simdata$paintedTree, "phylo")) simdata$paintedTree
  else ape::as.phylo(simdata$paintedTree)

  # Subsample 100 tips (for speed and determinism)
  set.seed(123)
  keep <- sample(tree0$tip.label, size = 100)
  tree0 <- ape::drop.tip(tree0, setdiff(tree0$tip.label, keep))

  # Paint a global baseline state "0" from the root
  root <- ape::Ntip(tree0) + 1L
  baseline <- phytools::paintSubTree(
    tree      = tree0,
    node      = root,
    state     = "0",
    anc.state = "0",
    stem      = FALSE
  )

  # Align trait data order to tree tips
  X <- simdata$simulatedData
  testthat::expect_true(is.matrix(X) || is.data.frame(X))
  testthat::expect_true(all(baseline$tip.label %in% rownames(X)))
  X <- X[baseline$tip.label, , drop = FALSE]

  list(tree = baseline, X = X)
}

# ---- tiny helpers for tolerant assertions ------------------------------------
expect_phylo_or_null <- function(x) {
  testthat::expect_true(is.null(x) || inherits(x, "phylo"))
}
expect_numeric_scalar <- function(x) {
  testthat::expect_true(is.numeric(x) && length(x) == 1L && is.finite(x))
}

# ---- Test 1: original end-to-end with parallel-capable path (num_cores = 1) --
test_that("searchOptimalConfiguration runs end-to-end on simulated data", {
  skip_if_missing_deps()
  simdata <- load_simdata_fixture()
  built <- build_baseline_and_data(simdata)
  baseline <- built$tree
  X <- built$X

  set.seed(123)
  res <- searchOptimalConfiguration(
    baseline_tree              = baseline,
    trait_data                 = X,
    formula                    = "trait_data ~ 1",
    min_descendant_tips        = 10,
    num_cores                  = 1,
    shift_acceptance_threshold = 5,
    plot                       = FALSE,
    IC                         = "GIC",
    store_model_fit_history    = FALSE,
    method                     = "LL"
  )

  # Core structure checks (present names)
  testthat::expect_type(res, "list")
  testthat::expect_true(all(c(
    "tree_no_uncertainty_transformed",
    "tree_no_uncertainty_untransformed",
    "model_no_uncertainty",
    "shift_nodes_no_uncertainty",
    "optimal_ic",
    "baseline_ic",
    "IC_used",
    "num_candidates"
  ) %in% names(res)))

  # Types/values
  expect_numeric_scalar(res$baseline_ic)
  expect_numeric_scalar(res$optimal_ic)
  testthat::expect_true(res$IC_used %in% c("GIC", "BIC"))
  # These may be NULL if no shifts are accepted; otherwise phylo
  expect_phylo_or_null(res$tree_no_uncertainty_untransformed)
  expect_phylo_or_null(res$tree_no_uncertainty_transformed)
  testthat::expect_true(is.list(res$VCVs) || is.null(res$VCVs))

  # If shifts were accepted, best should improve or match baseline
  if (length(res$shift_nodes_no_uncertainty) > 0) {
    testthat::expect_lte(res$optimal_ic, res$baseline_ic)
  } else {
    # No shifts: optimal equals baseline (within tiny tolerance)
    testthat::expect_equal(res$optimal_ic, res$baseline_ic, tolerance = 1e-8)
  }
})

# ---- Test 2: explicitly non-parallel (force sequential plan) -----------------
test_that("searchOptimalConfiguration also runs in purely sequential mode", {
  skip_if_missing_deps()
  simdata <- load_simdata_fixture()
  built <- build_baseline_and_data(simdata)
  baseline <- built$tree
  X <- built$X

  # Force a sequential future plan for the duration of this test
  old_plan <- NULL
  if (requireNamespace("future", quietly = TRUE)) {
    old_plan <- future::plan()
    future::plan(future::sequential)
  }
  on.exit({
    if (!is.null(old_plan)) future::plan(old_plan)
  }, add = TRUE)

  set.seed(456)
  res <- searchOptimalConfiguration(
    baseline_tree              = baseline,
    trait_data                 = X,
    formula                    = "trait_data ~ 1",
    min_descendant_tips        = 10,
    num_cores                  = 1,  # hint to not parallelize
    shift_acceptance_threshold = 5,
    plot                       = FALSE,
    IC                         = "GIC",
    store_model_fit_history    = FALSE,
    method                     = "LL"
  )

  # Minimal sanity checks
  testthat::expect_type(res, "list")
  expect_numeric_scalar(res$baseline_ic)
  expect_numeric_scalar(res$optimal_ic)
  expect_phylo_or_null(res$tree_no_uncertainty_untransformed)
})

# ---- Test 3: forced no-shift scenario (very high acceptance threshold) -------
test_that("searchOptimalConfiguration returns sensible output when no shifts are accepted", {
  skip_if_missing_deps()
  simdata <- load_simdata_fixture()
  built <- build_baseline_and_data(simdata)
  baseline <- built$tree
  X <- built$X

  set.seed(789)
  res <- searchOptimalConfiguration(
    baseline_tree              = baseline,
    trait_data                 = X,
    formula                    = "trait_data ~ 1",
    min_descendant_tips        = 10,
    num_cores                  = 1,
    shift_acceptance_threshold = 1e9,  # effectively forbids accepting any shift
    plot                       = FALSE,
    IC                         = "GIC",
    store_model_fit_history    = TRUE,
    method                     = "LL"
  )

  # No shifts detected
  testthat::expect_equal(length(res$shift_nodes_no_uncertainty), 0L)

  # Optimal IC equals baseline IC
  testthat::expect_equal(res$optimal_ic, res$baseline_ic, tolerance = 1e-8)

  # Trees for "no-uncertainty" path may be NULL (since no shift model exists); allow NULL or phylo
  expect_phylo_or_null(res$tree_no_uncertainty_untransformed)
  expect_phylo_or_null(res$tree_no_uncertainty_transformed)

  # Model may be NULL in this path
  testthat::expect_true(is.null(res$model_no_uncertainty) || is.list(res$model_no_uncertainty))

  # History exists and contains only rejected entries — be robust to type (logical/num/char/factor)
  if (!is.null(res$model_fit_history) && is.list(res$model_fit_history)) {
    hist <- res$model_fit_history
    if (!is.null(hist$ic_acceptance_matrix)) {
      acc <- hist$ic_acceptance_matrix

      testthat::expect_true(is.matrix(acc) || is.data.frame(acc))

      if (NCOL(acc) >= 2) {
        vals <- acc[, 2]

        # flatten to a simple vector
        vals <- if (is.data.frame(vals)) unlist(vals, use.names = FALSE) else vals
        vals <- if (is.list(vals)) unlist(vals, use.names = FALSE) else vals
        if (is.factor(vals)) vals <- as.character(vals)

        # Coerce robustly and ensure no TRUE
        is_true <- rep(FALSE, length(vals))
        if (is.logical(vals)) {
          is_true <- vals
        } else if (is.numeric(vals)) {
          is_true <- vals != 0
        } else if (is.character(vals)) {
          v <- trimws(tolower(vals))
          is_true <- v %in% c("true", "t", "1")
        }
        testthat::expect_false(any(is_true, na.rm = TRUE))
      }
    }
  }
})
