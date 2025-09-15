# tests/testthat/test-searchOptimalConfiguration.R

testthat::skip_on_cran()

# ---- dependency guard --------------------------------------------------------
skip_if_missing_deps <- function() {
  testthat::skip_if_not_installed("ape")
  testthat::skip_if_not_installed("phytools")
  testthat::skip_if_not_installed("mvMORPH")
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

# ---- the test ---------------------------------------------------------------
test_that("searchOptimalConfiguration runs end-to-end on simulated data", {
  skip_if_missing_deps()

  simdata <- load_simdata_fixture()

  # Coerce paintedTree to phylo if needed
  tree0 <-
    if (inherits(simdata$paintedTree, "phylo")) simdata$paintedTree
  else ape::as.phylo(simdata$paintedTree)

  # tree0 is your unpainted phylo object
  ntips <- ape::Ntip(tree0)

  set.seed(123)  # reproducible sampling
  tips_to_keep <- sample(tree0$tip.label, size = 100)

  # Drop all other tips
  tree0 <- ape::drop.tip(tree0, setdiff(tree0$tip.label, tips_to_keep))

  # Paint a global baseline state "0" from the root (stem = FALSE to avoid splitting)
  root <- ape::Ntip(tree0) + 1L
  baseline <- phytools::paintSubTree(
    tree      = tree0,
    node      = root,
    state     = "0",
    anc.state = "0",
    stem      = FALSE
  )

  # Trait matrix/data.frame; ensure row order matches tip.labels
  X <- simdata$simulatedData
  testthat::expect_true(is.matrix(X) || is.data.frame(X))
  testthat::expect_true(all(baseline$tip.label %in% rownames(X)))
  X <- X[baseline$tip.label, , drop = FALSE]

  # Run with modest resources; plotting off for tests
  set.seed(123)
  res <- searchOptimalConfiguration(
    baseline_tree              = baseline,
    trait_data                 = (X),
    formula                    = "trait_data ~ 1",
    min_descendant_tips        = 10,
    num_cores                  = 1,
    shift_acceptance_threshold = 5,
    plot                       = FALSE,
    IC                         = "GIC",
    store_model_fit_history    = FALSE,
    method                     = "LL"
  )

  # Core structure checks
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
  testthat::expect_true(is.finite(res$baseline_ic))
  testthat::expect_true(is.finite(res$optimal_ic))
  testthat::expect_true(res$IC_used %in% c("GIC", "BIC"))
  testthat::expect_true(inherits(res$tree_no_uncertainty_untransformed, "phylo"))
  testthat::expect_true(inherits(res$tree_no_uncertainty_transformed, "phylo"))
  testthat::expect_true(is.list(res$VCVs) || is.null(res$VCVs))

  # Best model should not be worse than baseline
  testthat::expect_lte(res$optimal_ic, res$baseline_ic)
})
