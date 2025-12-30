# tests/testthat/test-fitMvglsAndExtractBIC-formula.R

# testthat::local_edition(3)

# Dependency guard
skip_if_missing_deps <- function() {
  testthat::skip_if_not_installed("ape")
  testthat::skip_if_not_installed("phytools")
  testthat::skip_if_not_installed("mvMORPH")
}

# Inline helpers (scoped to this test file)

# Build a simple dataset with 1 predictor and p responses; rownames=tip labels
make_data_for_tree <- function(tree, p = 2, beta = 0.5, sd = 0.5, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  n <- ape::Ntip(tree)
  x <- stats::rnorm(n)
  Y <- replicate(p, beta * x + stats::rnorm(n, sd = sd))
  df <- data.frame(x = x, Y)
  colnames(df) <- c("x", paste0("y", seq_len(p)))
  rownames(df) <- tree$tip.label
  df
}

# Initialize a SIMMAP by painting every *edge* with global baseline "0"
# (paint child node of each edge; avoids painting the root which has no incoming edge)
init_global_map <- function(tr, global_state = "0") {
  child_nodes <- unique(tr$edge[, 2])
  phytools::paintBranches(
    tr,
    edge      = child_nodes,
    state     = global_state,
    anc.state = global_state
  )
}

# Single-regime SIMMAP: entire tree in state "0"
make_single_regime_tree <- function(n_tip = 30, global_state = "0", seed = 101, scale = 1) {
  set.seed(seed)
  tr  <- phytools::pbtree(n = n_tip, scale = scale)
  sim <- init_global_map(tr, global_state = global_state)
  sim
}

# Two-regime SIMMAP: global "0" + nested "1" on a deterministic internal clade
make_two_regime_tree <- function(n_tip = 30, seed = 202, scale = 1) {
  set.seed(seed)
  tr  <- phytools::pbtree(n = n_tip, scale = scale)
  sim <- init_global_map(tr, global_state = "0")

  internal_nodes <- (ape::Ntip(tr) + 1L):(ape::Ntip(tr) + ape::Nnode(tr))
  target_node <- internal_nodes[ceiling(length(internal_nodes) / 2)]

  sim <- phytools::paintSubTree(sim, node = target_node, state = "1", stem = TRUE)
  sim
}

# Extract a scalar BIC from either a numeric or a list-like object
# (defensive, though stats::BIC usually returns numeric)
get_bic_scalar <- function(bic_obj) {
  if (is.list(bic_obj) && !is.null(bic_obj$BIC)) {
    return(as.numeric(bic_obj$BIC))
  }
  as.numeric(bic_obj)
}

# Group: core fit behavior (mirrors GIC suite)

# Test: Single-regime tree fits and returns a scalar BIC (checks return value)
test_that("Single-regime tree fits and returns a scalar BIC", {
  skip_if_missing_deps()

  painted_tree <- make_single_regime_tree(20)
  set.seed(123)
  dat <- make_data_for_tree(painted_tree, p = 2)

  form_chr <- "cbind(y1, y2) ~ x"

  res <- suppressWarnings(
    fitMvglsAndExtractBIC.formula(
      formula      = form_chr,
      painted_tree = painted_tree,
      trait_data   = dat,
      data         = dat,
      method       = "LL"     # BIC is well-defined for LL fits
    )
  )

  expect_type(res, "list")
  expect_true("model" %in% names(res))
  expect_true("BIC"   %in% names(res))
  expect_true(inherits(res$model, "mvgls"))

  bic <- get_bic_scalar(res$BIC)
  expect_true(is.finite(bic))
  expect_length(bic, 1)
})

# Test: Multi-regime tree fits and returns a scalar BIC (BMM path) (checks return value)
test_that("Multi-regime tree fits and returns a scalar BIC (BMM path)", {
  skip_if_missing_deps()

  painted_tree <- make_two_regime_tree(25)
  set.seed(456)
  dat <- make_data_for_tree(painted_tree, p = 3)

  form_chr <- "cbind(y1, y2, y3) ~ x"

  res <- suppressWarnings(
    fitMvglsAndExtractBIC.formula(
      formula      = form_chr,
      painted_tree = painted_tree,
      trait_data   = dat,
      data         = dat,
      method       = "LL"
    )
  )

  expect_true(inherits(res$model, "mvgls"))

  bic <- get_bic_scalar(res$BIC)
  expect_true(is.finite(bic))
  expect_length(bic, 1)
})

# Group: input validation errors
# Test: Row names mismatch triggers informative error (expects error)
test_that("Row names mismatch triggers informative error", {
  skip_if_missing_deps()

  painted_tree <- make_two_regime_tree(15)
  set.seed(789)
  dat <- make_data_for_tree(painted_tree, p = 2)

  bad <- dat
  rownames(bad) <- paste0(rownames(bad), "_X")

  form_chr <- "cbind(y1, y2) ~ x"

  expect_error(
    fitMvglsAndExtractBIC.formula(
      formula      = form_chr,
      painted_tree = painted_tree,
      trait_data   = bad,
      data         = bad
    ),
    "Row names of trait_data must exactly match the tip labels of the tree."
  )
})

# Test: Formula must be provided as a character string
test_that("Formula must be provided as a character string", {
  skip_if_missing_deps()

  painted_tree <- make_single_regime_tree(10)
  set.seed(321)
  dat <- make_data_for_tree(painted_tree, p = 2)

  form_obj <- cbind(y1, y2) ~ x

  expect_error(
    fitMvglsAndExtractBIC.formula(
      formula      = form_obj,
      painted_tree = painted_tree,
      trait_data   = dat,
      data         = dat
    ),
    "A character formula must be provided."
  )

  expect_error(
    fitMvglsAndExtractBIC.formula(
      painted_tree = painted_tree,
      trait_data   = dat,
      data         = dat
    ),
    "A character formula must be provided."
  )
})

# Group: determinism, branch selection, and option passthrough

# Test: Determinism: same inputs yield identical BIC (re-runs with same inputs)
test_that("Determinism: same inputs yield identical BIC", {
  skip_if_missing_deps()
  painted_tree <- make_two_regime_tree(25)
  set.seed(999)
  dat <- make_data_for_tree(painted_tree, p = 2)
  form_chr <- "cbind(y1, y2) ~ x"

  r1 <- suppressWarnings(
    fitMvglsAndExtractBIC.formula(form_chr, painted_tree, dat, data = dat, method = "LL")
  )
  r2 <- suppressWarnings(
    fitMvglsAndExtractBIC.formula(form_chr, painted_tree, dat, data = dat, method = "LL")
  )

  b1 <- get_bic_scalar(r1$BIC)
  b2 <- get_bic_scalar(r2$BIC)
  expect_equal(b1, b2, tolerance = 1e-10)
})

# Test: BM/BMM branch selection aligns with regime count (BIC)
test_that("BM/BMM branch selection aligns with regime count (BIC)", {
  skip_if_missing_deps()
  tr_one <- make_single_regime_tree(20)
  tr_two <- make_two_regime_tree(20)

  expect_equal(length(unique(phytools::getStates(tr_one))), 1)
  expect_gt(length(unique(phytools::getStates(tr_two))), 1)

  set.seed(135)
  d1 <- make_data_for_tree(tr_one, p = 2)
  set.seed(136)
  d2 <- make_data_for_tree(tr_two, p = 2)
  f <- "cbind(y1, y2) ~ x"

  r1 <- suppressWarnings(fitMvglsAndExtractBIC.formula(f, tr_one, d1, data = d1, method = "LL"))
  r2 <- suppressWarnings(fitMvglsAndExtractBIC.formula(f, tr_two, d2, data = d2, method = "LL"))

  expect_true(inherits(r1$model, "mvgls"))
  expect_true(is.finite(get_bic_scalar(r1$BIC)))
  expect_true(inherits(r2$model, "mvgls"))
  expect_true(is.finite(get_bic_scalar(r2$BIC)))
})

# Test: mvgls options via ... are honoured for BIC (LL with REML toggles)
test_that("mvgls options via ... are honoured for BIC (LL with REML toggles)", {
  skip_if_missing_deps()
  painted_tree <- make_two_regime_tree(20)
  set.seed(42)
  dat <- make_data_for_tree(painted_tree, p = 2)
  f <- "cbind(y1, y2) ~ x"

  # Use method = "LL" for BIC; vary REML
  r_reml_true  <- suppressWarnings(
    fitMvglsAndExtractBIC.formula(f, painted_tree, dat, data = dat, method = "LL", REML = TRUE)
  )
  r_reml_false <- suppressWarnings(
    fitMvglsAndExtractBIC.formula(f, painted_tree, dat, data = dat, method = "LL", REML = FALSE)
  )

  expect_true(is.finite(get_bic_scalar(r_reml_true$BIC)))
  expect_true(is.finite(get_bic_scalar(r_reml_false$BIC)))
})

# Test: Rownames must match in name and order (not just set) for BIC
test_that("Rownames must match in name and order (not just set) for BIC", {
  skip_if_missing_deps()
  painted_tree <- make_single_regime_tree(15)
  set.seed(7)
  dat <- make_data_for_tree(painted_tree, p = 2)
  shuffled <- dat[sample(nrow(dat)), , drop = FALSE]  # same names, different order
  f <- "cbind(y1, y2) ~ x"

  expect_error(
    fitMvglsAndExtractBIC.formula(f, painted_tree, shuffled, data = shuffled),
    "Row names of trait_data must exactly match the tip labels of the tree\\."
  )
})
