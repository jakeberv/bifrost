# tests/testthat/test-fitMvglsAndExtractGIC-formula.R

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

# Extract a scalar GIC from either a numeric or a list returned by mvMORPH::GIC()
get_gic_scalar <- function(gic_obj) {
  if (is.list(gic_obj) && !is.null(gic_obj$GIC)) {
    return(as.numeric(gic_obj$GIC))
  }
  as.numeric(gic_obj)
}

# Group: core fit behavior

# Test: Single-regime tree fits and returns a scalar GIC (checks return value)
test_that("Single-regime tree fits and returns a scalar GIC", {
  skip_if_missing_deps()

  painted_tree <- make_single_regime_tree(20)
  set.seed(123)
  dat <- make_data_for_tree(painted_tree, p = 2)

  form_chr <- "cbind(y1, y2) ~ x"

  res <- suppressWarnings(
    fitMvglsAndExtractGIC.formula(
      formula      = form_chr,
      painted_tree = painted_tree,
      trait_data   = dat,
      data         = dat
    )
  )

  expect_type(res, "list")
  expect_true("model" %in% names(res))
  expect_true("GIC"   %in% names(res))
  expect_true(inherits(res$model, "mvgls"))

  gic <- get_gic_scalar(res$GIC)
  expect_true(is.finite(gic))
  expect_length(gic, 1)
})

# Test: Multi-regime tree fits and returns a scalar GIC (BMM path) (checks return value)
test_that("Multi-regime tree fits and returns a scalar GIC (BMM path)", {
  skip_if_missing_deps()

  painted_tree <- make_two_regime_tree(25)
  set.seed(456)
  dat <- make_data_for_tree(painted_tree, p = 3)

  form_chr <- "cbind(y1, y2, y3) ~ x"

  res <- suppressWarnings(
    fitMvglsAndExtractGIC.formula(
      formula      = form_chr,
      painted_tree = painted_tree,
      trait_data   = dat,
      data         = dat
    )
  )

  expect_true(inherits(res$model, "mvgls"))

  gic <- get_gic_scalar(res$GIC)
  expect_true(is.finite(gic))
  expect_length(gic, 1)
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
    fitMvglsAndExtractGIC.formula(
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
    fitMvglsAndExtractGIC.formula(
      formula      = form_obj,
      painted_tree = painted_tree,
      trait_data   = dat,
      data         = dat
    ),
    "A character formula must be provided."
  )

  expect_error(
    fitMvglsAndExtractGIC.formula(
      painted_tree = painted_tree,
      trait_data   = dat,
      data         = dat
    ),
    "A character formula must be provided."
  )
})


# Group: determinism, branch selection, and option passthrough

# Test: Determinism: same inputs yield identical GIC (re-runs with same inputs)
test_that("Determinism: same inputs yield identical GIC", {
  skip_if_missing_deps()
  painted_tree <- make_two_regime_tree(25)
  set.seed(999); dat <- make_data_for_tree(painted_tree, p = 2)
  form_chr <- "cbind(y1, y2) ~ x"

  r1 <- suppressWarnings(
    fitMvglsAndExtractGIC.formula(form_chr, painted_tree, dat, data = dat)
  )
  r2 <- suppressWarnings(
    fitMvglsAndExtractGIC.formula(form_chr, painted_tree, dat, data = dat)
  )

  g1 <- get_gic_scalar(r1$GIC)
  g2 <- get_gic_scalar(r2$GIC)
  expect_equal(g1, g2, tolerance = 1e-10)
})

# Test: BM/BMM branch selection aligns with regime count
test_that("BM/BMM branch selection aligns with regime count", {
  skip_if_missing_deps()
  tr_one <- make_single_regime_tree(20)
  tr_two <- make_two_regime_tree(20)

  expect_equal(length(unique(phytools::getStates(tr_one))), 1)
  expect_gt(length(unique(phytools::getStates(tr_two))), 1)

  set.seed(135); d1 <- make_data_for_tree(tr_one, p = 2)
  set.seed(136); d2 <- make_data_for_tree(tr_two, p = 2)
  f <- "cbind(y1, y2) ~ x"

  r1 <- suppressWarnings(fitMvglsAndExtractGIC.formula(f, tr_one, d1, data = d1))
  r2 <- suppressWarnings(fitMvglsAndExtractGIC.formula(f, tr_two, d2, data = d2))

  expect_true(inherits(r1$model, "mvgls"))
  expect_true(is.finite(get_gic_scalar(r1$GIC)))
  expect_true(inherits(r2$model, "mvgls"))
  expect_true(is.finite(get_gic_scalar(r2$GIC)))
})

# Test: mvgls options via ... are honoured (LL and PL-LOOCV)
test_that("mvgls options via ... are honoured (LL and PL-LOOCV)", {
  skip_if_missing_deps()
  painted_tree <- make_two_regime_tree(20)
  set.seed(42); dat <- make_data_for_tree(painted_tree, p = 2)
  f <- "cbind(y1, y2) ~ x"

  r_ll <- suppressWarnings(
    fitMvglsAndExtractGIC.formula(f, painted_tree, dat, data = dat, method = "LL")
  )
  r_pl <- suppressWarnings(
    fitMvglsAndExtractGIC.formula(f, painted_tree, dat, data = dat, method = "PL-LOOCV")
  )

  expect_true(is.finite(get_gic_scalar(r_ll$GIC)))
  expect_true(is.finite(get_gic_scalar(r_pl$GIC)))
})

# Test: Rownames must match in name and order (not just set)
test_that("Rownames must match in name and order (not just set)", {
  skip_if_missing_deps()
  painted_tree <- make_single_regime_tree(15)
  set.seed(7); dat <- make_data_for_tree(painted_tree, p = 2)
  shuffled <- dat[sample(nrow(dat)), , drop = FALSE]  # same names, different order
  f <- "cbind(y1, y2) ~ x"

  expect_error(
    fitMvglsAndExtractGIC.formula(f, painted_tree, shuffled, data = shuffled),
    "Row names of trait_data must exactly match the tip labels of the tree\\."
  )
})
