# tests/testthat/test-fitMvglsAndExtractIC.formula.R

skip_if_ic_formula_deps <- function() {
  testthat::skip_if_not_installed("ape")
  testthat::skip_if_not_installed("phytools")
  testthat::skip_if_not_installed("mvMORPH")
}

make_ic_formula_tree <- function(n_tip = 20, seed = 101, regimes = c("single", "multi")) {
  regimes <- match.arg(regimes)
  set.seed(seed)
  tree <- phytools::pbtree(n = n_tip, scale = 1)
  painted <- phytools::paintSubTree(
    tree,
    node = ape::Ntip(tree) + 1L,
    state = "0",
    anc.state = "0",
    stem = FALSE
  )
  if (identical(regimes, "multi")) {
    painted <- phytools::paintSubTree(
      painted,
      node = ape::Ntip(painted) + 2L,
      state = "1",
      stem = TRUE
    )
  }
  painted
}

make_ic_formula_data <- function(tree, n_response = 2, seed = 202) {
  set.seed(seed)
  n <- ape::Ntip(tree)
  size <- exp(stats::rnorm(n))
  group <- factor(rep(c("a", "b"), length.out = n))
  responses <- replicate(n_response, 0.5 * log(size) + stats::rnorm(n, sd = 0.4))
  out <- data.frame(size = size, group = group, responses)
  names(out) <- c("size", "group", paste0("y", seq_len(n_response)))
  rownames(out) <- tree$tip.label
  out
}

ic_formula_specs <- function() {
  list(
    GIC = list(
      fit = fitMvglsAndExtractGIC.formula,
      field = "GIC",
      scalar = function(x) as.numeric(x$GIC$GIC),
      args = list()
    ),
    BIC = list(
      fit = fitMvglsAndExtractBIC.formula,
      field = "BIC",
      scalar = function(x) as.numeric(x$BIC$BIC),
      args = list(method = "LL")
    )
  )
}

fit_ic_formula <- function(spec, formula, tree, data, ...) {
  do.call(
    spec$fit,
    c(
      list(formula = formula, painted_tree = tree, trait_data = data),
      spec$args,
      list(...)
    )
  )
}

for (ic_name in names(ic_formula_specs())) {
  local({
    label <- ic_name
    spec <- ic_formula_specs()[[label]]

    testthat::test_that(paste(label, "formula wrapper fits single- and multi-regime trees"), {
      skip_if_ic_formula_deps()

      single_tree <- make_ic_formula_tree(regimes = "single", seed = 11)
      multi_tree <- make_ic_formula_tree(regimes = "multi", seed = 12)
      single_data <- make_ic_formula_data(single_tree, n_response = 2, seed = 21)
      multi_data <- make_ic_formula_data(multi_tree, n_response = 3, seed = 22)

      single <- suppressWarnings(fit_ic_formula(
        spec,
        cbind(y1, y2) ~ size,
        single_tree,
        single_data
      ))
      multi <- suppressWarnings(fit_ic_formula(
        spec,
        "cbind(y1, y2, y3) ~ size",
        multi_tree,
        multi_data
      ))

      testthat::expect_named(single, c("model", spec$field))
      testthat::expect_s3_class(single$model, "mvgls")
      testthat::expect_true(is.finite(spec$scalar(single)))
      testthat::expect_named(multi, c("model", spec$field))
      testthat::expect_s3_class(multi$model, "mvgls")
      testthat::expect_true(is.finite(spec$scalar(multi)))
    })

    testthat::test_that(paste(label, "formula wrapper injects trait_data for named formula objects"), {
      skip_if_ic_formula_deps()

      tree <- make_ic_formula_tree(regimes = "single", seed = 13)
      data <- make_ic_formula_data(tree, n_response = 2, seed = 23)

      fit <- suppressWarnings(fit_ic_formula(
        spec,
        cbind(y1, y2) ~ log(size) * group,
        tree,
        data
      ))

      testthat::expect_s3_class(fit$model, "mvgls")
      testthat::expect_true(is.finite(spec$scalar(fit)))
    })

    testthat::test_that(paste(label, "formula wrapper validates formula and row alignment"), {
      skip_if_ic_formula_deps()

      tree <- make_ic_formula_tree(regimes = "multi", seed = 14)
      data <- make_ic_formula_data(tree, n_response = 2, seed = 24)
      shuffled <- data[rev(seq_len(nrow(data))), , drop = FALSE]
      missing_name <- data
      rownames(missing_name)[1] <- "not_a_tip"

      testthat::expect_error(
        fit_ic_formula(spec, 1, tree, data),
        "character formula or formula object"
      )
      testthat::expect_error(
        fit_ic_formula(spec, cbind(y1, y2) ~ size, tree, shuffled),
        "Row names of trait_data must exactly match the tip labels of the tree"
      )
      testthat::expect_error(
        fit_ic_formula(spec, cbind(y1, y2) ~ size, tree, missing_name),
        "Row names of trait_data must exactly match the tip labels of the tree"
      )
    })
  })
}

testthat::test_that("GIC formula wrapper forwards explicit LL and PL-LOOCV methods", {
  skip_if_ic_formula_deps()

  # Test method forwarding with enough tips per regime for stable starting
  # values. A two-tip regime can saturate the local regression used by mvgls.
  tree <- ape::stree(32, type = "balanced")
  tree$edge.length <- rep(1, nrow(tree$edge))
  root <- ape::Ntip(tree) + 1L
  tree <- phytools::paintSubTree(
    tree, node = root, state = "0", anc.state = "0", stem = FALSE
  )
  root_child <- tree$edge[tree$edge[, 1] == root, 2][1]
  tree <- phytools::paintSubTree(tree, node = root_child, state = "1", stem = TRUE)
  testthat::expect_equal(
    as.integer(table(phytools::getStates(tree, "tips"))), c(16L, 16L)
  )

  set.seed(42)
  size <- rnorm(32)
  residuals <- t(chol(ape::vcv.phylo(tree))) %*% matrix(rnorm(64), ncol = 2)
  data <- data.frame(
    size = size,
    y1 = 0.5 * size + residuals[, 1],
    y2 = 0.5 * size + residuals[, 2],
    row.names = tree$tip.label
  )

  expected_methods <- c(LL = "LL", "PL-LOOCV" = "LOOCV")
  for (method in names(expected_methods)) {
    fit <- fitMvglsAndExtractGIC.formula(
      cbind(y1, y2) ~ size, tree, data, method = method
    )

    testthat::expect_s3_class(fit$model, "mvgls")
    testthat::expect_identical(fit$model$call$method, method)
    testthat::expect_identical(fit$model$method, unname(expected_methods[[method]]))
    testthat::expect_true(is.finite(as.numeric(fit$GIC$GIC)))
  }
})

testthat::test_that("GIC formula wrapper preserves the mvgls default method", {
  skip_if_ic_formula_deps()

  tree <- make_ic_formula_tree(regimes = "multi", seed = 19)
  data <- make_ic_formula_data(tree, n_response = 2, seed = 29)

  fit <- suppressWarnings(fitMvglsAndExtractGIC.formula(
    cbind(y1, y2) ~ size,
    tree,
    data
  ))

  testthat::expect_null(fit$model$call$method)
  testthat::expect_identical(fit$model$method, "LOOCV")
  testthat::expect_true(is.finite(as.numeric(fit$GIC$GIC)))
})

testthat::test_that("BIC formula wrapper defaults to LL when method is omitted", {
  skip_if_ic_formula_deps()

  tree <- make_ic_formula_tree(regimes = "multi", seed = 16)
  data <- make_ic_formula_data(tree, n_response = 2, seed = 26)

  fit <- suppressWarnings(fitMvglsAndExtractBIC.formula(
    cbind(y1, y2) ~ size,
    tree,
    data
  ))

  testthat::expect_s3_class(fit$model, "mvgls")
  testthat::expect_identical(fit$model$call$method, "LL")
  testthat::expect_identical(fit$model$method, "LL")
  testthat::expect_true(is.finite(as.numeric(fit$BIC$BIC)))
})

testthat::test_that("BIC formula wrapper accepts LL REML toggles", {
  skip_if_ic_formula_deps()

  tree <- make_ic_formula_tree(regimes = "multi", seed = 18)
  data <- make_ic_formula_data(tree, n_response = 2, seed = 28)

  reml <- suppressWarnings(fitMvglsAndExtractBIC.formula(
    cbind(y1, y2) ~ size,
    tree,
    data,
    method = "LL",
    REML = TRUE
  ))
  ml <- suppressWarnings(fitMvglsAndExtractBIC.formula(
    cbind(y1, y2) ~ size,
    tree,
    data,
    method = "LL",
    REML = FALSE
  ))

  testthat::expect_true(is.finite(as.numeric(reml$BIC$BIC)))
  testthat::expect_true(is.finite(as.numeric(ml$BIC$BIC)))
})

testthat::test_that("IC formula wrappers are deterministic for repeated fits", {
  skip_if_ic_formula_deps()

  tree <- make_ic_formula_tree(regimes = "multi", seed = 19)
  data <- make_ic_formula_data(tree, n_response = 2, seed = 29)
  formula <- cbind(y1, y2) ~ size

  for (spec in ic_formula_specs()) {
    first <- suppressWarnings(fit_ic_formula(spec, formula, tree, data))
    second <- suppressWarnings(fit_ic_formula(spec, formula, tree, data))
    testthat::expect_equal(
      spec$scalar(first),
      spec$scalar(second),
      tolerance = 1e-10
    )
  }
})
