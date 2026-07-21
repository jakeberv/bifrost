# tests/testthat/test-mvgls-functions.R

skip_if_mvgls_function_deps <- function() {
  testthat::skip_if_not_installed("ape")
  testthat::skip_if_not_installed("phytools")
  testthat::skip_if_not_installed("mvMORPH")
}

make_mvgls_function_tree <- function(n_tip = 18, seed = 301) {
  set.seed(seed)
  tree <- phytools::pbtree(n = n_tip, scale = 1)
  baseline <- phytools::paintSubTree(
    tree,
    node = ape::Ntip(tree) + 1L,
    state = "0",
    anc.state = "0",
    stem = FALSE
  )
  phytools::paintSubTree(
    baseline,
    node = ape::Ntip(baseline) + 2L,
    state = "1",
    stem = TRUE
  )
}

make_mvgls_trait_matrix <- function(tree, n_traits = 2, seed = 302) {
  set.seed(seed)
  n <- ape::Ntip(tree)
  x <- stats::rnorm(n)
  out <- replicate(n_traits, x + stats::rnorm(n, sd = 0.5))
  colnames(out) <- paste0("trait", seq_len(n_traits))
  rownames(out) <- tree$tip.label
  out
}

mvgls_ic_specs <- function() {
  list(
    GIC = list(
      fit = fitMvglsAndExtractGIC,
      field = "GIC",
      scalar = function(x) as.numeric(x$GIC$GIC),
      expected_length = 7L
    ),
    BIC = list(
      fit = fitMvglsAndExtractBIC,
      field = "BIC",
      scalar = function(x) as.numeric(x$BIC$BIC),
      expected_length = 4L
    )
  )
}

for (ic_name in names(mvgls_ic_specs())) {
  local({
    label <- ic_name
    spec <- mvgls_ic_specs()[[label]]

    testthat::test_that(paste(label, "intercept wrapper fits a BMM model and returns a scalar criterion"), {
      skip_if_mvgls_function_deps()

      tree <- make_mvgls_function_tree(n_tip = 18, seed = 31)
      trait_data <- make_mvgls_trait_matrix(tree, n_traits = 2, seed = 32)

      fit <- suppressWarnings(spec$fit(tree, trait_data))

      testthat::expect_named(fit, c("model", spec$field))
      testthat::expect_s3_class(fit$model, "mvgls")
      testthat::expect_length(fit[[spec$field]], spec$expected_length)
      testthat::expect_true(is.finite(spec$scalar(fit)))
    })

    testthat::test_that(paste(label, "intercept wrapper supports higher-dimensional response matrices"), {
      skip_if_mvgls_function_deps()

      tree <- make_mvgls_function_tree(n_tip = 24, seed = 33)
      trait_data <- make_mvgls_trait_matrix(tree, n_traits = 5, seed = 34)

      fit <- suppressWarnings(spec$fit(tree, trait_data))

      testthat::expect_s3_class(fit$model, "mvgls")
      testthat::expect_true(is.finite(spec$scalar(fit)))
    })

    testthat::test_that(paste(label, "intercept wrapper validates matrix shape and row alignment"), {
      skip_if_mvgls_function_deps()

      tree <- make_mvgls_function_tree(n_tip = 18, seed = 35)
      trait_data <- make_mvgls_trait_matrix(tree, n_traits = 2, seed = 36)
      data_frame <- as.data.frame(trait_data)
      univariate <- trait_data[, 1, drop = FALSE]
      mismatched <- trait_data
      rownames(mismatched)[1] <- "not_a_tip"
      unnamed <- trait_data
      rownames(unnamed) <- NULL

      testthat::expect_error(
        spec$fit(tree, data_frame),
        "trait_data must be a matrix"
      )
      testthat::expect_error(
        suppressWarnings(spec$fit(tree, univariate)),
        "multivariate datasets"
      )
      testthat::expect_error(
        spec$fit(tree, mismatched),
        "Row names of trait_data must exactly match the tip labels of the tree"
      )
      testthat::expect_error(
        spec$fit(tree, unnamed),
        "Row names of trait_data must exactly match the tip labels of the tree"
      )
    })
  })
}

testthat::test_that("GIC and BIC intercept wrappers fit the same underlying model", {
  skip_if_mvgls_function_deps()

  tree <- make_mvgls_function_tree(n_tip = 18, seed = 37)
  trait_data <- make_mvgls_trait_matrix(tree, n_traits = 2, seed = 38)

  gic <- suppressWarnings(fitMvglsAndExtractGIC(tree, trait_data))
  bic <- suppressWarnings(fitMvglsAndExtractBIC(tree, trait_data))

  testthat::expect_equal(stats::logLik(gic$model), stats::logLik(bic$model))
  testthat::expect_true(is.finite(as.numeric(gic$GIC$GIC)))
  testthat::expect_true(is.finite(as.numeric(bic$BIC$BIC)))
})
