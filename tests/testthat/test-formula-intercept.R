test_that("formula normalization preserves the requested design without an intercept", {
  trait_data <- data.frame(
    y1 = c(1, 3, 2, 5, 4, 6),
    y2 = c(3, 1, 4, 2, 6, 5),
    x = c(2, 4, 1, 6, 3, 5),
    group = factor(rep(c("a", "b", "c"), 2))
  )
  formulas <- list(
    cbind(y1, y2) ~ 0 + x,
    cbind(y1, y2) ~ x - 1,
    cbind(y1, y2) ~ 0 + group,
    cbind(y1, y2) ~ group - 1,
    cbind(y1, y2) ~ 0 + x:group,
    trait_data[, 1:2] ~ 0 + trait_data[, 3],
    trait_data[, 1:2] ~ trait_data[, 3] - 1,
    cbind(y1, y2) ~ 0
  )
  reference_formulas <- formulas
  reference_formulas[[6L]] <- ~ 0 + x
  reference_formulas[[7L]] <- ~ x - 1
  for (i in seq_along(formulas)) {
    normalized <- normalizeMvglsFormulaCall(formulas[[i]], trait_data, list())
    expected <- stats::model.matrix(reference_formulas[[i]], trait_data)
    actual <- stats::model.matrix(normalized$formula, normalized$args_list$data)
    expect_equal(dim(actual), dim(expected))
    expect_equal(as.numeric(actual), as.numeric(expected))
  }
})

test_that("search fits preserve no-intercept numeric and factor models", {
  withr::local_seed(12)
  tree <- ape::rtree(30)
  dat <- data.frame(
    y1 = rnorm(30), y2 = rnorm(30), x = rnorm(30),
    group = factor(rep(c("a", "b", "c"), 10)),
    row.names = tree$tip.label
  )
  formulas <- list(
    cbind(y1, y2) ~ 0 + x,
    cbind(y1, y2) ~ x - 1,
    cbind(y1, y2) ~ 0 + group,
    cbind(y1, y2) ~ x + group
  )
  for (formula in formulas) {
    reference <- mvMORPH::mvgls(
      formula, data = dat, tree = tree, model = "BM",
      method = "LL", REML = FALSE
    )
    for (ic in c("BIC", "GIC")) {
      expect_warning(
        result <- searchOptimalConfiguration(
          tree, dat, formula = formula, min_descendant_tips = 30,
          IC = ic, method = "LL", REML = FALSE,
          num_cores = 1, progress = FALSE
        ),
        "No non-root internal nodes"
      )
      expect_equal(result$model_no_uncertainty$coefficients,
                   reference$coefficients, tolerance = 1e-8)
      expect_equal(result$model_no_uncertainty$sigma$Pinv,
                   reference$sigma$Pinv, tolerance = 1e-8)
      expected_ic <- if (ic == "BIC") {
        stats::BIC(reference)$BIC
      } else {
        mvMORPH::GIC(reference)$GIC
      }
      expect_equal(as.numeric(result$baseline_ic), as.numeric(expected_ic),
                   tolerance = 1e-8)
    }
  }
})
