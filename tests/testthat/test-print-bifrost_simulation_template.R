testthat::skip_on_cran()

skip_if_print_template_deps <- function() {
  testthat::skip_if_not_installed("ape")
  testthat::skip_if_not_installed("phytools")
  testthat::skip_if_not_installed("mvMORPH")
}

test_that("print.bifrost_simulation_template prints the key sections", {
  skip_if_print_template_deps()

  set.seed(60)
  tr <- ape::rtree(18)
  X <- matrix(rnorm(18 * 3), ncol = 3)
  rownames(X) <- tr$tip.label
  tmpl <- createSimulationTemplate(tr, X, formula = "trait_data ~ 1", method = "LL")

  out <- paste(capture.output(print(tmpl)), collapse = "\n")

  testthat::expect_match(out, "Bifrost Simulation Template")
  testthat::expect_match(out, "Empirical Covariance Summaries")
  testthat::expect_match(out, "Calibration formula:")
  testthat::expect_match(out, "Simulation search:")
  testthat::expect_match(out, "Tips:")
})

test_that("print.bifrost_simulation_template prints the error-term branch", {
  skip_if_print_template_deps()

  set.seed(61)
  tr <- ape::rtree(16)
  X <- matrix(rnorm(16 * 2), ncol = 2)
  rownames(X) <- tr$tip.label
  tmpl <- createSimulationTemplate(tr, X, formula = "trait_data ~ 1", method = "LL", error = TRUE)

  out <- paste(capture.output(print(tmpl)), collapse = "\n")

  testthat::expect_match(out, "Error term:")
})

test_that("print.bifrost_simulation_template uses stored evaluated fit settings", {
  skip_if_print_template_deps()

  set.seed(62)
  tr <- ape::rtree(16)
  X <- matrix(rnorm(16 * 2), ncol = 2)
  rownames(X) <- tr$tip.label
  method_value <- "LL"
  error_value <- TRUE
  tmpl <- createSimulationTemplate(
    tr,
    X,
    formula = "trait_data ~ 1",
    method = method_value,
    error = error_value
  )

  out <- paste(capture.output(print(tmpl)), collapse = "\n")

  testthat::expect_match(out, "Method: LL")
  testthat::expect_match(out, "Error term: TRUE")
})

test_that("print.bifrost_simulation_template falls back to model call settings", {
  skip_if_print_template_deps()

  set.seed(63)
  tr <- ape::rtree(16)
  X <- matrix(rnorm(16 * 2), ncol = 2)
  rownames(X) <- tr$tip.label
  tmpl <- createSimulationTemplate(tr, X, formula = "trait_data ~ 1", method = "LL", error = TRUE)
  tmpl$fit_method <- NULL
  tmpl$fit_error <- NULL

  out <- paste(capture.output(print(tmpl)), collapse = "\n")

  testthat::expect_match(out, "Method: LL")
  testthat::expect_match(out, "Error term: TRUE")
})

test_that("print.bifrost_simulation_template ignores unevaluable model call errors", {
  skip_if_print_template_deps()

  set.seed(64)
  tr <- ape::rtree(16)
  X <- matrix(rnorm(16 * 2), ncol = 2)
  rownames(X) <- tr$tip.label
  tmpl <- createSimulationTemplate(tr, X, formula = "trait_data ~ 1", method = "LL")
  tmpl$fit_method <- NULL
  tmpl$fit_error <- NULL
  tmpl$global_model$call$error <- quote(missing_error_flag)

  out <- paste(capture.output(print(tmpl)), collapse = "\n")

  testthat::expect_match(out, "Method: LL")
  testthat::expect_no_match(out, "Error term:")
})
