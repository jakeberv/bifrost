testthat::skip_on_cran()

skip_if_null_sim_deps <- function() {
  testthat::skip_if_not_installed("ape")
  testthat::skip_if_not_installed("phytools")
  testthat::skip_if_not_installed("mvMORPH")
}

make_null_template <- function() {
  set.seed(10)
  tr <- ape::rtree(20)
  X <- matrix(rnorm(20 * 3), ncol = 3)
  rownames(X) <- tr$tip.label
  createSimulationTemplate(
    baseline_tree = tr,
    trait_data = X,
    formula = "trait_data ~ 1",
    method = "LL"
  )
}

make_formula_null_template <- function() {
  set.seed(11)
  tr <- ape::rtree(20)
  X <- cbind(y1 = rnorm(20), y2 = rnorm(20), mass = rnorm(20))
  rownames(X) <- tr$tip.label
  createSimulationTemplate(
    baseline_tree = tr,
    trait_data = X,
    formula = "trait_data[, 1:2] ~ trait_data[, 3]",
    response_columns = 1:2,
    predictor_columns = 3,
    method = "LL"
  )
}

make_single_trait_null_template <- function() {
  set.seed(12)
  tr <- ape::rtree(14)
  X <- matrix(rnorm(14), ncol = 1, dimnames = list(tr$tip.label, "y1"))
  structure(
    list(
      baseline_tree = tr,
      trait_data = X,
      fitted_values = matrix(0, nrow = 14, ncol = 1, dimnames = list(tr$tip.label, "y1")),
      response_column_names = "y1",
      response_columns = 1L,
      predictor_columns = integer(0),
      trait_data_is_matrix = TRUE,
      variance_mean = 0.2,
      variance_sd = 0,
      covariance_mean = 0,
      covariance_sd = 0,
      n_response_traits = 1L,
      n_tips = ape::Ntip(tr)
    ),
    class = c("bifrost_simulation_template", "list")
  )
}

make_manual_null_template_df <- function(with_predictor = FALSE) {
  set.seed(13)
  tr <- ape::rtree(10)
  if (with_predictor) {
    dat <- data.frame(
      y1 = rnorm(10),
      y2 = rnorm(10),
      mass = rnorm(10),
      row.names = tr$tip.label
    )
    fitted <- matrix(0, nrow = 10, ncol = 2, dimnames = list(tr$tip.label, c("y1", "y2")))
    response_columns <- 1:2
    predictor_columns <- 3L
    response_names <- c("y1", "y2")
    n_response_traits <- 2L
  } else {
    dat <- data.frame(
      y1 = rnorm(10),
      y2 = rnorm(10),
      row.names = tr$tip.label
    )
    fitted <- matrix(0, nrow = 10, ncol = 2, dimnames = list(tr$tip.label, c("y1", "y2")))
    response_columns <- 1:2
    predictor_columns <- integer(0)
    response_names <- c("y1", "y2")
    n_response_traits <- 2L
  }

  structure(
    list(
      baseline_tree = tr,
      trait_data = dat,
      fitted_values = fitted,
      response_column_names = response_names,
      response_columns = response_columns,
      predictor_columns = predictor_columns,
      trait_data_is_matrix = FALSE,
      variance_mean = 0.2,
      variance_sd = 0.01,
      covariance_mean = 0.05,
      covariance_sd = 0.01,
      n_response_traits = n_response_traits,
      n_tips = ape::Ntip(tr)
    ),
    class = c("bifrost_simulation_template", "list")
  )
}

make_zero_cov_null_template <- function() {
  set.seed(18)
  tr <- ape::rtree(10)
  X <- matrix(0, nrow = 10, ncol = 2, dimnames = list(tr$tip.label, c("y1", "y2")))
  structure(
    list(
      baseline_tree = tr,
      trait_data = X,
      fitted_values = X,
      response_column_names = c("y1", "y2"),
      response_columns = 1:2,
      predictor_columns = integer(0),
      trait_data_is_matrix = TRUE,
      variance_mean = 0,
      variance_sd = 0,
      covariance_mean = 0,
      covariance_sd = 0,
      n_response_traits = 2L,
      n_tips = ape::Ntip(tr)
    ),
    class = c("bifrost_simulation_template", "list")
  )
}

test_that("simulateNullDataset returns an aligned null replicate", {
  skip_if_null_sim_deps()

  tmpl <- make_null_template()
  sim <- simulateNullDataset(tmpl, tree_tip_count = 15, seed = 2)

  testthat::expect_s3_class(sim, "bifrost_simulation_replicate_null")
  testthat::expect_equal(length(sim$tree$tip.label), 15)
  testthat::expect_identical(rownames(sim$trait_data), sim$tree$tip.label)
  testthat::expect_equal(dim(sim$simulatedData), c(15, 3))
  testthat::expect_true(all(eigen(sim$covariance_matrix, symmetric = TRUE)$values > 0))
  testthat::expect_identical(sim$generating_scenario, "null")
})

test_that("simulateNullDataset preserves predictor columns for formula templates", {
  skip_if_null_sim_deps()

  tmpl <- make_formula_null_template()
  sim <- simulateNullDataset(tmpl, tree_tip_count = 16, seed = 3)

  testthat::expect_identical(
    sim$trait_data[, 3],
    tmpl$trait_data[sim$tree$tip.label, 3]
  )
  testthat::expect_equal(dim(sim$simulatedData), c(16, 2))
})

test_that("simulateNullDataset requires predictors for formula templates", {
  skip_if_null_sim_deps()

  tmpl <- make_formula_null_template()

  testthat::expect_error(
    simulateNullDataset(tmpl, tree_tip_count = 16, seed = 3, preserve_predictors = FALSE),
    "preserve_predictors must be TRUE"
  )
})

test_that("simulateNullDataset validates tree_tip_count", {
  skip_if_null_sim_deps()

  tmpl <- make_null_template()
  testthat::expect_error(
    simulateNullDataset(tmpl, tree_tip_count = 100),
    "cannot exceed"
  )
})

test_that("simulateNullDataset validates template and argument types", {
  skip_if_null_sim_deps()

  tmpl <- make_null_template()

  testthat::expect_error(
    simulateNullDataset(list()),
    "bifrost_simulation_template"
  )
  testthat::expect_error(
    simulateNullDataset(tmpl, tree_tip_count = 1),
    "single integer >= 2"
  )
  testthat::expect_error(
    simulateNullDataset(tmpl, preserve_predictors = NA),
    "TRUE or FALSE"
  )
})

test_that("simulateNullDataset supports full-tree and single-trait templates", {
  skip_if_null_sim_deps()

  tmpl <- make_single_trait_null_template()
  sim <- simulateNullDataset(tmpl, seed = 14)

  testthat::expect_equal(length(sim$tree$tip.label), tmpl$n_tips)
  testthat::expect_equal(dim(sim$simulatedData), c(tmpl$n_tips, 1))
  testthat::expect_equal(dim(sim$covariance_matrix), c(1, 1))
})

test_that("simulateNullDataset preserves data.frame outputs when requested", {
  skip_if_null_sim_deps()

  sim_no_predictor <- simulateNullDataset(make_manual_null_template_df(FALSE), seed = 15)
  sim_with_predictor <- simulateNullDataset(make_manual_null_template_df(TRUE), seed = 16)

  testthat::expect_true(is.data.frame(sim_no_predictor$trait_data))
  testthat::expect_true(is.data.frame(sim_with_predictor$trait_data))
  testthat::expect_identical(
    sim_with_predictor$trait_data$mass,
    make_manual_null_template_df(TRUE)$trait_data[sim_with_predictor$tree$tip.label, "mass"]
  )
})

test_that("simulateNullDataset errors when covariance generation never succeeds", {
  skip_if_null_sim_deps()

  tmpl <- make_zero_cov_null_template()

  testthat::expect_error(
    simulateNullDataset(tmpl, seed = 17),
    "positive-definite covariance matrix"
  )
})
