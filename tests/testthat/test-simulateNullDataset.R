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
  testthat::skip_on_cran()
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

test_that("simulateNullDataset supports an explicit empirical Wishart design", {
  testthat::skip_on_cran()
  skip_if_null_sim_deps()

  tmpl <- make_null_template()
  sim <- simulateNullDataset(
    tmpl,
    tree_tip_count = 15,
    simulation_generator = "empirical",
    seed = 2
  )

  testthat::expect_identical(sim$simulation_generator, "empirical")
  testthat::expect_false("simulation_design" %in% names(sim))
  testthat::expect_identical(sim$covariance_df, tmpl$residual_df)
  testthat::expect_named(
    sim$covariance_diagnostics,
    c(
      "mean_absolute_correlation",
      "effective_dimensionality",
      "condition_number"
    )
  )
  testthat::expect_true(all(is.finite(unlist(sim$covariance_diagnostics))))
  testthat::expect_gte(tmpl$residual_df, tmpl$n_response_traits + 1L)
})

test_that("the direct simulator rejects generator vectors before RNG setup", {
  skip_if_null_sim_deps()
  rng_setup_called <- FALSE
  testthat::local_mocked_bindings(
    .simulation_set_seed = function(...) {
      rng_setup_called <<- TRUE
      stats::runif(1)
      stop("RNG setup reached", call. = FALSE)
    }
  )
  tmpl <- make_null_template()
  set.seed(20260709)
  before <- .Random.seed
  testthat::expect_error(
    simulateNullDataset(
      tmpl,
      seed = 2,
      simulation_generator = c("original", "empirical")
    ),
    "simulation_generator.*original.*empirical"
  )
  testthat::expect_false(rng_setup_called)
  testthat::expect_identical(.Random.seed, before)
})

test_that("simulateNullDataset defaults to the original covariance draw", {
  testthat::skip_on_cran()
  skip_if_null_sim_deps()

  expected_sigma <- matrix(
    c(
      2.42880177391557, 0.295943316837243, 0.438489181274454,
      0.295943316837243, 1.32591624061694, 0.300309524421243,
      0.438489181274454, 0.300309524421243, 1.84819988452061
    ),
    nrow = 3L
  )
  sim <- simulateNullDataset(
    make_null_template(),
    tree_tip_count = 15,
    seed = 2
  )

  testthat::expect_identical(sim$simulation_generator, "original")
  testthat::expect_true(is.na(sim$covariance_df))
  testthat::expect_equal(sim$covariance_matrix, expected_sigma, tolerance = 1e-14)
})

test_that("simulateNullDataset supports its legacy positional seed", {
  testthat::skip_on_cran()
  skip_if_null_sim_deps()

  positional <- simulateNullDataset(make_null_template(), 15, 2)
  named <- simulateNullDataset(
    make_null_template(),
    tree_tip_count = 15,
    seed = 2
  )
  positional$user_input <- named$user_input <- NULL

  testthat::expect_identical(positional, named)
})

test_that("simulateNullDataset returns response-only trait blocks for calibration templates", {
  testthat::skip_on_cran()
  skip_if_null_sim_deps()

  tmpl <- make_formula_null_template()
  sim <- simulateNullDataset(tmpl, tree_tip_count = 16, seed = 3)

  testthat::expect_equal(colnames(sim$trait_data), c("y1", "y2"))
  testthat::expect_false("mass" %in% colnames(sim$trait_data))
  testthat::expect_equal(dim(sim$simulatedData), c(16, 2))
})

test_that("simulateNullDataset validates tree_tip_count", {
  testthat::skip_on_cran()
  skip_if_null_sim_deps()

  tmpl <- make_null_template()
  testthat::expect_error(
    simulateNullDataset(tmpl, tree_tip_count = 3.5),
    "tree_tip_count"
  )
  testthat::expect_error(
    simulateNullDataset(tmpl, tree_tip_count = 100),
    "cannot exceed"
  )
})

test_that("simulateNullDataset validates template and argument types", {
  testthat::skip_on_cran()
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
    simulateNullDataset(tmpl, seed = 1.5),
    "seed"
  )
  testthat::expect_error(
    simulateNullDataset(tmpl, unexpected_arg = NA),
    "unused argument"
  )
  testthat::expect_error(
    simulateNullDataset(tmpl, simulation_generator = "unknown"),
    "simulation_generator"
  )
  testthat::expect_error(
    simulateNullDataset(tmpl, covariance_df = 2),
    "at least the number of response traits"
  )
})

test_that("simulateNullDataset restores the caller RNG state when seeded", {
  testthat::skip_on_cran()
  skip_if_null_sim_deps()

  tmpl <- make_null_template()
  expect_global_rng_restored(
    simulateNullDataset(tmpl, tree_tip_count = 15, seed = 2)
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

test_that("simulateNullDataset returns response-only data.frames for compatibility", {
  skip_if_null_sim_deps()

  sim_no_predictor <- simulateNullDataset(make_manual_null_template_df(FALSE), seed = 15)
  sim_with_predictor <- simulateNullDataset(make_manual_null_template_df(TRUE), seed = 16)

  testthat::expect_true(is.data.frame(sim_no_predictor$trait_data))
  testthat::expect_true(is.data.frame(sim_with_predictor$trait_data))
  testthat::expect_false("mass" %in% colnames(sim_with_predictor$trait_data))
  testthat::expect_equal(colnames(sim_with_predictor$trait_data), c("y1", "y2"))
})

test_that("simulateNullDataset errors when covariance generation never succeeds", {
  skip_if_null_sim_deps()

  tmpl <- make_zero_cov_null_template()

  testthat::expect_error(
    simulateNullDataset(tmpl, seed = 17),
    "positive-definite covariance matrix"
  )
})
