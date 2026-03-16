testthat::skip_on_cran()

skip_if_shift_sim_deps <- function() {
  testthat::skip_if_not_installed("ape")
  testthat::skip_if_not_installed("phytools")
  testthat::skip_if_not_installed("mvMORPH")
  testthat::skip_if_not_installed("RRphylo")
}

make_shift_template <- function() {
  set.seed(20)
  tr <- ape::rtree(28)
  X <- matrix(rnorm(28 * 3), ncol = 3)
  rownames(X) <- tr$tip.label
  createSimulationTemplate(
    baseline_tree = tr,
    trait_data = X,
    formula = "trait_data ~ 1",
    method = "LL"
  )
}

make_formula_shift_template <- function() {
  set.seed(21)
  tr <- ape::rtree(28)
  X <- cbind(y1 = rnorm(28), y2 = rnorm(28), mass = rnorm(28))
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

make_single_trait_shift_template <- function() {
  set.seed(22)
  tr <- ape::rtree(18)
  X <- matrix(rnorm(18), ncol = 1, dimnames = list(tr$tip.label, "y1"))
  structure(
    list(
      baseline_tree = tr,
      trait_data = X,
      fitted_values = matrix(0, nrow = 18, ncol = 1, dimnames = list(tr$tip.label, "y1")),
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

make_manual_shift_template_df <- function(with_predictor = FALSE) {
  set.seed(23)
  tr <- ape::rtree(16)
  if (with_predictor) {
    dat <- data.frame(
      y1 = rnorm(16),
      y2 = rnorm(16),
      mass = rnorm(16),
      row.names = tr$tip.label
    )
    fitted <- matrix(0, nrow = 16, ncol = 2, dimnames = list(tr$tip.label, c("y1", "y2")))
    response_columns <- 1:2
    predictor_columns <- 3L
    response_names <- c("y1", "y2")
    n_response_traits <- 2L
  } else {
    dat <- data.frame(
      y1 = rnorm(16),
      y2 = rnorm(16),
      row.names = tr$tip.label
    )
    fitted <- matrix(0, nrow = 16, ncol = 2, dimnames = list(tr$tip.label, c("y1", "y2")))
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

make_zero_cov_shift_template <- function() {
  set.seed(33)
  tr <- ape::rtree(16)
  X <- matrix(0, nrow = 16, ncol = 2, dimnames = list(tr$tip.label, c("y1", "y2")))
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

test_that("simulateShiftedDataset returns a proportional shifted replicate", {
  skip_if_shift_sim_deps()

  tmpl <- make_shift_template()
  sim <- simulateShiftedDataset(
    tmpl,
    tree_tip_count = 24,
    num_shifts = 2,
    min_shift_tips = 3,
    max_shift_tips = 8,
    scale_mode = "proportional",
    buffer = 2,
    seed = 4
  )

  testthat::expect_s3_class(sim, "bifrost_simulation_replicate_shifted")
  testthat::expect_length(sim$shiftNodes, 2)
  testthat::expect_length(sim$VCVs, 3)
  testthat::expect_identical(sim$generating_scenario, "proportional")
  testthat::expect_identical(rownames(sim$trait_data), sim$paintedTree$tip.label)
  testthat::expect_true(all(vapply(sim$VCVs, function(x) {
    all(eigen(x, symmetric = TRUE)$values > 0)
  }, logical(1))))

  clade_sizes <- vapply(sim$shiftNodes, function(node) {
    sum(getDescendants(sim$paintedTree, node) <= ape::Ntip(sim$paintedTree))
  }, integer(1))
  testthat::expect_true(all(clade_sizes >= 3 & clade_sizes <= 8))

  node_distance <- RRphylo::distNodes(
    ape::reorder.phylo(ape::as.phylo(sim$paintedTree), order = "cladewise"),
    node = sim$shiftNodes,
    clus = 0
  )$node
  testthat::expect_gte(node_distance, 2)
})

test_that("simulateShiftedDataset supports the correlation scenario and preserves predictors", {
  skip_if_shift_sim_deps()

  tmpl <- make_formula_shift_template()
  sim <- simulateShiftedDataset(
    tmpl,
    tree_tip_count = 24,
    num_shifts = 2,
    min_shift_tips = 3,
    max_shift_tips = 7,
    scale_mode = "correlation",
    seed = 5
  )

  testthat::expect_identical(sim$generating_scenario, "correlation")
  testthat::expect_identical(
    sim$trait_data[, 3],
    tmpl$trait_data[sim$sampled_tree$tip.label, 3]
  )
  testthat::expect_equal(dim(sim$simulatedData), c(24, 2))
})

test_that("simulateShiftedDataset requires predictors for formula templates", {
  skip_if_shift_sim_deps()

  tmpl <- make_formula_shift_template()

  testthat::expect_error(
    simulateShiftedDataset(
      tmpl,
      tree_tip_count = 24,
      num_shifts = 2,
      min_shift_tips = 3,
      max_shift_tips = 7,
      scale_mode = "correlation",
      preserve_predictors = FALSE,
      seed = 5
    ),
    "preserve_predictors must be TRUE"
  )
})

test_that("simulateShiftedDataset validates argument types and edge cases", {
  skip_if_shift_sim_deps()

  tmpl <- make_shift_template()
  one_trait <- make_single_trait_shift_template()

  testthat::expect_error(simulateShiftedDataset(list(), num_shifts = 1, min_shift_tips = 1, max_shift_tips = 2), "bifrost_simulation_template")
  testthat::expect_error(simulateShiftedDataset(tmpl, num_shifts = 0, min_shift_tips = 1, max_shift_tips = 2), "num_shifts")
  testthat::expect_error(simulateShiftedDataset(tmpl, num_shifts = 1, min_shift_tips = 0, max_shift_tips = 2), "min_shift_tips")
  testthat::expect_error(simulateShiftedDataset(tmpl, num_shifts = 1, min_shift_tips = 2, max_shift_tips = 1), "max_shift_tips")
  testthat::expect_error(simulateShiftedDataset(tmpl, num_shifts = 1, min_shift_tips = 1, max_shift_tips = 2, buffer = -1), "buffer")
  testthat::expect_error(simulateShiftedDataset(tmpl, num_shifts = 1, min_shift_tips = 1, max_shift_tips = 2, scale_factor_range = c(1, 1)), "scale_factor_range")
  testthat::expect_error(simulateShiftedDataset(tmpl, num_shifts = 1, min_shift_tips = 1, max_shift_tips = 2, exclude_range = c(2, 1)), "exclude_range")
  testthat::expect_error(simulateShiftedDataset(tmpl, num_shifts = 1, min_shift_tips = 1, max_shift_tips = 2, exclude_range = c(0, 2)), "strictly inside")
  testthat::expect_error(simulateShiftedDataset(one_trait, num_shifts = 1, min_shift_tips = 1, max_shift_tips = 2, scale_mode = "correlation"), "requires at least two response traits")
  testthat::expect_error(simulateShiftedDataset(tmpl, tree_tip_count = 1, num_shifts = 1, min_shift_tips = 1, max_shift_tips = 2), "single integer >= 2")
  testthat::expect_error(simulateShiftedDataset(tmpl, tree_tip_count = 100, num_shifts = 1, min_shift_tips = 1, max_shift_tips = 2), "cannot exceed")
  testthat::expect_error(simulateShiftedDataset(tmpl, num_shifts = 1, min_shift_tips = 1, max_shift_tips = 2, preserve_predictors = NA), "TRUE or FALSE")
})

test_that("simulateShiftedDataset supports full-tree and single-trait proportional runs", {
  skip_if_shift_sim_deps()

  tmpl <- make_single_trait_shift_template()
  sim <- simulateShiftedDataset(
    tmpl,
    num_shifts = 1,
    min_shift_tips = 2,
    max_shift_tips = 6,
    scale_mode = "proportional",
    buffer = 1,
    seed = 24
  )

  testthat::expect_equal(length(sim$paintedTree$tip.label), tmpl$n_tips)
  testthat::expect_equal(dim(sim$simulatedData), c(tmpl$n_tips, 1))
})

test_that("simulateShiftedDataset preserves data.frame outputs when requested", {
  skip_if_shift_sim_deps()

  sim_no_predictor <- simulateShiftedDataset(
    make_manual_shift_template_df(FALSE),
    num_shifts = 1,
    min_shift_tips = 2,
    max_shift_tips = 6,
    scale_mode = "proportional",
    seed = 25
  )
  template_with_predictor <- make_manual_shift_template_df(TRUE)
  sim_with_predictor <- simulateShiftedDataset(
    template_with_predictor,
    num_shifts = 1,
    min_shift_tips = 2,
    max_shift_tips = 6,
    scale_mode = "proportional",
    seed = 26
  )

  testthat::expect_true(is.data.frame(sim_no_predictor$trait_data))
  testthat::expect_true(is.data.frame(sim_with_predictor$trait_data))
  testthat::expect_identical(
    sim_with_predictor$trait_data$mass,
    template_with_predictor$trait_data[sim_with_predictor$sampled_tree$tip.label, "mass"]
  )
})

test_that("simulateShiftedDataset fails cleanly when no valid shift configuration exists", {
  skip_if_shift_sim_deps()

  tmpl <- make_shift_template()

  testthat::expect_error(
    simulateShiftedDataset(
      tmpl,
      num_shifts = 2,
      min_shift_tips = 50,
      max_shift_tips = 60,
      scale_mode = "proportional",
      seed = 27
    ),
    "Failed to place all requested shifts"
  )
})

test_that("simulateShiftedDataset covers placement failure after candidate removal", {
  skip_if_shift_sim_deps()

  tree_text <- "(((t1:1,t2:1):1,t3:1):1,t4:1);"
  pectinate <- ape::read.tree(text = tree_text)
  X <- matrix(rnorm(4 * 2), ncol = 2)
  rownames(X) <- pectinate$tip.label
  tmpl <- createSimulationTemplate(pectinate, X, formula = "trait_data ~ 1", method = "LL")

  testthat::expect_error(
    simulateShiftedDataset(
      tmpl,
      num_shifts = 2,
      min_shift_tips = 2,
      max_shift_tips = 4,
      scale_mode = "proportional",
      seed = 28
    ),
    "Failed to place all requested shifts"
  )
})

test_that("simulateShiftedDataset handles singleton candidate sets without resampling 1:n", {
  skip_if_shift_sim_deps()

  tree_text <- "((((t1:1,t2:1):1,t3:2):1,t4:3):1,t5:4);"
  pectinate <- ape::read.tree(text = tree_text)
  X <- matrix(rnorm(5 * 2), ncol = 2)
  rownames(X) <- pectinate$tip.label
  tmpl <- createSimulationTemplate(pectinate, X, formula = "trait_data ~ 1", method = "LL")

  sim <- simulateShiftedDataset(
    tmpl,
    num_shifts = 1,
    min_shift_tips = 4,
    max_shift_tips = 4,
    scale_mode = "proportional",
    seed = 29
  )

  testthat::expect_identical(sim$shiftNodes, 7L)
  testthat::expect_true(all(rownames(sim$trait_data) == sim$paintedTree$tip.label))
})

test_that("simulateShiftedDataset covers buffer failures", {
  skip_if_shift_sim_deps()

  testthat::local_mocked_bindings(
    distNodes = function(...) list(node = 0),
    .package = "RRphylo"
  )

  testthat::expect_error(
    simulateShiftedDataset(
      make_shift_template(),
      num_shifts = 2,
      min_shift_tips = 3,
      max_shift_tips = 8,
      scale_mode = "proportional",
      buffer = 3,
      seed = 29
    ),
    "Failed to place all requested shifts"
  )
})

test_that("simulateShiftedDataset covers ancestral covariance failures", {
  skip_if_shift_sim_deps()

  testthat::expect_error(
    simulateShiftedDataset(
      make_zero_cov_shift_template(),
      num_shifts = 1,
      min_shift_tips = 2,
      max_shift_tips = 5,
      scale_mode = "proportional",
      seed = 30
    ),
    "Failed to place all requested shifts"
  )
})

test_that("simulateShiftedDataset covers derived covariance failures", {
  skip_if_shift_sim_deps()

  testthat::expect_error(
    simulateShiftedDataset(
      make_shift_template(),
      num_shifts = 1,
      min_shift_tips = 3,
      max_shift_tips = 8,
      scale_mode = "proportional",
      scale_factor_range = c(-2, -0.1),
      exclude_range = c(-1.5, -0.5),
      seed = 31
    ),
    "Failed to place all requested shifts"
  )
})

test_that("simulateShiftedDataset covers state validation failures", {
  skip_if_shift_sim_deps()

  testthat::local_mocked_bindings(
    getStates = function(...) "ancestral",
    .package = "phytools"
  )

  testthat::expect_error(
    simulateShiftedDataset(
      make_shift_template(),
      num_shifts = 1,
      min_shift_tips = 3,
      max_shift_tips = 8,
      scale_mode = "proportional",
      seed = 32
    ),
    "Failed to place all requested shifts"
  )
})
