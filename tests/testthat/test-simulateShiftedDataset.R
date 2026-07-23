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
  testthat::skip_on_cran()
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
  testthat::expect_false(any(vapply(seq_along(sim$shiftNodes), function(i) {
    any(vapply(seq_along(sim$shiftNodes), function(j) {
      if (i == j) {
        return(FALSE)
      }
      sim$shiftNodes[j] %in% getDescendants(sim$paintedTree, sim$shiftNodes[i])
    }, logical(1)))
  }, logical(1))))

  node_distance <- RRphylo::distNodes(
    ape::reorder.phylo(ape::as.phylo(sim$paintedTree), order = "cladewise"),
    node = sim$shiftNodes,
    clus = 0
  )$node
  testthat::expect_gte(node_distance, 2)
})

test_that("simulateShiftedDataset returns response-only trait blocks in the correlation scenario", {
  testthat::skip_on_cran()
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
  testthat::expect_equal(colnames(sim$trait_data), c("y1", "y2"))
  testthat::expect_false("mass" %in% colnames(sim$trait_data))
  testthat::expect_equal(dim(sim$simulatedData), c(24, 2))
})

test_that("empirical correlation shifts impose the integration-rate trade-off", {
  testthat::skip_on_cran()
  skip_if_shift_sim_deps()

  tmpl <- make_shift_template()
  sim <- simulateShiftedDataset(
    tmpl,
    tree_tip_count = 24,
    num_shifts = 2,
    min_shift_tips = 3,
    max_shift_tips = 7,
    scale_mode = "correlation",
    simulation_generator = "empirical",
    seed = 1
  )

  ancestral_sigma <- sim$VCVs[["ancestral"]]
  derived_vcvs <- sim$VCVs[setdiff(names(sim$VCVs), "ancestral")]

  testthat::expect_identical(sim$simulation_generator, "empirical")
  testthat::expect_false("simulation_design" %in% names(sim))
  testthat::expect_identical(sim$covariance_df, tmpl$residual_df)
  testthat::expect_length(sim$sampledIntegrationPowers, 2L)
  testthat::expect_true(all(
    sim$sampledIntegrationPowers < 0.8 |
      sim$sampledIntegrationPowers > 1.1
  ))
  testthat::expect_equal(
    sim$sampledScaleFactors,
    sim$sampledIntegrationPowers,
    tolerance = 0
  )
  testthat::expect_equal(
    sim$sampledVarianceScaleFactors,
    1 / sim$sampledIntegrationPowers,
    tolerance = 1e-14
  )
  for (i in seq_along(derived_vcvs)) {
    power <- sim$sampledIntegrationPowers[i]
    testthat::expect_equal(
      diag(derived_vcvs[[i]]),
      diag(ancestral_sigma) / power,
      tolerance = 1e-12
    )
    diagnostics <- sim$covariance_diagnostics[[i + 1L]]
    if (power < 1) {
      testthat::expect_gt(
        diagnostics$effective_dimensionality_after,
        diagnostics$effective_dimensionality_before
      )
      testthat::expect_gt(sim$sampledVarianceScaleFactors[i], 1)
    } else {
      testthat::expect_lt(
        diagnostics$effective_dimensionality_after,
        diagnostics$effective_dimensionality_before
      )
      testthat::expect_lt(sim$sampledVarianceScaleFactors[i], 1)
    }
  }
  testthat::expect_named(
    sim$covariance_diagnostics,
    names(sim$VCVs)
  )
})

test_that("the direct simulator rejects generator vectors before RNG setup", {
  skip_if_shift_sim_deps()
  rng_setup_called <- FALSE
  testthat::local_mocked_bindings(
    .simulation_set_seed = function(...) {
      rng_setup_called <<- TRUE
      stats::runif(1)
      stop("RNG setup reached", call. = FALSE)
    }
  )
  tmpl <- make_shift_template()
  set.seed(20260709)
  before <- .Random.seed
  testthat::expect_error(
    simulateShiftedDataset(
      tmpl,
      num_shifts = 1,
      min_shift_tips = 3,
      max_shift_tips = 8,
      seed = 4,
      simulation_generator = c("original", "empirical")
    ),
    "simulation_generator.*original.*empirical"
  )
  testthat::expect_false(rng_setup_called)
  testthat::expect_identical(.Random.seed, before)
})

test_that("simulateShiftedDataset defaults to the original correlation matrices", {
  testthat::skip_on_cran()
  skip_if_shift_sim_deps()

  expected_ancestral <- matrix(
    c(
      1.32011756719237, 0.0544924885729728, -0.000401143385120667,
      0.0544924885729728, 1.40986503043954, -0.000423077927204905,
      -0.000401143385120667, -0.000423077927204905, 1.21737523566526
    ),
    nrow = 3L
  )
  expected_derived_1 <- matrix(
    c(
      8.54111332430542, 0.0544924885729728, -0.000401143385120667,
      0.0544924885729728, 9.1217761934416, -0.000423077927204905,
      -0.000401143385120667, -0.000423077927204905, 7.87637412335476
    ),
    nrow = 3L
  )
  expected_derived_2 <- matrix(
    c(
      2.95610148716854, 0.0544924885729728, -0.000401143385120667,
      0.0544924885729728, 3.15707041309443, -0.000423077927204905,
      -0.000401143385120667, -0.000423077927204905, 2.72603352460942
    ),
    nrow = 3L
  )
  sim <- simulateShiftedDataset(
    make_shift_template(),
    tree_tip_count = 24,
    num_shifts = 2,
    min_shift_tips = 3,
    max_shift_tips = 7,
    scale_mode = "correlation",
    seed = 5
  )

  testthat::expect_identical(sim$simulation_generator, "original")
  testthat::expect_identical(sim$shiftNodes, c(26L, 37L))
  testthat::expect_equal(
    sim$sampledScaleFactors,
    c(0.154560361988842, 0.446573831420392),
    tolerance = 1e-14
  )
  testthat::expect_equal(sim$VCVs$ancestral, expected_ancestral, tolerance = 1e-14)
  testthat::expect_equal(sim$VCVs$derived_1, expected_derived_1, tolerance = 1e-14)
  testthat::expect_equal(sim$VCVs$derived_2, expected_derived_2, tolerance = 1e-14)
})

test_that("simulateShiftedDataset supports its legacy positional buffer and seed", {
  testthat::skip_on_cran()
  skip_if_shift_sim_deps()

  positional <- simulateShiftedDataset(
    make_shift_template(), 24, 2, 3, 7, "correlation",
    c(0.1, 2), c(0.5, 1.5), 2, 5
  )
  named <- simulateShiftedDataset(
    make_shift_template(),
    tree_tip_count = 24,
    num_shifts = 2,
    min_shift_tips = 3,
    max_shift_tips = 7,
    scale_mode = "correlation",
    scale_factor_range = c(0.1, 2),
    exclude_range = c(0.5, 1.5),
    buffer = 2,
    seed = 5
  )
  positional$user_input <- named$user_input <- NULL

  testthat::expect_identical(positional, named)
})

test_that("simulateShiftedDataset correlation mode is not just proportional rescaling", {
  testthat::skip_on_cran()
  skip_if_shift_sim_deps()

  tmpl <- make_shift_template()
  sim <- simulateShiftedDataset(
    tmpl,
    tree_tip_count = 24,
    num_shifts = 2,
    min_shift_tips = 3,
    max_shift_tips = 7,
    scale_mode = "correlation",
    simulation_generator = "empirical",
    seed = 5
  )

  ancestral_sigma <- sim$VCVs[["ancestral"]]
  derived_sigma <- sim$VCVs[[setdiff(names(sim$VCVs), "ancestral")[[1L]]]]
  scaling_factor <- derived_sigma[1, 1] / ancestral_sigma[1, 1]

  testthat::expect_gt(max(abs(derived_sigma - (scaling_factor * ancestral_sigma))), 1e-8)
  testthat::expect_gt(
    max(abs(stats::cov2cor(derived_sigma) - stats::cov2cor(ancestral_sigma))[lower.tri(ancestral_sigma)]),
    1e-8
  )
})

test_that("simulateShiftedDataset validates argument types and edge cases", {
  testthat::skip_on_cran()
  skip_if_shift_sim_deps()

  tmpl <- make_shift_template()
  one_trait <- make_single_trait_shift_template()

  testthat::expect_error(simulateShiftedDataset(list(), num_shifts = 1, min_shift_tips = 1, max_shift_tips = 2), "bifrost_simulation_template")
  testthat::expect_error(simulateShiftedDataset(tmpl, num_shifts = 0, min_shift_tips = 1, max_shift_tips = 2), "num_shifts")
  testthat::expect_error(simulateShiftedDataset(tmpl, num_shifts = 1.5, min_shift_tips = 1, max_shift_tips = 2), "num_shifts")
  testthat::expect_error(simulateShiftedDataset(tmpl, num_shifts = 1, min_shift_tips = 0, max_shift_tips = 2), "min_shift_tips")
  testthat::expect_error(simulateShiftedDataset(tmpl, num_shifts = 1, min_shift_tips = 1.5, max_shift_tips = 2), "min_shift_tips")
  testthat::expect_error(simulateShiftedDataset(tmpl, num_shifts = 1, min_shift_tips = 2, max_shift_tips = 1), "max_shift_tips")
  testthat::expect_error(simulateShiftedDataset(tmpl, num_shifts = 1, min_shift_tips = 1, max_shift_tips = 2.5), "max_shift_tips")
  testthat::expect_error(simulateShiftedDataset(tmpl, num_shifts = 1, min_shift_tips = 1, max_shift_tips = 2, buffer = -1), "buffer")
  testthat::expect_error(simulateShiftedDataset(tmpl, num_shifts = 1, min_shift_tips = 1, max_shift_tips = 2, buffer = 1.5), "buffer")
  testthat::expect_error(simulateShiftedDataset(tmpl, num_shifts = 1, min_shift_tips = 1, max_shift_tips = 2, scale_factor_range = c(1, 1)), "scale_factor_range")
  testthat::expect_error(simulateShiftedDataset(tmpl, num_shifts = 1, min_shift_tips = 1, max_shift_tips = 2, scale_factor_range = c(0, 2)), "finite positive")
  testthat::expect_error(simulateShiftedDataset(tmpl, num_shifts = 1, min_shift_tips = 1, max_shift_tips = 2, scale_factor_range = c(0.1, Inf)), "finite positive")
  testthat::expect_error(simulateShiftedDataset(tmpl, num_shifts = 1, min_shift_tips = 1, max_shift_tips = 2, exclude_range = c(2, 1)), "exclude_range")
  testthat::expect_error(simulateShiftedDataset(tmpl, num_shifts = 1, min_shift_tips = 1, max_shift_tips = 2, exclude_range = c(0, 2)), "strictly inside")
  testthat::expect_error(simulateShiftedDataset(one_trait, num_shifts = 1, min_shift_tips = 1, max_shift_tips = 2, scale_mode = "correlation"), "requires at least two response traits")
  testthat::expect_error(simulateShiftedDataset(tmpl, tree_tip_count = 1, num_shifts = 1, min_shift_tips = 1, max_shift_tips = 2), "single integer >= 2")
  testthat::expect_error(simulateShiftedDataset(tmpl, tree_tip_count = 3.5, num_shifts = 1, min_shift_tips = 1, max_shift_tips = 2), "tree_tip_count")
  testthat::expect_error(simulateShiftedDataset(tmpl, tree_tip_count = 100, num_shifts = 1, min_shift_tips = 1, max_shift_tips = 2), "cannot exceed")
  testthat::expect_error(simulateShiftedDataset(tmpl, num_shifts = 1, min_shift_tips = 1, max_shift_tips = 2, seed = 1.5), "seed")
  testthat::expect_error(simulateShiftedDataset(tmpl, num_shifts = 1, min_shift_tips = 1, max_shift_tips = 2, unexpected_arg = NA), "unused argument")
  testthat::expect_error(simulateShiftedDataset(tmpl, num_shifts = 1, min_shift_tips = 1, max_shift_tips = 2, simulation_generator = "unknown"), "simulation_generator")
  testthat::expect_error(simulateShiftedDataset(tmpl, num_shifts = 1, min_shift_tips = 1, max_shift_tips = 2, covariance_df = 2), "at least the number of response traits")
  testthat::expect_error(simulateShiftedDataset(tmpl, num_shifts = 1, min_shift_tips = 1, max_shift_tips = 2, integration_power_range = c(1, 0.5)), "integration_power_range")
  testthat::expect_error(simulateShiftedDataset(tmpl, num_shifts = 1, min_shift_tips = 1, max_shift_tips = 2, integration_exclude_range = c(0.2, 0.9)), "contains 1")
  testthat::expect_error(simulateShiftedDataset(tmpl, num_shifts = 1, min_shift_tips = 1, max_shift_tips = 2, integration_exclude_range = c(0.2, 1.3)), "strictly inside")
  testthat::expect_error(simulateShiftedDataset(tmpl, num_shifts = 1, min_shift_tips = 1, max_shift_tips = 2, eigen_floor = 1), "eigen_floor")
})

test_that("simulateShiftedDataset restores the caller RNG state when seeded", {
  testthat::skip_on_cran()
  skip_if_shift_sim_deps()

  tmpl <- make_shift_template()
  expect_global_rng_restored(
    simulateShiftedDataset(
      tmpl,
      tree_tip_count = 24,
      num_shifts = 1,
      min_shift_tips = 3,
      max_shift_tips = 8,
      scale_mode = "proportional",
      buffer = 2,
      seed = 4
    )
  )
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

test_that("simulateShiftedDataset returns response-only data.frames for compatibility", {
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
  testthat::expect_false("mass" %in% colnames(sim_with_predictor$trait_data))
  testthat::expect_equal(colnames(sim_with_predictor$trait_data), c("y1", "y2"))
})

test_that("simulateShiftedDataset fails cleanly when no valid shift configuration exists", {
  testthat::skip_on_cran()
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

test_that("simulateShiftedDataset fails cleanly when candidate removal prevents placement", {
  testthat::skip_on_cran()
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
  testthat::skip_on_cran()
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

test_that("simulateShiftedDataset fails cleanly when buffer constraints remove all candidates", {
  testthat::skip_on_cran()
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

test_that("simulateShiftedDataset fails cleanly when ancestral covariance generation fails", {
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

test_that("simulateShiftedDataset fails cleanly when derived covariance generation fails", {
  skip_if_shift_sim_deps()

  tmpl <- make_manual_shift_template_df()
  tmpl$variance_mean <- 0.01
  tmpl$variance_sd <- 0
  tmpl$covariance_mean <- 1
  tmpl$covariance_sd <- 0

  testthat::expect_error(
    simulateShiftedDataset(
      tmpl,
      num_shifts = 1,
      min_shift_tips = 3,
      max_shift_tips = 8,
      scale_mode = "correlation",
      scale_factor_range = c(3, 4),
      exclude_range = c(3.2, 3.8),
      seed = 31
    ),
    "Failed to place all requested shifts"
  )
})

test_that("simulateShiftedDataset validates painted-tree states before simulation", {
  testthat::skip_on_cran()
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

test_that("simulateShiftedDataset evaluates each deterministic derived covariance once", {
  skip_if_shift_sim_deps()

  transform_count <- 0L
  real_transform <- .simulation_transform_integration
  testthat::local_mocked_bindings(
    .simulation_transform_integration = function(...) {
      transform_count <<- transform_count + 1L
      transformed <- real_transform(...)
      transformed$sigma[,] <- 0
      transformed
    }
  )

  testthat::expect_error(
    simulateShiftedDataset(
      make_shift_template(),
      num_shifts = 1,
      min_shift_tips = 3,
      max_shift_tips = 8,
      scale_mode = "correlation",
      simulation_generator = "empirical",
      seed = 34
    ),
    "Failed to place all requested shifts"
  )
  testthat::expect_identical(transform_count, 100L)
})
