make_covariance_helper_template <- function() {
  sigma <- matrix(
    c(
      4.0, 1.2, 0.4,
      1.2, 2.0, 0.3,
      0.4, 0.3, 1.0
    ),
    nrow = 3L,
    byrow = TRUE,
    dimnames = list(paste0("trait", 1:3), paste0("trait", 1:3))
  )

  list(
    residual_covariance = sigma,
    residual_df = 40L,
    n_response_traits = 3L,
    variance_mean = mean(diag(sigma)),
    variance_sd = stats::sd(diag(sigma)),
    covariance_mean = mean(sigma[lower.tri(sigma)]),
    covariance_sd = stats::sd(sigma[lower.tri(sigma)])
  )
}

test_that("simulation generator names and public defaults are explicit", {
  testthat::expect_identical(.simulation_check_generator("original"), "original")
  testthat::expect_identical(.simulation_check_generator("empirical"), "empirical")
  testthat::expect_error(
    .simulation_check_generator("corrected"),
    "simulation_generator.*original.*empirical"
  )
  testthat::expect_identical(
    eval(formals(simulateNullDataset)$simulation_generator),
    c("original", "empirical")
  )
  testthat::expect_identical(
    eval(formals(simulateShiftedDataset)$simulation_generator),
    c("original", "empirical")
  )
  testthat::expect_identical(.simulation_generator_from_options(list()), "original")
  testthat::expect_error(
    .simulation_generator_from_options(list(
      simulation_generator = c("original", "empirical")
    )),
    "simulation_generator.*original.*empirical"
  )
})

test_that("simulation generators retain legacy positional argument order", {
  testthat::expect_identical(
    names(formals(simulateNullDataset))[seq_len(3L)],
    c("template", "tree_tip_count", "seed")
  )
  testthat::expect_identical(
    names(formals(simulateShiftedDataset))[seq_len(10L)],
    c(
      "template", "tree_tip_count", "num_shifts", "min_shift_tips",
      "max_shift_tips", "scale_mode", "scale_factor_range", "exclude_range",
      "buffer", "seed"
    )
  )
})

test_that("internal covariance draws default to the original generator", {
  template <- make_covariance_helper_template()

  set.seed(102)
  default <- .simulation_draw_covariance(template)
  set.seed(102)
  explicit <- .simulation_draw_covariance(template, "original")

  testthat::expect_identical(default, explicit)
})

test_that("empirical covariance draws are reproducible and positive definite", {
  template <- make_covariance_helper_template()

  set.seed(100)
  first <- .simulation_draw_covariance(template, "empirical", NULL)
  set.seed(100)
  second <- .simulation_draw_covariance(template, "empirical", NULL)

  testthat::expect_identical(first, second)
  testthat::expect_identical(first$covariance_df, 40L)
  testthat::expect_equal(first$sigma, t(first$sigma), tolerance = 1e-14)
  testthat::expect_true(all(eigen(first$sigma, symmetric = TRUE)$values > 0))
})

test_that("empirical covariance draws are centered on the empirical covariance", {
  template <- make_covariance_helper_template()
  set.seed(101)
  draws <- replicate(
    3000L,
    .simulation_draw_covariance(template, "empirical", NULL)$sigma,
    simplify = "array"
  )
  draw_mean <- apply(draws, c(1L, 2L), mean)

  testthat::expect_lt(
    max(abs(draw_mean - template$residual_covariance)) /
      max(abs(template$residual_covariance)),
    0.08
  )
})

test_that("empirical covariance degrees of freedom are validated", {
  template <- make_covariance_helper_template()

  testthat::expect_identical(
    .simulation_draw_covariance(template, "empirical", 12)$covariance_df,
    12L
  )
  testthat::expect_error(
    .simulation_draw_covariance(template, "empirical", 2),
    "at least the number of response traits"
  )
  testthat::expect_error(
    .simulation_draw_covariance(template, "unknown", NULL),
    "simulation_generator"
  )
})

test_that("integration transforms preserve marginal variances and trait units", {
  sigma <- make_covariance_helper_template()$residual_covariance
  units <- diag(c(10, 0.1, 3))

  transformed <- .simulation_transform_integration(sigma, power = 1.2)
  transformed_units <- .simulation_transform_integration(
    units %*% sigma %*% units,
    power = 1.2
  )

  testthat::expect_equal(diag(transformed$sigma), diag(sigma), tolerance = 1e-12)
  testthat::expect_equal(
    transformed_units$sigma,
    units %*% transformed$sigma %*% units,
    tolerance = 1e-10
  )
  testthat::expect_true(all(eigen(transformed$sigma, symmetric = TRUE)$values > 0))
})

test_that("integration powers move effective dimensionality in both directions", {
  sigma <- make_covariance_helper_template()$residual_covariance
  reduced <- .simulation_transform_integration(sigma, power = 0.6)
  unchanged <- .simulation_transform_integration(sigma, power = 1)
  increased <- .simulation_transform_integration(sigma, power = 1.2)

  testthat::expect_gt(
    reduced$diagnostics$effective_dimensionality_after,
    reduced$diagnostics$effective_dimensionality_before
  )
  testthat::expect_lt(
    increased$diagnostics$effective_dimensionality_after,
    increased$diagnostics$effective_dimensionality_before
  )
  testthat::expect_lt(
    reduced$diagnostics$mean_absolute_correlation_after,
    reduced$diagnostics$mean_absolute_correlation_before
  )
  testthat::expect_gt(
    increased$diagnostics$mean_absolute_correlation_after,
    increased$diagnostics$mean_absolute_correlation_before
  )
  testthat::expect_equal(unchanged$sigma, sigma, tolerance = 1e-12)
})

test_that("integration transform rejects invalid covariance inputs and powers", {
  sigma <- make_covariance_helper_template()$residual_covariance

  testthat::expect_error(
    .simulation_transform_integration(sigma, power = 0),
    "power"
  )
  testthat::expect_error(
    .simulation_transform_integration(sigma[, 1:2], power = 1.2),
    "square"
  )
  testthat::expect_error(
    .simulation_transform_integration(diag(c(1, 0, 1)), power = 1.2),
    "positive diagonal"
  )
})
