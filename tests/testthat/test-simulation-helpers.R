test_that("simulation integer validators accept and reject scalar/vector controls", {
  testthat::expect_identical(.simulation_check_integer_scalar(NULL, "x", allow_null = TRUE), NULL)
  testthat::expect_identical(.simulation_check_integer_scalar(3, "x", minimum = 1, maximum = 5), 3L)
  testthat::expect_error(
    .simulation_check_integer_scalar(6, "x", maximum = 5),
    "`x` must be a single finite integer"
  )

  testthat::expect_identical(.simulation_check_integer_vector(c(1, 3), "x", minimum = 1), c(1L, 3L))
  testthat::expect_error(
    .simulation_check_integer_vector(c(1, NA), "x"),
    "`x` must be a finite integer vector"
  )
  testthat::expect_error(
    .simulation_check_integer_vector(0, "x", minimum = 1),
    "`x` must be a finite integer vector"
  )
})

test_that("simulation seed helpers restore caller RNG state", {
  old_seed <- if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
    get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  } else {
    NULL
  }
  on.exit({
    if (is.null(old_seed)) {
      if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
        rm(".Random.seed", envir = .GlobalEnv)
      }
    } else {
      assign(".Random.seed", old_seed, envir = .GlobalEnv)
    }
  }, add = TRUE)

  if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
    rm(".Random.seed", envir = .GlobalEnv)
  }
  testthat::expect_null(.simulation_restore_seed(NULL))
  testthat::expect_null(.simulation_set_seed(NULL))
  seed_state <- .simulation_set_seed(42)
  testthat::expect_false(seed_state$had_seed)
  testthat::expect_true(exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
  .simulation_restore_seed(seed_state)
  testthat::expect_false(exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))

  set.seed(100)
  existing_seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  seed_state <- .simulation_set_seed(101)
  .simulation_restore_seed(seed_state)
  testthat::expect_identical(get(".Random.seed", envir = .GlobalEnv, inherits = FALSE), existing_seed)

  kind_before <- RNGkind()
  seed_state <- .simulation_set_seed(102, kind = "L'Ecuyer-CMRG")
  testthat::expect_identical(RNGkind()[[1L]], "L'Ecuyer-CMRG")
  .simulation_restore_seed(seed_state)
  testthat::expect_identical(RNGkind(), kind_before)
  testthat::expect_identical(get(".Random.seed", envir = .GlobalEnv, inherits = FALSE), existing_seed)
})

test_that("simulation replicate seeds preserve unseeded execution", {
  seeds <- .simulation_replicate_seeds(3L, seeded = FALSE)

  testthat::expect_length(seeds, 3L)
  testthat::expect_true(all(vapply(seeds, is.null, logical(1))))
})

test_that("simulation seed runner restores absent and existing caller state on errors", {
  old_seed <- if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
    get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  } else {
    NULL
  }
  on.exit({
    if (is.null(old_seed)) {
      if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
        rm(".Random.seed", envir = .GlobalEnv)
      }
    } else {
      assign(".Random.seed", old_seed, envir = .GlobalEnv)
    }
  }, add = TRUE)

  runner <- .simulation_seed_runner()
  if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
    rm(".Random.seed", envir = .GlobalEnv)
  }
  testthat::expect_identical(runner(NULL, 42L), 42L)
  testthat::expect_false(exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
  existing_kind <- RNGkind()
  testthat::expect_type(runner(102, stats::runif(1L)), "double")
  testthat::expect_identical(RNGkind(), existing_kind)
  testthat::expect_false(exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
  testthat::expect_error(runner(103, stop("seeded failure")), "seeded failure")
  testthat::expect_identical(RNGkind(), existing_kind)
  testthat::expect_false(exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))

  set.seed(104)
  existing_seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  testthat::expect_type(runner(105, stats::runif(1L)), "double")
  testthat::expect_identical(
    get(".Random.seed", envir = .GlobalEnv, inherits = FALSE),
    existing_seed
  )
  testthat::expect_error(runner(105, stop("seeded failure")), "seeded failure")
  testthat::expect_identical(
    get(".Random.seed", envir = .GlobalEnv, inherits = FALSE),
    existing_seed
  )
})

test_that("simulation covariance validation characterizes malformed matrices", {
  testthat::expect_error(
    .simulation_validate_covariance(matrix("x", nrow = 1L)),
    "finite numeric values"
  )
  testthat::expect_error(
    .simulation_validate_covariance(matrix(c(1, 0.1, 0.2, 1), nrow = 2L)),
    "symmetric"
  )
  testthat::expect_error(
    .simulation_validate_covariance(matrix(c(1, 2, 2, 1), nrow = 2L)),
    "positive definite"
  )
})

test_that("empirical covariance draws use template fallbacks and validate dimensions", {
  template <- list(
    n_response_traits = 2L,
    residual_covariance = NULL,
    covariance_mean = 0,
    variance_mean = 1,
    residual_df = NULL,
    n_tips = 8L
  )

  withr::local_seed(106)
  draw <- .simulation_draw_covariance(template, simulation_generator = "empirical")
  testthat::expect_equal(dim(draw$sigma), c(2L, 2L))
  testthat::expect_identical(draw$covariance_df, 7L)
  testthat::expect_true(all(eigen(draw$sigma, symmetric = TRUE, only.values = TRUE)$values > 0))

  mismatched <- template
  mismatched$n_response_traits <- 3L
  mismatched$residual_covariance <- diag(2)
  mismatched$residual_df <- 4L
  testthat::expect_error(
    .simulation_draw_covariance(mismatched, simulation_generator = "empirical"),
    "dimensions must match the number of response traits"
  )
})

test_that("integration covariance transforms validate floors and preserve marginal variances", {
  sigma <- matrix(
    c(2, 0.6, 0.6, 1),
    nrow = 2L,
    dimnames = list(c("a", "b"), c("a", "b"))
  )

  testthat::expect_error(
    .simulation_transform_integration(sigma, power = 2, eigen_floor = 1),
    "strictly between zero and one"
  )

  identity_transform <- .simulation_transform_integration(sigma, power = 1)
  stronger_transform <- .simulation_transform_integration(sigma, power = 2)
  testthat::expect_identical(identity_transform$sigma, sigma)
  testthat::expect_equal(diag(stronger_transform$sigma), diag(sigma))
  testthat::expect_identical(dimnames(stronger_transform$sigma), dimnames(sigma))
  testthat::expect_false(isTRUE(all.equal(stronger_transform$sigma, sigma)))
})

test_that("simulation unused-dot checker names offending arguments", {
  testthat::expect_true(.simulation_check_unused_dots(list()))
  testthat::expect_error(
    .simulation_check_unused_dots(list(extra = 1, 2)),
    "Unused argument\\(s\\): extra, <unnamed>"
  )
})
