# tests/testthat/test-extractRegimeVCVs.R

# testthat::local_edition(3)

# -------------------------------------------------------------------
# Helpers
# -------------------------------------------------------------------
mock_model_output <- function(Pinv, params_named) {
  stopifnot(is.matrix(Pinv), is.numeric(params_named), !is.null(names(params_named)))
  list(
    param = params_named,
    sigma = list(Pinv = Pinv)
  )
}

# -------------------------------------------------------------------
# Tests
# -------------------------------------------------------------------

test_that("Returns NULL when required components are missing", {
  # Missing 'param'
  m1 <- list(sigma = list(Pinv = diag(2)))
  expect_null(extractRegimeVCVs(m1))

  # Missing 'sigma'
  m2 <- list(param = c(r1 = 1))
  expect_null(extractRegimeVCVs(m2))

  # Missing 'Pinv' inside sigma
  m3 <- list(param = c(r1 = 1), sigma = list())
  expect_null(extractRegimeVCVs(m3))
})

test_that("Single-regime: returns list of one matrix equal to base Pinv", {
  Pinv <- matrix(c(2, 0.3,
                   0.3, 1), nrow = 2, byrow = TRUE)
  params <- c(r1 = 1.5)

  model_out <- mock_model_output(Pinv, params)
  vcv_list  <- extractRegimeVCVs(model_out)

  expect_type(vcv_list, "list")
  expect_equal(names(vcv_list), "r1")
  expect_true(is.matrix(vcv_list[[1]]))
  expect_equal(dim(vcv_list[[1]]), dim(Pinv))
  expect_equal(vcv_list[[1]], Pinv)
})

test_that("Multi-regime: later regimes are scaled by param/base_param ratio", {
  Pinv <- matrix(c(1, 0.2,
                   0.2, 0.5), nrow = 2, byrow = TRUE)
  # base_param = 2; others scaled vs 2
  params <- c(r1 = 2, r2 = 3, r3 = 0.5)

  model_out <- mock_model_output(Pinv, params)
  vcv_list  <- extractRegimeVCVs(model_out)

  expect_equal(names(vcv_list), c("r1", "r2", "r3"))
  expect_true(all(vapply(vcv_list, is.matrix, logical(1))))
  expect_true(all(vapply(vcv_list, function(m) all(dim(m) == dim(Pinv)), logical(1))))

  # r1: unscaled
  expect_equal(vcv_list$r1, Pinv)

  # r2: scaled by 3/2
  expect_equal(vcv_list$r2, Pinv * (3/2))

  # r3: scaled by 0.5/2 = 0.25
  expect_equal(vcv_list$r3, Pinv * 0.25)
})

test_that("Works with non-identity, non-diagonal base matrix and preserves numeric type", {
  Pinv <- matrix(c(4, 1, 0.5,
                   1, 3, 0.2,
                   0.5, 0.2, 2),
                 nrow = 3, byrow = TRUE)
  params <- c(A = 1, B = 4)

  model_out <- mock_model_output(Pinv, params)
  vcv_list  <- extractRegimeVCVs(model_out)

  expect_equal(names(vcv_list), c("A", "B"))
  expect_equal(vcv_list$A, Pinv)              # base unchanged
  expect_equal(vcv_list$B, Pinv * (4/1))      # scaled by 4
  expect_true(is.numeric(vcv_list$B[1, 1]))
})

test_that("Handles tiny floating ratios without precision blow-ups", {
  Pinv <- diag(2)
  params <- c(baseline = 1e6, slow = 1e-6)  # ratio = 1e-12

  model_out <- mock_model_output(Pinv, params)
  vcv_list  <- extractRegimeVCVs(model_out)

  expect_equal(vcv_list$baseline, Pinv)
  expect_equal(vcv_list$slow, Pinv * 1e-12, tolerance = 1e-18)
})

# ---- Test XX (NEW): returns NULL when required components are missing ----------
test_that("extractRegimeVCVs returns NULL when required components are missing", {
  testthat::expect_null(extractRegimeVCVs(list()))

  # has param but missing sigma$Pinv
  testthat::expect_null(extractRegimeVCVs(list(
    param = c(r1 = 1),
    sigma = list()
  )))

  # has sigma$Pinv but missing param
  testthat::expect_null(extractRegimeVCVs(list(
    sigma = list(Pinv = diag(2))
  )))
})

# ---- Test XX (NEW): returns scaled list of matrices when components exist -----
test_that("extractRegimeVCVs returns per-regime matrices and scales by param ratio", {
  model_output <- list(
    param = c(r1 = 1, r2 = 3),
    sigma = list(Pinv = matrix(c(1, 2,
                                 3, 4), nrow = 2))
  )

  out <- extractRegimeVCVs(model_output)

  testthat::expect_true(is.list(out))
  testthat::expect_true(all(c("r1", "r2") %in% names(out)))
  testthat::expect_true(is.matrix(out[["r1"]]) && is.matrix(out[["r2"]]))

  # First regime is base matrix
  testthat::expect_equal(out[["r1"]], model_output$sigma$Pinv)

  # Second regime scaled by param ratio (r2/r1 = 3)
  testthat::expect_equal(out[["r2"]], model_output$sigma$Pinv * 3)
})
