testthat::skip_on_cran()

skip_if_simulation_template_deps <- function() {
  testthat::skip_if_not_installed("ape")
  testthat::skip_if_not_installed("phytools")
  testthat::skip_if_not_installed("mvMORPH")
}

test_that("createSimulationTemplate builds an intercept-only template", {
  skip_if_simulation_template_deps()

  set.seed(1)
  tr <- ape::rtree(18)
  X <- matrix(rnorm(18 * 3), ncol = 3)
  rownames(X) <- tr$tip.label

  tmpl <- createSimulationTemplate(
    baseline_tree = tr,
    trait_data = X,
    formula = "trait_data ~ 1",
    method = "LL"
  )

  testthat::expect_s3_class(tmpl, "bifrost_simulation_template")
  testthat::expect_identical(rownames(tmpl$trait_data), tr$tip.label)
  testthat::expect_identical(tmpl$response_columns, 1:3)
  testthat::expect_length(tmpl$predictor_columns, 0)
  testthat::expect_identical(tmpl$fit_method, "LL")
  testthat::expect_null(tmpl$fit_error)
  testthat::expect_equal(tmpl$n_response_traits, 3)
  testthat::expect_true(all(eigen(tmpl$residual_covariance, symmetric = TRUE)$values > 0))
})

test_that("createSimulationTemplate supports numeric formula-based inputs", {
  skip_if_simulation_template_deps()

  set.seed(2)
  tr <- ape::rtree(20)
  dat <- data.frame(
    y1 = rnorm(20),
    y2 = rnorm(20),
    mass = rnorm(20)
  )
  rownames(dat) <- tr$tip.label

  tmpl <- createSimulationTemplate(
    baseline_tree = tr,
    trait_data = dat,
    formula = "trait_data[, 1:2] ~ trait_data[, 3]",
    response_columns = 1:2,
    predictor_columns = 3,
    method = "LL"
  )

  testthat::expect_s3_class(tmpl, "bifrost_simulation_template")
  testthat::expect_false(tmpl$trait_data_is_matrix)
  testthat::expect_identical(tmpl$response_columns, 1:2)
  testthat::expect_identical(tmpl$predictor_columns, 3L)
  testthat::expect_identical(tmpl$formula_mode, "legacy_indexed")
  testthat::expect_identical(tmpl$formula_normalized, "cbind(y1, y2) ~ mass")
  testthat::expect_equal(dim(tmpl$fitted_values), c(20, 2))
})

test_that("createSimulationTemplate supports subset-response intercept formulas", {
  skip_if_simulation_template_deps()

  set.seed(22)
  tr <- ape::rtree(16)
  dat <- cbind(y1 = rnorm(16), y2 = rnorm(16), y3 = rnorm(16))
  rownames(dat) <- tr$tip.label

  tmpl <- createSimulationTemplate(
    baseline_tree = tr,
    trait_data = dat,
    formula = "trait_data[, 1:2] ~ 1",
    response_columns = 1:2,
    method = "LL"
  )

  testthat::expect_s3_class(tmpl, "bifrost_simulation_template")
  testthat::expect_identical(tmpl$response_columns, 1:2)
  testthat::expect_length(tmpl$predictor_columns, 0)
  testthat::expect_equal(ncol(tmpl$trait_data), 2)
  testthat::expect_equal(dim(tmpl$fitted_values), c(16, 2))
})

test_that("createSimulationTemplate validates alignment and column inputs", {
  skip_if_simulation_template_deps()

  set.seed(3)
  tr <- ape::rtree(12)
  X <- matrix(rnorm(12 * 3), ncol = 3)
  rownames(X) <- paste0("sp", seq_len(12))

  testthat::expect_error(
    createSimulationTemplate(
      baseline_tree = tr,
      trait_data = X,
      formula = "trait_data ~ 1",
      method = "LL"
    ),
    "tree tip labels"
  )

  rownames(X) <- tr$tip.label
  tmpl <- createSimulationTemplate(
    baseline_tree = tr,
    trait_data = X,
    formula = "trait_data[, 1:2] ~ trait_data[, 3]",
    method = "LL"
  )
  testthat::expect_identical(tmpl$response_columns, 1:2)
  testthat::expect_identical(tmpl$predictor_columns, 3L)
})

test_that("createSimulationTemplate validates basic inputs and column specifications", {
  skip_if_simulation_template_deps()

  set.seed(4)
  tr <- ape::rtree(10)
  X <- matrix(rnorm(10 * 3), ncol = 3)
  rownames(X) <- tr$tip.label
  colnames(X) <- c("y1", "y2", "mass")

  testthat::expect_error(
    createSimulationTemplate(tr, X, formula = "trait_data ~ 1", model = "OU"),
    "global fit is always BM"
  )
  testthat::expect_error(
    createSimulationTemplate(tr, 1:3, formula = "trait_data ~ 1"),
    "matrix or data.frame"
  )
  testthat::expect_error(
    createSimulationTemplate(tr, unname(X), formula = "trait_data ~ 1"),
    "must have row names"
  )
  testthat::expect_error(
    createSimulationTemplate(tr, X, formula = NA_character_),
    "single character string"
  )

  X_extra <- rbind(X, rogue = c(0, 0, 0))
  testthat::expect_error(
    createSimulationTemplate(tr, X_extra, formula = "trait_data ~ 1"),
    "row names that are not present"
  )

  testthat::expect_error(
    createSimulationTemplate(tr, X, formula = "trait_data ~ 1", predictor_columns = 3),
    "should not be supplied"
  )
  testthat::expect_error(
    createSimulationTemplate(
      tr,
      X,
      formula = "trait_data[, 1:2] ~ trait_data[, 3]",
      response_columns = integer(0),
      predictor_columns = 3
    ),
    "at least one response column"
  )
  testthat::expect_error(
    createSimulationTemplate(
      tr,
      X,
      formula = "trait_data[, 1:2] ~ trait_data[, 3]",
      response_columns = c(1, 99),
      predictor_columns = 3
    ),
    "invalid column positions"
  )
  testthat::expect_error(
    createSimulationTemplate(
      tr,
      X,
      formula = "trait_data[, 1:2] ~ trait_data[, 3]",
      response_columns = c("y1", "missing"),
      predictor_columns = "mass"
    ),
    "names not found"
  )
  testthat::expect_error(
    createSimulationTemplate(
      tr,
      X,
      formula = "trait_data[, 1:2] ~ trait_data[, 3]",
      response_columns = list(1, 2),
      predictor_columns = 3
    ),
    "must be NULL, numeric positions, or character column names"
  )
  testthat::expect_error(
    createSimulationTemplate(
      tr,
      X,
      formula = "trait_data[, 1:2] ~ trait_data[, 3]",
      response_columns = 1:2,
      predictor_columns = 2:3
    ),
    "must not overlap"
  )
  testthat::expect_error(
    createSimulationTemplate(
      tr,
      X,
      formula = "trait_data[, 1:3] ~ trait_data[, 3]",
      response_columns = integer(0)
    ),
    "at least one response column"
  )
  testthat::expect_error(
    createSimulationTemplate(
      tr,
      X,
      formula = "trait_data[, 1:3] ~ trait_data[, 3]"
    ),
    "at least one predictor column"
  )
})

test_that("createSimulationTemplate validates phylo coercion failures", {
  skip_if_simulation_template_deps()

  set.seed(40)
  tr <- ape::rtree(10)
  X <- matrix(rnorm(10 * 3), ncol = 3)
  rownames(X) <- tr$tip.label

  testthat::local_mocked_bindings(
    as.phylo = function(x) structure(list(), class = "notphylo"),
    .package = "ape"
  )

  testthat::expect_error(
    createSimulationTemplate(tr, X, formula = "trait_data ~ 1"),
    "coercible to class 'phylo'"
  )
})

test_that("createSimulationTemplate assigns synthetic names to unnamed matrices", {
  skip_if_simulation_template_deps()

  set.seed(41)
  tr <- ape::rtree(8)
  X <- matrix(rnorm(8 * 2), ncol = 2)
  rownames(X) <- tr$tip.label

  tmpl <- createSimulationTemplate(
    baseline_tree = tr,
    trait_data = X,
    formula = "trait_data[, 1:2] ~ 1"
  )

  testthat::expect_identical(colnames(tmpl$trait_data), c("V1", "V2"))
  testthat::expect_identical(tmpl$response_column_names, c("V1", "V2"))
})

test_that("createSimulationTemplate resolves response and predictor columns by name", {
  skip_if_simulation_template_deps()

  set.seed(5)
  tr <- ape::rtree(12)
  dat <- data.frame(
    y1 = rnorm(12),
    y2 = rnorm(12),
    mass = rnorm(12)
  )
  rownames(dat) <- tr$tip.label

  tmpl <- createSimulationTemplate(
    baseline_tree = tr,
    trait_data = dat,
    formula = "trait_data[, 1:2] ~ trait_data[, 3]",
    response_columns = c("y1", "y2"),
    method = "LL"
  )

  testthat::expect_identical(tmpl$response_columns, 1:2)
  testthat::expect_identical(tmpl$predictor_columns, 3L)
})

test_that("createSimulationTemplate rejects character predictors", {
  skip_if_simulation_template_deps()

  set.seed(6)
  tr <- ape::rtree(10)
  dat <- data.frame(
    y1 = rnorm(10),
    y2 = rnorm(10),
    grp = rep(letters[1:2], each = 5),
    stringsAsFactors = FALSE
  )
  rownames(dat) <- tr$tip.label

  testthat::expect_error(
    createSimulationTemplate(
      baseline_tree = tr,
      trait_data = dat,
      formula = "trait_data[, 1:2] ~ trait_data[, 3]",
      response_columns = 1:2,
      predictor_columns = 3,
      method = "LL"
    ),
    "Character predictors are not supported"
  )
})

test_that("createSimulationTemplate supports named formulas with factor predictors", {
  skip_if_simulation_template_deps()

  set.seed(61)
  tr <- ape::rtree(18)
  dat <- data.frame(
    y1 = rnorm(18),
    y2 = rnorm(18),
    size = exp(rnorm(18)),
    grp = factor(rep(c("a", "b"), length.out = 18))
  )
  rownames(dat) <- tr$tip.label

  tmpl <- createSimulationTemplate(
    baseline_tree = tr,
    trait_data = dat,
    formula = cbind(y1, y2) ~ log(size) + grp,
    method = "LL"
  )

  testthat::expect_identical(tmpl$formula_mode, "named_formula")
  testthat::expect_identical(tmpl$formula_normalized, "cbind(y1, y2) ~ log(size) + grp")
  testthat::expect_identical(tmpl$response_column_names, c("y1", "y2"))
  testthat::expect_identical(tmpl$predictor_column_names, c("size", "grp"))
  testthat::expect_identical(tmpl$predictor_schema$size$type, "numeric")
  testthat::expect_identical(tmpl$predictor_schema$grp$type, "factor")
  testthat::expect_identical(tmpl$predictor_schema$grp$levels, c("a", "b"))
  testthat::expect_s3_class(tmpl$trait_data, "data.frame")
  testthat::expect_identical(colnames(tmpl$data_prototype), colnames(dat))
})

test_that("createSimulationTemplate accepts formula objects", {
  skip_if_simulation_template_deps()

  set.seed(62)
  tr <- ape::rtree(14)
  dat <- data.frame(
    y1 = rnorm(14),
    y2 = rnorm(14),
    size = rnorm(14)
  )
  rownames(dat) <- tr$tip.label

  tmpl <- createSimulationTemplate(
    baseline_tree = tr,
    trait_data = dat,
    formula = cbind(y1, y2) ~ size,
    method = "LL"
  )

  testthat::expect_s3_class(tmpl$formula_original, "formula")
  testthat::expect_identical(tmpl$formula_normalized, "cbind(y1, y2) ~ size")
})

test_that("createSimulationTemplate rejects unsupported formula structures", {
  skip_if_simulation_template_deps()

  set.seed(63)
  tr <- ape::rtree(12)
  dat <- data.frame(
    y1 = rnorm(12),
    y2 = rnorm(12),
    size = exp(rnorm(12))
  )
  rownames(dat) <- tr$tip.label

  testthat::expect_error(
    createSimulationTemplate(
      baseline_tree = tr,
      trait_data = dat,
      formula = log(y1) ~ size,
      method = "LL"
    ),
    "Transformed responses are not supported"
  )

  testthat::expect_error(
    createSimulationTemplate(
      baseline_tree = tr,
      trait_data = dat,
      formula = cbind(y1, y2) ~ .,
      method = "LL"
    ),
    "'.' shorthand is not supported"
  )
})

test_that("createSimulationTemplate handles sigma fallbacks and mocked edge cases", {
  skip_if_simulation_template_deps()

  set.seed(7)
  tr <- ape::rtree(4)
  X1 <- matrix(rnorm(4), ncol = 1)
  rownames(X1) <- tr$tip.label

  testthat::local_mocked_bindings(
    mvgls = function(...) {
      list(
        sigma = list(Pinv = NULL, S = matrix(2, nrow = 1, ncol = 1)),
        fitted = matrix(0, nrow = 4, ncol = 1),
        call = list()
      )
    },
    .package = "mvMORPH"
  )

  tmpl <- createSimulationTemplate(
    baseline_tree = tr,
    trait_data = X1,
    formula = "trait_data ~ 1"
  )

  testthat::expect_equal(tmpl$variance_sd, 0)
  testthat::expect_equal(tmpl$covariance_mean, 0)
  testthat::expect_equal(tmpl$covariance_sd, 0)
  testthat::expect_identical(colnames(tmpl$fitted_values), "V1")
})

test_that("createSimulationTemplate errors when no covariance matrix can be found", {
  skip_if_simulation_template_deps()

  set.seed(8)
  tr <- ape::rtree(4)
  X1 <- matrix(rnorm(4), ncol = 1)
  rownames(X1) <- tr$tip.label

  testthat::local_mocked_bindings(
    mvgls = function(...) {
      list(
        sigma = list(Pinv = NULL, S = NULL),
        fitted = matrix(0, nrow = 4, ncol = 1),
        call = list()
      )
    },
    .package = "mvMORPH"
  )

  testthat::expect_error(
    createSimulationTemplate(
      baseline_tree = tr,
      trait_data = X1,
      formula = "trait_data ~ 1"
    ),
    "Could not locate a covariance matrix"
  )
})

test_that("createSimulationTemplate stores evaluated fit settings from variables", {
  skip_if_simulation_template_deps()

  set.seed(42)
  tr <- ape::rtree(10)
  X <- matrix(rnorm(10 * 2), ncol = 2)
  rownames(X) <- tr$tip.label

  method_value <- "LL"
  error_value <- TRUE

  tmpl <- createSimulationTemplate(
    baseline_tree = tr,
    trait_data = X,
    formula = "trait_data ~ 1",
    method = method_value,
    error = error_value
  )

  testthat::expect_identical(tmpl$fit_method, "LL")
  testthat::expect_identical(tmpl$fit_error, TRUE)
})
