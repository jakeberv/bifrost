test_that("simulation formula helpers validate and normalize direct inputs", {
  dat_names <- c("y1", "y2", "mass", "grp")
  dat <- data.frame(
    y1 = rnorm(4),
    y2 = rnorm(4),
    mass = rnorm(4),
    grp = factor(c("a", "b", "a", "b"))
  )

  testthat::expect_error(
    bifrost:::simulationFormulaToString(NA_character_),
    "single character string or a formula object"
  )
  testthat::expect_identical(
    bifrost:::validateSimulationStudyFormula(stats::as.formula("trait_data ~ 1")),
    stats::as.formula("trait_data ~ 1")
  )
  testthat::expect_error(
    bifrost:::validateSimulationStudyFormula("cbind(y1, y2) ~ mass"),
    "intercept-only search formulas only"
  )

  testthat::expect_true(bifrost:::isSupportedNamedResponse(quote(y1)))
  testthat::expect_true(bifrost:::isSupportedNamedResponse(quote(cbind(y1, y2))))
  testthat::expect_false(bifrost:::isSupportedNamedResponse(quote(trait_data)))
  testthat::expect_false(bifrost:::isSupportedNamedResponse(quote(log(y1))))

  testthat::expect_error(
    bifrost:::extractNamedResponseNames(quote(y3), dat_names),
    "was not found in trait_data"
  )
  testthat::expect_error(
    bifrost:::extractNamedResponseNames(quote(cbind(y1, y3)), dat_names),
    "were not found in trait_data"
  )
  testthat::expect_error(
    bifrost:::extractNamedResponseNames(quote(cbind(y1, y1)), dat_names),
    "must be unique"
  )
  testthat::expect_error(
    bifrost:::extractNamedResponseNames(quote(log(y1)), dat_names),
    "Unsupported response form"
  )

  testthat::expect_error(
    bifrost:::resolveLegacyTraitDataColumns(quote(trait_data[, missing_col]), dat_names),
    "could not be resolved from the supplied formula"
  )
  testthat::expect_error(
    bifrost:::resolveLegacyTraitDataColumns(quote(data.frame(other = 1)), dat_names),
    "could not be resolved to named columns"
  )
  testthat::expect_error(
    bifrost:::resolveLegacyTraitDataColumns(quote(c(99, 100)), dat_names),
    "could not be resolved to valid columns"
  )

  testthat::expect_identical(
    bifrost:::collectLegacyTraitDataColumns(quote(mass), dat_names),
    integer(0)
  )
  testthat::expect_identical(
    bifrost:::collectLegacyTraitDataColumns(
      quote(log(trait_data[, 3]) + trait_data[, 4]),
      dat_names
    ),
    c(3L, 4L)
  )

  testthat::expect_identical(
    as.character(bifrost:::rewriteLegacyTraitDataExpr(quote(trait_data), dat_names, side = "lhs"))[1],
    "cbind"
  )
  testthat::expect_identical(
    bifrost:::rewriteLegacyTraitDataExpr(quote(trait_data), "y1", side = "lhs"),
    quote(y1)
  )
  testthat::expect_identical(
    bifrost:::rewriteLegacyTraitDataExpr(quote(trait_data), dat_names, side = "rhs"),
    quote(trait_data)
  )
  testthat::expect_identical(
    as.character(bifrost:::rewriteLegacyTraitDataExpr(quote(trait_data[, 1:2]), dat_names, side = "lhs"))[1],
    "cbind"
  )
  testthat::expect_identical(
    bifrost:::rewriteLegacyTraitDataExpr(quote(log(trait_data[, 3])), dat_names, side = "rhs"),
    quote(log(mass))
  )
  testthat::expect_error(
    bifrost:::rewriteLegacyTraitDataExpr(quote(trait_data[, 1:2]), dat_names, side = "rhs"),
    "must resolve to single raw columns"
  )

  dat$flag <- c(TRUE, FALSE, TRUE, FALSE)
  schema <- bifrost:::predictorSchemaForData(dat, c("mass", "grp", "flag"))
  testthat::expect_identical(schema$mass$type, "numeric")
  testthat::expect_identical(schema$grp$type, "factor")
  testthat::expect_identical(schema$flag$type, "logical")
  testthat::expect_identical(bifrost:::predictorSchemaForData(dat, character(0)), list())
  testthat::expect_error(
    bifrost:::predictorSchemaForData(
      data.frame(y1 = I(list(1:2, 1:2)), row.names = c("a", "b")),
      "y1"
    ),
    "must be numeric, logical, factor, or ordered factor"
  )
})

test_that("normalizeSimulationFormulaSpec covers legacy, named, and validation branches", {
  dat <- data.frame(
    y1 = rnorm(4),
    y2 = rnorm(4),
    mass = rnorm(4),
    grp = factor(c("a", "b", "a", "b")),
    row.names = paste0("sp", 1:4)
  )

  spec_single <- bifrost:::normalizeSimulationFormulaSpec(
    formula = y1 ~ 1,
    trait_data = dat
  )
  testthat::expect_true(spec_single$intercept_only)
  testthat::expect_identical(spec_single$response_column_names, "y1")

  spec_trait_data <- bifrost:::normalizeSimulationFormulaSpec(
    formula = trait_data ~ mass,
    trait_data = dat,
    response_columns = c("y1", "y2")
  )
  testthat::expect_identical(spec_trait_data$formula_mode, "legacy_indexed")
  testthat::expect_identical(spec_trait_data$predictor_column_names, c("mass", "grp"))

  spec_named_single <- bifrost:::normalizeSimulationFormulaSpec(
    formula = y1 ~ mass,
    trait_data = dat
  )
  testthat::expect_identical(spec_named_single$formula_normalized_chr, "y1 ~ mass")

  testthat::expect_error(
    bifrost:::normalizeSimulationFormulaSpec(
      formula = y1 ~ 1,
      trait_data = dat,
      predictor_columns = integer(0)
    ),
    "at least one predictor column when supplied"
  )
  testthat::expect_error(
    bifrost:::normalizeSimulationFormulaSpec(
      formula = stats::as.formula("trait_data[, 0] ~ 1"),
      trait_data = dat
    ),
    "response_columns must identify at least one response column"
  )
  testthat::expect_error(
    bifrost:::normalizeSimulationFormulaSpec(
      formula = log(y1) ~ mass,
      trait_data = dat
    ),
    "Transformed responses are not supported"
  )
  testthat::expect_error(
    bifrost:::normalizeSimulationFormulaSpec(
      formula = I(y1) ~ 1,
      trait_data = dat
    ),
    "Unsupported response form in intercept-only simulation formula"
  )
  testthat::expect_error(
    bifrost:::normalizeSimulationFormulaSpec(
      formula = trait_data ~ mass,
      trait_data = dat
    ),
    "response_columns must be supplied"
  )
  testthat::expect_error(
    bifrost:::normalizeSimulationFormulaSpec(
      formula = trait_data[, 1:2] ~ mass,
      trait_data = dat,
      response_columns = 1L
    ),
    "response_columns conflict"
  )
  testthat::expect_error(
    bifrost:::normalizeSimulationFormulaSpec(
      formula = cbind(y1, y2) ~ missing_mass,
      trait_data = dat
    ),
    "may only reference raw columns present in trait_data"
  )
  testthat::expect_error(
    bifrost:::normalizeSimulationFormulaSpec(
      formula = cbind(y1, y2) ~ mass,
      trait_data = dat,
      response_columns = "y1"
    ),
    "response_columns conflict"
  )
  testthat::expect_error(
    bifrost:::normalizeSimulationFormulaSpec(
      formula = cbind(y1, y2) ~ mass,
      trait_data = dat,
      predictor_columns = "grp"
    ),
    "predictor_columns conflict"
  )
  testthat::expect_error(
    bifrost:::normalizeSimulationFormulaSpec(
      formula = stats::as.formula("trait_data[, 0] ~ mass"),
      trait_data = dat,
      predictor_columns = "mass"
    ),
    "must identify at least one response column"
  )
  testthat::expect_error(
    bifrost:::normalizeSimulationFormulaSpec(
      formula = 1 ~ mass,
      trait_data = dat
    ),
    "Unsupported response form in simulation formula"
  )
})
