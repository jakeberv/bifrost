test_that("normalizeMvglsFormulaCall rewrites response-only data.frames for trait_data ~ 1", {
  dat <- data.frame(
    y1 = c(1, 2, 3),
    y2 = c(4, 5, 6),
    row.names = c("sp1", "sp2", "sp3")
  )

  normalized <- bifrost:::normalizeMvglsFormulaCall(
    formula = "trait_data ~ 1",
    trait_data = dat,
    args_list = list()
  )

  testthat::expect_identical(
    paste(deparse(normalized$formula), collapse = " "),
    "cbind(y1, y2) ~ 1"
  )
  testthat::expect_true(is.data.frame(normalized$args_list$data))
  testthat::expect_identical(colnames(normalized$args_list$data), c("y1", "y2"))
  testthat::expect_identical(rownames(normalized$args_list$data), rownames(dat))
})

test_that("normalizeMvglsFormulaCall validates inputs and leaves unsupported legacy forms unchanged", {
  dat <- data.frame(
    y1 = c(1, 2, 3),
    grp = factor(c("a", "b", "a")),
    row.names = c("sp1", "sp2", "sp3")
  )

  testthat::expect_error(
    bifrost:::normalizeMvglsFormulaCall(
      trait_data = dat,
      args_list = list()
    ),
    "must be provided"
  )
  testthat::expect_error(
    bifrost:::normalizeMvglsFormulaCall(
      formula = 1,
      trait_data = dat,
      args_list = list()
    ),
    "must be provided"
  )

  unchanged <- bifrost:::normalizeMvglsFormulaCall(
    formula = "trait_data ~ 1",
    trait_data = dat,
    args_list = list(data = dat)
  )

  testthat::expect_identical(
    paste(deparse(unchanged$formula), collapse = " "),
    "trait_data ~ 1"
  )
  testthat::expect_identical(unchanged$args_list$data, dat)
})

test_that("normalizeMvglsFormulaCall assigns synthetic names and enforces multivariate responses by default", {
  dat <- data.frame(
    a = c(1, 2, 3),
    b = c(4, 5, 6),
    row.names = c("sp1", "sp2", "sp3")
  )
  names(dat) <- c("", "")

  normalized <- bifrost:::normalizeMvglsFormulaCall(
    formula = "trait_data ~ 1",
    trait_data = dat,
    args_list = list(data = dat)
  )

  testthat::expect_identical(
    paste(deparse(normalized$formula), collapse = " "),
    "cbind(Y1, Y2) ~ 1"
  )
  testthat::expect_identical(colnames(normalized$args_list$data), c("Y1", "Y2"))

  testthat::expect_error(
    bifrost:::normalizeMvglsFormulaCall(
      formula = y1 ~ mass,
      trait_data = data.frame(y1 = c(1, 2, 3), mass = c(4, 5, 6)),
      args_list = list(data = data.frame(y1 = c(1, 2, 3), mass = c(4, 5, 6)))
    ),
    "multivariate response"
  )
})

test_that("normalizeMvglsFormulaCall can synthesize multivariate names from matrix responses", {
  resp <- matrix(c(1, 2, 3, 4, 5, 6), ncol = 2)
  dat <- data.frame(
    resp = I(resp),
    mass = c(1, 2, 3),
    row.names = c("sp1", "sp2", "sp3")
  )

  normalized <- bifrost:::normalizeMvglsFormulaCall(
    formula = resp ~ mass,
    trait_data = dat,
    args_list = list(data = dat)
  )

  testthat::expect_identical(
    paste(deparse(normalized$formula), collapse = " "),
    "cbind(Y1, Y2) ~ mass"
  )
  testthat::expect_identical(colnames(normalized$args_list$data), c("Y1", "Y2", "mass"))
})

test_that("normalizeMvglsFormulaCall rewrites legacy indexed formulas onto named columns", {
  dat <- data.frame(
    y1 = c(1, 2, 3),
    y2 = c(4, 5, 6),
    y3 = c(7, 8, 9),
    mass = c(10, 11, 12),
    row.names = c("sp1", "sp2", "sp3")
  )

  normalized <- bifrost:::normalizeMvglsFormulaCall(
    formula = "trait_data[, 2:3] ~ trait_data[, 1]",
    trait_data = dat,
    args_list = list(data = dat)
  )

  testthat::expect_identical(
    paste(deparse(normalized$formula), collapse = " "),
    "cbind(y2, y3) ~ y1"
  )
  testthat::expect_identical(rownames(normalized$args_list$data), rownames(dat))

  intercept_only <- bifrost:::normalizeMvglsFormulaCall(
    formula = "trait_data[, 1:2] ~ 1",
    trait_data = dat,
    args_list = list(data = dat)
  )

  testthat::expect_identical(
    paste(deparse(intercept_only$formula), collapse = " "),
    "cbind(y1, y2) ~ 1"
  )
  testthat::expect_identical(colnames(intercept_only$args_list$data), c("y1", "y2"))
})

test_that("normalizeMvglsFormulaCall handles legacy edge cases without changing formulas unnecessarily", {
  dat <- data.frame(
    y1 = c(1, 2, 3),
    y2 = c(4, 5, 6),
    trait_data = c(7, 8, 9),
    row.names = c("sp1", "sp2", "sp3")
  )

  unchanged <- bifrost:::normalizeMvglsFormulaCall(
    formula = "cbind(y1, y2) ~ trait_data",
    trait_data = dat,
    args_list = list(data = dat)
  )

  testthat::expect_identical(
    paste(deparse(unchanged$formula), collapse = " "),
    "cbind(y1, y2) ~ trait_data"
  )
  testthat::expect_identical(unchanged$args_list$data, dat)

  unnamed_dat <- structure(
    list(c(1, 2, 3), c(4, 5, 6)),
    class = "data.frame",
    row.names = c("sp1", "sp2", "sp3")
  )

  normalized_unnamed <- bifrost:::normalizeMvglsFormulaCall(
    formula = "trait_data[, 2] ~ trait_data[, 1]",
    trait_data = unnamed_dat,
    args_list = list(data = unnamed_dat),
    allow_single_response = TRUE
  )

  testthat::expect_identical(
    paste(deparse(normalized_unnamed$formula), collapse = " "),
    "cbind(V2) ~ V1"
  )
  testthat::expect_identical(colnames(normalized_unnamed$args_list$data), c("V2", "V1"))

  drop_names_method <- function(x, i, j, ..., drop = FALSE) {
    out <- NextMethod("[")
    if (is.data.frame(out)) {
      names(out) <- rep("", ncol(out))
    }
    out
  }
  assign("[.bifrost_drop_names_df", drop_names_method, envir = globalenv())
  on.exit(rm("[.bifrost_drop_names_df", envir = globalenv()), add = TRUE)

  drop_names_dat <- data.frame(
    y1 = c(1, 2, 3),
    y2 = c(4, 5, 6),
    row.names = c("sp1", "sp2", "sp3")
  )
  class(drop_names_dat) <- c("bifrost_drop_names_df", class(drop_names_dat))

  normalized_drop_names <- bifrost:::normalizeMvglsFormulaCall(
    formula = "trait_data[, 1:2] ~ 1",
    trait_data = drop_names_dat,
    args_list = list(data = drop_names_dat)
  )

  testthat::expect_identical(
    paste(deparse(normalized_drop_names$formula), collapse = " "),
    "cbind(Y1, Y2) ~ 1"
  )
  testthat::expect_identical(colnames(normalized_drop_names$args_list$data), c("Y1", "Y2"))
})

test_that("normalizeMvglsFormulaCall can rewrite single-response data.frames when allowed", {
  dat <- data.frame(
    y1 = c(1, 2, 3),
    row.names = c("sp1", "sp2", "sp3")
  )

  normalized <- bifrost:::normalizeMvglsFormulaCall(
    formula = "trait_data ~ 1",
    trait_data = dat,
    args_list = list(),
    allow_single_response = TRUE
  )

  testthat::expect_identical(
    paste(deparse(normalized$formula), collapse = " "),
    "cbind(y1) ~ 1"
  )
  testthat::expect_true(is.data.frame(normalized$args_list$data))
  testthat::expect_identical(colnames(normalized$args_list$data), "y1")
  testthat::expect_identical(rownames(normalized$args_list$data), rownames(dat))

  indexed <- bifrost:::normalizeMvglsFormulaCall(
    formula = "trait_data[, 1] ~ 1",
    trait_data = dat,
    args_list = list(data = dat),
    allow_single_response = TRUE
  )

  testthat::expect_identical(
    paste(deparse(indexed$formula), collapse = " "),
    "cbind(y1) ~ 1"
  )
  testthat::expect_identical(colnames(indexed$args_list$data), "y1")
})

test_that("normalizeMvglsFormulaCall can synthesize names for unnamed single responses", {
  dat <- data.frame(
    y1 = c(1, 2, 3),
    mass = c(4, 5, 6),
    row.names = c("sp1", "sp2", "sp3")
  )

  normalized <- bifrost:::normalizeMvglsFormulaCall(
    formula = I(y1) ~ mass,
    trait_data = dat,
    args_list = list(data = dat),
    allow_single_response = TRUE
  )

  testthat::expect_identical(
    paste(deparse(normalized$formula), collapse = " "),
    "cbind(Y1) ~ mass"
  )
  testthat::expect_true(is.data.frame(normalized$args_list$data))
  testthat::expect_identical(colnames(normalized$args_list$data), c("Y1", "mass"))
})
