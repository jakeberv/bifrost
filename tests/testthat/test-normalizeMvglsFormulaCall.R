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
})
