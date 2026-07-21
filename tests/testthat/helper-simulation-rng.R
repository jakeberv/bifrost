expect_global_rng_restored <- function(code) {
  set.seed(20260709)
  before <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  force(code)
  testthat::expect_true(exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
  if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
    testthat::expect_identical(
      get(".Random.seed", envir = .GlobalEnv, inherits = FALSE),
      before
    )
  }
}
