# tests/testthat/test-search-characterization.R

testthat::skip_on_cran()

if (!exists("bifrost_build_search_characterization", mode = "function")) {
  source(testthat::test_path("helper-search-characterization.R"))
}

test_that("search behavior matches the characterization baseline", {
  bifrost_skip_search_characterization_deps()

  baseline_path <- testthat::test_path(
    "fixtures",
    "search-characterization-baseline.rds"
  )
  if (!file.exists(baseline_path)) {
    testthat::skip(paste("Characterization baseline not found:", baseline_path))
  }

  expected <- readRDS(baseline_path)
  observed <- bifrost_build_search_characterization()

  testthat::expect_equal(
    observed,
    expected,
    tolerance = 1e-4
  )
})
