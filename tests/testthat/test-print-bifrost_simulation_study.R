testthat::skip_on_cran()

test_that("print.bifrost_simulation_study prints false-positive summaries", {
  study <- structure(
    list(
      study_type = "false_positive",
      generating_scenario = "null",
      study_summary = list(
        n_replicates = 2L,
        n_completed = 2L,
        n_failed = 0L,
        mean_false_positive_rate = 0.05,
        median_false_positive_rate = 0.04
      )
    ),
    class = c("bifrost_simulation_study", "list")
  )

  out <- paste(capture.output(print(study)), collapse = "\n")

  testthat::expect_match(out, "Bifrost Simulation Study")
  testthat::expect_match(out, "False-Positive Summary")
  testthat::expect_match(out, "Mean FP rate")
})

test_that("print.bifrost_simulation_study prints recovery summaries", {
  study <- structure(
    list(
      study_type = "shift_recovery",
      generating_scenario = "correlation",
      study_summary = list(
        n_replicates = 2L,
        n_completed = 2L,
        n_failed = 0L
      ),
      evaluation = list(
        strict = list(precision = 0.7, recall = 0.6),
        fuzzy = list(precision = 0.8, recall = 0.75)
      )
    ),
    class = c("bifrost_simulation_study", "list")
  )

  out <- paste(capture.output(print(study)), collapse = "\n")

  testthat::expect_match(out, "Recovery Summary")
  testthat::expect_match(out, "Scenario: correlation")
  testthat::expect_match(out, "Strict precision")
})

test_that("print.bifrost_simulation_study collapses multiple scenario labels", {
  study <- structure(
    list(
      study_type = "shift_recovery",
      generating_scenario = c("proportional", "correlation"),
      study_summary = list(
        n_replicates = 2L,
        n_completed = 2L,
        n_failed = 0L
      ),
      evaluation = list(
        strict = list(precision = 0.7, recall = 0.6),
        fuzzy = list(precision = 0.8, recall = 0.75)
      )
    ),
    class = c("bifrost_simulation_study", "list")
  )

  out <- paste(capture.output(print(study)), collapse = "\n")

  testthat::expect_match(out, "Scenario: proportional, correlation")
})
