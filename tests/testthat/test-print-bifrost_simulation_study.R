test_that("print.bifrost_simulation_study prints false-positive summaries", {
  study <- structure(
    list(
      study_type = "false_positive",
      generating_scenario = "null",
      study_summary = list(
        n_replicates = 2L,
        n_completed = 2L,
        n_failed = 0L,
        n_evaluable_replicates = 2L,
        mean_false_positive_rate = 0.05,
        median_false_positive_rate = 0.04
      ),
      per_replicate = data.frame(replicate = 1:2, false_positive_rate = c(0.04, 0.06))
    ),
    class = c("bifrost_simulation_study", "list")
  )

  out <- paste(capture.output(print.bifrost_simulation_study(study)), collapse = "\n")

  testthat::expect_match(out, "Bifrost Simulation Study")
  testthat::expect_match(out, "False-Positive Summary")
  testthat::expect_match(out, "Mean FP rate")
  testthat::expect_identical(
    as.data.frame.bifrost_simulation_study(study),
    data.frame(replicate = 1:2, false_positive_rate = c(0.04, 0.06))
  )
  summary <- as.data.frame.bifrost_simulation_study(study, component = "summary")
  testthat::expect_identical(summary$study_type, "false_positive")
  testthat::expect_equal(summary$n_evaluable_replicates, 2L)
})

test_that("print.bifrost_simulation_study prints recovery summaries", {
  study <- structure(
    list(
      study_type = "shift_recovery",
      generating_scenario = "correlation",
      study_summary = list(
        n_replicates = 2L,
        n_completed = 2L,
        n_failed = 0L,
        strict = list(
          precision = 0.7,
          recall = 0.6,
          f1 = 0.646,
          specificity = 0.65,
          fpr = 0.35,
          balanced_accuracy = 0.625
        ),
        fuzzy = list(
          precision = 0.8,
          recall = 0.75,
          f1 = 0.774,
          specificity = 0.85,
          fpr = 0.15,
          balanced_accuracy = 0.8
        )
      ),
      evaluation = list(
        strict = list(
          precision = 0.7,
          recall = 0.6,
          f1 = 0.646,
          specificity = 0.65,
          fpr = 0.35,
          balanced_accuracy = 0.625
        ),
        fuzzy = list(
          precision = 0.8,
          recall = 0.75,
          f1 = 0.774,
          specificity = 0.85,
          fpr = 0.15,
          balanced_accuracy = 0.8
        )
      ),
      per_replicate = data.frame(replicate = 1:2, status = "ok")
    ),
    class = c("bifrost_simulation_study", "list")
  )

  out <- paste(capture.output(print.bifrost_simulation_study(study)), collapse = "\n")

  testthat::expect_match(out, "Recovery Summary")
  testthat::expect_match(out, "Scenario: correlation")
  testthat::expect_match(out, "Strict precision")
  testthat::expect_match(out, "Fuzzy balanced accuracy")
  summary <- as.data.frame.bifrost_simulation_study(study, component = "summary")
  testthat::expect_identical(summary$study_type, "shift_recovery")
  testthat::expect_equal(summary$fuzzy_f1, 0.774)
  testthat::expect_equal(summary$strict_specificity, 0.65)
  testthat::expect_equal(summary$strict_fpr, 0.35)
  testthat::expect_equal(summary$strict_balanced_accuracy, 0.625)
  testthat::expect_equal(summary$fuzzy_specificity, 0.85)
  testthat::expect_equal(summary$fuzzy_fpr, 0.15)
  testthat::expect_equal(summary$fuzzy_balanced_accuracy, 0.8)
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

  out <- paste(capture.output(print.bifrost_simulation_study(study)), collapse = "\n")

  testthat::expect_match(out, "Scenario: proportional, correlation")
})

test_that("as.data.frame.bifrost_simulation_study collapses false-positive scenario labels", {
  study <- structure(
    list(
      study_type = "false_positive",
      generating_scenario = c("null_a", "null_b"),
      study_summary = list(
        n_replicates = 2L,
        n_completed = 2L,
        n_failed = 0L,
        n_evaluable_replicates = 2L,
        mean_false_positive_rate = 0.05,
        median_false_positive_rate = 0.04
      ),
      per_replicate = data.frame(replicate = 1:2)
    ),
    class = c("bifrost_simulation_study", "list")
  )

  summary <- as.data.frame.bifrost_simulation_study(study, component = "summary")

  testthat::expect_identical(summary$generating_scenario, "null_a, null_b")
})
