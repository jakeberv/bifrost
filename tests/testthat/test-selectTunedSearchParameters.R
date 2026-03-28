make_tuning_grid_stub <- function(summary_table, IC = "GIC") {
  structure(
    list(
      IC = IC,
      summary_table = summary_table,
      base_search_options = list(
        formula = "trait_data ~ 1",
        method = "LL",
        error = TRUE
      )
    ),
    class = c("bifrost_search_tuning_grid", "list")
  )
}

test_that("selectTunedSearchParameters picks a conservative winner within one IC family", {
  summary_table <- data.frame(
    setting_id = 1:3,
    IC = rep("GIC", 3),
    shift_acceptance_threshold = c(5, 10, 2),
    min_descendant_tips = c(5, 10, 3),
    null_mean_false_positive_rate = c(0.02, 0.02, 0.10),
    null_fraction_any_false_positive = c(0.10, 0.10, 0.30),
    null_evaluable_fraction = c(1, 1, 1),
    proportional_evaluable_fraction = c(1, 1, 1),
    correlation_evaluable_fraction = c(1, 1, 1),
    proportional_fuzzy_f1 = c(0.80, 0.80, 0.95),
    correlation_fuzzy_f1 = c(0.70, 0.70, 0.95),
    proportional_fuzzy_recall = c(0.85, 0.85, 0.98),
    correlation_fuzzy_recall = c(0.75, 0.75, 0.98),
    proportional_strict_f1 = c(0.60, 0.60, 0.80),
    correlation_strict_f1 = c(0.50, 0.50, 0.78),
    stringsAsFactors = FALSE
  )
  tuning_grid <- make_tuning_grid_stub(summary_table)

  selected <- selectTunedSearchParameters(
    tuning_grid,
    max_false_positive_rate = 0.05,
    max_any_false_positive = 0.2,
    min_evaluable_fraction = 0.9,
    primary_metric = "fuzzy_f1",
    tie_break = "conservative"
  )

  testthat::expect_s3_class(selected, "bifrost_search_tuning_selection")
  testthat::expect_identical(selected$selected_row$setting_id, 2L)
  testthat::expect_identical(selected$recommended_search_options$IC, "GIC")
  testthat::expect_identical(selected$recommended_search_options$shift_acceptance_threshold, 10)
  testthat::expect_equal(selected$recommended_search_options$min_descendant_tips, 10)
  testthat::expect_identical(selected$n_feasible_settings, 2L)
  testthat::expect_false(selected$used_all_settings)
})

test_that("selectTunedSearchParameters can rank weighted metrics and fall back to the full grid", {
  summary_table <- data.frame(
    setting_id = 1:2,
    IC = rep("BIC", 2),
    shift_acceptance_threshold = c(2, 10),
    min_descendant_tips = c(3, 8),
    null_mean_false_positive_rate = c(0.20, 0.15),
    null_fraction_any_false_positive = c(0.50, 0.40),
    null_evaluable_fraction = c(1, 1),
    proportional_evaluable_fraction = c(1, 1),
    correlation_evaluable_fraction = c(1, 1),
    proportional_fuzzy_f1 = c(0.70, 0.75),
    correlation_fuzzy_f1 = c(0.60, 0.70),
    proportional_fuzzy_recall = c(0.80, 0.82),
    correlation_fuzzy_recall = c(0.65, 0.74),
    proportional_strict_f1 = c(0.50, 0.60),
    correlation_strict_f1 = c(0.45, 0.55),
    proportional_weighted_fuzzy_f1 = c(0.72, 0.90),
    correlation_weighted_fuzzy_f1 = c(0.62, 0.88),
    stringsAsFactors = FALSE
  )
  tuning_grid <- make_tuning_grid_stub(summary_table, IC = "BIC")

  selected <- selectTunedSearchParameters(
    tuning_grid,
    max_false_positive_rate = 0.01,
    max_any_false_positive = 0.05,
    primary_metric = "weighted_fuzzy_f1",
    scenario_weights = c(0.25, 0.75)
  )

  testthat::expect_true(selected$used_all_settings)
  testthat::expect_identical(selected$n_feasible_settings, 0L)
  testthat::expect_identical(selected$selected_row$setting_id, 2L)
  testthat::expect_identical(selected$recommended_search_options$IC, "BIC")
})

test_that("selectTunedSearchParameters validates its inputs", {
  summary_table <- data.frame(
    setting_id = 1L,
    IC = "GIC",
    shift_acceptance_threshold = 5,
    min_descendant_tips = 5L,
    null_mean_false_positive_rate = 0.02,
    null_fraction_any_false_positive = 0.10,
    null_evaluable_fraction = 1,
    proportional_evaluable_fraction = 1,
    correlation_evaluable_fraction = 1,
    proportional_fuzzy_f1 = 0.8,
    correlation_fuzzy_f1 = 0.7,
    proportional_fuzzy_recall = 0.85,
    correlation_fuzzy_recall = 0.75,
    proportional_strict_f1 = 0.6,
    correlation_strict_f1 = 0.5,
    stringsAsFactors = FALSE
  )
  tuning_grid <- make_tuning_grid_stub(summary_table)

  testthat::expect_error(
    selectTunedSearchParameters(list()),
    "bifrost_search_tuning_grid"
  )
  testthat::expect_error(
    selectTunedSearchParameters(tuning_grid, max_false_positive_rate = -1),
    "max_false_positive_rate"
  )
  testthat::expect_error(
    selectTunedSearchParameters(tuning_grid, max_any_false_positive = 2),
    "between 0 and 1"
  )
  testthat::expect_error(
    selectTunedSearchParameters(tuning_grid, min_evaluable_fraction = 2),
    "between 0 and 1"
  )
  testthat::expect_error(
    selectTunedSearchParameters(tuning_grid, scenario_weights = c(1, -1)),
    "scenario_weights"
  )
  testthat::expect_error(
    selectTunedSearchParameters(tuning_grid, primary_metric = "weighted_fuzzy_f1"),
    "missing required columns"
  )
})
