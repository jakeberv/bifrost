test_that("tuning recommendations preserve inherited fitting settings and explicit overrides", {
  withr::local_seed(60)
  tree <- ape::rtree(20)
  traits <- matrix(rnorm(40), 20, 2, dimnames = list(tree$tip.label, NULL))
  template <- createSimulationTemplate(tree, traits, method = "LL", error = TRUE)
  cases <- list(
    list(options = list(), method = "LL", error = TRUE),
    list(options = list(method = "PL-LOOCV", error = FALSE),
         method = "PL-LOOCV", error = FALSE),
    list(options = list(method = NULL, error = NULL),
         method = NULL, error = NULL)
  )
  for (case in cases) {
    grid <- runSearchTuningGrid(
      template, IC = "GIC", shift_acceptance_thresholds = 1e6,
      min_descendant_tips_values = 2, null_replicates = 1,
      recovery_replicates = 1,
      proportional_simulation_options = list(
        num_shifts = 1, min_shift_tips = 2, max_shift_tips = 5, buffer = 0
      ),
      base_search_options = case$options, weighted = FALSE,
      seed = 17, store_studies = TRUE
    )
    expect_equal(unname(unlist(grid$summary_table[c(
      "null_completion_rate", "proportional_completion_rate",
      "correlation_completion_rate"
    )])), rep(1, 3))
    selected <- selectTunedSearchParameters(grid)
    expect_identical(selected$recommended_search_options$method, case$method)
    expect_identical(selected$recommended_search_options$error, case$error)
    for (study in grid$studies[[1L]]) {
      expect_identical(study$search_options$method, case$method)
      expect_identical(study$search_options$error, case$error)
    }
    # Recommendations must retain the settings even when raw studies are not saved.
    grid$studies <- NULL
    compact_selected <- selectTunedSearchParameters(grid)
    expect_identical(compact_selected$recommended_search_options,
                     selected$recommended_search_options)
  }
})
