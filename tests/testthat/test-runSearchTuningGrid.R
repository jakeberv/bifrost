testthat::skip_on_cran()

skip_if_tuning_grid_deps <- function() {
  testthat::skip_if_not_installed("ape")
  testthat::skip_if_not_installed("phytools")
  testthat::skip_if_not_installed("mvMORPH")
  testthat::skip_if_not_installed("RRphylo")
  testthat::skip_if_not_installed("future")
  testthat::skip_if_not_installed("future.apply")
  testthat::skip_if_not_installed("progressr")
}

local_rebind <- function(name, value, env) {
  old_value <- get(name, envir = env, inherits = FALSE)
  was_locked <- bindingIsLocked(name, env)
  if (was_locked) {
    unlockBinding(name, env)
  }
  assign(name, value, envir = env)
  withr::defer({
    assign(name, old_value, envir = env)
    if (was_locked) {
      lockBinding(name, env)
    }
  }, envir = parent.frame())
}

make_tuning_template <- function() {
  set.seed(60)
  tr <- ape::rtree(20)
  X <- matrix(rnorm(20 * 2), ncol = 2)
  rownames(X) <- tr$tip.label
  createSimulationTemplate(
    baseline_tree = tr,
    trait_data = X,
    formula = "trait_data ~ 1",
    method = "LL"
  )
}

make_formula_tuning_template <- function() {
  set.seed(61)
  tr <- ape::rtree(20)
  grp <- factor(rep(c("a", "b"), length.out = 20))
  size <- rnorm(20)
  X <- data.frame(
    y1 = 0.3 * size + ifelse(grp == "b", 0.7, 0) + rnorm(20),
    y2 = -0.25 * size + ifelse(grp == "b", -0.4, 0) + rnorm(20),
    size = size,
    grp = grp,
    row.names = tr$tip.label
  )
  createSimulationTemplate(
    baseline_tree = tr,
    trait_data = X,
    formula = cbind(y1, y2) ~ size + grp,
    method = "LL"
  )
}

test_that("runSearchTuningGrid returns a fixed-IC tuning grid on a small real run", {
  skip_if_tuning_grid_deps()

  tmpl <- make_tuning_template()
  grid <- suppressWarnings(
    runSearchTuningGrid(
      template = tmpl,
      IC = "GIC",
      shift_acceptance_thresholds = 5,
      min_descendant_tips_values = 3,
      tree_tip_count = 16,
      null_replicates = 1,
      recovery_replicates = 1,
      proportional_simulation_options = list(
        num_shifts = 1,
        min_shift_tips = 2,
        max_shift_tips = 5
      ),
      base_search_options = list(
        formula = "trait_data ~ 1",
        method = "LL",
        num_cores = 1
      ),
      weighted = FALSE,
      num_cores = 1,
      seed = 61,
      store_studies = TRUE
    )
  )

  testthat::expect_s3_class(grid, "bifrost_search_tuning_grid")
  testthat::expect_identical(grid$IC, "GIC")
  testthat::expect_equal(nrow(grid$summary_table), 1)
  testthat::expect_length(grid$studies, 1)
  testthat::expect_true(all(c(
    "null_mean_false_positive_rate",
    "proportional_fuzzy_f1",
    "correlation_fuzzy_f1",
    "proportional_evaluable_fraction",
    "correlation_evaluable_fraction"
  ) %in% names(grid$summary_table)))
  testthat::expect_identical(grid$studies[[1]]$null$generating_scenario, "null")
  testthat::expect_identical(grid$studies[[1]]$proportional$generating_scenario, "proportional")
  testthat::expect_identical(grid$studies[[1]]$correlation$generating_scenario, "correlation")
})

test_that("runSearchTuningGrid works in parallel for formula-calibrated templates", {
  skip_if_tuning_grid_deps()

  tmpl <- make_formula_tuning_template()
  grid <- suppressWarnings(
    runSearchTuningGrid(
      template = tmpl,
      IC = "GIC",
      shift_acceptance_thresholds = c(5, 7),
      min_descendant_tips_values = 3,
      tree_tip_count = 16,
      null_replicates = 1,
      recovery_replicates = 1,
      proportional_simulation_options = list(
        num_shifts = 1,
        min_shift_tips = 2,
        max_shift_tips = 5
      ),
      base_search_options = list(
        formula = "trait_data ~ 1",
        method = "LL",
        num_cores = 1
      ),
      weighted = FALSE,
      num_cores = 2,
      seed = 62,
      store_studies = TRUE
    )
  )

  testthat::expect_true(all(vapply(grid$studies, function(x) {
    all(x$proportional$per_replicate$status == "ok")
  }, logical(1))))
  testthat::expect_true(all(vapply(grid$studies, function(x) {
    identical(colnames(x$null$simdata[[1]]$data), c("y1", "y2"))
  }, logical(1))))
})

test_that("runSearchTuningGrid propagates fixed IC and scenario-specific options", {
  skip_if_tuning_grid_deps()

  tmpl <- make_tuning_template()
  mock_env <- new.env(parent = emptyenv())
  mock_env$null_calls <- list()
  mock_env$shift_calls <- list()

  local_rebind(
    "runFalsePositiveSimulationStudy",
    function(template, search_options, seed, num_cores, ...) {
      mock_env$null_calls[[length(mock_env$null_calls) + 1L]] <- list(
        IC = search_options$IC,
        threshold = search_options$shift_acceptance_threshold,
        min_tips = search_options$min_descendant_tips,
        search_num_cores = search_options$num_cores,
        num_cores = num_cores,
        seed = seed
      )
      list(
        generating_scenario = "null",
        per_replicate = data.frame(
          n_inferred_shifts = 1,
          false_positive_rate = 0.05,
          n_candidates = 5,
          status = "ok"
        ),
        study_summary = list(mean_false_positive_rate = 0.05)
      )
    },
    environment(runSearchTuningGrid)
  )
  local_rebind(
    "runShiftRecoverySimulationStudy",
    function(template, simulation_options, search_options, weighted, seed, num_cores, ...) {
      mock_env$shift_calls[[length(mock_env$shift_calls) + 1L]] <- list(
        IC = search_options$IC,
        threshold = search_options$shift_acceptance_threshold,
        min_tips = search_options$min_descendant_tips,
        search_num_cores = search_options$num_cores,
        num_cores = num_cores,
        scale_mode = simulation_options$scale_mode,
        weighted = weighted,
        seed = seed
      )
      scenario <- simulation_options$scale_mode
      list(
        generating_scenario = scenario,
        per_replicate = data.frame(
          n_candidates = 4,
          n_inferred_shifts = 1,
          status = "ok"
        ),
        evaluation = list(
          strict = list(f1 = if (scenario == "proportional") 0.5 else 0.4),
          fuzzy = list(
            f1 = if (scenario == "proportional") 0.8 else 0.6,
            recall = if (scenario == "proportional") 0.9 else 0.7
          ),
          weighted = if (weighted) {
            list(fuzzy = list(f1 = if (scenario == "proportional") 0.75 else 0.55))
          } else {
            NULL
          }
        )
      )
    },
    environment(runSearchTuningGrid)
  )

  grid <- runSearchTuningGrid(
    template = tmpl,
    IC = "BIC",
    shift_acceptance_thresholds = c(2, 10),
    min_descendant_tips_values = 5,
    null_replicates = 2,
    recovery_replicates = 3,
    proportional_simulation_options = list(
      num_shifts = 1,
      min_shift_tips = 2,
      max_shift_tips = 5
    ),
    base_search_options = list(formula = "trait_data ~ 1", method = "LL", num_cores = 7),
    weighted = TRUE,
    num_cores = 1,
    seed = 62,
    store_studies = FALSE
  )

  testthat::expect_null(grid$studies)
  testthat::expect_true(all(grid$summary_table$IC == "BIC"))
  testthat::expect_equal(nrow(grid$summary_table), 2)
  testthat::expect_true(all(vapply(mock_env$null_calls, `[[`, character(1), "IC") == "BIC"))
  testthat::expect_true(all(vapply(mock_env$shift_calls, `[[`, character(1), "IC") == "BIC"))
  testthat::expect_true(all(vapply(mock_env$null_calls, `[[`, numeric(1), "num_cores") == 1))
  testthat::expect_true(all(vapply(mock_env$shift_calls, `[[`, numeric(1), "num_cores") == 1))
  testthat::expect_true(all(vapply(mock_env$null_calls, `[[`, numeric(1), "search_num_cores") == 1))
  testthat::expect_true(all(vapply(mock_env$shift_calls, `[[`, numeric(1), "search_num_cores") == 1))
  testthat::expect_equal(
    vapply(mock_env$shift_calls[c(FALSE, TRUE, FALSE, TRUE)], `[[`, character(1), "scale_mode"),
    c("correlation", "correlation")
  )
  testthat::expect_true(all(c(
    "proportional_weighted_fuzzy_f1",
    "correlation_weighted_fuzzy_f1"
  ) %in% names(grid$summary_table)))
})

test_that("runSearchTuningGrid validates inputs", {
  skip_if_tuning_grid_deps()

  tmpl <- make_tuning_template()

  testthat::expect_error(
    runSearchTuningGrid(
      template = list(),
      IC = "GIC",
      shift_acceptance_thresholds = 5,
      min_descendant_tips_values = 3,
      proportional_simulation_options = list(num_shifts = 1, min_shift_tips = 2, max_shift_tips = 5)
    ),
    "bifrost_simulation_template"
  )
  testthat::expect_error(
    runSearchTuningGrid(
      template = tmpl,
      IC = "GIC",
      shift_acceptance_thresholds = numeric(0),
      min_descendant_tips_values = 3,
      proportional_simulation_options = list(num_shifts = 1, min_shift_tips = 2, max_shift_tips = 5)
    ),
    "shift_acceptance_thresholds"
  )
  testthat::expect_error(
    runSearchTuningGrid(
      template = tmpl,
      IC = "GIC",
      shift_acceptance_thresholds = 5,
      min_descendant_tips_values = 3,
      proportional_simulation_options = list(num_shifts = 1, min_shift_tips = 2, max_shift_tips = 5),
      base_search_options = list(formula = "trait_data[, 1] ~ x")
    ),
    "intercept-only search formulas only"
  )
  testthat::expect_error(
    runSearchTuningGrid(
      template = tmpl,
      IC = "GIC",
      shift_acceptance_thresholds = 5,
      min_descendant_tips_values = 0,
      proportional_simulation_options = list(num_shifts = 1, min_shift_tips = 2, max_shift_tips = 5)
    ),
    "min_descendant_tips_values"
  )
  testthat::expect_error(
    runSearchTuningGrid(
      template = tmpl,
      IC = "GIC",
      shift_acceptance_thresholds = 5,
      min_descendant_tips_values = 3
    ),
    "must be lists"
  )
  testthat::expect_error(
    runSearchTuningGrid(
      template = tmpl,
      IC = "GIC",
      shift_acceptance_thresholds = 5,
      min_descendant_tips_values = 3,
      proportional_simulation_options = list(num_shifts = 1)
    ),
    "must include num_shifts, min_shift_tips, and max_shift_tips"
  )
})

test_that("runSearchTuningGrid can parallelize over settings in a small real run", {
  skip_if_tuning_grid_deps()

  tmpl <- make_tuning_template()
  grid <- suppressWarnings(
    runSearchTuningGrid(
      template = tmpl,
      IC = "GIC",
      shift_acceptance_thresholds = c(2, 5),
      min_descendant_tips_values = 3,
      tree_tip_count = 16,
      null_replicates = 1,
      recovery_replicates = 1,
      proportional_simulation_options = list(
        num_shifts = 1,
        min_shift_tips = 2,
        max_shift_tips = 5
      ),
      base_search_options = list(
        formula = "trait_data ~ 1",
        method = "LL",
        num_cores = 4
      ),
      weighted = FALSE,
      num_cores = 2,
      seed = 63,
      store_studies = FALSE
    )
  )

  testthat::expect_s3_class(grid, "bifrost_search_tuning_grid")
  testthat::expect_equal(nrow(grid$summary_table), 2)
  testthat::expect_true(all(grid$summary_table$IC == "GIC"))
})
