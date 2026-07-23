skip_if_tuning_grid_deps <- function() {
  testthat::skip_if_not_installed("ape")
  testthat::skip_if_not_installed("phytools")
  testthat::skip_if_not_installed("mvMORPH")
  testthat::skip_if_not_installed("RRphylo")
  testthat::skip_if_not_installed("future")
  testthat::skip_if_not_installed("future.apply")
  testthat::skip_if_not_installed("progressr")
  withr::local_options(list(progressr.enable = FALSE))
}

.tuning_template_cache <- new.env(parent = emptyenv())

get_cached_tuning_template <- function(key, builder) {
  if (!exists(key, envir = .tuning_template_cache, inherits = FALSE)) {
    assign(key, builder(), envir = .tuning_template_cache)
  }
  unserialize(serialize(get(key, envir = .tuning_template_cache, inherits = FALSE), NULL))
}

make_tuning_template <- function() {
  get_cached_tuning_template("tuning_template", function() {
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
  })
}

make_formula_tuning_template <- function() {
  get_cached_tuning_template("formula_tuning_template", function() {
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
  })
}

test_that("bifrost_search_tuning_grid has print and data-frame methods", {
  summary_table <- data.frame(
    setting_id = 1L,
    IC = "GIC",
    shift_acceptance_threshold = 5,
    min_descendant_tips = 3L,
    null_mean_false_positive_rate = 0.01,
    proportional_fuzzy_f1 = 0.8,
    correlation_fuzzy_f1 = 0.7,
    stringsAsFactors = FALSE
  )
  grid <- structure(
    list(
      IC = "GIC",
      grid = summary_table[, c("setting_id", "shift_acceptance_threshold", "min_descendant_tips")],
      summary_table = summary_table,
      studies = NULL,
      null_replicates = 2L,
      recovery_replicates = 3L,
      weighted = FALSE,
      store_studies = FALSE
    ),
    class = c("bifrost_search_tuning_grid", "list")
  )

  testthat::expect_output(print(grid), "Bifrost Search Tuning Grid")
  testthat::expect_identical(as.data.frame(grid), summary_table)
})

test_that("runSearchTuningGrid returns a fixed-IC tuning grid on a small real run", {
  testthat::skip_on_cran()
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
      null_simulation_options = list(),
      proportional_simulation_options = list(
        num_shifts = 1,
        min_shift_tips = 2,
        max_shift_tips = 5
      ),
      correlation_simulation_options = list(),
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
  testthat::expect_identical(
    grid$simulation_generators,
    c(
      null = "original",
      proportional = "original",
      correlation = "original"
    )
  )
  testthat::expect_false("simulation_designs" %in% names(grid))
  testthat::expect_identical(
    grid$studies[[1]]$correlation$simulation_generator,
    "original"
  )
})

test_that("runSearchTuningGrid works in parallel for formula-calibrated templates", {
  testthat::skip_on_cran()
  skip_if_tuning_grid_deps()
  testthat::skip_on_covr()

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
  testthat::skip_on_cran()
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
        search_progress = search_options$progress,
        num_cores = num_cores,
        seed = seed
      )
      list(
        generating_scenario = "null",
        per_replicate = data.frame(
          n_inferred_shifts = c(1, NA_real_),
          false_positive_rate = c(0.05, NA_real_),
          n_candidates = c(5, 5),
          status = c("ok", "error")
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
        search_progress = search_options$progress,
        num_cores = num_cores,
        scale_mode = simulation_options$scale_mode,
        weighted = weighted,
        seed = seed
      )
      scenario <- simulation_options$scale_mode
      list(
        generating_scenario = scenario,
        per_replicate = data.frame(
          n_candidates = c(4, 4),
          n_inferred_shifts = c(1, NA_real_),
          status = c("ok", "error")
        ),
        evaluation = list(
          strict = if (scenario == "proportional") {
            list(
              precision = 0.11,
              recall = 0.12,
              f1 = 0.13,
              specificity = 0.14,
              fpr = 0.15,
              balanced_accuracy = 0.16
            )
          } else {
            list(
              precision = 0.21,
              recall = 0.22,
              f1 = 0.23,
              specificity = 0.24,
              fpr = 0.25,
              balanced_accuracy = 0.26
            )
          },
          fuzzy = if (scenario == "proportional") {
            list(
              precision = 0.31,
              recall = 0.32,
              f1 = 0.33,
              specificity = 0.34,
              fpr = 0.35,
              balanced_accuracy = 0.36
            )
          } else {
            list(
              precision = 0.41,
              recall = 0.42,
              f1 = 0.43,
              specificity = 0.44,
              fpr = 0.45,
              balanced_accuracy = 0.46
            )
          },
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
    base_search_options = list(
      formula = "trait_data ~ 1",
      method = "LL",
      num_cores = 7,
      progress = TRUE
    ),
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
  testthat::expect_false(any(vapply(mock_env$null_calls, `[[`, logical(1), "search_progress")))
  testthat::expect_false(any(vapply(mock_env$shift_calls, `[[`, logical(1), "search_progress")))
  testthat::expect_true(all(c("null_seed", "proportional_seed", "correlation_seed") %in% names(grid$summary_table)))
  testthat::expect_equal(
    grid$summary_table$null_seed,
    vapply(mock_env$null_calls, `[[`, numeric(1), "seed")
  )
  testthat::expect_equal(
    grid$summary_table$proportional_seed,
    vapply(mock_env$shift_calls[c(TRUE, FALSE, TRUE, FALSE)], `[[`, numeric(1), "seed")
  )
  testthat::expect_equal(
    grid$summary_table$correlation_seed,
    vapply(mock_env$shift_calls[c(FALSE, TRUE, FALSE, TRUE)], `[[`, numeric(1), "seed")
  )
  testthat::expect_equal(anyDuplicated(unlist(grid$summary_table[
    ,
    c("null_seed", "proportional_seed", "correlation_seed")
  ])), 0L)
  testthat::expect_equal(
    vapply(mock_env$shift_calls[c(FALSE, TRUE, FALSE, TRUE)], `[[`, character(1), "scale_mode"),
    c("correlation", "correlation")
  )
  testthat::expect_true(all(c(
    "proportional_weighted_fuzzy_f1",
    "correlation_weighted_fuzzy_f1"
  ) %in% names(grid$summary_table)))
  recovery_metric_columns <- unlist(lapply(
    c("proportional", "correlation"),
    function(scenario) {
      paste(
        scenario,
        rep(c("strict", "fuzzy"), each = 6L),
        rep(c("precision", "recall", "f1", "specificity", "fpr", "balanced_accuracy"), 2L),
        sep = "_"
      )
    }
  ))
  testthat::expect_true(all(recovery_metric_columns %in% names(grid$summary_table)))
  testthat::expect_equal(
    unname(unlist(grid$summary_table[1L, paste("proportional", c(
      "strict_precision", "strict_recall", "strict_f1", "strict_specificity", "strict_fpr", "strict_balanced_accuracy",
      "fuzzy_precision", "fuzzy_recall", "fuzzy_f1", "fuzzy_specificity", "fuzzy_fpr", "fuzzy_balanced_accuracy"
    ), sep = "_")])),
    c(0.11, 0.12, 0.13, 0.14, 0.15, 0.16, 0.31, 0.32, 0.33, 0.34, 0.35, 0.36)
  )
  testthat::expect_equal(
    unname(unlist(grid$summary_table[1L, paste("correlation", c(
      "strict_precision", "strict_recall", "strict_f1", "strict_specificity", "strict_fpr", "strict_balanced_accuracy",
      "fuzzy_precision", "fuzzy_recall", "fuzzy_f1", "fuzzy_specificity", "fuzzy_fpr", "fuzzy_balanced_accuracy"
    ), sep = "_")])),
    c(0.21, 0.22, 0.23, 0.24, 0.25, 0.26, 0.41, 0.42, 0.43, 0.44, 0.45, 0.46)
  )
  testthat::expect_equal(grid$summary_table$null_evaluable_fraction, rep(1, 2))
  testthat::expect_equal(grid$summary_table$proportional_evaluable_fraction, rep(1, 2))
  testthat::expect_equal(grid$summary_table$correlation_evaluable_fraction, rep(1, 2))
  testthat::expect_equal(grid$summary_table$null_completion_rate, rep(0.5, 2))
  testthat::expect_equal(grid$summary_table$proportional_completion_rate, rep(0.5, 2))
  testthat::expect_equal(grid$summary_table$correlation_completion_rate, rep(0.5, 2))
  testthat::expect_equal(grid$summary_table$null_failure_rate, rep(0.5, 2))
  testthat::expect_equal(grid$summary_table$proportional_failure_rate, rep(0.5, 2))
  testthat::expect_equal(grid$summary_table$correlation_failure_rate, rep(0.5, 2))
})

test_that("runSearchTuningGrid reports missing means when every replicate fails", {
  testthat::skip_on_cran()
  skip_if_tuning_grid_deps()

  tmpl <- make_tuning_template()
  local_rebind(
    "runFalsePositiveSimulationStudy",
    function(...) {
      list(
        generating_scenario = "null",
        per_replicate = data.frame(
          n_inferred_shifts = NA_real_,
          false_positive_rate = NA_real_,
          n_candidates = NA_real_,
          status = "error"
        ),
        study_summary = list(mean_false_positive_rate = NA_real_)
      )
    },
    environment(runSearchTuningGrid)
  )
  local_rebind(
    "runShiftRecoverySimulationStudy",
    function(template, simulation_options, ...) {
      list(
        generating_scenario = simulation_options$scale_mode,
        per_replicate = data.frame(
          n_candidates = NA_real_,
          n_inferred_shifts = NA_real_,
          status = "error"
        ),
        evaluation = list(
          strict = list(f1 = NA_real_),
          fuzzy = list(f1 = NA_real_, recall = NA_real_),
          weighted = NULL
        )
      )
    },
    environment(runSearchTuningGrid)
  )

  grid <- runSearchTuningGrid(
    template = tmpl,
    IC = "GIC",
    shift_acceptance_thresholds = 5,
    min_descendant_tips_values = 3,
    null_replicates = 1,
    recovery_replicates = 1,
    proportional_simulation_options = list(
      num_shifts = 1,
      min_shift_tips = 2,
      max_shift_tips = 5
    ),
    base_search_options = list(formula = "trait_data ~ 1", method = "LL"),
    weighted = FALSE,
    num_cores = 1,
    seed = 63,
    store_studies = FALSE
  )

  testthat::expect_true(is.na(grid$summary_table$null_fraction_any_false_positive))
  testthat::expect_true(is.na(grid$summary_table$null_mean_inferred_shifts))
  testthat::expect_true(is.na(grid$summary_table$proportional_mean_inferred_shifts))
  testthat::expect_true(is.na(grid$summary_table$correlation_mean_inferred_shifts))
  testthat::expect_equal(grid$summary_table$null_failure_rate, 1)
  testthat::expect_equal(grid$summary_table$proportional_failure_rate, 1)
  testthat::expect_equal(grid$summary_table$correlation_failure_rate, 1)
})

test_that("runSearchTuningGrid inherits template defaults and can choose multisession plans", {
  testthat::skip_on_cran()
  skip_if_tuning_grid_deps()

  tmpl <- make_tuning_template()
  tmpl$search_formula <- NULL
  mock_env <- new.env(parent = emptyenv())
  mock_env$null_calls <- list()
  mock_env$shift_calls <- list()
  plans <- list()

  old_rstudio <- Sys.getenv("RSTUDIO", unset = NA_character_)
  old_rstudio_init <- Sys.getenv("RSTUDIO_SESSION_INITIALIZED", unset = NA_character_)
  Sys.setenv(RSTUDIO = "1", RSTUDIO_SESSION_INITIALIZED = "1")
  on.exit({
    if (is.na(old_rstudio)) Sys.unsetenv("RSTUDIO") else Sys.setenv(RSTUDIO = old_rstudio)
    if (is.na(old_rstudio_init)) Sys.unsetenv("RSTUDIO_SESSION_INITIALIZED") else Sys.setenv(RSTUDIO_SESSION_INITIALIZED = old_rstudio_init)
  }, add = TRUE)

  testthat::local_mocked_bindings(
    plan = function(strategy, workers = NULL, ...) {
      if (missing(strategy)) {
        return("old_plan")
      }
      plans[[length(plans) + 1L]] <<- list(strategy = strategy, workers = workers)
      invisible("mock_plan")
    },
    .package = "future"
  )
  testthat::local_mocked_bindings(
    future_lapply = function(X, FUN, ...) lapply(X, FUN),
    .package = "future.apply"
  )

  local_rebind(
    "runFalsePositiveSimulationStudy",
    function(template, search_options, seed, ...) {
      mock_env$null_calls[[length(mock_env$null_calls) + 1L]] <- list(
        formula = search_options$formula,
        seed = seed
      )
      list(
        generating_scenario = "null",
        per_replicate = data.frame(
          n_inferred_shifts = 0,
          false_positive_rate = 0,
          n_candidates = 5,
          status = "ok"
        ),
        study_summary = list(mean_false_positive_rate = 0)
      )
    },
    environment(runSearchTuningGrid)
  )
  local_rebind(
    "runShiftRecoverySimulationStudy",
    function(template, simulation_options, search_options, seed, ...) {
      mock_env$shift_calls[[length(mock_env$shift_calls) + 1L]] <- list(
        formula = search_options$formula,
        scale_mode = simulation_options$scale_mode,
        max_shift_tips = simulation_options$max_shift_tips,
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
          strict = list(f1 = 0.4),
          fuzzy = list(f1 = 0.6, recall = 0.7),
          weighted = NULL
        )
      )
    },
    environment(runSearchTuningGrid)
  )

  grid <- runSearchTuningGrid(
    template = tmpl,
    IC = "GIC",
    shift_acceptance_thresholds = 5,
    min_descendant_tips_values = 3,
    null_replicates = 1,
    recovery_replicates = 1,
    proportional_simulation_options = list(
      num_shifts = 1,
      min_shift_tips = 2,
      max_shift_tips = 5
    ),
    correlation_simulation_options = list(max_shift_tips = 9),
    base_search_options = list(method = "LL", num_cores = 1),
    weighted = FALSE,
    num_cores = 2,
    seed = NULL,
    store_studies = FALSE
  )

  testthat::expect_identical(mock_env$null_calls[[1]]$formula, "trait_data ~ 1")
  testthat::expect_identical(mock_env$shift_calls[[1]]$formula, "trait_data ~ 1")
  testthat::expect_identical(mock_env$shift_calls[[2]]$scale_mode, "correlation")
  testthat::expect_identical(mock_env$shift_calls[[2]]$max_shift_tips, 9)
  testthat::expect_null(mock_env$null_calls[[1]]$seed)
  testthat::expect_true(any(vapply(plans, function(x) identical(x$strategy, future::multisession), logical(1))))
  testthat::expect_equal(nrow(grid$summary_table), 1)

  plans <- list()
  Sys.unsetenv(c("RSTUDIO", "RSTUDIO_SESSION_INITIALIZED"))
  grid <- runSearchTuningGrid(
    template = tmpl,
    IC = "GIC",
    shift_acceptance_thresholds = 5,
    min_descendant_tips_values = 3,
    null_replicates = 1,
    recovery_replicates = 1,
    proportional_simulation_options = list(
      num_shifts = 1,
      min_shift_tips = 2,
      max_shift_tips = 5
    ),
    correlation_simulation_options = list(max_shift_tips = 9),
    base_search_options = list(method = "LL", num_cores = 1),
    weighted = FALSE,
    num_cores = 2,
    seed = NULL,
    store_studies = FALSE
  )

  expected_strategy <- if (.Platform$OS.type == "windows") {
    future::multisession
  } else {
    future::multicore
  }
  testthat::expect_true(any(vapply(plans, function(x) identical(x$strategy, expected_strategy), logical(1))))
  testthat::expect_equal(nrow(grid$summary_table), 1)
})

test_that("runSearchTuningGrid inherits template search formulas when base options omit one", {
  testthat::skip_on_cran()
  skip_if_tuning_grid_deps()

  tmpl <- make_tuning_template()
  mock_env <- new.env(parent = emptyenv())
  mock_env$null_formula <- NULL

  local_rebind(
    "runFalsePositiveSimulationStudy",
    function(template, search_options, ...) {
      mock_env$null_formula <- search_options$formula
      list(
        generating_scenario = "null",
        per_replicate = data.frame(
          n_inferred_shifts = 0,
          false_positive_rate = 0,
          n_candidates = 5,
          status = "ok"
        ),
        study_summary = list(mean_false_positive_rate = 0)
      )
    },
    environment(runSearchTuningGrid)
  )
  local_rebind(
    "runShiftRecoverySimulationStudy",
    function(template, simulation_options, ...) {
      list(
        generating_scenario = simulation_options$scale_mode,
        per_replicate = data.frame(
          n_candidates = 4,
          n_inferred_shifts = 1,
          status = "ok"
        ),
        evaluation = list(
          strict = list(f1 = 0.4),
          fuzzy = list(f1 = 0.6, recall = 0.7),
          weighted = NULL
        )
      )
    },
    environment(runSearchTuningGrid)
  )

  runSearchTuningGrid(
    template = tmpl,
    IC = "GIC",
    shift_acceptance_thresholds = 5,
    min_descendant_tips_values = 3,
    null_replicates = 1,
    recovery_replicates = 1,
    proportional_simulation_options = list(
      num_shifts = 1,
      min_shift_tips = 2,
      max_shift_tips = 5
    ),
    base_search_options = list(method = "LL", num_cores = 1),
    weighted = FALSE,
    num_cores = 1,
    seed = 64,
    store_studies = FALSE
  )

  testthat::expect_identical(mock_env$null_formula, "trait_data ~ 1")
})

test_that("the tuning grid rejects generator vectors before RNG setup", {
  testthat::skip_on_cran()
  skip_if_tuning_grid_deps()
  tmpl <- make_tuning_template()
  rng_setup_called <- FALSE
  local_rebind(
    ".simulation_set_seed",
    function(...) {
      rng_setup_called <<- TRUE
      stats::runif(1)
      stop("RNG setup reached", call. = FALSE)
    },
    environment(runSearchTuningGrid)
  )
  set.seed(20260709)
  before <- .Random.seed
  testthat::expect_error(
    runSearchTuningGrid(
      template = tmpl,
      IC = "GIC",
      shift_acceptance_thresholds = 5,
      min_descendant_tips_values = 3,
      null_replicates = 1,
      recovery_replicates = 1,
      null_simulation_options = list(
        simulation_generator = c("original", "empirical")
      ),
      proportional_simulation_options = list(
        num_shifts = 1,
        min_shift_tips = 2,
        max_shift_tips = 5
      ),
      base_search_options = list(formula = "trait_data ~ 1", method = "LL"),
      weighted = FALSE,
      num_cores = 1,
      seed = 64,
      store_studies = FALSE
    ),
    "simulation_generator.*original.*empirical"
  )
  testthat::expect_false(rng_setup_called)
  testthat::expect_identical(.Random.seed, before)
})

test_that("runSearchTuningGrid validates integer controls and restores the caller RNG", {
  testthat::skip_on_cran()
  skip_if_tuning_grid_deps()

  tmpl <- make_tuning_template()
  local_rebind(
    "runFalsePositiveSimulationStudy",
    function(...) {
      stats::rnorm(1)
      list(
        generating_scenario = "null",
        per_replicate = data.frame(
          n_inferred_shifts = 0,
          false_positive_rate = 0,
          n_candidates = 5,
          status = "ok"
        ),
        study_summary = list(mean_false_positive_rate = 0)
      )
    },
    environment(runSearchTuningGrid)
  )
  local_rebind(
    "runShiftRecoverySimulationStudy",
    function(template, simulation_options, weighted, ...) {
      stats::rnorm(1)
      list(
        generating_scenario = simulation_options$scale_mode,
        per_replicate = data.frame(
          n_candidates = 4,
          n_inferred_shifts = 1,
          status = "ok"
        ),
        evaluation = list(
          strict = list(f1 = 0.4),
          fuzzy = list(f1 = 0.6, recall = 0.7),
          weighted = if (weighted) list(fuzzy = list(f1 = 0.5)) else NULL
        )
      )
    },
    environment(runSearchTuningGrid)
  )

  run_grid <- function(...) {
    args <- utils::modifyList(list(
      template = tmpl,
      IC = "GIC",
      shift_acceptance_thresholds = 5,
      min_descendant_tips_values = 3,
      null_replicates = 1,
      recovery_replicates = 1,
      proportional_simulation_options = list(
        num_shifts = 1,
        min_shift_tips = 2,
        max_shift_tips = 5
      ),
      base_search_options = list(formula = "trait_data ~ 1", method = "LL"),
      weighted = FALSE,
      num_cores = 1,
      seed = 64,
      store_studies = FALSE
    ), list(...))
    do.call(runSearchTuningGrid, args)
  }

  testthat::expect_error(run_grid(min_descendant_tips_values = 3.5), "min_descendant_tips_values")
  testthat::expect_error(run_grid(null_replicates = 1.5), "null_replicates")
  testthat::expect_error(run_grid(recovery_replicates = 1.5), "recovery_replicates")
  testthat::expect_error(run_grid(tree_tip_count = 3.5), "tree_tip_count")
  testthat::expect_error(run_grid(fuzzy_distance = 1.5), "fuzzy_distance")
  testthat::expect_error(run_grid(num_cores = 1.5), "num_cores")
  testthat::expect_error(run_grid(seed = 1.5), "seed")
  expect_global_rng_restored(run_grid(seed = 65))
})

test_that("runSearchTuningGrid validates inputs", {
  testthat::skip_on_cran()
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
      shift_acceptance_thresholds = c(5, -1),
      min_descendant_tips_values = 3,
      proportional_simulation_options = list(num_shifts = 1, min_shift_tips = 2, max_shift_tips = 5)
    ),
    "finite non-negative"
  )
  testthat::expect_error(
    runSearchTuningGrid(
      template = tmpl,
      IC = "GIC",
      shift_acceptance_thresholds = Inf,
      min_descendant_tips_values = 3,
      proportional_simulation_options = list(num_shifts = 1, min_shift_tips = 2, max_shift_tips = 5)
    ),
    "finite non-negative"
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
      min_descendant_tips_values = NULL,
      proportional_simulation_options = list(num_shifts = 1, min_shift_tips = 2, max_shift_tips = 5)
    ),
    "min_descendant_tips_values"
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
      min_descendant_tips_values = 3,
      tree_tip_count = tmpl$n_tips + 1,
      proportional_simulation_options = list(num_shifts = 1, min_shift_tips = 2, max_shift_tips = 5)
    ),
    "tree_tip_count cannot exceed"
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
  testthat::expect_error(
    runSearchTuningGrid(
      template = tmpl,
      IC = "GIC",
      shift_acceptance_thresholds = 5,
      min_descendant_tips_values = 3,
      null_replicates = 0,
      proportional_simulation_options = list(num_shifts = 1, min_shift_tips = 2, max_shift_tips = 5)
    ),
    "null_replicates"
  )
  testthat::expect_error(
    runSearchTuningGrid(
      template = tmpl,
      IC = "GIC",
      shift_acceptance_thresholds = 5,
      min_descendant_tips_values = 3,
      recovery_replicates = 0,
      proportional_simulation_options = list(num_shifts = 1, min_shift_tips = 2, max_shift_tips = 5)
    ),
    "recovery_replicates"
  )
  testthat::expect_error(
    runSearchTuningGrid(
      template = tmpl,
      IC = "GIC",
      shift_acceptance_thresholds = 5,
      min_descendant_tips_values = 3,
      proportional_simulation_options = list(num_shifts = 1, min_shift_tips = 2, max_shift_tips = 5),
      fuzzy_distance = -1
    ),
    "fuzzy_distance"
  )
  testthat::expect_error(
    runSearchTuningGrid(
      template = tmpl,
      IC = "GIC",
      shift_acceptance_thresholds = 5,
      min_descendant_tips_values = 3,
      proportional_simulation_options = list(num_shifts = 1, min_shift_tips = 2, max_shift_tips = 5),
      weighted = NA
    ),
    "weighted"
  )
  testthat::expect_error(
    runSearchTuningGrid(
      template = tmpl,
      IC = "GIC",
      shift_acceptance_thresholds = 5,
      min_descendant_tips_values = 3,
      proportional_simulation_options = list(num_shifts = 1, min_shift_tips = 2, max_shift_tips = 5),
      num_cores = 0
    ),
    "num_cores"
  )
  testthat::expect_error(
    runSearchTuningGrid(
      template = tmpl,
      IC = "GIC",
      shift_acceptance_thresholds = 5,
      min_descendant_tips_values = 3,
      proportional_simulation_options = list(num_shifts = 1, min_shift_tips = 2, max_shift_tips = 5),
      store_studies = NA
    ),
    "store_studies"
  )
})

test_that("runSearchTuningGrid can parallelize over settings in a small real run", {
  testthat::skip_on_cran()
  skip_if_tuning_grid_deps()
  testthat::skip_on_covr()

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
