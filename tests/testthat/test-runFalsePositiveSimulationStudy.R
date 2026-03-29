testthat::skip_on_cran()

skip_if_fp_study_deps <- function() {
  testthat::skip_if_not_installed("ape")
  testthat::skip_if_not_installed("phytools")
  testthat::skip_if_not_installed("mvMORPH")
  testthat::skip_if_not_installed("future")
  testthat::skip_if_not_installed("future.apply")
  testthat::skip_if_not_installed("progressr")
}

make_fp_template <- function() {
  set.seed(30)
  tr <- ape::rtree(22)
  X <- matrix(rnorm(22 * 2), ncol = 2)
  rownames(X) <- tr$tip.label
  createSimulationTemplate(
    baseline_tree = tr,
    trait_data = X,
    formula = "trait_data ~ 1",
    method = "LL"
  )
}

make_formula_fp_template <- function() {
  set.seed(31)
  tr <- ape::rtree(22)
  grp <- factor(rep(c("a", "b"), length.out = 22))
  size <- rnorm(22)
  X <- data.frame(
    y1 = 0.4 * size + ifelse(grp == "b", 0.8, 0) + rnorm(22),
    y2 = -0.2 * size + ifelse(grp == "b", -0.5, 0) + rnorm(22),
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

make_fp_template_with_error <- function() {
  set.seed(32)
  tr <- ape::rtree(18)
  X <- matrix(rnorm(18 * 2), ncol = 2)
  rownames(X) <- tr$tip.label
  createSimulationTemplate(
    baseline_tree = tr,
    trait_data = X,
    formula = "trait_data ~ 1",
    method = "LL",
    error = TRUE
  )
}

test_that("runFalsePositiveSimulationStudy returns study summaries", {
  skip_if_fp_study_deps()

  tmpl <- make_fp_template()
  study <- runFalsePositiveSimulationStudy(
    tmpl,
    n_replicates = 2,
    tree_tip_count = 18,
    search_options = list(
      formula = "trait_data ~ 1",
      min_descendant_tips = 3,
      shift_acceptance_threshold = 5,
      num_cores = 1,
      IC = "GIC",
      method = "LL"
    ),
    num_cores = 1,
    seed = 6
  )

  testthat::expect_s3_class(study, "bifrost_simulation_study")
  testthat::expect_identical(study$study_type, "false_positive")
  testthat::expect_identical(study$generating_scenario, "null")
  testthat::expect_length(study$simdata, 2)
  testthat::expect_length(study$results, 2)
  testthat::expect_true(all(c(
    "replicate", "n_candidates", "n_inferred_shifts", "false_positive_rate",
    "status", "error"
  ) %in% names(study$per_replicate)))
  testthat::expect_equal(
    study$per_replicate$false_positive_rate,
    study$per_replicate$n_inferred_shifts / study$per_replicate$n_candidates
  )
})

test_that("runFalsePositiveSimulationStudy rejects per-replicate seeds", {
  skip_if_fp_study_deps()

  tmpl <- make_fp_template()

  testthat::expect_error(
    runFalsePositiveSimulationStudy(
      tmpl,
      n_replicates = 2,
      simulation_options = list(seed = 99),
      search_options = list(
        formula = "trait_data ~ 1",
        min_descendant_tips = 3,
        shift_acceptance_threshold = 5,
        num_cores = 1,
        IC = "GIC",
        method = "LL"
      ),
      num_cores = 1,
      seed = 6
    ),
    "simulation_options\\$seed"
  )
})

test_that("runFalsePositiveSimulationStudy uses intercept-only searches for formula-calibrated templates", {
  skip_if_fp_study_deps()

  tmpl <- make_formula_fp_template()
  study <- runFalsePositiveSimulationStudy(
    tmpl,
    n_replicates = 1,
    tree_tip_count = 20,
    search_options = list(
      formula = "trait_data ~ 1",
      min_descendant_tips = 3,
      shift_acceptance_threshold = 5,
      num_cores = 1,
      IC = "GIC",
      method = "LL"
    ),
    num_cores = 1,
    seed = 123
  )

  testthat::expect_identical(study$search_options$formula, "trait_data ~ 1")
  testthat::expect_identical(study$per_replicate$status, "ok")
  testthat::expect_equal(colnames(study$simdata[[1]]$data), c("y1", "y2"))
  testthat::expect_false("mass" %in% colnames(study$simdata[[1]]$data))
})

test_that("runFalsePositiveSimulationStudy inherits the default intercept-only search for data.frame templates", {
  skip_if_fp_study_deps()

  tmpl <- make_formula_fp_template()
  study <- runFalsePositiveSimulationStudy(
    tmpl,
    n_replicates = 1,
    tree_tip_count = 20,
    search_options = list(
      min_descendant_tips = 3,
      shift_acceptance_threshold = 5,
      num_cores = 1,
      IC = "GIC",
      method = "LL"
    ),
    num_cores = 1,
    seed = 124
  )

  testthat::expect_identical(study$search_options$formula, "trait_data ~ 1")
  testthat::expect_identical(study$per_replicate$status, "ok")
  testthat::expect_true(all(is.na(study$per_replicate$error)))
})

test_that("runFalsePositiveSimulationStudy rejects non-intercept search formulas", {
  skip_if_fp_study_deps()

  tmpl <- make_formula_fp_template()

  testthat::expect_error(
    runFalsePositiveSimulationStudy(
      tmpl,
      n_replicates = 1,
      tree_tip_count = 20,
      search_options = list(
        formula = "trait_data[, 1:2] ~ trait_data[, 3]",
        min_descendant_tips = 3,
        shift_acceptance_threshold = 5,
        num_cores = 1,
        IC = "GIC",
        method = "LL"
      ),
      num_cores = 1,
      seed = 123
    ),
    "intercept-only search formulas only"
  )
})

test_that("runFalsePositiveSimulationStudy validates inputs and inherits template error settings", {
  skip_if_fp_study_deps()

  tmpl <- make_fp_template_with_error()

  testthat::expect_error(
    runFalsePositiveSimulationStudy(list(), n_replicates = 1),
    "bifrost_simulation_template"
  )
  testthat::expect_error(
    runFalsePositiveSimulationStudy(tmpl, n_replicates = 0),
    "n_replicates"
  )
  testthat::expect_error(
    runFalsePositiveSimulationStudy(tmpl, n_replicates = 1, num_cores = 0),
    "num_cores"
  )
  testthat::expect_error(
    runFalsePositiveSimulationStudy(tmpl, n_replicates = 1, simulation_options = 1),
    "must both be lists"
  )

  study <- runFalsePositiveSimulationStudy(
    tmpl,
    n_replicates = 1,
    search_options = list(
      formula = "trait_data ~ 1",
      min_descendant_tips = 3,
      shift_acceptance_threshold = 5,
      num_cores = 1,
      IC = "GIC",
      method = "LL"
    ),
    num_cores = 1,
    seed = 33
  )

  testthat::expect_true(isTRUE(study$search_options$error))
})

test_that("runFalsePositiveSimulationStudy inherits evaluated template settings", {
  skip_if_fp_study_deps()

  set.seed(36)
  tr <- ape::rtree(18)
  X <- matrix(rnorm(18 * 2), ncol = 2)
  rownames(X) <- tr$tip.label
  method_value <- "LL"
  error_value <- TRUE

  tmpl <- createSimulationTemplate(
    baseline_tree = tr,
    trait_data = X,
    formula = "trait_data ~ 1",
    method = method_value,
    error = error_value
  )

  study <- runFalsePositiveSimulationStudy(
    tmpl,
    n_replicates = 1,
    search_options = list(
      formula = "trait_data ~ 1",
      min_descendant_tips = 3,
      shift_acceptance_threshold = 5,
      num_cores = 1,
      IC = "GIC"
    ),
    num_cores = 1,
    seed = 37
  )

  testthat::expect_identical(study$search_options$method, "LL")
  testthat::expect_true(isTRUE(study$search_options$error))
  testthat::expect_identical(study$per_replicate$status, "ok")
})

test_that("runFalsePositiveSimulationStudy falls back to template call settings when stored values are missing", {
  skip_if_fp_study_deps()

  method_value <- "LL"
  tmpl <- make_fp_template_with_error()
  tmpl$fit_method <- NULL
  tmpl$fit_error <- NULL
  tmpl$global_model$call$method <- method_value
  tmpl$global_model$call$error <- quote(TRUE)
  tmpl$search_formula <- NULL

  study <- runFalsePositiveSimulationStudy(
    tmpl,
    n_replicates = 1,
    search_options = list(
      min_descendant_tips = 3,
      shift_acceptance_threshold = 5,
      num_cores = 1,
      IC = "GIC"
    ),
    num_cores = 1,
    seed = 39
  )

  testthat::expect_identical(study$search_options$formula, "trait_data ~ 1")
  testthat::expect_identical(study$search_options$method, "LL")
  testthat::expect_true(isTRUE(study$search_options$error))
})

test_that("runFalsePositiveSimulationStudy tolerates unevaluable template error calls", {
  skip_if_fp_study_deps()

  tmpl <- make_fp_template()
  tmpl$fit_method <- NULL
  tmpl$fit_error <- NULL
  tmpl$global_model$call$method <- "LL"
  tmpl$global_model$call$error <- quote(missing_error_flag)
  tmpl$search_formula <- NULL

  study <- runFalsePositiveSimulationStudy(
    tmpl,
    n_replicates = 1,
    search_options = list(
      min_descendant_tips = 3,
      shift_acceptance_threshold = 5,
      num_cores = 1,
      IC = "GIC"
    ),
    num_cores = 1,
    seed = 40
  )

  testthat::expect_false("error" %in% names(study$search_options))
  testthat::expect_identical(study$per_replicate$status, "ok")
})

test_that("runFalsePositiveSimulationStudy is reproducible with the same seed and workers", {
  skip_if_fp_study_deps()

  tmpl <- make_fp_template()
  run_once <- function() {
    runFalsePositiveSimulationStudy(
      tmpl,
      n_replicates = 2,
      tree_tip_count = 18,
      search_options = list(
        formula = "trait_data ~ 1",
        min_descendant_tips = 3,
        shift_acceptance_threshold = 5,
        num_cores = 1,
        IC = "GIC",
        method = "LL"
      ),
      num_cores = 2,
      seed = 6
    )
  }

  study_a <- run_once()
  study_b <- run_once()

  testthat::expect_identical(study_a$per_replicate, study_b$per_replicate)
  testthat::expect_true(all(vapply(seq_along(study_a$simdata), function(i) {
    isTRUE(all.equal(study_a$simdata[[i]]$data, study_b$simdata[[i]]$data))
  }, logical(1))))
  testthat::expect_identical(
    lapply(study_a$results, `[[`, "shift_nodes_no_uncertainty"),
    lapply(study_b$results, `[[`, "shift_nodes_no_uncertainty")
  )
})

test_that("runFalsePositiveSimulationStudy warns on nested parallelism", {
  skip_if_fp_study_deps()

  tmpl <- make_fp_template()

  testthat::local_mocked_bindings(
    future_lapply = function(X, FUN, ...) lapply(X, FUN),
    .package = "future.apply"
  )
  testthat::local_mocked_bindings(
    plan = function(...) list(),
    .package = "future"
  )

  testthat::expect_warning(
    runFalsePositiveSimulationStudy(
      tmpl,
      n_replicates = 1,
      search_options = list(
        formula = "trait_data ~ 1",
        min_descendant_tips = 3,
        shift_acceptance_threshold = 5,
        num_cores = 2,
        IC = "GIC",
        method = "LL"
      ),
      num_cores = 2,
      seed = 34
    ),
    "nested parallelism"
  )
})

test_that("runFalsePositiveSimulationStudy validates formula types and can choose multisession plans", {
  skip_if_fp_study_deps()

  tmpl <- make_fp_template()

  testthat::expect_error(
    runFalsePositiveSimulationStudy(
      tmpl,
      n_replicates = 1,
      search_options = list(
        formula = 1,
        min_descendant_tips = 3,
        shift_acceptance_threshold = 5,
        num_cores = 1,
        IC = "GIC"
      ),
      num_cores = 1,
      seed = 41
    ),
    "single character string or a formula object"
  )

  plans <- list()
  old_rstudio <- Sys.getenv("RSTUDIO", unset = NA_character_)
  old_rstudio_init <- Sys.getenv("RSTUDIO_SESSION_INITIALIZED", unset = NA_character_)
  Sys.setenv(RSTUDIO = "1", RSTUDIO_SESSION_INITIALIZED = "1")
  on.exit({
    if (is.na(old_rstudio)) Sys.unsetenv("RSTUDIO") else Sys.setenv(RSTUDIO = old_rstudio)
    if (is.na(old_rstudio_init)) Sys.unsetenv("RSTUDIO_SESSION_INITIALIZED") else Sys.setenv(RSTUDIO_SESSION_INITIALIZED = old_rstudio_init)
  }, add = TRUE)

  testthat::local_mocked_bindings(
    future_lapply = function(X, FUN, ...) lapply(X, FUN),
    .package = "future.apply"
  )
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

  runFalsePositiveSimulationStudy(
    tmpl,
    n_replicates = 1,
    search_options = list(
      formula = "trait_data ~ 1",
      min_descendant_tips = 3,
      shift_acceptance_threshold = 5,
      num_cores = 1,
      IC = "GIC",
      method = "LL"
    ),
    num_cores = 2,
    seed = 42
  )

  testthat::expect_true(any(vapply(plans, function(x) identical(x$strategy, future::multisession), logical(1))))
})

test_that("runFalsePositiveSimulationStudy accepts formula objects in search options", {
  skip_if_fp_study_deps()

  tmpl <- make_fp_template()
  study <- runFalsePositiveSimulationStudy(
    tmpl,
    n_replicates = 1,
    search_options = list(
      formula = stats::as.formula("trait_data ~ 1"),
      min_descendant_tips = 3,
      shift_acceptance_threshold = 5,
      num_cores = 1,
      IC = "GIC",
      method = "LL"
    ),
    num_cores = 1,
    seed = 43
  )

  testthat::expect_s3_class(study$search_options$formula, "formula")
  testthat::expect_identical(study$per_replicate$status, "ok")
})

test_that("runFalsePositiveSimulationStudy records search fallback errors", {
  skip_if_fp_study_deps()

  tmpl <- make_fp_template()

  study <- runFalsePositiveSimulationStudy(
    tmpl,
    n_replicates = 1,
    search_options = list(
      formula = "trait_data ~ 1",
      min_descendant_tips = 3,
      shift_acceptance_threshold = 5,
      num_cores = 1,
      IC = "BAD",
      method = "LL"
    ),
    num_cores = 1,
    seed = 35
  )

  testthat::expect_identical(study$per_replicate$status, "error")
  testthat::expect_match(study$per_replicate$error, "IC must be GIC or BIC")
  testthat::expect_equal(
    study$results[[1]]$num_candidates,
    max(length(generatePaintedTrees(ape::as.phylo(study$simdata[[1]]$tree), min_tips = 3)) - 1L, 0L)
  )
})

test_that("runFalsePositiveSimulationStudy reports NA summaries when no replicate is evaluable", {
  skip_if_fp_study_deps()

  tmpl <- make_fp_template()

  study <- runFalsePositiveSimulationStudy(
    tmpl,
    n_replicates = 1,
    search_options = list(
      formula = "trait_data ~ 1",
      min_descendant_tips = 100,
      shift_acceptance_threshold = 5,
      num_cores = 1,
      IC = "GIC",
      method = "LL"
    ),
    num_cores = 1,
    seed = 38
  )

  testthat::expect_equal(study$per_replicate$n_candidates, 0)
  testthat::expect_true(is.na(study$per_replicate$false_positive_rate))
  testthat::expect_identical(study$study_summary$n_evaluable_replicates, 0L)
  testthat::expect_true(is.na(study$study_summary$mean_false_positive_rate))
  testthat::expect_true(is.na(study$study_summary$median_false_positive_rate))
})
