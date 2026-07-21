skip_if_fp_study_deps <- function() {
  testthat::skip_if_not_installed("ape")
  testthat::skip_if_not_installed("phytools")
  testthat::skip_if_not_installed("mvMORPH")
  testthat::skip_if_not_installed("future")
  testthat::skip_if_not_installed("future.apply")
  testthat::skip_if_not_installed("progressr")
  withr::local_options(list(progressr.enable = FALSE))
}

.fp_template_cache <- new.env(parent = emptyenv())

get_cached_fp_template <- function(key, builder) {
  if (!exists(key, envir = .fp_template_cache, inherits = FALSE)) {
    assign(key, builder(), envir = .fp_template_cache)
  }
  unserialize(serialize(get(key, envir = .fp_template_cache, inherits = FALSE), NULL))
}

make_fp_template <- function() {
  get_cached_fp_template("fp_template", function() {
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
  })
}

make_formula_fp_template <- function() {
  get_cached_fp_template("formula_fp_template", function() {
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
  })
}

make_fp_template_with_error <- function() {
  get_cached_fp_template("fp_template_with_error", function() {
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
  })
}

test_that("runFalsePositiveSimulationStudy returns study summaries", {
  testthat::skip_on_cran()
  skip_if_fp_study_deps()

  tmpl <- make_fp_template()
  study <- runFalsePositiveSimulationStudy(
    tmpl,
    n_replicates = 1,
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
  testthat::expect_identical(study$simulation_generator, "original")
  testthat::expect_identical(
    unique(study$per_replicate$simulation_generator),
    "original"
  )
  testthat::expect_false("simulation_design" %in% names(study))
  testthat::expect_false("simulation_design" %in% names(study$per_replicate))
  testthat::expect_false(study$search_options$progress)
  testthat::expect_length(study$simdata, 1)
  testthat::expect_length(study$results, 1)
  testthat::expect_true(all(c(
    "replicate", "n_candidates", "n_inferred_shifts", "false_positive_rate",
    "status", "error"
  ) %in% names(study$per_replicate)))
  testthat::expect_equal(
    study$per_replicate$false_positive_rate,
    study$per_replicate$n_inferred_shifts / study$per_replicate$n_candidates
  )
  testthat::expect_identical(as.data.frame(study), study$per_replicate)
})

test_that("runFalsePositiveSimulationStudy forwards the empirical generator", {
  testthat::skip_on_cran()
  skip_if_fp_study_deps()

  study <- runFalsePositiveSimulationStudy(
    make_fp_template(),
    n_replicates = 1,
    tree_tip_count = 18,
    simulation_options = list(simulation_generator = "empirical"),
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

  testthat::expect_identical(study$simulation_generator, "empirical")
  testthat::expect_identical(
    study$simdata[[1L]]$simulation_generator,
    "empirical"
  )
})

test_that("runFalsePositiveSimulationStudy stores compact replicate calls", {
  testthat::skip_on_cran()
  skip_if_fp_study_deps()

  study <- runFalsePositiveSimulationStudy(
    make_fp_template(),
    n_replicates = 1,
    tree_tip_count = 18,
    search_options = list(
      formula = "trait_data ~ 1",
      min_descendant_tips = 3,
      shift_acceptance_threshold = 5,
      num_cores = 1,
      IC = "BAD",
      method = "LL"
    ),
    num_cores = 1,
    seed = 56
  )

  replicate_call <- study$simdata[[1L]]$user_input
  testthat::expect_identical(replicate_call[[1L]], quote(simulateNullDataset))
  testthat::expect_identical(
    replicate_call$template,
    "<bifrost_simulation_template>"
  )
  testthat::expect_false(is.symbol(replicate_call$template))
  testthat::expect_lt(as.numeric(object.size(replicate_call)), 10000)
})

test_that("the false-positive study rejects generator vectors before RNG setup", {
  skip_if_fp_study_deps()
  tmpl <- make_fp_template()
  rng_setup_called <- FALSE
  local_rebind(
    ".simulation_set_seed",
    function(...) {
      rng_setup_called <<- TRUE
      stats::runif(1)
      stop("RNG setup reached", call. = FALSE)
    },
    environment(runFalsePositiveSimulationStudy)
  )
  set.seed(20260709)
  before <- .Random.seed
  testthat::expect_error(
    runFalsePositiveSimulationStudy(
      tmpl,
      n_replicates = 1,
      simulation_options = list(
        simulation_generator = c("original", "empirical")
      ),
      seed = 6
    ),
    "simulation_generator.*original.*empirical"
  )
  testthat::expect_false(rng_setup_called)
  testthat::expect_identical(.Random.seed, before)
})

test_that("runFalsePositiveSimulationStudy rejects per-replicate seeds", {
  testthat::skip_on_cran()
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
  testthat::skip_on_cran()
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

  testthat::expect_identical(
    study$search_options$formula,
    simulationStudyFormulaCases$string
  )
  testthat::expect_identical(study$per_replicate$status, "ok")
  testthat::expect_equal(colnames(study$simdata[[1]]$data), c("y1", "y2"))
  testthat::expect_false("mass" %in% colnames(study$simdata[[1]]$data))
})

test_that("runFalsePositiveSimulationStudy inherits the default intercept-only search for data.frame templates", {
  testthat::skip_on_cran()
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
  testthat::skip_on_cran()
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
    simulationStudyFormulaErrors$non_intercept
  )
})

test_that("runFalsePositiveSimulationStudy validates inputs and inherits template error settings", {
  testthat::skip_on_cran()
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
    runFalsePositiveSimulationStudy(tmpl, n_replicates = 1.5),
    "n_replicates"
  )
  testthat::expect_error(
    runFalsePositiveSimulationStudy(tmpl, n_replicates = 1, tree_tip_count = 3.5),
    "tree_tip_count"
  )
  testthat::expect_error(
    runFalsePositiveSimulationStudy(tmpl, n_replicates = 1, tree_tip_count = tmpl$n_tips + 1L),
    "tree_tip_count cannot exceed"
  )
  testthat::expect_error(
    runFalsePositiveSimulationStudy(tmpl, n_replicates = 1, num_cores = 0),
    "num_cores"
  )
  testthat::expect_error(
    runFalsePositiveSimulationStudy(tmpl, n_replicates = 1, num_cores = 1.5),
    "num_cores"
  )
  testthat::expect_error(
    runFalsePositiveSimulationStudy(tmpl, n_replicates = 1, seed = 1.5),
    "seed"
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

test_that("runFalsePositiveSimulationStudy restores the caller RNG state when seeded", {
  testthat::skip_on_cran()
  skip_if_fp_study_deps()

  tmpl <- make_fp_template()
  local_rebind(
    "simulateNullDataset",
    function(template, tree_tip_count = NULL, ...) {
      stats::rnorm(1)
      list(
        tree = template$baseline_tree,
        data = template$trait_data
      )
    },
    environment(runFalsePositiveSimulationStudy)
  )
  local_rebind(
    "searchOptimalConfiguration",
    function(...) {
      stats::rnorm(1)
      list(
        shift_nodes_no_uncertainty = integer(0),
        num_candidates = 3L,
        ic_weights = data.frame()
      )
    },
    environment(runFalsePositiveSimulationStudy)
  )

  expect_global_rng_restored(
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
      num_cores = 1,
      seed = 6
    )
  )
})

test_that("runFalsePositiveSimulationStudy inherits evaluated template settings", {
  testthat::skip_on_cran()
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
  testthat::skip_on_cran()
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
  testthat::skip_on_cran()
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
  testthat::skip_on_cran()
  skip_if_fp_study_deps()
  testthat::skip_on_covr()

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

test_that("runFalsePositiveSimulationStudy gives seeded draws independent of worker count", {
  testthat::skip_on_cran()
  skip_if_fp_study_deps()
  testthat::skip_on_covr()

  tmpl <- make_fp_template()
  local_rebind(
    "simulateNullDataset",
    function(template, tree_tip_count = NULL, ...) {
      draw <- stats::rnorm(1)
      list(
        tree = template$baseline_tree,
        data = template$trait_data,
        random_draw = draw
      )
    },
    environment(runFalsePositiveSimulationStudy)
  )
  local_rebind(
    "searchOptimalConfiguration",
    function(...) {
      list(
        shift_nodes_no_uncertainty = integer(0),
        num_candidates = 3L,
        ic_weights = data.frame(),
        random_draw = stats::rnorm(1)
      )
    },
    environment(runFalsePositiveSimulationStudy)
  )
  withr::local_envvar(c(
    RSTUDIO = "1",
    RSTUDIO_SESSION_INITIALIZED = "1"
  ))

  run_once <- function(num_cores) {
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
      num_cores = num_cores,
      seed = 606
    )
  }

  serial <- run_once(1L)
  multisession <- run_once(2L)

  testthat::expect_identical(
    vapply(serial$simdata, `[[`, numeric(1), "random_draw"),
    vapply(multisession$simdata, `[[`, numeric(1), "random_draw")
  )
  testthat::expect_identical(
    vapply(serial$results, `[[`, numeric(1), "random_draw"),
    vapply(multisession$results, `[[`, numeric(1), "random_draw")
  )
})

test_that("runFalsePositiveSimulationStudy makes outer parallelism override inner cores", {
  testthat::skip_on_cran()
  skip_if_fp_study_deps()

  tmpl <- make_fp_template()

  testthat::local_mocked_bindings(
    future_lapply = serial_future_lapply,
    .package = "future.apply"
  )
  testthat::local_mocked_bindings(
    plan = function(...) list(),
    .package = "future"
  )

  testthat::expect_warning(
    study <- runFalsePositiveSimulationStudy(
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
    "forcing each search to use one core"
  )
  testthat::expect_identical(study$search_options$num_cores, 1L)
})

test_that("runFalsePositiveSimulationStudy does not capture simdata in future workers", {
  testthat::skip_on_cran()
  skip_if_fp_study_deps()

  worker_environments <- list()
  testthat::local_mocked_bindings(
    future_lapply = function(X, FUN, ...) {
      worker_environments[[length(worker_environments) + 1L]] <<- environment(FUN)
      serial_future_lapply(X, FUN, ...)
    },
    .package = "future.apply"
  )
  testthat::local_mocked_bindings(
    plan = function(...) list(),
    .package = "future"
  )
  local_rebind(
    "searchOptimalConfiguration",
    function(...) {
      list(
        shift_nodes_no_uncertainty = integer(0),
        num_candidates = 3L,
        ic_weights = data.frame()
      )
    },
    environment(runFalsePositiveSimulationStudy)
  )

  runFalsePositiveSimulationStudy(
    make_fp_template(),
    n_replicates = 1,
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
    seed = 54
  )

  testthat::expect_length(worker_environments, 1L)
  testthat::expect_false(any(vapply(worker_environments, function(worker_env) {
    exists("simdata", envir = worker_env, inherits = FALSE)
  }, logical(1L))))
})

test_that("runFalsePositiveSimulationStudy validates formula types and can choose multisession plans", {
  testthat::skip_on_cran()
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
    simulationStudyFormulaErrors$invalid_type
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
    future_lapply = serial_future_lapply,
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
  testthat::skip_on_cran()
  skip_if_fp_study_deps()

  tmpl <- make_fp_template()
  study <- runFalsePositiveSimulationStudy(
    tmpl,
    n_replicates = 1,
    search_options = list(
      formula = simulationStudyFormulaCases$formula,
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
  testthat::expect_identical(
    study$search_options$formula,
    simulationStudyFormulaCases$formula
  )
  testthat::expect_identical(study$per_replicate$status, "ok")
})

test_that("runFalsePositiveSimulationStudy records search fallback errors", {
  testthat::skip_on_cran()
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
  testthat::expect_true(is.na(study$per_replicate$n_inferred_shifts))
  testthat::expect_true(is.na(study$per_replicate$false_positive_rate))
  testthat::expect_identical(study$study_summary$n_completed, 0L)
  testthat::expect_identical(study$study_summary$n_failed, 1L)
  testthat::expect_identical(study$study_summary$n_evaluable_replicates, 0L)
  testthat::expect_equal(study$study_summary$completion_rate, 0)
  testthat::expect_equal(study$study_summary$failure_rate, 1)
  testthat::expect_true(is.na(study$study_summary$mean_false_positive_rate))
  testthat::expect_true(is.na(study$study_summary$median_false_positive_rate))
  testthat::expect_equal(
    study$results[[1]]$num_candidates,
    max(length(generatePaintedTrees(ape::as.phylo(study$simdata[[1]]$tree), min_tips = 3)) - 1L, 0L)
  )
})

test_that("runFalsePositiveSimulationStudy reports NA summaries when no replicate is evaluable", {
  testthat::skip_on_cran()
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
