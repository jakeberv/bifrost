skip_if_recovery_study_deps <- function() {
  testthat::skip_if_not_installed("ape")
  testthat::skip_if_not_installed("phytools")
  testthat::skip_if_not_installed("mvMORPH")
  testthat::skip_if_not_installed("RRphylo")
  testthat::skip_if_not_installed("future")
  testthat::skip_if_not_installed("future.apply")
  testthat::skip_if_not_installed("progressr")
  withr::local_options(list(progressr.enable = FALSE))
}

.recovery_template_cache <- new.env(parent = emptyenv())

get_cached_recovery_template <- function(key, builder) {
  if (!exists(key, envir = .recovery_template_cache, inherits = FALSE)) {
    assign(key, builder(), envir = .recovery_template_cache)
  }
  unserialize(serialize(get(key, envir = .recovery_template_cache, inherits = FALSE), NULL))
}

make_recovery_template <- function() {
  get_cached_recovery_template("recovery_template", function() {
    set.seed(40)
    tr <- ape::rtree(24)
    X <- matrix(rnorm(24 * 2), ncol = 2)
    rownames(X) <- tr$tip.label
    createSimulationTemplate(
      baseline_tree = tr,
      trait_data = X,
      formula = "trait_data ~ 1",
      method = "LL"
    )
  })
}

make_formula_recovery_template <- function() {
  get_cached_recovery_template("formula_recovery_template", function() {
    set.seed(41)
    tr <- ape::rtree(24)
    grp <- factor(rep(c("a", "b"), each = 12))
    size <- rnorm(24)
    X <- data.frame(
      y1 = 0.5 * size + ifelse(grp == "b", 1, 0) + rnorm(24),
      y2 = -0.3 * size + ifelse(grp == "b", -0.5, 0) + rnorm(24),
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

make_recovery_template_with_error <- function() {
  get_cached_recovery_template("recovery_template_with_error", function() {
    set.seed(42)
    tr <- ape::rtree(20)
    X <- matrix(rnorm(20 * 2), ncol = 2)
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

test_that("runShiftRecoverySimulationStudy supports the proportional scenario", {
  testthat::skip_on_cran()
  skip_if_recovery_study_deps()

  tmpl <- make_recovery_template()
  study <- suppressWarnings(
    runShiftRecoverySimulationStudy(
      tmpl,
      n_replicates = 1,
      tree_tip_count = 20,
      simulation_options = list(
        num_shifts = 2,
        min_shift_tips = 3,
        max_shift_tips = 7,
        scale_mode = "proportional"
      ),
      search_options = list(
        formula = "trait_data ~ 1",
        min_descendant_tips = 3,
        shift_acceptance_threshold = 5,
        num_cores = 1,
        IC = "GIC",
        method = "LL",
        uncertaintyweights_par = FALSE
      ),
      weighted = FALSE,
      num_cores = 1,
      seed = 7
    )
  )

  testthat::expect_s3_class(study, "bifrost_simulation_study")
  testthat::expect_identical(study$study_type, "shift_recovery")
  testthat::expect_identical(study$generating_scenario, "proportional")
  testthat::expect_identical(study$simulation_generator, "original")
  testthat::expect_identical(
    unique(study$per_replicate$simulation_generator),
    "original"
  )
  testthat::expect_false("simulation_design" %in% names(study))
  testthat::expect_false("simulation_design" %in% names(study$per_replicate))
  testthat::expect_false(study$search_options$progress)
  testthat::expect_length(study$simdata, 1)
  testthat::expect_true(all(vapply(study$simdata, function(x) length(x$shiftNodes), integer(1)) == 2L))
  testthat::expect_false(any(vapply(study$simdata, function(sim) {
    any(vapply(seq_along(sim$shiftNodes), function(i) {
      any(vapply(seq_along(sim$shiftNodes), function(j) {
        if (i == j) {
          return(FALSE)
        }
        sim$shiftNodes[j] %in% getDescendants(sim$paintedTree, sim$shiftNodes[i])
      }, logical(1)))
    }, logical(1)))
  }, logical(1))))
  testthat::expect_true(all(c("strict", "fuzzy", "weighted", "counts") %in% names(study$evaluation)))
  testthat::expect_identical(as.data.frame(study), study$per_replicate)
})

test_that("runShiftRecoverySimulationStudy supports the correlation scenario", {
  testthat::skip_on_cran()
  skip_if_recovery_study_deps()

  tmpl <- make_recovery_template()
  study <- suppressWarnings(
    runShiftRecoverySimulationStudy(
      tmpl,
      n_replicates = 1,
      tree_tip_count = 20,
      simulation_options = list(
        num_shifts = 2,
        min_shift_tips = 3,
        max_shift_tips = 7,
        scale_mode = "correlation",
        simulation_generator = "empirical"
      ),
      search_options = list(
        formula = "trait_data ~ 1",
        min_descendant_tips = 3,
        shift_acceptance_threshold = 5,
        num_cores = 1,
        IC = "GIC",
        method = "LL",
        uncertaintyweights_par = FALSE
      ),
      weighted = FALSE,
      num_cores = 1,
      seed = 8
    )
  )

  testthat::expect_identical(study$generating_scenario, "correlation")
  testthat::expect_equal(study$per_replicate$generating_scenario, "correlation")
  testthat::expect_identical(study$simulation_generator, "empirical")
  testthat::expect_identical(
    study$per_replicate$simulation_generator,
    "empirical"
  )
})

test_that("runShiftRecoverySimulationStudy stores compact replicate calls", {
  testthat::skip_on_cran()
  skip_if_recovery_study_deps()

  study <- runShiftRecoverySimulationStudy(
    make_recovery_template(),
    n_replicates = 1,
    tree_tip_count = 20,
    simulation_options = list(
      num_shifts = 1,
      min_shift_tips = 3,
      max_shift_tips = 7,
      scale_mode = "proportional"
    ),
    search_options = list(
      formula = "trait_data ~ 1",
      min_descendant_tips = 3,
      shift_acceptance_threshold = 5,
      num_cores = 1,
      IC = "BAD",
      method = "LL",
      uncertaintyweights_par = FALSE
    ),
    weighted = FALSE,
    num_cores = 1,
    seed = 56
  )

  replicate_call <- study$simdata[[1L]]$user_input
  testthat::expect_identical(replicate_call[[1L]], quote(simulateShiftedDataset))
  testthat::expect_identical(names(replicate_call), c("", "num_shifts", "min_shift_tips", "max_shift_tips", "scale_mode", "template", "tree_tip_count"))
  testthat::expect_identical(
    replicate_call$template,
    "<bifrost_simulation_template>"
  )
  testthat::expect_false(is.symbol(replicate_call$template))
  testthat::expect_lt(as.numeric(object.size(replicate_call)), 10000)
})

test_that("the recovery study rejects generator vectors before RNG setup", {
  skip_if_recovery_study_deps()
  tmpl <- make_recovery_template()
  rng_setup_called <- FALSE
  local_rebind(
    ".simulation_set_seed",
    function(...) {
      rng_setup_called <<- TRUE
      stats::runif(1)
      stop("RNG setup reached", call. = FALSE)
    },
    environment(runShiftRecoverySimulationStudy)
  )
  set.seed(20260709)
  before <- .Random.seed
  testthat::expect_error(
    runShiftRecoverySimulationStudy(
      tmpl,
      n_replicates = 1,
      simulation_options = list(
        num_shifts = 1,
        min_shift_tips = 2,
        max_shift_tips = 5,
        simulation_generator = c("original", "empirical")
      ),
      seed = 43
    ),
    "simulation_generator.*original.*empirical"
  )
  testthat::expect_false(rng_setup_called)
  testthat::expect_identical(.Random.seed, before)
})

test_that("runShiftRecoverySimulationStudy rejects per-replicate seeds", {
  testthat::skip_on_cran()
  skip_if_recovery_study_deps()

  tmpl <- make_recovery_template()

  testthat::expect_error(
    runShiftRecoverySimulationStudy(
      tmpl,
      n_replicates = 1,
      simulation_options = list(
        num_shifts = 2,
        min_shift_tips = 3,
        max_shift_tips = 7,
        scale_mode = "proportional",
        seed = 99
      ),
      search_options = list(
        formula = "trait_data ~ 1",
        min_descendant_tips = 3,
        shift_acceptance_threshold = 5,
        num_cores = 1,
        IC = "GIC",
        method = "LL",
        uncertaintyweights_par = FALSE
      ),
      weighted = FALSE,
      num_cores = 1,
      seed = 8
    ),
    "simulation_options\\$seed"
  )
})

test_that("runShiftRecoverySimulationStudy uses intercept-only searches for formula-calibrated templates", {
  testthat::skip_on_cran()
  skip_if_recovery_study_deps()

  tmpl <- make_formula_recovery_template()
  study <- suppressWarnings(
    runShiftRecoverySimulationStudy(
      tmpl,
      n_replicates = 1,
      tree_tip_count = 20,
      simulation_options = list(
        num_shifts = 1,
        min_shift_tips = 3,
        max_shift_tips = 7,
        scale_mode = "proportional"
      ),
      search_options = list(
        formula = "trait_data ~ 1",
        min_descendant_tips = 3,
        shift_acceptance_threshold = 5,
        num_cores = 1,
        IC = "GIC",
        method = "LL",
        uncertaintyweights_par = FALSE
      ),
      weighted = FALSE,
      num_cores = 1,
      seed = 123
    )
  )

  testthat::expect_identical(
    study$search_options$formula,
    simulationStudyFormulaCases$string
  )
  testthat::expect_identical(study$per_replicate$status, "ok")
  testthat::expect_equal(study$per_replicate$n_true_shifts, 1L)
  testthat::expect_equal(colnames(study$simdata[[1]]$trait_data), c("y1", "y2"))
  testthat::expect_false("mass" %in% colnames(study$simdata[[1]]$trait_data))
})

test_that("runShiftRecoverySimulationStudy works in parallel for formula-calibrated templates", {
  testthat::skip_on_cran()
  skip_if_recovery_study_deps()
  testthat::skip_on_covr()

  tmpl <- make_formula_recovery_template()
  study <- suppressWarnings(
    runShiftRecoverySimulationStudy(
      tmpl,
      n_replicates = 2,
      tree_tip_count = 20,
      simulation_options = list(
        num_shifts = 1,
        min_shift_tips = 3,
        max_shift_tips = 7,
        scale_mode = "proportional"
      ),
      search_options = list(
        formula = "trait_data ~ 1",
        min_descendant_tips = 3,
        shift_acceptance_threshold = 5,
        num_cores = 1,
        IC = "GIC",
        method = "LL",
        uncertaintyweights_par = FALSE
      ),
      weighted = FALSE,
      num_cores = 2,
      seed = 124
    )
  )

  testthat::expect_true(all(study$per_replicate$status == "ok"))
  testthat::expect_equal(colnames(study$simdata[[1]]$trait_data), c("y1", "y2"))
  testthat::expect_true(all(is.na(study$per_replicate$error)))
})

test_that("runShiftRecoverySimulationStudy rejects non-intercept search formulas", {
  testthat::skip_on_cran()
  skip_if_recovery_study_deps()

  tmpl <- make_formula_recovery_template()

  testthat::expect_error(
    runShiftRecoverySimulationStudy(
      tmpl,
      n_replicates = 1,
      simulation_options = list(
        num_shifts = 2,
        min_shift_tips = 3,
        max_shift_tips = 7,
        scale_mode = "correlation"
      ),
      search_options = list(
        formula = "trait_data[, 1:2] ~ trait_data[, 3]",
        min_descendant_tips = 3,
        shift_acceptance_threshold = 5,
        num_cores = 1,
        IC = "GIC",
        method = "LL",
        uncertaintyweights_par = FALSE
      ),
      weighted = FALSE,
      num_cores = 1,
      seed = 8
    ),
    simulationStudyFormulaErrors$non_intercept
  )
})

test_that("runShiftRecoverySimulationStudy validates inputs and inherits template error settings", {
  testthat::skip_on_cran()
  skip_if_recovery_study_deps()

  tmpl <- make_recovery_template_with_error()

  testthat::expect_error(
    runShiftRecoverySimulationStudy(list(), n_replicates = 1, simulation_options = list(num_shifts = 1, min_shift_tips = 2, max_shift_tips = 3)),
    "bifrost_simulation_template"
  )
  testthat::expect_error(
    runShiftRecoverySimulationStudy(tmpl, n_replicates = 0, simulation_options = list(num_shifts = 1, min_shift_tips = 2, max_shift_tips = 3)),
    "n_replicates"
  )
  testthat::expect_error(
    runShiftRecoverySimulationStudy(tmpl, n_replicates = 1.5, simulation_options = list(num_shifts = 1, min_shift_tips = 2, max_shift_tips = 3)),
    "n_replicates"
  )
  testthat::expect_error(
    runShiftRecoverySimulationStudy(tmpl, n_replicates = 1, tree_tip_count = 3.5, simulation_options = list(num_shifts = 1, min_shift_tips = 2, max_shift_tips = 3)),
    "tree_tip_count"
  )
  testthat::expect_error(
    runShiftRecoverySimulationStudy(tmpl, n_replicates = 1, tree_tip_count = tmpl$n_tips + 1L, simulation_options = list(num_shifts = 1, min_shift_tips = 2, max_shift_tips = 3)),
    "tree_tip_count cannot exceed"
  )
  testthat::expect_error(
    runShiftRecoverySimulationStudy(tmpl, n_replicates = 1, simulation_options = list(num_shifts = 1, min_shift_tips = 2, max_shift_tips = 3), fuzzy_distance = -1),
    "fuzzy_distance"
  )
  testthat::expect_error(
    runShiftRecoverySimulationStudy(tmpl, n_replicates = 1, simulation_options = list(num_shifts = 1, min_shift_tips = 2, max_shift_tips = 3), fuzzy_distance = 1.5),
    "fuzzy_distance"
  )
  testthat::expect_error(
    runShiftRecoverySimulationStudy(tmpl, n_replicates = 1, simulation_options = list(num_shifts = 1, min_shift_tips = 2, max_shift_tips = 3), weighted = NA),
    "weighted"
  )
  testthat::expect_error(
    runShiftRecoverySimulationStudy(tmpl, n_replicates = 1, simulation_options = list(num_shifts = 1, min_shift_tips = 2, max_shift_tips = 3), num_cores = 0),
    "num_cores"
  )
  testthat::expect_error(
    runShiftRecoverySimulationStudy(tmpl, n_replicates = 1, simulation_options = list(num_shifts = 1, min_shift_tips = 2, max_shift_tips = 3), num_cores = 1.5),
    "num_cores"
  )
  testthat::expect_error(
    runShiftRecoverySimulationStudy(tmpl, n_replicates = 1, simulation_options = list(num_shifts = 1, min_shift_tips = 2, max_shift_tips = 3), seed = 1.5),
    "seed"
  )
  testthat::expect_error(
    runShiftRecoverySimulationStudy(tmpl, n_replicates = 1, simulation_options = 1),
    "must both be lists"
  )
  testthat::expect_error(
    runShiftRecoverySimulationStudy(tmpl, n_replicates = 1, simulation_options = list()),
    "must include num_shifts, min_shift_tips, and max_shift_tips"
  )

  study <- runShiftRecoverySimulationStudy(
    tmpl,
    n_replicates = 1,
    simulation_options = list(
      num_shifts = 1,
      min_shift_tips = 2,
      max_shift_tips = 5,
      scale_mode = "proportional"
    ),
    search_options = list(
      formula = "trait_data ~ 1",
      min_descendant_tips = 2,
      shift_acceptance_threshold = 5,
      num_cores = 1,
      IC = "GIC",
      method = "LL"
    ),
    weighted = FALSE,
    num_cores = 1,
    seed = 43
  )

  testthat::expect_true(isTRUE(study$search_options$error))
})

test_that("runShiftRecoverySimulationStudy restores the caller RNG state when seeded", {
  testthat::skip_on_cran()
  skip_if_recovery_study_deps()

  tmpl <- make_recovery_template()
  local_rebind(
    "simulateShiftedDataset",
    function(template, tree_tip_count = NULL, simulation_options = NULL, ...) {
      stats::rnorm(1)
      list(
        paintedTree = template$baseline_tree,
        trait_data = template$trait_data,
        simulatedData = template$trait_data,
        shiftNodes = ape::Ntip(template$baseline_tree) + 1L,
        generating_scenario = "proportional"
      )
    },
    environment(runShiftRecoverySimulationStudy)
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
    environment(runShiftRecoverySimulationStudy)
  )

  expect_global_rng_restored(
    runShiftRecoverySimulationStudy(
      tmpl,
      n_replicates = 2,
      tree_tip_count = 20,
      simulation_options = list(
        num_shifts = 1,
        min_shift_tips = 2,
        max_shift_tips = 5,
        scale_mode = "proportional"
      ),
      search_options = list(
        formula = "trait_data ~ 1",
        min_descendant_tips = 2,
        shift_acceptance_threshold = 5,
        num_cores = 1,
        IC = "GIC",
        method = "LL"
      ),
      weighted = FALSE,
      num_cores = 1,
      seed = 43
    )
  )
})

test_that("runShiftRecoverySimulationStudy inherits evaluated template settings", {
  testthat::skip_on_cran()
  skip_if_recovery_study_deps()

  set.seed(47)
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

  study <- runShiftRecoverySimulationStudy(
    tmpl,
    n_replicates = 1,
    simulation_options = list(
      num_shifts = 1,
      min_shift_tips = 2,
      max_shift_tips = 5,
      scale_mode = "proportional"
    ),
    search_options = list(
      formula = "trait_data ~ 1",
      min_descendant_tips = 2,
      shift_acceptance_threshold = 5,
      num_cores = 1,
      IC = "GIC"
    ),
    weighted = FALSE,
    num_cores = 1,
    seed = 48
  )

  testthat::expect_identical(study$search_options$method, "LL")
  testthat::expect_true(isTRUE(study$search_options$error))
  testthat::expect_identical(study$per_replicate$status, "ok")
})

test_that("runShiftRecoverySimulationStudy falls back to template call settings when stored values are missing", {
  testthat::skip_on_cran()
  skip_if_recovery_study_deps()

  tmpl <- make_recovery_template_with_error()
  tmpl$fit_method <- NULL
  tmpl$fit_error <- NULL
  tmpl$global_model$call$method <- "LL"
  tmpl$global_model$call$error <- quote(TRUE)
  tmpl$search_formula <- NULL

  study <- runShiftRecoverySimulationStudy(
    tmpl,
    n_replicates = 1,
    simulation_options = list(
      num_shifts = 1,
      min_shift_tips = 2,
      max_shift_tips = 5,
      scale_mode = "proportional"
    ),
    search_options = list(
      min_descendant_tips = 2,
      shift_acceptance_threshold = 5,
      num_cores = 1,
      IC = "GIC"
    ),
    weighted = FALSE,
    num_cores = 1,
    seed = 49
  )

  testthat::expect_identical(study$search_options$formula, "trait_data ~ 1")
  testthat::expect_identical(study$search_options$method, "LL")
  testthat::expect_true(isTRUE(study$search_options$error))
})

test_that("runShiftRecoverySimulationStudy tolerates unevaluable template error calls", {
  testthat::skip_on_cran()
  skip_if_recovery_study_deps()

  tmpl <- make_recovery_template()
  tmpl$fit_method <- NULL
  tmpl$fit_error <- NULL
  tmpl$global_model$call$method <- "LL"
  tmpl$global_model$call$error <- quote(missing_error_flag)
  tmpl$search_formula <- NULL

  study <- runShiftRecoverySimulationStudy(
    tmpl,
    n_replicates = 1,
    simulation_options = list(
      num_shifts = 1,
      min_shift_tips = 2,
      max_shift_tips = 5,
      scale_mode = "proportional"
    ),
    search_options = list(
      min_descendant_tips = 2,
      shift_acceptance_threshold = 5,
      num_cores = 1,
      IC = "GIC"
    ),
    weighted = FALSE,
    num_cores = 1,
    seed = 50
  )

  testthat::expect_false("error" %in% names(study$search_options))
  testthat::expect_identical(study$per_replicate$status, "ok")
})

test_that("runShiftRecoverySimulationStudy is reproducible with the same seed and workers", {
  testthat::skip_on_cran()
  skip_if_recovery_study_deps()
  testthat::skip_on_covr()

  tmpl <- make_recovery_template()
  run_once <- function() {
    suppressWarnings(
      runShiftRecoverySimulationStudy(
        tmpl,
        n_replicates = 2,
        tree_tip_count = 20,
        simulation_options = list(
          num_shifts = 2,
          min_shift_tips = 3,
          max_shift_tips = 7,
          scale_mode = "proportional"
        ),
        search_options = list(
          formula = "trait_data ~ 1",
          min_descendant_tips = 3,
          shift_acceptance_threshold = 5,
          num_cores = 1,
          IC = "GIC",
          method = "LL",
          uncertaintyweights_par = FALSE
        ),
        weighted = FALSE,
        num_cores = 2,
        seed = 7
      )
    )
  }

  study_a <- run_once()
  study_b <- run_once()

  testthat::expect_identical(study_a$per_replicate, study_b$per_replicate)
  testthat::expect_true(all(vapply(seq_along(study_a$simdata), function(i) {
    isTRUE(all.equal(study_a$simdata[[i]]$simulatedData, study_b$simdata[[i]]$simulatedData))
  }, logical(1))))
  testthat::expect_identical(
    lapply(study_a$simdata, `[[`, "shiftNodes"),
    lapply(study_b$simdata, `[[`, "shiftNodes")
  )
  testthat::expect_identical(
    lapply(study_a$results, `[[`, "shift_nodes_no_uncertainty"),
    lapply(study_b$results, `[[`, "shift_nodes_no_uncertainty")
  )
})

test_that("runShiftRecoverySimulationStudy gives seeded draws independent of worker count", {
  testthat::skip_on_cran()
  skip_if_recovery_study_deps()
  testthat::skip_on_covr()

  tmpl <- make_recovery_template()
  local_rebind(
    "simulateShiftedDataset",
    function(template, tree_tip_count = NULL, ...) {
      list(
        paintedTree = template$baseline_tree,
        trait_data = template$trait_data,
        simulatedData = template$trait_data,
        shiftNodes = ape::Ntip(template$baseline_tree) + 1L,
        generating_scenario = "proportional",
        random_draw = stats::rnorm(1)
      )
    },
    environment(runShiftRecoverySimulationStudy)
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
    environment(runShiftRecoverySimulationStudy)
  )
  withr::local_envvar(c(
    RSTUDIO = "1",
    RSTUDIO_SESSION_INITIALIZED = "1"
  ))

  run_once <- function(num_cores) {
    runShiftRecoverySimulationStudy(
      tmpl,
      n_replicates = 2,
      tree_tip_count = 20,
      simulation_options = list(
        num_shifts = 1,
        min_shift_tips = 2,
        max_shift_tips = 5,
        scale_mode = "proportional"
      ),
      search_options = list(
        formula = "trait_data ~ 1",
        min_descendant_tips = 2,
        shift_acceptance_threshold = 5,
        num_cores = 1,
        IC = "GIC",
        method = "LL"
      ),
      weighted = FALSE,
      num_cores = num_cores,
      seed = 707
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

test_that("runShiftRecoverySimulationStudy makes outer parallelism override inner cores", {
  testthat::skip_on_cran()
  skip_if_recovery_study_deps()

  tmpl <- make_recovery_template()

  testthat::local_mocked_bindings(
    future_lapply = serial_future_lapply,
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
        ic_weights = data.frame(
          node = integer(0),
          ic_with_shift = numeric(0),
          ic_without_shift = numeric(0),
          delta_ic = numeric(0),
          ic_weight_withshift = numeric(0),
          ic_weight_withoutshift = numeric(0),
          evidence_ratio = numeric(0)
        )
      )
    },
    environment(runShiftRecoverySimulationStudy)
  )

  testthat::expect_warning(
    study <- runShiftRecoverySimulationStudy(
      tmpl,
      n_replicates = 1,
      simulation_options = list(
        num_shifts = 1,
        min_shift_tips = 2,
        max_shift_tips = 5,
        scale_mode = "proportional"
      ),
      search_options = list(
        formula = "trait_data ~ 1",
        min_descendant_tips = 2,
        shift_acceptance_threshold = 5,
        num_cores = 2,
        IC = "GIC",
        method = "LL"
      ),
      weighted = FALSE,
      num_cores = 2,
      seed = 44
    ),
    "forcing each search to use one core"
  )
  testthat::expect_identical(study$search_options$num_cores, 1L)
})

test_that("runShiftRecoverySimulationStudy does not capture simdata in future workers", {
  testthat::skip_on_cran()
  skip_if_recovery_study_deps()

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
    environment(runShiftRecoverySimulationStudy)
  )

  runShiftRecoverySimulationStudy(
    make_recovery_template(),
    n_replicates = 1,
    tree_tip_count = 20,
    simulation_options = list(
      num_shifts = 1,
      min_shift_tips = 3,
      max_shift_tips = 7,
      scale_mode = "proportional"
    ),
    search_options = list(
      formula = "trait_data ~ 1",
      min_descendant_tips = 3,
      shift_acceptance_threshold = 5,
      num_cores = 1,
      IC = "GIC",
      method = "LL",
      uncertaintyweights_par = FALSE
    ),
    weighted = FALSE,
    num_cores = 2,
    seed = 55
  )

  testthat::expect_length(worker_environments, 2L)
  testthat::expect_false(any(vapply(worker_environments, function(worker_env) {
    exists("simdata", envir = worker_env, inherits = FALSE)
  }, logical(1L))))
})

test_that("runShiftRecoverySimulationStudy validates formula types and can choose multisession plans", {
  testthat::skip_on_cran()
  skip_if_recovery_study_deps()

  tmpl <- make_recovery_template()

  testthat::expect_error(
    runShiftRecoverySimulationStudy(
      tmpl,
      n_replicates = 1,
      simulation_options = list(
        num_shifts = 1,
        min_shift_tips = 2,
        max_shift_tips = 5,
        scale_mode = "proportional"
      ),
      search_options = list(
        formula = 1,
        min_descendant_tips = 2,
        shift_acceptance_threshold = 5,
        num_cores = 1,
        IC = "GIC"
      ),
      weighted = FALSE,
      num_cores = 1,
      seed = 51
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
  local_rebind(
    "searchOptimalConfiguration",
    function(...) {
      list(
        shift_nodes_no_uncertainty = integer(0),
        num_candidates = 3L,
        ic_weights = data.frame(
          node = integer(0),
          ic_with_shift = numeric(0),
          ic_without_shift = numeric(0),
          delta_ic = numeric(0),
          ic_weight_withshift = numeric(0),
          ic_weight_withoutshift = numeric(0),
          evidence_ratio = numeric(0)
        )
      )
    },
    environment(runShiftRecoverySimulationStudy)
  )

  runShiftRecoverySimulationStudy(
    tmpl,
    n_replicates = 1,
    simulation_options = list(
      num_shifts = 1,
      min_shift_tips = 2,
      max_shift_tips = 5,
      scale_mode = "proportional"
    ),
    search_options = list(
      formula = "trait_data ~ 1",
      min_descendant_tips = 2,
      shift_acceptance_threshold = 5,
      num_cores = 1,
      IC = "GIC",
      method = "LL"
    ),
    weighted = FALSE,
    num_cores = 2,
    seed = 52
  )

  testthat::expect_true(any(vapply(plans, function(x) identical(x$strategy, future::multisession), logical(1))))
})

test_that("runShiftRecoverySimulationStudy accepts formula objects in search options", {
  testthat::skip_on_cran()
  skip_if_recovery_study_deps()

  tmpl <- make_recovery_template()
  study <- runShiftRecoverySimulationStudy(
    tmpl,
    n_replicates = 1,
    simulation_options = list(
      num_shifts = 1,
      min_shift_tips = 2,
      max_shift_tips = 5,
      scale_mode = "proportional"
    ),
    search_options = list(
      formula = simulationStudyFormulaCases$formula,
      min_descendant_tips = 2,
      shift_acceptance_threshold = 5,
      num_cores = 1,
      IC = "GIC",
      method = "LL"
    ),
    weighted = FALSE,
    num_cores = 1,
    seed = 53
  )

  testthat::expect_s3_class(study$search_options$formula, "formula")
  testthat::expect_identical(
    study$search_options$formula,
    simulationStudyFormulaCases$formula
  )
  testthat::expect_identical(study$per_replicate$status, "ok")
})

test_that("runShiftRecoverySimulationStudy records search fallback errors", {
  testthat::skip_on_cran()
  skip_if_recovery_study_deps()

  tmpl <- make_recovery_template()

  study <- runShiftRecoverySimulationStudy(
    tmpl,
    n_replicates = 1,
    simulation_options = list(
      num_shifts = 1,
      min_shift_tips = 2,
      max_shift_tips = 5,
      scale_mode = "proportional"
    ),
    search_options = list(
      formula = "trait_data ~ 1",
      min_descendant_tips = 2,
      shift_acceptance_threshold = 5,
      num_cores = 1,
      IC = "BAD",
      method = "LL"
    ),
    weighted = FALSE,
    num_cores = 1,
    seed = 45
  )

  testthat::expect_identical(study$per_replicate$status, "error")
  testthat::expect_match(study$per_replicate$error, "IC must be GIC or BIC")
  testthat::expect_true(is.na(study$per_replicate$n_inferred_shifts))
  testthat::expect_identical(study$study_summary$n_completed, 0L)
  testthat::expect_identical(study$study_summary$n_failed, 1L)
  testthat::expect_identical(study$study_summary$n_evaluable_replicates, 0L)
  testthat::expect_equal(study$study_summary$completion_rate, 0)
  testthat::expect_equal(study$study_summary$failure_rate, 1)
  testthat::expect_equal(unname(study$evaluation$counts$strict), rep(0, 4))
  testthat::expect_equal(unname(study$evaluation$counts$fuzzy), rep(0, 4))
  testthat::expect_true(all(is.na(unlist(study$evaluation$strict))))
  testthat::expect_true(all(is.na(unlist(study$evaluation$fuzzy))))
  testthat::expect_equal(
    study$results[[1]]$num_candidates,
    max(length(generatePaintedTrees(ape::as.phylo(study$simdata[[1]]$paintedTree), min_tips = 2)) - 1L, 0L)
  )
})

test_that("runShiftRecoverySimulationStudy marks zero-candidate recovery metrics as unevaluable", {
  testthat::skip_on_cran()
  skip_if_recovery_study_deps()

  tmpl <- make_recovery_template()

  study <- runShiftRecoverySimulationStudy(
    tmpl,
    n_replicates = 1,
    simulation_options = list(
      num_shifts = 1,
      min_shift_tips = 2,
      max_shift_tips = 5,
      scale_mode = "proportional"
    ),
    search_options = list(
      formula = "trait_data ~ 1",
      min_descendant_tips = 100,
      shift_acceptance_threshold = 5,
      num_cores = 1,
      IC = "GIC",
      method = "LL"
    ),
    weighted = FALSE,
    num_cores = 1,
    seed = 49
  )

  testthat::expect_equal(study$per_replicate$n_candidates, 0)
  testthat::expect_true(is.na(study$evaluation$strict$specificity))
  testthat::expect_true(is.na(study$evaluation$strict$fpr))
  testthat::expect_true(is.na(study$evaluation$strict$balanced_accuracy))
  testthat::expect_true(is.na(study$evaluation$fuzzy$specificity))
  testthat::expect_true(is.na(study$evaluation$fuzzy$fpr))
  testthat::expect_true(is.na(study$evaluation$fuzzy$balanced_accuracy))
})

test_that("runShiftRecoverySimulationStudy handles missing trait_data and mixed scenarios", {
  testthat::skip_on_cran()
  skip_if_recovery_study_deps()

  tmpl <- make_recovery_template()
  tr <- phytools::paintSubTree(ape::rtree(6), node = 7, state = "ancestral", anc.state = "ancestral")
  sim_stub <- list(
    list(
      paintedTree = tr,
      shiftNodes = 7L,
      simulatedData = matrix(rnorm(12), ncol = 2, dimnames = list(tr$tip.label, c("y1", "y2"))),
      trait_data = NULL,
      generating_scenario = "proportional"
    ),
    list(
      paintedTree = tr,
      shiftNodes = 7L,
      simulatedData = matrix(rnorm(12), ncol = 2, dimnames = list(tr$tip.label, c("y1", "y2"))),
      trait_data = NULL,
      generating_scenario = "correlation"
    )
  )

  mock_env <- new.env(parent = emptyenv())
  mock_env$sim_i <- 0L
  mock_env$search_i <- 0L
  mock_env$seen_trait_data <- vector("list", length(sim_stub))
  local_rebind(
    "simulateShiftedDataset",
    function(...) {
      mock_env$sim_i <- mock_env$sim_i + 1L
      sim_stub[[mock_env$sim_i]]
    },
    environment(runShiftRecoverySimulationStudy)
  )
  local_rebind(
    "searchOptimalConfiguration",
    function(baseline_tree, trait_data, ...) {
      mock_env$search_i <- mock_env$search_i + 1L
      search_i <- mock_env$search_i
      mock_env$seen_trait_data[[search_i]] <- trait_data
      list(
        shift_nodes_no_uncertainty = integer(0),
        num_candidates = 3L,
        ic_weights = data.frame(
          node = integer(0),
          ic_with_shift = numeric(0),
          ic_without_shift = numeric(0),
          delta_ic = numeric(0),
          ic_weight_withshift = numeric(0),
          ic_weight_withoutshift = numeric(0),
          evidence_ratio = numeric(0)
        )
      )
    },
    environment(runShiftRecoverySimulationStudy)
  )
  local_rebind(
    "evaluateShiftRecovery",
    function(...) {
      list(
        strict = list(precision = 0, recall = 0),
        fuzzy = list(precision = 0, recall = 0),
        weighted = NULL,
        counts = list(strict = c(TP = 0, FP = 0, FN = 0, TN = 0), fuzzy = c(TP = 0, FP = 0, FN = 0, TN = 0))
      )
    },
    environment(runShiftRecoverySimulationStudy)
  )
  study <- runShiftRecoverySimulationStudy(
    tmpl,
    n_replicates = 2,
    simulation_options = list(
      num_shifts = 1,
      min_shift_tips = 2,
      max_shift_tips = 5,
      scale_mode = "proportional"
    ),
    search_options = list(
      formula = "trait_data ~ 1",
      min_descendant_tips = 2,
      shift_acceptance_threshold = 5,
      num_cores = 1,
      IC = "GIC",
      method = "LL"
    ),
    weighted = FALSE,
    num_cores = 1,
    seed = 46
  )

  testthat::expect_identical(study$generating_scenario, c("proportional", "correlation"))
  testthat::expect_equal(mock_env$seen_trait_data, lapply(sim_stub, `[[`, "simulatedData"))
})

test_that("runShiftRecoverySimulationStudy integrates weighted metrics across mixed evaluable replicates", {
  testthat::skip_on_cran()
  skip_if_recovery_study_deps()

  tmpl <- make_recovery_template()
  tr <- phytools::paintSubTree(ape::rtree(6), node = 7, state = "ancestral", anc.state = "ancestral")
  sim_stub <- list(
    list(
      paintedTree = tr,
      shiftNodes = 7L,
      simulatedData = matrix(rnorm(12), ncol = 2, dimnames = list(tr$tip.label, c("y1", "y2"))),
      trait_data = NULL,
      generating_scenario = "proportional"
    ),
    list(
      paintedTree = tr,
      shiftNodes = 7L,
      simulatedData = matrix(rnorm(12), ncol = 2, dimnames = list(tr$tip.label, c("y1", "y2"))),
      trait_data = NULL,
      generating_scenario = "proportional"
    )
  )
  result_stub <- list(
    list(
      shift_nodes_no_uncertainty = integer(0),
      num_candidates = 0L,
      ic_weights = data.frame(
        node = integer(0),
        ic_with_shift = numeric(0),
        ic_without_shift = numeric(0),
        delta_ic = numeric(0),
        ic_weight_withshift = numeric(0),
        ic_weight_withoutshift = numeric(0),
        evidence_ratio = numeric(0)
      )
    ),
    list(
      shift_nodes_no_uncertainty = 7L,
      num_candidates = 4L,
      ic_weights = data.frame(
        node = 7L,
        ic_with_shift = 1,
        ic_without_shift = 2,
        delta_ic = 1,
        ic_weight_withshift = 0.8,
        ic_weight_withoutshift = 0.2,
        evidence_ratio = 4
      )
    )
  )

  mock_env <- new.env(parent = emptyenv())
  mock_env$sim_i <- 0L
  mock_env$search_i <- 0L
  mock_env$seen_trait_data <- vector("list", length(sim_stub))
  local_rebind(
    "simulateShiftedDataset",
    function(...) {
      mock_env$sim_i <- mock_env$sim_i + 1L
      sim_stub[[mock_env$sim_i]]
    },
    environment(runShiftRecoverySimulationStudy)
  )
  local_rebind(
    "searchOptimalConfiguration",
    function(baseline_tree, trait_data, ...) {
      mock_env$search_i <- mock_env$search_i + 1L
      search_i <- mock_env$search_i
      mock_env$seen_trait_data[[search_i]] <- trait_data
      result_stub[[search_i]]
    },
    environment(runShiftRecoverySimulationStudy)
  )

  study <- runShiftRecoverySimulationStudy(
    tmpl,
    n_replicates = 2,
    simulation_options = list(
      num_shifts = 1,
      min_shift_tips = 2,
      max_shift_tips = 5,
      scale_mode = "proportional"
    ),
    search_options = list(
      formula = "trait_data ~ 1",
      min_descendant_tips = 2,
      shift_acceptance_threshold = 5,
      num_cores = 1,
      IC = "GIC",
      method = "LL"
    ),
    weighted = TRUE,
    num_cores = 1,
    seed = 50
  )

  testthat::expect_equal(study$per_replicate$n_candidates, c(0, 4))
  testthat::expect_equal(mock_env$seen_trait_data, lapply(sim_stub, `[[`, "simulatedData"))
  testthat::expect_false(is.null(study$evaluation$weighted))
  testthat::expect_equal(study$evaluation$strict$precision, 1)
  testthat::expect_equal(study$evaluation$strict$recall, 0.5)
  testthat::expect_equal(study$evaluation$strict$specificity, 1)
  testthat::expect_equal(study$evaluation$weighted$strict$precision, 1)
  testthat::expect_equal(study$evaluation$weighted$strict$recall, 0.4)
})
