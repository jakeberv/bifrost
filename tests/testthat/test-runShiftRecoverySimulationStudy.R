testthat::skip_on_cran()

skip_if_recovery_study_deps <- function() {
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

make_recovery_template <- function() {
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
}

make_formula_recovery_template <- function() {
  set.seed(41)
  tr <- ape::rtree(24)
  X <- cbind(y1 = rnorm(24), y2 = rnorm(24), mass = rnorm(24))
  rownames(X) <- tr$tip.label
  createSimulationTemplate(
    baseline_tree = tr,
    trait_data = X,
    formula = "trait_data[, 1:2] ~ trait_data[, 3]",
    response_columns = 1:2,
    predictor_columns = 3,
    method = "LL"
  )
}

make_recovery_template_with_error <- function() {
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
}

test_that("runShiftRecoverySimulationStudy supports the proportional scenario", {
  skip_if_recovery_study_deps()

  tmpl <- make_recovery_template()
  study <- suppressWarnings(
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
      num_cores = 1,
      seed = 7
    )
  )

  testthat::expect_s3_class(study, "bifrost_simulation_study")
  testthat::expect_identical(study$study_type, "shift_recovery")
  testthat::expect_identical(study$generating_scenario, "proportional")
  testthat::expect_length(study$simdata, 2)
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
})

test_that("runShiftRecoverySimulationStudy supports the correlation scenario", {
  skip_if_recovery_study_deps()

  tmpl <- make_recovery_template()
  study <- runShiftRecoverySimulationStudy(
    tmpl,
    n_replicates = 1,
    tree_tip_count = 20,
    simulation_options = list(
      num_shifts = 2,
      min_shift_tips = 3,
      max_shift_tips = 7,
      scale_mode = "correlation"
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

  testthat::expect_identical(study$generating_scenario, "correlation")
  testthat::expect_equal(study$per_replicate$generating_scenario, "correlation")
})

test_that("runShiftRecoverySimulationStudy rejects per-replicate seeds", {
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

  testthat::expect_identical(study$search_options$formula, "trait_data ~ 1")
  testthat::expect_identical(study$per_replicate$status, "ok")
  testthat::expect_equal(study$per_replicate$n_true_shifts, 1L)
  testthat::expect_equal(colnames(study$simdata[[1]]$trait_data), c("y1", "y2"))
  testthat::expect_false("mass" %in% colnames(study$simdata[[1]]$trait_data))
})

test_that("runShiftRecoverySimulationStudy rejects non-intercept search formulas", {
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
    "intercept-only search formulas only"
  )
})

test_that("runShiftRecoverySimulationStudy validates inputs and inherits template error settings", {
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
    runShiftRecoverySimulationStudy(tmpl, n_replicates = 1, simulation_options = list(num_shifts = 1, min_shift_tips = 2, max_shift_tips = 3), fuzzy_distance = -1),
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

test_that("runShiftRecoverySimulationStudy inherits evaluated template settings", {
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

test_that("runShiftRecoverySimulationStudy is reproducible with the same seed and workers", {
  skip_if_recovery_study_deps()

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

test_that("runShiftRecoverySimulationStudy warns on nested parallelism", {
  skip_if_recovery_study_deps()

  tmpl <- make_recovery_template()

  testthat::local_mocked_bindings(
    future_lapply = function(X, FUN, ...) lapply(X, FUN),
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
        num_cores = 2,
        IC = "GIC",
        method = "LL"
      ),
      weighted = FALSE,
      num_cores = 2,
      seed = 44
    ),
    "nested parallelism"
  )
})

test_that("runShiftRecoverySimulationStudy records search fallback errors", {
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
  testthat::expect_equal(
    study$results[[1]]$num_candidates,
    max(length(generatePaintedTrees(ape::as.phylo(study$simdata[[1]]$paintedTree), min_tips = 2)) - 1L, 0L)
  )
})

test_that("runShiftRecoverySimulationStudy marks zero-candidate recovery metrics as unevaluable", {
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
  mock_env$mock_i <- 0L
  local_rebind(
    "simulateShiftedDataset",
    function(...) {
      sim_stub[[get("mock_i", envir = mock_env)]]
    },
    environment(runShiftRecoverySimulationStudy)
  )
  local_rebind(
    "searchOptimalConfiguration",
    function(baseline_tree, trait_data, ...) {
      testthat::expect_equal(trait_data, sim_stub[[get("mock_i", envir = mock_env)]]$simulatedData)
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
  testthat::local_mocked_bindings(
    future_lapply = function(X, FUN, ...) {
      lapply(X, function(i) {
        mock_env$mock_i <- i
        FUN(i)
      })
    },
    .package = "future.apply"
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

  testthat::expect_length(study$generating_scenario, 2)
  testthat::expect_true(all(c("proportional", "correlation") %in% study$generating_scenario))
})

test_that("runShiftRecoverySimulationStudy integrates weighted metrics across mixed evaluable replicates", {
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
  mock_env$mock_i <- 0L
  local_rebind(
    "simulateShiftedDataset",
    function(...) {
      sim_stub[[get("mock_i", envir = mock_env)]]
    },
    environment(runShiftRecoverySimulationStudy)
  )
  local_rebind(
    "searchOptimalConfiguration",
    function(baseline_tree, trait_data, ...) {
      testthat::expect_equal(trait_data, sim_stub[[get("mock_i", envir = mock_env)]]$simulatedData)
      result_stub[[get("mock_i", envir = mock_env)]]
    },
    environment(runShiftRecoverySimulationStudy)
  )
  testthat::local_mocked_bindings(
    future_lapply = function(X, FUN, ...) {
      lapply(X, function(i) {
        mock_env$mock_i <- i
        FUN(i)
      })
    },
    .package = "future.apply"
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
  testthat::expect_false(is.null(study$evaluation$weighted))
  testthat::expect_equal(study$evaluation$strict$precision, 1)
  testthat::expect_equal(study$evaluation$strict$recall, 0.5)
  testthat::expect_equal(study$evaluation$strict$specificity, 1)
  testthat::expect_equal(study$evaluation$weighted$strict$precision, 1)
  testthat::expect_equal(study$evaluation$weighted$strict$recall, 0.4)
})
