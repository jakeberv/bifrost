.simulation_study_compact_call <- function(function_name, arguments) {
  arguments$template <- "<bifrost_simulation_template>"
  c(list(as.name(function_name)), arguments)
}

.simulation_study_simulate_null <- function(replicate,
                                             template,
                                             tree_tip_count,
                                             simulation_options,
                                             simulation_seeds,
                                             simulate_dataset,
                                             with_replicate_seed) {
  force(replicate)
  sim_args <- utils::modifyList(
    simulation_options,
    list(
      template = template,
      tree_tip_count = tree_tip_count
    )
  )
  result <- with_replicate_seed(
    simulation_seeds[[replicate]],
    do.call(simulate_dataset, sim_args)
  )
  result$user_input <- .simulation_study_compact_call("simulateNullDataset", sim_args)
  result
}

.simulation_study_simulate_shifted <- function(replicate,
                                                template,
                                                tree_tip_count,
                                                simulation_options,
                                                simulation_seeds,
                                                simulate_dataset,
                                                with_replicate_seed,
                                                compact_call,
                                                progress = NULL) {
  force(replicate)
  on.exit({
    if (!is.null(progress)) progress()
  }, add = TRUE)
  sim_args <- utils::modifyList(
    simulation_options,
    list(
      template = template,
      tree_tip_count = tree_tip_count
    )
  )
  result <- with_replicate_seed(
    simulation_seeds[[replicate]],
    do.call(simulate_dataset, sim_args)
  )
  result$user_input <- compact_call("simulateShiftedDataset", sim_args)
  result
}

.simulation_study_search <- function(replicate,
                                     search_options,
                                     tree_component,
                                     data_components,
                                     search_optimal_configuration,
                                     generate_painted_trees,
                                     with_replicate_seed,
                                     progress = NULL) {
  on.exit({
    if (!is.null(progress)) progress()
  }, add = TRUE)
  sim <- replicate$sim

  trait_data <- NULL
  for (component in data_components) {
    if (!is.null(sim[[component]])) {
      trait_data <- sim[[component]]
      break
    }
  }
  search_args <- utils::modifyList(
    search_options,
    list(
      baseline_tree = ape::as.phylo(sim[[tree_component]]),
      trait_data = trait_data
    )
  )

  with_replicate_seed(
    replicate$seed,
    tryCatch(
      do.call(search_optimal_configuration, search_args),
      error = function(e) {
        candidate_count <- max(length(generate_painted_trees(
          ape::as.phylo(sim[[tree_component]]),
          min_tips = search_options$min_descendant_tips
        )) - 1L, 0L)
        list(
          shift_nodes_no_uncertainty = integer(0),
          num_candidates = candidate_count,
          ic_weights = data.frame(
            node = integer(0),
            ic_with_shift = numeric(0),
            ic_without_shift = numeric(0),
            delta_ic = numeric(0),
            ic_weight_withshift = numeric(0),
            ic_weight_withoutshift = numeric(0),
            evidence_ratio = numeric(0)
          ),
          error = conditionMessage(e)
        )
      }
    )
  )
}

#' Run a False-Positive Simulation Study
#'
#' @description
#' Repeatedly simulate null datasets from a
#' [`createSimulationTemplate()`][createSimulationTemplate] object and analyze
#' each replicate with [searchOptimalConfiguration()] to estimate the expected
#' false-positive behavior of a candidate `bifrost` search setup.
#'
#' @param template A `bifrost_simulation_template` returned by
#'   [createSimulationTemplate()].
#' @param n_replicates Integer number of simulation replicates to run.
#' @param tree_tip_count Optional integer tip count for random subtree sampling.
#'   If `NULL`, each replicate uses the full empirical tree.
#' @param simulation_options Named list of additional arguments passed to
#'   [simulateNullDataset()], including `simulation_generator` and `covariance_df`.
#'   The manuscript-compatible `"original"` generator is used by default;
#'   request `simulation_generator = "empirical"` for full-covariance draws.
#'   The `simulation_generator` entry must contain exactly one supported string.
#'   For replicated studies, `simulation_options$seed` is not allowed; use the
#'   wrapper-level `seed` argument instead.
#' @param search_options Named list of arguments passed to
#'   [searchOptimalConfiguration()]. These override the default response-only
#'   intercept-only search settings assembled from `template`. If `formula` is
#'   omitted, the wrapper inherits `template$search_formula`, which is
#'   intercept-only in the dataset-matched simulation workflow.
#' @param num_cores Integer number of workers used across replicate analyses.
#'   When greater than one, outer replicate parallelism takes precedence and
#'   each inner search is forced to `search_options$num_cores = 1`.
#' @param seed Optional integer random seed. When supplied, the wrapper restores
#'   the caller's previous RNG state before returning.
#'
#' @details
#' Simulated datasets, search results, a per-replicate summary table, and a
#' compact study summary are returned in a single object of class
#' `bifrost_simulation_study`.
#'
#' Regardless of the global calibration model used to build `template`, the
#' downstream search is intentionally restricted to intercept-only formulas.
#' This keeps the simulation study focused on branch-shift detection in the
#' response block, matching the residual-calibration logic. The
#' simulated `trait_data` supplied to each search result therefore contains only
#' the regenerated response variables.
#' Templates calibrated with covariates are supported; those covariates influence
#' the fitted means and residual covariance used to generate each response block,
#' but they are not re-fit within each simulated search replicate.
#'
#' Failed search replicates are retained in the output, but their inferred-shift
#' and false-positive metrics are recorded as `NA` and excluded from scientific
#' summaries. Their error messages and candidate counts remain available for
#' diagnosis, and completion and failure rates are reported separately.
#' Reproducibility for replicated studies is controlled by the wrapper-level
#' `seed` argument. A supplied seed deterministically creates separate
#' per-replicate simulation and search streams, so serial and parallel runs use
#' the same replicate-level random draws. Parallel execution also uses
#' `future.seed = TRUE`; per-replicate simulator seeds are intentionally
#' disallowed here.
#' If no replicate yields an evaluable false-positive rate, the study-level mean
#' and median false-positive summaries are returned as `NA`.
#'
#' @return A list of class `bifrost_simulation_study` containing the simulated
#'   datasets, raw search results, a per-replicate summary table, and a compact
#'   study-level summary for the null scenario.
#'
#' @seealso [simulateNullDataset()], [runShiftRecoverySimulationStudy()],
#'   [searchOptimalConfiguration()]
#'
#' @examples
#' \dontrun{
#' set.seed(1)
#' tr <- ape::rtree(30)
#' X <- matrix(rnorm(30 * 3), ncol = 3)
#' rownames(X) <- tr$tip.label
#' tmpl <- createSimulationTemplate(tr, X, formula = "trait_data ~ 1", method = "LL")
#'
#' fp_study <- runFalsePositiveSimulationStudy(
#'   tmpl,
#'   n_replicates = 2,
#'   tree_tip_count = 20,
#'   search_options = list(
#'     formula = "trait_data ~ 1",
#'     min_descendant_tips = 3,
#'     shift_acceptance_threshold = 5,
#'     num_cores = 1,
#'     IC = "GIC",
#'     method = "LL"
#'   ),
#'   num_cores = 1,
#'   seed = 2
#' )
#'
#' fp_study
#' }
#'
#' @export
runFalsePositiveSimulationStudy <- function(template,
                                            n_replicates,
                                            tree_tip_count = NULL,
                                            simulation_options = list(),
                                            search_options = list(),
                                            num_cores = 1,
                                            seed = NULL) {
  if (!inherits(template, "bifrost_simulation_template")) {
    stop("template must be a 'bifrost_simulation_template' object.")
  }
  n_replicates <- .simulation_check_integer_scalar(
    n_replicates,
    name = "n_replicates",
    minimum = 1L,
    message = "n_replicates must be a single integer >= 1."
  )
  tree_tip_count <- .simulation_check_integer_scalar(
    tree_tip_count,
    name = "tree_tip_count",
    minimum = 2L,
    allow_null = TRUE,
    message = "tree_tip_count must be NULL or a single integer >= 2."
  )
  if (!is.null(tree_tip_count) && tree_tip_count > template$n_tips) {
    stop("tree_tip_count cannot exceed the number of tips in the template.")
  }
  num_cores <- .simulation_check_integer_scalar(
    num_cores,
    name = "num_cores",
    minimum = 1L,
    message = "num_cores must be a single integer >= 1."
  )
  if (!is.list(simulation_options) || !is.list(search_options)) {
    stop("simulation_options and search_options must both be lists.")
  }
  if (!is.null(simulation_options$seed)) {
    stop(
      "Do not supply simulation_options$seed to runFalsePositiveSimulationStudy(); ",
      "use the wrapper-level seed argument instead."
    )
  }
  requested_generator <- .simulation_generator_from_options(simulation_options)
  seed_state <- .simulation_set_seed(seed)
  on.exit(.simulation_restore_seed(seed_state), add = TRUE)
  simulation_seeds <- .simulation_replicate_seeds(
    n_replicates,
    seeded = !is.null(seed)
  )
  search_seeds <- .simulation_replicate_seeds(
    n_replicates,
    seeded = !is.null(seed)
  )

  template_search_formula <- if (!is.null(template$search_formula)) {
    template$search_formula
  } else {
    "trait_data ~ 1"
  }

  search_defaults <- list(
    formula = template_search_formula,
    min_descendant_tips = 10,
    num_cores = 1,
    shift_acceptance_threshold = 10,
    IC = "GIC",
    plot = FALSE,
    store_model_fit_history = FALSE,
    verbose = FALSE,
    progress = FALSE
  )
  fit_settings <- .simulation_template_fit_settings(template)
  if (!is.null(fit_settings$method)) {
    search_defaults$method <- fit_settings$method
  }
  if (!is.null(fit_settings$error)) {
    search_defaults$error <- fit_settings$error
  }
  search_opts <- utils::modifyList(search_defaults, search_options)
  search_opts$formula <- validateSimulationStudyFormula(search_opts$formula)
  if (isTRUE(num_cores > 1L) && isTRUE(search_opts$num_cores > 1L)) {
    warning(
      "Both wrapper-level and search-level parallelism are > 1; ",
      "forcing each search to use one core so outer replicate parallelism ",
      "takes precedence.",
      call. = FALSE
    )
    search_opts$num_cores <- 1L
  }

  simulate_null_dataset_fn <- simulateNullDataset
  search_optimal_configuration_fn <- searchOptimalConfiguration
  generate_painted_trees_fn <- generatePaintedTrees
  with_replicate_seed_fn <- .simulation_seed_runner()

  simdata <- lapply(
    seq_len(n_replicates),
    .simulation_study_simulate_null,
    template = template,
    tree_tip_count = tree_tip_count,
    simulation_options = simulation_options,
    simulation_seeds = simulation_seeds,
    simulate_dataset = simulate_null_dataset_fn,
    with_replicate_seed = with_replicate_seed_fn
  )

  search_replicates <- Map(
    function(sim, replicate_seed) list(sim = sim, seed = replicate_seed),
    simdata,
    search_seeds
  )

  if (num_cores > 1L) {
    old_plan <- future::plan()
    on.exit(future::plan(old_plan), add = TRUE)
    use_multicore <- .Platform$OS.type != "windows" &&
      !identical(Sys.getenv("RSTUDIO"), "1") &&
      !identical(Sys.getenv("RSTUDIO_SESSION_INITIALIZED"), "1")
    if (use_multicore) {
      future::plan(future::multicore, workers = num_cores)
    } else {
      future::plan(future::multisession, workers = num_cores)
    }

    results <- progressr::with_progress({
      progress <- progressr::progressor(along = search_replicates)
      future.apply::future_lapply(
        search_replicates,
        .simulation_study_search,
        search_options = search_opts,
        tree_component = "tree",
        data_components = "data",
        search_optimal_configuration = search_optimal_configuration_fn,
        generate_painted_trees = generate_painted_trees_fn,
        with_replicate_seed = with_replicate_seed_fn,
        progress = progress,
        future.seed = TRUE
      )
    })
  } else {
    results <- lapply(
      search_replicates,
      .simulation_study_search,
      search_options = search_opts,
      tree_component = "tree",
      data_components = "data",
      search_optimal_configuration = search_optimal_configuration_fn,
      generate_painted_trees = generate_painted_trees_fn,
      with_replicate_seed = with_replicate_seed_fn
    )
  }

  replicate_generators <- vapply(simdata, function(sim) {
    if (is.null(sim$simulation_generator)) requested_generator else sim$simulation_generator
  }, character(1L))
  simulation_generators <- unique(replicate_generators)

  per_replicate <- do.call(rbind, lapply(seq_along(results), function(i) {
    n_candidates <- results[[i]]$num_candidates
    status <- if (is.null(results[[i]]$error)) "ok" else "error"
    n_inferred <- if (identical(status, "ok")) {
      length(results[[i]]$shift_nodes_no_uncertainty)
    } else {
      NA_integer_
    }
    data.frame(
      replicate = i,
      generating_scenario = "null",
      simulation_generator = replicate_generators[i],
      n_candidates = n_candidates,
      n_inferred_shifts = n_inferred,
      false_positive_rate = if (identical(status, "ok") && is.numeric(n_candidates) &&
          !is.na(n_candidates) && n_candidates > 0) {
        n_inferred / n_candidates
      } else {
        NA_real_
      },
      status = status,
      error = if (is.null(results[[i]]$error)) NA_character_ else as.character(results[[i]]$error),
      stringsAsFactors = FALSE
    )
  }))
  rownames(per_replicate) <- NULL

  evaluable_fp_rates <- per_replicate$false_positive_rate[!is.na(per_replicate$false_positive_rate)]
  study_summary <- list(
    n_replicates = n_replicates,
    n_completed = sum(per_replicate$status == "ok"),
    n_failed = sum(per_replicate$status == "error"),
    simulation_generator = simulation_generators,
    completion_rate = mean(per_replicate$status == "ok"),
    failure_rate = mean(per_replicate$status == "error"),
    n_evaluable_replicates = length(evaluable_fp_rates),
    mean_false_positive_rate = if (length(evaluable_fp_rates) > 0L) {
      mean(evaluable_fp_rates)
    } else {
      NA_real_
    },
    median_false_positive_rate = if (length(evaluable_fp_rates) > 0L) {
      stats::median(evaluable_fp_rates)
    } else {
      NA_real_
    }
  )

  out <- list(
    user_input = as.list(match.call()),
    study_type = "false_positive",
    generating_scenario = "null",
    simulation_generator = simulation_generators,
    simdata = simdata,
    results = results,
    per_replicate = per_replicate,
    study_summary = study_summary,
    simulation_options = simulation_options,
    search_options = search_opts
  )
  class(out) <- c("bifrost_simulation_study", class(out))
  out
}

#' Run a Shift-Recovery Simulation Study
#'
#' @description
#' Repeatedly simulate known-shift datasets from a
#' [`createSimulationTemplate()`][createSimulationTemplate] object and analyze
#' each replicate with [searchOptimalConfiguration()] to estimate expected shift
#' recovery performance. The default `"original"` generator reproduces the
#' published simulation operations. The `"empirical"` generator remains
#' available explicitly for full-covariance draws and an empirically calibrated
#' integration-rate robustness scenario.
#'
#' @param template A `bifrost_simulation_template` returned by
#'   [createSimulationTemplate()].
#' @param n_replicates Integer number of simulation replicates to run.
#' @param tree_tip_count Optional integer tip count for random subtree sampling.
#'   If `NULL`, each replicate uses the full empirical tree.
#' @param simulation_options Named list of additional arguments passed to
#'   [simulateShiftedDataset()], including covariance and integration-rate
#'   trade-off controls. The manuscript-compatible `"original"` generator is
#'   used by default; request `simulation_generator = "empirical"` for the
#'   full-covariance option. The `simulation_generator` entry must contain
#'   exactly one supported string. For replicated studies,
#'   `simulation_options$seed` is not allowed; use the wrapper-level `seed`
#'   argument instead.
#' @param search_options Named list of arguments passed to
#'   [searchOptimalConfiguration()]. These override the default response-only
#'   intercept-only search settings assembled from `template`. If `formula` is
#'   omitted, the wrapper inherits `template$search_formula`, which is
#'   intercept-only in the dataset-matched simulation workflow.
#' @param fuzzy_distance Integer node distance used for fuzzy matching in
#'   [evaluateShiftRecovery()].
#' @param weighted Logical; if `TRUE`, compute weighted recovery summaries using
#'   IC weights when available.
#' @param num_cores Integer number of workers used across replicate analyses.
#'   When greater than one, outer replicate parallelism takes precedence and
#'   each inner search is forced to `search_options$num_cores = 1`.
#' @param seed Optional integer random seed. When supplied, the wrapper restores
#'   the caller's previous RNG state before returning.
#'
#' @details
#' By default, the search configuration uses `uncertaintyweights_par = TRUE`
#' when `weighted = TRUE`, so shift-recovery summaries can report both
#' unweighted and IC-weighted metrics.
#'
#' Regardless of the global calibration model used to build `template`, the
#' downstream search is intentionally restricted to intercept-only formulas.
#' This keeps the study focused on shift recovery in the simulated response
#' block rather than re-estimating predictor effects within each replicate. The
#' simulated `trait_data` supplied to each search result therefore contains only
#' the regenerated response variables.
#' Templates calibrated with covariates are supported; those covariates influence
#' the fitted means and residual covariance used to generate each response block,
#' but they are not re-fit within each simulated search replicate.
#'
#' Failed search replicates are retained in the output, but their inferred-shift
#' metrics are recorded as `NA` and excluded from recovery evaluation. Their
#' error messages and candidate counts remain available for diagnosis, and
#' completion and failure rates are reported separately. Reproducibility for
#' replicated studies is controlled by the wrapper-level `seed` argument. A
#' supplied seed deterministically creates separate per-replicate simulation
#' and search streams, so serial and parallel runs use the same replicate-level
#' random draws. Parallel execution also uses `future.seed = TRUE`;
#' per-replicate simulator seeds are intentionally disallowed here.
#' When `num_cores > 1`, namespace-level workers receive only their current
#' replicate and compact shared settings; replicate call records use a literal
#' template placeholder so the complete calibration template is not serialized
#' again inside every simulated object.
#'
#' @return A list of class `bifrost_simulation_study` containing the simulated
#'   datasets, raw search results, per-replicate summaries, a study-level
#'   recovery summary, and the output of [evaluateShiftRecovery()].
#'
#' @seealso [simulateShiftedDataset()], [evaluateShiftRecovery()],
#'   [runFalsePositiveSimulationStudy()]
#'
#' @examples
#' \dontrun{
#' set.seed(1)
#' tr <- ape::rtree(40)
#' X <- matrix(rnorm(40 * 3), ncol = 3)
#' rownames(X) <- tr$tip.label
#' tmpl <- createSimulationTemplate(tr, X, formula = "trait_data ~ 1", method = "LL")
#'
#' recovery_study <- runShiftRecoverySimulationStudy(
#'   tmpl,
#'   n_replicates = 2,
#'   tree_tip_count = 25,
#'   simulation_options = list(
#'     num_shifts = 2,
#'     min_shift_tips = 3,
#'     max_shift_tips = 8,
#'     scale_mode = "proportional"
#'   ),
#'   search_options = list(
#'     formula = "trait_data ~ 1",
#'     min_descendant_tips = 3,
#'     shift_acceptance_threshold = 5,
#'     num_cores = 1,
#'     IC = "GIC",
#'     method = "LL"
#'   ),
#'   num_cores = 1,
#'   seed = 2
#' )
#'
#' recovery_study
#' }
#'
#' @export
runShiftRecoverySimulationStudy <- function(template,
                                            n_replicates,
                                            tree_tip_count = NULL,
                                            simulation_options = list(),
                                            search_options = list(),
                                            fuzzy_distance = 2,
                                            weighted = TRUE,
                                            num_cores = 1,
                                            seed = NULL) {
  if (!inherits(template, "bifrost_simulation_template")) {
    stop("template must be a 'bifrost_simulation_template' object.")
  }
  n_replicates <- .simulation_check_integer_scalar(
    n_replicates,
    name = "n_replicates",
    minimum = 1L,
    message = "n_replicates must be a single integer >= 1."
  )
  tree_tip_count <- .simulation_check_integer_scalar(
    tree_tip_count,
    name = "tree_tip_count",
    minimum = 2L,
    allow_null = TRUE,
    message = "tree_tip_count must be NULL or a single integer >= 2."
  )
  if (!is.null(tree_tip_count) && tree_tip_count > template$n_tips) {
    stop("tree_tip_count cannot exceed the number of tips in the template.")
  }
  fuzzy_distance <- .simulation_check_integer_scalar(
    fuzzy_distance,
    name = "fuzzy_distance",
    minimum = 0L,
    message = "fuzzy_distance must be a single non-negative integer."
  )
  if (!is.logical(weighted) || length(weighted) != 1L || is.na(weighted)) {
    stop("weighted must be TRUE or FALSE.")
  }
  num_cores <- .simulation_check_integer_scalar(
    num_cores,
    name = "num_cores",
    minimum = 1L,
    message = "num_cores must be a single integer >= 1."
  )
  if (!is.list(simulation_options) || !is.list(search_options)) {
    stop("simulation_options and search_options must both be lists.")
  }
  if (!is.null(simulation_options$seed)) {
    stop(
      "Do not supply simulation_options$seed to runShiftRecoverySimulationStudy(); ",
      "use the wrapper-level seed argument instead."
    )
  }
  if (is.null(simulation_options$num_shifts) ||
      is.null(simulation_options$min_shift_tips) ||
      is.null(simulation_options$max_shift_tips)) {
    stop("simulation_options must include num_shifts, min_shift_tips, and max_shift_tips.")
  }
  requested_generator <- .simulation_generator_from_options(simulation_options)
  seed_state <- .simulation_set_seed(seed)
  on.exit(.simulation_restore_seed(seed_state), add = TRUE)
  simulation_seeds <- .simulation_replicate_seeds(
    n_replicates,
    seeded = !is.null(seed)
  )
  search_seeds <- .simulation_replicate_seeds(
    n_replicates,
    seeded = !is.null(seed)
  )

  template_search_formula <- if (!is.null(template$search_formula)) {
    template$search_formula
  } else {
    "trait_data ~ 1"
  }

  search_defaults <- list(
    formula = template_search_formula,
    min_descendant_tips = 10,
    num_cores = 1,
    shift_acceptance_threshold = 10,
    IC = "GIC",
    plot = FALSE,
    store_model_fit_history = FALSE,
    verbose = FALSE,
    progress = FALSE,
    uncertaintyweights = FALSE,
    uncertaintyweights_par = isTRUE(weighted)
  )
  fit_settings <- .simulation_template_fit_settings(template)
  if (!is.null(fit_settings$method)) {
    search_defaults$method <- fit_settings$method
  }
  if (!is.null(fit_settings$error)) {
    search_defaults$error <- fit_settings$error
  }
  search_opts <- utils::modifyList(search_defaults, search_options)
  search_opts$formula <- validateSimulationStudyFormula(search_opts$formula)
  if (isTRUE(num_cores > 1L) && isTRUE(search_opts$num_cores > 1L)) {
    warning(
      "Both wrapper-level and search-level parallelism are > 1; ",
      "forcing each search to use one core so outer replicate parallelism ",
      "takes precedence.",
      call. = FALSE
    )
    search_opts$num_cores <- 1L
  }

  simulate_shifted_dataset_fn <- simulateShiftedDataset
  search_optimal_configuration_fn <- searchOptimalConfiguration
  generate_painted_trees_fn <- generatePaintedTrees
  with_replicate_seed_fn <- .simulation_seed_runner()
  compact_call_fn <- .simulation_study_compact_call

  if (num_cores > 1L) {
    old_plan <- future::plan()
    on.exit(future::plan(old_plan), add = TRUE)
    use_multicore <- .Platform$OS.type != "windows" &&
      !identical(Sys.getenv("RSTUDIO"), "1") &&
      !identical(Sys.getenv("RSTUDIO_SESSION_INITIALIZED"), "1")
    if (use_multicore) {
      future::plan(future::multicore, workers = num_cores)
    } else {
      future::plan(future::multisession, workers = num_cores)
    }

    simdata <- progressr::with_progress({
      progress <- progressr::progressor(steps = n_replicates)
      future.apply::future_lapply(
        seq_len(n_replicates),
        .simulation_study_simulate_shifted,
        template = template,
        tree_tip_count = tree_tip_count,
        simulation_options = simulation_options,
        simulation_seeds = simulation_seeds,
        simulate_dataset = simulate_shifted_dataset_fn,
        with_replicate_seed = with_replicate_seed_fn,
        compact_call = compact_call_fn,
        progress = progress,
        future.seed = TRUE
      )
    })

    search_replicates <- Map(
      function(sim, replicate_seed) list(sim = sim, seed = replicate_seed),
      simdata,
      search_seeds
    )

    results <- progressr::with_progress({
      progress <- progressr::progressor(along = search_replicates)
      future.apply::future_lapply(
        search_replicates,
        .simulation_study_search,
        search_options = search_opts,
        tree_component = "paintedTree",
        data_components = c("trait_data", "simulatedData"),
        search_optimal_configuration = search_optimal_configuration_fn,
        generate_painted_trees = generate_painted_trees_fn,
        with_replicate_seed = with_replicate_seed_fn,
        progress = progress,
        future.seed = TRUE
      )
    })
  } else {
    simdata <- lapply(
      seq_len(n_replicates),
      .simulation_study_simulate_shifted,
      template = template,
      tree_tip_count = tree_tip_count,
      simulation_options = simulation_options,
      simulation_seeds = simulation_seeds,
      simulate_dataset = simulate_shifted_dataset_fn,
      with_replicate_seed = with_replicate_seed_fn,
      compact_call = compact_call_fn
    )

    search_replicates <- Map(
      function(sim, replicate_seed) list(sim = sim, seed = replicate_seed),
      simdata,
      search_seeds
    )

    results <- lapply(
      search_replicates,
      .simulation_study_search,
      search_options = search_opts,
      tree_component = "paintedTree",
      data_components = c("trait_data", "simulatedData"),
      search_optimal_configuration = search_optimal_configuration_fn,
      generate_painted_trees = generate_painted_trees_fn,
      with_replicate_seed = with_replicate_seed_fn
    )
  }

  replicate_generators <- vapply(simdata, function(sim) {
    if (is.null(sim$simulation_generator)) requested_generator else sim$simulation_generator
  }, character(1L))
  simulation_generators <- unique(replicate_generators)

  per_replicate <- do.call(rbind, lapply(seq_along(results), function(i) {
    status <- if (is.null(results[[i]]$error)) "ok" else "error"
    data.frame(
      replicate = i,
      generating_scenario = simdata[[i]]$generating_scenario,
      simulation_generator = replicate_generators[i],
      n_true_shifts = length(simdata[[i]]$shiftNodes),
      n_inferred_shifts = if (identical(status, "ok")) {
        length(results[[i]]$shift_nodes_no_uncertainty)
      } else {
        NA_integer_
      },
      n_candidates = results[[i]]$num_candidates,
      status = status,
      error = if (is.null(results[[i]]$error)) NA_character_ else as.character(results[[i]]$error),
      stringsAsFactors = FALSE
    )
  }))
  rownames(per_replicate) <- NULL

  evaluation <- evaluateShiftRecovery(
    simdata = simdata,
    simresults = results,
    fuzzy_distance = fuzzy_distance,
    weighted = weighted,
    verbose = FALSE
  )

  generating_scenarios <- unique(vapply(simdata, `[[`, character(1), "generating_scenario"))
  study_summary <- list(
    n_replicates = n_replicates,
    n_completed = sum(per_replicate$status == "ok"),
    n_failed = sum(per_replicate$status == "error"),
    simulation_generator = simulation_generators,
    completion_rate = mean(per_replicate$status == "ok"),
    failure_rate = mean(per_replicate$status == "error"),
    n_evaluable_replicates = evaluation$n_evaluable_replicates,
    generating_scenario = generating_scenarios,
    strict = evaluation$strict,
    fuzzy = evaluation$fuzzy,
    weighted = evaluation$weighted
  )

  out <- list(
    user_input = as.list(match.call()),
    study_type = "shift_recovery",
    generating_scenario = study_summary$generating_scenario,
    simulation_generator = simulation_generators,
    simdata = simdata,
    results = results,
    per_replicate = per_replicate,
    evaluation = evaluation,
    study_summary = study_summary,
    simulation_options = simulation_options,
    search_options = search_opts
  )
  class(out) <- c("bifrost_simulation_study", class(out))
  out
}

#' Evaluate Shift Recovery Against Simulated Ground Truth
#'
#' @description
#' Compare inferred shift locations against known simulated shift locations across
#' multiple datasets using both strict node matching and fuzzy matching within a
#' configurable node distance. This is a package-clean export of the manuscript's
#' `evaluate_shift_recovery()` function.
#'
#' @param simdata A list of simulated shifted datasets, typically the `simdata`
#'   component returned by [runShiftRecoverySimulationStudy()]. Each element must
#'   contain at least `paintedTree` and `shiftNodes`.
#' @param simresults A list of search results corresponding to `simdata`, usually
#'   the `results` component returned by [runShiftRecoverySimulationStudy()].
#'   Each element should contain `shift_nodes_no_uncertainty` and
#'   `num_candidates`; if `weighted = TRUE`, `ic_weights` is also used when
#'   available.
#' @param fuzzy_distance Integer node distance threshold used for fuzzy matching.
#' @param weighted Logical; if `TRUE`, compute weighted precision/recall/F1 using
#'   the inferred-node IC weights from each search result when available.
#' @param verbose Logical; if `TRUE`, print a compact summary of the strict,
#'   fuzzy, and weighted metrics.
#'
#' @details
#' The evaluation uses three complementary summaries:
#' \itemize{
#'   \item strict metrics count only exact node matches,
#'   \item fuzzy metrics allow inferred shifts within `fuzzy_distance` nodes of a
#'         true shift, using a greedy one-to-one assignment, and
#'   \item weighted metrics apply IC weights only to inferred nodes, following
#'         the simulation-study implementation.
#' }
#'
#' Search results carrying a non-empty `error` field are retained by the study
#' wrappers for diagnosis but excluded from recovery counts and metrics.
#'
#' When a replicate has no evaluable candidate shifts (`num_candidates == 0`),
#' recall-style quantities can still be computed from the true and inferred
#' shifts, but specificity, false-positive rate, and balanced accuracy are
#' returned as `NA` because there are no evaluable negatives.
#'
#' @return A list of class `bifrost_shift_recovery_evaluation` with five
#'   components:
#' \describe{
#'   \item{`strict`}{Strict precision, recall, F1, specificity, false-positive
#'   rate, and balanced accuracy.}
#'   \item{`fuzzy`}{The same metrics under fuzzy matching.}
#'   \item{`weighted`}{Weighted strict and fuzzy precision/recall/F1 summaries
#'   when `weighted = TRUE`; otherwise `NULL`.}
#'   \item{`counts`}{Aggregated contingency-table counts for the strict and fuzzy
#'   matching schemes.}
#'   \item{`n_evaluable_replicates`}{Number of replicate pairs included in the
#'   aggregated counts after failed or structurally incomplete records are
#'   excluded.}
#' }
#'
#' @seealso [runShiftRecoverySimulationStudy()], [simulateShiftedDataset()]
#'
#' @examples
#' \dontrun{
#' # Usually called on the output of runShiftRecoverySimulationStudy():
#' # evaluateShiftRecovery(study$simdata, study$results, fuzzy_distance = 2)
#' }
#'
#' @export
evaluateShiftRecovery <- function(simdata,
                                  simresults,
                                  fuzzy_distance = 2,
                                  weighted = TRUE,
                                  verbose = TRUE) {
  if (!is.list(simdata) || !is.list(simresults)) {
    stop("simdata and simresults must both be lists.")
  }
  if (length(simdata) != length(simresults)) {
    stop("simdata and simresults must have the same length.")
  }
  fuzzy_distance <- .simulation_check_integer_scalar(
    fuzzy_distance,
    name = "fuzzy_distance",
    minimum = 0L,
    message = "fuzzy_distance must be a single non-negative integer."
  )
  if (!is.logical(weighted) || length(weighted) != 1L || is.na(weighted)) {
    stop("weighted must be TRUE or FALSE.")
  }
  if (!is.logical(verbose) || length(verbose) != 1L || is.na(verbose)) {
    stop("verbose must be TRUE or FALSE.")
  }

  safe_divide <- function(num, den) {
    ifelse(den == 0, NA_real_, num / den)
  }
  harmonic_mean <- function(p, r) {
    ifelse(is.na(p + r) || (p + r) == 0, NA_real_, (2 * p * r) / (p + r))
  }

  strict_counts <- c(TP = 0, FP = 0, FN = 0, TN = 0)
  fuzzy_counts <- c(TP = 0, FP = 0, FN = 0, TN = 0)
  weighted_strict <- c(TP = 0, FP = 0)
  weighted_fuzzy <- c(TP = 0, FP = 0)
  n_evaluable_replicates <- 0L

  for (k in seq_along(simdata)) {
    result_error <- simresults[[k]]$error
    result_failed <- !is.null(result_error) && length(result_error) > 0L &&
      any(!is.na(result_error) & nzchar(as.character(result_error)))
    if (result_failed) {
      next
    }

    true_nodes <- simdata[[k]]$shiftNodes
    inferred_nodes <- simresults[[k]]$shift_nodes_no_uncertainty
    candidate_count <- simresults[[k]]$num_candidates
    tree_k <- simdata[[k]]$paintedTree

    if (is.null(true_nodes) || is.null(inferred_nodes) || is.null(candidate_count) ||
        is.null(tree_k)) {
      next
    }
    n_evaluable_replicates <- n_evaluable_replicates + 1L

    weights <- NULL
    if (weighted && !is.null(simresults[[k]]$ic_weights) &&
        nrow(simresults[[k]]$ic_weights) > 0L) {
      weights <- setNames(
        simresults[[k]]$ic_weights$ic_weight_withshift,
        simresults[[k]]$ic_weights$node
      )
    }

    strict_tp_nodes <- intersect(true_nodes, inferred_nodes)
    strict_fp_nodes <- setdiff(inferred_nodes, true_nodes)
    strict_fn_nodes <- setdiff(true_nodes, inferred_nodes)
    strict_tn <- max(
      candidate_count - length(strict_tp_nodes) -
        length(strict_fp_nodes) - length(strict_fn_nodes),
      0L
    )

    strict_counts <- strict_counts + c(
      TP = length(strict_tp_nodes),
      FP = length(strict_fp_nodes),
      FN = length(strict_fn_nodes),
      TN = strict_tn
    )

    if (weighted && !is.null(weights)) {
      weighted_strict["TP"] <- weighted_strict["TP"] +
        sum(weights[as.character(strict_tp_nodes)], na.rm = TRUE)
      weighted_strict["FP"] <- weighted_strict["FP"] +
        sum(weights[as.character(strict_fp_nodes)], na.rm = TRUE)
    }

    if (length(true_nodes) == 0L || length(inferred_nodes) == 0L) {
      fuzzy_tp <- 0
      fuzzy_fp <- length(inferred_nodes)
      fuzzy_fn <- length(true_nodes)
      matched_inferred <- rep(FALSE, length(inferred_nodes))
    } else {
      distance_matrix <- matrix(
        Inf,
        nrow = length(inferred_nodes),
        ncol = length(true_nodes)
      )

      for (i in seq_along(inferred_nodes)) {
        for (j in seq_along(true_nodes)) {
          node_distance <- tryCatch(
            length(ape::nodepath(tree_k, inferred_nodes[i], true_nodes[j])) - 1L,
            error = function(e) Inf
          )
          if (node_distance <= fuzzy_distance) {
            distance_matrix[i, j] <- node_distance
          }
        }
      }

      matched_inferred <- rep(FALSE, length(inferred_nodes))
      matched_true <- rep(FALSE, length(true_nodes))
      # Manuscript compatibility: consume the nearest available pair greedily,
      # including the original row/column tie-breaking order. Replacing this
      # with maximum-cardinality matching can change recovery counts.
      while (TRUE) {
        min_distance <- min(distance_matrix, na.rm = TRUE)
        if (!is.finite(min_distance)) {
          break
        }
        match_index <- which(distance_matrix == min_distance, arr.ind = TRUE)[1, ]
        matched_inferred[match_index[1]] <- TRUE
        matched_true[match_index[2]] <- TRUE
        distance_matrix[match_index[1], ] <- Inf
        distance_matrix[, match_index[2]] <- Inf
      }

      fuzzy_tp <- sum(matched_inferred)
      fuzzy_fp <- length(inferred_nodes) - fuzzy_tp
      fuzzy_fn <- length(true_nodes) - sum(matched_true)
    }

    fuzzy_tn <- max(candidate_count - fuzzy_tp - fuzzy_fp - fuzzy_fn, 0L)
    fuzzy_counts <- fuzzy_counts + c(TP = fuzzy_tp, FP = fuzzy_fp, FN = fuzzy_fn, TN = fuzzy_tn)

    if (weighted && !is.null(weights)) {
      weighted_fuzzy["TP"] <- weighted_fuzzy["TP"] +
        sum(weights[as.character(inferred_nodes[matched_inferred])], na.rm = TRUE)
      weighted_fuzzy["FP"] <- weighted_fuzzy["FP"] +
        sum(weights[as.character(inferred_nodes[!matched_inferred])], na.rm = TRUE)
    }
  }

  strict_precision <- safe_divide(strict_counts["TP"], strict_counts["TP"] + strict_counts["FP"])
  strict_recall <- safe_divide(strict_counts["TP"], strict_counts["TP"] + strict_counts["FN"])
  strict_specificity <- safe_divide(strict_counts["TN"], strict_counts["TN"] + strict_counts["FP"])
  strict_f1 <- harmonic_mean(strict_precision, strict_recall)
  strict_fpr <- safe_divide(strict_counts["FP"], strict_counts["FP"] + strict_counts["TN"])
  strict_balanced <- if (is.na(strict_recall) || is.na(strict_specificity)) {
    NA_real_
  } else {
    mean(c(strict_recall, strict_specificity))
  }

  fuzzy_precision <- safe_divide(fuzzy_counts["TP"], fuzzy_counts["TP"] + fuzzy_counts["FP"])
  fuzzy_recall <- safe_divide(fuzzy_counts["TP"], fuzzy_counts["TP"] + fuzzy_counts["FN"])
  fuzzy_specificity <- safe_divide(fuzzy_counts["TN"], fuzzy_counts["TN"] + fuzzy_counts["FP"])
  fuzzy_f1 <- harmonic_mean(fuzzy_precision, fuzzy_recall)
  fuzzy_fpr <- safe_divide(fuzzy_counts["FP"], fuzzy_counts["FP"] + fuzzy_counts["TN"])
  fuzzy_balanced <- if (is.na(fuzzy_recall) || is.na(fuzzy_specificity)) {
    NA_real_
  } else {
    mean(c(fuzzy_recall, fuzzy_specificity))
  }

  scalar_metric <- function(x) {
    unname(as.numeric(x))
  }

  strict_precision <- scalar_metric(strict_precision)
  strict_recall <- scalar_metric(strict_recall)
  strict_specificity <- scalar_metric(strict_specificity)
  strict_f1 <- scalar_metric(strict_f1)
  strict_fpr <- scalar_metric(strict_fpr)
  strict_balanced <- scalar_metric(strict_balanced)

  fuzzy_precision <- scalar_metric(fuzzy_precision)
  fuzzy_recall <- scalar_metric(fuzzy_recall)
  fuzzy_specificity <- scalar_metric(fuzzy_specificity)
  fuzzy_f1 <- scalar_metric(fuzzy_f1)
  fuzzy_fpr <- scalar_metric(fuzzy_fpr)
  fuzzy_balanced <- scalar_metric(fuzzy_balanced)

  weighted_metrics <- NULL
  if (weighted) {
    weighted_precision_strict <- safe_divide(
      weighted_strict["TP"],
      weighted_strict["TP"] + weighted_strict["FP"]
    )
    weighted_recall_strict <- safe_divide(
      weighted_strict["TP"],
      strict_counts["TP"] + strict_counts["FN"]
    )
    weighted_f1_strict <- harmonic_mean(weighted_precision_strict, weighted_recall_strict)

    weighted_precision_fuzzy <- safe_divide(
      weighted_fuzzy["TP"],
      weighted_fuzzy["TP"] + weighted_fuzzy["FP"]
    )
    weighted_recall_fuzzy <- safe_divide(
      weighted_fuzzy["TP"],
      fuzzy_counts["TP"] + fuzzy_counts["FN"]
    )
    weighted_f1_fuzzy <- harmonic_mean(weighted_precision_fuzzy, weighted_recall_fuzzy)

    weighted_precision_strict <- scalar_metric(weighted_precision_strict)
    weighted_recall_strict <- scalar_metric(weighted_recall_strict)
    weighted_f1_strict <- scalar_metric(weighted_f1_strict)
    weighted_precision_fuzzy <- scalar_metric(weighted_precision_fuzzy)
    weighted_recall_fuzzy <- scalar_metric(weighted_recall_fuzzy)
    weighted_f1_fuzzy <- scalar_metric(weighted_f1_fuzzy)

    weighted_metrics <- list(
      strict = list(
        precision = weighted_precision_strict,
        recall = weighted_recall_strict,
        f1 = weighted_f1_strict
      ),
      fuzzy = list(
        precision = weighted_precision_fuzzy,
        recall = weighted_recall_fuzzy,
        f1 = weighted_f1_fuzzy
      )
    )
  }

  if (verbose) {
    cat("\nStrict Performance Metrics\n-------------------------\n")
    cat(sprintf("Precision        : %.3f\n", strict_precision))
    cat(sprintf("Recall           : %.3f\n", strict_recall))
    cat(sprintf("F1 Score         : %.3f\n", strict_f1))
    cat(sprintf("Specificity      : %.3f\n", strict_specificity))
    cat(sprintf("False Pos. Rate  : %.3f\n", strict_fpr))
    cat(sprintf("Balanced Accuracy: %.3f\n\n", strict_balanced))

    cat(sprintf("Fuzzy Matching (distance <= %d)\n----------------------------------------\n",
                as.integer(fuzzy_distance)))
    cat(sprintf("Fuzzy Precision        : %.3f\n", fuzzy_precision))
    cat(sprintf("Fuzzy Recall           : %.3f\n", fuzzy_recall))
    cat(sprintf("Fuzzy F1 Score         : %.3f\n", fuzzy_f1))
    cat(sprintf("Fuzzy Specificity      : %.3f\n", fuzzy_specificity))
    cat(sprintf("Fuzzy False Pos. Rate  : %.3f\n", fuzzy_fpr))
    cat(sprintf("Fuzzy Balanced Accuracy: %.3f\n\n", fuzzy_balanced))

    if (weighted && !is.null(weighted_metrics)) {
      cat("Weighted Metrics (weights on inferred nodes)\n")
      cat("-------------------------------------------\n")
      cat(sprintf("Weighted Precision (strict): %.3f\n", weighted_metrics$strict$precision))
      cat(sprintf("Weighted Recall    (strict): %.3f\n", weighted_metrics$strict$recall))
      cat(sprintf("Weighted F1 Score  (strict): %.3f\n", weighted_metrics$strict$f1))
      cat(sprintf("Weighted Precision (fuzzy) : %.3f\n", weighted_metrics$fuzzy$precision))
      cat(sprintf("Weighted Recall    (fuzzy) : %.3f\n", weighted_metrics$fuzzy$recall))
      cat(sprintf("Weighted F1 Score  (fuzzy) : %.3f\n\n", weighted_metrics$fuzzy$f1))
    }
  }

  out <- list(
    strict = list(
      precision = strict_precision,
      recall = strict_recall,
      f1 = strict_f1,
      specificity = strict_specificity,
      fpr = strict_fpr,
      balanced_accuracy = strict_balanced
    ),
    fuzzy = list(
      precision = fuzzy_precision,
      recall = fuzzy_recall,
      f1 = fuzzy_f1,
      specificity = fuzzy_specificity,
      fpr = fuzzy_fpr,
      balanced_accuracy = fuzzy_balanced
    ),
    weighted = weighted_metrics,
    counts = list(
      strict = strict_counts,
      fuzzy = fuzzy_counts
    ),
    n_evaluable_replicates = n_evaluable_replicates
  )
  class(out) <- c("bifrost_shift_recovery_evaluation", class(out))
  out
}
