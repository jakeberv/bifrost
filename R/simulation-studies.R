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
#'   [simulateNullDataset()]. For replicated studies, `simulation_options$seed`
#'   is not allowed; use the wrapper-level `seed` argument instead.
#' @param search_options Named list of arguments passed to
#'   [searchOptimalConfiguration()]. These override the default manuscript-style
#'   intercept-only search settings assembled from `template`.
#' @param num_cores Integer number of workers used across replicate analyses.
#' @param seed Optional integer random seed.
#'
#' @details
#' This function is the package-native port of the manuscript's
#' `run_FP_shift_inference()` workflow. It preserves the main structure of the
#' original analysis while making the outputs stable and package-friendly:
#' simulated datasets, search results, a per-replicate summary table, and a
#' compact study summary are all returned in a single object of class
#' `bifrost_simulation_study`.
#'
#' Regardless of the global calibration model used to build `template`, the
#' downstream search is intentionally restricted to intercept-only formulas.
#' This keeps the simulation study focused on branch-shift detection in the
#' response block, matching the manuscript-style residual-calibration logic.
#'
#' Failed search replicates are retained in the output and counted as zero
#' inferred shifts; the corresponding error messages are recorded in the
#' per-replicate summary. Reproducibility for replicated studies is controlled
#' by the wrapper-level `seed` argument together with `future.seed = TRUE`;
#' per-replicate simulator seeds are intentionally disallowed here.
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
  if (!is.numeric(n_replicates) || length(n_replicates) != 1L ||
      is.na(n_replicates) || n_replicates < 1L) {
    stop("n_replicates must be a single integer >= 1.")
  }
  if (!is.numeric(num_cores) || length(num_cores) != 1L ||
      is.na(num_cores) || num_cores < 1L) {
    stop("num_cores must be a single integer >= 1.")
  }
  if (!is.list(simulation_options) || !is.list(search_options)) {
    stop("simulation_options and search_options must both be lists.")
  }
  if (!is.null(simulation_options$seed)) {
    stop(
      "Do not supply simulation_options$seed to runFalsePositiveSimulationStudy(); ",
      "use the wrapper-level seed argument instead."
    )
  }
  if (!is.null(seed)) {
    set.seed(seed)
  }

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
    verbose = FALSE
  )
  method_setting <- template$fit_method
  if (is.null(method_setting) && !is.null(template$global_model$call$method)) {
    method_setting <- as.character(template$global_model$call$method)
  }
  if (!is.null(method_setting)) {
    search_defaults$method <- method_setting
  }
  error_setting <- template$fit_error
  if (is.null(error_setting) && !is.null(template$global_model$call$error)) {
    error_setting <- tryCatch(
      eval(template$global_model$call$error),
      error = function(e) NULL
    )
  }
  if (!is.null(error_setting)) {
    search_defaults$error <- error_setting
  }
  search_opts <- utils::modifyList(search_defaults, search_options)
  validate_simulation_study_formula_fn <- function(formula) {
    formula_obj <- if (inherits(formula, "formula")) {
      formula
    } else if (is.character(formula) && length(formula) == 1L && !is.na(formula)) {
      stats::as.formula(formula)
    } else {
      stop("formula must be a single character string or a formula object.")
    }
    formula_terms <- stats::terms(formula_obj)
    intercept_only <- length(attr(formula_terms, "term.labels")) == 0L &&
      identical(attr(formula_terms, "intercept"), 1L)
    if (!intercept_only) {
      stop(
        "Simulation studies currently support intercept-only search formulas only. ",
        "Use formula = \"trait_data ~ 1\" (or an equivalent intercept-only response formula)."
      )
    }
    formula
  }
  search_opts$formula <- validate_simulation_study_formula_fn(search_opts$formula)
  if (isTRUE(num_cores > 1L) && isTRUE(search_opts$num_cores > 1L)) {
    warning("Both wrapper-level and search-level parallelism are > 1; nested parallelism may be inefficient.")
  }

  simulate_null_dataset_fn <- simulateNullDataset
  search_optimal_configuration_fn <- searchOptimalConfiguration
  generate_painted_trees_fn <- generatePaintedTrees

  simdata <- lapply(seq_len(n_replicates), function(i) {
    sim_args <- utils::modifyList(
      simulation_options,
      list(
        template = template,
        tree_tip_count = tree_tip_count
      )
    )
    do.call(simulate_null_dataset_fn, sim_args)
  })

  old_plan <- future::plan()
  on.exit(future::plan(old_plan), add = TRUE)
  use_multicore <- .Platform$OS.type != "windows" &&
    !identical(Sys.getenv("RSTUDIO"), "1") &&
    !identical(Sys.getenv("RSTUDIO_SESSION_INITIALIZED"), "1")
  if (num_cores > 1L) {
    if (use_multicore) {
      future::plan(future::multicore, workers = num_cores)
    } else {
      future::plan(future::multisession, workers = num_cores)
    }
  } else {
    future::plan(future::sequential)
  }

  results <- progressr::with_progress({
    progress <- progressr::progressor(along = simdata)
    future.apply::future_lapply(seq_along(simdata), function(i) {
      search_args <- utils::modifyList(
        search_opts,
        list(
          baseline_tree = ape::as.phylo(simdata[[i]]$tree),
          trait_data = simdata[[i]]$data
        )
      )
      out <- tryCatch(
        do.call(search_optimal_configuration_fn, search_args),
        error = function(e) {
          candidate_count <- max(length(generate_painted_trees_fn(
            ape::as.phylo(simdata[[i]]$tree),
            min_tips = search_opts$min_descendant_tips
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
      progress()
      out
    }, future.seed = TRUE)
  })

  per_replicate <- do.call(rbind, lapply(seq_along(results), function(i) {
    n_candidates <- results[[i]]$num_candidates
    n_inferred <- length(results[[i]]$shift_nodes_no_uncertainty)
    data.frame(
      replicate = i,
      generating_scenario = "null",
      n_candidates = n_candidates,
      n_inferred_shifts = n_inferred,
      false_positive_rate = if (is.numeric(n_candidates) && !is.na(n_candidates) && n_candidates > 0) {
        n_inferred / n_candidates
      } else {
        NA_real_
      },
      status = if (is.null(results[[i]]$error)) "ok" else "error",
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
#' recovery performance. The simulation stage supports both the manuscript's
#' proportional generating-model scenario and its correlation-changing
#' non-generating robustness scenario.
#'
#' @param template A `bifrost_simulation_template` returned by
#'   [createSimulationTemplate()].
#' @param n_replicates Integer number of simulation replicates to run.
#' @param tree_tip_count Optional integer tip count for random subtree sampling.
#'   If `NULL`, each replicate uses the full empirical tree.
#' @param simulation_options Named list of additional arguments passed to
#'   [simulateShiftedDataset()]. For replicated studies,
#'   `simulation_options$seed` is not allowed; use the wrapper-level `seed`
#'   argument instead.
#' @param search_options Named list of arguments passed to
#'   [searchOptimalConfiguration()]. These override the default manuscript-style
#'   intercept-only search settings assembled from `template`.
#' @param fuzzy_distance Integer node distance used for fuzzy matching in
#'   [evaluateShiftRecovery()].
#' @param weighted Logical; if `TRUE`, compute weighted recovery summaries using
#'   IC weights when available.
#' @param num_cores Integer number of workers used across replicate analyses.
#' @param seed Optional integer random seed.
#'
#' @details
#' This function is the package-native port of the manuscript's
#' `run_FN_shift_inference()` workflow. To match the manuscript as closely as
#' possible, the default search configuration uses `uncertaintyweights_par = TRUE`
#' when `weighted = TRUE`, so that shift-recovery summaries can report both
#' unweighted and IC-weighted metrics.
#'
#' Regardless of the global calibration model used to build `template`, the
#' downstream search is intentionally restricted to intercept-only formulas.
#' This keeps the study focused on shift recovery in the simulated response
#' block rather than re-estimating predictor effects within each replicate.
#'
#' Failed search replicates are retained in the output and treated as zero-shift
#' recoveries in the per-replicate summary; the corresponding error messages are
#' recorded. Reproducibility for replicated studies is controlled by the
#' wrapper-level `seed` argument together with `future.seed = TRUE`;
#' per-replicate simulator seeds are intentionally disallowed here.
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
  if (!is.numeric(n_replicates) || length(n_replicates) != 1L ||
      is.na(n_replicates) || n_replicates < 1L) {
    stop("n_replicates must be a single integer >= 1.")
  }
  if (!is.numeric(fuzzy_distance) || length(fuzzy_distance) != 1L ||
      is.na(fuzzy_distance) || fuzzy_distance < 0) {
    stop("fuzzy_distance must be a single non-negative number.")
  }
  if (!is.logical(weighted) || length(weighted) != 1L || is.na(weighted)) {
    stop("weighted must be TRUE or FALSE.")
  }
  if (!is.numeric(num_cores) || length(num_cores) != 1L ||
      is.na(num_cores) || num_cores < 1L) {
    stop("num_cores must be a single integer >= 1.")
  }
  if (!is.list(simulation_options) || !is.list(search_options)) {
    stop("simulation_options and search_options must both be lists.")
  }
  if (!is.null(simulation_options$seed)) {
    stop(
      "Do not supply simulation_options$seed to runShiftRecoverySimulationStudy(); ",
      "use the wrapper-level seed argument instead."
    )
  }
  if (!is.null(seed)) {
    set.seed(seed)
  }
  if (is.null(simulation_options$num_shifts) ||
      is.null(simulation_options$min_shift_tips) ||
      is.null(simulation_options$max_shift_tips)) {
    stop("simulation_options must include num_shifts, min_shift_tips, and max_shift_tips.")
  }

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
    uncertaintyweights = FALSE,
    uncertaintyweights_par = isTRUE(weighted)
  )
  method_setting <- template$fit_method
  if (is.null(method_setting) && !is.null(template$global_model$call$method)) {
    method_setting <- as.character(template$global_model$call$method)
  }
  if (!is.null(method_setting)) {
    search_defaults$method <- method_setting
  }
  error_setting <- template$fit_error
  if (is.null(error_setting) && !is.null(template$global_model$call$error)) {
    error_setting <- tryCatch(
      eval(template$global_model$call$error),
      error = function(e) NULL
    )
  }
  if (!is.null(error_setting)) {
    search_defaults$error <- error_setting
  }
  search_opts <- utils::modifyList(search_defaults, search_options)
  validate_simulation_study_formula_fn <- function(formula) {
    formula_obj <- if (inherits(formula, "formula")) {
      formula
    } else if (is.character(formula) && length(formula) == 1L && !is.na(formula)) {
      stats::as.formula(formula)
    } else {
      stop("formula must be a single character string or a formula object.")
    }
    formula_terms <- stats::terms(formula_obj)
    intercept_only <- length(attr(formula_terms, "term.labels")) == 0L &&
      identical(attr(formula_terms, "intercept"), 1L)
    if (!intercept_only) {
      stop(
        "Simulation studies currently support intercept-only search formulas only. ",
        "Use formula = \"trait_data ~ 1\" (or an equivalent intercept-only response formula)."
      )
    }
    formula
  }
  search_opts$formula <- validate_simulation_study_formula_fn(search_opts$formula)
  if (isTRUE(num_cores > 1L) && isTRUE(search_opts$num_cores > 1L)) {
    warning("Both wrapper-level and search-level parallelism are > 1; nested parallelism may be inefficient.")
  }

  simulate_shifted_dataset_fn <- simulateShiftedDataset
  search_optimal_configuration_fn <- searchOptimalConfiguration
  generate_painted_trees_fn <- generatePaintedTrees

  old_plan <- future::plan()
  on.exit(future::plan(old_plan), add = TRUE)
  use_multicore <- .Platform$OS.type != "windows" &&
    !identical(Sys.getenv("RSTUDIO"), "1") &&
    !identical(Sys.getenv("RSTUDIO_SESSION_INITIALIZED"), "1")
  if (num_cores > 1L) {
    if (use_multicore) {
      future::plan(future::multicore, workers = num_cores)
    } else {
      future::plan(future::multisession, workers = num_cores)
    }
  } else {
    future::plan(future::sequential)
  }

  simdata <- progressr::with_progress({
    progress <- progressr::progressor(steps = n_replicates)
    future.apply::future_lapply(seq_len(n_replicates), function(i) {
      sim_args <- utils::modifyList(
        simulation_options,
        list(
          template = template,
          tree_tip_count = tree_tip_count
        )
      )
      out <- do.call(simulate_shifted_dataset_fn, sim_args)
      progress()
      out
    }, future.seed = TRUE)
  })

  results <- progressr::with_progress({
    progress <- progressr::progressor(along = simdata)
    future.apply::future_lapply(seq_along(simdata), function(i) {
      search_args <- utils::modifyList(
        search_opts,
        list(
          baseline_tree = ape::as.phylo(simdata[[i]]$paintedTree),
          trait_data = if (!is.null(simdata[[i]]$trait_data)) {
            simdata[[i]]$trait_data
          } else {
            simdata[[i]]$simulatedData
          }
        )
      )
      out <- tryCatch(
        do.call(search_optimal_configuration_fn, search_args),
        error = function(e) {
          candidate_count <- max(length(generate_painted_trees_fn(
            ape::as.phylo(simdata[[i]]$paintedTree),
            min_tips = search_opts$min_descendant_tips
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
      progress()
      out
    }, future.seed = TRUE)
  })

  per_replicate <- do.call(rbind, lapply(seq_along(results), function(i) {
    data.frame(
      replicate = i,
      generating_scenario = simdata[[i]]$generating_scenario,
      n_true_shifts = length(simdata[[i]]$shiftNodes),
      n_inferred_shifts = length(results[[i]]$shift_nodes_no_uncertainty),
      n_candidates = results[[i]]$num_candidates,
      status = if (is.null(results[[i]]$error)) "ok" else "error",
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
    generating_scenario = if (length(generating_scenarios) == 1L) {
      generating_scenarios
    } else {
      generating_scenarios
    },
    strict = evaluation$strict,
    fuzzy = evaluation$fuzzy,
    weighted = evaluation$weighted
  )

  out <- list(
    user_input = as.list(match.call()),
    study_type = "shift_recovery",
    generating_scenario = study_summary$generating_scenario,
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
#' This function preserves the manuscript evaluation framework as closely as
#' possible:
#' \itemize{
#'   \item strict metrics count only exact node matches,
#'   \item fuzzy metrics allow inferred shifts within `fuzzy_distance` nodes of a
#'         true shift, using a greedy one-to-one assignment, and
#'   \item weighted metrics apply IC weights only to inferred nodes, following
#'         the manuscript implementation.
#' }
#'
#' When a replicate has no evaluable candidate shifts (`num_candidates == 0`),
#' recall-style quantities can still be computed from the true and inferred
#' shifts, but specificity, false-positive rate, and balanced accuracy are
#' returned as `NA` because there are no evaluable negatives.
#'
#' @return An invisible list with four components:
#' \describe{
#'   \item{`strict`}{Strict precision, recall, F1, specificity, false-positive
#'   rate, and balanced accuracy.}
#'   \item{`fuzzy`}{The same metrics under fuzzy matching.}
#'   \item{`weighted`}{Weighted strict and fuzzy precision/recall/F1 summaries
#'   when `weighted = TRUE`; otherwise `NULL`.}
#'   \item{`counts`}{Aggregated contingency-table counts for the strict and fuzzy
#'   matching schemes.}
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
  if (!is.numeric(fuzzy_distance) || length(fuzzy_distance) != 1L ||
      is.na(fuzzy_distance) || fuzzy_distance < 0) {
    stop("fuzzy_distance must be a single non-negative number.")
  }
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

  for (k in seq_along(simdata)) {
    true_nodes <- simdata[[k]]$shiftNodes
    inferred_nodes <- simresults[[k]]$shift_nodes_no_uncertainty
    candidate_count <- simresults[[k]]$num_candidates
    tree_k <- simdata[[k]]$paintedTree

    if (is.null(true_nodes) || is.null(inferred_nodes) || is.null(candidate_count) ||
        is.null(tree_k)) {
      next
    }

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

  invisible(list(
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
    )
  ))
}
