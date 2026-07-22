#' Run a Fixed-IC Search Tuning Grid
#'
#' @description
#' Run a small simulation-based tuning grid for one information-criterion
#' workflow at a time. This helper evaluates combinations of
#' `shift_acceptance_threshold` and `min_descendant_tips` under a fixed `IC`
#' value by combining:
#'
#' \itemize{
#'   \item a null false-positive study,
#'   \item a proportional shift-recovery study, and
#'   \item a non-proportional correlation robustness study.
#' }
#'
#' The goal is to support practical tuning workflows such as "find the best
#' search settings under GIC" or "find the best search settings under BIC"
#' without treating `IC` as just another free parameter in one large
#' optimization.
#'
#' @param template A `bifrost_simulation_template` returned by
#'   [createSimulationTemplate()].
#' @param IC Character scalar, either `"GIC"` or `"BIC"`. The tuning grid is
#'   run for this IC family only.
#' @param shift_acceptance_thresholds Numeric vector of candidate
#'   `shift_acceptance_threshold` values to evaluate.
#' @param min_descendant_tips_values Integer vector of candidate
#'   `min_descendant_tips` values to evaluate.
#' @param tree_tip_count Optional integer tip count passed to the simulation
#'   study wrappers. If `NULL`, each study uses the full empirical tree.
#' @param null_replicates Integer number of null replicates run for each grid
#'   row.
#' @param recovery_replicates Integer number of shifted replicates run for each
#'   grid row in both the proportional and correlation scenarios.
#' @param null_simulation_options Named list forwarded to
#'   [runFalsePositiveSimulationStudy()]. The `"original"` manuscript generator
#'   is used when `simulation_generator` is omitted. Its `simulation_generator`
#'   entry must contain exactly one supported string.
#' @param proportional_simulation_options Named list forwarded to
#'   [runShiftRecoverySimulationStudy()] for the proportional generating-model
#'   scenario. Must include `num_shifts`, `min_shift_tips`, and
#'   `max_shift_tips`. Its `simulation_generator` entry must contain exactly one
#'   supported string.
#' @param correlation_simulation_options Optional named list forwarded to
#'   [runShiftRecoverySimulationStudy()] for the non-proportional robustness
#'   scenario. The argument name is retained for compatibility with the
#'   `scale_mode = "correlation"` label. Under the default `"original"`
#'   generator it controls the published transform; under
#'   `simulation_generator = "empirical"` it controls the integration-rate
#'   trade-off. Its `simulation_generator` entry must contain exactly one
#'   supported string. If
#'   `NULL`, the proportional options are reused and only `scale_mode` is
#'   changed to `"correlation"`.
#' @param base_search_options Named list of search options shared across all
#'   grid rows. The tuning helper overwrites `IC`, `shift_acceptance_threshold`,
#'   and `min_descendant_tips` for each row. Simulation-grid searches are
#'   intentionally intercept-only, so `formula` should remain `"trait_data ~ 1"`
#'   (or an equivalent intercept-only response formula). If omitted, the helper
#'   inherits `template$search_formula`.
#' @param fuzzy_distance Integer node distance passed to
#'   [evaluateShiftRecovery()] inside the shifted-study wrappers.
#' @param weighted Logical; if `TRUE`, request weighted recovery summaries from
#'   [runShiftRecoverySimulationStudy()] and include weighted fuzzy F1 columns
#'   in the output table.
#' @param num_cores Integer number of workers used across grid settings. When
#'   `num_cores > 1`, [runSearchTuningGrid()] parallelizes over settings and
#'   forces the dependent study wrappers and search calls to run serially within
#'   each setting.
#' @param seed Optional integer seed used to derive deterministic per-study
#'   seeds across the entire grid. When supplied, the wrapper restores the
#'   caller's previous RNG state before returning.
#' @param store_studies Logical; if `TRUE`, retain the raw study objects for
#'   every grid row and scenario. If `FALSE`, return only the summary table and
#'   metadata.
#'
#' @details
#' This function is intentionally conservative in scope. It does not try to
#' optimize over `IC`; instead, it assumes that GIC and BIC should be tuned as
#' separate workflows. In typical use, you call it twice, once with
#' `IC = "GIC"` and once with `IC = "BIC"`, then select one recommended setting
#' from each grid with [selectTunedSearchParameters()].
#'
#' The template may come from a richer global calibration model, but each grid
#' row is still evaluated with an intercept-only shift search on the simulated
#' response block. This follows the residual-calibration workflow used by the
#' simulation-study wrappers. Covariates in the calibration model influence the
#' fitted means and residual covariance used for simulation, but they are not
#' re-fit inside each grid replicate.
#'
#' The returned `summary_table` contains one row per grid combination. Null
#' summaries emphasize false-positive behavior, while proportional and
#' correlation-scenario summaries emphasize recovery, including fuzzy balanced
#' accuracy. The corresponding output columns retain their `correlation_`
#' prefixes for backward compatibility.
#' When `seed` is supplied, the
#' `null_seed`, `proportional_seed`, and `correlation_seed` columns record the
#' deterministic per-study seeds used for each setting; otherwise those columns
#' are `NA`. Candidate-set availability among completed searches is tracked via
#' evaluable fractions so that overly strict `min_descendant_tips` settings can
#' be screened out before choosing a final workflow. Completion and failure
#' rates are reported separately and failed searches are excluded from
#' scientific means and proportions.
#'
#' Parallelism is owned by the top-level function that the user calls. In this
#' helper, `num_cores` is interpreted at the setting level, so dependent study
#' wrappers are always called with `num_cores = 1` and their embedded search
#' calls are forced to `num_cores = 1` as well. This avoids nested parallelism
#' while keeping the user-facing API simple.
#'
#' @return A list of class `bifrost_search_tuning_grid` with components:
#' \describe{
#'   \item{`IC`}{The fixed IC family used across the grid.}
#'   \item{`grid`}{The evaluated combinations of thresholds and minimum clade
#'   sizes.}
#'   \item{`summary_table`}{A data frame with one row per setting, deterministic
#'   per-study seed columns, and summary metrics from the null, proportional,
#'   and correlation-scenario studies.}
#'   \item{`simulation_generators`}{A named character vector recording the null,
#'   proportional, and non-proportional simulation generators.}
#'   \item{`studies`}{Either `NULL` or a list of raw study objects for each grid
#'   row and scenario, depending on `store_studies`.}
#'   \item{`base_search_options`}{The shared search options used before the
#'   tuned grid parameters were injected.}
#' }
#'
#' @seealso [runFalsePositiveSimulationStudy()],
#'   [runShiftRecoverySimulationStudy()], [selectTunedSearchParameters()]
#'
#' @examples
#' \dontrun{
#' set.seed(1)
#' tr <- ape::rtree(40)
#' X <- matrix(rnorm(40 * 3), ncol = 3)
#' rownames(X) <- tr$tip.label
#' tmpl <- createSimulationTemplate(tr, X, formula = "trait_data ~ 1", method = "LL")
#'
#' gic_grid <- runSearchTuningGrid(
#'   template = tmpl,
#'   IC = "GIC",
#'   shift_acceptance_thresholds = c(2, 10),
#'   min_descendant_tips_values = c(5, 10),
#'   null_replicates = 10,
#'   recovery_replicates = 10,
#'   proportional_simulation_options = list(
#'     num_shifts = 2,
#'     min_shift_tips = 5,
#'     max_shift_tips = 10
#'   ),
#'   base_search_options = list(
#'     formula = "trait_data ~ 1",
#'     method = "LL"
#'   ),
#'   num_cores = 1,
#'   seed = 1
#' )
#'
#' head(gic_grid$summary_table)
#' }
#'
#' @export
runSearchTuningGrid <- function(template,
                                IC = c("GIC", "BIC"),
                                shift_acceptance_thresholds,
                                min_descendant_tips_values,
                                tree_tip_count = NULL,
                                null_replicates = 50,
                                recovery_replicates = 50,
                                null_simulation_options = list(),
                                proportional_simulation_options,
                                correlation_simulation_options = NULL,
                                base_search_options = list(),
                                fuzzy_distance = 2,
                                weighted = TRUE,
                                num_cores = 1,
                                seed = NULL,
                                store_studies = FALSE) {
  IC <- match.arg(IC)

  if (!inherits(template, "bifrost_simulation_template")) {
    stop("template must be a 'bifrost_simulation_template' object.")
  }
  if (missing(shift_acceptance_thresholds) ||
      !is.numeric(shift_acceptance_thresholds) ||
      length(shift_acceptance_thresholds) < 1L ||
      anyNA(shift_acceptance_thresholds)) {
    stop("shift_acceptance_thresholds must be a numeric vector with at least one value.")
  }
  if (missing(min_descendant_tips_values) ||
      is.null(min_descendant_tips_values)) {
    stop("min_descendant_tips_values must be a numeric vector of integers >= 1.")
  }
  min_descendant_tips_values <- .simulation_check_integer_vector(
    min_descendant_tips_values,
    name = "min_descendant_tips_values",
    minimum = 1L,
    message = "min_descendant_tips_values must be a numeric vector of integers >= 1."
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
  null_replicates <- .simulation_check_integer_scalar(
    null_replicates,
    name = "null_replicates",
    minimum = 1L,
    message = "null_replicates must be a single integer >= 1."
  )
  recovery_replicates <- .simulation_check_integer_scalar(
    recovery_replicates,
    name = "recovery_replicates",
    minimum = 1L,
    message = "recovery_replicates must be a single integer >= 1."
  )
  if (!is.list(null_simulation_options) ||
      missing(proportional_simulation_options) ||
      !is.list(proportional_simulation_options) ||
      !is.null(correlation_simulation_options) && !is.list(correlation_simulation_options) ||
      !is.list(base_search_options)) {
    stop(
      "null_simulation_options, proportional_simulation_options, ",
      "correlation_simulation_options, and base_search_options must be lists."
    )
  }
  if (is.null(proportional_simulation_options$num_shifts) ||
      is.null(proportional_simulation_options$min_shift_tips) ||
      is.null(proportional_simulation_options$max_shift_tips)) {
    stop(
      "proportional_simulation_options must include num_shifts, ",
      "min_shift_tips, and max_shift_tips."
    )
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
  if (!is.logical(store_studies) || length(store_studies) != 1L || is.na(store_studies)) {
    stop("store_studies must be TRUE or FALSE.")
  }
  proportional_options <- proportional_simulation_options
  proportional_options$scale_mode <- "proportional"

  correlation_options <- if (is.null(correlation_simulation_options)) {
    proportional_options
  } else {
    utils::modifyList(proportional_simulation_options, correlation_simulation_options)
  }
  correlation_options$scale_mode <- "correlation"
  simulation_generators <- c(
    null = .simulation_generator_from_options(null_simulation_options),
    proportional = .simulation_generator_from_options(proportional_options),
    correlation = .simulation_generator_from_options(correlation_options)
  )
  seed_state <- .simulation_set_seed(seed)
  on.exit(.simulation_restore_seed(seed_state), add = TRUE)

  if (is.null(base_search_options$formula)) {
    base_search_options$formula <- if (!is.null(template$search_formula)) {
      template$search_formula
    } else {
      "trait_data ~ 1"
    }
  }
  base_search_options$formula <- validateSimulationStudyFormula(base_search_options$formula)

  threshold_values <- unique(as.numeric(shift_acceptance_thresholds))
  min_tip_values <- unique(min_descendant_tips_values)
  grid <- expand.grid(
    shift_acceptance_threshold = threshold_values,
    min_descendant_tips = min_tip_values,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  grid$setting_id <- seq_len(nrow(grid))
  grid <- grid[, c("setting_id", "shift_acceptance_threshold", "min_descendant_tips")]

  setting_seeds <- if (!is.null(seed_state)) {
    matrix(
      sample.int(.Machine$integer.max, size = nrow(grid) * 3L),
      ncol = 3L,
      byrow = TRUE,
      dimnames = list(NULL, c("null", "proportional", "correlation"))
    )
  } else {
    NULL
  }

  run_false_positive_study_fn <- runFalsePositiveSimulationStudy
  run_shift_recovery_study_fn <- runShiftRecoverySimulationStudy

  evaluate_setting <- function(i) {
    tuning_search_options <- utils::modifyList(
      base_search_options,
      list(
        IC = IC,
        shift_acceptance_threshold = grid$shift_acceptance_threshold[i],
        min_descendant_tips = grid$min_descendant_tips[i],
        num_cores = 1L,
        progress = FALSE
      )
    )

    null_study <- run_false_positive_study_fn(
      template = template,
      n_replicates = null_replicates,
      tree_tip_count = tree_tip_count,
      simulation_options = null_simulation_options,
      search_options = tuning_search_options,
      num_cores = 1L,
      seed = if (is.null(setting_seeds)) NULL else setting_seeds[i, "null"]
    )

    proportional_study <- run_shift_recovery_study_fn(
      template = template,
      n_replicates = recovery_replicates,
      tree_tip_count = tree_tip_count,
      simulation_options = proportional_options,
      search_options = tuning_search_options,
      fuzzy_distance = fuzzy_distance,
      weighted = weighted,
      num_cores = 1L,
      seed = if (is.null(setting_seeds)) NULL else setting_seeds[i, "proportional"]
    )

    correlation_study <- run_shift_recovery_study_fn(
      template = template,
      n_replicates = recovery_replicates,
      tree_tip_count = tree_tip_count,
      simulation_options = correlation_options,
      search_options = tuning_search_options,
      fuzzy_distance = fuzzy_distance,
      weighted = weighted,
      num_cores = 1L,
      seed = if (is.null(setting_seeds)) NULL else setting_seeds[i, "correlation"]
    )

    completed_mean <- function(x, status) {
      completed <- status == "ok"
      if (!any(completed)) {
        return(NA_real_)
      }
      mean(x[completed])
    }
    completion_rate <- function(status) mean(status == "ok")
    failure_rate <- function(status) mean(status == "error")

    null_status <- null_study$per_replicate$status
    proportional_status <- proportional_study$per_replicate$status
    correlation_status <- correlation_study$per_replicate$status
    null_any_fp <- completed_mean(
      null_study$per_replicate$n_inferred_shifts > 0,
      null_status
    )
    null_evaluable_fraction <- completed_mean(
      !is.na(null_study$per_replicate$false_positive_rate),
      null_status
    )
    proportional_evaluable_fraction <- completed_mean(
      proportional_study$per_replicate$n_candidates > 0,
      proportional_status
    )
    correlation_evaluable_fraction <- completed_mean(
      correlation_study$per_replicate$n_candidates > 0,
      correlation_status
    )
    recovery_metric <- function(study, match, metric) {
      value <- study$evaluation[[match]][[metric]]
      if (is.null(value)) {
        return(NA_real_)
      }
      unname(as.numeric(value))
    }

    summary_row <- data.frame(
      setting_id = grid$setting_id[i],
      IC = IC,
      shift_acceptance_threshold = grid$shift_acceptance_threshold[i],
      min_descendant_tips = grid$min_descendant_tips[i],
      null_seed = if (is.null(setting_seeds)) NA_integer_ else setting_seeds[i, "null"],
      proportional_seed = if (is.null(setting_seeds)) NA_integer_ else setting_seeds[i, "proportional"],
      correlation_seed = if (is.null(setting_seeds)) NA_integer_ else setting_seeds[i, "correlation"],
      null_mean_false_positive_rate = unname(as.numeric(null_study$study_summary$mean_false_positive_rate)),
      null_fraction_any_false_positive = null_any_fp,
      null_evaluable_fraction = null_evaluable_fraction,
      null_completion_rate = completion_rate(null_status),
      null_failure_rate = failure_rate(null_status),
      null_mean_inferred_shifts = completed_mean(
        null_study$per_replicate$n_inferred_shifts,
        null_status
      ),
      proportional_strict_precision = recovery_metric(proportional_study, "strict", "precision"),
      proportional_strict_recall = recovery_metric(proportional_study, "strict", "recall"),
      proportional_strict_f1 = recovery_metric(proportional_study, "strict", "f1"),
      proportional_strict_specificity = recovery_metric(proportional_study, "strict", "specificity"),
      proportional_strict_fpr = recovery_metric(proportional_study, "strict", "fpr"),
      proportional_strict_balanced_accuracy = recovery_metric(proportional_study, "strict", "balanced_accuracy"),
      proportional_fuzzy_precision = recovery_metric(proportional_study, "fuzzy", "precision"),
      proportional_fuzzy_recall = recovery_metric(proportional_study, "fuzzy", "recall"),
      proportional_fuzzy_f1 = recovery_metric(proportional_study, "fuzzy", "f1"),
      proportional_fuzzy_specificity = recovery_metric(proportional_study, "fuzzy", "specificity"),
      proportional_fuzzy_fpr = recovery_metric(proportional_study, "fuzzy", "fpr"),
      proportional_fuzzy_balanced_accuracy = recovery_metric(proportional_study, "fuzzy", "balanced_accuracy"),
      proportional_evaluable_fraction = proportional_evaluable_fraction,
      proportional_completion_rate = completion_rate(proportional_status),
      proportional_failure_rate = failure_rate(proportional_status),
      proportional_mean_inferred_shifts = completed_mean(
        proportional_study$per_replicate$n_inferred_shifts,
        proportional_status
      ),
      correlation_strict_precision = recovery_metric(correlation_study, "strict", "precision"),
      correlation_strict_recall = recovery_metric(correlation_study, "strict", "recall"),
      correlation_strict_f1 = recovery_metric(correlation_study, "strict", "f1"),
      correlation_strict_specificity = recovery_metric(correlation_study, "strict", "specificity"),
      correlation_strict_fpr = recovery_metric(correlation_study, "strict", "fpr"),
      correlation_strict_balanced_accuracy = recovery_metric(correlation_study, "strict", "balanced_accuracy"),
      correlation_fuzzy_precision = recovery_metric(correlation_study, "fuzzy", "precision"),
      correlation_fuzzy_recall = recovery_metric(correlation_study, "fuzzy", "recall"),
      correlation_fuzzy_f1 = recovery_metric(correlation_study, "fuzzy", "f1"),
      correlation_fuzzy_specificity = recovery_metric(correlation_study, "fuzzy", "specificity"),
      correlation_fuzzy_fpr = recovery_metric(correlation_study, "fuzzy", "fpr"),
      correlation_fuzzy_balanced_accuracy = recovery_metric(correlation_study, "fuzzy", "balanced_accuracy"),
      correlation_evaluable_fraction = correlation_evaluable_fraction,
      correlation_completion_rate = completion_rate(correlation_status),
      correlation_failure_rate = failure_rate(correlation_status),
      correlation_mean_inferred_shifts = completed_mean(
        correlation_study$per_replicate$n_inferred_shifts,
        correlation_status
      ),
      stringsAsFactors = FALSE
    )

    if (weighted) {
      summary_row$proportional_weighted_fuzzy_f1 <- if (!is.null(proportional_study$evaluation$weighted)) {
        proportional_study$evaluation$weighted$fuzzy$f1
      } else {
        NA_real_
      }
      summary_row$correlation_weighted_fuzzy_f1 <- if (!is.null(correlation_study$evaluation$weighted)) {
        correlation_study$evaluation$weighted$fuzzy$f1
      } else {
        NA_real_
      }
    }

    list(
      summary_row = summary_row,
      studies = if (store_studies) {
        list(
          null = null_study,
          proportional = proportional_study,
          correlation = correlation_study
        )
      } else {
        NULL
      }
    )
  }

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

    setting_results <- progressr::with_progress({
      progress <- progressr::progressor(along = seq_len(nrow(grid)))
      future.apply::future_lapply(
        seq_len(nrow(grid)),
        function(i) {
          on.exit(progress(), add = TRUE)
          evaluate_setting(i)
        },
        future.seed = TRUE
      )
    })
  } else {
    setting_results <- lapply(seq_len(nrow(grid)), evaluate_setting)
  }

  studies <- if (store_studies) lapply(setting_results, `[[`, "studies") else NULL
  summary_table <- do.call(rbind, lapply(setting_results, `[[`, "summary_row"))
  rownames(summary_table) <- NULL

  out <- list(
    user_input = as.list(match.call()),
    IC = IC,
    grid = grid,
    summary_table = summary_table,
    studies = studies,
    simulation_generators = simulation_generators,
    base_search_options = base_search_options,
    null_replicates = as.integer(null_replicates),
    recovery_replicates = as.integer(recovery_replicates),
    fuzzy_distance = fuzzy_distance,
    weighted = weighted,
    store_studies = store_studies
  )
  class(out) <- c("bifrost_search_tuning_grid", class(out))
  out
}

#' Select Tuned Search Parameters from a Fixed-IC Grid
#'
#' @description
#' Choose a recommended `shift_acceptance_threshold` and `min_descendant_tips`
#' combination from a `bifrost_search_tuning_grid` produced by
#' [runSearchTuningGrid()]. Selection happens within a single IC family, making
#' it easy to report one tuned workflow for GIC and another for BIC.
#'
#' @param tuning_grid A `bifrost_search_tuning_grid` returned by
#'   [runSearchTuningGrid()].
#' @param max_false_positive_rate Maximum acceptable mean false-positive rate
#'   under the null study.
#' @param max_any_false_positive Maximum acceptable fraction of null replicates
#'   that infer at least one shift.
#' @param min_evaluable_fraction Minimum acceptable fraction of replicates with
#'   at least one candidate shift. This filter is applied to the null,
#'   proportional, and integration-rate summaries.
#' @param primary_metric Character scalar indicating the recovery metric used to
#'   rank feasible settings. Defaults to `"fuzzy_balanced_accuracy"`.
#'   Supported legacy values are `"fuzzy_f1"`, `"fuzzy_recall"`, `"strict_f1"`,
#'   and `"weighted_fuzzy_f1"`.
#' @param scenario_weights Numeric length-2 vector giving the weights applied to
#'   the proportional and integration-rate scenarios when forming the ranking
#'   score. The internal weight name `correlation` is retained for compatibility.
#'   If unnamed, the first value is used for proportional and the second for
#'   the integration-rate scenario.
#' @param tie_break Character scalar. `"conservative"` prefers larger thresholds
#'   and larger `min_descendant_tips` when the primary ranking score is tied;
#'   `"liberal"` prefers the smaller values.
#'
#' @details
#' The selector uses a two-step decision rule. First, it filters out settings
#' that violate the null false-positive or evaluability constraints. Second, it
#' ranks the remaining settings using a weighted average of the chosen recovery
#' metric across the proportional and integration-rate scenarios. By default,
#' this ranking metric is fuzzy balanced accuracy.
#'
#' Only settings with finite recovery metrics for every scenario having positive
#' weight are rankable. If no rankable settings pass the filters, the function
#' warns, falls back to the full rankable grid, and records that no feasible
#' settings were available under the supplied constraints. If the grid contains
#' no rankable settings, the function stops rather than returning an unsupported
#' recommendation.
#'
#' @return A list of class `bifrost_search_tuning_selection` with components:
#' \describe{
#'   \item{`selected_row`}{The chosen row from `tuning_grid$summary_table`,
#'   augmented with a ranking score.}
#'   \item{`recommended_search_options`}{A search-options list containing the
#'   tuned controls from the simulation grid. For response-only empirical
#'   searches, the inherited intercept-only formula can usually be used directly.
#'   For empirical analyses with covariates, carry over the tuned controls but set
#'   `formula` to the empirical analysis formula.}
#'   \item{`feasible_table`}{The filtered and ranked table used for selection.}
#'   \item{`n_feasible_settings`}{The number of settings that passed the
#'   supplied constraints.}
#'   \item{`used_all_settings`}{Logical indicating whether the selector had to
#'   fall back to ranking the full grid because no feasible settings remained.}
#' }
#'
#' @seealso [runSearchTuningGrid()], [searchOptimalConfiguration()]
#'
#' @examples
#' \dontrun{
#' # Usually called after runSearchTuningGrid():
#' # tuned <- selectTunedSearchParameters(gic_grid)
#' # tuned$recommended_search_options
#' }
#'
#' @export
selectTunedSearchParameters <- function(tuning_grid,
                                        max_false_positive_rate = 0.05,
                                        max_any_false_positive = 0.2,
                                        min_evaluable_fraction = 0.9,
                                        primary_metric = c(
                                          "fuzzy_balanced_accuracy",
                                          "fuzzy_f1",
                                          "fuzzy_recall",
                                          "strict_f1",
                                          "weighted_fuzzy_f1"
                                        ),
                                        scenario_weights = c(proportional = 0.5, correlation = 0.5),
                                        tie_break = c("conservative", "liberal")) {
  primary_metric <- match.arg(primary_metric)
  tie_break <- match.arg(tie_break)

  if (!inherits(tuning_grid, "bifrost_search_tuning_grid")) {
    stop("tuning_grid must be a 'bifrost_search_tuning_grid' object.")
  }
  if (!is.numeric(max_false_positive_rate) || length(max_false_positive_rate) != 1L ||
      is.na(max_false_positive_rate) || max_false_positive_rate < 0) {
    stop("max_false_positive_rate must be a single non-negative number.")
  }
  if (!is.numeric(max_any_false_positive) || length(max_any_false_positive) != 1L ||
      is.na(max_any_false_positive) || max_any_false_positive < 0 ||
      max_any_false_positive > 1) {
    stop("max_any_false_positive must be a single number between 0 and 1.")
  }
  if (!is.numeric(min_evaluable_fraction) || length(min_evaluable_fraction) != 1L ||
      is.na(min_evaluable_fraction) || min_evaluable_fraction < 0 ||
      min_evaluable_fraction > 1) {
    stop("min_evaluable_fraction must be a single number between 0 and 1.")
  }
  if (!is.numeric(scenario_weights) || length(scenario_weights) != 2L ||
      anyNA(scenario_weights) || any(!is.finite(scenario_weights)) ||
      any(scenario_weights < 0) ||
      sum(scenario_weights) <= 0) {
    stop("scenario_weights must be a non-negative numeric vector of length 2.")
  }

  metric_suffix <- primary_metric

  proportional_metric <- paste0("proportional_", metric_suffix)
  correlation_metric <- paste0("correlation_", metric_suffix)
  required_columns <- c(
    "shift_acceptance_threshold",
    "min_descendant_tips",
    "null_mean_false_positive_rate",
    "null_fraction_any_false_positive",
    "null_evaluable_fraction",
    "proportional_evaluable_fraction",
    "correlation_evaluable_fraction",
    proportional_metric,
    correlation_metric
  )
  missing_columns <- setdiff(required_columns, names(tuning_grid$summary_table))
  if (length(missing_columns) > 0L) {
    stop(
      "tuning_grid$summary_table is missing required columns: ",
      paste(missing_columns, collapse = ", "),
      "."
    )
  }

  weights_named <- if (is.null(names(scenario_weights)) || any(names(scenario_weights) == "")) {
    stats::setNames(as.numeric(scenario_weights), c("proportional", "correlation"))
  } else {
    scenario_weights
  }
  if (!all(c("proportional", "correlation") %in% names(weights_named))) {
    stop("scenario_weights must be named 'proportional' and 'correlation', or unnamed length-2.")
  }
  weights_named <- weights_named[c("proportional", "correlation")]
  weights_named <- weights_named / sum(weights_named)

  summary_table <- tuning_grid$summary_table
  score <- rep(0, nrow(summary_table))
  rankable <- rep(TRUE, nrow(summary_table))
  if (weights_named[["proportional"]] > 0) {
    proportional_values <- summary_table[[proportional_metric]]
    rankable <- rankable & is.finite(proportional_values)
    score <- score + weights_named[["proportional"]] * proportional_values
  }
  if (weights_named[["correlation"]] > 0) {
    correlation_values <- summary_table[[correlation_metric]]
    rankable <- rankable & is.finite(correlation_values)
    score <- score + weights_named[["correlation"]] * correlation_values
  }
  rankable <- rankable & is.finite(score)
  rankable_table <- summary_table[rankable, , drop = FALSE]
  rankable_table$score <- score[rankable]
  if (nrow(rankable_table) == 0L) {
    stop(
      "No tuning settings have finite recovery evidence for every scenario ",
      "with positive weight; no recommendation can be made.",
      call. = FALSE
    )
  }

  finite_filter_values <- is.finite(rankable_table$null_mean_false_positive_rate) &
    is.finite(rankable_table$null_fraction_any_false_positive) &
    is.finite(rankable_table$null_evaluable_fraction) &
    is.finite(rankable_table$proportional_evaluable_fraction) &
    is.finite(rankable_table$correlation_evaluable_fraction)
  feasible <- rankable_table[
    finite_filter_values &
      rankable_table$null_mean_false_positive_rate <= max_false_positive_rate &
      rankable_table$null_fraction_any_false_positive <= max_any_false_positive &
      rankable_table$null_evaluable_fraction >= min_evaluable_fraction &
      rankable_table$proportional_evaluable_fraction >= min_evaluable_fraction &
      rankable_table$correlation_evaluable_fraction >= min_evaluable_fraction,
    ,
    drop = FALSE
  ]

  n_feasible_settings <- nrow(feasible)
  used_all_settings <- FALSE
  if (nrow(feasible) == 0L) {
    warning(
      "No rankable tuning settings met the false-positive and evaluability ",
      "constraints; ranking the full rankable grid instead.",
      call. = FALSE
    )
    feasible <- rankable_table
    used_all_settings <- TRUE
  }

  order_score <- feasible$score
  order_null_fp <- ifelse(is.na(feasible$null_mean_false_positive_rate), Inf, feasible$null_mean_false_positive_rate)
  order_null_any <- ifelse(is.na(feasible$null_fraction_any_false_positive), Inf, feasible$null_fraction_any_false_positive)

  if (tie_break == "conservative") {
    ord <- order(
      -order_score,
      order_null_fp,
      order_null_any,
      -feasible$shift_acceptance_threshold,
      -feasible$min_descendant_tips
    )
  } else {
    ord <- order(
      -order_score,
      order_null_fp,
      order_null_any,
      feasible$shift_acceptance_threshold,
      feasible$min_descendant_tips
    )
  }
  feasible <- feasible[ord, , drop = FALSE]
  rownames(feasible) <- NULL

  selected_row <- feasible[1, , drop = FALSE]
  recommended_search_options <- utils::modifyList(
    tuning_grid$base_search_options,
    list(
      IC = tuning_grid$IC,
      shift_acceptance_threshold = selected_row$shift_acceptance_threshold[[1L]],
      min_descendant_tips = selected_row$min_descendant_tips[[1L]]
    )
  )

  out <- list(
    IC = tuning_grid$IC,
    primary_metric = primary_metric,
    scenario_weights = weights_named,
    selected_row = selected_row,
    recommended_search_options = recommended_search_options,
    feasible_table = feasible,
    n_feasible_settings = n_feasible_settings,
    used_all_settings = used_all_settings,
    filters = list(
      max_false_positive_rate = max_false_positive_rate,
      max_any_false_positive = max_any_false_positive,
      min_evaluable_fraction = min_evaluable_fraction,
      tie_break = tie_break
    )
  )
  class(out) <- c("bifrost_search_tuning_selection", class(out))
  out
}
