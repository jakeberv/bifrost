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
#'   \item a correlation-based robustness study.
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
#'   [runFalsePositiveSimulationStudy()].
#' @param proportional_simulation_options Named list forwarded to
#'   [runShiftRecoverySimulationStudy()] for the proportional generating-model
#'   scenario. Must include `num_shifts`, `min_shift_tips`, and
#'   `max_shift_tips`.
#' @param correlation_simulation_options Optional named list forwarded to
#'   [runShiftRecoverySimulationStudy()] for the correlation robustness
#'   scenario. If `NULL`, the proportional options are reused and only
#'   `scale_mode` is changed to `"correlation"`.
#' @param base_search_options Named list of search options shared across all
#'   grid rows. The tuning helper overwrites `IC`, `shift_acceptance_threshold`,
#'   and `min_descendant_tips` for each row.
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
#'   seeds across the entire grid.
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
#' The returned `summary_table` contains one row per grid combination. Null
#' summaries emphasize false-positive behavior, while proportional and
#' correlation summaries emphasize recovery. Candidate-set availability is
#' tracked separately via evaluable fractions so that overly strict
#' `min_descendant_tips` settings can be screened out before choosing a final
#' workflow.
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
#'   \item{`summary_table`}{A data frame with one row per setting and summary
#'   metrics from the null, proportional, and correlation studies.}
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
      !is.numeric(min_descendant_tips_values) ||
      length(min_descendant_tips_values) < 1L ||
      anyNA(min_descendant_tips_values) ||
      any(min_descendant_tips_values < 1L)) {
    stop("min_descendant_tips_values must be a numeric vector of integers >= 1.")
  }
  if (!is.numeric(null_replicates) || length(null_replicates) != 1L ||
      is.na(null_replicates) || null_replicates < 1L) {
    stop("null_replicates must be a single integer >= 1.")
  }
  if (!is.numeric(recovery_replicates) || length(recovery_replicates) != 1L ||
      is.na(recovery_replicates) || recovery_replicates < 1L) {
    stop("recovery_replicates must be a single integer >= 1.")
  }
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
  if (!is.logical(store_studies) || length(store_studies) != 1L || is.na(store_studies)) {
    stop("store_studies must be TRUE or FALSE.")
  }

  threshold_values <- unique(as.numeric(shift_acceptance_thresholds))
  min_tip_values <- unique(as.integer(min_descendant_tips_values))
  grid <- expand.grid(
    shift_acceptance_threshold = threshold_values,
    min_descendant_tips = min_tip_values,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  grid$setting_id <- seq_len(nrow(grid))
  grid <- grid[, c("setting_id", "shift_acceptance_threshold", "min_descendant_tips")]

  proportional_options <- proportional_simulation_options
  proportional_options$scale_mode <- "proportional"

  correlation_options <- if (is.null(correlation_simulation_options)) {
    proportional_options
  } else {
    utils::modifyList(proportional_simulation_options, correlation_simulation_options)
  }
  correlation_options$scale_mode <- "correlation"

  setting_seeds <- if (!is.null(seed)) {
    set.seed(seed)
    matrix(
      sample.int(.Machine$integer.max, size = nrow(grid) * 3L),
      ncol = 3L,
      byrow = TRUE,
      dimnames = list(NULL, c("null", "proportional", "correlation"))
    )
  } else {
    NULL
  }

  evaluate_setting <- function(i) {
    tuning_search_options <- utils::modifyList(
      base_search_options,
      list(
        IC = IC,
        shift_acceptance_threshold = grid$shift_acceptance_threshold[i],
        min_descendant_tips = grid$min_descendant_tips[i],
        num_cores = 1L
      )
    )

    null_study <- runFalsePositiveSimulationStudy(
      template = template,
      n_replicates = null_replicates,
      tree_tip_count = tree_tip_count,
      simulation_options = null_simulation_options,
      search_options = tuning_search_options,
      num_cores = 1L,
      seed = if (is.null(setting_seeds)) NULL else setting_seeds[i, "null"]
    )

    proportional_study <- runShiftRecoverySimulationStudy(
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

    correlation_study <- runShiftRecoverySimulationStudy(
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

    null_any_fp <- mean(null_study$per_replicate$n_inferred_shifts > 0)
    null_evaluable_fraction <- mean(!is.na(null_study$per_replicate$false_positive_rate))
    proportional_evaluable_fraction <- mean(proportional_study$per_replicate$n_candidates > 0)
    correlation_evaluable_fraction <- mean(correlation_study$per_replicate$n_candidates > 0)

    summary_row <- data.frame(
      setting_id = grid$setting_id[i],
      IC = IC,
      shift_acceptance_threshold = grid$shift_acceptance_threshold[i],
      min_descendant_tips = grid$min_descendant_tips[i],
      null_mean_false_positive_rate = unname(as.numeric(null_study$study_summary$mean_false_positive_rate)),
      null_fraction_any_false_positive = null_any_fp,
      null_evaluable_fraction = null_evaluable_fraction,
      null_mean_inferred_shifts = mean(null_study$per_replicate$n_inferred_shifts),
      proportional_strict_f1 = proportional_study$evaluation$strict$f1,
      proportional_fuzzy_f1 = proportional_study$evaluation$fuzzy$f1,
      proportional_fuzzy_recall = proportional_study$evaluation$fuzzy$recall,
      proportional_evaluable_fraction = proportional_evaluable_fraction,
      proportional_mean_inferred_shifts = mean(proportional_study$per_replicate$n_inferred_shifts),
      correlation_strict_f1 = correlation_study$evaluation$strict$f1,
      correlation_fuzzy_f1 = correlation_study$evaluation$fuzzy$f1,
      correlation_fuzzy_recall = correlation_study$evaluation$fuzzy$recall,
      correlation_evaluable_fraction = correlation_evaluable_fraction,
      correlation_mean_inferred_shifts = mean(correlation_study$per_replicate$n_inferred_shifts),
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

  old_plan <- future::plan()
  on.exit(future::plan(old_plan), add = TRUE)
  if (num_cores > 1L) {
    future::plan(future::multisession, workers = num_cores)
  } else {
    future::plan(future::sequential)
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

  studies <- if (store_studies) lapply(setting_results, `[[`, "studies") else NULL
  summary_table <- do.call(rbind, lapply(setting_results, `[[`, "summary_row"))
  rownames(summary_table) <- NULL

  out <- list(
    user_input = as.list(match.call()),
    IC = IC,
    grid = grid,
    summary_table = summary_table,
    studies = studies,
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
#'   proportional, and correlation summaries.
#' @param primary_metric Character scalar indicating the recovery metric used to
#'   rank feasible settings. Supported values are `"fuzzy_f1"`,
#'   `"fuzzy_recall"`, `"strict_f1"`, and `"weighted_fuzzy_f1"`.
#' @param scenario_weights Numeric length-2 vector giving the weights applied to
#'   the proportional and correlation scenarios when forming the ranking score.
#'   If unnamed, the first value is used for proportional and the second for
#'   correlation.
#' @param tie_break Character scalar. `"conservative"` prefers larger thresholds
#'   and larger `min_descendant_tips` when the primary ranking score is tied;
#'   `"liberal"` prefers the smaller values.
#'
#' @details
#' The selector uses a two-step decision rule. First, it filters out settings
#' that violate the null false-positive or evaluability constraints. Second, it
#' ranks the remaining settings using a weighted average of the chosen recovery
#' metric across the proportional and correlation scenarios.
#'
#' If no settings pass the filters, the function falls back to ranking the full
#' grid and records that no feasible settings were available under the supplied
#' constraints. This makes it easier to diagnose over-strict tuning criteria
#' without discarding the grid output entirely.
#'
#' @return A list of class `bifrost_search_tuning_selection` with components:
#' \describe{
#'   \item{`selected_row`}{The chosen row from `tuning_grid$summary_table`,
#'   augmented with a ranking score.}
#'   \item{`recommended_search_options`}{A search-options list ready to combine
#'   with empirical calls to [searchOptimalConfiguration()].}
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
      anyNA(scenario_weights) || any(scenario_weights < 0) ||
      sum(scenario_weights) <= 0) {
    stop("scenario_weights must be a non-negative numeric vector of length 2.")
  }

  metric_suffix <- switch(
    primary_metric,
    fuzzy_f1 = "fuzzy_f1",
    fuzzy_recall = "fuzzy_recall",
    strict_f1 = "strict_f1",
    weighted_fuzzy_f1 = "weighted_fuzzy_f1"
  )

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
  feasible <- summary_table[
    !is.na(summary_table$null_mean_false_positive_rate) &
      summary_table$null_mean_false_positive_rate <= max_false_positive_rate &
      summary_table$null_fraction_any_false_positive <= max_any_false_positive &
      summary_table$null_evaluable_fraction >= min_evaluable_fraction &
      summary_table$proportional_evaluable_fraction >= min_evaluable_fraction &
      summary_table$correlation_evaluable_fraction >= min_evaluable_fraction,
    ,
    drop = FALSE
  ]

  used_all_settings <- FALSE
  if (nrow(feasible) == 0L) {
    feasible <- summary_table
    used_all_settings <- TRUE
  }

  feasible$score <- (
    weights_named["proportional"] * feasible[[proportional_metric]]
  ) + (
    weights_named["correlation"] * feasible[[correlation_metric]]
  )
  order_score <- ifelse(is.na(feasible$score), -Inf, feasible$score)
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
    n_feasible_settings = if (used_all_settings) 0L else nrow(feasible),
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
