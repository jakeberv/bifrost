.simulation_template_fit_settings <- function(template) {
  method_setting <- template$fit_method
  if (is.null(method_setting) && !is.null(template$global_model$call$method)) {
    method_setting <- as.character(template$global_model$call$method)
  }

  error_setting <- template$fit_error
  if (is.null(error_setting) && !is.null(template$global_model$call$error)) {
    error_setting <- tryCatch(
      eval(template$global_model$call$error),
      error = function(e) NULL
    )
  }

  list(method = method_setting, error = error_setting)
}

#' Print method for bifrost simulation templates
#'
#' Print a compact summary of a `bifrost_simulation_template`, including the
#' aligned tree size, response/predictor structure, global calibration formula,
#' downstream simulation-search formula, and the empirical covariance summaries
#' used to generate simulation replicates.
#'
#' @param x A `bifrost_simulation_template` object returned by
#'   [createSimulationTemplate()].
#' @param ... Unused (S3 compatibility).
#'
#' @return Invisibly returns `x`. Called for its printing side effects.
#'
#' @seealso [createSimulationTemplate()]
#'
#' @examples
#' \dontrun{
#' set.seed(1)
#' tr <- ape::rtree(20)
#' X <- matrix(rnorm(20 * 3), ncol = 3)
#' rownames(X) <- tr$tip.label
#' tmpl <- createSimulationTemplate(tr, X, formula = "trait_data ~ 1", method = "LL")
#' tmpl
#' }
#'
#' @export
print.bifrost_simulation_template <- function(x, ...) {
  cat("Bifrost Simulation Template\n")
  cat("===========================\n\n")
  cat("Tree\n")
  cat("  Tips: ", x$n_tips, "\n", sep = "")
  cat("  Response traits: ", x$n_response_traits, "\n", sep = "")
  cat("  Predictor columns: ", length(x$predictor_columns), "\n\n", sep = "")

  cat("Model\n")
  cat("  Calibration formula: ", x$formula, "\n", sep = "")
  if (!is.null(x$search_formula)) {
    cat("  Simulation search:   ", x$search_formula, "\n", sep = "")
  }
  fit_settings <- .simulation_template_fit_settings(x)
  if (!is.null(fit_settings$method)) {
    cat("  Method: ", as.character(fit_settings$method), "\n", sep = "")
  }
  if (!is.null(fit_settings$error)) {
    cat("  Error term: ", as.character(fit_settings$error), "\n", sep = "")
  }
  cat("\n")

  cat("Empirical Covariance Summaries\n")
  cat("  Variance mean: ", formatC(x$variance_mean, format = "f", digits = 6), "\n", sep = "")
  cat("  Variance sd:   ", formatC(x$variance_sd, format = "f", digits = 6), "\n", sep = "")
  cat("  Covariance mean: ", formatC(x$covariance_mean, format = "f", digits = 6), "\n", sep = "")
  cat("  Covariance sd:   ", formatC(x$covariance_sd, format = "f", digits = 6), "\n", sep = "")

  invisible(x)
}

#' Print method for bifrost simulation studies
#'
#' Print a compact summary of a `bifrost_simulation_study`, including the study
#' type, generating scenario, replicate counts, and the main false-positive or
#' recovery summaries attached to the study object.
#'
#' @param x A `bifrost_simulation_study` object returned by
#'   [runFalsePositiveSimulationStudy()] or
#'   [runShiftRecoverySimulationStudy()].
#' @param ... Unused (S3 compatibility).
#'
#' @return Invisibly returns `x`. Called for its printing side effects.
#'
#' @seealso [runFalsePositiveSimulationStudy()],
#'   [runShiftRecoverySimulationStudy()]
#'
#' @examples
#' \dontrun{
#' # Usually called on the output of a simulation-study wrapper:
#' # study
#' }
#'
#' @export
print.bifrost_simulation_study <- function(x, ...) {
  cat("Bifrost Simulation Study\n")
  cat("=======================\n\n")
  cat("Study\n")
  cat("  Type: ", x$study_type, "\n", sep = "")
  scenario <- x$generating_scenario
  if (length(scenario) > 1L) {
    scenario <- paste(scenario, collapse = ", ")
  }
  cat("  Scenario: ", scenario, "\n", sep = "")
  cat("  Replicates: ", x$study_summary$n_replicates, "\n", sep = "")
  cat("  Completed: ", x$study_summary$n_completed, "\n", sep = "")
  cat("  Failed: ", x$study_summary$n_failed, "\n\n", sep = "")

  if (identical(x$study_type, "false_positive")) {
    cat("False-Positive Summary\n")
    cat(
      "  Mean FP rate:   ",
      formatC(x$study_summary$mean_false_positive_rate, format = "f", digits = 4),
      "\n",
      sep = ""
    )
    cat(
      "  Median FP rate: ",
      formatC(x$study_summary$median_false_positive_rate, format = "f", digits = 4),
      "\n",
      sep = ""
    )
  } else if (identical(x$study_type, "shift_recovery") && !is.null(x$evaluation)) {
    cat("Recovery Summary\n")
    cat(
      "  Strict precision: ",
      formatC(x$evaluation$strict$precision, format = "f", digits = 3),
      "\n",
      sep = ""
    )
    cat(
      "  Strict recall:    ",
      formatC(x$evaluation$strict$recall, format = "f", digits = 3),
      "\n",
      sep = ""
    )
    cat(
      "  Strict F1:        ",
      formatC(x$evaluation$strict$f1, format = "f", digits = 3),
      "\n",
      sep = ""
    )
    cat(
      "  Strict specificity: ",
      formatC(x$evaluation$strict$specificity, format = "f", digits = 3),
      "\n",
      sep = ""
    )
    cat(
      "  Strict FPR:       ",
      formatC(x$evaluation$strict$fpr, format = "f", digits = 3),
      "\n",
      sep = ""
    )
    cat(
      "  Strict balanced accuracy: ",
      formatC(x$evaluation$strict$balanced_accuracy, format = "f", digits = 3),
      "\n",
      sep = ""
    )
    cat(
      "  Fuzzy precision:  ",
      formatC(x$evaluation$fuzzy$precision, format = "f", digits = 3),
      "\n",
      sep = ""
    )
    cat(
      "  Fuzzy recall:     ",
      formatC(x$evaluation$fuzzy$recall, format = "f", digits = 3),
      "\n",
      sep = ""
    )
    cat(
      "  Fuzzy F1:         ",
      formatC(x$evaluation$fuzzy$f1, format = "f", digits = 3),
      "\n",
      sep = ""
    )
    cat(
      "  Fuzzy specificity: ",
      formatC(x$evaluation$fuzzy$specificity, format = "f", digits = 3),
      "\n",
      sep = ""
    )
    cat(
      "  Fuzzy FPR:        ",
      formatC(x$evaluation$fuzzy$fpr, format = "f", digits = 3),
      "\n",
      sep = ""
    )
    cat(
      "  Fuzzy balanced accuracy: ",
      formatC(x$evaluation$fuzzy$balanced_accuracy, format = "f", digits = 3),
      "\n",
      sep = ""
    )
  }

  invisible(x)
}

#' Convert a bifrost simulation study to a data frame
#'
#' @param x A `bifrost_simulation_study` object.
#' @param row.names,optional Arguments required by the S3 generic; ignored.
#' @param component Which table to return: `"per_replicate"` for the main
#'   replicate table or `"summary"` for a one-row study summary.
#' @param ... Unused.
#'
#' @return A data frame for the requested study component.
#'
#' @export
as.data.frame.bifrost_simulation_study <- function(x,
                                                   row.names = NULL,
                                                   optional = FALSE,
                                                   component = c("per_replicate", "summary"),
                                                   ...) {
  .simulation_check_unused_dots(list(...))
  component <- match.arg(component)

  if (identical(component, "per_replicate")) {
    out <- x$per_replicate
    rownames(out) <- NULL
    return(out)
  }

  scenario <- x$generating_scenario
  if (length(scenario) > 1L) {
    scenario <- paste(scenario, collapse = ", ")
  }

  if (identical(x$study_type, "false_positive")) {
    return(data.frame(
      study_type = x$study_type,
      generating_scenario = scenario,
      n_replicates = x$study_summary$n_replicates,
      n_completed = x$study_summary$n_completed,
      n_failed = x$study_summary$n_failed,
      n_evaluable_replicates = x$study_summary$n_evaluable_replicates,
      mean_false_positive_rate = x$study_summary$mean_false_positive_rate,
      median_false_positive_rate = x$study_summary$median_false_positive_rate,
      row.names = NULL,
      check.names = FALSE
    ))
  }

  recovery_summary_metric <- function(match, metric) {
    value <- x$study_summary[[match]][[metric]]
    if (is.null(value)) NA_real_ else unname(as.numeric(value))
  }

  data.frame(
    study_type = x$study_type,
    generating_scenario = scenario,
    n_replicates = x$study_summary$n_replicates,
    n_completed = x$study_summary$n_completed,
    n_failed = x$study_summary$n_failed,
    strict_precision = recovery_summary_metric("strict", "precision"),
    strict_recall = recovery_summary_metric("strict", "recall"),
    strict_f1 = recovery_summary_metric("strict", "f1"),
    strict_specificity = recovery_summary_metric("strict", "specificity"),
    strict_fpr = recovery_summary_metric("strict", "fpr"),
    strict_balanced_accuracy = recovery_summary_metric("strict", "balanced_accuracy"),
    fuzzy_precision = recovery_summary_metric("fuzzy", "precision"),
    fuzzy_recall = recovery_summary_metric("fuzzy", "recall"),
    fuzzy_f1 = recovery_summary_metric("fuzzy", "f1"),
    fuzzy_specificity = recovery_summary_metric("fuzzy", "specificity"),
    fuzzy_fpr = recovery_summary_metric("fuzzy", "fpr"),
    fuzzy_balanced_accuracy = recovery_summary_metric("fuzzy", "balanced_accuracy"),
    row.names = NULL,
    check.names = FALSE
  )
}

#' Print method for fixed-IC search tuning grids
#'
#' @param x A `bifrost_search_tuning_grid` object returned by
#'   [runSearchTuningGrid()].
#' @param ... Unused (S3 compatibility).
#'
#' @return Invisibly returns `x`. Called for its printing side effects.
#'
#' @seealso [runSearchTuningGrid()], [selectTunedSearchParameters()]
#'
#' @export
print.bifrost_search_tuning_grid <- function(x, ...) {
  cat("Bifrost Search Tuning Grid\n")
  cat("==========================\n\n")
  cat("Grid\n")
  cat("  IC: ", x$IC, "\n", sep = "")
  cat("  Settings: ", nrow(x$summary_table), "\n", sep = "")
  cat("  Null replicates: ", x$null_replicates, "\n", sep = "")
  cat("  Recovery replicates: ", x$recovery_replicates, "\n", sep = "")
  cat("  Weighted summaries: ", isTRUE(x$weighted), "\n", sep = "")
  cat("  Stored studies: ", isTRUE(x$store_studies), "\n\n", sep = "")

  if (!is.null(x$summary_table) && nrow(x$summary_table) > 0L) {
    cat("Ranges\n")
    cat(
      "  shift_acceptance_threshold: ",
      paste(range(x$summary_table$shift_acceptance_threshold, na.rm = TRUE), collapse = " to "),
      "\n",
      sep = ""
    )
    cat(
      "  min_descendant_tips: ",
      paste(range(x$summary_table$min_descendant_tips, na.rm = TRUE), collapse = " to "),
      "\n",
      sep = ""
    )
  }

  invisible(x)
}

#' Convert a fixed-IC search tuning grid to a data frame
#'
#' @param x A `bifrost_search_tuning_grid` object.
#' @param row.names,optional Arguments required by the S3 generic; ignored.
#' @param ... Unused.
#'
#' @return The `summary_table` component as a data frame.
#'
#' @export
as.data.frame.bifrost_search_tuning_grid <- function(x,
                                                     row.names = NULL,
                                                     optional = FALSE,
                                                     ...) {
  .simulation_check_unused_dots(list(...))
  out <- x$summary_table
  rownames(out) <- NULL
  out
}

#' Print method for tuned search-parameter selections
#'
#' @param x A `bifrost_search_tuning_selection` object returned by
#'   [selectTunedSearchParameters()].
#' @param ... Unused (S3 compatibility).
#'
#' @return Invisibly returns `x`. Called for its printing side effects.
#'
#' @seealso [selectTunedSearchParameters()], [runSearchTuningGrid()]
#'
#' @export
print.bifrost_search_tuning_selection <- function(x, ...) {
  selected <- x$selected_row
  cat("Bifrost Tuned Search Selection\n")
  cat("==============================\n\n")
  cat("Selection\n")
  cat("  IC: ", x$IC, "\n", sep = "")
  cat("  Primary metric: ", x$primary_metric, "\n", sep = "")
  cat("  Feasible settings: ", x$n_feasible_settings, "\n", sep = "")
  cat("  Used all settings: ", isTRUE(x$used_all_settings), "\n\n", sep = "")
  cat("Recommended Search Options\n")
  cat("  shift_acceptance_threshold: ", selected$shift_acceptance_threshold[[1L]], "\n", sep = "")
  cat("  min_descendant_tips: ", selected$min_descendant_tips[[1L]], "\n", sep = "")
  if ("score" %in% names(selected)) {
    cat("  score: ", formatC(selected$score[[1L]], format = "f", digits = 3), "\n", sep = "")
  }

  invisible(x)
}

#' Convert a tuned search-parameter selection to a data frame
#'
#' @param x A `bifrost_search_tuning_selection` object.
#' @param row.names,optional Arguments required by the S3 generic; ignored.
#' @param component Which table to return: `"selected"` for the chosen row or
#'   `"feasible"` for the filtered/ranked table.
#' @param ... Unused.
#'
#' @return A data frame for the requested selection component.
#'
#' @export
as.data.frame.bifrost_search_tuning_selection <- function(x,
                                                          row.names = NULL,
                                                          optional = FALSE,
                                                          component = c("selected", "feasible"),
                                                          ...) {
  .simulation_check_unused_dots(list(...))
  component <- match.arg(component)
  out <- if (identical(component, "selected")) {
    x$selected_row
  } else {
    x$feasible_table
  }
  rownames(out) <- NULL
  out
}

#' Print method for shift-recovery evaluations
#'
#' @param x A `bifrost_shift_recovery_evaluation` object returned by
#'   [evaluateShiftRecovery()].
#' @param ... Unused (S3 compatibility).
#'
#' @return Invisibly returns `x`. Called for its printing side effects.
#'
#' @seealso [evaluateShiftRecovery()]
#'
#' @export
print.bifrost_shift_recovery_evaluation <- function(x, ...) {
  fmt <- function(value, digits = 3L) {
    if (is.null(value) || length(value) == 0L || is.na(value)) {
      return("NA")
    }
    formatC(value, format = "f", digits = digits)
  }

  cat("Bifrost Shift-Recovery Evaluation\n")
  cat("=================================\n\n")
  cat("Strict\n")
  cat("  Precision: ", fmt(x$strict$precision), "\n", sep = "")
  cat("  Recall:    ", fmt(x$strict$recall), "\n", sep = "")
  cat("  F1:        ", fmt(x$strict$f1), "\n\n", sep = "")
  cat("Fuzzy\n")
  cat("  Precision: ", fmt(x$fuzzy$precision), "\n", sep = "")
  cat("  Recall:    ", fmt(x$fuzzy$recall), "\n", sep = "")
  cat("  F1:        ", fmt(x$fuzzy$f1), "\n", sep = "")

  if (!is.null(x$weighted)) {
    cat("\nWeighted Fuzzy\n")
    cat("  Precision: ", fmt(x$weighted$fuzzy$precision), "\n", sep = "")
    cat("  Recall:    ", fmt(x$weighted$fuzzy$recall), "\n", sep = "")
    cat("  F1:        ", fmt(x$weighted$fuzzy$f1), "\n", sep = "")
  }

  invisible(x)
}

#' Convert a shift-recovery evaluation to a data frame
#'
#' @param x A `bifrost_shift_recovery_evaluation` object.
#' @param row.names,optional Arguments required by the S3 generic; ignored.
#' @param component Which table to return: `"metrics"` for performance metrics
#'   or `"counts"` for aggregated contingency counts.
#' @param ... Unused.
#'
#' @return A data frame for the requested evaluation component.
#'
#' @export
as.data.frame.bifrost_shift_recovery_evaluation <- function(x,
                                                            row.names = NULL,
                                                            optional = FALSE,
                                                            component = c("metrics", "counts"),
                                                            ...) {
  .simulation_check_unused_dots(list(...))
  component <- match.arg(component)

  if (identical(component, "counts")) {
    out <- data.frame(
      match = c("strict", "fuzzy"),
      rbind(x$counts$strict, x$counts$fuzzy),
      row.names = NULL,
      check.names = FALSE
    )
    rownames(out) <- NULL
    return(out)
  }

  metric_row <- function(label, values) {
    data.frame(
      match = label,
      precision = values$precision,
      recall = values$recall,
      f1 = values$f1,
      specificity = if (is.null(values$specificity)) NA_real_ else values$specificity,
      fpr = if (is.null(values$fpr)) NA_real_ else values$fpr,
      balanced_accuracy = if (is.null(values$balanced_accuracy)) NA_real_ else values$balanced_accuracy,
      row.names = NULL,
      check.names = FALSE
    )
  }

  rows <- list(
    metric_row("strict", x$strict),
    metric_row("fuzzy", x$fuzzy)
  )
  if (!is.null(x$weighted)) {
    rows <- c(
      rows,
      list(
        metric_row("weighted_strict", x$weighted$strict),
        metric_row("weighted_fuzzy", x$weighted$fuzzy)
      )
    )
  }

  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}
