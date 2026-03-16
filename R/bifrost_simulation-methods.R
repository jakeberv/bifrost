#' Print method for bifrost simulation templates
#'
#' @description
#' Print a compact summary of a `bifrost_simulation_template`, including the
#' aligned tree size, response/predictor structure, fitted formula, and the
#' empirical covariance summaries used to generate simulation replicates.
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
  cat("  Formula: ", x$formula, "\n", sep = "")
  method_setting <- x$fit_method
  if (is.null(method_setting) && !is.null(x$global_model$call$method)) {
    method_setting <- as.character(x$global_model$call$method)
  }
  if (!is.null(method_setting)) {
    cat("  Method: ", as.character(method_setting), "\n", sep = "")
  }
  error_setting <- x$fit_error
  if (is.null(error_setting) && !is.null(x$global_model$call$error)) {
    error_setting <- tryCatch(
      eval(x$global_model$call$error),
      error = function(e) NULL
    )
  }
  if (!is.null(error_setting)) {
    cat("  Error term: ", as.character(error_setting), "\n", sep = "")
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
#' @description
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
    cat("  Mean FP rate:   ",
        formatC(x$study_summary$mean_false_positive_rate, format = "f", digits = 4),
        "\n", sep = "")
    cat("  Median FP rate: ",
        formatC(x$study_summary$median_false_positive_rate, format = "f", digits = 4),
        "\n", sep = "")
  } else if (identical(x$study_type, "shift_recovery") && !is.null(x$evaluation)) {
    cat("Recovery Summary\n")
    cat("  Strict precision: ",
        formatC(x$evaluation$strict$precision, format = "f", digits = 3),
        "\n", sep = "")
    cat("  Strict recall:    ",
        formatC(x$evaluation$strict$recall, format = "f", digits = 3),
        "\n", sep = "")
    cat("  Fuzzy precision:  ",
        formatC(x$evaluation$fuzzy$precision, format = "f", digits = 3),
        "\n", sep = "")
    cat("  Fuzzy recall:     ",
        formatC(x$evaluation$fuzzy$recall, format = "f", digits = 3),
        "\n", sep = "")
  }

  invisible(x)
}
