#' Create a Dataset-Matched Simulation Template for `bifrost`
#'
#' @description
#' Fit a global single-regime multivariate GLS model to an empirical dataset and
#' convert the fitted object into a reusable simulation template. This function is
#' the package-native entry point for the manuscript workflow in which simulation
#' parameters were estimated from a global empirical fit before running null and
#' shift-recovery studies.
#'
#' @param baseline_tree A rooted phylogenetic tree of class `phylo` (or coercible
#'   via [ape::as.phylo()]) whose tip labels match the row names of `trait_data`.
#' @param trait_data A numeric matrix or data frame containing the empirical data
#'   used to parameterize the simulations. Row names must match
#'   `baseline_tree$tip.label`.
#' @param formula Formula specification passed to [mvMORPH::mvgls()]. May be a
#'   single character string or a formula object. Legacy forms that reference
#'   `trait_data` directly remain supported (for example, `"trait_data ~ 1"` or
#'   `"trait_data[, 1:12] ~ trait_data[, 13]"`), but named-column formulas such
#'   as `cbind(y1, y2) ~ size + grp` are also accepted.
#' @param response_columns Optional column specification identifying the
#'   multivariate response. May be integer positions, character column names, or
#'   a mix resolvable against `trait_data`. For intercept-only workflows this
#'   defaults to all columns. For any workflow where the response is a subset of
#'   `trait_data`, including intercept-only formulas such as
#'   `"trait_data[, 1:12] ~ 1"`, this should be supplied explicitly.
#' @param predictor_columns Optional column specification identifying predictor
#'   columns used by the global calibration model. For named-column formulas
#'   these are usually inferred from the formula itself. For non-intercept
#'   workflows where `predictor_columns` is omitted, the complement of
#'   `response_columns` is used unless the formula identifies raw predictor
#'   columns directly.
#' @param ... Additional arguments passed to [mvMORPH::mvgls()] for the global
#'   empirical fit (for example `method = "LL"` or `error = TRUE`).
#'
#' @details
#' This function ports the manuscript parameterization step into the package as
#' closely as possible while cleaning it up for package use. The returned template
#' stores:
#' \itemize{
#'   \item the aligned empirical tree and calibration data,
#'   \item the fitted global `mvgls` model,
#'   \item evaluated copies of key fit settings (`fit_method`, `fit_error`) that
#'         can be reused safely by downstream simulation-study wrappers,
#'   \item the fitted response mean structure,
#'   \item the covariance matrix used to summarize empirical variance and
#'         covariance magnitudes, and
#'   \item manuscript-style summaries of the diagonal and off-diagonal elements
#'         (`variance_mean`, `variance_sd`, `covariance_mean`, `covariance_sd`)
#'         used by downstream simulation functions to generate a fresh ancestral
#'         covariance matrix for each replicate.
#' }
#'
#' The template is a calibration object, not a full simulation-study design. Its
#' formula defines the one global mean model used to estimate fitted values and
#' residual covariance structure. Downstream simulation studies then regenerate
#' the response block around those fitted means and evaluate intercept-only
#' shift-search behavior on the simulated responses. In other words, a richer
#' calibration model can still feed the manuscript-style `trait_data ~ 1`
#' simulation workflow.
#'
#' Internally, all supported formulas are normalized into a single simulation
#' formula specification. Legacy indexed formulas are rewritten onto named
#' columns, formula objects are accepted alongside character strings, predictor
#' schemas record numeric/logical/factor/ordered predictors, and unsupported
#' forms such as transformed responses, `.` shorthand, and character predictors
#' are rejected early.
#'
#' @return A list of class `bifrost_simulation_template` containing the aligned
#'   inputs, fitted global model, empirical mean structure, response/predictor
#'   column metadata, stored fit settings, and manuscript-style covariance
#'   summaries used for simulation.
#'
#' @seealso [simulateNullDataset()], [simulateShiftedDataset()],
#'   [runFalsePositiveSimulationStudy()], [runShiftRecoverySimulationStudy()],
#'   [searchOptimalConfiguration()]
#'
#' @examples
#' \dontrun{
#' set.seed(1)
#' tr <- ape::rtree(20)
#' X <- matrix(rnorm(20 * 3), ncol = 3)
#' rownames(X) <- tr$tip.label
#'
#' tmpl <- createSimulationTemplate(
#'   baseline_tree = tr,
#'   trait_data = X,
#'   formula = "trait_data ~ 1",
#'   method = "LL"
#' )
#'
#' tmpl
#' }
#'
#' @export
createSimulationTemplate <- function(baseline_tree,
                                     trait_data,
                                     formula = "trait_data ~ 1",
                                     response_columns = NULL,
                                     predictor_columns = NULL,
                                     ...) {
  extra_args <- list(...)

  if ("model" %in% names(extra_args)) {
    stop("Do not supply 'model' to createSimulationTemplate(); the global fit is always BM.")
  }
  if ("data" %in% names(extra_args)) {
    stop("Do not supply 'data' to createSimulationTemplate(); trait_data is always the empirical data source.")
  }

  if (!inherits(ape::as.phylo(baseline_tree), "phylo")) {
    stop("baseline_tree must be coercible to class 'phylo'.")
  }
  if (!is.matrix(trait_data) && !is.data.frame(trait_data)) {
    stop("trait_data must be a matrix or data.frame.")
  }
  if (is.null(rownames(trait_data))) {
    stop("trait_data must have row names matching the tree tip labels.")
  }
  formula_obj <- simulationFormulaToObject(formula)

  aligned_tree <- ape::as.phylo(baseline_tree)
  if (!all(aligned_tree$tip.label %in% rownames(trait_data))) {
    stop("All tree tip labels must be present in rownames(trait_data).")
  }
  if (!all(rownames(trait_data) %in% aligned_tree$tip.label)) {
    stop("trait_data contains row names that are not present in baseline_tree$tip.label.")
  }

  aligned_trait_data <- trait_data[aligned_tree$tip.label, , drop = FALSE]
  if (is.matrix(aligned_trait_data) && is.null(colnames(aligned_trait_data))) {
    colnames(aligned_trait_data) <- paste0("V", seq_len(ncol(aligned_trait_data)))
  }
  aligned_data_for_spec <- if (is.data.frame(aligned_trait_data)) {
    aligned_trait_data
  } else {
    as.data.frame(aligned_trait_data, check.names = FALSE, stringsAsFactors = FALSE)
  }
  formula_spec <- normalizeSimulationFormulaSpec(
    formula = formula_obj,
    trait_data = aligned_data_for_spec,
    response_columns = response_columns,
    predictor_columns = predictor_columns
  )

  response_idx <- formula_spec$response_columns
  predictor_idx <- formula_spec$predictor_columns
  intercept_only <- formula_spec$intercept_only

  if (intercept_only) {
    analysis_trait_data <- as.matrix(aligned_data_for_spec[, response_idx, drop = FALSE])
    colnames(analysis_trait_data) <- formula_spec$response_column_names
    fit_formula_obj <- formula_spec$formula_normalized_obj
    fit_call_args <- c(
      list(
        formula = fit_formula_obj,
        tree = phytools::paintSubTree(
          ape::reorder.phylo(aligned_tree, order = "postorder"),
          node = ape::Ntip(aligned_tree) + 1L,
          state = 0,
          anc.state = 0
        ),
        model = "BM",
        data = list(trait_data = analysis_trait_data)
      ),
      extra_args
    )
  } else {
    analysis_trait_data <- aligned_data_for_spec
    response_block <- analysis_trait_data[, formula_spec$response_column_names, drop = FALSE]
    if (!all(vapply(response_block, is.numeric, logical(1)))) {
      stop("Response columns in simulation templates must be numeric.")
    }
    formula_spec$data_prototype <- analysis_trait_data[0, , drop = FALSE]
    normalized_fit <- normalizeMvglsFormulaCall(
      formula = formula_spec$formula_normalized_obj,
      trait_data = analysis_trait_data,
      args_list = list(data = analysis_trait_data),
      allow_single_response = TRUE
    )
    fit_formula_obj <- normalized_fit$formula
    fit_call_args <- c(
      list(
        formula = fit_formula_obj,
        tree = phytools::paintSubTree(
          ape::reorder.phylo(aligned_tree, order = "postorder"),
          node = ape::Ntip(aligned_tree) + 1L,
          state = 0,
          anc.state = 0
        ),
        model = "BM"
      ),
      normalized_fit$args_list,
      extra_args
    )
  }

  global_model <- do.call(mvMORPH::mvgls, fit_call_args)

  sigma_source <- global_model$sigma$Pinv
  if (is.null(sigma_source)) {
    sigma_source <- global_model$sigma$S
  }
  if (is.null(sigma_source)) {
    stop("Could not locate a covariance matrix in the fitted global model.")
  }
  sigma_source <- as.matrix(sigma_source)

  variances <- diag(sigma_source)
  variance_mean <- mean(variances)
  variance_sd <- stats::sd(variances)
  if (length(variances) == 1L || is.na(variance_sd)) {
    variance_sd <- 0
  }

  if (ncol(sigma_source) > 1L) {
    covariance_values <- sigma_source[upper.tri(sigma_source)]
    covariance_mean <- mean(covariance_values)
    covariance_sd <- stats::sd(covariance_values)
    if (is.na(covariance_sd)) {
      covariance_sd <- 0
    }
  } else {
    covariance_mean <- 0
    covariance_sd <- 0
  }

  fitted_values <- as.matrix(global_model$fitted)
  rownames(fitted_values) <- rownames(analysis_trait_data)
  if (is.null(colnames(fitted_values))) {
    colnames(fitted_values) <- formula_spec$response_column_names
  }

  result <- list(
    user_input = as.list(match.call()),
    baseline_tree = aligned_tree,
    trait_data = analysis_trait_data,
    calibration_formula = formula_spec$formula_original_chr,
    formula = formula_spec$formula_original_chr,
    formula_original = formula_spec$formula_original,
    formula_normalized = formula_spec$formula_normalized_chr,
    formula_normalized_obj = formula_spec$formula_normalized_obj,
    search_formula = "trait_data ~ 1",
    formula_mode = formula_spec$formula_mode,
    response_columns = response_idx,
    response_column_names = formula_spec$response_column_names,
    predictor_columns = predictor_idx,
    predictor_column_names = formula_spec$predictor_column_names,
    predictor_schema = formula_spec$predictor_schema,
    data_prototype = if (intercept_only) NULL else formula_spec$data_prototype,
    trait_data_is_matrix = is.matrix(aligned_trait_data),
    fit_method = if ("method" %in% names(extra_args)) extra_args$method else NULL,
    fit_error = if ("error" %in% names(extra_args)) extra_args$error else NULL,
    global_model = global_model,
    fitted_values = fitted_values,
    residual_covariance = sigma_source,
    variance_mean = variance_mean,
    variance_sd = variance_sd,
    covariance_mean = covariance_mean,
    covariance_sd = covariance_sd,
    n_response_traits = ncol(fitted_values),
    n_tips = ape::Ntip(aligned_tree)
  )

  class(result) <- c("bifrost_simulation_template", class(result))
  result
}
