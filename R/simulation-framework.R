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

#' Simulate a Dataset-Matched Null Replicate
#'
#' @description
#' Generate a single no-shift simulation replicate under a uniform Brownian
#' motion process using the manuscript-style covariance-generation framework
#' ported into `bifrost`. Each replicate draws a fresh positive-definite
#' covariance matrix using empirical variance and covariance summaries from a
#' [`createSimulationTemplate()`] object.
#'
#' @param template A `bifrost_simulation_template` returned by
#'   [createSimulationTemplate()].
#' @param tree_tip_count Optional integer tip count for random subtree sampling.
#'   If `NULL`, the full empirical tree from `template` is used.
#' @param seed Optional integer random seed.
#' @param preserve_predictors Retained for backward compatibility. Simulation
#'   studies generated from a template always operate on the response block
#'   alone, so predictor columns are not carried through into simulated
#'   datasets.
#' @param ... Reserved for future extensions. Currently ignored.
#'
#' @details
#' This function is the package-clean analogue of the manuscript's
#' `simulate_traits_BM1()` workflow. The core scientific behavior is preserved:
#' a new covariance matrix is generated for each replicate from the empirical
#' variance/covariance summaries in the template, rather than reusing the fitted
#' covariance matrix directly. For templates calibrated from richer global
#' models, the empirical fitted mean structure is added back to the simulated
#' residual process before the downstream intercept-only search is run on the
#' regenerated response block.
#'
#' @return A list of class `bifrost_simulation_replicate_null` containing the
#'   sampled tree, a single-regime baseline tree, the simulated response matrix,
#'   the combined `trait_data` object used for downstream model fitting, and the
#'   generated covariance matrix.
#'
#' @seealso [createSimulationTemplate()], [simulateShiftedDataset()],
#'   [runFalsePositiveSimulationStudy()]
#'
#' @examples
#' \dontrun{
#' set.seed(1)
#' tr <- ape::rtree(20)
#' X <- matrix(rnorm(20 * 3), ncol = 3)
#' rownames(X) <- tr$tip.label
#' tmpl <- createSimulationTemplate(tr, X, formula = "trait_data ~ 1", method = "LL")
#'
#' sim_null <- simulateNullDataset(tmpl, tree_tip_count = 15, seed = 2)
#' str(sim_null$trait_data)
#' }
#'
#' @export
simulateNullDataset <- function(template,
                                tree_tip_count = NULL,
                                seed = NULL,
                                preserve_predictors = TRUE,
                                ...) {
  if (!inherits(template, "bifrost_simulation_template")) {
    stop("template must be a 'bifrost_simulation_template' object.")
  }
  if (!is.null(seed)) {
    set.seed(seed)
  }
  if (!is.null(tree_tip_count)) {
    if (!is.numeric(tree_tip_count) || length(tree_tip_count) != 1L ||
        is.na(tree_tip_count) || tree_tip_count < 2L) {
      stop("tree_tip_count must be NULL or a single integer >= 2.")
    }
    tree_tip_count <- as.integer(tree_tip_count)
    if (tree_tip_count > template$n_tips) {
      stop("tree_tip_count cannot exceed the number of tips in the template.")
    }
  }
  if (!is.logical(preserve_predictors) || length(preserve_predictors) != 1L ||
      is.na(preserve_predictors)) {
    stop("preserve_predictors must be TRUE or FALSE.")
  }
  sampled_tree <- if (is.null(tree_tip_count) || tree_tip_count == template$n_tips) {
    ape::reorder.phylo(ape::as.phylo(template$baseline_tree), order = "postorder")
  } else {
    ape::keep.tip(
      ape::reorder.phylo(ape::as.phylo(template$baseline_tree), order = "postorder"),
      tip = sample(template$baseline_tree$tip.label, size = tree_tip_count)
    )
  }

  sampled_fitted <- template$fitted_values[sampled_tree$tip.label, , drop = FALSE]
  n_traits <- template$n_response_traits

  sigma <- NULL
  attempt <- 1L
  while (attempt <= 100L && is.null(sigma)) {
    variances <- abs(stats::rnorm(
      n_traits,
      mean = template$variance_mean,
      sd = template$variance_sd
    ))

    if (n_traits == 1L) {
      candidate_sigma <- matrix(variances[1], nrow = 1L, ncol = 1L)
    } else {
      lower_tri <- matrix(
        stats::rnorm(
          n_traits * n_traits,
          mean = template$covariance_mean,
          sd = template$covariance_sd
        ),
        ncol = n_traits
      )
      lower_tri[upper.tri(lower_tri)] <- 0
      candidate_sigma <- t(lower_tri) %*% lower_tri
      diag(candidate_sigma) <- variances + diag(candidate_sigma)
    }

    eigenvalues <- eigen(candidate_sigma, symmetric = TRUE, only.values = TRUE)$values
    if (all(eigenvalues > 0)) {
      sigma <- candidate_sigma
    } else {
      attempt <- attempt + 1L
    }
  }
  if (is.null(sigma)) {
    stop("Failed to generate a positive-definite covariance matrix for the null simulation.")
  }

  simulated <- mvMORPH::mvSIM(
    tree = sampled_tree,
    nsim = 1,
    model = "BM1",
    param = list(
      ntraits = n_traits,
      sigma = sigma,
      theta = rep(0, n_traits)
    )
  )

  simulated_residuals <- if (is.list(simulated)) simulated[[1]] else simulated
  simulated_residuals <- as.matrix(simulated_residuals)
  rownames(simulated_residuals) <- sampled_tree$tip.label
  colnames(simulated_residuals) <- template$response_column_names

  simulated_response <- sampled_fitted + simulated_residuals
  colnames(simulated_response) <- template$response_column_names

  trait_data <- if (template$trait_data_is_matrix) {
    simulated_response
  } else {
    as.data.frame(simulated_response, stringsAsFactors = FALSE)
  }

  baseline_tree <- phytools::paintSubTree(
    ape::reorder.phylo(ape::as.phylo(sampled_tree), order = "postorder"),
    node = ape::Ntip(sampled_tree) + 1L,
    state = 0,
    anc.state = 0
  )

  result <- list(
    user_input = as.list(match.call()),
    tree = sampled_tree,
    baseline_tree = baseline_tree,
    data = trait_data,
    trait_data = trait_data,
    simulatedData = simulated_response,
    covariance_matrix = sigma,
    generating_scenario = "null"
  )

  class(result) <- c("bifrost_simulation_replicate_null", class(result))
  result
}

#' Simulate a Dataset-Matched Shifted Replicate
#'
#' @description
#' Generate a single known-shift simulation replicate using the manuscript
#' multi-shift framework ported into `bifrost`. The function supports both the
#' proportional generating-model scenario and the manuscript's non-generating
#' correlation-changing robustness scenario.
#'
#' @param template A `bifrost_simulation_template` returned by
#'   [createSimulationTemplate()].
#' @param tree_tip_count Optional integer tip count for random subtree sampling.
#'   If `NULL`, the full empirical tree from `template` is used.
#' @param num_shifts Integer number of shifts to simulate.
#' @param min_shift_tips Integer minimum size of any shifted clade.
#' @param max_shift_tips Integer maximum size of any shifted clade.
#' @param scale_mode Character string specifying the generating scenario:
#'   `"proportional"` for the manuscript generating model or `"correlation"` for
#'   the manuscript non-generating robustness scenario in which trait
#'   correlations evolve.
#' @param scale_factor_range Numeric length-2 vector giving the full range of
#'   possible shift scalars.
#' @param exclude_range Numeric length-2 vector specifying the central interval
#'   excluded from sampled shift scalars.
#' @param buffer Integer minimum node distance between simulated shifts.
#' @param seed Optional integer random seed.
#' @param preserve_predictors Retained for backward compatibility. Simulation
#'   studies generated from a template always operate on the response block
#'   alone, so predictor columns are not carried through into simulated
#'   datasets.
#' @param ... Reserved for future extensions. Currently ignored.
#'
#' @details
#' This function is a package-clean port of the manuscript's
#' `simulateAndPaintMultipleShifts()` workflow. It preserves the core behavior:
#' \itemize{
#'   \item random subtree sampling from the empirical tree,
#'   \item shift placement subject to clade-size and non-overlap constraints,
#'   \item manuscript-style buffer enforcement using node distances,
#'   \item generation of a fresh ancestral covariance matrix from empirical
#'         variance/covariance summaries for each replicate,
#'   \item construction of derived regimes by either proportional scaling or the
#'         non-generating correlation-scaling scenario, and
#'   \item simulation of multivariate BMM residuals that are added to the
#'         empirical fitted mean structure from the global calibration model.
#' }
#'
#' Under `scale_mode = "proportional"`, derived regimes are scalar multiples of
#' the ancestral covariance matrix and therefore match the generating assumptions
#' of the current `bifrost` BMM search. Under `scale_mode = "correlation"`, the
#' off-diagonal correlation structure is changed while variances are adjusted,
#' mirroring the manuscript's deliberate model-misspecification robustness test.
#' Downstream simulation studies still evaluate intercept-only shift searches on
#' the regenerated response block.
#'
#' @return A list of class `bifrost_simulation_replicate_shifted` containing the
#'   painted generating tree, the true shift nodes, the simulated response
#'   matrix, the combined `trait_data` object used for downstream model fitting,
#'   regime covariance matrices, sampled scale factors, and a scenario label.
#'
#' @seealso [simulateNullDataset()], [createSimulationTemplate()],
#'   [runShiftRecoverySimulationStudy()]
#'
#' @examples
#' \dontrun{
#' set.seed(1)
#' tr <- ape::rtree(40)
#' X <- matrix(rnorm(40 * 3), ncol = 3)
#' rownames(X) <- tr$tip.label
#' tmpl <- createSimulationTemplate(tr, X, formula = "trait_data ~ 1", method = "LL")
#'
#' sim_shift <- simulateShiftedDataset(
#'   tmpl,
#'   tree_tip_count = 25,
#'   num_shifts = 2,
#'   min_shift_tips = 3,
#'   max_shift_tips = 8,
#'   scale_mode = "proportional",
#'   seed = 3
#' )
#'
#' sim_shift$shiftNodes
#' }
#'
#' @export
simulateShiftedDataset <- function(template,
                                   tree_tip_count = NULL,
                                   num_shifts,
                                   min_shift_tips,
                                   max_shift_tips,
                                   scale_mode = c("proportional", "correlation"),
                                   scale_factor_range = c(0.1, 2.0),
                                   exclude_range = c(0.5, 1.5),
                                   buffer = 3,
                                   seed = NULL,
                                   preserve_predictors = TRUE,
                                   ...) {
  scale_mode <- match.arg(scale_mode)

  if (!inherits(template, "bifrost_simulation_template")) {
    stop("template must be a 'bifrost_simulation_template' object.")
  }
  if (!is.null(seed)) {
    set.seed(seed)
  }
  if (!is.numeric(num_shifts) || length(num_shifts) != 1L || is.na(num_shifts) ||
      num_shifts < 1L) {
    stop("num_shifts must be a single integer >= 1.")
  }
  if (!is.numeric(min_shift_tips) || length(min_shift_tips) != 1L ||
      is.na(min_shift_tips) || min_shift_tips < 1L) {
    stop("min_shift_tips must be a single integer >= 1.")
  }
  if (!is.numeric(max_shift_tips) || length(max_shift_tips) != 1L ||
      is.na(max_shift_tips) || max_shift_tips < min_shift_tips) {
    stop("max_shift_tips must be >= min_shift_tips.")
  }
  if (!is.numeric(buffer) || length(buffer) != 1L || is.na(buffer) || buffer < 0) {
    stop("buffer must be a single non-negative number.")
  }
  if (!is.numeric(scale_factor_range) || length(scale_factor_range) != 2L ||
      anyNA(scale_factor_range) || diff(scale_factor_range) <= 0) {
    stop("scale_factor_range must be a numeric length-2 vector in increasing order.")
  }
  if (!is.numeric(exclude_range) || length(exclude_range) != 2L ||
      anyNA(exclude_range) || diff(exclude_range) <= 0) {
    stop("exclude_range must be a numeric length-2 vector in increasing order.")
  }
  if (exclude_range[1] <= scale_factor_range[1] || exclude_range[2] >= scale_factor_range[2]) {
    stop("exclude_range must lie strictly inside scale_factor_range.")
  }
  if (scale_mode == "correlation" && template$n_response_traits < 2L) {
    stop("scale_mode = 'correlation' requires at least two response traits.")
  }
  if (!is.null(tree_tip_count)) {
    if (!is.numeric(tree_tip_count) || length(tree_tip_count) != 1L ||
        is.na(tree_tip_count) || tree_tip_count < 2L) {
      stop("tree_tip_count must be NULL or a single integer >= 2.")
    }
    tree_tip_count <- as.integer(tree_tip_count)
    if (tree_tip_count > template$n_tips) {
      stop("tree_tip_count cannot exceed the number of tips in the template.")
    }
  }
  if (!is.logical(preserve_predictors) || length(preserve_predictors) != 1L ||
      is.na(preserve_predictors)) {
    stop("preserve_predictors must be TRUE or FALSE.")
  }

  sampled_tree <- if (is.null(tree_tip_count) || tree_tip_count == template$n_tips) {
    ape::reorder.phylo(ape::as.phylo(template$baseline_tree), order = "postorder")
  } else {
    ape::keep.tip(
      ape::reorder.phylo(ape::as.phylo(template$baseline_tree), order = "postorder"),
      tip = sample(template$baseline_tree$tip.label, size = tree_tip_count)
    )
  }

  sampled_fitted <- template$fitted_values[sampled_tree$tip.label, , drop = FALSE]
  n_traits <- template$n_response_traits

  result <- NULL
  outer_attempt <- 0L

  while (outer_attempt < 100L && is.null(result)) {
    outer_attempt <- outer_attempt + 1L
    generating_tree <- phytools::paintSubTree(
      ape::reorder.phylo(ape::as.phylo(sampled_tree), order = "postorder"),
      node = ape::Ntip(sampled_tree) + 1L,
      state = "ancestral",
      anc.state = "ancestral"
    )

    root_node <- ape::Ntip(generating_tree) + 1L
    valid_candidates <- setdiff(
      unique(generating_tree$edge[, 1]),
      c(seq_len(ape::Ntip(generating_tree)), root_node)
    )
    valid_candidates <- Filter(function(node) {
      descendants <- getDescendants(generating_tree, node)
      n_desc_tips <- sum(descendants <= ape::Ntip(generating_tree))
      n_desc_tips >= min_shift_tips && n_desc_tips <= max_shift_tips
    }, valid_candidates)

    if (length(valid_candidates) < num_shifts) {
      next
    }

    shift_nodes <- integer(0)
    remaining_candidates <- valid_candidates
    placement_failed <- FALSE

    for (shift_index in seq_len(num_shifts)) {
      if (length(remaining_candidates) == 0L) {
        placement_failed <- TRUE
        break
      }

      sampled_node <- if (length(remaining_candidates) == 1L) {
        remaining_candidates[[1L]]
      } else {
        remaining_candidates[[sample.int(length(remaining_candidates), 1L)]]
      }
      shift_nodes <- c(shift_nodes, sampled_node)

      descendants <- getDescendants(generating_tree, sampled_node)
      ancestors <- integer(0)
      current_parent <- generating_tree$edge[generating_tree$edge[, 2] == sampled_node, 1]
      while (length(current_parent) == 1L &&
             !is.na(current_parent) &&
             current_parent >= (ape::Ntip(generating_tree) + 1L)) {
        ancestors <- c(ancestors, current_parent)
        current_parent <- generating_tree$edge[generating_tree$edge[, 2] == current_parent, 1]
      }
      remaining_candidates <- setdiff(
        remaining_candidates,
        c(descendants, ancestors, sampled_node)
      )
    }

    if (placement_failed) {
      next
    }

    if (length(shift_nodes) > 1L) {
      buffer_ok <- TRUE
      for (i in seq_len(length(shift_nodes) - 1L)) {
        for (j in seq.int(i + 1L, length(shift_nodes))) {
          distance_info <- RRphylo::distNodes(
            ape::reorder.phylo(ape::as.phylo(generating_tree), order = "cladewise"),
            node = c(shift_nodes[i], shift_nodes[j]),
            clus = 0
          )
          if (distance_info$node < buffer) {
            buffer_ok <- FALSE
            break
          }
        }
        if (!buffer_ok) {
          break
        }
      }
      if (!buffer_ok) {
        next
      }
    }

    ancestral_sigma <- NULL
    sigma_attempt <- 1L
    while (sigma_attempt <= 100L && is.null(ancestral_sigma)) {
      sampled_variances <- abs(stats::rnorm(
        n_traits,
        mean = template$variance_mean,
        sd = template$variance_sd
      ))

      if (n_traits == 1L) {
        candidate_sigma <- matrix(sampled_variances[1], nrow = 1L, ncol = 1L)
      } else {
        lower_tri <- matrix(
          stats::rnorm(
            n_traits * n_traits,
            mean = template$covariance_mean,
            sd = template$covariance_sd
          ),
          ncol = n_traits
        )
        lower_tri[upper.tri(lower_tri)] <- 0
        candidate_sigma <- t(lower_tri) %*% lower_tri
        diag(candidate_sigma) <- sampled_variances + diag(candidate_sigma)
      }

      if (all(eigen(candidate_sigma, symmetric = TRUE, only.values = TRUE)$values > 0)) {
        ancestral_sigma <- candidate_sigma
      } else {
        sigma_attempt <- sigma_attempt + 1L
      }
    }
    if (is.null(ancestral_sigma)) {
      next
    }

    sigma_list <- list(ancestral = ancestral_sigma)
    vcvs <- list(ancestral = ancestral_sigma)
    sampled_scale_factors <- numeric(num_shifts)
    tree_with_shifts <- generating_tree
    derived_failed <- FALSE

    for (shift_count in seq_len(num_shifts)) {
      derived_state <- paste("derived", shift_count, sep = "_")
      part1 <- stats::runif(1, min = scale_factor_range[1], max = exclude_range[1])
      part2 <- stats::runif(1, min = exclude_range[2], max = scale_factor_range[2])
      current_scale_factor <- sample(c(part1, part2), 1L)
      sampled_scale_factors[shift_count] <- current_scale_factor

      derived_sigma <- NULL
      derived_attempt <- 1L
      while (derived_attempt <= 100L && is.null(derived_sigma)) {
        if (scale_mode == "correlation") {
          correlation_matrix <- stats::cov2cor(ancestral_sigma)
          scaled_cor <- correlation_matrix
          scaled_cor[lower.tri(scaled_cor)] <-
            correlation_matrix[lower.tri(correlation_matrix)] * current_scale_factor
          scaled_cor[upper.tri(scaled_cor)] <- t(scaled_cor)[upper.tri(scaled_cor)]
          sd_vector <- sqrt(diag(ancestral_sigma))
          new_variances <- (sd_vector^2) / current_scale_factor
          candidate_sigma <- diag(sqrt(new_variances)) %*%
            scaled_cor %*%
            diag(sqrt(new_variances))
        } else {
          candidate_sigma <- current_scale_factor * ancestral_sigma
        }

        if (all(eigen(candidate_sigma, symmetric = TRUE, only.values = TRUE)$values > 0)) {
          derived_sigma <- candidate_sigma
        } else {
          derived_attempt <- derived_attempt + 1L
        }
      }

      if (is.null(derived_sigma)) {
        derived_failed <- TRUE
        break
      }

      sigma_list[[derived_state]] <- derived_sigma
      vcvs[[derived_state]] <- derived_sigma
      tree_with_shifts <- phytools::paintSubTree(
        tree_with_shifts,
        node = shift_nodes[shift_count],
        state = derived_state,
        anc.state = "ancestral"
      )
    }

    if (derived_failed) {
      next
    }
    if (length(unique(phytools::getStates(tree_with_shifts))) != (num_shifts + 1L)) {
      next
    }

    simulated <- mvMORPH::mvSIM(
      tree = tree_with_shifts,
      nsim = 1,
      model = "BMM",
      param = list(
        ntraits = n_traits,
        sigma = sigma_list,
        theta = rep(0, n_traits)
      )
    )
    simulated_residuals <- if (is.list(simulated)) simulated[[1]] else simulated
    simulated_residuals <- as.matrix(simulated_residuals)
    rownames(simulated_residuals) <- tree_with_shifts$tip.label
    colnames(simulated_residuals) <- template$response_column_names

    simulated_response <- sampled_fitted + simulated_residuals
    colnames(simulated_response) <- template$response_column_names

    trait_data <- if (template$trait_data_is_matrix) {
      simulated_response
    } else {
      as.data.frame(simulated_response, stringsAsFactors = FALSE)
    }

    baseline_tree <- phytools::paintSubTree(
      ape::reorder.phylo(ape::as.phylo(sampled_tree), order = "postorder"),
      node = ape::Ntip(sampled_tree) + 1L,
      state = 0,
      anc.state = 0
    )

    result <- list(
      user_input = as.list(match.call()),
      sampled_tree = sampled_tree,
      baseline_tree = baseline_tree,
      paintedTree = tree_with_shifts,
      shiftNodes = shift_nodes,
      simulatedData = simulated_response,
      trait_data = trait_data,
      VCVs = vcvs,
      sampledScaleFactors = sampled_scale_factors,
      generating_scenario = scale_mode
    )
  }

  if (is.null(result)) {
    stop("Failed to place all requested shifts after 100 attempts.")
  }

  class(result) <- c("bifrost_simulation_replicate_shifted", class(result))
  result
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
  search_opts$formula <- validateSimulationStudyFormula(search_opts$formula)
  if (isTRUE(num_cores > 1L) && isTRUE(search_opts$num_cores > 1L)) {
    warning("Both wrapper-level and search-level parallelism are > 1; nested parallelism may be inefficient.")
  }

  simdata <- lapply(seq_len(n_replicates), function(i) {
    sim_args <- utils::modifyList(
      simulation_options,
      list(
        template = template,
        tree_tip_count = tree_tip_count
      )
    )
    do.call(simulateNullDataset, sim_args)
  })

  old_plan <- future::plan()
  on.exit(future::plan(old_plan), add = TRUE)
  if (num_cores > 1L) {
    future::plan(future::multisession, workers = num_cores)
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
        do.call(searchOptimalConfiguration, search_args),
        error = function(e) {
          candidate_count <- max(length(generatePaintedTrees(
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
  search_opts$formula <- validateSimulationStudyFormula(search_opts$formula)
  if (isTRUE(num_cores > 1L) && isTRUE(search_opts$num_cores > 1L)) {
    warning("Both wrapper-level and search-level parallelism are > 1; nested parallelism may be inefficient.")
  }

  old_plan <- future::plan()
  on.exit(future::plan(old_plan), add = TRUE)
  if (num_cores > 1L) {
    future::plan(future::multisession, workers = num_cores)
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
      out <- do.call(simulateShiftedDataset, sim_args)
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
        do.call(searchOptimalConfiguration, search_args),
        error = function(e) {
          candidate_count <- max(length(generatePaintedTrees(
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
