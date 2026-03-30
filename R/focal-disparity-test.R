.bifrost_resolve_focal_tips <- function(tree, focal) {
  if (is.character(focal)) {
    return(unique(focal))
  }

  if (is.logical(focal)) {
    if (length(focal) != ape::Ntip(tree)) {
      stop("Logical focal selectors must have length equal to ape::Ntip(tree).")
    }
    return(tree$tip.label[which(focal)])
  }

  if (is.numeric(focal) && length(focal) == 1L && !is.na(focal)) {
    node <- as.integer(focal)
    if (node <= ape::Ntip(tree)) {
      stop("Numeric focal values must be internal node numbers, not tip indices.")
    }
    if (node > (ape::Ntip(tree) + ape::Nnode(tree))) {
      stop("Numeric focal values must be internal node numbers present in tree.")
    }
    return(ape::extract.clade(tree, node = node)$tip.label)
  }

  stop("focal must be a character vector of tip labels, a logical tip selector, or a single internal node number.")
}

.bifrost_align_focal_test_inputs <- function(tree, trait_data, focal, min_focal_tips) {
  if (!inherits(tree, "phylo")) {
    tree <- ape::as.phylo(tree)
  }
  if (!inherits(tree, "phylo")) {
    stop("tree must be coercible to class 'phylo'.")
  }
  if (!is.matrix(trait_data) && !is.data.frame(trait_data)) {
    stop("trait_data must be a matrix or data.frame.")
  }
  if (is.null(rownames(trait_data))) {
    stop("trait_data must have row names matching tree tip labels.")
  }
  if (!is.numeric(min_focal_tips) || length(min_focal_tips) != 1L ||
      is.na(min_focal_tips) || min_focal_tips < 2) {
    stop("min_focal_tips must be a single integer >= 2.")
  }
  min_focal_tips <- as.integer(min_focal_tips)

  focal_tips_full <- .bifrost_resolve_focal_tips(tree, focal)
  keep <- intersect(tree$tip.label, rownames(trait_data))
  if (length(keep) < 3L) {
    stop("Fewer than 3 overlapping taxa between tree and trait_data.")
  }

  aligned_tree <- if (length(keep) == ape::Ntip(tree)) {
    tree
  } else {
    ape::drop.tip(tree, setdiff(tree$tip.label, keep))
  }

  aligned_tree <- ape::reorder.phylo(aligned_tree, order = "postorder")
  aligned_trait_data <- trait_data[aligned_tree$tip.label, , drop = FALSE]
  focal_tips <- intersect(focal_tips_full, aligned_tree$tip.label)

  if (length(focal_tips) < min_focal_tips) {
    stop(sprintf(
      "Need at least %d focal tips after alignment; found %d.",
      min_focal_tips,
      length(focal_tips)
    ))
  }

  list(
    tree = aligned_tree,
    trait_data = aligned_trait_data,
    focal_tips = focal_tips
  )
}

.bifrost_make_disparity_metric <- function(metric) {
  if (is.function(metric)) {
    return(list(fun = metric, label = "user_function"))
  }

  metric <- match.arg(
    metric,
    choices = c("centroid_msd", "trace_cov", "mean_pairwise_sqdist")
  )

  if (identical(metric, "centroid_msd")) {
    return(list(
      fun = function(X) {
        ctr <- colMeans(X)
        mean(rowSums((X - matrix(ctr, nrow(X), ncol(X), byrow = TRUE))^2))
      },
      label = metric
    ))
  }

  if (identical(metric, "trace_cov")) {
    return(list(
      fun = function(X) sum(diag(stats::cov(X))),
      label = metric
    ))
  }

  list(
    fun = function(X) {
      dmat <- as.matrix(stats::dist(X))
      mean(dmat[upper.tri(dmat)]^2)
    },
    label = metric
  )
}

.bifrost_eval_disparity_metric <- function(metric_fun, X) {
  value <- metric_fun(as.matrix(X))
  if (!is.numeric(value) || length(value) != 1L || !is.finite(value)) {
    stop("metric must return a single finite numeric value.")
  }
  as.numeric(value)
}

.bifrost_extract_sigma_matrix <- function(fit) {
  sigma_hat <- fit$sigma

  if (is.null(sigma_hat)) {
    stop("Could not locate a sigma estimate in the fitted mvgls object.")
  }
  if (is.matrix(sigma_hat)) {
    return(as.matrix(sigma_hat))
  }
  if (is.list(sigma_hat) && !is.null(sigma_hat$S)) {
    return(as.matrix(sigma_hat$S))
  }
  if (is.list(sigma_hat) && !is.null(sigma_hat$Pinv)) {
    return(as.matrix(sigma_hat$Pinv))
  }

  stop("Could not extract a usable sigma matrix from the fitted mvgls object.")
}

.bifrost_resolve_simulation_model <- function(fit_model) {
  switch(
    fit_model,
    BM = "BM1",
    OU = "OU1",
    EB = "EB",
    stop("model must be one of 'BM', 'OU', or 'EB'.")
  )
}

.bifrost_expand_ou_alpha <- function(alpha, n_traits) {
  if (is.matrix(alpha)) {
    if (!identical(dim(alpha), c(n_traits, n_traits))) {
      stop("OU null simulation requires an alpha matrix with dimensions n_traits x n_traits.")
    }
    return(alpha)
  }

  alpha <- as.numeric(alpha)
  if (length(alpha) == 1L) {
    return(diag(alpha, n_traits))
  }
  if (length(alpha) == n_traits) {
    return(diag(alpha, n_traits))
  }

  stop("OU null simulation currently requires a scalar, length-n_traits vector, or square alpha matrix.")
}

.bifrost_expand_eb_beta <- function(beta, n_traits) {
  if (is.matrix(beta)) {
    if (!identical(dim(beta), c(n_traits, n_traits))) {
      stop("EB null simulation requires a beta matrix with dimensions n_traits x n_traits.")
    }
    return(beta)
  }

  beta <- as.numeric(beta)
  if (length(beta) == 1L) {
    return(beta)
  }
  if (length(beta) == n_traits) {
    return(diag(beta, n_traits))
  }

  stop("EB null simulation currently requires a scalar, length-n_traits vector, or square beta matrix.")
}

.bifrost_build_null_simulation_params <- function(fit, simulation_model, n_traits, sigma_hat) {
  params <- list(
    ntraits = n_traits,
    sigma = sigma_hat,
    theta = rep(0, n_traits)
  )

  fit_param <- fit$param
  if (identical(simulation_model, "OU1")) {
    if (is.null(fit_param) || anyNA(fit_param)) {
      stop("Could not extract OU alpha from the fitted mvgls object.")
    }
    params$alpha <- .bifrost_expand_ou_alpha(fit_param, n_traits)
    return(params)
  }

  if (identical(simulation_model, "EB")) {
    if (is.null(fit_param) || anyNA(fit_param)) {
      stop("Could not extract EB beta from the fitted mvgls object.")
    }
    params$beta <- .bifrost_expand_eb_beta(fit_param, n_traits)
    return(params)
  }

  params
}

.bifrost_extract_simulated_matrix <- function(simulated, tip_labels, trait_names) {
  Ysim <- simulated
  if (is.list(Ysim)) {
    Ysim <- Ysim[[1L]]
  }
  if (length(dim(Ysim)) == 3L) {
    Ysim <- Ysim[, , 1L, drop = TRUE]
  }

  Ysim <- as.matrix(Ysim)
  rownames(Ysim) <- tip_labels
  colnames(Ysim) <- trait_names
  Ysim
}

.bifrost_get_response_names <- function(fitted_values, residuals_matrix) {
  fitted_names <- colnames(fitted_values)
  if (!is.null(fitted_names) && length(fitted_names) > 0L) {
    return(as.character(fitted_names))
  }

  residual_names <- colnames(residuals_matrix)
  if (!is.null(residual_names) && length(residual_names) > 0L) {
    return(as.character(residual_names))
  }

  stop("Could not determine response column names from the fitted mvgls object.")
}

.bifrost_replace_response_block <- function(data, response_names, response_matrix) {
  out <- as.data.frame(data, check.names = FALSE, stringsAsFactors = FALSE)
  response_df <- as.data.frame(response_matrix, check.names = FALSE, stringsAsFactors = FALSE)
  colnames(response_df) <- response_names
  out[, response_names] <- response_df[, response_names, drop = FALSE]
  out
}

.bifrost_build_mvgls_call_args <- function(formula,
                                           data,
                                           tree,
                                           model,
                                           method,
                                           normalized_extra_args,
                                           extra_args) {
  c(
    list(
      formula = formula,
      data = data,
      tree = tree,
      model = model,
      method = method
    ),
    normalized_extra_args,
    extra_args
  )
}

.bifrost_formula_from_string <- function(formula_chr) {
  stats::as.formula(formula_chr, env = baseenv())
}

.bifrost_call_mvgls <- function(call_args) {
  call <- as.call(c(list(quote(mvMORPH::mvgls)), unname(call_args)))
  names(call)[-1] <- names(call_args)
  eval(call, envir = parent.frame())
}

.bifrost_run_bootstrap <- function(nsim,
                                   sim_one,
                                   seed,
                                   show_progress,
                                   use_future) {
  if (isTRUE(use_future)) {
    if (isTRUE(show_progress)) {
      return(progressr::with_progress({
        progress <- progressr::progressor(steps = nsim)
        future.apply::future_sapply(
          X = seq_len(nsim),
          FUN = function(i) sim_one(i, progress = progress),
          future.seed = seed
        )
      }))
    }

    return(future.apply::future_sapply(
      X = seq_len(nsim),
      FUN = function(i) sim_one(i, progress = NULL),
      future.seed = seed
    ))
  }

  old_seed_exists <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (old_seed_exists) {
    old_seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  }
  set.seed(seed)
  on.exit({
    if (old_seed_exists) {
      assign(".Random.seed", old_seed, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)

  if (isTRUE(show_progress)) {
    return(progressr::with_progress({
      progress <- progressr::progressor(steps = nsim)
      vapply(
        X = seq_len(nsim),
        FUN = function(i) sim_one(i, progress = progress),
        FUN.VALUE = numeric(1)
      )
    }))
  }

  vapply(
    X = seq_len(nsim),
    FUN = function(i) sim_one(i, progress = NULL),
    FUN.VALUE = numeric(1)
  )
}

.bifrost_summarize_null <- function(observed, null_statistics) {
  null_mean <- mean(null_statistics)
  null_sd <- stats::sd(null_statistics)
  ses <- if (is.finite(null_sd) && null_sd > 0) {
    (observed - null_mean) / null_sd
  } else {
    NA_real_
  }

  p_upper <- (1 + sum(null_statistics >= observed)) / (length(null_statistics) + 1)
  p_lower <- (1 + sum(null_statistics <= observed)) / (length(null_statistics) + 1)

  list(
    observed = observed,
    null_mean = null_mean,
    null_sd = null_sd,
    SES = ses,
    p_upper = p_upper,
    p_lower = p_lower,
    p_two_tailed = min(1, 2 * min(p_upper, p_lower))
  )
}

.bifrost_find_example_focal_node <- function(tree, min_tips = 8L, max_tips = 15L) {
  nodes <- seq.int(ape::Ntip(tree) + 1L, ape::Ntip(tree) + ape::Nnode(tree))
  tip_counts <- vapply(
    nodes,
    function(node) length(ape::extract.clade(tree, node = node)$tip.label),
    integer(1)
  )

  eligible <- nodes[tip_counts >= min_tips & tip_counts <= max_tips]
  if (length(eligible) > 0L) {
    return(eligible[1L])
  }

  nodes[order(abs(tip_counts - round((min_tips + max_tips) / 2)))[1L]]
}

.bifrost_make_significant_focal_disparity_example <- function(seed = 2L) {
  set.seed(seed)

  tree <- ape::rtree(40)
  focal_node <- .bifrost_find_example_focal_node(tree, min_tips = 8L, max_tips = 14L)
  focal_tips <- ape::extract.clade(tree, node = focal_node)$tip.label

  n_tips <- ape::Ntip(tree)
  size <- stats::rnorm(n_tips)

  residual_block <- matrix(stats::rnorm(n_tips * 3L, sd = 0.2), ncol = 3L)
  focal_idx <- match(focal_tips, tree$tip.label)
  residual_block[focal_idx, ] <- cbind(
    stats::rnorm(length(focal_tips), sd = 0.95),
    stats::rnorm(length(focal_tips), sd = 0.85),
    stats::rnorm(length(focal_tips), sd = 0.9)
  )

  trait_data <- data.frame(
    y1 = 0.8 * size + residual_block[, 1L],
    y2 = -0.5 * size + residual_block[, 2L],
    y3 = 0.35 * size + residual_block[, 3L],
    size = size
  )
  rownames(trait_data) <- tree$tip.label

  list(
    tree = tree,
    trait_data = trait_data,
    focal_node = focal_node,
    focal_tips = focal_tips
  )
}

#' Create Reproducible Example Data for `testFocalDisparity()`
#'
#' @description
#' Generate a small synthetic example designed to yield elevated focal disparity
#' for a chosen subclade after accounting for a size covariate.
#'
#' @param seed Integer random seed controlling the tree, focal node, and data.
#'
#' @return A list with components `tree`, `trait_data`, `focal_node`, and
#'   `focal_tips`, suitable for passing directly to [testFocalDisparity()].
#'
#' @examples
#' ex <- focalDisparityExampleData(seed = 2)
#' names(ex)
#' length(ex$focal_tips)
#'
#' @export
focalDisparityExampleData <- function(seed = 2L) {
  .bifrost_make_significant_focal_disparity_example(seed = seed)
}

#' Test Focal Disparity by Parametric Bootstrap
#'
#' @description
#' Fit a multivariate phylogenetic GLS model with [mvMORPH::mvgls()], compute a
#' focal disparity statistic on the residuals of a focal set of tips, and
#' compare that observed value to a parametric bootstrap null distribution.
#'
#' @param tree A rooted phylogeny of class `phylo` or coercible with
#'   [ape::as.phylo()].
#' @param trait_data A numeric matrix or data frame with row names matching tip
#'   labels in `tree`. Extra rows or tree tips are dropped by intersection.
#' @param focal Focal definition. Supply either a character vector of tip labels,
#'   a logical selector of length `ape::Ntip(tree)`, or a single internal node
#'   number whose descendant tips define the focal clade.
#' @param formula Formula passed to [mvMORPH::mvgls()]. Named-column formulas
#'   such as `cbind(y1, y2) ~ size + group` are supported, as are response-only
#'   matrix workflows such as `"trait_data ~ 1"`.
#' @param metric Either one of `"centroid_msd"`, `"trace_cov"`, or
#'   `"mean_pairwise_sqdist"`, or a function that accepts the focal residual
#'   matrix and returns a single numeric value. The built-in options are all
#'   disparity/dispersion summaries.
#' @param statistic Deprecated alias for `metric`. Kept for backward
#'   compatibility.
#' @param nsim Integer number of bootstrap replicates.
#' @param model Character `mvgls` fit-model code used for the observed fit,
#'   bootstrap simulations, and each bootstrap refit. Supported values are
#'   `"BM"`, `"OU"`, and `"EB"`.
#' @param method Character fitting method passed to [mvMORPH::mvgls()].
#' @param workers Integer number of workers used when `future_plan` is
#'   `"multisession"` or `"multicore"`.
#' @param seed Integer seed passed to `future.apply` for reproducible bootstrap
#'   streams.
#' @param future_plan Future backend used for bootstrap refits. One of
#'   `"multisession"`, `"multicore"`, or `"sequential"`. If `NULL`, the
#'   function uses `"sequential"` when `workers = 1` and `"multisession"`
#'   otherwise.
#' @param show_progress Logical; if `TRUE`, wrap bootstrap evaluation in
#'   [progressr::with_progress()]. Users can set their preferred progress
#'   handler with `progressr::handlers()`.
#' @param min_focal_tips Minimum number of focal tips required after alignment.
#' @param ... Additional arguments passed to [mvMORPH::mvgls()] for both the
#'   observed fit and bootstrap refits.
#'
#' @details
#' The bootstrap null simulates a residual process under the same single-regime
#' model used for fitting: `"BM"` uses `mvSIM(model = "BM1")`, `"OU"` uses
#' `mvSIM(model = "OU1")`, and `"EB"` uses `mvSIM(model = "EB")`. The fitted
#' residual covariance matrix and fitted mean structure from the observed
#' `mvgls` model are reused in each bootstrap replicate.
#'
#' If `focal` is supplied as an internal node number, descendant tips are
#' resolved on the original tree before overlap pruning, so the focal definition
#' remains stable when non-overlapping taxa are dropped during alignment.
#'
#' @return A list of class `bifrost_focal_disparity_test` containing the aligned
#'   data, the fitted `mvgls` object, the focal tip set, the observed statistic,
#'   the bootstrap null distribution, and standardized effect-size and p-value
#'   summaries.
#'
#' @seealso [createSimulationTemplate()]
#'
#' @examples
#' \dontrun{
#' ex <- focalDisparityExampleData(seed = 2)
#'
#' res <- testFocalDisparity(
#'   tree = ex$tree,
#'   trait_data = ex$trait_data,
#'   focal = ex$focal_tips,
#'   formula = cbind(y1, y2, y3) ~ size,
#'   metric = "centroid_msd",
#'   nsim = 100,
#'   method = "LL",
#'   workers = 1,
#'   future_plan = "sequential",
#'   show_progress = FALSE,
#'   seed = 2
#' )
#'
#' ex$focal_node
#' length(ex$focal_tips)
#' print(res)
#' plot(res)
#' }
#'
#' @export
testFocalDisparity <- function(tree,
                               trait_data,
                               focal,
                               formula = "trait_data ~ 1",
                               metric = c("centroid_msd", "trace_cov", "mean_pairwise_sqdist"),
                               statistic = NULL,
                               nsim = 999,
                               model = "BM",
                               method = "PL-LOOCV",
                               workers = future::availableCores(),
                               seed = 1,
                               future_plan = NULL,
                               show_progress = TRUE,
                               min_focal_tips = 3L,
                               ...) {
  metric_supplied <- !missing(metric)
  extra_args <- list(...)
  reserved_args <- intersect(
    c("formula", "tree", "model", "data", "method"),
    names(extra_args)
  )
  if (length(reserved_args) > 0L) {
    stop(
      "Do not pass reserved mvgls arguments via ...: ",
      paste(reserved_args, collapse = ", "),
      "."
    )
  }

  if (!is.null(statistic)) {
    if (metric_supplied) {
      stop("Specify only one of 'metric' and 'statistic'.")
    }
    metric <- statistic
  }

  if (!is.numeric(nsim) || length(nsim) != 1L || is.na(nsim) || nsim < 1) {
    stop("nsim must be a single integer >= 1.")
  }
  if (!is.numeric(workers) || length(workers) != 1L || is.na(workers) || workers < 1) {
    stop("workers must be a single integer >= 1.")
  }

  nsim <- as.integer(nsim)
  workers <- as.integer(workers)
  model <- match.arg(model, choices = c("BM", "OU", "EB"))

  if (is.null(future_plan)) {
    future_plan <- if (workers > 1L) "multisession" else "sequential"
  } else {
    future_plan <- match.arg(future_plan, choices = c("multisession", "multicore", "sequential"))
  }

  aligned <- .bifrost_align_focal_test_inputs(
    tree = tree,
    trait_data = trait_data,
    focal = focal,
    min_focal_tips = min_focal_tips
  )
  tree <- aligned$tree
  trait_data <- aligned$trait_data
  focal_tips <- aligned$focal_tips

  normalized <- normalizeMvglsFormulaCall(
    formula = formula,
    trait_data = trait_data,
    args_list = list(data = trait_data),
    allow_single_response = FALSE
  )

  fit_formula <- normalized$formula
  fit_data <- normalized$args_list$data
  if (is.null(fit_data)) {
    stop("Could not construct a data argument for mvgls.")
  }
  fit_data <- as.data.frame(fit_data, check.names = FALSE, stringsAsFactors = FALSE)
  rownames(fit_data) <- tree$tip.label
  formula_input_chr <- if (inherits(formula, "formula")) {
    paste(deparse(formula, width.cutoff = 500), collapse = " ")
  } else {
    as.character(formula)[1L]
  }
  fit_formula_chr <- paste(deparse(fit_formula, width.cutoff = 500), collapse = " ")

  normalized_extra_args <- normalized$args_list[setdiff(names(normalized$args_list), "data")]
  fit_call_args <- .bifrost_build_mvgls_call_args(
    formula = .bifrost_formula_from_string(fit_formula_chr),
    data = fit_data,
    tree = tree,
    model = model,
    method = method,
    normalized_extra_args = normalized_extra_args,
    extra_args = extra_args
  )

  fit <- .bifrost_call_mvgls(fit_call_args)

  residuals_observed <- as.matrix(residuals(fit))
  rownames(residuals_observed) <- tree$tip.label
  fitted_values <- as.matrix(fit$fitted)
  rownames(fitted_values) <- tree$tip.label
  response_names <- .bifrost_get_response_names(
    fitted_values = fitted_values,
    residuals_matrix = residuals_observed
  )
  colnames(residuals_observed) <- response_names
  colnames(fitted_values) <- response_names

  metric_spec <- .bifrost_make_disparity_metric(metric)
  observed_statistic <- .bifrost_eval_disparity_metric(
    metric_fun = metric_spec$fun,
    X = residuals_observed[focal_tips, , drop = FALSE]
  )

  simulation_model <- .bifrost_resolve_simulation_model(model)
  sigma_hat <- .bifrost_extract_sigma_matrix(fit)
  null_simulation_params <- .bifrost_build_null_simulation_params(
    fit = fit,
    simulation_model = simulation_model,
    n_traits = ncol(fitted_values),
    sigma_hat = sigma_hat
  )
  sim_tree <- ape::reorder.phylo(ape::as.phylo(tree), order = "postorder")
  extract_simulated_matrix <- .bifrost_extract_simulated_matrix
  replace_response_block <- .bifrost_replace_response_block
  build_mvgls_call_args <- .bifrost_build_mvgls_call_args
  formula_from_string <- .bifrost_formula_from_string
  call_mvgls <- .bifrost_call_mvgls
  eval_disparity_metric <- .bifrost_eval_disparity_metric
  metric_fun <- metric_spec$fun

  sim_one <- function(i, progress = NULL) {
    if (!is.null(progress)) {
      progress(sprintf("simulation %d/%d", i, nsim))
    }

    simulated <- mvMORPH::mvSIM(
      tree = sim_tree,
      nsim = 1,
      model = simulation_model,
      param = null_simulation_params
    )

    simulated_residuals <- extract_simulated_matrix(
      simulated = simulated,
      tip_labels = sim_tree$tip.label,
      trait_names = response_names
    )
    simulated_response <- fitted_values + simulated_residuals
    rownames(simulated_response) <- tree$tip.label
    colnames(simulated_response) <- response_names
    simulated_data <- replace_response_block(
      data = fit_data,
      response_names = response_names,
      response_matrix = simulated_response
    )
    simulated_data <- as.data.frame(
      simulated_data,
      check.names = FALSE,
      stringsAsFactors = FALSE
    )
    rownames(simulated_data) <- tree$tip.label

    sim_fit <- call_mvgls(
      build_mvgls_call_args(
        formula = formula_from_string(fit_formula_chr),
        data = simulated_data,
        tree = tree,
        model = model,
        method = method,
        normalized_extra_args = normalized_extra_args,
        extra_args = extra_args
      )
    )

    residuals_simulated <- as.matrix(residuals(sim_fit))
    rownames(residuals_simulated) <- tree$tip.label
    eval_disparity_metric(
      metric_fun = metric_fun,
      X = residuals_simulated[focal_tips, , drop = FALSE]
    )
  }

  old_plan <- future::plan()
  on.exit(future::plan(old_plan), add = TRUE)

  if (identical(future_plan, "multisession")) {
    future::plan(future::multisession, workers = workers)
  } else if (identical(future_plan, "multicore")) {
    future::plan(future::multicore, workers = workers)
  } else {
    future::plan(future::sequential)
  }

  null_statistics <- .bifrost_run_bootstrap(
    nsim = nsim,
    sim_one = sim_one,
    seed = seed,
    show_progress = show_progress,
    use_future = !identical(future_plan, "sequential")
  )

  null_statistics <- unname(as.numeric(null_statistics))
  summary_stats <- .bifrost_summarize_null(
    observed = observed_statistic,
    null_statistics = null_statistics
  )

  out <- list(
    call = match.call(),
    fit = fit,
    tree = tree,
    focal_input = focal,
    trait_data = trait_data,
    aligned_trait_data = trait_data,
    fit_data = fit_data,
    focal_tips = focal_tips,
    formula = formula_input_chr,
    formula_input = formula_input_chr,
    formula_normalized = fit_formula_chr,
    fit_formula = fit_formula_chr,
    metric = metric_spec$label,
    statistic = metric_spec$label,
    observed_statistic = summary_stats$observed,
    null_statistics = null_statistics,
    null_mean = summary_stats$null_mean,
    null_sd = summary_stats$null_sd,
    SES = summary_stats$SES,
    p_upper = summary_stats$p_upper,
    p_lower = summary_stats$p_lower,
    p_two_tailed = summary_stats$p_two_tailed,
    nsim = nsim,
    fit_model = model,
    model = model,
    null_model = model,
    simulation_model = simulation_model,
    method = method,
    workers = workers,
    seed = seed,
    bootstrap_model = simulation_model
  )

  class(out) <- c("bifrost_focal_disparity_test", class(out))
  out
}

# Backward-compatible alias kept unexported while this API settles.
testFocalResidualStatistic <- function(...) {
  testFocalDisparity(...)
}

#' Print method for focal disparity bootstrap tests
#'
#' @param x A `bifrost_focal_disparity_test` object returned by
#'   `testFocalDisparity()`.
#' @param ... Unused (S3 compatibility).
#'
#' @return Invisibly returns `x`.
#'
#' @export
print.bifrost_focal_disparity_test <- function(x, ...) {
  cat("Bifrost Focal Disparity Test\n")
  cat("============================\n\n")
  cat("Model\n")
  cat("  Input formula: ", x$formula_input, "\n", sep = "")
  if (!identical(x$formula_input, x$fit_formula)) {
    cat("  Fit formula: ", x$fit_formula, "\n", sep = "")
  }
  cat("  Fit model: ", x$fit_model, "\n", sep = "")
  cat("  Fit method: ", x$method, "\n", sep = "")
  cat("  Bootstrap model: ", x$null_model, "\n", sep = "")
  cat("  Simulation code: ", x$simulation_model, "\n\n", sep = "")

  cat("Focal Set\n")
  cat("  Tips: ", length(x$focal_tips), "\n", sep = "")
  cat("  Metric: ", x$metric, "\n", sep = "")
  cat("  Simulations: ", x$nsim, "\n\n", sep = "")

  cat("Results\n")
  cat("  Observed: ", formatC(x$observed_statistic, format = "f", digits = 6), "\n", sep = "")
  cat("  Null mean: ", formatC(x$null_mean, format = "f", digits = 6), "\n", sep = "")
  cat("  Null sd: ", formatC(x$null_sd, format = "f", digits = 6), "\n", sep = "")
  cat("  SES: ", formatC(x$SES, format = "f", digits = 6), "\n", sep = "")
  cat("  Upper-tail p: ", formatC(x$p_upper, format = "f", digits = 6), "\n", sep = "")
  cat("  Lower-tail p: ", formatC(x$p_lower, format = "f", digits = 6), "\n", sep = "")
  cat("  Two-tailed p: ", formatC(x$p_two_tailed, format = "f", digits = 6), "\n", sep = "")

  invisible(x)
}

#' Plot method for focal disparity bootstrap tests
#'
#' @param x A `bifrost_focal_disparity_test` object returned by
#'   `testFocalDisparity()`.
#' @param breaks Histogram break specification passed to [graphics::hist()].
#' @param main Plot title. Defaults to `"Null Distribution of Focal Disparity"`.
#' @param xlab X-axis label. Defaults to a metric-aware label.
#' @param ... Additional graphical parameters passed to [graphics::hist()].
#'
#' @return Invisibly returns `x`.
#'
#' @export
plot.bifrost_focal_disparity_test <- function(x,
                                              breaks = 30,
                                              main = "Null Distribution of Focal Disparity",
                                              xlab = NULL,
                                              ...) {
  if (is.null(xlab)) {
    xlab <- if (!is.null(x$metric) && !identical(x$metric, "user_function")) {
      paste0("Bootstrap metric: ", x$metric)
    } else {
      "Bootstrap statistic"
    }
  }

  graphics::hist(
    x$null_statistics,
    breaks = breaks,
    main = main,
    xlab = xlab,
    ...
  )
  graphics::abline(v = x$observed_statistic, lwd = 2)
  invisible(x)
}
