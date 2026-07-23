#' Fit Independent Post-hoc Regime Covariance Models
#'
#' Fit separate `mvMORPH::mvgls()` Brownian-motion models to the tips assigned
#' to each contemporary mapped regime. This is a post-hoc workflow: it estimates
#' independent regime covariance matrices after a mapped-regime analysis has
#' already been chosen.
#'
#' @param x Optional `bifrost_search` object, compatible list with a
#'   `tree_no_uncertainty_untransformed` component, or SIMMAP-style `phylo`
#'   tree. If `tree` is supplied, `x` is only used as an optional source of
#'   regime-rate parameters.
#' @param tree Optional SIMMAP-style `phylo` tree with mapped states, non-empty
#'   regime-state identifiers, and unique, non-empty tip labels. Required when
#'   `x` is `NULL`.
#' @param trait_data Matrix or data frame with unique, non-empty row names
#'   matching the tree tip labels. Column names may be omitted; when present,
#'   they must also be unique and non-empty.
#' @param formula Formula used for each regime-specific `mvgls()` fit. The
#'   default treats `trait_data` as the multivariate response. Formulae may
#'   refer to `trait_data` directly, as in `trait_data[, 1:5] ~ trait_data[, 6]`,
#'   or to columns by name.
#' @param model Evolutionary model passed to [mvMORPH::mvgls()]. The default
#'   `"BM"` matches the independent post-hoc covariance workflow.
#' @param min_tips Minimum number of tips required before a regime is fitted.
#'   Regimes are fitted when `tip_count >= min_tips`; regimes with fewer tips
#'   are returned with status `"skipped"`. The default of two tips matches the
#'   Berv et al. post-hoc fitting script.
#' @param cores Number of workers. Values greater than one use
#'   [future.apply::future_lapply()] with a temporary multisession plan.
#' @param error Logical passed to [mvMORPH::mvgls()].
#' @param ... Additional arguments passed to [mvMORPH::mvgls()].
#'
#' @return An object of class `regime_covariances`, a list with:
#' \describe{
#'   \item{fits}{Named list of regime-specific `mvgls` fits or `NULL`.}
#'   \item{covariances}{Named list of extracted post-hoc covariance matrices or
#'     `NULL` for skipped/failed regimes.}
#'   \item{status}{Data frame with regime ID, tip count, regime age, status, and
#'     message.}
#'   \item{rates}{Named regime-rate vector if available from `x`; otherwise
#'     `NULL`.}
#' }
#'
#' @details
#' The covariance matrices returned here are independent post-hoc estimates.
#' They are not the same object as `search$VCVs` in a `bifrost_search` result:
#' the current `bifrost` search model uses scalar transformations of a shared
#' phenotypic covariance matrix, so `search$VCVs` are proportional joint-model
#' matrices rather than independently refit regime matrices.
#' Extracted matrices must be symmetric and positive semidefinite, contain
#' finite entries and strictly positive diagonal variances, and, when named,
#' have unique matching row and column trait names. Fits that do not satisfy
#' those requirements are retained with status `"failed"` and a diagnostic
#' message.
#'
#' @export
fit_regime_covariances <- function(x = NULL,
                                   tree = NULL,
                                   trait_data,
                                   formula = trait_data ~ 1,
                                   model = "BM",
                                   min_tips = 2L,
                                   cores = 1L,
                                   error = TRUE,
                                   ...) {
  resolved_tree <- .regime_resolve_tree(x = x, tree = tree)
  trait_data <- .regime_validate_trait_data(trait_data, resolved_tree)
  min_tips <- .regime_validate_min_tips(min_tips)
  cores <- .regime_validate_cores(cores)

  tip_states <- .regime_tip_states(resolved_tree)
  tips_by_regime <- split(names(tip_states), unname(tip_states))
  regimes <- names(tips_by_regime)
  subtree_age <- .regime_subtree_age
  formula_with_data <- .regime_formula_with_data
  extract_covariance <- .regime_extract_covariance
  validate_covariance <- .regime_validate_covariance_matrix

  fit_one <- function(regime) {
    tips <- tips_by_regime[[regime]]
    tip_count <- length(tips)
    regime_age <- subtree_age(resolved_tree, tips)

    if (tip_count < min_tips) {
      return(list(
        regime = regime,
        tip_count = tip_count,
        regime_age = regime_age,
        fit = NULL,
        covariance = NULL,
        status = "skipped",
        message = sprintf(
          "Regime has %d tip%s; min_tips is %d.",
          tip_count,
          if (identical(tip_count, 1L)) "" else "s",
          min_tips
        )
      ))
    }

    subtree <- tryCatch(
      ape::keep.tip(resolved_tree, tips),
      error = function(e) e
    )
    if (inherits(subtree, "error")) { # nocov start
      return(list(
        regime = regime,
        tip_count = tip_count,
        regime_age = regime_age,
        fit = NULL,
        covariance = NULL,
        status = "failed",
        message = paste("Failed to extract regime subtree:", subtree$message)
      ))
    } # nocov end

    subset_data <- trait_data[
      match(subtree$tip.label, rownames(trait_data)),
      ,
      drop = FALSE
    ]

    formula_obj <- formula_with_data(formula, subset_data)
    fit <- tryCatch(
      mvMORPH::mvgls(
        formula_obj,
        tree = subtree,
        model = model,
        error = error,
        ...
      ),
      error = function(e) e
    )

    if (inherits(fit, "error")) {
      return(list(
        regime = regime,
        tip_count = tip_count,
        regime_age = regime_age,
        fit = NULL,
        covariance = NULL,
        status = "failed",
        message = paste(
          "Formula/model fit failed for",
          paste(deparse(formula_obj), collapse = " "),
          ":",
          fit$message
        )
      ))
    }

    covariance <- extract_covariance(fit)
    if (is.null(covariance)) { # nocov start
      return(list(
        regime = regime,
        tip_count = tip_count,
        regime_age = regime_age,
        fit = fit,
        covariance = NULL,
        status = "failed",
        message = "Fitted model did not contain a usable covariance matrix."
      ))
    } # nocov end

    covariance_error <- tryCatch(
      {
        validate_covariance(covariance, regime)
        NULL
      },
      error = function(e) e
    )
    if (inherits(covariance_error, "error")) {
      return(list(
        regime = regime,
        tip_count = tip_count,
        regime_age = regime_age,
        fit = fit,
        covariance = NULL,
        status = "failed",
        message = paste(
          "Fitted model returned an unusable covariance matrix:",
          conditionMessage(covariance_error)
        )
      ))
    }

    list(
      regime = regime,
      tip_count = tip_count,
      regime_age = regime_age,
      fit = fit,
      covariance = covariance,
      status = "ok",
      message = NA_character_
    )
  }

  results <- if (cores > 1L) {
    old_plan <- future::plan()
    on.exit(future::plan(old_plan), add = TRUE)
    future::plan(future::multisession, workers = cores)
    future.apply::future_lapply(regimes, fit_one, future.seed = TRUE)
  } else {
    lapply(regimes, fit_one)
  }
  names(results) <- regimes
  parallel_strategy <- if (cores > 1L) {
    "future.apply::future_lapply"
  } else {
    "base::lapply"
  }

  fits <- setNames(lapply(results, `[[`, "fit"), regimes)
  covariances <- setNames(lapply(results, `[[`, "covariance"), regimes)
  status <- data.frame(
    regime = vapply(results, `[[`, character(1), "regime"),
    tip_count = vapply(results, `[[`, integer(1), "tip_count"),
    regime_age = vapply(results, `[[`, numeric(1), "regime_age"),
    status = vapply(results, `[[`, character(1), "status"),
    message = vapply(results, `[[`, character(1), "message"),
    row.names = NULL,
    check.names = FALSE
  )

  structure(
    list(
      fits = fits,
      covariances = covariances,
      status = status,
      rates = .regime_extract_rates(x),
      formula = formula,
      model = model,
      min_tips = min_tips,
      cores = cores,
      parallel_strategy = parallel_strategy
    ),
    class = c("regime_covariances", "list")
  )
}

#' Fit Post-hoc Regime Covariance Models Across Runs
#'
#' Apply [fit_regime_covariances()] to each search/run in a named list. This is
#' the package analogue of the manuscript `fitPosthocModels()` step: the
#' top-level run iteration is explicit, while each run can still fit its
#' regime subtrees in parallel through `cores`.
#'
#' @param x Named list of `bifrost_search` objects, compatible search-like
#'   lists, or SIMMAP-style `phylo` trees. Missing or blank names are generated
#'   and names must be unique after normalization.
#' @inheritParams fit_regime_covariances
#' @param tree_element For search-like list inputs, the element containing the
#'   mapped tree to use for post-hoc refits. Set to `NULL` to let
#'   [fit_regime_covariances()] resolve the tree from each run directly.
#' @param verbose Logical; if `TRUE`, print progress messages by run name.
#'
#' @return A named list of `regime_covariances` objects with class
#'   `regime_covariance_runs`.
#' @export
fit_regime_covariance_runs <- function(x,
                                       trait_data,
                                       formula = trait_data ~ 1,
                                       model = "BM",
                                       min_tips = 2L,
                                       cores = 1L,
                                       error = TRUE,
                                       tree_element = "tree_no_uncertainty_untransformed",
                                       verbose = FALSE,
                                       ...) {
  runs <- .regime_validate_run_list(x, "x")
  verbose <- .regime_validate_logical(verbose, "verbose")
  out <- stats::setNames(vector("list", length(runs)), names(runs))

  for (run_name in names(runs)) {
    if (verbose) {
      message("Fitting post-hoc regime covariances for: ", run_name)
    }
    run <- runs[[run_name]]
    tree <- NULL
    if (!is.null(tree_element) &&
        !inherits(run, "phylo") &&
        is.list(run) &&
        tree_element %in% names(run)) {
      tree <- run[[tree_element]]
    }
    out[[run_name]] <- fit_regime_covariances(
      x = run,
      tree = tree,
      trait_data = trait_data,
      formula = formula,
      model = model,
      min_tips = min_tips,
      cores = cores,
      error = error,
      ...
    )
  }

  structure(
    out,
    class = c("regime_covariance_runs", "list"),
    settings = list(
      formula = formula,
      model = model,
      min_tips = min_tips,
      cores = cores,
      tree_element = tree_element
    )
  )
}

#' Summarize Post-hoc Regime Covariance Matrices
#'
#' Compute regime-rate, mean variance, mean absolute trait correlation,
#' Fisher-Z transformed mean absolute correlation, tip count, and regime age
#' from independent post-hoc covariance matrices.
#'
#' @param x A `regime_covariances` object returned by
#'   [fit_regime_covariances()] or a named list of covariance/correlation
#'   matrices.
#' @param search Optional `bifrost_search` object used as a source of named
#'   regime-rate parameters and, when possible, mapped-regime tip counts/ages.
#' @param rates Optional named numeric vector of regime rates. Names must be
#'   unique and non-empty and are matched to regime IDs; numeric equality is
#'   never used for matching.
#' @param tree Optional SIMMAP-style tree used to compute tip counts and regime
#'   ages when they are not already supplied.
#' @param tip_counts Optional named numeric vector of tip counts. When names are
#'   supplied, they must be unique and non-empty.
#' @param regime_ages Optional named numeric vector of regime ages. When names
#'   are supplied, they must be unique and non-empty.
#' @param fisher_boundary How to handle correlations with `abs(r) >= 1`, where
#'   Fisher-Z is undefined. `"NA"` returns `NA` for those regimes; `"error"`
#'   stops with an error.
#' @param remove_high_corr Logical; if `TRUE`, drop rows whose mean absolute
#'   correlation exceeds `corr_threshold`. This reproduces the manuscript
#'   `generateVarsCorsList(remove_high_corr = TRUE, corr_threshold = 0.95)`
#'   filtering step.
#' @param corr_threshold Correlation threshold used when `remove_high_corr` is
#'   `TRUE`.
#'
#' @return A data frame with one row per regime and columns `regime`, `rate`,
#'   `mean_variance`, `mean_abs_correlation`,
#'   `fisher_z_mean_abs_correlation`, `tip_count`, `regime_age`, `status`, and
#'   `message`.
#'
#' @details
#' This summary is designed for independent post-hoc matrices. Do not pass
#' proportional `search$VCVs` from a `bifrost_search` object when the question is
#' whether regimes differ in phenotypic integration or correlation structure.
#' If `x` is detectably the same object as `search$VCVs`, the function warns and
#' still returns the requested descriptive summary.
#' Raw matrix-list inputs must have unique, non-empty regime names. Their
#' matrices must be symmetric and positive semidefinite, contain finite entries
#' and strictly positive diagonal variances, and, when named, have unique
#' matching row and column trait names; violations produce a regime-specific
#' error. Invalid matrices
#' encountered in a `regime_covariances` fit object are instead represented as
#' failed rows with missing summaries and a diagnostic message.
#' A one-trait matrix has no pairwise correlations, so its mean absolute
#' correlation and Fisher-Z summary are returned as `NA`.
#'
#' @export
summarize_regime_covariances <- function(x,
                                         search = NULL,
                                         rates = NULL,
                                         tree = NULL,
                                         tip_counts = NULL,
                                         regime_ages = NULL,
                                         fisher_boundary = c("NA", "error"),
                                         remove_high_corr = FALSE,
                                         corr_threshold = 0.95) {
  fisher_boundary <- match.arg(fisher_boundary)
  remove_high_corr <- .regime_validate_logical(remove_high_corr, "remove_high_corr")
  corr_threshold <- .regime_validate_corr_threshold(corr_threshold)
  extracted <- .regime_extract_covariance_list(x)
  .regime_warn_if_search_vcvs(x, search)
  covariances <- extracted$covariances
  status <- extracted$status
  regimes <- names(covariances)

  if (is.null(rates)) {
    rates <- extracted$rates
  }
  if (is.null(rates)) {
    rates <- .regime_extract_rates(search)
  }
  rates <- .regime_match_rates(rates, regimes)

  if (is.null(tree) && !is.null(search)) {
    tree <- tryCatch(.regime_resolve_tree(x = search), error = function(e) NULL)
  }
  if (is.null(tip_counts) &&
      !is.null(status$tip_count) &&
      any(!is.na(status$tip_count))) {
    tip_counts <- stats::setNames(status$tip_count, status$regime)
  }
  if (is.null(regime_ages) &&
      !is.null(status$regime_age) &&
      any(!is.na(status$regime_age))) {
    regime_ages <- stats::setNames(status$regime_age, status$regime)
  }
  if ((is.null(tip_counts) || is.null(regime_ages)) && !is.null(tree)) {
    tree_stats <- .regime_tree_stats(tree)
    if (is.null(tip_counts)) {
      tip_counts <- tree_stats$tip_counts
    }
    if (is.null(regime_ages)) {
      regime_ages <- tree_stats$regime_ages
    }
  }

  tip_counts <- .regime_match_optional(tip_counts, regimes)
  regime_ages <- .regime_match_optional(regime_ages, regimes)

  out <- lapply(regimes, function(regime) {
    mat <- covariances[[regime]]
    status_i <- status$status[match(regime, status$regime)]
    message_i <- status$message[match(regime, status$regime)]
    invalid_message <- extracted$invalid_messages[[regime]]
    if (length(invalid_message) == 1L && !is.na(invalid_message)) {
      status_i <- "failed"
      message_i <- invalid_message
    }
    if (length(status_i) == 0L || is.na(status_i)) {
      status_i <- if (is.matrix(mat)) "ok" else "failed"
    }
    if (length(message_i) == 0L) {
      message_i <- NA_character_
    }

    if (!is.matrix(mat)) {
      return(data.frame(
        regime = regime,
        rate = rates[[regime]],
        mean_variance = NA_real_,
        mean_abs_correlation = NA_real_,
        fisher_z_mean_abs_correlation = NA_real_,
        tip_count = tip_counts[[regime]],
        regime_age = regime_ages[[regime]],
        status = status_i,
        message = message_i,
        row.names = NULL,
        check.names = FALSE
      ))
    }

    .regime_validate_covariance_matrix(mat, regime)
    mean_variance <- mean(diag(mat), na.rm = TRUE)
    cor_mat <- stats::cov2cor(mat)
    pairwise_correlations <- cor_mat[upper.tri(cor_mat, diag = FALSE)]
    mean_abs_correlation <- if (length(pairwise_correlations) == 0L) {
      NA_real_
    } else {
      mean(abs(pairwise_correlations), na.rm = TRUE)
    }
    fisher_value <- .regime_fisher_z(mean_abs_correlation, fisher_boundary)

    data.frame(
      regime = regime,
      rate = rates[[regime]],
      mean_variance = mean_variance,
      mean_abs_correlation = mean_abs_correlation,
      fisher_z_mean_abs_correlation = fisher_value,
      tip_count = tip_counts[[regime]],
      regime_age = regime_ages[[regime]],
      status = status_i,
      message = message_i,
      row.names = NULL,
      check.names = FALSE
    )
  })

  out <- do.call(rbind, out)
  if (remove_high_corr) {
    out <- out[
      !is.na(out$mean_abs_correlation) &
        abs(out$mean_abs_correlation) <= corr_threshold,
      ,
      drop = FALSE
    ]
    rownames(out) <- NULL
  }
  out
}

#' Summarize Post-hoc Regime Covariance Models Across Runs
#'
#' Apply [summarize_regime_covariances()] to each element returned by
#' [fit_regime_covariance_runs()]. This mirrors the manuscript
#' `generateVarsCorsList()` step while keeping rate matching, tip counts, and
#' high-correlation filtering explicit and reusable.
#'
#' @param x Named list of `regime_covariances` objects, usually from
#'   [fit_regime_covariance_runs()]. Missing or blank names are generated and
#'   names must be unique after normalization.
#' @param searches Optional named list of corresponding `bifrost_search`
#'   objects. Run names are matched to `x`.
#' @param rates Optional named list of regime-rate vectors, or one named rate
#'   vector to reuse for every run.
#' @param tree,tip_counts,regime_ages Optional named lists of per-run values,
#'   or a single value to reuse for every run.
#' @inheritParams summarize_regime_covariances
#'
#' @return A named list of summary data frames with class
#'   `regime_covariance_run_summaries`.
#' @export
summarize_regime_covariance_runs <- function(x,
                                             searches = NULL,
                                             rates = NULL,
                                             tree = NULL,
                                             tip_counts = NULL,
                                             regime_ages = NULL,
                                             fisher_boundary = c("NA", "error"),
                                             remove_high_corr = FALSE,
                                             corr_threshold = 0.95) {
  fits <- .regime_validate_run_list(x, "x")
  run_names <- names(fits)
  searches <- .regime_match_run_argument(searches, run_names, "searches")
  rates <- .regime_match_run_argument(rates, run_names, "rates")
  tree <- .regime_match_run_argument(tree, run_names, "tree")
  tip_counts <- .regime_match_run_argument(tip_counts, run_names, "tip_counts")
  regime_ages <- .regime_match_run_argument(regime_ages, run_names, "regime_ages")

  out <- stats::setNames(vector("list", length(fits)), run_names)
  for (run_name in run_names) {
    out[[run_name]] <- summarize_regime_covariances(
      fits[[run_name]],
      search = searches[[run_name]],
      rates = rates[[run_name]],
      tree = tree[[run_name]],
      tip_counts = tip_counts[[run_name]],
      regime_ages = regime_ages[[run_name]],
      fisher_boundary = fisher_boundary,
      remove_high_corr = remove_high_corr,
      corr_threshold = corr_threshold
    )
  }

  structure(
    out,
    class = c("regime_covariance_run_summaries", "list"),
    settings = list(
      remove_high_corr = remove_high_corr,
      corr_threshold = corr_threshold
    )
  )
}

.regime_check_unused_dots <- function(dots) {
  if (length(dots) > 0L) {
    dot_names <- names(dots)
    if (is.null(dot_names)) {
      dot_names <- rep.int("", length(dots))
    }
    missing_names <- is.na(dot_names) | !nzchar(trimws(dot_names))
    dot_names[missing_names] <- paste0("..", which(missing_names))
    stop(
      "Unused argument(s): ",
      paste(dot_names, collapse = ", "),
      call. = FALSE
    )
  }
}

#' Run PCA on Regime Correlation Structures
#'
#' Convert each regime covariance matrix to a correlation matrix, vectorize the
#' upper-triangle off-diagonal entries, and run [stats::prcomp()] on the
#' resulting regime-by-correlation matrix. When `use_correlation = TRUE`, the
#' function warns if all retained covariance matrices are scalar-proportional to
#' one another. That pattern is expected for proportional `search$VCVs` from the
#' scalar `bifrost` search model and should not be interpreted as independent
#' post-hoc covariance reconfiguration evidence.
#'
#' @param x A `regime_covariances` object returned by
#'   [fit_regime_covariances()] or a named list of covariance/correlation
#'   matrices. Matrices must contain finite entries and strictly positive
#'   diagonal variances.
#' @param use_correlation Logical; if `TRUE`, convert each matrix with
#'   [stats::cov2cor()] before vectorizing.
#' @param center,scale. Passed to [stats::prcomp()]. Defaults match the
#'   manuscript workflow.
#' @param tip_counts Optional named numeric vector of tip counts.
#' @param regime_ages Optional named numeric vector of regime ages.
#' @param trait_labels Optional labels for the matrix rows/columns. Supply a
#'   named vector to map internal trait names to display labels, or an unnamed
#'   vector in matrix order. Resolved labels must be unique and not blank.
#' @param min_tips Optional positive whole-number minimum tip count for
#'   inclusion when `tip_counts` are available. Matching the manuscript `min_n`
#'   convention, regimes are retained only when `tip_count > min_tips`. At
#'   least one non-missing tip count must be available when this filter is used.
#' @param ... Reserved for future extensions.
#'
#' @return An object of class `regime_correlation_pca` containing the `prcomp`
#'   object, vectorized input matrix, scores, loadings, variance explained,
#'   regime IDs, trait labels, optional diagnostics, and a `settings` list that
#'   records the PCA scaling, tip-count filter, and whether retained covariance
#'   matrices were scalar-proportional.
#'
#' @export
regime_correlation_pca <- function(x,
                                   use_correlation = TRUE,
                                   center = TRUE,
                                   scale. = TRUE,
                                   tip_counts = NULL,
                                   regime_ages = NULL,
                                   trait_labels = NULL,
                                   min_tips = NULL,
                                   ...) {
  dots <- list(...)
  .regime_check_unused_dots(dots)
  if (!is.null(min_tips)) {
    min_tips <- .regime_validate_min_tips(min_tips)
  }

  extracted <- .regime_extract_covariance_list(x)
  matrices <- extracted$covariances
  if (is.null(tip_counts) && !is.null(extracted$status$tip_count)) {
    tip_counts <- stats::setNames(extracted$status$tip_count, extracted$status$regime)
  }
  if (is.null(regime_ages) && !is.null(extracted$status$regime_age)) {
    regime_ages <- stats::setNames(extracted$status$regime_age, extracted$status$regime)
  }

  regimes <- names(matrices)
  tip_counts <- .regime_match_optional(tip_counts, regimes)
  regime_ages <- .regime_match_optional(regime_ages, regimes)
  if (!is.null(min_tips) && all(is.na(tip_counts))) {
    stop(
      "`min_tips` requires available `tip_counts`; supply at least one ",
      "non-missing tip count.",
      call. = FALSE
    )
  }

  keep <- vapply(matrices, is.matrix, logical(1))
  if (!is.null(min_tips)) {
    keep <- keep & !is.na(tip_counts) & tip_counts > min_tips
  }
  matrices <- matrices[keep]
  regimes <- names(matrices)
  tip_counts <- tip_counts[regimes]
  regime_ages <- regime_ages[regimes]

  if (length(matrices) < 2L) {
    stop("At least two valid regime matrices are required for PCA.", call. = FALSE)
  }

  template <- matrices[[1L]]
  .regime_validate_covariance_matrix(template, regimes[[1L]])
  internal_traits <- rownames(template)
  if (is.null(internal_traits)) {
    internal_traits <- colnames(template)
  }
  if (is.null(internal_traits)) {
    internal_traits <- paste0("trait", seq_len(nrow(template)))
  }
  labels <- .regime_resolve_trait_labels(trait_labels, internal_traits)
  .regime_validate_identifiers(labels, "trait labels")
  upper_names <- .regime_upper_pair_names(labels)

  aligned_matrices <- lapply(regimes, function(regime) {
    mat <- matrices[[regime]]
    .regime_validate_covariance_matrix(mat, regime)
    if (!identical(dim(mat), dim(template))) {
      stop("All matrices must have the same dimensions.", call. = FALSE)
    }
    .regime_align_matrix_traits(
      mat = mat,
      template_traits = internal_traits,
      regime = regime
    )
  })
  names(aligned_matrices) <- regimes

  proportional_covariances <- .regime_covariance_list_is_scalar_proportional(aligned_matrices)
  if (use_correlation && proportional_covariances) {
    warning(
      "`x` contains scalar-proportional covariance matrices. These may be ",
      "proportional `search$VCVs` from a scalar `bifrost` search and should ",
      "not be interpreted as independent post-hoc covariance reconfiguration ",
      "evidence.",
      call. = FALSE
    )
  }

  vectorized <- lapply(aligned_matrices, function(mat) {
    if (use_correlation) {
      mat <- stats::cov2cor(mat)
    }
    mat[upper.tri(mat, diag = FALSE)]
  })

  vectorized_matrix <- do.call(rbind, vectorized)
  rownames(vectorized_matrix) <- regimes
  colnames(vectorized_matrix) <- upper_names

  pca <- stats::prcomp(vectorized_matrix, center = center, scale. = scale.)
  variance_explained <- pca$sdev^2 / sum(pca$sdev^2)

  structure(
    list(
      pca = pca,
      vectorized_matrix = vectorized_matrix,
      scores = pca$x,
      loadings = pca$rotation,
      variance_explained = variance_explained,
      regime_ids = regimes,
      trait_labels = labels,
      tip_counts = tip_counts,
      regime_ages = regime_ages,
      use_correlation = use_correlation,
      center = center,
      scale = scale.,
      settings = list(
        use_correlation = use_correlation,
        center = center,
        scale = scale.,
        min_tips = min_tips,
        tip_filter = if (is.null(min_tips)) NULL else "tip_count > min_tips",
        proportional_covariances = proportional_covariances
      )
    ),
    class = c("regime_correlation_pca", "list")
  )
}

#' Diagnose PCA Axes with User-defined Trait Modules
#'
#' Summarize within- and between-module correlations for each post-hoc regime
#' and correlate those module summaries with regime correlation PCA scores.
#'
#' @param pca A `regime_correlation_pca` object from
#'   [regime_correlation_pca()] computed with `use_correlation = TRUE`.
#' @param modules Named list of character vectors. Module names must be unique
#'   and not blank. Each element names the unique trait labels belonging to one
#'   anatomical, developmental, or functional module. Singleton modules are
#'   allowed; their within-module scores and PC correlations are `NA` because
#'   no pairwise correlation is defined.
#' @param comparisons Optional named list defining module comparisons to score.
#'   Each element must be a length-two character vector naming entries in
#'   `modules`. When both names are the same, the score is the mean
#'   upper-triangle within-module correlation; otherwise it is the mean
#'   between-module correlation. If `NULL`, all within-module and pairwise
#'   between-module comparisons are generated. Missing or blank comparison
#'   names are generated from their module pairs; the resulting names must be
#'   unique and not blank.
#' @param pcs Principal components to correlate with module summaries. Supply
#'   numeric indices or names such as `"PC1"`. Defaults to all PCA score
#'   columns.
#'
#' @return An object of class `regime_module_diagnostics`, containing
#'   per-regime module scores, PC/module correlations, comparison definitions,
#'   and settings.
#' @export
regime_module_diagnostics <- function(pca,
                                      modules,
                                      comparisons = NULL,
                                      pcs = NULL) {
  if (!inherits(pca, "regime_correlation_pca")) {
    stop("`pca` must be a `regime_correlation_pca` object.", call. = FALSE)
  }
  if (!isTRUE(pca$use_correlation)) {
    stop(
      "`regime_module_diagnostics()` requires a correlation-mode PCA. ",
      "Recompute with `regime_correlation_pca(..., use_correlation = TRUE)`.",
      call. = FALSE
    )
  }
  modules <- .regime_validate_modules(modules, pca$trait_labels)
  comparisons <- .regime_validate_module_comparisons(comparisons, modules)

  pc_names <- colnames(pca$scores)
  pc_idx <- if (is.null(pcs)) {
    seq_along(pc_names)
  } else {
    .regime_pca_component_indices(pcs, pc_names)
  }

  regime_matrices <- lapply(seq_len(nrow(pca$vectorized_matrix)), function(i) {
    .regime_vector_to_symmetric_matrix(
      pca$vectorized_matrix[i, ],
      trait_labels = pca$trait_labels
    )
  })
  names(regime_matrices) <- pca$regime_ids

  module_scores <- data.frame(
    regime = pca$regime_ids,
    row.names = NULL,
    check.names = FALSE
  )
  for (comparison_name in names(comparisons)) {
    modules_i <- comparisons[[comparison_name]]
    module_scores[[comparison_name]] <- vapply(
      regime_matrices,
      .regime_module_comparison_score,
      numeric(1),
      module_a = modules[[modules_i[[1L]]]],
      module_b = modules[[modules_i[[2L]]]]
    )
  }

  correlation_rows <- lapply(pc_names[pc_idx], function(pc) {
    data.frame(
      PC = pc,
      module_comparison = names(comparisons),
      correlation = vapply(
        names(comparisons),
        function(comparison_name) {
          .regime_complete_correlation(
            pca$scores[, pc],
            module_scores[[comparison_name]]
          )
        },
        numeric(1)
      ),
      row.names = NULL,
      check.names = FALSE
    )
  })
  correlations <- do.call(rbind, correlation_rows)

  plot_data <- data.frame(
    regime = pca$regime_ids,
    pca$scores[, pc_idx, drop = FALSE],
    module_scores[, names(comparisons), drop = FALSE],
    row.names = NULL,
    check.names = FALSE
  )

  structure(
    list(
      module_scores = module_scores,
      correlations = correlations,
      plot_data = plot_data,
      modules = modules,
      comparisons = comparisons,
      pcs = pc_names[pc_idx],
      settings = list(
        comparisons = comparisons,
        pcs = pc_names[pc_idx]
      )
    ),
    class = c("regime_module_diagnostics", "list")
  )
}

#' Fit the Manuscript-style Regime Integration pGLS
#'
#' Reproduce the representative pGLS test used in the post-hoc integration
#' analysis: log regime rate is modeled as a function of log post-hoc mean
#' variance and Fisher-Z transformed mean absolute trait correlation on a
#' collapsed regime phylogeny.
#'
#' @param summary_data A data frame from [summarize_regime_covariances()] or a
#'   manuscript-compatible `vars_cors` table with columns `rate`, `vars`,
#'   `corrs`, and `State`. Regime IDs must be unique and non-empty.
#' @param search A `bifrost_search` object containing the mapped regime tree.
#' @param tree Optional SIMMAP-style mapped tree. Ignored when `search` is
#'   supplied.
#' @param model Evolutionary model passed to [phylolm::phylolm()]. Defaults to
#'   `"BM"`, matching the manuscript.
#' @param min_tips Optional minimum tip count for downstream inclusion when
#'   `summary_data` contains a `tip_count` column. Downstream summaries are
#'   retained when `tip_count >= min_tips`; this differs from the strict PCA
#'   manuscript filter used by [regime_correlation_pca()].
#' @param ... Additional arguments passed to [phylolm::phylolm()].
#'
#' @return A `phylolm` fit.
#' @details The collapse step follows the manuscript implementation: a
#'   monophyletic regime is represented by one collapsed tip, whereas all tips
#'   assigned to a nonmonophyletic regime are removed. When removals occur, the
#'   function emits one warning listing every dropped regime ID. The collapse
#'   stops if regime relabeling would create duplicated output tip labels.
#' @export
regime_integration_pgls <- function(summary_data,
                                    search = NULL,
                                    tree = NULL,
                                    model = "BM",
                                    min_tips = NULL,
                                    ...) {
  if (!requireNamespace("phylolm", quietly = TRUE)) {
    stop("Package `phylolm` is required for `regime_integration_pgls()`.", call. = FALSE)
  }
  if (is.null(tree)) {
    tree <- .regime_resolve_tree(x = search)
  } else {
    tree <- .regime_resolve_tree(tree = tree)
  }

  data <- .regime_standardize_summary_data(summary_data)
  data <- .regime_filter_summary_min_tips(
    data,
    min_tips = min_tips,
    caller = "regime_integration_pgls"
  )
  phy <- .collapse_regime_phylogeny(tree)
  data$z_log_rate <- as.numeric(scale(log(data$rate)))
  data$z_log_vars <- as.numeric(scale(log(data$vars)))
  data$z_fisher_corr <- as.numeric(scale(fisher_z_transform(data$corrs)))
  data <- data[data$regime %in% phy$tip.label, , drop = FALSE]
  rownames(data) <- data$regime
  if (length(setdiff(phy$tip.label, data$regime)) > 0L) {
    phy <- ape::drop.tip(phy, setdiff(phy$tip.label, data$regime))
  }

  if (nrow(data) < 3L) {
    stop("At least three collapsed-regime rows are required for pGLS.", call. = FALSE)
  }

  fit <- phylolm::phylolm(
    z_log_rate ~ z_log_vars + z_fisher_corr,
    data = data,
    phy = phy,
    model = model,
    ...
  )
  old_terms <- c("z_log_vars", "z_fisher_corr")
  new_terms <- c("scale(log(vars))", "scale(fisher_z_transform(corrs))")
  names(fit$coefficients) <- sub(old_terms[1L], new_terms[1L], names(fit$coefficients), fixed = TRUE)
  names(fit$coefficients) <- sub(old_terms[2L], new_terms[2L], names(fit$coefficients), fixed = TRUE)
  if (!is.null(dimnames(fit$vcov))) {
    dimnames(fit$vcov) <- lapply(dimnames(fit$vcov), function(nm) {
      nm <- sub(old_terms[1L], new_terms[1L], nm, fixed = TRUE)
      sub(old_terms[2L], new_terms[2L], nm, fixed = TRUE)
    })
  }
  fit
}

#' Analyze Rate-Integration Relationships Across Regimes
#'
#' Prepare the point sets and bootstrap confidence curves used for the
#' Supplementary Figure 4A-style rate-vs-variance and rate-vs-integration
#' panels. Use `plot()` on the returned object to draw the panels.
#'
#' @param summaries A list of manuscript-compatible `vars_cors` tables or a data
#'   frame from [summarize_regime_covariances()]. Regime IDs within each input
#'   table must be unique and non-empty. List inputs are pooled with a `run`
#'   column; missing or blank list names become `run1`, `run2`, and so on based
#'   on their positions, and the normalized names must be unique.
#' @param resid_sd_threshold_vars,resid_sd_threshold_corrs Non-negative finite
#'   studentized-residual thresholds used to filter the variance and
#'   correlation panels.
#' @param n_boot Number of bootstrap replicates for confidence curves.
#' @param ci_level Confidence level for bootstrap ribbons.
#' @param seed Optional random seed for reproducible bootstrap curves.
#' @param min_tips Optional minimum tip count for inclusion when summaries
#'   contain a `tip_count` column. Relationship summaries are retained when
#'   `tip_count >= min_tips`; this differs from the strict PCA manuscript filter
#'   used by [regime_correlation_pca()].
#'
#' @return An object of class `regime_integration_relationships`, containing
#'   plot-ready data frames, fitted `lm` objects, bootstrap curves, and settings.
#'
#' @details
#' High-correlation filtering is applied before this step with
#' [summarize_regime_covariances()] via `remove_high_corr` and
#' `corr_threshold`. The relationship object records residual filters and
#' bootstrap settings, but it does not reapply matrix-level correlation filters.
#' Rows with incomplete data for a panel remain in `combined` with an `NA`
#' studentized residual and are omitted from that panel's point and removed-row
#' data frames. Non-missing rates and variances must be finite and strictly
#' positive, and non-missing correlations must be finite values in `[-1, 1]`.
#' Boundary correlations at `-1` or `1` are retained in `combined`, but their
#' undefined Fisher-Z transforms and correlation-panel residuals are `NA`.
#' @export
regime_integration_relationships <- function(summaries,
                                             resid_sd_threshold_vars = 2,
                                             resid_sd_threshold_corrs = 2,
                                             n_boot = 1000,
                                             ci_level = 0.99,
                                             seed = NULL,
                                             min_tips = NULL) {
  resid_sd_threshold_vars <- .regime_validate_resid_sd_threshold(
    resid_sd_threshold_vars,
    "resid_sd_threshold_vars"
  )
  resid_sd_threshold_corrs <- .regime_validate_resid_sd_threshold(
    resid_sd_threshold_corrs,
    "resid_sd_threshold_corrs"
  )
  validated_seed <- .regime_validate_module_plot_seed(seed)
  correlation_seed <- if (is.null(validated_seed)) {
    NULL
  } else if (validated_seed == .Machine$integer.max) {
    0L
  } else {
    validated_seed + 1L
  }
  combined <- .regime_bind_summary_runs(summaries)
  combined <- .regime_filter_summary_min_tips(
    combined,
    min_tips = min_tips,
    caller = "regime_integration_relationships"
  )
  .regime_validate_relationship_values(combined)
  combined$fisher_z_corr <- vapply(
    combined$corrs,
    .regime_fisher_z,
    numeric(1),
    boundary = "NA"
  )
  combined$log_rate <- log(combined$rate)
  combined$log_vars <- log(combined$vars)

  lm_vars <- stats::lm(
    log_rate ~ log_vars,
    data = combined,
    na.action = stats::na.exclude
  )
  combined$vars_resid <- stats::rstudent(lm_vars)
  variance_evaluable <- !is.na(combined$vars_resid)
  variance_points <- combined[
    variance_evaluable &
      abs(combined$vars_resid) <= resid_sd_threshold_vars,
    ,
    drop = FALSE
  ]
  variance_removed <- combined[
    variance_evaluable &
      abs(combined$vars_resid) > resid_sd_threshold_vars,
    ,
    drop = FALSE
  ]

  lm_corrs <- stats::lm(
    log_rate ~ fisher_z_corr,
    data = combined,
    na.action = stats::na.exclude
  )
  combined$corrs_resid <- stats::rstudent(lm_corrs)
  correlation_evaluable <- !is.na(combined$corrs_resid)
  correlation_points <- combined[
    correlation_evaluable &
      abs(combined$corrs_resid) <= resid_sd_threshold_corrs,
    ,
    drop = FALSE
  ]
  correlation_removed <- combined[
    correlation_evaluable &
      abs(combined$corrs_resid) > resid_sd_threshold_corrs,
    ,
    drop = FALSE
  ]

  variance_lm <- stats::lm(log_rate ~ log_vars, data = variance_points)
  correlation_lm <- stats::lm(log_rate ~ fisher_z_corr, data = correlation_points)

  variance_curve <- .regime_bootstrap_curve(
    data = variance_points,
    x_col = "log_vars",
    y_col = "log_rate",
    formula = log_rate ~ log_vars,
    n_boot = n_boot,
    ci_level = ci_level,
    seed = seed
  )
  names(variance_curve)[names(variance_curve) == "x"] <- "log_vars"

  correlation_curve <- .regime_bootstrap_curve(
    data = correlation_points,
    x_col = "fisher_z_corr",
    y_col = "log_rate",
    formula = log_rate ~ fisher_z_corr,
    n_boot = n_boot,
    ci_level = ci_level,
    seed = correlation_seed
  )
  names(correlation_curve)[names(correlation_curve) == "x"] <- "fisher_z_corr"

  structure(
    list(
      combined = combined,
      variance_points = variance_points,
      correlation_points = correlation_points,
      variance_removed = variance_removed,
      correlation_removed = correlation_removed,
      variance_curve = variance_curve,
      correlation_curve = correlation_curve,
      variance_lm = variance_lm,
      correlation_lm = correlation_lm,
      settings = list(
        resid_sd_threshold_vars = resid_sd_threshold_vars,
        resid_sd_threshold_corrs = resid_sd_threshold_corrs,
        n_boot = as.integer(n_boot),
        ci_level = ci_level,
        seed = seed,
        min_tips = min_tips
      )
    ),
    class = c("regime_integration_relationships", "list")
  )
}

#' Plot Regime Rate-Integration Relationships
#'
#' Draw Supplementary Figure 4A-style panels from
#' [regime_integration_relationships()].
#'
#' @param x A `regime_integration_relationships` object.
#' @param panel Which panel to draw: `"both"`, `"variance"`, or
#'   `"correlation"`.
#' @param main Optional panel title. For `panel = "both"`, supply a length-two
#'   vector to title the variance and correlation panels separately.
#' @param xlab,ylab Optional axis labels.
#' @param point_alpha Point fill transparency.
#' @param variance_col,correlation_col Point fill colors for the two panels.
#' @param variance_line_col,correlation_line_col Curve colors for the two
#'   panels.
#' @param ribbon_col Ribbon color.
#' @param pch,lwd Base plotting symbol and line width.
#' @param ... Additional graphical parameters passed to [graphics::plot()].
#'
#' @return Invisibly returns `x`.
#' @export
plot.regime_integration_relationships <- function(x,
                                                  panel = c(
                                                    "both",
                                                    "variance",
                                                    "correlation"
                                                  ),
                                                  main = NULL,
                                                  xlab = NULL,
                                                  ylab = "Log inferred regime rate",
                                                  point_alpha = 0.45,
                                                  variance_col = "#4C78A8",
                                                  correlation_col = "#D55E00",
                                                  variance_line_col = "#1B4F72",
                                                  correlation_line_col = "#8A3A00",
                                                  ribbon_col = "grey70",
                                                  pch = 21,
                                                  lwd = 2,
                                                  ...) {
  panel <- match.arg(panel)
  if (identical(panel, "both")) {
    old_par <- graphics::par(mfrow = c(1, 2))
    on.exit(graphics::par(old_par), add = TRUE)
    main <- .regime_integration_plot_main(main, c("Rate vs. variance", "Rate vs. integration"))
    .regime_integration_plot_panel(
      x = x,
      panel = "variance",
      main = main[[1L]],
      xlab = xlab,
      ylab = ylab,
      point_alpha = point_alpha,
      point_col = variance_col,
      line_col = variance_line_col,
      ribbon_col = ribbon_col,
      pch = pch,
      lwd = lwd,
      ...
    )
    .regime_integration_plot_panel(
      x = x,
      panel = "correlation",
      main = main[[2L]],
      xlab = xlab,
      ylab = ylab,
      point_alpha = point_alpha,
      point_col = correlation_col,
      line_col = correlation_line_col,
      ribbon_col = ribbon_col,
      pch = pch,
      lwd = lwd,
      ...
    )
    return(invisible(x))
  }

  default_main <- if (identical(panel, "variance")) {
    "Rate vs. variance"
  } else {
    "Rate vs. integration"
  }
  .regime_integration_plot_panel(
    x = x,
    panel = panel,
    main = .regime_integration_plot_main(main, default_main)[[1L]],
    xlab = xlab,
    ylab = ylab,
    point_alpha = point_alpha,
    point_col = if (identical(panel, "variance")) variance_col else correlation_col,
    line_col = if (identical(panel, "variance")) variance_line_col else correlation_line_col,
    ribbon_col = ribbon_col,
    pch = pch,
    lwd = lwd,
    ...
  )
  invisible(x)
}

#' Fisher's Z Transformation
#'
#' Transform correlation coefficients with `atanh(r)` while rejecting values
#' outside the open interval `(-1, 1)`.
#'
#' @param correlations Numeric correlation vector.
#'
#' @return Numeric vector of transformed correlations.
#' @export
fisher_z_transform <- function(correlations) {
  if (any(correlations <= -1 | correlations >= 1, na.rm = TRUE)) {
    stop("Fisher Z is undefined for correlations with abs(r) >= 1.", call. = FALSE)
  }
  atanh(correlations)
}

#' Convert Regime Correlation PCA Results to Data Frames
#'
#' @param x A `regime_correlation_pca` object from
#'   [regime_correlation_pca()].
#' @param row.names,optional Arguments required by the S3 generic; ignored.
#' @param component Which component to return: `"scores"`, `"loadings"`,
#'   `"variance"`, or `"diagnostics"`.
#' @param ... Unused.
#'
#' @return A data frame for the requested PCA component.
#' @export
as.data.frame.regime_correlation_pca <- function(x,
                                                 row.names = NULL,
                                                 optional = FALSE,
                                                 component = c(
                                                   "scores",
                                                   "loadings",
                                                   "variance",
                                                   "diagnostics"
                                                 ),
                                                 ...) {
  dots <- list(...)
  .regime_check_unused_dots(dots)
  component <- match.arg(component)

  if (identical(component, "scores")) {
    df <- data.frame(
      regime = x$regime_ids,
      x$scores,
      row.names = NULL,
      check.names = FALSE
    )
    df$tip_count <- unname(x$tip_counts[x$regime_ids])
    df$regime_age <- unname(x$regime_ages[x$regime_ids])
    return(df)
  }

  if (identical(component, "loadings")) {
    return(data.frame(
      loading = rownames(x$loadings),
      x$loadings,
      row.names = NULL,
      check.names = FALSE
    ))
  }

  if (identical(component, "variance")) {
    return(data.frame(
      PC = paste0("PC", seq_along(x$variance_explained)),
      eigenvalue = x$pca$sdev^2,
      variance_explained = x$variance_explained,
      cumulative_variance = cumsum(x$variance_explained),
      row.names = NULL,
      check.names = FALSE
    ))
  }

  data.frame(
    regime = x$regime_ids,
    tip_count = unname(x$tip_counts[x$regime_ids]),
    regime_age = unname(x$regime_ages[x$regime_ids]),
    row.names = NULL,
    check.names = FALSE
  )
}

#' Plot Regime Correlation PCA Results
#'
#' Draw compact PCA views from [regime_correlation_pca()]. The default
#' `"loadings"` view reconstructs upper-triangle loading vectors as symmetric
#' trait-by-trait heatmaps, matching the Supplementary Figure 4B-style
#' interpretation of correlation-structure axes.
#'
#' @param x A `regime_correlation_pca` object.
#' @param type Plot type: `"loadings"`, `"scores"`, or `"variance"`.
#' @param components Principal components to draw for loading or variance
#'   plots. Supply numeric indices or names such as `"PC1"`. By default, the
#'   first four available components, or all components when fewer exist, are
#'   drawn.
#' @param pc_x,pc_y Components for the score scatter plot. The default uses
#'   PC1 and PC2 when available, or PC1 for both axes when only one component
#'   exists.
#' @param palette Optional color vector for loading heatmaps.
#' @param main Optional plot title. For loading plots with multiple
#'   components, supply one title per component.
#' @param cluster Loading-heatmap clustering style: `"none"` preserves trait
#'   order, `"global"` reuses the first selected PC's clustered row/column
#'   order for all panels, and `"local"` reclusters each PC panel.
#' @param heatmap_engine Loading-heatmap renderer. `"auto"` uses
#'   `ComplexHeatmap` for clustered heatmaps when that optional package is
#'   installed and otherwise falls back to base graphics.
#' @param show_dendrogram,show_legend Logical controls used by the optional
#'   `ComplexHeatmap` renderer.
#' @param cex.axis Axis label size for loading heatmaps.
#' @param ... Additional graphical parameters passed to base plotting
#'   functions.
#'
#' @return Invisibly returns `x`.
#' @export
plot.regime_correlation_pca <- function(x,
                                        type = c("loadings", "scores", "variance"),
                                        components = 1:4,
                                        pc_x = 1L,
                                        pc_y = 2L,
                                        palette = NULL,
                                        main = NULL,
                                        cluster = c("none", "global", "local"),
                                        heatmap_engine = c("auto", "base", "ComplexHeatmap"),
                                        show_dendrogram = TRUE,
                                        show_legend = TRUE,
                                        cex.axis = 0.62,
                                        ...) {
  components_missing <- missing(components)
  pc_y_missing <- missing(pc_y)
  type <- match.arg(type)
  cluster <- match.arg(cluster)
  heatmap_engine <- match.arg(heatmap_engine)
  pc_names <- colnames(x$scores)
  if (components_missing) {
    components <- seq_len(min(4L, length(pc_names)))
  }
  if (pc_y_missing) {
    pc_y <- min(2L, length(pc_names))
  }

  if (identical(type, "loadings")) {
    component_idx <- .regime_pca_component_indices(components, pc_names)
    if (.regime_should_use_complex_heatmap(heatmap_engine, cluster)) {
      drawn <- .regime_pca_plot_loadings_complex_heatmap(
        x = x,
        component_idx = component_idx,
        main = main,
        palette = palette,
        cluster = cluster,
        show_dendrogram = show_dendrogram,
        show_legend = show_legend,
        cex.axis = cex.axis
      )
      if (drawn) {
        return(invisible(x))
      }
      if (identical(heatmap_engine, "ComplexHeatmap")) {
        stop(
          "`heatmap_engine = \"ComplexHeatmap\"` requires packages ",
          "`ComplexHeatmap` and `circlize`.",
          call. = FALSE
        )
      }
    }

    if (length(component_idx) > 1L) {
      old_par <- graphics::par(mfrow = grDevices::n2mfrow(length(component_idx)))
      on.exit(graphics::par(old_par), add = TRUE)
    }
    main <- .regime_integration_plot_main(
      main,
      paste0(pc_names[component_idx], " loadings")
    )
    global_order <- NULL
    if (identical(cluster, "global")) {
      first_mat <- .regime_pca_loading_matrix(
        loadings = x$loadings[, component_idx[[1L]]],
        trait_labels = x$trait_labels
      )
      global_order <- .regime_heatmap_orders(first_mat)
    }
    for (i in seq_along(component_idx)) {
      pc <- component_idx[[i]]
      mat <- .regime_pca_loading_matrix(
        loadings = x$loadings[, pc],
        trait_labels = x$trait_labels
      )
      orders <- if (identical(cluster, "global")) {
        global_order
      } else if (identical(cluster, "local")) {
        .regime_heatmap_orders(mat)
      } else {
        NULL
      }
      .regime_pca_plot_loading_matrix(
        mat = mat,
        main = main[[i]],
        palette = palette,
        orders = orders,
        cex.axis = cex.axis,
        ...
      )
    }
    return(invisible(x))
  }

  if (identical(type, "scores")) {
    pc_x <- .regime_pca_component_indices(pc_x, pc_names)[[1L]]
    pc_y <- .regime_pca_component_indices(pc_y, pc_names)[[1L]]
    graphics::plot(
      x$scores[, pc_x],
      x$scores[, pc_y],
      xlab = pc_names[[pc_x]],
      ylab = pc_names[[pc_y]],
      main = if (is.null(main)) "Regime correlation PCA scores" else main[[1L]],
      ...
    )
    return(invisible(x))
  }

  component_idx <- .regime_pca_component_indices(components, pc_names)
  heights <- 100 * x$variance_explained[component_idx]
  graphics::barplot(
    heights,
    names.arg = pc_names[component_idx],
    ylab = "Variance explained (%)",
    main = if (is.null(main)) "Regime correlation PCA variance" else main[[1L]],
    ...
  )
  invisible(x)
}

#' Convert Regime Module Diagnostics to Data Frames
#'
#' @param x A `regime_module_diagnostics` object.
#' @param row.names,optional Arguments required by the S3 generic; ignored.
#' @param component Which component to return: `"scores"`, `"correlations"`,
#'   `"plot_data"`, or `"comparisons"`.
#' @param ... Unused.
#'
#' @return A data frame for the requested component.
#' @export
as.data.frame.regime_module_diagnostics <- function(x,
                                                    row.names = NULL,
                                                    optional = FALSE,
                                                    component = c(
                                                      "scores",
                                                      "correlations",
                                                      "plot_data",
                                                      "comparisons"
                                                    ),
                                                    ...) {
  dots <- list(...)
  .regime_check_unused_dots(dots)
  component <- match.arg(component)
  if (identical(component, "scores")) {
    return(x$module_scores)
  }
  if (identical(component, "correlations")) {
    return(x$correlations)
  }
  if (identical(component, "plot_data")) {
    return(x$plot_data)
  }
  data.frame(
    module_comparison = names(x$comparisons),
    module_a = vapply(x$comparisons, `[[`, character(1), 1L),
    module_b = vapply(x$comparisons, `[[`, character(1), 2L),
    row.names = NULL,
    check.names = FALSE
  )
}

#' Plot Regime Module Diagnostics
#'
#' Draw scatterplots comparing selected PCA scores with selected module
#' summaries from [regime_module_diagnostics()].
#'
#' @param x A `regime_module_diagnostics` object.
#' @param pc Principal component name(s), such as `"PC1"`.
#' @param comparison Module comparison name(s).
#' @param main Optional panel title(s).
#' @param xlab,ylab Optional axis labels. Recycled across panels.
#' @param point_col,line_col Point fill and fitted-line colors.
#' @param ribbon_col Ribbon color for bootstrap confidence intervals.
#' @param point_alpha,pch,lwd Graphical controls.
#' @param n_boot Number of bootstrap replicates for confidence ribbons. Set to
#'   `0` to omit the ribbon and draw only the fitted line.
#' @param ci_level Confidence level for bootstrap ribbons.
#' @param seed Optional random seed for reproducible bootstrap curves.
#' @param ... Additional graphical parameters passed to [graphics::plot()].
#'
#' @return Invisibly returns `x`.
#' @export
plot.regime_module_diagnostics <- function(x,
                                           pc,
                                           comparison,
                                           main = NULL,
                                           xlab = NULL,
                                           ylab = NULL,
                                           point_col = "#4C78A8",
                                           line_col = "#1B4F72",
                                           ribbon_col = "grey70",
                                           point_alpha = 0.5,
                                           pch = 21,
                                           lwd = 2,
                                           n_boot = 1000,
                                           ci_level = 0.99,
                                           seed = 1L,
                                           ...) {
  if (!inherits(x, "regime_module_diagnostics")) {
    stop("`x` must be a `regime_module_diagnostics` object.", call. = FALSE)
  }
  if (missing(pc) || missing(comparison)) {
    stop("Supply both `pc` and `comparison`.", call. = FALSE)
  }
  selections <- .regime_recycle_plot_selection(pc, comparison, "pc", "comparison")
  pc <- selections$pc
  comparison <- selections$comparison
  missing_pc <- setdiff(pc, names(x$plot_data))
  if (length(missing_pc) > 0L) {
    stop("Unknown PC column: ", paste(missing_pc, collapse = ", "), call. = FALSE)
  }
  missing_comparison <- setdiff(comparison, names(x$module_scores))
  if (length(missing_comparison) > 0L) {
    stop(
      "Unknown module comparison: ",
      paste(missing_comparison, collapse = ", "),
      call. = FALSE
    )
  }
  n_boot <- .regime_validate_module_plot_bootstrap(n_boot)
  .regime_validate_module_plot_ci_level(ci_level)
  seed <- .regime_validate_module_plot_seed(seed)

  n_panels <- length(pc)
  if (!is.null(seed) && n_boot > 0L) {
    bootstrap_panels <- vapply(seq_len(n_panels), function(i) {
      sum(stats::complete.cases(
        x$plot_data[[pc[[i]]]],
        x$module_scores[[comparison[[i]]]]
      )) >= 3L
    }, logical(1))
    if (any(bootstrap_panels)) {
      max_offset <- as.double(max(which(bootstrap_panels)) - 1L)
      .regime_validate_module_plot_seed(as.double(seed) + max_offset)
    }
  }
  if (n_panels > 1L) {
    old_par <- graphics::par(mfrow = grDevices::n2mfrow(n_panels))
    on.exit(graphics::par(old_par), add = TRUE)
  }
  main <- .regime_integration_plot_main(
    main,
    paste(pc, "vs", .regime_pretty_module_name(comparison))
  )
  xlab <- .regime_recycle_label(xlab, paste0(pc, " score"), n_panels)
  ylab <- .regime_recycle_label(
    ylab,
    .regime_pretty_module_name(comparison),
    n_panels
  )
  point_col <- rep(point_col, length.out = n_panels)
  line_col <- rep(line_col, length.out = n_panels)
  ribbon_col <- rep(ribbon_col, length.out = n_panels)

  for (i in seq_len(n_panels)) {
    x_vals <- x$plot_data[[pc[[i]]]]
    y_vals <- x$module_scores[[comparison[[i]]]]
    fit_data <- stats::na.omit(data.frame(x_vals = x_vals, y_vals = y_vals))
    lm_fit <- NULL
    curve <- NULL
    if (nrow(fit_data) >= 2L) {
      lm_fit <- stats::lm(y_vals ~ x_vals, data = fit_data)
    }
    if (nrow(fit_data) >= 3L && n_boot > 0L) {
      panel_seed <- if (is.null(seed)) {
        NULL
      } else {
        as.double(seed) + as.double(i - 1L)
      }
      curve <- .regime_bootstrap_curve(
        data = fit_data,
        x_col = "x_vals",
        y_col = "y_vals",
        formula = y_vals ~ x_vals,
        n_boot = n_boot,
        ci_level = ci_level,
        seed = panel_seed
      )
    }

    plot_args <- list(...)
    if (is.null(plot_args$type)) {
      plot_args$type <- "n"
    }
    if (is.null(plot_args$ylim)) {
      if (!is.null(curve)) {
        plot_args$ylim <- range(
          fit_data$y_vals,
          curve$ymin,
          curve$ymax,
          na.rm = TRUE
        )
      } else if (!any(is.finite(y_vals))) {
        plot_args$ylim <- c(-1, 1)
      }
    }
    do.call(
      graphics::plot,
      c(
        list(
          x = x_vals,
          y = y_vals,
          xlab = xlab[[i]],
          ylab = ylab[[i]],
          main = main[[i]]
        ),
        plot_args
      )
    )
    if (!is.null(curve)) {
      graphics::polygon(
        c(curve$x, rev(curve$x)),
        c(curve$ymin, rev(curve$ymax)),
        col = grDevices::adjustcolor(ribbon_col[[i]], alpha.f = 0.45),
        border = NA
      )
      graphics::lines(curve$x, curve$y, col = line_col[[i]], lwd = lwd)
    } else if (!is.null(lm_fit)) {
      x_seq <- seq(min(fit_data$x_vals), max(fit_data$x_vals), length.out = 100)
      graphics::lines(
        x_seq,
        stats::predict(lm_fit, newdata = data.frame(x_vals = x_seq)),
        col = line_col[[i]],
        lwd = lwd
      )
    }
    graphics::points(
      x_vals,
      y_vals,
      pch = pch,
      bg = grDevices::adjustcolor(point_col[[i]], alpha.f = point_alpha),
      col = "white"
    )
    r <- x$correlations$correlation[
      x$correlations$PC == pc[[i]] &
        x$correlations$module_comparison == comparison[[i]]
    ]
    if (length(r) == 1L && is.finite(r)) {
      graphics::legend(
        "topleft",
        legend = paste0("r = ", round(r, 3)),
        bty = "n"
      )
    }
  }
  invisible(x)
}

.regime_validate_module_plot_bootstrap <- function(n_boot) {
  if (!is.numeric(n_boot) || length(n_boot) != 1L ||
      !is.finite(n_boot) || n_boot < 0 ||
      n_boot != floor(n_boot) || n_boot > .Machine$integer.max) {
    stop(
      "`n_boot` must be a non-negative integer (a whole number) in R integer range.",
      call. = FALSE
    )
  }
  as.integer(n_boot)
}

.regime_validate_module_plot_ci_level <- function(ci_level) {
  if (!is.numeric(ci_level) || length(ci_level) != 1L ||
      !is.finite(ci_level) || ci_level <= 0 || ci_level >= 1) {
    stop("`ci_level` must be between 0 and 1.", call. = FALSE)
  }
  invisible(ci_level)
}

.regime_validate_resid_sd_threshold <- function(threshold, name) {
  if (!is.numeric(threshold) || length(threshold) != 1L ||
      !is.finite(threshold) || threshold < 0) {
    stop("`", name, "` must be a non-negative finite number.", call. = FALSE)
  }
  as.numeric(threshold)
}

.regime_validate_module_plot_seed <- function(seed) {
  if (is.null(seed)) {
    return(NULL)
  }
  if (!is.numeric(seed) || length(seed) != 1L || !is.finite(seed) ||
      seed != floor(seed) ||
      seed < -.Machine$integer.max || seed > .Machine$integer.max) {
    stop(
      "`seed` must be `NULL` or a whole number in R integer range.",
      call. = FALSE
    )
  }
  as.integer(seed)
}

#' Convert Regime Integration Relationship Results to Data Frames
#'
#' @param x A `regime_integration_relationships` object.
#' @param row.names,optional Arguments required by the S3 generic; ignored.
#' @param component Which component to return.
#' @param ... Unused.
#'
#' @return A data frame for the requested component.
#' @export
as.data.frame.regime_integration_relationships <- function(x,
                                                           row.names = NULL,
                                                           optional = FALSE,
                                                           component = c(
                                                             "combined",
                                                             "variance_points",
                                                             "correlation_points",
                                                             "variance_removed",
                                                             "correlation_removed",
                                                             "variance_curve",
                                                             "correlation_curve",
                                                             "models"
                                                           ),
                                                           ...) {
  dots <- list(...)
  .regime_check_unused_dots(dots)
  component <- match.arg(component)
  if (identical(component, "models")) {
    return(.regime_integration_model_summary(x))
  }
  x[[component]]
}

.regime_resolve_tree <- function(x = NULL, tree = NULL) {
  if (!is.null(tree)) {
    out <- tree
  } else if (inherits(x, "phylo")) {
    out <- x
  } else if (is.list(x) && !is.null(x$tree_no_uncertainty_untransformed)) {
    out <- x$tree_no_uncertainty_untransformed
  } else {
    stop(
      "Supply a SIMMAP-style tree via `tree` or a `bifrost_search`/tree via `x`.",
      call. = FALSE
    )
  }

  if (!inherits(out, "phylo") || is.null(out$maps)) {
    stop("Regime integration requires a SIMMAP-style `phylo` tree with `$maps`.", call. = FALSE)
  }
  .regime_validate_identifiers(out$tip.label, "tree tip labels")
  map_state_names <- unlist(lapply(out$maps, names), use.names = FALSE)
  if (length(map_state_names) != sum(lengths(out$maps))) {
    stop(
      "`mapped regime-state identifiers` must contain non-empty identifiers.",
      call. = FALSE
    )
  }
  .regime_validate_identifiers(
    map_state_names,
    "mapped regime-state identifiers",
    require_unique = FALSE
  )
  out
}

.regime_validate_trait_data <- function(trait_data, tree) {
  if (!is.matrix(trait_data) && !is.data.frame(trait_data)) {
    stop("`trait_data` must be a matrix or data frame.", call. = FALSE)
  }
  if (is.null(rownames(trait_data))) {
    stop("`trait_data` must have row names matching tree tip labels.", call. = FALSE)
  }
  .regime_validate_identifiers(rownames(trait_data), "trait-data row names")
  if (!is.null(colnames(trait_data))) {
    .regime_validate_identifiers(
      colnames(trait_data),
      "trait-data column names"
    )
  }
  missing_tips <- setdiff(tree$tip.label, rownames(trait_data))
  if (length(missing_tips) > 0L) {
    stop(
      "`trait_data` is missing rows for tree tips: ",
      paste(utils::head(missing_tips, 5L), collapse = ", "),
      if (length(missing_tips) > 5L) ", ..." else "",
      call. = FALSE
    )
  }
  trait_data
}

.regime_validate_identifiers <- function(
    identifiers,
    label,
    require_unique = TRUE) {
  if (is.null(identifiers)) {
    stop("`", label, "` must contain non-empty identifiers.", call. = FALSE)
  }
  identifiers <- as.character(identifiers)
  invalid <- is.na(identifiers) | !nzchar(trimws(identifiers))
  if (any(invalid)) {
    stop("`", label, "` must contain non-empty identifiers.", call. = FALSE)
  }
  if (isTRUE(require_unique)) {
    duplicated_identifiers <- unique(identifiers[
      duplicated(identifiers) | duplicated(identifiers, fromLast = TRUE)
    ])
    if (length(duplicated_identifiers) > 0L) {
      stop(
        "`", label, "` contains duplicated identifier(s): ",
        paste(duplicated_identifiers, collapse = ", "),
        ".",
        call. = FALSE
      )
    }
  }
  invisible(identifiers)
}

.regime_validate_min_tips <- function(min_tips) {
  if (!is.numeric(min_tips) || length(min_tips) != 1L ||
      !is.finite(min_tips) || min_tips < 1 ||
      min_tips != floor(min_tips) || min_tips > .Machine$integer.max) {
    stop(
      "`min_tips` must be a positive integer (a whole number) in R integer range.",
      call. = FALSE
    )
  }
  as.integer(min_tips)
}

.regime_validate_cores <- function(cores) {
  if (!is.numeric(cores) || length(cores) != 1L ||
      !is.finite(cores) || cores < 1 ||
      cores != floor(cores) || cores > .Machine$integer.max) {
    stop(
      "`cores` must be a positive integer (a whole number) in R integer range.",
      call. = FALSE
    )
  }
  as.integer(cores)
}

.regime_validate_logical <- function(x, name) {
  if (!is.logical(x) || length(x) != 1L || is.na(x)) {
    stop("`", name, "` must be TRUE or FALSE.", call. = FALSE)
  }
  x
}

.regime_normalize_run_names <- function(run_names, n_runs, name) {
  if (is.null(run_names)) {
    run_names <- rep(NA_character_, n_runs)
  }
  missing_names <- is.na(run_names) | !nzchar(trimws(run_names))
  run_names[missing_names] <- paste0("run", which(missing_names))
  duplicated_run_names <- unique(run_names[duplicated(run_names)])
  if (length(duplicated_run_names) > 0L) {
    stop(
      "`", name, "` contains duplicated run name(s): ",
      paste(duplicated_run_names, collapse = ", "),
      ".",
      call. = FALSE
    )
  }
  run_names
}

.regime_validate_run_list <- function(x, name) {
  if (inherits(x, "phylo") ||
      inherits(x, "bifrost_search") ||
      inherits(x, "regime_covariances")) {
    stop(
      "`",
      name,
      "` must be a list of runs. Use the singular function for one run.",
      call. = FALSE
    )
  }
  if (!is.list(x) || length(x) == 0L) {
    stop("`", name, "` must be a non-empty list.", call. = FALSE)
  }
  run_names <- .regime_normalize_run_names(names(x), length(x), name)
  names(x) <- run_names
  x
}

.regime_match_run_argument <- function(value, run_names, argument) {
  out <- stats::setNames(vector("list", length(run_names)), run_names)
  if (is.null(value)) {
    return(out)
  }

  if (is.list(value) &&
      !inherits(value, "phylo") &&
      !inherits(value, "bifrost_search") &&
      !inherits(value, "regime_covariances")) {
    value_names <- names(value)
    if (!is.null(value_names) && all(run_names %in% value_names)) {
      return(value[run_names])
    }
    if (length(value) == length(run_names)) {
      names(value) <- run_names
      return(value)
    }
    stop(
      "`",
      argument,
      "` must be named with the same run names as `x`, have one value per run, ",
      "or be a single value reused for every run.",
      call. = FALSE
    )
  }

  for (run_name in run_names) {
    out[[run_name]] <- value
  }
  out
}

.regime_validate_corr_threshold <- function(corr_threshold) {
  if (!is.numeric(corr_threshold) || length(corr_threshold) != 1L ||
      !is.finite(corr_threshold) || corr_threshold <= 0 ||
      corr_threshold >= 1) {
    stop("`corr_threshold` must be a finite number between 0 and 1.", call. = FALSE)
  }
  corr_threshold
}

.regime_tip_states <- function(tree) {
  states <- phytools::getStates(tree, type = "tips")
  states <- as.character(states)
  .regime_validate_identifiers(
    states,
    "mapped tip-state identifiers",
    require_unique = FALSE
  )
  if (is.null(names(states))) {
    names(states) <- tree$tip.label
  }
  states
}

.regime_subtree_age <- function(tree, tips) {
  if (length(tips) < 2L) {
    return(NA_real_)
  }
  subtree <- tryCatch(ape::keep.tip(tree, tips), error = function(e) NULL)
  if (is.null(subtree) || is.null(subtree$edge.length)) {
    return(NA_real_)
  }
  max(phytools::nodeHeights(subtree), na.rm = TRUE)
}

.regime_formula_with_data <- function(formula, trait_data) {
  formula_obj <- if (inherits(formula, "formula")) {
    formula
  } else {
    stats::as.formula(formula)
  }
  parent_env <- environment(formula_obj)
  if (is.null(parent_env)) {
    parent_env <- parent.frame()
  }
  eval_env <- new.env(parent = parent_env)
  assign("trait_data", trait_data, envir = eval_env)

  col_names <- colnames(trait_data)
  if (!is.null(col_names)) {
    for (name in setdiff(col_names, "trait_data")) {
      assign(name, trait_data[, name], envir = eval_env)
    }
  }
  environment(formula_obj) <- eval_env
  formula_obj
}

.regime_extract_covariance <- function(fit) {
  if (is.list(fit$sigma) && !is.null(fit$sigma$Pinv) && is.matrix(fit$sigma$Pinv)) {
    return(fit$sigma$Pinv)
  }
  if (!is.null(fit$sigma) && is.matrix(fit$sigma)) {
    return(fit$sigma)
  }
  NULL
}

.regime_extract_rates <- function(x) {
  if (is.null(x)) {
    return(NULL)
  }
  rates <- NULL
  if (is.list(x) && !is.null(x$model_no_uncertainty$param)) {
    rates <- x$model_no_uncertainty$param
  } else if (is.list(x) && !is.null(x$param)) {
    rates <- x$param
  }
  if (is.null(rates)) {
    return(NULL)
  }
  if (!is.numeric(rates) || is.null(names(rates))) {
    return(NULL)
  }
  .regime_validate_identifiers(names(rates), "rate names")
  rates
}

.regime_extract_covariance_list <- function(x) {
  if (inherits(x, "regime_covariances")) {
    covariances <- x$covariances
    if (!is.list(covariances) || length(covariances) == 0L) {
      stop("`x$covariances` must be a non-empty named list.", call. = FALSE)
    }
    .regime_validate_identifiers(
      names(covariances),
      "covariance-list regime names"
    )
    if (!is.data.frame(x$status) || !"regime" %in% names(x$status)) {
      stop("`x$status` must be a data frame with a `regime` column.", call. = FALSE)
    }
    if (nrow(x$status) > 0L) {
      .regime_validate_identifiers(
        x$status$regime,
        "status-table regime IDs"
      )
    }
    invalid_messages <- stats::setNames(
      rep(NA_character_, length(covariances)),
      names(covariances)
    )
    for (regime in names(covariances)) {
      mat <- covariances[[regime]]
      if (!is.matrix(mat)) {
        next
      }
      validation_error <- tryCatch(
        {
          .regime_validate_covariance_matrix(mat, regime)
          NULL
        },
        error = function(e) e
      )
      if (inherits(validation_error, "error")) {
        covariances[regime] <- list(NULL)
        invalid_messages[[regime]] <- conditionMessage(validation_error)
      }
    }
    return(list(
      covariances = covariances,
      status = x$status,
      rates = x$rates,
      invalid_messages = invalid_messages
    ))
  }

  if (!is.list(x) || length(x) == 0L) {
    stop("`x` must be a `regime_covariances` object or a named list of matrices.", call. = FALSE)
  }
  .regime_validate_identifiers(names(x), "covariance-list regime names")
  status <- data.frame(
    regime = names(x),
    tip_count = rep(NA_real_, length(x)),
    regime_age = rep(NA_real_, length(x)),
    status = ifelse(vapply(x, is.matrix, logical(1)), "ok", "failed"),
    message = NA_character_,
    row.names = NULL,
    check.names = FALSE
  )
  list(
    covariances = x,
    status = status,
    rates = NULL,
    invalid_messages = stats::setNames(rep(NA_character_, length(x)), names(x))
  )
}

.regime_warn_if_search_vcvs <- function(x, search) {
  if (is.null(search) ||
      is.null(search$VCVs) ||
      !identical(x, search$VCVs)) {
    return(invisible(FALSE))
  }

  warning(
    "`x` appears to be proportional `search$VCVs` from a scalar `bifrost` ",
    "search. These matrices summarize the joint model and should not be used ",
    "as independent post-hoc covariance refits for integration or correlation ",
    "reconfiguration questions.",
    call. = FALSE
  )
  invisible(TRUE)
}

.regime_match_rates <- function(rates, regimes) {
  if (is.null(rates)) {
    return(stats::setNames(rep(NA_real_, length(regimes)), regimes))
  }
  if (!is.numeric(rates) || is.null(names(rates))) {
    stop("`rates` must be a named numeric vector.", call. = FALSE)
  }
  .regime_validate_identifiers(names(rates), "rate names")
  out <- stats::setNames(rep(NA_real_, length(regimes)), regimes)
  matched <- intersect(regimes, names(rates))
  out[matched] <- rates[matched]
  out
}

.regime_match_optional <- function(values, regimes) {
  if (is.null(values)) {
    return(stats::setNames(rep(NA_real_, length(regimes)), regimes))
  }
  if (is.null(names(values))) {
    if (length(values) != length(regimes)) {
      stop("Unnamed optional vectors must have one value per regime.", call. = FALSE)
    }
    return(stats::setNames(as.numeric(values), regimes))
  }
  .regime_validate_identifiers(
    names(values),
    "optional-vector regime names"
  )
  out <- stats::setNames(rep(NA_real_, length(regimes)), regimes)
  matched <- intersect(regimes, names(values))
  out[matched] <- as.numeric(values[matched])
  out
}

.regime_tree_stats <- function(tree) {
  tree <- .regime_resolve_tree(tree = tree)
  tip_states <- .regime_tip_states(tree)
  tips_by_regime <- split(names(tip_states), unname(tip_states))
  regimes <- names(tips_by_regime)
  list(
    tip_counts = stats::setNames(
      vapply(tips_by_regime, length, integer(1)),
      regimes
    ),
    regime_ages = stats::setNames(
      vapply(tips_by_regime, .regime_subtree_age, numeric(1), tree = tree),
      regimes
    )
  )
}

.regime_validate_square_matrix <- function(mat, regime) {
  if (!is.matrix(mat) || !is.numeric(mat) || nrow(mat) != ncol(mat)) {
    stop("Regime `", regime, "` must contain a square numeric matrix.", call. = FALSE)
  }
  invisible(TRUE)
}

.regime_validate_covariance_matrix <- function(
    mat,
    regime,
    tolerance = sqrt(.Machine$double.eps)) {
  .regime_validate_square_matrix(mat, regime)
  if (any(!is.finite(mat))) {
    stop(
      "Regime `", regime, "` covariance entries must all be finite.",
      call. = FALSE
    )
  }
  symmetry_scale <- max(1, max(abs(mat)))
  if (max(abs(mat - t(mat))) > tolerance * symmetry_scale) {
    stop(
      "Regime `", regime, "` covariance matrix must be symmetric within ",
      "numeric tolerance.",
      call. = FALSE
    )
  }
  row_traits <- rownames(mat)
  col_traits <- colnames(mat)
  if (xor(is.null(row_traits), is.null(col_traits))) {
    stop(
      "Regime `", regime, "` must have matching row and column trait names.",
      call. = FALSE
    )
  }
  if (!is.null(row_traits)) {
    .regime_validate_identifiers(
      row_traits,
      paste0("Regime `", regime, "` row trait names")
    )
    .regime_validate_identifiers(
      col_traits,
      paste0("Regime `", regime, "` column trait names")
    )
    if (!identical(row_traits, col_traits)) {
      stop(
        "Regime `", regime, "` must have matching row and column trait names.",
        call. = FALSE
      )
    }
  }
  if (any(diag(mat) <= 0)) {
    stop(
      "Regime `", regime, "` must have strictly positive diagonal variances.",
      call. = FALSE
    )
  }
  symmetric_mat <- (mat + t(mat)) / 2
  eigenvalues <- eigen(
    symmetric_mat,
    symmetric = TRUE,
    only.values = TRUE
  )$values
  if (min(eigenvalues) < -tolerance * symmetry_scale) {
    stop(
      "Regime `", regime, "` covariance matrix must be positive semidefinite.",
      call. = FALSE
    )
  }
  invisible(TRUE)
}

.regime_align_matrix_traits <- function(mat, template_traits, regime) {
  row_traits <- rownames(mat)
  col_traits <- colnames(mat)
  if (is.null(row_traits) && is.null(col_traits)) {
    return(mat)
  }
  if (is.null(row_traits) || is.null(col_traits) || !identical(row_traits, col_traits)) {
    stop("Regime `", regime, "` must have matching row and column trait names.", call. = FALSE)
  }
  if (!setequal(row_traits, template_traits)) {
    stop("Regime `", regime, "` has trait names that do not match the first matrix.", call. = FALSE)
  }
  mat[template_traits, template_traits, drop = FALSE]
}

.regime_covariance_list_is_scalar_proportional <- function(matrices, tolerance = sqrt(.Machine$double.eps)) {
  template <- matrices[[1L]]
  template_values <- as.numeric(template)
  template_norm <- sqrt(sum(template_values^2))
  if (!is.finite(template_norm) || template_norm <= tolerance) {
    return(FALSE)
  }

  vapply(matrices[-1L], function(mat) {
    values <- as.numeric(mat)
    scalar <- sum(values * template_values) / sum(template_values^2)
    is.finite(scalar) &&
      max(abs(values - scalar * template_values), na.rm = TRUE) <=
        tolerance * max(1, sqrt(sum(values^2)))
  }, logical(1)) |>
    all()
}

.regime_fisher_z <- function(r, boundary) {
  if (is.na(r)) {
    return(NA_real_)
  }
  if (abs(r) >= 1) {
    if (identical(boundary, "error")) {
      stop("Fisher Z is undefined for correlations with abs(r) >= 1.", call. = FALSE)
    }
    return(NA_real_)
  }
  atanh(r)
}

.regime_standardize_summary_data <- function(summary_data) {
  if (!is.data.frame(summary_data)) {
    stop("`summary_data` must be a data frame.", call. = FALSE)
  }

  if (all(c("rate", "vars", "corrs") %in% names(summary_data))) {
    out <- data.frame(
      regime = if ("State" %in% names(summary_data)) as.character(summary_data$State) else rownames(summary_data),
      rate = summary_data$rate,
      vars = summary_data$vars,
      corrs = summary_data$corrs,
      row.names = NULL,
      check.names = FALSE
    )
    .regime_validate_identifiers(out$regime, "summary-data regime IDs")
    return(.regime_append_optional_summary_columns(out, summary_data))
  }

  required <- c("regime", "rate", "mean_variance", "mean_abs_correlation")
  missing <- setdiff(required, names(summary_data))
  if (length(missing) > 0L) {
    stop(
      "`summary_data` is missing required columns: ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }
  out <- data.frame(
    regime = as.character(summary_data$regime),
    rate = summary_data$rate,
    vars = summary_data$mean_variance,
    corrs = summary_data$mean_abs_correlation,
    row.names = NULL,
    check.names = FALSE
  )
  .regime_validate_identifiers(out$regime, "summary-data regime IDs")
  .regime_append_optional_summary_columns(out, summary_data)
}

.regime_validate_relationship_values <- function(data) {
  validate_positive <- function(values, column) {
    if (!is.numeric(values)) {
      stop("`", column, "` must be numeric.", call. = FALSE)
    }
    invalid <- is.nan(values) |
      (!is.na(values) & (!is.finite(values) | values <= 0))
    if (any(invalid)) {
      stop(
        "`", column,
        "` values must be finite and strictly positive when present; ",
        "invalid regime(s): ",
        paste(data$regime[invalid], collapse = ", "),
        ".",
        call. = FALSE
      )
    }
  }

  validate_positive(data$rate, "rate")
  validate_positive(data$vars, "vars")

  if (!is.numeric(data$corrs)) {
    stop("`corrs` must be numeric.", call. = FALSE)
  }
  invalid_correlation <- is.nan(data$corrs) |
    (!is.na(data$corrs) &
      (!is.finite(data$corrs) | abs(data$corrs) > 1))
  if (any(invalid_correlation)) {
    stop(
      "`corrs` values must be finite numbers between -1 and 1 when present; ",
      "invalid regime(s): ",
      paste(data$regime[invalid_correlation], collapse = ", "),
      ".",
      call. = FALSE
    )
  }

  invisible(data)
}

.regime_append_optional_summary_columns <- function(out, summary_data) {
  optional <- c("tip_count", "regime_age", "status", "message")
  for (column in intersect(optional, names(summary_data))) {
    out[[column]] <- summary_data[[column]]
  }
  out
}

.regime_filter_summary_min_tips <- function(data, min_tips, caller) {
  if (is.null(min_tips)) {
    return(data)
  }
  min_tips <- .regime_validate_min_tips(min_tips)
  if (!"tip_count" %in% names(data)) {
    stop(
      "`min_tips` was supplied to `",
      caller,
      "()` but `summary_data` does not contain a `tip_count` column.",
      call. = FALSE
    )
  }
  data[!is.na(data$tip_count) & data$tip_count >= min_tips, , drop = FALSE]
}

.regime_bind_summary_runs <- function(vars_cors) {
  if (is.data.frame(vars_cors)) {
    out <- .regime_standardize_summary_data(vars_cors)
    out$run <- if ("run" %in% names(vars_cors)) vars_cors$run else "run1"
    return(out[, c("run", setdiff(names(out), "run")), drop = FALSE])
  }

  if (!is.list(vars_cors) || length(vars_cors) == 0L) {
    stop("`vars_cors` must be a data frame or a non-empty list of data frames.", call. = FALSE)
  }
  run_names <- .regime_normalize_run_names(
    names(vars_cors),
    length(vars_cors),
    "summaries"
  )
  out <- lapply(seq_along(vars_cors), function(i) {
    df <- .regime_standardize_summary_data(vars_cors[[i]])
    df$run <- run_names[[i]]
    df[, c("run", setdiff(names(df), "run")), drop = FALSE]
  })
  do.call(rbind, out)
}

.regime_bootstrap_curve <- function(data,
                                    x_col,
                                    y_col,
                                    formula,
                                    n_boot,
                                    ci_level,
                                    seed) {
  if (nrow(data) < 3L) {
    stop("At least three rows are required to build bootstrap plot curves.", call. = FALSE)
  }
  if (!is.numeric(n_boot) || length(n_boot) != 1L ||
      !is.finite(n_boot) || n_boot < 1 ||
      n_boot != floor(n_boot) || n_boot > .Machine$integer.max) {
    stop(
      "`n_boot` must be a positive integer (a whole number) in R integer range.",
      call. = FALSE
    )
  }
  n_boot <- as.integer(n_boot)
  seed <- .regime_validate_module_plot_seed(seed)
  if (!is.numeric(ci_level) || length(ci_level) != 1L ||
      !is.finite(ci_level) || ci_level <= 0 || ci_level >= 1) {
    stop("`ci_level` must be between 0 and 1.", call. = FALSE)
  }

  x_seq <- seq(min(data[[x_col]], na.rm = TRUE), max(data[[x_col]], na.rm = TRUE), length.out = 100)
  newdata <- stats::setNames(data.frame(x_seq), x_col)

  if (!is.null(seed)) {
    old_seed <- if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      get(".Random.seed", envir = .GlobalEnv)
    } else {
      NULL
    }
    on.exit({
      if (is.null(old_seed)) {
        rm(".Random.seed", envir = .GlobalEnv)
      } else {
        assign(".Random.seed", old_seed, envir = .GlobalEnv)
      }
    }, add = TRUE)
    set.seed(seed)
  }

  preds <- replicate(n_boot, {
    boot_sample <- data[sample.int(nrow(data), replace = TRUE), , drop = FALSE]
    stats::predict(stats::lm(formula, data = boot_sample), newdata = newdata)
  })
  probs <- c((1 - ci_level) / 2, 1 - (1 - ci_level) / 2)
  ci <- apply(preds, 1L, stats::quantile, probs = probs, na.rm = TRUE)
  data.frame(
    x = x_seq,
    y = rowMeans(preds, na.rm = TRUE),
    ymin = ci[1L, ],
    ymax = ci[2L, ],
    row.names = NULL
  )
}

.regime_integration_model_summary <- function(x) {
  variance_summary <- summary(x$variance_lm)
  correlation_summary <- summary(x$correlation_lm)
  data.frame(
    panel = c("variance", "correlation"),
    n = c(nrow(x$variance_points), nrow(x$correlation_points)),
    intercept = c(
      unname(stats::coef(x$variance_lm)[[1L]]),
      unname(stats::coef(x$correlation_lm)[[1L]])
    ),
    slope = c(
      unname(stats::coef(x$variance_lm)[[2L]]),
      unname(stats::coef(x$correlation_lm)[[2L]])
    ),
    r_squared = c(
      variance_summary$r.squared,
      correlation_summary$r.squared
    ),
    p_value = c(
      stats::coef(variance_summary)[2L, 4L],
      stats::coef(correlation_summary)[2L, 4L]
    ),
    row.names = NULL,
    check.names = FALSE
  )
}

.regime_integration_plot_main <- function(main, default) {
  if (is.null(main)) {
    return(as.list(default))
  }
  if (length(main) == 1L && length(default) > 1L) {
    return(as.list(rep(main, length(default))))
  }
  if (length(main) < length(default)) {
    stop("`main` must have one value per requested panel.", call. = FALSE)
  }
  as.list(main[seq_along(default)])
}

.regime_integration_plot_panel <- function(x,
                                           panel,
                                           main,
                                           xlab,
                                           ylab,
                                           point_alpha,
                                           point_col,
                                           line_col,
                                           ribbon_col,
                                           pch,
                                           lwd,
                                           ...) {
  if (identical(panel, "variance")) {
    points <- x$variance_points
    curve <- x$variance_curve
    x_col <- "log_vars"
    if (is.null(xlab)) {
      xlab <- "Log mean post-hoc variance"
    }
  } else {
    points <- x$correlation_points
    curve <- x$correlation_curve
    x_col <- "fisher_z_corr"
    if (is.null(xlab)) {
      xlab <- "Fisher-Z mean absolute correlation"
    }
  }

  graphics::plot(
    points[[x_col]],
    points$log_rate,
    type = "n",
    xlab = xlab,
    ylab = ylab,
    main = main,
    ...
  )
  graphics::polygon(
    c(curve[[x_col]], rev(curve[[x_col]])),
    c(curve$ymin, rev(curve$ymax)),
    col = grDevices::adjustcolor(ribbon_col, alpha.f = 0.45),
    border = NA
  )
  graphics::lines(curve[[x_col]], curve$y, col = line_col, lwd = lwd)
  graphics::points(
    points[[x_col]],
    points$log_rate,
    pch = pch,
    bg = grDevices::adjustcolor(point_col, alpha.f = point_alpha),
    col = "white"
  )
}

.regime_pca_component_indices <- function(components, pc_names) {
  if (length(components) < 1L) {
    stop("`components` must identify at least one principal component.", call. = FALSE)
  }
  if (is.character(components)) {
    idx <- match(components, pc_names)
    if (anyNA(idx)) {
      stop("Unknown PCA component: ", paste(components[is.na(idx)], collapse = ", "), call. = FALSE)
    }
    return(as.integer(idx))
  }
  if (!is.numeric(components) || is.complex(components) ||
      any(!is.finite(components))) {
    stop("`components` must identify valid principal components.", call. = FALSE)
  }
  if (any(components != floor(components))) {
    stop("Numeric `components` indices must be whole-number values.", call. = FALSE)
  }
  if (any(components < 1 | components > length(pc_names))) {
    stop("`components` must identify valid principal components.", call. = FALSE)
  }
  as.integer(components)
}

.regime_pca_loading_matrix <- function(loadings, trait_labels) {
  p <- length(trait_labels)
  mat <- matrix(
    0,
    nrow = p,
    ncol = p,
    dimnames = list(trait_labels, trait_labels)
  )
  idx <- which(upper.tri(mat, diag = FALSE), arr.ind = TRUE)
  mat[idx] <- loadings
  mat <- mat + t(mat)
  diag(mat) <- 0
  mat
}

.regime_vector_to_symmetric_matrix <- function(values, trait_labels) {
  p <- length(trait_labels)
  mat <- matrix(
    0,
    nrow = p,
    ncol = p,
    dimnames = list(trait_labels, trait_labels)
  )
  idx <- which(upper.tri(mat, diag = FALSE), arr.ind = TRUE)
  mat[idx] <- values
  mat <- mat + t(mat)
  diag(mat) <- NA_real_
  mat
}

.regime_pca_plot_loading_matrix <- function(mat,
                                            main,
                                            palette,
                                            orders,
                                            cex.axis,
                                            ...) {
  if (is.null(palette)) {
    palette <- .regime_pca_loading_palette()
  }
  if (!is.null(orders)) {
    mat <- mat[orders$row_order, orders$column_order, drop = FALSE]
  }
  lim <- max(abs(mat), na.rm = TRUE)
  if (!is.finite(lim) || lim == 0) {
    lim <- 1
  }
  breaks <- seq(-lim, lim, length.out = length(palette) + 1L)
  z <- t(mat[nrow(mat):1L, , drop = FALSE])
  graphics::image(
    seq_len(ncol(mat)),
    seq_len(nrow(mat)),
    z,
    col = palette,
    breaks = breaks,
    axes = FALSE,
    xlab = "",
    ylab = "",
    main = main,
    ...
  )
  graphics::axis(
    1,
    at = seq_len(ncol(mat)),
    labels = colnames(mat),
    las = 2,
    cex.axis = cex.axis
  )
  graphics::axis(
    2,
    at = seq_len(nrow(mat)),
    labels = rev(rownames(mat)),
    las = 2,
    cex.axis = cex.axis
  )
  graphics::box()
}

.regime_pca_loading_palette <- function(n = 101L) {
  if (requireNamespace("RColorBrewer", quietly = TRUE)) {
    return(grDevices::colorRampPalette(
      rev(RColorBrewer::brewer.pal(11, "RdYlBu"))
    )(n))
  }
  grDevices::colorRampPalette(c("#3B4CC0", "#F7F7F7", "#B40426"))(n)
}

.regime_should_use_complex_heatmap <- function(heatmap_engine, cluster) {
  if (identical(heatmap_engine, "base")) {
    return(FALSE)
  }
  if (identical(heatmap_engine, "ComplexHeatmap")) {
    return(TRUE)
  }
  !identical(cluster, "none") &&
    requireNamespace("ComplexHeatmap", quietly = TRUE) &&
    requireNamespace("circlize", quietly = TRUE)
}

.regime_pca_plot_loadings_complex_heatmap <- function(x,
                                                      component_idx,
                                                      main,
                                                      palette,
                                                      cluster,
                                                      show_dendrogram,
                                                      show_legend,
                                                      cex.axis) {
  if (!requireNamespace("ComplexHeatmap", quietly = TRUE) ||
      !requireNamespace("circlize", quietly = TRUE)) {
    return(FALSE)
  }
  if (is.null(palette)) {
    palette <- .regime_pca_loading_palette()
  }

  pc_names <- colnames(x$scores)
  main <- .regime_integration_plot_main(
    main,
    paste0(pc_names[component_idx], " loadings")
  )
  mats <- lapply(component_idx, function(pc) {
    .regime_pca_loading_matrix(
      loadings = x$loadings[, pc],
      trait_labels = x$trait_labels
    )
  })
  lim <- max(abs(unlist(mats, use.names = FALSE)), na.rm = TRUE)
  if (!is.finite(lim) || lim == 0) {
    lim <- 1
  }
  col_fun <- circlize::colorRamp2(
    seq(-lim, lim, length.out = length(palette)),
    palette
  )
  global_orders <- if (identical(cluster, "global")) {
    .regime_heatmap_orders(mats[[1L]])
  } else {
    NULL
  }

  ht_list <- NULL
  for (i in seq_along(mats)) {
    cluster_rows <- FALSE
    cluster_columns <- FALSE
    if (identical(cluster, "local")) {
      cluster_rows <- TRUE
      cluster_columns <- TRUE
    } else if (identical(cluster, "global")) {
      cluster_rows <- global_orders$row_dend
      cluster_columns <- global_orders$column_dend
    }
    body_size <- .regime_heatmap_body_size(mats[[i]])

    ht <- ComplexHeatmap::Heatmap(
      mats[[i]],
      name = pc_names[[component_idx[[i]]]],
      col = col_fun,
      width = body_size$width,
      height = body_size$height,
      cluster_rows = cluster_rows,
      cluster_columns = cluster_columns,
      show_row_dend = show_dendrogram && !identical(cluster, "none"),
      show_column_dend = show_dendrogram && !identical(cluster, "none"),
      row_names_side = "right",
      column_names_side = "bottom",
      row_names_gp = grid::gpar(fontsize = 9 * cex.axis / 0.62),
      column_names_gp = grid::gpar(fontsize = 9 * cex.axis / 0.62),
      column_names_rot = 90,
      column_title = main[[i]],
      heatmap_legend_param = list(
        title = "Loading",
        at = c(-lim, 0, lim),
        labels = round(c(-lim, 0, lim), 2)
      ),
      show_heatmap_legend = show_legend && identical(i, 1L)
    )
    ht_list <- if (is.null(ht_list)) ht else ht_list + ht
  }
  ComplexHeatmap::draw(ht_list, merge_legend = TRUE)
  TRUE
}

.regime_heatmap_body_size <- function(mat, cell_size_mm = 4.5) {
  side_mm <- max(dim(mat)) * cell_size_mm
  list(
    width = grid::unit(side_mm, "mm"),
    height = grid::unit(side_mm, "mm")
  )
}

.regime_heatmap_orders <- function(mat) {
  row_hc <- if (nrow(mat) > 1L) {
    stats::hclust(stats::dist(mat))
  } else {
    NULL
  }
  column_hc <- if (ncol(mat) > 1L) {
    stats::hclust(stats::dist(t(mat)))
  } else {
    NULL
  }
  list(
    row_order = if (is.null(row_hc)) seq_len(nrow(mat)) else row_hc$order,
    column_order = if (is.null(column_hc)) seq_len(ncol(mat)) else column_hc$order,
    row_dend = if (is.null(row_hc)) FALSE else stats::as.dendrogram(row_hc),
    column_dend = if (is.null(column_hc)) FALSE else stats::as.dendrogram(column_hc)
  )
}

.regime_validate_modules <- function(modules, trait_labels) {
  if (!is.list(modules) || length(modules) == 0L) {
    stop("`modules` must be a non-empty named list.", call. = FALSE)
  }
  if (is.null(names(modules))) {
    stop("`modules` must be named.", call. = FALSE)
  }
  .regime_validate_identifiers(names(modules), "module names")
  modules <- lapply(modules, as.character)
  empty <- names(modules)[vapply(modules, length, integer(1)) == 0L]
  if (length(empty) > 0L) {
    stop("Empty module definition(s): ", paste(empty, collapse = ", "), call. = FALSE)
  }
  for (module_name in names(modules)) {
    duplicated_traits <- unique(
      modules[[module_name]][duplicated(modules[[module_name]])]
    )
    if (length(duplicated_traits) > 0L) {
      stop(
        "Module `", module_name, "` contains duplicated trait label(s): ",
        paste(duplicated_traits, collapse = ", "),
        ".",
        call. = FALSE
      )
    }
    .regime_validate_identifiers(
      modules[[module_name]],
      paste0("trait labels in module `", module_name, "`")
    )
  }
  module_traits <- unique(unlist(modules, use.names = FALSE))
  missing_traits <- setdiff(module_traits, trait_labels)
  if (length(missing_traits) > 0L) {
    stop(
      "`modules` contains trait labels not present in `pca`: ",
      paste(missing_traits, collapse = ", "),
      call. = FALSE
    )
  }
  modules
}

.regime_validate_module_comparisons <- function(comparisons, modules) {
  if (is.null(comparisons)) {
    comparisons <- .regime_default_module_comparisons(modules)
  }
  if (!is.list(comparisons) || length(comparisons) == 0L) {
    stop("`comparisons` must be a non-empty named list.", call. = FALSE)
  }
  comparison_names <- names(comparisons)
  generated_names <- vapply(
    comparisons,
    function(x) paste(as.character(x), collapse = "_vs_"),
    character(1)
  )
  if (is.null(comparison_names)) {
    comparison_names <- generated_names
  } else {
    missing_names <- is.na(comparison_names) |
      !nzchar(trimws(comparison_names))
    comparison_names[missing_names] <- generated_names[missing_names]
  }
  names(comparisons) <- comparison_names
  .regime_validate_identifiers(names(comparisons), "comparison names")
  comparisons <- lapply(comparisons, as.character)
  invalid_length <- names(comparisons)[
    vapply(comparisons, length, integer(1)) != 2L
  ]
  if (length(invalid_length) > 0L) {
    stop(
      "`comparisons` entries must name exactly two modules: ",
      paste(invalid_length, collapse = ", "),
      call. = FALSE
    )
  }
  missing_modules <- setdiff(
    unique(unlist(comparisons, use.names = FALSE)),
    names(modules)
  )
  if (length(missing_modules) > 0L) {
    stop(
      "`comparisons` names unknown module(s): ",
      paste(missing_modules, collapse = ", "),
      call. = FALSE
    )
  }
  comparisons
}

.regime_default_module_comparisons <- function(modules) {
  module_names <- names(modules)
  within <- lapply(module_names, function(module) c(module, module))
  names(within) <- paste0("within_", module_names)
  between <- if (length(module_names) >= 2L) {
    utils::combn(module_names, 2L, simplify = FALSE)
  } else {
    list()
  }
  names(between) <- vapply(
    between,
    function(x) paste(x, collapse = "_vs_"),
    character(1)
  )
  c(within, between)
}

.regime_module_comparison_score <- function(mat, module_a, module_b) {
  if (identical(module_a, module_b)) {
    if (length(module_a) < 2L) {
      return(NA_real_)
    }
    sub_mat <- mat[module_a, module_a, drop = FALSE]
    return(mean(sub_mat[upper.tri(sub_mat, diag = FALSE)], na.rm = TRUE))
  }
  mean(mat[module_a, module_b, drop = FALSE], na.rm = TRUE)
}

.regime_complete_correlation <- function(x, y) {
  complete <- stats::complete.cases(x, y)
  if (sum(complete) < 2L) {
    return(NA_real_)
  }
  x <- x[complete]
  y <- y[complete]
  if (all(x == x[[1L]]) || all(y == y[[1L]])) {
    return(NA_real_)
  }
  stats::cor(x, y)
}

.regime_recycle_plot_selection <- function(x, y, x_name, y_name) {
  if (length(x) == length(y)) {
    return(stats::setNames(list(x, y), c(x_name, y_name)))
  }
  if (length(x) == 1L) {
    x <- rep(x, length(y))
    return(stats::setNames(list(x, y), c(x_name, y_name)))
  }
  if (length(y) == 1L) {
    y <- rep(y, length(x))
    return(stats::setNames(list(x, y), c(x_name, y_name)))
  }
  stop(
    "`",
    x_name,
    "` and `",
    y_name,
    "` must have the same length, or one must have length one.",
    call. = FALSE
  )
}

.regime_recycle_label <- function(label, default, n) {
  if (is.null(label)) {
    return(as.list(default))
  }
  as.list(rep(label, length.out = n))
}

.regime_pretty_module_name <- function(x) {
  gsub("_", " ", x, fixed = TRUE)
}

.regime_resolve_trait_labels <- function(trait_labels, internal_traits) {
  if (is.null(trait_labels)) {
    return(internal_traits)
  }
  if (!is.null(names(trait_labels))) {
    missing <- setdiff(internal_traits, names(trait_labels))
    if (length(missing) > 0L) {
      stop(
        "`trait_labels` is missing labels for: ",
        paste(missing, collapse = ", "),
        call. = FALSE
      )
    }
    return(unname(trait_labels[internal_traits]))
  }
  if (length(trait_labels) != length(internal_traits)) {
    stop("Unnamed `trait_labels` must match the matrix dimension.", call. = FALSE)
  }
  as.character(trait_labels)
}

.regime_upper_pair_names <- function(trait_names) {
  idx <- which(
    upper.tri(matrix(NA, length(trait_names), length(trait_names)), diag = FALSE),
    arr.ind = TRUE
  )
  paste(trait_names[idx[, 1]], trait_names[idx[, 2]], sep = "__")
}

.regime_validate_collapsed_tip_labels <- function(tree, tips, state) {
  prospective_tip_labels <- tree$tip.label
  if (length(tips) > 1L) {
    prospective_tip_labels <- prospective_tip_labels[
      !prospective_tip_labels %in% tips[-1L]
    ]
  }
  prospective_tip_labels[prospective_tip_labels == tips[[1L]]] <- state
  .regime_validate_identifiers(
    prospective_tip_labels,
    "collapsed tree tip labels"
  )
}

.collapse_regime_phylogeny <- function(simmap_tree) {
  tree <- .regime_resolve_tree(tree = simmap_tree)
  tip_states <- .regime_tip_states(tree)
  out <- tree
  dropped_regimes <- character()

  for (state in unique(unname(tip_states))) {
    tips <- names(tip_states)[tip_states == state]
    tips <- intersect(tips, out$tip.label)
    if (length(tips) == 0L) { # nocov start
      next
    } # nocov end
    if (length(tips) == 1L) {
      .regime_validate_collapsed_tip_labels(out, tips, state)
      out$tip.label[out$tip.label == tips] <- state
      next
    }
    if (ape::is.monophyletic(ape::as.phylo(out), tips)) {
      .regime_validate_collapsed_tip_labels(out, tips, state)
      out <- ape::drop.tip(out, tips[-1L])
      out$tip.label[out$tip.label == tips[[1L]]] <- state
    } else {
      dropped_regimes <- c(dropped_regimes, state)
      out <- ape::drop.tip(out, tips)
    }
  }

  .regime_validate_identifiers(out$tip.label, "collapsed tree tip labels")

  if (length(dropped_regimes) > 0L) {
    warning(
      "The pGLS collapse dropped nonmonophyletic regime(s): ",
      paste(dropped_regimes, collapse = ", "),
      ". All tips assigned to these regimes were removed, matching the ",
      "manuscript workflow.",
      call. = FALSE
    )
  }

  ape::as.phylo(out)
}
