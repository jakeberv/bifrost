bifrost_skip_search_characterization_deps <- function() {
  testthat::skip_if_not_installed("ape")
  testthat::skip_if_not_installed("phytools")
  testthat::skip_if_not_installed("mvMORPH")
  testthat::skip_if_not_installed("future")
}

bifrost_characterization_formula_label <- function(formula) {
  paste(deparse(formula, width.cutoff = 500), collapse = " ")
}

bifrost_characterization_run_search <- function(...) {
  suppressWarnings(suppressMessages(searchOptimalConfiguration(...)))
}

bifrost_characterization_tree_summary <- function(tree) {
  if (is.null(tree)) {
    return(NULL)
  }

  state_counts <- function(type) {
    states <- tryCatch(
      phytools::getStates(tree, type = type),
      error = function(e) NULL
    )
    if (is.null(states)) {
      return(NULL)
    }
    as.integer(sort(table(as.character(states))))
  }

  state_count_names <- function(type) {
    states <- tryCatch(
      phytools::getStates(tree, type = type),
      error = function(e) NULL
    )
    if (is.null(states)) {
      return(NULL)
    }
    names(sort(table(as.character(states))))
  }

  list(
    class = class(tree),
    n_tips = ape::Ntip(tree),
    n_nodes = ape::Nnode(tree),
    edge = tree$edge,
    edge_length = tree$edge.length,
    mapped_edge = tree$mapped.edge,
    tip_state_names = state_count_names("tips"),
    tip_state_counts = state_counts("tips"),
    node_state_names = state_count_names("nodes"),
    node_state_counts = state_counts("nodes"),
    both_state_names = state_count_names("both"),
    both_state_counts = state_counts("both")
  )
}

bifrost_characterization_model_summary <- function(model) {
  if (is.null(model)) {
    return(NULL)
  }

  model_call <- model$call
  call_value <- function(name) {
    if (is.null(model_call[[name]])) {
      return(NULL)
    }
    as.character(model_call[[name]])
  }

  list(
    class = class(model),
    call_model = call_value("model"),
    call_method = call_value("method"),
    formula = if (!is.null(model$formula)) {
      bifrost_characterization_formula_label(model$formula)
    } else {
      NULL
    },
    y_dim = dim(model$Y),
    param = model$param,
    log_lik = as.numeric(stats::logLik(model))
  )
}

bifrost_characterization_vcv_summary <- function(vcv) {
  if (is.null(vcv)) {
    return(NULL)
  }
  lapply(vcv, function(mat) {
    storage.mode(mat) <- "double"
    mat
  })
}

bifrost_characterization_history_summary <- function(history, ic_used) {
  if (is.null(history)) {
    return(NULL)
  }

  fit_summary <- NULL
  if (is.list(history$fits)) {
    fit_summary <- lapply(history$fits, function(entry) {
      model_ic <- NA_real_
      if (!is.null(entry$model)) {
        model_ic <- if (identical(ic_used, "GIC")) {
          entry$model$GIC$GIC
        } else {
          entry$model$BIC$BIC
        }
      }

      list(
        has_model = !is.null(entry$model),
        accepted = isTRUE(entry$accepted),
        delta_ic = entry$delta_ic,
        error = if (!is.null(entry$error)) entry$error else NULL,
        model_ic = model_ic
      )
    })
  }

  list(
    names = names(history),
    fits_length = if (is.list(history$fits)) length(history$fits) else NULL,
    fits = fit_summary,
    ic_acceptance_matrix = history$ic_acceptance_matrix
  )
}

bifrost_characterization_weights_summary <- function(weights) {
  if (is.null(weights)) {
    return(NULL)
  }
  rownames(weights) <- NULL
  weights
}

bifrost_characterization_result_summary <- function(result) {
  list(
    class = class(result),
    names = names(result),
    user_input_names = names(result$user_input),
    IC_used = result$IC_used,
    baseline_ic = result$baseline_ic,
    optimal_ic = result$optimal_ic,
    global_delta_ic = result$baseline_ic - result$optimal_ic,
    num_candidates = result$num_candidates,
    shift_nodes_no_uncertainty = result$shift_nodes_no_uncertainty,
    transformed_tree = bifrost_characterization_tree_summary(
      result$tree_no_uncertainty_transformed
    ),
    untransformed_tree = bifrost_characterization_tree_summary(
      result$tree_no_uncertainty_untransformed
    ),
    model_no_uncertainty = bifrost_characterization_model_summary(
      result$model_no_uncertainty
    ),
    model_fit_history = bifrost_characterization_history_summary(
      result$model_fit_history,
      result$IC_used
    ),
    VCVs = bifrost_characterization_vcv_summary(result$VCVs),
    ic_weights = bifrost_characterization_weights_summary(result$ic_weights),
    warnings = if (!is.null(result$warnings)) unlist(result$warnings) else character()
  )
}

bifrost_characterization_matrix_data <- function(seed, n_tips, n_traits) {
  set.seed(seed)
  tree <- ape::rtree(n_tips)
  trait_data <- matrix(rnorm(n_tips * n_traits), ncol = n_traits)
  colnames(trait_data) <- paste0("y", seq_len(n_traits))
  rownames(trait_data) <- tree$tip.label
  list(tree = tree, trait_data = trait_data)
}

bifrost_characterization_named_data <- function(seed, n_tips) {
  set.seed(seed)
  tree <- ape::rtree(n_tips)
  data <- data.frame(
    y1 = rnorm(n_tips),
    y2 = rnorm(n_tips),
    size = exp(rnorm(n_tips)),
    grp = factor(rep(c("a", "b"), length.out = n_tips))
  )
  rownames(data) <- tree$tip.label
  list(tree = tree, trait_data = data)
}

bifrost_build_search_characterization <- function() {
  gic_no_shift <- bifrost_characterization_matrix_data(
    seed = 4101,
    n_tips = 24,
    n_traits = 2
  )
  bic_no_shift <- bifrost_characterization_matrix_data(
    seed = 4102,
    n_tips = 24,
    n_traits = 2
  )
  gic_forced_serial <- bifrost_characterization_matrix_data(
    seed = 4103,
    n_tips = 20,
    n_traits = 2
  )
  gic_forced_parallel <- bifrost_characterization_matrix_data(
    seed = 4104,
    n_tips = 20,
    n_traits = 2
  )
  named_formula <- bifrost_characterization_named_data(
    seed = 4105,
    n_tips = 20
  )
  legacy_formula <- bifrost_characterization_named_data(
    seed = 4106,
    n_tips = 20
  )

  list(
    gic_no_shift_matrix = bifrost_characterization_result_summary(
      bifrost_characterization_run_search(
        baseline_tree = gic_no_shift$tree,
        trait_data = gic_no_shift$trait_data,
        formula = "trait_data ~ 1",
        min_descendant_tips = 8,
        num_cores = 1,
        shift_acceptance_threshold = 1e9,
        plot = FALSE,
        IC = "GIC",
        store_model_fit_history = TRUE,
        method = "LL",
        uncertaintyweights_par = TRUE
      )
    ),
    bic_no_shift_matrix = bifrost_characterization_result_summary(
      bifrost_characterization_run_search(
        baseline_tree = bic_no_shift$tree,
        trait_data = bic_no_shift$trait_data,
        formula = "trait_data ~ 1",
        min_descendant_tips = 8,
        num_cores = 1,
        shift_acceptance_threshold = 1e9,
        plot = FALSE,
        IC = "BIC",
        store_model_fit_history = TRUE,
        method = "LL",
        uncertaintyweights_par = TRUE
      )
    ),
    gic_forced_shift_serial_weights = bifrost_characterization_result_summary(
      bifrost_characterization_run_search(
        baseline_tree = gic_forced_serial$tree,
        trait_data = gic_forced_serial$trait_data,
        formula = "trait_data ~ 1",
        min_descendant_tips = 7,
        num_cores = 1,
        shift_acceptance_threshold = -Inf,
        plot = FALSE,
        IC = "GIC",
        store_model_fit_history = TRUE,
        method = "LL",
        uncertaintyweights = TRUE
      )
    ),
    gic_forced_shift_parallel_weights = bifrost_characterization_result_summary(
      bifrost_characterization_run_search(
        baseline_tree = gic_forced_parallel$tree,
        trait_data = gic_forced_parallel$trait_data,
        formula = "trait_data ~ 1",
        min_descendant_tips = 7,
        num_cores = 1,
        shift_acceptance_threshold = -Inf,
        plot = FALSE,
        IC = "GIC",
        store_model_fit_history = TRUE,
        method = "LL",
        uncertaintyweights_par = TRUE
      )
    ),
    named_formula_gic = bifrost_characterization_result_summary(
      bifrost_characterization_run_search(
        baseline_tree = named_formula$tree,
        trait_data = named_formula$trait_data,
        formula = cbind(y1, y2) ~ log(size) * grp,
        min_descendant_tips = 8,
        num_cores = 1,
        shift_acceptance_threshold = 1e9,
        plot = FALSE,
        IC = "GIC",
        store_model_fit_history = FALSE,
        method = "LL"
      )
    ),
    legacy_indexed_formula_bic = bifrost_characterization_result_summary(
      bifrost_characterization_run_search(
        baseline_tree = legacy_formula$tree,
        trait_data = legacy_formula$trait_data,
        formula = "trait_data[, 1:2] ~ trait_data[, 3]",
        min_descendant_tips = 8,
        num_cores = 1,
        shift_acceptance_threshold = 1e9,
        plot = FALSE,
        IC = "BIC",
        store_model_fit_history = FALSE,
        method = "LL"
      )
    )
  )
}
