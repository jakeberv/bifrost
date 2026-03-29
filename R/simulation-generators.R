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
                                ...) {
  extra_args <- list(...)
  if ("preserve_predictors" %in% names(extra_args)) {
    stop(
      "preserve_predictors has been removed. ",
      "Simulation replicates now always contain the response block only."
    )
  }
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
                                   ...) {
  extra_args <- list(...)
  if ("preserve_predictors" %in% names(extra_args)) {
    stop(
      "preserve_predictors has been removed. ",
      "Simulation replicates now always contain the response block only."
    )
  }
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
