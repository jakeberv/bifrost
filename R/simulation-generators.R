#' Simulate a Dataset-Matched Null Replicate
#'
#' @description
#' Generate a single no-shift simulation replicate under a uniform Brownian
#' motion process using the dataset-matched covariance-generation framework.
#' By default, each replicate uses the element-sampling operations from the
#' published simulation code. The empirical generator is available explicitly
#' for full-covariance draws centered on a [`createSimulationTemplate()`] object.
#'
#' @param template A `bifrost_simulation_template` returned by
#'   [createSimulationTemplate()].
#' @param tree_tip_count Optional integer tip count for random subtree sampling.
#'   If `NULL`, the full empirical tree from `template` is used.
#' @param seed Optional integer random seed. When supplied, the function
#'   restores the caller's previous RNG state before returning.
#' @param simulation_generator Simulation generator to use. `"original"`, the
#'   default, reproduces the element-sampling operations used in the published
#'   simulation code. `"empirical"` draws around the full fitted residual
#'   covariance. Supply exactly one supported string; explicitly supplied
#'   vectors are rejected. When omitted, `"original"` is used.
#' @param covariance_df Optional integer degrees of freedom for empirical
#'   Wishart draws. It must be at least the number of response traits. The
#'   default uses the residual degrees of freedom stored in `template`.
#' @details
#' The default `"original"` generator samples marginal variance and covariance
#' summaries and retains the published mathematical operations exactly. Under
#' the explicit `"empirical"` generator, a new covariance matrix \eqn{W} is drawn as
#' \eqn{W \sim Wishart(\nu, S / \nu)}, where \eqn{S} is the full empirical
#' residual covariance and \eqn{\nu} is `covariance_df`. Thus
#' \eqn{E[W] = S}: empirical integration is retained on average while individual
#' replicates vary.
#'
#' For templates calibrated from richer global
#' models, the empirical fitted mean structure is added back to the simulated
#' residual process before the downstream intercept-only search is run on the
#' regenerated response block. The returned `trait_data` therefore contains only
#' the simulated response variables, even when the global calibration model
#' included predictors. This is the lower-level replicate generator used by
#' [runFalsePositiveSimulationStudy()]; use that wrapper for replicated
#' false-positive simulation workflows.
#'
#' @return A list of class `bifrost_simulation_replicate_null` containing the
#'   sampled tree, a single-regime baseline tree, the simulated response matrix,
#'   a response-only `trait_data` object used for downstream model fitting, and
#'   the generated covariance matrix, generator metadata, resolved
#'   covariance degrees of freedom, and covariance diagnostics.
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
                                simulation_generator = c(
                                  "original",
                                  "empirical"
                                ),
                                covariance_df = NULL) {
  simulation_generator_missing <- missing(simulation_generator)

  if (!inherits(template, "bifrost_simulation_template")) {
    stop("template must be a 'bifrost_simulation_template' object.")
  }
  if (simulation_generator_missing) {
    simulation_generator <- "original"
  }
  simulation_generator <- .simulation_check_generator(simulation_generator)
  if (!is.null(covariance_df)) {
    covariance_df <- .simulation_check_integer_scalar(
      covariance_df,
      name = "covariance_df",
      minimum = template$n_response_traits,
      message = paste0(
        "`covariance_df` must be a single finite integer at least the number ",
        "of response traits (", template$n_response_traits, ")."
      )
    )
  }
  seed_state <- .simulation_set_seed(seed)
  on.exit(.simulation_restore_seed(seed_state), add = TRUE)

  tree_tip_count <- .simulation_check_integer_scalar(
    tree_tip_count,
    name = "tree_tip_count",
    minimum = 2L,
    allow_null = TRUE,
    message = "tree_tip_count must be NULL or a single integer >= 2."
  )
  if (!is.null(tree_tip_count)) {
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

  covariance_draw <- tryCatch(
    .simulation_draw_covariance(
      template = template,
      simulation_generator = simulation_generator,
      covariance_df = covariance_df
    ),
    error = function(e) NULL
  )
  if (is.null(covariance_draw)) {
    stop("Failed to generate a positive-definite covariance matrix for the null simulation.")
  }
  sigma <- covariance_draw$sigma

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
    simulation_generator = simulation_generator,
    covariance_df = covariance_draw$covariance_df,
    covariance_diagnostics = covariance_draw$diagnostics,
    generating_scenario = "null"
  )

  class(result) <- c("bifrost_simulation_replicate_null", class(result))
  result
}

#' Simulate a Dataset-Matched Shifted Replicate
#'
#' @description
#' Generate a single known-shift simulation replicate using the dataset-matched
#' multi-shift framework. The function supports both the proportional
#' generating-model scenario and a generator-specific non-proportional robustness
#' scenario.
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
#'   the generator-specific non-proportional robustness scenario. With the
#'   default `"original"` generator, the latter reproduces the published
#'   correlation transform; with `simulation_generator = "empirical"`, it is an
#'   integration-rate trade-off.
#' @param scale_factor_range Strictly positive numeric length-2 vector giving
#'   the full range of possible shift scalars.
#' @param exclude_range Numeric length-2 vector specifying the central interval
#'   excluded from sampled shift scalars.
#' @param buffer Integer minimum node distance between simulated shifts.
#' @param seed Optional integer random seed. When supplied, the function
#'   restores the caller's previous RNG state before returning.
#' @param simulation_generator Covariance and correlation-shift generator to
#'   use. `"original"`, the default, reproduces the published simulation
#'   operations exactly. `"empirical"` draws around the full fitted covariance
#'   and changes integration and marginal evolutionary variance in opposing
#'   directions. Supply exactly one supported string; explicitly supplied
#'   vectors are rejected. When omitted, `"original"` is used.
#' @param covariance_df Optional integer degrees of freedom for empirical
#'   Wishart draws. It must be at least the number of response traits. The
#'   default uses the residual degrees of freedom stored in `template`.
#' @param integration_power_range Strictly positive increasing numeric length-2
#'   vector spanning one. Corrected correlation shifts draw spectral powers from
#'   the two tails of this range.
#' @param integration_exclude_range Increasing numeric length-2 vector that
#'   contains one and lies strictly inside `integration_power_range`. Powers in
#'   this central interval are excluded so shifts have a material effect.
#' @param eigen_floor Numeric stability floor, as a fraction of the largest
#'   eigenvalue, used by empirical correlation transforms.
#' @details
#' This function generates shifted residual processes around the fitted mean
#' structure stored in the template:
#' \itemize{
#'   \item random subtree sampling from the empirical tree,
#'   \item shift placement subject to clade-size and non-overlap constraints,
#'   \item buffer enforcement using node distances,
#'   \item generation of a fresh ancestral covariance matrix centered on the
#'         full empirical residual covariance for each empirical replicate,
#'   \item construction of derived regimes by either proportional scaling or the
#'         generator-specific non-proportional robustness scenario, and
#'   \item simulation of multivariate BMM residuals that are added to the
#'         empirical fitted mean structure from the global calibration model.
#' }
#'
#' Under `scale_mode = "proportional"`, derived regimes are scalar multiples of
#' the ancestral covariance matrix and therefore match the generating assumptions
#' of the current `bifrost` BMM search. Under `scale_mode = "correlation"`, the
#' empirical generator raises the eigenvalues of the ancestral correlation matrix
#' to a sampled positive power and renormalizes the result to a correlation
#' matrix. Powers below one reduce integration; powers above one increase it.
#' The derived covariance is then divided by the same power. Consequently,
#' powers below one pair reduced integration with increased marginal
#' evolutionary variance, whereas powers above one pair increased integration
#' with decreased marginal variance. More precisely,
#' `diag(Sigma_derived) = diag(Sigma_ancestral) / power`. This deliberate
#' integration-rate trade-off remains a model-misspecification robustness test
#' because the search fits proportional covariance shifts.
#'
#' The `"original"` generator is supplied for exact reproduction. Its
#' correlation scenario scales correlations and inversely scales marginal
#' variances according to the original manuscript operations; it should not be
#' interpreted as holding marginal variances fixed.
#' Downstream simulation studies still evaluate intercept-only shift searches on
#' the regenerated response block, so the returned `trait_data` contains only
#' the simulated response variables. This is the lower-level shifted-replicate
#' generator used by [runShiftRecoverySimulationStudy()] and
#' [runSearchTuningGrid()]; use those wrappers for replicated shift-recovery or
#' tuning-grid workflows.
#'
#' @return A list of class `bifrost_simulation_replicate_shifted` containing the
#'   painted generating tree, the true shift nodes, the simulated response
#'   matrix, a response-only `trait_data` object used for downstream model
#'   fitting, regime covariance matrices, sampled scale factors or integration
#'   powers, reciprocal marginal-variance scale factors, generator
#'   metadata, covariance diagnostics, and a scenario label. For schema
#'   compatibility, `sampledScaleFactors` contains the sampled integration
#'   powers under the empirical correlation generator; the explicit
#'   `sampledIntegrationPowers` and `sampledVarianceScaleFactors` fields should
#'   be preferred for interpretation.
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
                                   simulation_generator = c(
                                     "original",
                                     "empirical"
                                   ),
                                   covariance_df = NULL,
                                   integration_power_range = c(0.5, 1.25),
                                   integration_exclude_range = c(0.8, 1.1),
                                   eigen_floor = 1e-8) {
  simulation_generator_missing <- missing(simulation_generator)

  scale_mode <- match.arg(scale_mode)

  if (!inherits(template, "bifrost_simulation_template")) {
    stop("template must be a 'bifrost_simulation_template' object.")
  }
  if (simulation_generator_missing) {
    simulation_generator <- "original"
  }
  simulation_generator <- .simulation_check_generator(simulation_generator)
  if (!is.null(covariance_df)) {
    covariance_df <- .simulation_check_integer_scalar(
      covariance_df,
      name = "covariance_df",
      minimum = template$n_response_traits,
      message = paste0(
        "`covariance_df` must be a single finite integer at least the number ",
        "of response traits (", template$n_response_traits, ")."
      )
    )
  }
  seed_state <- .simulation_set_seed(seed)
  on.exit(.simulation_restore_seed(seed_state), add = TRUE)

  num_shifts <- .simulation_check_integer_scalar(
    num_shifts,
    name = "num_shifts",
    minimum = 1L,
    message = "num_shifts must be a single integer >= 1."
  )
  min_shift_tips <- .simulation_check_integer_scalar(
    min_shift_tips,
    name = "min_shift_tips",
    minimum = 1L,
    message = "min_shift_tips must be a single integer >= 1."
  )
  max_shift_tips <- .simulation_check_integer_scalar(
    max_shift_tips,
    name = "max_shift_tips",
    minimum = 1L,
    message = "max_shift_tips must be a single integer >= min_shift_tips."
  )
  if (max_shift_tips < min_shift_tips) {
    stop("max_shift_tips must be >= min_shift_tips.")
  }
  buffer <- .simulation_check_integer_scalar(
    buffer,
    name = "buffer",
    minimum = 0L,
    message = "buffer must be a single non-negative integer."
  )
  if (!is.numeric(scale_factor_range) || length(scale_factor_range) != 2L ||
      anyNA(scale_factor_range) || any(!is.finite(scale_factor_range)) ||
      any(scale_factor_range <= 0) || diff(scale_factor_range) <= 0) {
    stop("scale_factor_range must be a finite positive numeric length-2 vector in increasing order.")
  }
  if (!is.numeric(exclude_range) || length(exclude_range) != 2L ||
      anyNA(exclude_range) || diff(exclude_range) <= 0) {
    stop("exclude_range must be a numeric length-2 vector in increasing order.")
  }
  if (exclude_range[1] <= scale_factor_range[1] || exclude_range[2] >= scale_factor_range[2]) {
    stop("exclude_range must lie strictly inside scale_factor_range.")
  }
  if (!is.numeric(integration_power_range) ||
      length(integration_power_range) != 2L ||
      anyNA(integration_power_range) ||
      any(!is.finite(integration_power_range)) ||
      any(integration_power_range <= 0) ||
      diff(integration_power_range) <= 0 ||
      integration_power_range[1L] >= 1 ||
      integration_power_range[2L] <= 1) {
    stop(
      paste0(
        "integration_power_range must be a finite positive numeric length-2 ",
        "vector in increasing order that spans 1."
      )
    )
  }
  if (!is.numeric(integration_exclude_range) ||
      length(integration_exclude_range) != 2L ||
      anyNA(integration_exclude_range) ||
      any(!is.finite(integration_exclude_range)) ||
      diff(integration_exclude_range) <= 0 ||
      integration_exclude_range[1L] >= 1 ||
      integration_exclude_range[2L] <= 1) {
    stop(
      paste0(
        "integration_exclude_range must be a finite numeric length-2 vector ",
        "in increasing order that contains 1."
      )
    )
  }
  if (integration_exclude_range[1L] <= integration_power_range[1L] ||
      integration_exclude_range[2L] >= integration_power_range[2L]) {
    stop(
      "integration_exclude_range must lie strictly inside integration_power_range."
    )
  }
  valid_eigen_floor <- is.numeric(eigen_floor) &&
    length(eigen_floor) == 1L &&
    !is.na(eigen_floor) &&
    is.finite(eigen_floor) &&
    eigen_floor > 0 &&
    eigen_floor < 1
  if (!valid_eigen_floor) {
    stop("eigen_floor must be a single number strictly between zero and one.")
  }
  if (scale_mode == "correlation" && template$n_response_traits < 2L) {
    stop("scale_mode = 'correlation' requires at least two response traits.")
  }
  tree_tip_count <- .simulation_check_integer_scalar(
    tree_tip_count,
    name = "tree_tip_count",
    minimum = 2L,
    allow_null = TRUE,
    message = "tree_tip_count must be NULL or a single integer >= 2."
  )
  if (!is.null(tree_tip_count)) {
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

    covariance_draw <- tryCatch(
      .simulation_draw_covariance(
        template = template,
        simulation_generator = simulation_generator,
        covariance_df = covariance_df
      ),
      error = function(e) NULL
    )
    if (is.null(covariance_draw)) {
      next
    }
    ancestral_sigma <- covariance_draw$sigma

    sigma_list <- list(ancestral = ancestral_sigma)
    vcvs <- list(ancestral = ancestral_sigma)
    covariance_diagnostics <- list(
      ancestral = covariance_draw$diagnostics
    )
    sampled_scale_factors <- numeric(num_shifts)
    sampled_integration_powers <- rep(NA_real_, num_shifts)
    sampled_variance_scale_factors <- numeric(num_shifts)
    tree_with_shifts <- generating_tree
    derived_failed <- FALSE

    for (shift_count in seq_len(num_shifts)) {
      derived_state <- paste("derived", shift_count, sep = "_")
      use_integration_power <-
        scale_mode == "correlation" && simulation_generator == "empirical"
      active_range <- if (use_integration_power) {
        integration_power_range
      } else {
        scale_factor_range
      }
      active_exclude <- if (use_integration_power) {
        integration_exclude_range
      } else {
        exclude_range
      }
      part1 <- stats::runif(1, min = active_range[1], max = active_exclude[1])
      part2 <- stats::runif(1, min = active_exclude[2], max = active_range[2])
      current_scale_factor <- sample(c(part1, part2), 1L)
      sampled_scale_factors[shift_count] <- current_scale_factor
      if (use_integration_power) {
        sampled_integration_powers[shift_count] <- current_scale_factor
      }
      sampled_variance_scale_factors[shift_count] <- if (
        scale_mode == "correlation"
      ) {
        1 / current_scale_factor
      } else {
        current_scale_factor
      }

      derived_sigma <- NULL
      if (use_integration_power) {
        transformed <- .simulation_transform_integration(
          ancestral_sigma,
          power = current_scale_factor,
          eigen_floor = eigen_floor
        )
        candidate_sigma <-
          sampled_variance_scale_factors[shift_count] * transformed$sigma
        transformed$diagnostics$marginal_variance_scale <-
          sampled_variance_scale_factors[shift_count]
        transformed$diagnostics$mean_marginal_variance_before <-
          mean(diag(ancestral_sigma))
        transformed$diagnostics$mean_marginal_variance_after <-
          mean(diag(candidate_sigma))
      } else if (scale_mode == "correlation") {
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
      }

      if (is.null(derived_sigma)) {
        derived_failed <- TRUE
        break
      }

      sigma_list[[derived_state]] <- derived_sigma
      vcvs[[derived_state]] <- derived_sigma
      covariance_diagnostics[[derived_state]] <- if (use_integration_power) {
        transformed$diagnostics
      } else {
        .simulation_covariance_diagnostics(derived_sigma)
      }
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
      sampledIntegrationPowers = sampled_integration_powers,
      sampledVarianceScaleFactors = sampled_variance_scale_factors,
      simulation_generator = simulation_generator,
      covariance_df = covariance_draw$covariance_df,
      covariance_diagnostics = covariance_diagnostics,
      generating_scenario = scale_mode
    )
  }

  if (is.null(result)) {
    stop("Failed to place all requested shifts after 100 attempts.")
  }

  class(result) <- c("bifrost_simulation_replicate_shifted", class(result))
  result
}
