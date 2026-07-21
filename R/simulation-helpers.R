.simulation_check_integer_scalar <- function(x,
                                             name,
                                             minimum = NULL,
                                             maximum = NULL,
                                             allow_null = FALSE,
                                             message = NULL) {
  if (is.null(x) && isTRUE(allow_null)) {
    return(NULL)
  }

  valid <- is.numeric(x) &&
    length(x) == 1L &&
    !is.na(x) &&
    is.finite(x) &&
    x == floor(x) &&
    x >= -.Machine$integer.max &&
    x <= .Machine$integer.max

  if (valid && !is.null(minimum)) {
    valid <- x >= minimum
  }
  if (valid && !is.null(maximum)) {
    valid <- x <= maximum
  }

  if (!valid) {
    if (is.null(message)) {
      message <- paste0("`", name, "` must be a single finite integer.")
    }
    stop(message, call. = FALSE)
  }

  as.integer(x)
}

.simulation_check_integer_vector <- function(x,
                                             name,
                                             minimum = NULL,
                                             message = NULL) {
  valid <- is.numeric(x) &&
    length(x) >= 1L &&
    !anyNA(x) &&
    all(is.finite(x)) &&
    all(x == floor(x)) &&
    all(x >= -.Machine$integer.max) &&
    all(x <= .Machine$integer.max)

  if (valid && !is.null(minimum)) {
    valid <- all(x >= minimum)
  }

  if (!valid) {
    if (is.null(message)) {
      message <- paste0("`", name, "` must be a finite integer vector.")
    }
    stop(message, call. = FALSE)
  }

  as.integer(x)
}

.simulation_set_seed <- function(seed, kind = NULL) {
  seed <- .simulation_check_integer_scalar(
    seed,
    name = "seed",
    allow_null = TRUE,
    message = "`seed` must be NULL or a single finite integer."
  )
  if (is.null(seed)) {
    return(NULL)
  }

  old_kind <- RNGkind()
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  old_seed <- if (had_seed) {
    get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  } else {
    NULL
  }

  if (is.null(kind)) {
    set.seed(seed)
  } else {
    set.seed(seed, kind = kind)
  }
  list(had_seed = had_seed, old_seed = old_seed, old_kind = old_kind)
}

.simulation_restore_seed <- function(seed_state) {
  if (is.null(seed_state)) {
    return(invisible(NULL))
  }

  if (!is.null(seed_state$old_kind)) {
    do.call(RNGkind, as.list(seed_state$old_kind))
  }
  if (isTRUE(seed_state$had_seed)) {
    assign(".Random.seed", seed_state$old_seed, envir = .GlobalEnv)
  } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
    rm(".Random.seed", envir = .GlobalEnv)
  }

  invisible(NULL)
}

.simulation_replicate_seeds <- function(n_replicates, seeded) {
  if (!isTRUE(seeded)) {
    return(rep(list(NULL), n_replicates))
  }

  as.list(sample.int(.Machine$integer.max, n_replicates, replace = FALSE))
}

.simulation_seed_runner <- function() {
  set_seed <- .simulation_set_seed
  restore_seed <- .simulation_restore_seed
  function(seed, code) {
    seed_state <- set_seed(seed, kind = "L'Ecuyer-CMRG")
    if (is.null(seed_state)) {
      return(force(code))
    }

    on.exit(restore_seed(seed_state), add = TRUE)
    force(code)
  }
}

.simulation_check_unused_dots <- function(dots) {
  if (length(dots) == 0L) {
    return(invisible(TRUE))
  }

  dot_names <- names(dots)
  dot_names[is.na(dot_names) | dot_names == ""] <- "<unnamed>"
  stop("Unused argument(s): ", paste(dot_names, collapse = ", "), call. = FALSE)
}

.simulation_validate_covariance <- function(sigma, name = "sigma") {
  sigma <- as.matrix(sigma)
  valid_square <- length(dim(sigma)) == 2L &&
    nrow(sigma) >= 1L &&
    nrow(sigma) == ncol(sigma)
  if (!valid_square) {
    stop("`", name, "` must be a non-empty square matrix.", call. = FALSE)
  }
  if (!is.numeric(sigma) || anyNA(sigma) || any(!is.finite(sigma))) {
    stop("`", name, "` must contain only finite numeric values.", call. = FALSE)
  }
  if (!isTRUE(all.equal(sigma, t(sigma), tolerance = 1e-10))) {
    stop("`", name, "` must be symmetric.", call. = FALSE)
  }
  if (any(diag(sigma) <= 0)) {
    stop("`", name, "` must have a positive diagonal.", call. = FALSE)
  }
  if (any(eigen(sigma, symmetric = TRUE, only.values = TRUE)$values <= 0)) {
    stop("`", name, "` must be positive definite.", call. = FALSE)
  }

  sigma
}

.simulation_covariance_diagnostics <- function(sigma) {
  correlation <- stats::cov2cor(sigma)
  eigenvalues <- eigen(
    correlation,
    symmetric = TRUE,
    only.values = TRUE
  )$values
  mean_absolute_correlation <- if (nrow(correlation) == 1L) {
    0
  } else {
    mean(abs(correlation[lower.tri(correlation)]))
  }

  list(
    mean_absolute_correlation = mean_absolute_correlation,
    effective_dimensionality = sum(eigenvalues)^2 / sum(eigenvalues^2),
    condition_number = max(eigenvalues) / min(eigenvalues)
  )
}

.simulation_check_generator <- function(simulation_generator) {
  valid_generators <- c("original", "empirical")
  if (!is.character(simulation_generator) ||
      length(simulation_generator) != 1L ||
      is.na(simulation_generator) ||
      !simulation_generator %in% valid_generators) {
    stop(
      "`simulation_generator` must be one of: ",
      paste(valid_generators, collapse = ", "),
      ".",
      call. = FALSE
    )
  }

  simulation_generator
}

.simulation_generator_from_options <- function(options) {
  if (is.null(options$simulation_generator)) return("original")

  .simulation_check_generator(options$simulation_generator)
}

.simulation_draw_covariance <- function(template,
                                        simulation_generator = "original",
                                        covariance_df = NULL) {
  simulation_generator <- .simulation_check_generator(simulation_generator)
  n_traits <- .simulation_check_integer_scalar(
    template$n_response_traits,
    name = "template$n_response_traits",
    minimum = 1L
  )

  if (simulation_generator == "empirical") {
    sigma_source <- template$residual_covariance
    if (is.null(sigma_source)) {
      sigma_source <- matrix(
        template$covariance_mean,
        nrow = n_traits,
        ncol = n_traits
      )
      diag(sigma_source) <- template$variance_mean
    }
    sigma_source <- .simulation_validate_covariance(
      sigma_source,
      name = "template$residual_covariance"
    )
    if (nrow(sigma_source) != n_traits) {
      stop(
        "`template$residual_covariance` dimensions must match the number of response traits.",
        call. = FALSE
      )
    }

    if (is.null(covariance_df)) {
      covariance_df <- template$residual_df
      if (is.null(covariance_df)) {
        covariance_df <- max(
          n_traits + 1L,
          as.integer(template$n_tips) - 1L
        )
      }
    }
    covariance_df <- .simulation_check_integer_scalar(
      covariance_df,
      name = "covariance_df",
      minimum = n_traits,
      message = paste0(
        "`covariance_df` must be a single finite integer at least the number ",
        "of response traits (", n_traits, ")."
      )
    )

    sigma <- stats::rWishart(
      1L,
      df = covariance_df,
      Sigma = sigma_source / covariance_df
    )[, , 1L]
    sigma <- (sigma + t(sigma)) / 2
    dimnames(sigma) <- dimnames(sigma_source)

    return(list(
      sigma = sigma,
      covariance_df = covariance_df,
      diagnostics = .simulation_covariance_diagnostics(sigma)
    ))
  }

  sigma <- NULL
  attempt <- 1L
  while (attempt <= 100L && is.null(sigma)) {
    variances <- abs(stats::rnorm(
      n_traits,
      mean = template$variance_mean,
      sd = template$variance_sd
    ))

    if (n_traits == 1L) {
      candidate_sigma <- matrix(variances[1L], nrow = 1L, ncol = 1L)
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

    eigenvalues <- eigen(
      candidate_sigma,
      symmetric = TRUE,
      only.values = TRUE
    )$values
    if (all(eigenvalues > 0)) {
      sigma <- candidate_sigma
    } else {
      attempt <- attempt + 1L
    }
  }
  if (is.null(sigma)) {
    stop("Failed to generate a positive-definite covariance matrix.", call. = FALSE)
  }

  list(
    sigma = sigma,
    covariance_df = NA_integer_,
    diagnostics = .simulation_covariance_diagnostics(sigma)
  )
}

.simulation_transform_integration <- function(sigma,
                                              power,
                                              eigen_floor = 1e-8) {
  sigma <- .simulation_validate_covariance(sigma)
  valid_power <- is.numeric(power) &&
    length(power) == 1L &&
    !is.na(power) &&
    is.finite(power) &&
    power > 0
  if (!valid_power) {
    stop("`power` must be a single finite positive number.", call. = FALSE)
  }
  valid_floor <- is.numeric(eigen_floor) &&
    length(eigen_floor) == 1L &&
    !is.na(eigen_floor) &&
    is.finite(eigen_floor) &&
    eigen_floor > 0 &&
    eigen_floor < 1
  if (!valid_floor) {
    stop("`eigen_floor` must be a single number strictly between zero and one.", call. = FALSE)
  }

  diagnostics_before <- .simulation_covariance_diagnostics(sigma)
  if (power == 1) {
    return(list(
      sigma = sigma,
      diagnostics = list(
        power = power,
        mean_absolute_correlation_before = diagnostics_before$mean_absolute_correlation,
        mean_absolute_correlation_after = diagnostics_before$mean_absolute_correlation,
        effective_dimensionality_before = diagnostics_before$effective_dimensionality,
        effective_dimensionality_after = diagnostics_before$effective_dimensionality,
        condition_number_before = diagnostics_before$condition_number,
        condition_number_after = diagnostics_before$condition_number
      )
    ))
  }

  marginal_sd <- sqrt(diag(sigma))
  correlation <- stats::cov2cor(sigma)
  decomposition <- eigen(correlation, symmetric = TRUE)
  minimum_eigenvalue <- max(decomposition$values) * eigen_floor
  stable_eigenvalues <- pmax(decomposition$values, minimum_eigenvalue)
  powered_eigenvalues <- stable_eigenvalues^power
  powered_eigenvalues <- pmax(
    powered_eigenvalues,
    max(powered_eigenvalues) * eigen_floor
  )
  powered_correlation <-
    sweep(decomposition$vectors, 2L, powered_eigenvalues, `*`) %*%
    t(decomposition$vectors)
  powered_correlation <- (powered_correlation + t(powered_correlation)) / 2
  powered_correlation <- stats::cov2cor(powered_correlation)
  powered_correlation <- (powered_correlation + t(powered_correlation)) / 2
  diag(powered_correlation) <- 1

  transformed <- (marginal_sd * powered_correlation) *
    rep(marginal_sd, each = nrow(sigma))
  transformed <- (transformed + t(transformed)) / 2
  diag(transformed) <- diag(sigma)
  dimnames(transformed) <- dimnames(sigma)
  diagnostics_after <- .simulation_covariance_diagnostics(transformed)

  list(
    sigma = transformed,
    diagnostics = list(
      power = power,
      mean_absolute_correlation_before = diagnostics_before$mean_absolute_correlation,
      mean_absolute_correlation_after = diagnostics_after$mean_absolute_correlation,
      effective_dimensionality_before = diagnostics_before$effective_dimensionality,
      effective_dimensionality_after = diagnostics_after$effective_dimensionality,
      condition_number_before = diagnostics_before$condition_number,
      condition_number_after = diagnostics_after$condition_number
    )
  )
}
