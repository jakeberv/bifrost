mean_cov_diagonals <- function(object) {
  # Extract covariance matrices from object$VCV
  vcv_matrices <- object$VCVs

  # Get matrix names (use indices if no names available)
  matrix_names <- names(vcv_matrices)
  if (is.null(matrix_names)) {
    matrix_names <- paste0("matrix_", seq_along(vcv_matrices))
  }

  # Extract diagonals and compute means vectorized
  diagonal_means <- vapply(vcv_matrices, function(mat) {
    mean(diag(mat))
  }, numeric(1))

  # Set names for the means vector
  names(diagonal_means) <- matrix_names

  # Create pairwise difference matrix (absolute differences)
  n_matrices <- length(diagonal_means)
  diff_matrix <- outer(diagonal_means, diagonal_means, function(x, y) abs(x - y))

  # Set row and column names for the difference matrix
  rownames(diff_matrix) <- matrix_names
  colnames(diff_matrix) <- matrix_names

  # Return results as a list
  return(list(
    diagonal_means = diagonal_means,
    pairwise_differences = diff_matrix
  ))
}



plot_heatmap_means <- function(analysis_result) {
  heatmap(analysis_result$pairwise_differences,
          symm = TRUE,
          main = "Pairwise Differences in Diagonal Means",
          xlab = "Matrix", ylab = "Matrix")
}


wasserstein_cov_diagonals <- function(object) {
  # Check if transport package is available
  if (!requireNamespace("transport", quietly = TRUE)) {
    stop("Package 'transport' is required for Wasserstein distance calculations")
  }

  # Extract covariance matrices from object$VCVs
  vcv_matrices <- object$VCVs

  # Get matrix names (use indices if no names available)
  matrix_names <- names(vcv_matrices)
  if (is.null(matrix_names)) {
    matrix_names <- paste0("matrix_", seq_along(vcv_matrices))
  }

  # Extract all diagonals first (vectorized)
  diagonals <- lapply(vcv_matrices, diag)
  names(diagonals) <- matrix_names

  # Create indices for vectorized computation
  n_matrices <- length(diagonals)
  indices <- seq_len(n_matrices)

  # Vectorized pairwise distance computation using outer
  wasserstein_matrix <- outer(indices, indices, function(i, j) {
    mapply(function(idx_i, idx_j) {
      if (idx_i == idx_j) {
        return(0)
      } else {
        return(transport::wasserstein1d(diagonals[[idx_i]], diagonals[[idx_j]]))
      }
    }, i, j)
  })

  # Set row and column names
  rownames(wasserstein_matrix) <- matrix_names
  colnames(wasserstein_matrix) <- matrix_names

  # Return results as a list
  return(list(
    diagonals = diagonals,
    wasserstein_distances = wasserstein_matrix
  ))
}

# Plot Wasserstein distances
plot_heatmap_wasserstein <- function(analysis_result) {
  heatmap(analysis_result$wasserstein_distances,
          symm = TRUE,
          main = "Pairwise Wasserstein Distances of Diagonal Values",
          xlab = "Matrix", ylab = "Matrix")
}


# Kolmorogov-Smirnov distances with p-values and q-values
ks_cov_diagonals <- function(object) {
  # Check if qvalue package is available
  if (!requireNamespace("qvalue", quietly = TRUE)) {
    stop("Package 'qvalue' is required for q-value calculations")
  }

  # Extract covariance matrices from object$VCVs
  vcv_matrices <- object$VCVs

  # Get matrix names (use indices if no names available)
  matrix_names <- names(vcv_matrices)
  if (is.null(matrix_names)) {
    matrix_names <- paste0("matrix_", seq_along(vcv_matrices))
  }

  # Extract all diagonals first (vectorized)
  diagonals <- lapply(vcv_matrices, diag)
  names(diagonals) <- matrix_names

  # Create indices for vectorized computation
  n_matrices <- length(diagonals)
  indices <- seq_len(n_matrices)

  # Matrices to store results
  ks_matrix <- matrix(0, nrow = n_matrices, ncol = n_matrices)
  pval_matrix <- matrix(1, nrow = n_matrices, ncol = n_matrices)  # Diagonal = 1 (no difference)

  # Vectorized pairwise KS test computation
  results <- outer(indices, indices, function(i, j) {
    mapply(function(idx_i, idx_j) {
      if (idx_i == idx_j) {
        return(list(statistic = 0, p.value = 1))
      } else {
        ks_result <- ks.test(diagonals[[idx_i]], diagonals[[idx_j]])
        return(list(statistic = as.numeric(ks_result$statistic),
                    p.value = ks_result$p.value))
      }
    }, i, j, SIMPLIFY = FALSE)
  })

  # Extract statistics and p-values
  for (i in 1:n_matrices) {
    for (j in 1:n_matrices) {
      ks_matrix[i, j] <- results[[i, j]]$statistic
      pval_matrix[i, j] <- results[[i, j]]$p.value
    }
  }

  # Calculate q-values from upper triangle p-values (excluding diagonal)
  upper_tri_indices <- which(upper.tri(pval_matrix), arr.ind = TRUE)
  upper_tri_pvals <- pval_matrix[upper_tri_indices]

  # Only calculate q-values if we have valid p-values
  if (length(upper_tri_pvals) > 0 && any(upper_tri_pvals < 1)) {
    tryCatch({
      qval_result <- qvalue::qvalue(upper_tri_pvals)
      upper_tri_qvals <- qval_result$qvalues

      # Create q-value matrix (symmetric)
      qval_matrix <- matrix(1, nrow = n_matrices, ncol = n_matrices)  # Diagonal = 1
      qval_matrix[upper_tri_indices] <- upper_tri_qvals
      qval_matrix[lower.tri(qval_matrix)] <- t(qval_matrix)[lower.tri(qval_matrix)]

    }, error = function(e) {
      warning("Q-value calculation failed: ", e$message, ". Returning NULL.")
      qval_matrix <- NULL
    })
  } else {
    qval_matrix <- pval_matrix
  }

  # Set row and column names
  rownames(ks_matrix) <- matrix_names
  colnames(ks_matrix) <- matrix_names
  rownames(pval_matrix) <- matrix_names
  colnames(pval_matrix) <- matrix_names
  rownames(qval_matrix) <- matrix_names
  colnames(qval_matrix) <- matrix_names

  # Return results as a list
  return(list(
    diagonals = diagonals,
    ks_distances = ks_matrix,
    pvalues = pval_matrix,
    qvalues = qval_matrix
  ))
}

# Energy distance version with p-values and q-values
energy_cov_diagonals <- function(object, R = 99) {
  # Check if energy package is available
  if (!requireNamespace("energy", quietly = TRUE)) {
    stop("Package 'energy' is required for energy distance calculations")
  }

  # Check if qvalue package is available
  if (!requireNamespace("qvalue", quietly = TRUE)) {
    stop("Package 'qvalue' is required for q-value calculations")
  }

  # Extract covariance matrices from object$VCVs
  vcv_matrices <- object$VCVs

  # Get matrix names (use indices if no names available)
  matrix_names <- names(vcv_matrices)
  if (is.null(matrix_names)) {
    matrix_names <- paste0("matrix_", seq_along(vcv_matrices))
  }

  # Extract all diagonals first (vectorized)
  diagonals <- lapply(vcv_matrices, diag)
  names(diagonals) <- matrix_names

  # Create indices for vectorized computation
  n_matrices <- length(diagonals)

  # Matrices to store results
  energy_matrix <- matrix(0, nrow = n_matrices, ncol = n_matrices)
  pval_matrix <- matrix(1, nrow = n_matrices, ncol = n_matrices)

  # Note: This version will be slower due to bootstrap resampling
  cat("Computing energy distances with", R, "bootstrap replicates...\n")

  for (i in 1:n_matrices) {
    for (j in i:n_matrices) {
      if (i == j) {
        energy_matrix[i, j] <- 0
        pval_matrix[i, j] <- 1
      } else {
        energy_result <- energy::eqdist.etest(rbind(matrix(diagonals[[i]]), matrix(diagonals[[j]])), sizes = c(c(length(diagonals[[i]]), length(diagonals[[j]]))), distance = FALSE, R = R)
        # When R > 0, returns a list with $statistic and $p.value
        energy_matrix[i, j] <- energy_result$statistic
        energy_matrix[j, i] <- energy_result$statistic  # Symmetric
        pval_matrix[i, j] <- energy_result$p.value
        pval_matrix[j, i] <- energy_result$p.value
      }
    }
  }

  # Calculate q-values from upper triangle p-values (excluding diagonal)
  upper_tri_indices <- which(upper.tri(pval_matrix), arr.ind = TRUE)
  upper_tri_pvals <- pval_matrix[upper_tri_indices]

  # Only calculate q-values if we have valid p-values
  if (length(upper_tri_pvals) > 0 && any(upper_tri_pvals < 1)) {
    tryCatch({
      qval_result <- qvalue::qvalue(upper_tri_pvals)
      upper_tri_qvals <- qval_result$qvalues

      # Create q-value matrix (symmetric)
      qval_matrix <- matrix(1, nrow = n_matrices, ncol = n_matrices)  # Diagonal = 1
      qval_matrix[upper_tri_indices] <- upper_tri_qvals
      qval_matrix[lower.tri(qval_matrix)] <- t(qval_matrix)[lower.tri(qval_matrix)]

    }, error = function(e) {
      warning("Q-value calculation failed: ", e$message, ". Returning NULL.")
      qval_matrix <- NULL
    })
  } else {
    qval_matrix <- pval_matrix
  }

  # Set row and column names
  rownames(energy_matrix) <- matrix_names
  colnames(energy_matrix) <- matrix_names
  rownames(pval_matrix) <- matrix_names
  colnames(pval_matrix) <- matrix_names
  rownames(qval_matrix) <- matrix_names
  colnames(qval_matrix) <- matrix_names

  # Return results as a list
  return(list(
    diagonals = diagonals,
    energy_distances = energy_matrix,
    pvalues = pval_matrix,
    qvalues = qval_matrix
  ))
}
