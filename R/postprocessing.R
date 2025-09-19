# =============================================================================
# Covariance Matrix Diagonal Analysis Functions
# =============================================================================
#
# This file provides functions to analyze and compare the diagonal elements
# of covariance matrices from searchOptimalConfiguration() outputs.
#
# OVERVIEW:
# These functions extract diagonal values from multiple covariance matrices and
# compare them using various distance/similarity metrics:
#
# DISTANCE METRICS:
# - Mean differences: Simple arithmetic differences between diagonal means
# - Variance differences: Differences in diagonal variance patterns
# - Wasserstein distance: Optimal transport distance between distributions
# - Kolmogorov-Smirnov: Non-parametric test for distribution differences
# - Energy distance: Another distribution comparison method with bootstrap tests
#
# TYPICAL WORKFLOW:
# 1. Run searchOptimalConfiguration() to get an object with $VCVs
# 2. Apply one of the distance functions (e.g., mean_cov_diagonals())
# 3. Visualize results with plot_cov_heatmap() or plot_pvalue_heatmap()
#
# DEPENDENCIES:
# - transport package (for Wasserstein distances)
# - energy package (for energy distances)
# - qvalue package (for multiple testing correction)
# - viridis package (optional, for better color schemes)
#
# =============================================================================

# Basic mean-based comparison of diagonal elements
mean_cov_diagonals <- function(object) {
  # Extract covariance matrices from object$VCVs
  vcv_matrices <- object$VCVs

  # Get matrix names (use indices if no names available)
  matrix_names <- names(vcv_matrices)
  if (is.null(matrix_names)) {
    matrix_names <- paste0("matrix_", seq_along(vcv_matrices))
  }

  # Extract diagonals and compute means vectorized
  diagonal_means <- vapply(vcv_matrices, function(mat) {
    mean(diag(mat))  # Mean of diagonal elements for each matrix
  }, numeric(1))

  # Set names for the means vector
  names(diagonal_means) <- matrix_names

  # Create pairwise difference matrix (absolute differences)
  # This shows how different the average variances are between matrices
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

# Variance-based comparison of diagonal variability
var_cov_diagonals <- function(object) {
  # Extract covariance matrices from object$VCVs
  vcv_matrices <- object$VCVs

  # Get matrix names (use indices if no names available)
  matrix_names <- names(vcv_matrices)
  if (is.null(matrix_names)) {
    matrix_names <- paste0("matrix_", seq_along(vcv_matrices))
  }

  # Extract diagonals and compute variance of diagonal elements
  # This measures how variable the diagonal elements are within each matrix
  diagonal_vars <- vapply(vcv_matrices, function(mat) {
    var(diag(mat))  # Variance of diagonal elements within each matrix
  }, numeric(1))

  # Set names for the variance vector
  names(diagonal_vars) <- matrix_names

  # Create pairwise difference matrix (absolute differences)
  # Shows differences in diagonal variability patterns between matrices
  n_matrices <- length(diagonal_vars)
  diff_matrix <- outer(diagonal_vars, diagonal_vars, function(x, y) abs(x - y))

  # Set row and column names for the difference matrix
  rownames(diff_matrix) <- matrix_names
  colnames(diff_matrix) <- matrix_names

  # Return results as a list
  return(list(
    diagonal_vars = diagonal_vars,
    variance_differences = diff_matrix
  ))
}


# Wasserstein (optimal transport) distance between diagonal distributions
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

  # Extract all diagonal vectors first (vectorized)
  # Wasserstein distance compares entire distributions, not just summary stats
  diagonals <- lapply(vcv_matrices, diag)
  names(diagonals) <- matrix_names

  # Create indices for vectorized computation
  n_matrices <- length(diagonals)
  indices <- seq_len(n_matrices)

  # Vectorized pairwise distance computation using outer
  # Wasserstein distance measures the "cost" of transforming one distribution to another
  wasserstein_matrix <- outer(indices, indices, function(i, j) {
    mapply(function(idx_i, idx_j) {
      if (idx_i == idx_j) {
        return(0)  # Distance from a distribution to itself is 0
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

# Kolmogorov-Smirnov test with statistical significance testing
ks_cov_diagonals <- function(object, qval = TRUE) {
  # Check if qvalue package is available (only if qval = TRUE)
  if (qval && !requireNamespace("qvalue", quietly = TRUE)) {
    stop("Package 'qvalue' is required for q-value calculations")
  }

  # Extract covariance matrices from object$VCVs
  vcv_matrices <- object$VCVs

  # Get matrix names (use indices if no names available)
  matrix_names <- names(vcv_matrices)
  if (is.null(matrix_names)) {
    matrix_names <- paste0("matrix_", seq_along(vcv_matrices))
  }

  # Extract all diagonal vectors first (vectorized)
  # KS test compares cumulative distribution functions
  diagonals <- lapply(vcv_matrices, diag)
  names(diagonals) <- matrix_names

  # Create indices for vectorized computation
  n_matrices <- length(diagonals)
  indices <- seq_len(n_matrices)

  # Matrices to store KS statistics and p-values
  ks_matrix <- matrix(0, nrow = n_matrices, ncol = n_matrices)
  pval_matrix <- matrix(1, nrow = n_matrices, ncol = n_matrices)  # Diagonal = 1 (no difference)

  # Vectorized pairwise KS test computation
  # KS test null hypothesis: both samples come from the same distribution
  results <- outer(indices, indices, function(i, j) {
    mapply(function(idx_i, idx_j) {
      if (idx_i == idx_j) {
        return(list(statistic = 0, p.value = 1))  # Same distribution
      } else {
        ks_result <- ks.test(diagonals[[idx_i]], diagonals[[idx_j]])
        return(list(statistic = as.numeric(ks_result$statistic),
                    p.value = ks_result$p.value))
      }
    }, i, j, SIMPLIFY = FALSE)
  })

  # Extract statistics and p-values from nested results
  for (i in 1:n_matrices) {
    for (j in 1:n_matrices) {
      ks_matrix[i, j] <- results[[i, j]]$statistic
      pval_matrix[i, j] <- results[[i, j]]$p.value
    }
  }

  # Calculate q-values for multiple testing correction (Benjamini-Hochberg FDR)
  if (qval) {
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
        warning("Q-value calculation failed: ", e$message, ". Returning placeholder.")
        qval_matrix <- "Q-value calculation failed - see warning message"
      })
    } else {
      qval_matrix <- "No valid p-values for q-value calculation"
    }
  } else {
    qval_matrix <- "Q-value calculation skipped (qval = FALSE)"
  }

  # Set row and column names
  rownames(ks_matrix) <- matrix_names
  colnames(ks_matrix) <- matrix_names
  rownames(pval_matrix) <- matrix_names
  colnames(pval_matrix) <- matrix_names

  # Only set names if qval_matrix is actually a matrix
  if (is.matrix(qval_matrix)) {
    rownames(qval_matrix) <- matrix_names
    colnames(qval_matrix) <- matrix_names
  }

  # Return results as a list
  return(list(
    diagonals = diagonals,
    ks_distances = ks_matrix,      # KS test statistics
    pvalues = pval_matrix,         # Raw p-values
    qvalues = qval_matrix          # Multiple testing corrected q-values
  ))
}

# Energy distance with bootstrap testing (slower but robust)
energy_cov_diagonals <- function(object, R = 99, qval = TRUE) {
  # Check if energy package is available
  if (!requireNamespace("energy", quietly = TRUE)) {
    stop("Package 'energy' is required for energy distance calculations")
  }

  # Check if qvalue package is available (only if qval = TRUE)
  if (qval && !requireNamespace("qvalue", quietly = TRUE)) {
    stop("Package 'qvalue' is required for q-value calculations")
  }

  # Extract covariance matrices from object$VCVs
  vcv_matrices <- object$VCVs

  # Get matrix names (use indices if no names available)
  matrix_names <- names(vcv_matrices)
  if (is.null(matrix_names)) {
    matrix_names <- paste0("matrix_", seq_along(vcv_matrices))
  }

  # Extract all diagonal vectors first (vectorized)
  # Energy distance is another way to compare distributions using E-statistics
  diagonals <- lapply(vcv_matrices, diag)
  names(diagonals) <- matrix_names

  # Create indices for computation
  n_matrices <- length(diagonals)

  # Matrices to store results
  energy_matrix <- matrix(0, nrow = n_matrices, ncol = n_matrices)
  pval_matrix <- matrix(1, nrow = n_matrices, ncol = n_matrices)

  # Note: This version will be slower due to bootstrap resampling
  cat("Computing energy distances with", R, "bootstrap replicates...\n")

  # Compute pairwise energy distances with bootstrap p-values
  for (i in 1:n_matrices) {
    for (j in i:n_matrices) {
      if (i == j) {
        energy_matrix[i, j] <- 0  # Distance from distribution to itself
        pval_matrix[i, j] <- 1
      } else {
        # Energy test with bootstrap resampling for p-value estimation
        energy_result <- energy::eqdist.etest(rbind(matrix(diagonals[[i]]), matrix(diagonals[[j]])), sizes = c(c(length(diagonals[[i]]), length(diagonals[[j]]))), distance = FALSE, R = R)
        # When R > 0, returns a list with $statistic and $p.value
        energy_matrix[i, j] <- energy_result$statistic
        energy_matrix[j, i] <- energy_result$statistic  # Symmetric
        pval_matrix[i, j] <- energy_result$p.value
        pval_matrix[j, i] <- energy_result$p.value
      }
    }
  }

  # Calculate q-values for multiple testing correction (same as KS function)
  if (qval) {
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
        warning("Q-value calculation failed: ", e$message, ". Returning placeholder.")
        qval_matrix <- "Q-value calculation failed - see warning message"
      })
    } else {
      qval_matrix <- "No valid p-values for q-value calculation"
    }
  } else {
    qval_matrix <- "Q-value calculation skipped (qval = FALSE)"
  }

  # Set row and column names
  rownames(energy_matrix) <- matrix_names
  colnames(energy_matrix) <- matrix_names
  rownames(pval_matrix) <- matrix_names
  colnames(pval_matrix) <- matrix_names

  # Only set names if qval_matrix is actually a matrix
  if (is.matrix(qval_matrix)) {
    rownames(qval_matrix) <- matrix_names
    colnames(qval_matrix) <- matrix_names
  }

  # Return results as a list
  return(list(
    diagonals = diagonals,
    energy_distances = energy_matrix,  # Energy test statistics
    pvalues = pval_matrix,            # Bootstrap p-values
    qvalues = qval_matrix             # Multiple testing corrected q-values
  ))
}

# Flexible heatmap plotting function for distance/difference matrices
plot_cov_heatmap <- function(analysis_result,
                             metric = "auto",           # Auto-detect or specify metric to plot
                             title = NULL,              # Custom plot title
                             color_scheme = "default",  # Color palette choice
                             show_values = FALSE,       # Show numeric values in cells
                             cex_values = 0.8,         # Text size for cell values
                             margins = c(5, 5),        # Plot margins
                             ...) {                    # Additional arguments to heatmap()

  # Auto-detect the appropriate matrix to plot based on available results
  if (metric == "auto") {
    if ("pairwise_differences" %in% names(analysis_result)) {
      plot_matrix <- analysis_result$pairwise_differences
      default_title <- "Pairwise Differences in Diagonal Means"
    } else if ("variance_differences" %in% names(analysis_result)) {
      plot_matrix <- analysis_result$variance_differences
      default_title <- "Pairwise Variance Distances"
    }else if ("wasserstein_distances" %in% names(analysis_result)) {
      plot_matrix <- analysis_result$wasserstein_distances
      default_title <- "Pairwise Wasserstein Distances"
    } else if ("ks_distances" %in% names(analysis_result)) {
      plot_matrix <- analysis_result$ks_distances
      default_title <- "Kolmogorov-Smirnov Test Statistics"
    } else if ("energy_distances" %in% names(analysis_result)) {
      plot_matrix <- analysis_result$energy_distances
      default_title <- "Energy Distance Statistics"
    } else {
      stop("Cannot auto-detect matrix type. Please specify 'metric' parameter.")
    }
  } else {
    # Manual metric selection
    valid_metrics <- c("pairwise_differences", "variance_differences", "wasserstein_distances",
                       "ks_distances", "energy_distances", "pvalues", "qvalues")

    if (!metric %in% valid_metrics) {
      stop("Invalid metric. Choose from: ", paste(valid_metrics, collapse = ", "))
    }

    if (!metric %in% names(analysis_result)) {
      stop("Metric '", metric, "' not found in analysis result.")
    }

    plot_matrix <- analysis_result[[metric]]

    # Handle non-matrix results (e.g., failed q-value calculations)
    if (!is.matrix(plot_matrix)) {
      stop("Selected metric is not a matrix: ", plot_matrix)
    }

    # Set appropriate default titles for each metric type
    default_title <- switch(metric,
                            "pairwise_differences" = "Pairwise Differences in Diagonal Means",
                            "variance_differences" = "Pairwise Differences in Diagonal Variances",
                            "wasserstein_distances" = "Pairwise Wasserstein Distances",
                            "ks_distances" = "Kolmogorov-Smirnov Test Statistics",
                            "energy_distances" = "Energy Distance Statistics",
                            "pvalues" = "P-values",
                            "qvalues" = "Q-values"
    )
  }

  # Use provided title or default
  plot_title <- if (is.null(title)) default_title else title

  # Set up color palette based on scheme choice
  if (color_scheme == "default") {
    col_palette <- heat.colors(50)  # Traditional heat colors (red-yellow-white)
  } else if (color_scheme == "blue_red") {
    col_palette <- colorRampPalette(c("blue", "white", "red"))(50)  # Diverging scheme
  } else if (color_scheme == "viridis") {
    if (requireNamespace("viridis", quietly = TRUE)) {
      col_palette <- viridis::viridis(50)  # Perceptually uniform colors
    } else {
      col_palette <- heat.colors(50)
      warning("viridis package not available, using default colors")
    }
  } else if (is.character(color_scheme) && length(color_scheme) > 1) {
    col_palette <- colorRampPalette(color_scheme)(50)  # Custom color vector
  } else {
    col_palette <- heat.colors(50)  # Fallback to default
  }

  # Create the heatmap with hierarchical clustering
  heatmap(plot_matrix,
          symm = TRUE,          # Matrix is symmetric
          main = plot_title,
          xlab = "Matrix",
          ylab = "Matrix",
          col = col_palette,
          margins = margins,
          ...)                  # Pass additional arguments to heatmap()

  # Add numeric values to cells if requested (approximate positioning)
  if (show_values) {
    # Get the reordered matrix from heatmap (this is tricky with base heatmap)
    # For simplicity, we'll add values to the original order
    # Note: This won't match the dendogram reordering, but gives an idea
    n <- nrow(plot_matrix)
    for (i in 1:n) {
      for (j in 1:n) {
        # Format small values in scientific notation
        if (plot_matrix[i,j] < 0.001) {
          text_val <- format(plot_matrix[i,j], scientific = TRUE, digits = 2)
        } else {
          text_val <- format(round(plot_matrix[i,j], 3), nsmall = 3)
        }
        # This is approximate positioning - exact positioning would require
        # extracting heatmap coordinates
        text(j/n, (n-i+1)/n, text_val, cex = cex_values, col = "black")
      }
    }
  }
}

# Specialized heatmap function for p-values and q-values with significance-based coloring
plot_pvalue_heatmap <- function(analysis_result,
                                value_type = "pvalues",    # "pvalues" or "qvalues"
                                title = NULL,              # Custom plot title
                                color_scheme = "pvalue",   # Color scheme optimized for p-values
                                show_values = TRUE,        # Show p-values/q-values in cells
                                cex_values = 0.8,         # Text size for cell values
                                margins = c(5, 5),        # Plot margins
                                ...) {                    # Additional arguments to heatmap()

  # Validate value_type parameter
  if (!value_type %in% c("pvalues", "qvalues")) {
    stop("value_type must be either 'pvalues' or 'qvalues'")
  }

  # Check if the requested values exist in the analysis result
  if (!value_type %in% names(analysis_result)) {
    stop("'", value_type, "' not found in analysis result")
  }

  plot_matrix <- analysis_result[[value_type]]

  # Handle cases where q-values couldn't be calculated
  if (!is.matrix(plot_matrix)) {
    if (is.character(plot_matrix)) {
      stop("Cannot plot ", value_type, ": ", plot_matrix)
    } else {
      stop("'", value_type, "' is not a matrix")
    }
  }

  # Set appropriate default title based on value type
  if (is.null(title)) {
    if (value_type == "pvalues") {
      title <- "Statistical Test P-values"
    } else {
      title <- "Multiple Testing Corrected Q-values"
    }

    # Add test type information if we can determine it from the result
    if ("ks_distances" %in% names(analysis_result)) {
      title <- paste("Kolmogorov-Smirnov", title)
    } else if ("energy_distances" %in% names(analysis_result)) {
      title <- paste("Energy Test", title)
    }
  }

  # Set up color schemes optimized for statistical significance visualization
  if (color_scheme == "pvalue") {
    # Red (significant) to white (non-significant) - emphasizes low p-values
    col_palette <- colorRampPalette(c("red", "orange", "yellow", "white", "lightblue", "blue"))(100)
  } else if (color_scheme == "significance") {
    # Simple red/white scheme based on significance levels
    col_palette <- colorRampPalette(c("red", "white"))(100)
  } else if (color_scheme == "grayscale") {
    col_palette <- colorRampPalette(c("black", "white"))(100)
  } else if (color_scheme == "viridis") {
    if (requireNamespace("viridis", quietly = TRUE)) {
      col_palette <- viridis::viridis(100, direction = -1)  # Reverse so low p-values are dark
    } else {
      col_palette <- colorRampPalette(c("red", "orange", "yellow", "white", "lightblue", "blue"))(100)
      warning("viridis package not available, using default p-value colors")
    }
  } else if (is.character(color_scheme) && length(color_scheme) > 1) {
    col_palette <- colorRampPalette(color_scheme)(100)  # Custom colors
  } else {
    col_palette <- colorRampPalette(c("red", "orange", "yellow", "white", "lightblue", "blue"))(100)
  }

  # Create the heatmap
  hm_result <- heatmap(plot_matrix,
                       symm = TRUE,
                       main = title,
                       xlab = "Matrix",
                       ylab = "Matrix",
                       col = col_palette,
                       margins = margins,
                       ...)

  # Add p-values/q-values to cells if requested
  if (show_values) {
    n <- nrow(plot_matrix)

    for (i in 1:n) {
      for (j in 1:n) {
        val <- plot_matrix[i, j]

        # Format the value appropriately (scientific notation for very small values)
        if (val < 0.001) {
          val_text <- format(val, scientific = TRUE, digits = 2)
        } else {
          val_text <- format(round(val, 3), nsmall = 3)
        }

        # Determine text color based on background intensity (white for dark backgrounds)
        text_color <- if (val < 0.1) "white" else "black"

        # Add the value to the plot (approximate positioning)
        text((j-0.5)/n, (n-i+0.5)/n, val_text, cex = cex_values, col = text_color)
      }
    }
  }

  # Return heatmap result invisibly for potential further manipulation
  invisible(hm_result)
}
