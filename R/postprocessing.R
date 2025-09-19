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

  # Calculate q-values only if requested
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
    ks_distances = ks_matrix,
    pvalues = pval_matrix,
    qvalues = qval_matrix
  ))
}

# Energy distance version with p-values and q-values
energy_cov_diagonals <- function(object, R = 99, qval = TRUE) {
  # Check if energy package is available
  if (!requireNamespace("energy", quietly = TRUE)) {
    stop("Package 'energy' is required for energy distance calculations")
  }

  # Check if qvalue package is available (only if qval = TRUE)
  if (qval && !requireNamespace("qvalue", quietly = TRUE)) {
    stop("Package 'qvalue' is required for q-value calculations")
  }

  # Extract covariance matrices from object$VCV
  vcv_matrices <- object$VCV

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

  # Calculate q-values only if requested
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
    energy_distances = energy_matrix,
    pvalues = pval_matrix,
    qvalues = qval_matrix
  ))
}


plot_cov_heatmap <- function(analysis_result,
                             metric = "auto",
                             title = NULL,
                             color_scheme = "default",
                             show_values = FALSE,
                             cex_values = 0.8,
                             margins = c(5, 5),
                             ...) {

  # Auto-detect the appropriate matrix to plot
  if (metric == "auto") {
    if ("pairwise_differences" %in% names(analysis_result)) {
      plot_matrix <- analysis_result$pairwise_differences
      default_title <- "Pairwise Differences in Diagonal Means"
    } else if ("wasserstein_distances" %in% names(analysis_result)) {
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
    # Manual selection
    valid_metrics <- c("pairwise_differences", "wasserstein_distances",
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

    default_title <- switch(metric,
                            "pairwise_differences" = "Pairwise Differences in Diagonal Means",
                            "wasserstein_distances" = "Pairwise Wasserstein Distances",
                            "ks_distances" = "Kolmogorov-Smirnov Test Statistics",
                            "energy_distances" = "Energy Distance Statistics",
                            "pvalues" = "P-values",
                            "qvalues" = "Q-values"
    )
  }

  # Use provided title or default
  plot_title <- if (is.null(title)) default_title else title

  # Set up color scheme
  if (color_scheme == "default") {
    col_palette <- heat.colors(50)
  } else if (color_scheme == "blue_red") {
    col_palette <- colorRampPalette(c("blue", "white", "red"))(50)
  } else if (color_scheme == "viridis") {
    if (requireNamespace("viridis", quietly = TRUE)) {
      col_palette <- viridis::viridis(50)
    } else {
      col_palette <- heat.colors(50)
      warning("viridis package not available, using default colors")
    }
  } else if (is.character(color_scheme) && length(color_scheme) > 1) {
    col_palette <- colorRampPalette(color_scheme)(50)
  } else {
    col_palette <- heat.colors(50)
  }

  # Create the heatmap
  heatmap(plot_matrix,
          symm = TRUE,
          main = plot_title,
          xlab = "Matrix",
          ylab = "Matrix",
          col = col_palette,
          margins = margins,
          ...)

  # Add values to cells if requested
  if (show_values) {
    # Get the reordered matrix from heatmap (this is tricky with base heatmap)
    # For simplicity, we'll add values to the original order
    # Note: This won't match the dendogram reordering, but gives an idea
    n <- nrow(plot_matrix)
    for (i in 1:n) {
      for (j in 1:n) {
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


plot_pvalue_heatmap <- function(analysis_result,
                                value_type = "pvalues",
                                title = NULL,
                                color_scheme = "pvalue",
                                show_values = TRUE,
                                cex_values = 0.8,
                                margins = c(5, 5),
                                ...) {

  # Validate value_type
  if (!value_type %in% c("pvalues", "qvalues")) {
    stop("value_type must be either 'pvalues' or 'qvalues'")
  }

  # Check if the requested values exist
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

  # Set default title
  if (is.null(title)) {
    if (value_type == "pvalues") {
      title <- "Statistical Test P-values"
    } else {
      title <- "Multiple Testing Corrected Q-values"
    }

    # Add test type info if we can determine it
    if ("ks_distances" %in% names(analysis_result)) {
      title <- paste("Kolmogorov-Smirnov", title)
    } else if ("energy_distances" %in% names(analysis_result)) {
      title <- paste("Energy Test", title)
    }
  }

  # Set up color scheme for p-values/q-values
  if (color_scheme == "pvalue") {
    # Red (significant) to white (non-significant) to blue (very non-significant)
    col_palette <- colorRampPalette(c("red", "orange", "yellow", "white", "lightblue", "blue"))(100)
  } else if (color_scheme == "significance") {
    # Simple red/white scheme based on significance
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
    col_palette <- colorRampPalette(color_scheme)(100)
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

  # Add values if requested
  if (show_values) {
    n <- nrow(plot_matrix)

    for (i in 1:n) {
      for (j in 1:n) {
        val <- plot_matrix[i, j]

        # Format the value appropriately
        if (val < 0.001) {
          val_text <- format(val, scientific = TRUE, digits = 2)
        } else {
          val_text <- format(round(val, 3), nsmall = 3)
        }

        # Determine text color based on background
        text_color <- if (val < 0.1) "white" else "black"

        # Add the value to the plot (approximate positioning)
        text((j-0.5)/n, (n-i+0.5)/n, val_text, cex = cex_values, col = text_color)
      }
    }
  }

  invisible(hm_result)
}
