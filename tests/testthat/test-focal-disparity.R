library(testthat)
library(ape)

make_focal_residual_test_data <- function() {
  set.seed(2026)
  tree <- ape::rtree(12)
  size <- stats::rnorm(12)

  dat <- data.frame(
    y1 = 0.8 * size + stats::rnorm(12, sd = 0.3),
    y2 = -0.4 * size + stats::rnorm(12, sd = 0.35),
    size = size
  )
  rownames(dat) <- tree$tip.label

  list(
    tree = tree,
    dat = dat,
    Y = as.matrix(dat[, c("y1", "y2"), drop = FALSE])
  )
}

pick_internal_node_with_min_tips <- function(tree, min_tips = 3L) {
  candidates <- seq.int(ape::Ntip(tree) + 1L, ape::Ntip(tree) + ape::Nnode(tree))
  for (node in candidates) {
    clade_tips <- ape::extract.clade(tree, node = node)$tip.label
    if (length(clade_tips) >= min_tips && length(clade_tips) < ape::Ntip(tree)) {
      return(node)
    }
  }
  candidates[1L]
}

test_that("focal helper validates logical, numeric, and unsupported selectors", {
  tree <- ape::rtree(5)

  expect_identical(
    .bifrost_resolve_focal_tips(tree, tree$tip.label[c(1, 1, 3)]),
    tree$tip.label[c(1, 3)]
  )
  expect_identical(
    .bifrost_resolve_focal_tips(tree, c(TRUE, FALSE, TRUE, FALSE, FALSE)),
    tree$tip.label[c(1, 3)]
  )
  expect_error(
    .bifrost_resolve_focal_tips(tree, c(TRUE, FALSE)),
    "length equal to ape::Ntip"
  )
  expect_error(
    .bifrost_resolve_focal_tips(tree, 1),
    "not tip indices"
  )
  expect_error(
    .bifrost_resolve_focal_tips(tree, list(tree$tip.label[1])),
    "focal must be"
  )
})

test_that("alignment helper covers coercion, pruning, and validation errors", {
  set.seed(7)
  hc <- stats::hclust(stats::dist(matrix(stats::rnorm(25), ncol = 5)))
  tree <- ape::as.phylo(hc)
  dat <- data.frame(y1 = stats::rnorm(5), y2 = stats::rnorm(5))
  rownames(dat) <- tree$tip.label

  aligned <- .bifrost_align_focal_test_inputs(
    tree = hc,
    trait_data = dat[1:4, , drop = FALSE],
    focal = tree$tip.label[1:3],
    min_focal_tips = 2
  )

  expect_s3_class(aligned$tree, "phylo")
  expect_identical(ape::Ntip(aligned$tree), 4L)
  expect_identical(rownames(aligned$trait_data), aligned$tree$tip.label)

  expect_error(
    .bifrost_align_focal_test_inputs(tree = tree, trait_data = 1:3, focal = tree$tip.label[1:2], min_focal_tips = 2),
    "matrix or data.frame"
  )

  dat_no_rownames <- as.matrix(dat)
  rownames(dat_no_rownames) <- NULL
  expect_error(
    .bifrost_align_focal_test_inputs(tree = tree, trait_data = dat_no_rownames, focal = tree$tip.label[1:2], min_focal_tips = 2),
    "row names"
  )

  expect_error(
    .bifrost_align_focal_test_inputs(tree = tree, trait_data = dat, focal = tree$tip.label[1:2], min_focal_tips = 1),
    "integer >= 2"
  )

  expect_error(
    .bifrost_align_focal_test_inputs(
      tree = tree,
      trait_data = dat[1:2, , drop = FALSE],
      focal = tree$tip.label[1:2],
      min_focal_tips = 2
    ),
    "Fewer than 3 overlapping taxa"
  )

  expect_error(
    .bifrost_align_focal_test_inputs(
      tree = tree,
      trait_data = dat[3:5, , drop = FALSE],
      focal = tree$tip.label[1:2],
      min_focal_tips = 2
    ),
    "Need at least 2 focal tips"
  )

  expect_error(
    .bifrost_align_focal_test_inputs(tree = "bad-tree", trait_data = dat, focal = rownames(dat)[1:2], min_focal_tips = 2),
    "coercible to class 'phylo'"
  )
})

test_that("metric helpers cover built-in variants and invalid return values", {
  X <- matrix(
    c(0, 0,
      1, 2,
      3, 1),
    ncol = 2,
    byrow = TRUE
  )

  trace_metric <- .bifrost_make_disparity_metric("trace_cov")
  pairwise_metric <- .bifrost_make_disparity_metric("mean_pairwise_sqdist")
  log_det_metric <- .bifrost_make_disparity_metric("log_det_cov")
  max_eigen_metric <- .bifrost_make_disparity_metric("max_eigenvalue")
  evenness_metric <- .bifrost_make_disparity_metric("eigenvalue_evenness")
  singular_X <- cbind(1:4, 2 * (1:4), 3 * (1:4))
  singular_eigenvalues <- .bifrost_cov_eigenvalues(singular_X)

  expect_identical(trace_metric$label, "trace_cov")
  expect_identical(pairwise_metric$label, "mean_pairwise_sqdist")
  expect_identical(log_det_metric$label, "log_det_cov")
  expect_identical(max_eigen_metric$label, "max_eigenvalue")
  expect_identical(evenness_metric$label, "eigenvalue_evenness")
  expect_true(trace_metric$fun(X) > 0)
  expect_true(pairwise_metric$fun(X) > 0)
  expect_true(is.finite(log_det_metric$fun(X)))
  expect_true(is.finite(log_det_metric$fun(singular_X)))
  expect_true(max_eigen_metric$fun(X) > 0)
  expect_gte(evenness_metric$fun(X), 0)
  expect_lte(evenness_metric$fun(X), 1)
  expect_identical(length(singular_eigenvalues), 3L)
  expect_true(all(singular_eigenvalues >= 0))

  expect_error(
    .bifrost_eval_disparity_metric(function(x) c(1, 2), X),
    "single finite numeric value"
  )
  expect_error(
    .bifrost_eval_disparity_metric(function(x) NA_real_, X),
    "single finite numeric value"
  )
})

test_that("testFocalDisparity supports additional covariance-based metrics", {
  skip_if_not_installed("mvMORPH")
  skip_if_not_installed("phytools")

  x <- make_focal_residual_test_data()
  metrics <- c("log_det_cov", "max_eigenvalue", "eigenvalue_evenness")

  for (metric_name in metrics) {
    res <- testFocalDisparity(
      tree = x$tree,
      trait_data = x$dat,
      focal = x$tree$tip.label[1:4],
      formula = cbind(y1, y2) ~ size,
      metric = metric_name,
      nsim = 2,
      method = "LL",
      workers = 1,
      future_plan = "sequential",
      show_progress = FALSE,
      seed = 2026
    )

    expect_identical(res$metric, metric_name)
    expect_true(is.finite(res$observed_statistic))
    expect_true(all(is.finite(res$null_statistics)))
  }
})

test_that("runFocalDisparityGrid evaluates combinations and optionally writes outputs", {
  skip_if_not_installed("mvMORPH")

  ex <- focalDisparityExampleData(seed = 2)
  out_dir <- file.path(tempdir(), "focal-disparity-grid-test")

  grid <- runFocalDisparityGrid(
    data_bundle = list(
      tree = ex$tree,
      data = ex$trait_data,
      focal_tips = ex$focal_tips
    ),
    predictor_col = 4,
    response_cols = 1:3,
    metrics = c("centroid_msd", "mean_pairwise_sqdist"),
    models = c("BM", "OU"),
    nsim = 1,
    method = "LL",
    workers = 1,
    seed = 2,
    future_plan = "sequential",
    show_progress = FALSE,
    output_dir = out_dir
  )

  expect_s3_class(grid$formula, "formula")
  expect_s3_class(grid, "bifrost_focal_disparity_grid")
  expect_identical(grid$predictor, "size")
  expect_identical(grid$responses, c("y1", "y2", "y3"))
  expect_equal(nrow(grid$summary), 4L)
  expect_length(grid$results, 4L)
  expect_true(all(c("model", "metric", "p_upper") %in% colnames(grid$summary)))
  expect_true(file.exists(file.path(out_dir, "focal-disparity-grid-summary.csv")))
  expect_true(file.exists(file.path(out_dir, "focal-disparity-grid-results.rds")))
  expect_true(file.exists(file.path(out_dir, "focal-disparity-grid-summary.png")))
  expect_length(list.files(out_dir, pattern = "\\.png$", full.names = TRUE), 5L)

  printed <- paste(capture.output(print(grid)), collapse = "\n")
  expect_match(printed, "Focal Disparity Grid")

  plot_file <- tempfile(fileext = ".pdf")
  grDevices::pdf(plot_file)
  plot(grid)
  grDevices::dev.off()
  expect_true(file.exists(plot_file))
})

test_that("plotFocalResidualPCA returns residual ordination data and draws a plot", {
  skip_if_not_installed("mvMORPH")

  ex <- focalDisparityExampleData(seed = 2)
  plot_file <- tempfile(fileext = ".pdf")

  grDevices::pdf(plot_file)
  pca_plot <- plotFocalResidualPCA(
    tree = ex$tree,
    trait_data = ex$trait_data,
    focal = ex$focal_tips,
    formula = cbind(y1, y2, y3) ~ size,
    model = "OU",
    method = "LL"
  )
  grDevices::dev.off()

  expect_true(file.exists(plot_file))
  expect_s3_class(pca_plot$pca, "prcomp")
  expect_equal(nrow(pca_plot$scores), ape::Ntip(ex$tree))
  expect_equal(length(pca_plot$focal_indicator), ape::Ntip(ex$tree))
  expect_true(all(pca_plot$focal_indicator == (rownames(pca_plot$scores) %in% ex$focal_tips)))
  expect_equal(length(pca_plot$variance_explained), ncol(ex$trait_data) - 1L)
  expect_identical(pca_plot$fit_model, "OU")
})

test_that("plotFocalResidualPCA supports multi-model panel plots", {
  skip_if_not_installed("mvMORPH")

  ex <- focalDisparityExampleData(seed = 2)
  plot_file <- tempfile(fileext = ".pdf")

  grDevices::pdf(plot_file)
  pca_panels <- plotFocalResidualPCA(
    tree = ex$tree,
    trait_data = ex$trait_data,
    focal = ex$focal_tips,
    formula = cbind(y1, y2, y3) ~ size,
    model = c("BM", "OU", "EB"),
    method = "LL",
    show_legend = FALSE
  )
  grDevices::dev.off()

  expect_true(file.exists(plot_file))
  expect_type(pca_panels, "list")
  expect_identical(names(pca_panels), c("BM", "OU", "EB"))
  expect_true(all(vapply(pca_panels, function(x) inherits(x$pca, "prcomp"), logical(1))))
  expect_identical(unname(vapply(pca_panels, `[[`, character(1), "fit_model")), c("BM", "OU", "EB"))
})

test_that("plotSubtreeCovarianceHeatmaps returns post-hoc subtree matrices", {
  skip_if_not_installed("mvMORPH")

  ex <- focalDisparityExampleData(seed = 2)
  plot_file <- tempfile(fileext = ".pdf")

  grDevices::pdf(plot_file)
  subtree_mats <- suppressWarnings(
    plotSubtreeCovarianceHeatmaps(
      tree = ex$tree,
      trait_data = ex$trait_data,
      focal = ex$focal_tips,
      formula = cbind(y1, y2, y3) ~ size,
      model = c("BM", "OU"),
      subsets = c("focal", "background"),
      matrix_type = c("covariance", "correlation"),
      method = "LL"
    )
  )
  grDevices::dev.off()

  expect_true(file.exists(plot_file))
  expect_s3_class(subtree_mats, "bifrost_subtree_covariance_heatmaps")
  expect_identical(names(subtree_mats), c("BM", "OU"))
  expect_identical(names(subtree_mats$BM), c("focal", "background"))
  expect_equal(colnames(subtree_mats$BM$focal$covariance), c("y1", "y2", "y3"))
  expect_equal(colnames(subtree_mats$BM$focal$correlation), c("y1", "y2", "y3"))
  expect_true(isSymmetric(subtree_mats$OU$background$covariance))
  expect_true(isSymmetric(subtree_mats$OU$background$correlation))
  expect_equal(unname(diag(subtree_mats$OU$background$correlation)), rep(1, 3))
  expect_true("summary" %in% names(subtree_mats$BM$focal))
  expect_equal(subtree_mats$BM$focal$summary$subset, "focal")

  summary_obj <- bifrost:::summary.bifrost_subtree_covariance_heatmaps(subtree_mats)
  expect_true(all(c("summary_table", "contrast_table") %in% names(summary_obj)))
  expect_equal(unique(summary_obj$summary_table$model), c("BM", "OU"))
  expect_equal(summary_obj$summary_table$subset, c("focal", "background", "focal", "background"))
  expect_equal(summary_obj$contrast_table$model, c("BM", "OU"))
  printed <- paste(capture.output(print(subtree_mats)), collapse = "\n")
  expect_match(printed, "Subtree Covariance Heatmaps")
  expect_match(printed, "Focal vs Background Contrasts")
})

test_that("sigma and model helpers cover alternative object shapes", {
  sigma <- diag(c(1, 2))

  expect_identical(.bifrost_extract_sigma_matrix(list(sigma = sigma)), sigma)
  expect_identical(.bifrost_extract_sigma_matrix(list(sigma = list(S = sigma))), sigma)
  expect_identical(.bifrost_extract_sigma_matrix(list(sigma = list(Pinv = sigma))), sigma)
  expect_error(
    .bifrost_extract_sigma_matrix(list(sigma = NULL)),
    "Could not locate a sigma estimate"
  )
  expect_error(
    .bifrost_extract_sigma_matrix(list(sigma = list(foo = sigma))),
    "usable sigma matrix"
  )

  expect_identical(.bifrost_resolve_simulation_model("BM"), "BM1")
  expect_error(
    .bifrost_resolve_simulation_model("BMM"),
    "model must be one of"
  )
})

test_that("OU and EB expansion helpers cover matrix, vector, scalar, and invalid inputs", {
  alpha_mat <- diag(c(0.2, 0.4))
  beta_mat <- diag(c(-0.3, -0.1))

  expect_identical(.bifrost_expand_ou_alpha(alpha_mat, 2L), alpha_mat)
  expect_identical(.bifrost_expand_ou_alpha(0.5, 2L), diag(0.5, 2))
  expect_identical(.bifrost_expand_ou_alpha(c(0.2, 0.4), 2L), alpha_mat)
  expect_error(
    .bifrost_expand_ou_alpha(diag(3), 2L),
    "alpha matrix with dimensions"
  )
  expect_error(
    .bifrost_expand_ou_alpha(c(0.1, 0.2, 0.3), 2L),
    "scalar, length-n_traits vector, or square alpha matrix"
  )

  expect_identical(.bifrost_expand_eb_beta(beta_mat, 2L), beta_mat)
  expect_identical(.bifrost_expand_eb_beta(-0.2, 2L), -0.2)
  expect_identical(.bifrost_expand_eb_beta(c(-0.3, -0.1), 2L), beta_mat)
  expect_error(
    .bifrost_expand_eb_beta(diag(3), 2L),
    "beta matrix with dimensions"
  )
  expect_error(
    .bifrost_expand_eb_beta(c(-0.2, -0.1, 0), 2L),
    "scalar, length-n_traits vector, or square beta matrix"
  )
})

test_that("null-simulation helpers cover OU and EB parameter errors and extraction variants", {
  sigma <- diag(c(1, 2))

  ou_params <- .bifrost_build_null_simulation_params(
    fit = list(param = c(0.2, 0.4)),
    simulation_model = "OU1",
    n_traits = 2,
    sigma_hat = sigma
  )
  eb_params <- .bifrost_build_null_simulation_params(
    fit = list(param = c(-0.3, -0.1)),
    simulation_model = "EB",
    n_traits = 2,
    sigma_hat = sigma
  )
  bm_params <- .bifrost_build_null_simulation_params(
    fit = list(param = NULL),
    simulation_model = "BM1",
    n_traits = 2,
    sigma_hat = sigma
  )

  expect_true(all(c("ntraits", "sigma", "theta", "alpha") %in% names(ou_params)))
  expect_true(all(c("ntraits", "sigma", "theta", "beta") %in% names(eb_params)))
  expect_identical(names(bm_params), c("ntraits", "sigma", "theta"))

  expect_error(
    .bifrost_build_null_simulation_params(list(param = NULL), "OU1", 2, sigma),
    "Could not extract OU alpha"
  )
  expect_error(
    .bifrost_build_null_simulation_params(list(param = NA_real_), "EB", 2, sigma),
    "Could not extract EB beta"
  )

  arr <- array(1:8, dim = c(2, 2, 2))
  sim_from_list <- .bifrost_extract_simulated_matrix(
    simulated = list(arr[, , 1, drop = FALSE]),
    tip_labels = c("t1", "t2"),
    trait_names = c("a", "b")
  )
  sim_from_array <- .bifrost_extract_simulated_matrix(
    simulated = arr[, , 1, drop = FALSE],
    tip_labels = c("t1", "t2"),
    trait_names = c("a", "b")
  )

  expect_identical(dim(sim_from_list), c(2L, 2L))
  expect_identical(dim(sim_from_array), c(2L, 2L))
  expect_identical(rownames(sim_from_array), c("t1", "t2"))
  expect_identical(colnames(sim_from_array), c("a", "b"))

  expect_identical(
    .bifrost_get_response_names(
      fitted_values = matrix(1:4, ncol = 2, dimnames = list(NULL, c("y1", "y2"))),
      residuals_matrix = matrix(1:4, ncol = 2)
    ),
    c("y1", "y2")
  )
  expect_identical(
    .bifrost_get_response_names(
      fitted_values = matrix(1:4, ncol = 2),
      residuals_matrix = matrix(1:4, ncol = 2, dimnames = list(NULL, c("r1", "r2")))
    ),
    c("r1", "r2")
  )
  expect_error(
    .bifrost_get_response_names(
      fitted_values = matrix(1:4, ncol = 2),
      residuals_matrix = matrix(1:4, ncol = 2)
    ),
    "Could not determine response column names"
  )
})

test_that("bootstrap helper covers sequential and future execution paths", {
  progressr::handlers("void")
  sim_one <- function(i, progress = NULL) {
    if (!is.null(progress)) {
      progress(sprintf("replicate %d", i))
    }
    i * 2
  }

  old_seed_exists <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (old_seed_exists) {
    old_seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  }
  if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
    rm(".Random.seed", envir = .GlobalEnv)
  }
  on.exit({
    if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
    if (old_seed_exists) {
      assign(".Random.seed", old_seed, envir = .GlobalEnv)
    }
  }, add = TRUE)

  expect_identical(
    .bifrost_run_bootstrap(3, sim_one, seed = 99, show_progress = FALSE, use_future = FALSE),
    c(2, 4, 6)
  )
  expect_false(exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))

  set.seed(123)
  current_seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  expect_identical(
    .bifrost_run_bootstrap(2, sim_one, seed = 99, show_progress = TRUE, use_future = FALSE),
    c(2, 4)
  )
  expect_identical(get(".Random.seed", envir = .GlobalEnv, inherits = FALSE), current_seed)

  old_plan <- future::plan()
  on.exit(future::plan(old_plan), add = TRUE)
  future::plan(future::sequential)

  expect_identical(
    .bifrost_run_bootstrap(2, sim_one, seed = 99, show_progress = FALSE, use_future = TRUE),
    c(2, 4)
  )
  expect_identical(
    .bifrost_run_bootstrap(2, sim_one, seed = 99, show_progress = TRUE, use_future = TRUE),
    c(2, 4)
  )
})

test_that("example helper falls back to the nearest clade size when needed", {
  tree <- ape::rtree(10)
  node <- .bifrost_find_example_focal_node(tree, min_tips = 20L, max_tips = 25L)

  expect_true(node > ape::Ntip(tree))
  expect_true(node <= ape::Ntip(tree) + ape::Nnode(tree))
})

test_that("testFocalDisparity fits a named-column formula workflow", {
  skip_if_not_installed("mvMORPH")
  skip_if_not_installed("phytools")

  x <- make_focal_residual_test_data()

  res <- testFocalDisparity(
    tree = x$tree,
    trait_data = x$dat,
    focal = x$tree$tip.label[1:4],
    formula = cbind(y1, y2) ~ size,
    nsim = 5,
    method = "LL",
    workers = 1,
    future_plan = "sequential",
    show_progress = FALSE
  )

  expect_s3_class(res, "bifrost_focal_disparity_test")
  expect_length(res$null_statistics, 5L)
  expect_true(is.numeric(res$observed_statistic))
  expect_true(all(res$null_statistics >= 0))
  expect_identical(res$metric, "centroid_msd")
  expect_identical(res$fit_model, "BM")
  expect_identical(res$simulation_model, "BM1")
  expect_true(res$p_upper >= 0 && res$p_upper <= 1)
  expect_true(res$p_lower >= 0 && res$p_lower <= 1)
  expect_true(res$p_two_tailed >= 0 && res$p_two_tailed <= 1)

  printed <- capture.output(print(res))
  expect_true(any(grepl("Bifrost Focal Disparity Test", printed, fixed = TRUE)))
  expect_true(any(grepl("Metric: centroid_msd", printed, fixed = TRUE)))
})

test_that("testFocalDisparity supports legacy indexed formulas over trait_data", {
  skip_if_not_installed("mvMORPH")
  skip_if_not_installed("phytools")

  x <- make_focal_residual_test_data()

  res <- testFocalDisparity(
    tree = x$tree,
    trait_data = x$dat,
    focal = x$tree$tip.label[1:4],
    formula = "trait_data[, 1:2] ~ trait_data[, 3]",
    nsim = 3,
    method = "LL",
    workers = 1,
    future_plan = "sequential",
    show_progress = FALSE,
    seed = 2026
  )

  expect_s3_class(res, "bifrost_focal_disparity_test")
  expect_identical(res$formula_input, "trait_data[, 1:2] ~ trait_data[, 3]")
  expect_identical(res$fit_formula, "cbind(y1, y2) ~ size")
  expect_true(all(is.finite(res$null_statistics)))
})

test_that("testFocalDisparity resolves focal clades from node numbers", {
  skip_if_not_installed("mvMORPH")
  skip_if_not_installed("phytools")

  x <- make_focal_residual_test_data()
  focal_node <- pick_internal_node_with_min_tips(x$tree, min_tips = 3L)
  expected_tips <- ape::extract.clade(x$tree, node = focal_node)$tip.label

  res <- testFocalDisparity(
    tree = x$tree,
    trait_data = x$dat,
    focal = focal_node,
    formula = cbind(y1, y2) ~ size,
    nsim = 4,
    method = "LL",
    workers = 1,
    future_plan = "sequential",
    show_progress = FALSE
  )

  expect_setequal(res$focal_tips, expected_tips)
})

test_that("testFocalDisparity accepts custom metric functions", {
  skip_if_not_installed("mvMORPH")
  skip_if_not_installed("phytools")

  x <- make_focal_residual_test_data()
  custom_stat <- function(X) max(rowSums(X^2))

  res <- testFocalDisparity(
    tree = x$tree,
    trait_data = x$Y,
    focal = c(TRUE, TRUE, TRUE, rep(FALSE, ape::Ntip(x$tree) - 3L)),
    formula = "trait_data ~ 1",
    metric = custom_stat,
    nsim = 3,
    method = "LL",
    workers = 1,
    future_plan = "sequential",
    show_progress = FALSE
  )

  expect_identical(res$metric, "user_function")
  expect_length(res$null_statistics, 3L)
  expect_true(is.numeric(res$observed_statistic))
})

test_that("testFocalDisparity handles OU fits with factor and numeric predictors", {
  skip_if_not_installed("mvMORPH")
  skip_if_not_installed("phytools")

  set.seed(11)
  tree <- ape::rtree(18)
  size <- stats::rnorm(18)
  grp <- factor(rep(c("a", "b", "c"), each = 6))
  dat <- data.frame(
    y1 = 0.5 * size + rep(c(-0.2, 0.1, 0.3), each = 6) + stats::rnorm(18, sd = 0.2),
    y2 = -0.3 * size + rep(c(0.1, -0.15, 0.2), each = 6) + stats::rnorm(18, sd = 0.25),
    size = size,
    grp = grp
  )
  rownames(dat) <- tree$tip.label

  res <- testFocalDisparity(
    tree = tree,
    trait_data = dat,
    focal = tree$tip.label[1:6],
    formula = cbind(y1, y2) ~ size + grp,
    nsim = 3,
    model = "OU",
    method = "LL",
    workers = 1,
    future_plan = "sequential",
    show_progress = FALSE,
    seed = 11
  )

  expect_s3_class(res, "bifrost_focal_disparity_test")
  expect_identical(res$fit_model, "OU")
  expect_identical(res$simulation_model, "OU1")
  expect_true(all(c("size", "grpb", "grpc") %in% colnames(res$fit_data)))
  expect_true(all(is.finite(res$null_statistics)))
})

test_that("testFocalDisparity stores fit-ready formula/data for factor predictors", {
  skip_if_not_installed("mvMORPH")
  skip_if_not_installed("phytools")

  set.seed(1)
  tree <- ape::rtree(12)
  dat <- data.frame(
    y1 = stats::rnorm(12),
    y2 = stats::rnorm(12),
    grp = factor(rep(c("a", "b"), each = 6))
  )
  rownames(dat) <- tree$tip.label

  res <- testFocalDisparity(
    tree = tree,
    trait_data = dat,
    focal = tree$tip.label[1:4],
    formula = cbind(y1, y2) ~ grp,
    nsim = 2,
    method = "LL",
    workers = 1,
    future_plan = "sequential",
    show_progress = FALSE
  )

  expect_identical(res$formula_input, "cbind(y1, y2) ~ grp")
  expect_identical(res$fit_formula, "cbind(y1, y2) ~ grpb")
  expect_true("grpb" %in% colnames(res$fit_data))
  expect_true("grp" %in% colnames(res$trait_data))
})

test_that("testFocalDisparity rejects invalid internal node ids", {
  skip_if_not_installed("mvMORPH")
  skip_if_not_installed("phytools")

  x <- make_focal_residual_test_data()
  bad_node <- ape::Ntip(x$tree) + ape::Nnode(x$tree) + 1L

  expect_error(
    testFocalDisparity(
      tree = x$tree,
      trait_data = x$dat,
      focal = bad_node,
      formula = cbind(y1, y2) ~ size,
      nsim = 2,
      method = "LL",
      workers = 1,
      future_plan = "sequential",
      show_progress = FALSE
    ),
    "internal node numbers present in tree"
  )
})

test_that("focalDisparityExampleData returns usable example inputs", {
  skip_if_not_installed("mvMORPH")
  skip_if_not_installed("phytools")

  ex <- focalDisparityExampleData(seed = 2)
  expect_true(inherits(ex$tree, "phylo"))
  expect_true(is.data.frame(ex$trait_data))
  expect_true(length(ex$focal_tips) >= 3L)
  expect_true(ex$focal_node > ape::Ntip(ex$tree))

  res <- testFocalDisparity(
    tree = ex$tree,
    trait_data = ex$trait_data,
    focal = ex$focal_tips,
    formula = cbind(y1, y2, y3) ~ size,
    nsim = 5,
    method = "LL",
    workers = 1,
    future_plan = "sequential",
    show_progress = FALSE,
    seed = 2
  )

  expect_s3_class(res, "bifrost_focal_disparity_test")
})

test_that("testFocalDisparity produces a non-degenerate OU null on example data", {
  skip_if_not_installed("mvMORPH")
  skip_if_not_installed("phytools")

  ex <- focalDisparityExampleData(seed = 2)
  res <- testFocalDisparity(
    tree = ex$tree,
    trait_data = ex$trait_data,
    focal = ex$focal_tips,
    formula = cbind(y1, y2, y3) ~ size,
    metric = "centroid_msd",
    nsim = 5,
    model = "OU",
    method = "LL",
    workers = 1,
    future_plan = "sequential",
    show_progress = FALSE,
    seed = 2
  )

  expect_identical(res$simulation_model, "OU1")
  expect_gt(length(unique(signif(res$null_statistics, 12))), 1L)
  expect_true(is.finite(res$null_sd))
  expect_gt(res$null_sd, 0)
})

test_that("testFocalDisparity supports EB fits on example data", {
  skip_if_not_installed("mvMORPH")
  skip_if_not_installed("phytools")

  ex <- focalDisparityExampleData(seed = 2)
  res <- testFocalDisparity(
    tree = ex$tree,
    trait_data = ex$trait_data,
    focal = ex$focal_tips,
    formula = cbind(y1, y2, y3) ~ size,
    metric = "centroid_msd",
    nsim = 5,
    model = "EB",
    method = "LL",
    workers = 1,
    future_plan = "sequential",
    show_progress = FALSE,
    seed = 2
  )

  expect_s3_class(res, "bifrost_focal_disparity_test")
  expect_identical(res$simulation_model, "EB")
  expect_length(res$null_statistics, 5L)
  expect_true(all(is.finite(res$null_statistics)))
})

test_that("testFocalDisparity rejects unsupported fit models", {
  skip_if_not_installed("mvMORPH")
  skip_if_not_installed("phytools")

  ex <- focalDisparityExampleData(seed = 2)
  expect_error(
    testFocalDisparity(
      tree = ex$tree,
      trait_data = ex$trait_data,
      focal = ex$focal_tips,
      formula = cbind(y1, y2, y3) ~ size,
      metric = "centroid_msd",
      nsim = 5,
      model = "BMM",
      method = "LL",
      workers = 1,
      future_plan = "sequential",
      show_progress = FALSE,
      seed = 2
    ),
    "one of"
  )
})

test_that("testFocalDisparity covers validation and progress-related branches", {
  skip_if_not_installed("mvMORPH")
  skip_if_not_installed("phytools")

  x <- make_focal_residual_test_data()

  expect_error(
    do.call(
      testFocalDisparity,
      c(
        list(
          tree = x$tree,
          trait_data = x$dat,
          focal = x$tree$tip.label[1:4],
          formula = cbind(y1, y2) ~ size,
          nsim = 2,
          method = "LL",
          workers = 1,
          future_plan = "sequential",
          show_progress = FALSE
        ),
        list(data = x$dat)
      )
    ),
    "reserved mvgls arguments"
  )

  expect_error(
    testFocalDisparity(
      tree = x$tree,
      trait_data = x$dat,
      focal = x$tree$tip.label[1:4],
      formula = cbind(y1, y2) ~ size,
      nsim = 0,
      method = "LL",
      workers = 1,
      future_plan = "sequential",
      show_progress = FALSE
    ),
    "nsim must be"
  )

  expect_error(
    testFocalDisparity(
      tree = x$tree,
      trait_data = x$dat,
      focal = x$tree$tip.label[1:4],
      formula = cbind(y1, y2) ~ size,
      nsim = 1,
      method = "LL",
      workers = 0,
      future_plan = "sequential",
      show_progress = FALSE
    ),
    "workers must be"
  )

  res_default_plan <- testFocalDisparity(
    tree = x$tree,
    trait_data = x$dat,
    focal = x$tree$tip.label[1:4],
    formula = cbind(y1, y2) ~ size,
    nsim = 2,
    method = "LL",
    workers = 1,
    show_progress = FALSE
  )
  expect_s3_class(res_default_plan, "bifrost_focal_disparity_test")

  progressr::handlers("void")
  res_with_progress <- testFocalDisparity(
    tree = x$tree,
    trait_data = x$dat,
    focal = x$tree$tip.label[1:4],
    formula = cbind(y1, y2) ~ size,
    nsim = 2,
    method = "LL",
    workers = 1,
    future_plan = "sequential",
    show_progress = TRUE
  )
  expect_length(res_with_progress$null_statistics, 2L)

  testthat::local_mocked_bindings(
    normalizeMvglsFormulaCall = function(...) {
      list(
        formula = stats::as.formula("cbind(y1, y2) ~ 1"),
        args_list = list(data = NULL)
      )
    },
    .package = "bifrost"
  )
  expect_error(
    testFocalDisparity(
      tree = x$tree,
      trait_data = x$dat,
      focal = x$tree$tip.label[1:4],
      formula = cbind(y1, y2) ~ size,
      nsim = 1,
      method = "LL",
      workers = 1,
      future_plan = "sequential",
      show_progress = FALSE
    ),
    "Could not construct a data argument"
  )
})

test_that("testFocalDisparity supports OU fits with multisession bootstrap", {
  skip_on_cran()
  skip_if_not_installed("mvMORPH")
  skip_if_not_installed("phytools")
  skip_if(future::availableCores() < 2L)

  ex <- focalDisparityExampleData(seed = 2)
  res <- testFocalDisparity(
    tree = ex$tree,
    trait_data = ex$trait_data,
    focal = ex$focal_tips,
    formula = cbind(y1, y2, y3) ~ size,
    metric = "centroid_msd",
    nsim = 3,
    model = "OU",
    method = "LL",
    workers = 2,
    future_plan = "multisession",
    show_progress = FALSE,
    seed = 2
  )

  expect_s3_class(res, "bifrost_focal_disparity_test")
  expect_length(res$null_statistics, 3L)
  expect_true(all(is.finite(res$null_statistics)))
  expect_gt(res$null_sd, 0)
})

test_that("testFocalDisparity supports BM fits with multicore bootstrap", {
  skip_if_not_installed("mvMORPH")
  skip_if_not_installed("phytools")
  skip_on_os("windows")
  skip_if(future::availableCores() < 2L)

  ex <- focalDisparityExampleData(seed = 2)
  res <- testFocalDisparity(
    tree = ex$tree,
    trait_data = ex$trait_data,
    focal = ex$focal_tips,
    formula = cbind(y1, y2, y3) ~ size,
    metric = "centroid_msd",
    nsim = 2,
    model = "BM",
    method = "LL",
    workers = 2,
    future_plan = "multicore",
    show_progress = FALSE,
    seed = 2
  )

  expect_s3_class(res, "bifrost_focal_disparity_test")
  expect_length(res$null_statistics, 2L)
  expect_true(all(is.finite(res$null_statistics)))
})

test_that("plot method expands x-axis so the observed line remains visible", {
  skip_if_not_installed("mvMORPH")
  skip_if_not_installed("phytools")

  ex <- focalDisparityExampleData(seed = 2)
  res <- testFocalDisparity(
    tree = ex$tree,
    trait_data = ex$trait_data,
    focal = ex$focal_tips,
    formula = cbind(y1, y2, y3) ~ size,
    metric = "centroid_msd",
    nsim = 5,
    model = "BM",
    method = "LL",
    workers = 1,
    future_plan = "sequential",
    show_progress = FALSE,
    seed = 2
  )

  res$observed_statistic <- max(res$null_statistics) + 1
  tmp <- tempfile(fileext = ".pdf")
  grDevices::pdf(tmp)
  on.exit({
    grDevices::dev.off()
    unlink(tmp)
  }, add = TRUE)

  plot(res)
  usr <- graphics::par("usr")
  expect_lte(usr[1], res$observed_statistic)
  expect_gte(usr[2], res$observed_statistic)
})

test_that("print and plot methods cover normalized-formula and xlim edge branches", {
  skip_if_not_installed("mvMORPH")
  skip_if_not_installed("phytools")

  set.seed(1)
  tree <- ape::rtree(12)
  dat <- data.frame(
    y1 = stats::rnorm(12),
    y2 = stats::rnorm(12),
    grp = factor(rep(c("a", "b"), each = 6))
  )
  rownames(dat) <- tree$tip.label

  res <- testFocalDisparity(
    tree = tree,
    trait_data = dat,
    focal = tree$tip.label[1:4],
    formula = cbind(y1, y2) ~ grp,
    nsim = 2,
    method = "LL",
    workers = 1,
    future_plan = "sequential",
    show_progress = FALSE
  )

  printed <- capture.output(print(res))
  expect_true(any(grepl("Fit formula:", printed, fixed = TRUE)))

  custom_plot_obj <- structure(
    list(
      null_statistics = c(1, 1.5, 2),
      observed_statistic = 1,
      metric = "user_function"
    ),
    class = "bifrost_focal_disparity_test"
  )

  tmp <- tempfile(fileext = ".pdf")
  grDevices::pdf(tmp)
  on.exit({
    grDevices::dev.off()
    unlink(tmp)
  }, add = TRUE)

  expect_invisible(plot(custom_plot_obj, xlim = c(1, 1)))
  usr <- graphics::par("usr")
  expect_lte(usr[1], custom_plot_obj$observed_statistic)
  expect_gte(usr[2], custom_plot_obj$observed_statistic)

  suppressWarnings(
    expect_error(
      plot(custom_plot_obj, xlim = c(NA_real_, NA_real_)),
      "xlim must contain two finite numeric values"
    )
  )
})

test_that("plot method widens a user-supplied xlim when needed", {
  skip_if_not_installed("mvMORPH")
  skip_if_not_installed("phytools")

  ex <- focalDisparityExampleData(seed = 2)
  res <- testFocalDisparity(
    tree = ex$tree,
    trait_data = ex$trait_data,
    focal = ex$focal_tips,
    formula = cbind(y1, y2, y3) ~ size,
    metric = "centroid_msd",
    nsim = 5,
    model = "BM",
    method = "LL",
    workers = 1,
    future_plan = "sequential",
    show_progress = FALSE,
    seed = 2
  )

  forced_xmax <- max(res$null_statistics)
  res$observed_statistic <- forced_xmax + 0.75
  tmp <- tempfile(fileext = ".pdf")
  grDevices::pdf(tmp)
  on.exit({
    grDevices::dev.off()
    unlink(tmp)
  }, add = TRUE)

  plot(res, xlim = c(min(res$null_statistics), forced_xmax))
  usr <- graphics::par("usr")
  expect_lte(usr[1], res$observed_statistic)
  expect_gte(usr[2], res$observed_statistic)
})
