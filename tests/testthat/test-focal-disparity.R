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
  expect_identical(res$null_model, "BM")
  expect_identical(res$simulation_model, "BM1")
  expect_true(res$p_upper >= 0 && res$p_upper <= 1)
  expect_true(res$p_lower >= 0 && res$p_lower <= 1)
  expect_true(res$p_two_tailed >= 0 && res$p_two_tailed <= 1)

  printed <- capture.output(print(res))
  expect_true(any(grepl("Bifrost Focal Disparity Test", printed, fixed = TRUE)))
  expect_true(any(grepl("Metric: centroid_msd", printed, fixed = TRUE)))
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

test_that("testFocalDisparity accepts custom statistic functions", {
  skip_if_not_installed("mvMORPH")
  skip_if_not_installed("phytools")

  x <- make_focal_residual_test_data()
  custom_stat <- function(X) max(rowSums(X^2))

  res <- testFocalDisparity(
    tree = x$tree,
    trait_data = x$Y,
    focal = c(TRUE, TRUE, TRUE, rep(FALSE, ape::Ntip(x$tree) - 3L)),
    formula = "trait_data ~ 1",
    statistic = custom_stat,
    nsim = 3,
    method = "LL",
    workers = 1,
    future_plan = "sequential",
    show_progress = FALSE
  )

  expect_identical(res$statistic, "user_function")
  expect_identical(res$metric, "user_function")
  expect_length(res$null_statistics, 3L)
  expect_true(is.numeric(res$observed_statistic))
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

  expect_identical(res$null_model, "OU")
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
  expect_identical(res$null_model, "EB")
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
