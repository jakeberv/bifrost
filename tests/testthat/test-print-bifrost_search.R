# tests/testthat/test-print-bifrost_search.R

testthat::skip_on_cran()

# Ensure helper exists regardless of test file order
if (!exists("skip_if_missing_deps", mode = "function")) {
  skip_if_missing_deps <- function() {
    testthat::skip_if_not_installed("ape")
    testthat::skip_if_not_installed("phytools")
    testthat::skip_if_not_installed("mvMORPH")
    testthat::skip_if_not_installed("future")
  }
}

# ---- Test XX: print method core output on a real no-shifts run ---------------
test_that("print.bifrost_search prints core sections on a real no-shifts run", {
  skip_if_missing_deps()

  set.seed(123)
  tr <- ape::rtree(20)
  X <- matrix(rnorm(20 * 2), ncol = 2)
  rownames(X) <- tr$tip.label

  res <- suppressWarnings(suppressMessages(searchOptimalConfiguration(
    baseline_tree              = tr,
    trait_data                 = X,
    formula                    = "trait_data ~ 1",
    min_descendant_tips        = 10,
    num_cores                  = 1,
    shift_acceptance_threshold = 1e9,
    plot                       = FALSE,
    IC                         = "GIC",
    store_model_fit_history    = FALSE,
    method                     = "LL",
    verbose                    = FALSE
  )))

  testthat::expect_s3_class(res, "bifrost_search")

  txt <- paste(testthat::capture_output(print(res)), collapse = "\n")

  testthat::expect_true(grepl("Bifrost Search Result", txt, fixed = TRUE))
  testthat::expect_true(grepl("IC (GIC)", txt, fixed = TRUE))
  testthat::expect_true(grepl("Search", txt, fixed = TRUE))
  testthat::expect_true(grepl("mvgls", txt, fixed = TRUE))
  testthat::expect_true(grepl("Shift Nodes", txt, fixed = TRUE))
  testthat::expect_true(grepl("Warnings", txt, fixed = TRUE))

  testthat::expect_false(grepl("IC History (Best IC by Iteration)", txt, fixed = TRUE))
  testthat::expect_false(grepl("Weights \\(Support\\)", txt))
})

# ---- Test XX: print method covers history plot + penalty/target + weights ----
test_that("print.bifrost_search prints history plot, penalty/target, and weights when present", {
  testthat::skip_if_not_installed("ape")
  testthat::skip_if_not_installed("phytools")
  testthat::skip_if_not_installed("txtplot")

  old_width <- getOption("width")
  old_tw <- getOption("bifrost.txtplot.width")
  old_th <- getOption("bifrost.txtplot.height")

  options(width = 80, bifrost.txtplot.width = 30L, bifrost.txtplot.height = 15L)
  on.exit({
    options(width = old_width)
    if (is.null(old_tw)) options(bifrost.txtplot.width = NULL) else options(bifrost.txtplot.width = old_tw)
    if (is.null(old_th)) options(bifrost.txtplot.height = NULL) else options(bifrost.txtplot.height = old_th)
  }, add = TRUE)

  tr <- ape::rtree(10)
  tr <- phytools::paintBranches(tr, edge = unique(tr$edge[, 2]), state = "0", anc.state = "0")
  class(tr) <- c("simmap", setdiff(class(tr), "simmap"))

  model_stub <- list(
    call = list(model = "BMM", method = "LL"),
    Y = matrix(0, nrow = ape::Ntip(tr), ncol = 3),
    formula = stats::as.formula("trait_data ~ 1")
  )

  ic_mat <- cbind(
    c(-100, -105, -103, -110, -115),
    c(1,    1,    0,    1,    1)
  )

  obj <- list(
    user_input = list(
      formula = "trait_data ~ 1",
      min_descendant_tips = 3,
      num_cores = 2,
      shift_acceptance_threshold = 20,
      plot = FALSE,
      verbose = TRUE,
      store_model_fit_history = TRUE,
      method = "LL",
      error = TRUE,
      penalty = "ridge",
      target = "CV"
    ),
    tree_no_uncertainty_untransformed = tr,
    model_no_uncertainty = model_stub,
    shift_nodes_no_uncertainty = c(12L, 15L),
    IC_used = "GIC",
    baseline_ic = -100,
    optimal_ic = -115,
    num_candidates = 5L,
    model_fit_history = list(ic_acceptance_matrix = ic_mat),
    ic_weights = data.frame(
      node = c(12L, 15L),
      ic_weight_withshift = c(0.9, 0.1)
    ),
    warnings = character(0)
  )
  class(obj) <- c("bifrost_search", "list")

  txt <- paste(testthat::capture_output(print(obj)), collapse = "\n")

  testthat::expect_true(grepl("IC History (Best IC by Iteration)", txt, fixed = TRUE))
  testthat::expect_true(grepl("\\*", txt))  # txtplot points

  testthat::expect_true(grepl("Penalty", txt, fixed = TRUE))
  testthat::expect_true(grepl("Target", txt, fixed = TRUE))

  testthat::expect_true(grepl("GIC Weights (Support)", txt, fixed = TRUE))
  testthat::expect_true(grepl("\\b12\\b", txt))
  testthat::expect_true(grepl("\\b15\\b", txt))
})

# ---- Test XX: print method edge-case branch coverage -------------------------
test_that("print.bifrost_search covers edge-case branches (NA types, fallback fields, schema issues)", {
  testthat::skip_if_not_installed("ape")
  testthat::skip_if_not_installed("phytools")
  testthat::skip_if_not_installed("txtplot")

  old_width <- getOption("width")
  old_tw <- getOption("bifrost.txtplot.width")
  old_th <- getOption("bifrost.txtplot.height")
  options(width = 80, bifrost.txtplot.width = 40L, bifrost.txtplot.height = 12L)
  on.exit({
    options(width = old_width)
    if (is.null(old_tw)) options(bifrost.txtplot.width = NULL) else options(bifrost.txtplot.width = old_tw)
    if (is.null(old_th)) options(bifrost.txtplot.height = NULL) else options(bifrost.txtplot.height = old_th)
  }, add = TRUE)

  # Case A
  obj_a <- list(
    user_input = list(
      min_descendant_tips = NA_integer_,
      num_cores = NA_integer_,
      shift_acceptance_threshold = NA_real_,
      plot = "maybe",
      verbose = character(0),
      store_model_fit_history = FALSE
    ),
    IC_used = "GIC",
    baseline_ic = "foo",
    optimal_ic = NA_real_,
    shift_nodes_no_uncertainty = NULL,
    model_no_uncertainty = NULL,
    ic_weights = data.frame(),
    warnings = character(0)
  )
  class(obj_a) <- c("bifrost_search", "list")
  txt_a <- paste(testthat::capture_output(print(obj_a)), collapse = "\n")
  testthat::expect_true(grepl("Requested, but no shifts detected", txt_a, fixed = TRUE))

  # Case B
  tr <- ape::rtree(10)
  tr0 <- phytools::paintBranches(tr, edge = unique(tr$edge[, 2]), state = "0", anc.state = "0")
  nd <- ape::Ntip(tr0) + 2L
  tr2 <- phytools::paintSubTree(tr0, node = nd, state = "1", anc.state = "0", stem = FALSE)

  long_sym <- as.name(paste(rep("x", 120), collapse = ""))

  model_stub <- list(
    call = list(method = "H&L"),
    formula = stats::as.formula("trait_data ~ 1"),
    residuals = rnorm(ape::Ntip(tr2))
  )

  obj_b <- list(
    user_input = list(
      store_model_fit_history = FALSE,
      error = 1.23,
      penalty = long_sym,
      target = as.name("CV")
    ),
    tree_no_uncertainty_transformed = tr2,
    model_no_uncertainty = model_stub,
    shift_nodes_no_uncertainty = c(12L, 15L),
    IC_used = "GIC",
    baseline_ic = -100,
    optimal_ic = -115,
    num_candidates = 5L,
    ic_weights = data.frame(foo = 1),
    warnings = character(0)
  )
  class(obj_b) <- c("bifrost_search", "list")
  txt_b <- paste(testthat::capture_output(print(obj_b)), collapse = "\n")
  testthat::expect_true(grepl("expected columns missing", txt_b, fixed = TRUE))

  # Case C (logical accept + ylab fallback)
  ic_mat_logical <- matrix(
    c(TRUE,  TRUE,
      TRUE,  FALSE,
      FALSE, TRUE,
      TRUE,  TRUE),
    ncol = 2, byrow = TRUE
  )
  obj_c <- list(
    user_input = list(store_model_fit_history = TRUE),
    IC_used = "(best)",
    baseline_ic = 1,
    optimal_ic  = 0,
    shift_nodes_no_uncertainty = integer(0),
    model_fit_history = list(ic_acceptance_matrix = ic_mat_logical),
    warnings = character(0)
  )
  class(obj_c) <- c("bifrost_search", "list")
  txt_c <- paste(testthat::capture_output(print(obj_c)), collapse = "\n")
  testthat::expect_true(grepl("IC History (Best IC by Iteration)", txt_c, fixed = TRUE))
})

# ---- Test XX: cover requireNamespace(FALSE) branches via isolated .libPaths ---
test_that("print.bifrost_search handles missing ape/phytools gracefully", {
  old_lib <- .libPaths()
  empty_lib <- tempfile("bifrost-empty-lib-")
  dir.create(empty_lib)
  .libPaths(empty_lib)
  on.exit(.libPaths(old_lib), add = TRUE)

  tree_stub <- structure(list(), class = "phylo")

  obj <- list(
    user_input = list(store_model_fit_history = FALSE),
    tree_no_uncertainty_untransformed = tree_stub,
    model_no_uncertainty = NULL,
    shift_nodes_no_uncertainty = integer(0),
    IC_used = "GIC",
    baseline_ic = -1,
    optimal_ic = -1,
    num_candidates = 0L,
    warnings = character(0)
  )
  class(obj) <- c("bifrost_search", "list")

  txt <- paste(testthat::capture_output(print(obj)), collapse = "\n")
  testthat::expect_true(grepl("Bifrost Search Result", txt, fixed = TRUE))
})

# ---- Test XX: penalty/target line is a no-op when both absent -----------------
test_that("print.bifrost_search does not print Penalty/Target when absent", {
  testthat::skip_if_not_installed("ape")
  testthat::skip_if_not_installed("phytools")

  tr <- ape::rtree(10)
  tr <- phytools::paintBranches(tr, edge = unique(tr$edge[, 2]), state = "0", anc.state = "0")
  class(tr) <- c("simmap", setdiff(class(tr), "simmap"))

  model_stub <- list(
    call = list(model = "BMM", method = "LL"),
    Y = matrix(0, nrow = ape::Ntip(tr), ncol = 2),
    formula = stats::as.formula("trait_data ~ 1")
  )

  obj <- list(
    user_input = list(
      formula = "trait_data ~ 1",
      method = "LL",
      store_model_fit_history = FALSE
      # no penalty, no target
    ),
    tree_no_uncertainty_untransformed = tr,
    model_no_uncertainty = model_stub,
    shift_nodes_no_uncertainty = integer(0),
    IC_used = "GIC",
    baseline_ic = -1,
    optimal_ic = -1,
    num_candidates = 0L,
    warnings = character(0)
  )
  class(obj) <- c("bifrost_search", "list")

  txt <- paste(testthat::capture_output(print(obj)), collapse = "\n")
  testthat::expect_false(grepl("Penalty:", txt, fixed = TRUE))
  testthat::expect_false(grepl("Target:", txt, fixed = TRUE))
})

