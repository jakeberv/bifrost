skip_if_missing_deps <- function() {
  testthat::skip_if_not_installed("ape")
  testthat::skip_if_not_installed("phytools")
  testthat::skip_if_not_installed("mvMORPH")
  testthat::skip_if_not_installed("future")
}

# ---- Test XX (NEW): print method core output on a real no-shifts run ----------
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

  out <- testthat::capture_output(print(res))
  txt <- paste(out, collapse = "\n")

  testthat::expect_true(grepl("Bifrost Search Result", txt, fixed = TRUE))
  testthat::expect_true(grepl("IC (GIC)", txt, fixed = TRUE))
  testthat::expect_true(grepl("Search", txt, fixed = TRUE))
  testthat::expect_true(grepl("mvgls", txt, fixed = TRUE))
  testthat::expect_true(grepl("Shift Nodes", txt, fixed = TRUE))
  testthat::expect_true(grepl("Warnings", txt, fixed = TRUE))

  # No history plot (history disabled)
  testthat::expect_false(grepl("IC History (Best IC by Iteration)", txt, fixed = TRUE))

  # No weights section unless ic_weights is present
  testthat::expect_false(grepl("Weights \\(Support\\)", txt))
})

# ---- Test XX (NEW): print method covers history plot + penalty/target + weights ---
test_that("print.bifrost_search prints history plot, penalty/target, and weights when present", {
  # Only the deps actually needed for this synthetic-object test
  testthat::skip_if_not_installed("ape")
  testthat::skip_if_not_installed("phytools")
  testthat::skip_if_not_installed("txtplot")

  # Replace testthat::local_options() with base options() + on.exit()
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

  out <- testthat::capture_output(print(obj))
  txt <- paste(out, collapse = "\n")

  testthat::expect_true(grepl("IC History (Best IC by Iteration)", txt, fixed = TRUE))
  testthat::expect_true(grepl("\\*", txt))  # txtplot uses '*' for points

  testthat::expect_true(grepl("Penalty", txt, fixed = TRUE))
  testthat::expect_true(grepl("Target", txt, fixed = TRUE))

  testthat::expect_true(grepl("GIC Weights (Support)", txt, fixed = TRUE))
  testthat::expect_true(grepl("\\b12\\b", txt))
  testthat::expect_true(grepl("\\b15\\b", txt))

  testthat::expect_true(grepl("Shift Nodes", txt, fixed = TRUE))
})

# ---- Test XX (NEW): print method edge-case branch coverage -------------------
test_that("print.bifrost_search covers edge-case branches (NA types, fallback fields, schema issues)", {
  testthat::skip_if_not_installed("ape")
  testthat::skip_if_not_installed("phytools")
  testthat::skip_if_not_installed("txtplot")

  # Keep output deterministic-ish for wrapping/plots
  old_width <- getOption("width")
  old_tw <- getOption("bifrost.txtplot.width")
  old_th <- getOption("bifrost.txtplot.height")
  options(width = 80, bifrost.txtplot.width = 40L, bifrost.txtplot.height = 12L)
  on.exit({
    options(width = old_width)
    if (is.null(old_tw)) options(bifrost.txtplot.width = NULL) else options(bifrost.txtplot.width = old_tw)
    if (is.null(old_th)) options(bifrost.txtplot.height = NULL) else options(bifrost.txtplot.height = old_th)
  }, add = TRUE)

  # --- Case A: missing tree/model, NA-ish user_input types, empty weights df ---
  obj_a <- list(
    user_input = list(
      min_descendant_tips = NA_integer_,
      num_cores = NA_integer_,
      shift_acceptance_threshold = NA_real_,
      plot = "maybe",                 # character -> .fmt_lgl as.character branch
      verbose = character(0),         # length==0 -> .fmt_lgl early return
      store_model_fit_history = FALSE
    ),
    IC_used = "GIC",
    baseline_ic = "foo",              # non-numeric -> .fmt_num non-numeric branch
    optimal_ic = NA_real_,            # numeric NA -> .fmt_num(all NA) branch for delta
    shift_nodes_no_uncertainty = NULL,
    model_no_uncertainty = NULL,
    ic_weights = data.frame(),        # triggers nrow(w)==0 message branch
    warnings = character(0)
  )
  class(obj_a) <- c("bifrost_search", "list")

  out_a <- testthat::capture_output(print(obj_a))
  txt_a <- paste(out_a, collapse = "\n")
  testthat::expect_true(grepl("Bifrost Search Result", txt_a, fixed = TRUE))
  testthat::expect_true(grepl("Requested, but no shifts detected", txt_a, fixed = TRUE))

  # --- Case B: transformed-tree selection + method/formula fallback + model_code fallback + bad weights schema ---
  tr <- ape::rtree(10)
  tr0 <- phytools::paintBranches(tr, edge = unique(tr$edge[, 2]), state = "0", anc.state = "0")
  # create a second regime so regimes > 1
  nd <- ape::Ntip(tr0) + 2L
  tr2 <- phytools::paintSubTree(tr0, node = nd, state = "1", anc.state = "0", stem = FALSE)

  # long symbol forces truncation branch in .fmt_val
  long_sym <- as.name(paste(rep("x", 120), collapse = ""))

  model_stub <- list(
    call = list(method = "H&L"),                   # method fallback path
    formula = stats::as.formula("trait_data ~ 1"), # formula fallback path
    residuals = rnorm(ape::Ntip(tr2))              # numeric residuals -> trait count returns 1
    # intentionally no call$model so model_code fallback uses regimes (>1)
  )

  obj_b <- list(
    user_input = list(
      store_model_fit_history = FALSE,
      # method/formula omitted so fallback from model_stub triggers
      error = 1.23,            # numeric -> covers .fmt_val numeric branch
      penalty = long_sym,      # long deparse -> truncation branch
      target = as.name("CV")   # short deparse -> "else s" (non-truncated) branch
    ),
    tree_no_uncertainty_transformed = tr2,
    model_no_uncertainty = model_stub,
    shift_nodes_no_uncertainty = c(12L, 15L),
    IC_used = "GIC",
    baseline_ic = -100,
    optimal_ic = -115,
    num_candidates = 5L,
    ic_weights = data.frame(foo = 1),              # wrong schema -> "expected columns missing" branch
    warnings = character(0)
  )
  class(obj_b) <- c("bifrost_search", "list")

  out_b <- testthat::capture_output(print(obj_b))
  txt_b <- paste(out_b, collapse = "\n")
  testthat::expect_true(grepl("Method:", txt_b, fixed = TRUE))      # fallback printed
  testthat::expect_true(grepl("Formula:", txt_b, fixed = TRUE))     # fallback printed
  testthat::expect_true(grepl("Penalty", txt_b, fixed = TRUE))
  testthat::expect_true(grepl("Target", txt_b, fixed = TRUE))
  testthat::expect_true(grepl("expected columns missing", txt_b, fixed = TRUE))

  # --- Case C: history plot with TRUE logical accept column + empty ylab fallback ---
  # Create a *logical* matrix so acc_raw is logical and hits as.integer(acc_raw)
  ic_mat_logical <- matrix(
    c(TRUE,  TRUE,
      TRUE,  FALSE,
      FALSE, TRUE,
      TRUE,  TRUE),
    ncol = 2, byrow = TRUE
  )

  obj_c <- list(
    user_input = list(store_model_fit_history = TRUE),
    IC_used = "(best)",  # gsub removes "(best)" -> empty -> ylab fallback to "IC"
    baseline_ic = 1,     # finite baseline
    optimal_ic  = 0,     # doesn't matter for plot
    shift_nodes_no_uncertainty = integer(0),
    model_fit_history = list(ic_acceptance_matrix = ic_mat_logical),
    warnings = character(0)
  )
  class(obj_c) <- c("bifrost_search", "list")

  out_c <- testthat::capture_output(print(obj_c))
  txt_c <- paste(out_c, collapse = "\n")
  testthat::expect_true(grepl("IC History (Best IC by Iteration)", txt_c, fixed = TRUE))
  testthat::expect_true(grepl("\\*", txt_c))  # plot points
})

# ---- Test XX (NEW): cover requireNamespace(FALSE) branches via isolated .libPaths ---
test_that("print.bifrost_search handles missing ape/phytools gracefully", {
  # We intentionally *hide* installed packages by changing .libPaths(),
  # so don't skip based on ape/phytools installation here.

  old_lib <- .libPaths()
  empty_lib <- tempfile("bifrost-empty-lib-")
  dir.create(empty_lib)
  .libPaths(empty_lib)
  on.exit(.libPaths(old_lib), add = TRUE)

  # Non-NULL tree object to ensure .safe_tip_count/.safe_regime_count reach requireNamespace() checks
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

  out <- testthat::capture_output(print(obj))
  txt <- paste(out, collapse = "\n")

  testthat::expect_true(grepl("Bifrost Search Result", txt, fixed = TRUE))
})

# ---- Test XX (NEW): cover requireNamespace(FALSE) guards for ape/phytools -----
test_that("print.bifrost_search covers ape/phytools requireNamespace guards", {
  # We intentionally hide installed packages by changing .libPaths()
  old_lib <- .libPaths()
  empty_lib <- tempfile("bifrost-empty-lib-")
  dir.create(empty_lib)
  .libPaths(empty_lib)
  on.exit(.libPaths(old_lib), add = TRUE)

  # Non-NULL tree so .safe_* functions reach requireNamespace() checks
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

  out <- testthat::capture_output(print(obj))
  txt <- paste(out, collapse = "\n")
  testthat::expect_true(grepl("Bifrost Search Result", txt, fixed = TRUE))
})

# ---- Test XX (NEW): print handles missing ape/phytools (covers requireNamespace guards) ----
test_that("print.bifrost_search covers ape/phytools requireNamespace guards", {
  # Hide installed packages from requireNamespace by pointing libPaths at an empty library
  old_lib <- .libPaths()
  empty_lib <- tempfile("bifrost-empty-lib-")
  dir.create(empty_lib)
  .libPaths(empty_lib)
  on.exit(.libPaths(old_lib), add = TRUE)

  # Provide a non-NULL tree so safe_tip_count/safe_regime_count reach requireNamespace()
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

  out <- testthat::capture_output(print(obj))
  txt <- paste(out, collapse = "\n")
  testthat::expect_true(grepl("Bifrost Search Result", txt, fixed = TRUE))
})
