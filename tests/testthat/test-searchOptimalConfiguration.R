# tests/testthat/test-searchOptimalConfiguration.R

testthat::skip_on_cran()

# ---- dependency guard --------------------------------------------------------
skip_if_missing_deps <- function() {
  testthat::skip_if_not_installed("ape")
  testthat::skip_if_not_installed("phytools")
  testthat::skip_if_not_installed("mvMORPH")
  testthat::skip_if_not_installed("future")
}

# ---- locate and load fixture -------------------------------------------------
load_simdata_fixture <- function() {
  # Expect the file at tests/testthat/fixtures/simdata.RDS
  rds_path <- testthat::test_path("fixtures", "simdata.RDS")
  if (!file.exists(rds_path)) {
    testthat::skip(paste("Fixture not found:", rds_path))
  }
  obj <- readRDS(rds_path)

  # Some fixtures are a list-of-lists; unwrap first element if needed
  bundle <- if (is.list(obj) && !all(c("paintedTree", "simulatedData") %in% names(obj))) obj[[1]] else obj

  # Basic sanity
  if (!all(c("paintedTree", "simulatedData") %in% names(bundle))) {
    testthat::skip("Fixture does not contain expected names: paintedTree, simulatedData")
  }
  bundle
}

# ---- helper to build 100-tip baseline + aligned data -------------------------
build_baseline_and_data <- function(simdata) {
  # Coerce paintedTree to phylo if needed
  tree0 <-
    if (inherits(simdata$paintedTree, "phylo")) simdata$paintedTree
  else ape::as.phylo(simdata$paintedTree)

  # Subsample 100 tips (for speed and determinism)
  set.seed(123)
  keep <- sample(tree0$tip.label, size = 100)
  tree0 <- ape::drop.tip(tree0, setdiff(tree0$tip.label, keep))

  # Paint a global baseline state "0" from the root
  # root <- ape::Ntip(tree0) + 1L
  # baseline <- phytools::paintSubTree(
  #   tree      = tree0,
  #   node      = root,
  #   state     = "0",
  #   anc.state = "0",
  #   stem      = FALSE
  # )
  baseline <- as.phylo(tree0)

  # Align trait data order to tree tips
  X <- simdata$simulatedData
  testthat::expect_true(is.matrix(X) || is.data.frame(X))
  testthat::expect_true(all(baseline$tip.label %in% rownames(X)))
  X <- X[baseline$tip.label, , drop = FALSE]

  list(tree = baseline, X = X)
}

# ---- tiny helpers for tolerant assertions ------------------------------------
expect_phylo_or_null <- function(x) {
  testthat::expect_true(is.null(x) || inherits(x, "phylo"))
}
expect_numeric_scalar <- function(x) {
  testthat::expect_true(is.numeric(x) && length(x) == 1L && is.finite(x))
}

# ---- Test 1: original end-to-end with parallel-capable path (num_cores = 1) --
test_that("searchOptimalConfiguration runs end-to-end on simulated data (GIC)", {
  skip_if_missing_deps()
  simdata <- load_simdata_fixture()
  built <- build_baseline_and_data(simdata)
  baseline <- built$tree
  X <- built$X

  set.seed(123)
  res <- searchOptimalConfiguration(
    baseline_tree              = baseline,
    trait_data                 = X,
    formula                    = "trait_data ~ 1",
    min_descendant_tips        = 10,
    num_cores                  = 1,
    shift_acceptance_threshold = 5,
    plot                       = FALSE,
    IC                         = "GIC",
    store_model_fit_history    = FALSE,
    method                     = "LL",
    uncertaintyweights = TRUE
  )

  # Core structure checks (present names)
  testthat::expect_type(res, "list")
  testthat::expect_true(all(c(
    "tree_no_uncertainty_transformed",
    "tree_no_uncertainty_untransformed",
    "model_no_uncertainty",
    "shift_nodes_no_uncertainty",
    "optimal_ic",
    "baseline_ic",
    "IC_used",
    "num_candidates"
  ) %in% names(res)))

  # Types/values
  expect_numeric_scalar(res$baseline_ic)
  expect_numeric_scalar(res$optimal_ic)
  testthat::expect_true(res$IC_used %in% c("GIC", "BIC"))
  # These may be NULL if no shifts are accepted; otherwise phylo
  expect_phylo_or_null(res$tree_no_uncertainty_untransformed)
  expect_phylo_or_null(res$tree_no_uncertainty_transformed)
  testthat::expect_true(is.list(res$VCVs) || is.null(res$VCVs))

  # If shifts were accepted, best should improve or match baseline
  if (length(res$shift_nodes_no_uncertainty) > 0) {
    testthat::expect_lte(res$optimal_ic, res$baseline_ic)
  } else {
    # No shifts: optimal equals baseline (within tiny tolerance)
    testthat::expect_equal(res$optimal_ic, res$baseline_ic, tolerance = 1e-8)
  }
})

# ---- Test 1a: original end-to-end with parallel-capable path (num_cores = 1) --
test_that("searchOptimalConfiguration runs end-to-end on simulated data (BIC)", {
  skip_if_missing_deps()
  simdata <- load_simdata_fixture()
  built <- build_baseline_and_data(simdata)
  baseline <- built$tree
  X <- built$X

  set.seed(123)
  res <- searchOptimalConfiguration(
    baseline_tree              = baseline,
    trait_data                 = X,
    formula                    = "trait_data ~ 1",
    min_descendant_tips        = 10,
    num_cores                  = 1,
    shift_acceptance_threshold = 5,
    plot                       = FALSE,
    IC                         = "BIC",
    store_model_fit_history    = TRUE,
    method                     = "LL",
    uncertaintyweights_par = TRUE
  )

  # Core structure checks (present names)
  testthat::expect_type(res, "list")
  testthat::expect_true(all(c(
    "tree_no_uncertainty_transformed",
    "tree_no_uncertainty_untransformed",
    "model_no_uncertainty",
    "shift_nodes_no_uncertainty",
    "optimal_ic",
    "baseline_ic",
    "IC_used",
    "num_candidates"
  ) %in% names(res)))

  # Types/values
  expect_numeric_scalar(res$baseline_ic)
  expect_numeric_scalar(res$optimal_ic)
  testthat::expect_true(res$IC_used %in% c("GIC", "BIC"))
  # These may be NULL if no shifts are accepted; otherwise phylo
  expect_phylo_or_null(res$tree_no_uncertainty_untransformed)
  expect_phylo_or_null(res$tree_no_uncertainty_transformed)
  testthat::expect_true(is.list(res$VCVs) || is.null(res$VCVs))

  # If shifts were accepted, best should improve or match baseline
  if (length(res$shift_nodes_no_uncertainty) > 0) {
    testthat::expect_lte(res$optimal_ic, res$baseline_ic)
  } else {
    # No shifts: optimal equals baseline (within tiny tolerance)
    testthat::expect_equal(res$optimal_ic, res$baseline_ic, tolerance = 1e-8)
  }
})

# ---- Test 1b: IC-weights invariants (non-brittle consistency check) ----------
test_that("ic_weights are internally consistent when present", {
  skip_if_missing_deps()

  set.seed(123)
  tr <- ape::rtree(40)
  X <- matrix(rnorm(40 * 2), ncol = 2)
  rownames(X) <- tr$tip.label

  res <- searchOptimalConfiguration(
    baseline_tree              = tr,
    trait_data                 = X,
    formula                    = "trait_data ~ 1",
    min_descendant_tips        = 5,
    num_cores                  = 1,
    shift_acceptance_threshold = -Inf,   # encourage accepting shifts
    plot                       = FALSE,
    store_model_fit_history    = FALSE,
    method                     = "LL",
    verbose                    = FALSE,
    IC                         = "GIC",
    uncertaintyweights_par     = TRUE
  )

  testthat::expect_true("ic_weights" %in% names(res))
  testthat::expect_true(is.data.frame(res$ic_weights))

  if (nrow(res$ic_weights) > 0) {
    # delta_ic should equal ic_with_shift - ic_without_shift
    testthat::expect_equal(
      res$ic_weights$delta_ic,
      res$ic_weights$ic_with_shift - res$ic_weights$ic_without_shift,
      tolerance = 1e-10
    )

    # evidence_ratio should equal weight_with / weight_without
    testthat::expect_equal(
      res$ic_weights$evidence_ratio,
      res$ic_weights$ic_weight_withshift / res$ic_weights$ic_weight_withoutshift,
      tolerance = 1e-10
    )
  }
})

# ---- Test 2: explicitly non-parallel (force sequential plan) -----------------
test_that("searchOptimalConfiguration also runs in purely sequential mode", {
  skip_if_missing_deps()
  simdata <- load_simdata_fixture()
  built <- build_baseline_and_data(simdata)
  baseline <- built$tree
  X <- built$X

  # Force a sequential future plan for the duration of this test
  old_plan <- NULL
  if (requireNamespace("future", quietly = TRUE)) {
    old_plan <- future::plan()
    future::plan(future::sequential)
  }
  on.exit({
    if (!is.null(old_plan)) future::plan(old_plan)
  }, add = TRUE)

  set.seed(456)
  res <- searchOptimalConfiguration(
    baseline_tree              = baseline,
    trait_data                 = X,
    formula                    = "trait_data ~ 1",
    min_descendant_tips        = 10,
    num_cores                  = 1,  # hint to not parallelize
    shift_acceptance_threshold = 5,
    plot                       = FALSE,
    IC                         = "GIC",
    store_model_fit_history    = FALSE,
    method                     = "LL"
  )

  # Minimal sanity checks
  testthat::expect_type(res, "list")
  expect_numeric_scalar(res$baseline_ic)
  expect_numeric_scalar(res$optimal_ic)
  expect_phylo_or_null(res$tree_no_uncertainty_untransformed)
})

# ---- Test 3: forced no-shift scenario (very high acceptance threshold) -------
test_that("searchOptimalConfiguration returns sensible output when no shifts are accepted", {
  skip_if_missing_deps()
  simdata <- load_simdata_fixture()
  built <- build_baseline_and_data(simdata)
  baseline <- built$tree
  X <- built$X

  set.seed(789)
  res <- searchOptimalConfiguration(
    baseline_tree              = baseline,
    trait_data                 = X,
    formula                    = "trait_data ~ 1",
    min_descendant_tips        = 10,
    num_cores                  = 1,
    shift_acceptance_threshold = 1e9,  # effectively forbids accepting any shift
    plot                       = FALSE,
    IC                         = "GIC",
    store_model_fit_history    = TRUE,
    method                     = "LL",
    uncertaintyweights = TRUE
  )

  # No shifts detected
  testthat::expect_equal(length(res$shift_nodes_no_uncertainty), 0L)

  # Optimal IC equals baseline IC
  testthat::expect_equal(res$optimal_ic, res$baseline_ic, tolerance = 1e-8)

  # Trees for "no-uncertainty" path may be NULL (since no shift model exists); allow NULL or phylo
  expect_phylo_or_null(res$tree_no_uncertainty_untransformed)
  expect_phylo_or_null(res$tree_no_uncertainty_transformed)

  # Model may be NULL in this path
  testthat::expect_true(is.null(res$model_no_uncertainty) || is.list(res$model_no_uncertainty))

  # History exists and contains only rejected entries — be robust to type (logical/num/char/factor)
  if (!is.null(res$model_fit_history) && is.list(res$model_fit_history)) {
    hist <- res$model_fit_history
    if (!is.null(hist$ic_acceptance_matrix)) {
      acc <- hist$ic_acceptance_matrix

      testthat::expect_true(is.matrix(acc) || is.data.frame(acc))

      if (NCOL(acc) >= 2) {
        vals <- acc[, 2]

        # flatten to a simple vector
        vals <- if (is.data.frame(vals)) unlist(vals, use.names = FALSE) else vals
        vals <- if (is.list(vals)) unlist(vals, use.names = FALSE) else vals
        if (is.factor(vals)) vals <- as.character(vals)

        # Coerce robustly and ensure no TRUE
        is_true <- rep(FALSE, length(vals))
        if (is.logical(vals)) {
          is_true <- vals
        } else if (is.numeric(vals)) {
          is_true <- vals != 0
        } else if (is.character(vals)) {
          v <- trimws(tolower(vals))
          is_true <- v %in% c("true", "t", "1")
        }
        testthat::expect_false(any(is_true, na.rm = TRUE))
      }
    }
  }
})

# ---- Test 4 (NEW): exercise acceptance + history (+plot) ---------
test_that("searchOptimalConfiguration records accepted steps with history (and covers plot/postorder)", {
  skip_if_missing_deps()
  simdata <- load_simdata_fixture()
  built <- build_baseline_and_data(simdata)
  baseline <- built$tree
  X <- built$X

  # Send plotting to a null device to cover the plotting branch safely in CI
  grDevices::pdf(NULL)
  on.exit({
    try(grDevices::dev.off(), silent = TRUE)
  }, add = TRUE)

  set.seed(10101)
  res <- searchOptimalConfiguration(
    baseline_tree              = baseline,
    trait_data                 = X,
    formula                    = "trait_data ~ 1",
    min_descendant_tips        = 10,      # broader candidate set
    num_cores                  = 1,
    shift_acceptance_threshold = -Inf,   # force acceptance of the first candidate evaluated
    plot                       = TRUE,   # hit plotSimmap/nodelabels branches
    #postorder_traversal        = TRUE,   # hit postorder switch
    IC                         = "GIC",
    store_model_fit_history    = TRUE,   # ensure history writer runs
    method                     = "LL"
  )

  # We expect at least one shift to be recorded/accepted under -Inf threshold
  testthat::expect_type(res, "list")
  testthat::expect_true(length(res$shift_nodes_no_uncertainty) >= 1L)

  # History object present and shows at least one accepted evaluation
  testthat::expect_true(!is.null(res$model_fit_history))
  hist <- res$model_fit_history
  # Flexible handling: some implementations capture acceptance in either a list of fits or a matrix
  any_accept <- FALSE
  if (is.list(hist$fits)) {
    flags <- vapply(hist$fits, function(e) isTRUE(e$accepted), logical(1))
    any_accept <- any(flags, na.rm = TRUE)
  }
  if (!any_accept && !is.null(hist$ic_acceptance_matrix)) {
    acc <- hist$ic_acceptance_matrix
    if (NCOL(acc) >= 2) {
      vals <- acc[, 2]
      vals <- if (is.data.frame(vals)) unlist(vals, use.names = FALSE) else vals
      vals <- if (is.list(vals)) unlist(vals, use.names = FALSE) else vals
      if (is.factor(vals)) vals <- as.character(vals)
      if (is.logical(vals)) {
        any_accept <- any(vals, na.rm = TRUE)
      } else if (is.numeric(vals)) {
        any_accept <- any(vals != 0, na.rm = TRUE)
      } else if (is.character(vals)) {
        v <- trimws(tolower(vals))
        any_accept <- any(v %in% c("true", "t", "1"))
      }
    }
  }
  testthat::expect_true(any_accept)

  # Final assembly fields present (ensures we traversed to the end)
  testthat::expect_true(all(c("user_input", "optimal_ic", "baseline_ic", "IC_used", "num_candidates") %in% names(res)))
})

# ---- Test 5 (NEW): testing verbose output --------
test_that("searchOptimalConfiguration emits progress output when verbose = TRUE", {
  skip_if_missing_deps()

  # Prevent helper-level verbosity (getOption("bifrost.verbose")) from interfering
  old_opt <- getOption("bifrost.verbose")
  options(bifrost.verbose = FALSE)
  on.exit(options(bifrost.verbose = old_opt), add = TRUE)

  set.seed(999)
  tr <- ape::rtree(25)
  X <- matrix(rnorm(25 * 2), ncol = 2)
  rownames(X) <- tr$tip.label

  # Helper: capture BOTH message() and stdout, then search across both
  capture_both <- function(expr) {
    out <- character()
    msgs <- testthat::capture_messages({
      out <<- testthat::capture_output(expr)
    })
    paste(c(out, msgs), collapse = "\n")
  }

  # Case A: plot = FALSE
  combined_a <- capture_both(
    searchOptimalConfiguration(
      baseline_tree              = tr,
      trait_data                 = X,
      formula                    = "trait_data ~ 1",
      min_descendant_tips        = 20,
      num_cores                  = 1,
      shift_acceptance_threshold = 1e9,
      plot                       = FALSE,
      store_model_fit_history    = FALSE,
      method                     = "LL",
      verbose                    = TRUE
    )
  )
  testthat::expect_true(grepl("Generating candidate shift models", combined_a))

  # Case B: plot = TRUE (use null device so plotting doesn’t pop windows)
  grDevices::pdf(NULL)
  on.exit(try(grDevices::dev.off(), silent = TRUE), add = TRUE)

  combined_b <- capture_both(
    searchOptimalConfiguration(
      baseline_tree              = tr,
      trait_data                 = X,
      formula                    = "trait_data ~ 1",
      min_descendant_tips        = 20,
      num_cores                  = 1,
      shift_acceptance_threshold = 1e9,
      plot                       = TRUE,
      store_model_fit_history    = FALSE,
      method                     = "LL",
      verbose                    = TRUE
    )
  )
  testthat::expect_true(grepl("Generating candidate shift models", combined_b))
})

# ---- Test 6 (NEW):  testing verbose output --------
test_that("searchOptimalConfiguration is quiet when verbose = FALSE", {
  skip_if_missing_deps()

  # Prevent helper-level verbosity (getOption("bifrost.verbose")) from interfering
  old_opt <- getOption("bifrost.verbose")
  options(bifrost.verbose = FALSE)
  on.exit(options(bifrost.verbose = old_opt), add = TRUE)

  set.seed(1000)
  tr <- ape::rtree(25)
  X <- matrix(rnorm(25 * 2), ncol = 2)
  rownames(X) <- tr$tip.label

  capture_both <- function(expr) {
    out <- character()
    msgs <- testthat::capture_messages({
      out <<- testthat::capture_output(expr)
    })
    list(out = paste(out, collapse = "\n"), msgs = paste(msgs, collapse = "\n"))
  }

  # plot = FALSE
  cap_a <- capture_both(
    searchOptimalConfiguration(
      baseline_tree              = tr,
      trait_data                 = X,
      formula                    = "trait_data ~ 1",
      min_descendant_tips        = 20,
      num_cores                  = 1,
      shift_acceptance_threshold = 1e9,
      plot                       = FALSE,
      store_model_fit_history    = FALSE,
      method                     = "LL",
      verbose                    = FALSE
    )
  )
  testthat::expect_equal(nchar(cap_a$msgs), 0)
  testthat::expect_equal(nchar(cap_a$out), 0)

  # plot = TRUE
  grDevices::pdf(NULL)
  on.exit(try(grDevices::dev.off(), silent = TRUE), add = TRUE)

  cap_b <- capture_both(
    searchOptimalConfiguration(
      baseline_tree              = tr,
      trait_data                 = X,
      formula                    = "trait_data ~ 1",
      min_descendant_tips        = 20,
      num_cores                  = 1,
      shift_acceptance_threshold = 1e9,
      plot                       = TRUE,
      store_model_fit_history    = FALSE,
      method                     = "LL",
      verbose                    = FALSE
    )
  )
  testthat::expect_equal(nchar(cap_b$msgs), 0)
  testthat::expect_equal(nchar(cap_b$out), 0)
})

# ---- Test 7 (NEW): invalid IC guard (covers stop("IC must be GIC or BIC")) ----
test_that("searchOptimalConfiguration errors on invalid IC", {
  skip_if_missing_deps()

  set.seed(1)
  tr <- ape::rtree(25)
  X <- matrix(rnorm(25 * 2), ncol = 2)
  rownames(X) <- tr$tip.label

  testthat::expect_error(
    searchOptimalConfiguration(
      baseline_tree              = tr,
      trait_data                 = X,
      formula                    = "trait_data ~ 1",
      min_descendant_tips        = 10,
      num_cores                  = 1,
      shift_acceptance_threshold = 1e9,
      plot                       = FALSE,
      store_model_fit_history    = FALSE,
      method                     = "LL",
      verbose                    = FALSE,
      IC                         = "AIC"
    ),
    "IC must be GIC or BIC"
  )
})

# ---- Test 8 (NEW): both IC-weight flags TRUE should error --------------------
test_that("searchOptimalConfiguration errors if both uncertaintyweights flags are TRUE", {
  skip_if_missing_deps()

  set.seed(2)
  tr <- ape::rtree(25)
  X <- matrix(rnorm(25 * 2), ncol = 2)
  rownames(X) <- tr$tip.label

  testthat::expect_error(
    searchOptimalConfiguration(
      baseline_tree              = tr,
      trait_data                 = X,
      formula                    = "trait_data ~ 1",
      min_descendant_tips        = 20,
      num_cores                  = 1,
      shift_acceptance_threshold = 1e9,
      plot                       = FALSE,
      store_model_fit_history    = FALSE,
      method                     = "LL",
      verbose                    = FALSE,
      IC                         = "GIC",
      uncertaintyweights         = TRUE,
      uncertaintyweights_par     = TRUE
    ),
    "Exactly one of uncertaintyweights or uncertaintyweights_par must be TRUE"
  )
})

# ---- Test 9 (NEW): no-shifts path for uncertaintyweights_par (ic_weights = NA) ----
test_that("searchOptimalConfiguration skips IC weights (parallel) when no shifts are detected", {
  skip_if_missing_deps()

  set.seed(3)
  tr <- ape::rtree(25)
  X <- matrix(rnorm(25 * 2), ncol = 2)
  rownames(X) <- tr$tip.label

  res <- searchOptimalConfiguration(
    baseline_tree              = tr,
    trait_data                 = X,
    formula                    = "trait_data ~ 1",
    min_descendant_tips        = 10,
    num_cores                  = 1,
    shift_acceptance_threshold = 1e9,
    plot                       = FALSE,
    store_model_fit_history    = FALSE,
    method                     = "LL",
    verbose                    = FALSE,
    IC                         = "GIC",
    uncertaintyweights_par     = TRUE
  )

  testthat::expect_true("ic_weights" %in% names(res))
  testthat::expect_true(is.data.frame(res$ic_weights))
  testthat::expect_equal(nrow(res$ic_weights), 0L)

  # Optional: check schema stays stable
  testthat::expect_true(all(c(
    "node","ic_with_shift","ic_without_shift","delta_ic",
    "ic_weight_withshift","ic_weight_withoutshift","evidence_ratio"
  ) %in% names(res$ic_weights)))
})

# ---- Test 10 (NEW): force multisession branch (RSTUDIO=1) --------------------
test_that("searchOptimalConfiguration takes multisession path when RSTUDIO=1", {
  skip_if_missing_deps()

  old_rstudio <- Sys.getenv("RSTUDIO", unset = NA_character_)
  Sys.setenv(RSTUDIO = "1")
  on.exit({
    if (is.na(old_rstudio)) Sys.unsetenv("RSTUDIO") else Sys.setenv(RSTUDIO = old_rstudio)
  }, add = TRUE)

  set.seed(4)
  tr <- ape::rtree(25)
  X <- matrix(rnorm(25 * 2), ncol = 2)
  rownames(X) <- tr$tip.label

  res <- searchOptimalConfiguration(
    baseline_tree              = tr,
    trait_data                 = X,
    formula                    = "trait_data ~ 1",
    min_descendant_tips        = 20,
    num_cores                  = 1,
    shift_acceptance_threshold = 1e9,
    plot                       = FALSE,
    store_model_fit_history    = FALSE,
    method                     = "LL",
    verbose                    = FALSE,
    IC                         = "GIC"
  )

  testthat::expect_type(res, "list")
})

# ---- Test 11 (NEW): restore_threads else-branch (pre-set env var restored) ----
test_that("searchOptimalConfiguration restores BLAS/OpenMP env vars after candidate scoring", {
  skip_if_missing_deps()

  thread_vars <- c(
    "OMP_NUM_THREADS","OPENBLAS_NUM_THREADS","MKL_NUM_THREADS",
    "VECLIB_MAXIMUM_THREADS","NUMEXPR_NUM_THREADS"
  )
  old <- Sys.getenv(thread_vars, unset = NA_character_)

  # Pre-set at least one var so restore hits the else branch
  Sys.setenv(OMP_NUM_THREADS = "3")
  on.exit({
    for (nm in names(old)) {
      val <- old[[nm]]
      if (is.na(val) || val == "") {
        Sys.unsetenv(nm)
      } else {
        do.call(Sys.setenv, setNames(list(val), nm))
      }
    }
  }, add = TRUE)

  set.seed(5)
  tr <- ape::rtree(25)
  X <- matrix(rnorm(25 * 2), ncol = 2)
  rownames(X) <- tr$tip.label

  searchOptimalConfiguration(
    baseline_tree              = tr,
    trait_data                 = X,
    formula                    = "trait_data ~ 1",
    min_descendant_tips        = 20,
    num_cores                  = 1,
    shift_acceptance_threshold = 1e9,
    plot                       = FALSE,
    store_model_fit_history    = FALSE,
    method                     = "LL",
    verbose                    = FALSE,
    IC                         = "GIC"
  )

  testthat::expect_identical(Sys.getenv("OMP_NUM_THREADS"), "3")
})

# ---- Test 12 (NEW): serial vs parallel IC-weights agree (num_cores = 1) ------
test_that("searchOptimalConfiguration returns consistent ic_weights for serial vs parallel modes", {
  skip_if_missing_deps()

  set.seed(1414)

  # Build a tree and pre-paint it (matches what searchOptimalConfiguration does internally)
  tr0 <- ape::rtree(40)
  tr  <- phytools::paintSubTree(ape::as.phylo(tr0), node = ape::Ntip(tr0) + 1L, state = 0)

  # Build trait matrix aligned to the painted tree tip order
  X <- matrix(rnorm(40 * 2), ncol = 2)
  rownames(X) <- tr$tip.label

  # Run 1: serial IC weights
  set.seed(1414)
  res_ser <- suppressWarnings(searchOptimalConfiguration(
    baseline_tree              = tr,
    trait_data                 = X,
    formula                    = "trait_data ~ 1",
    min_descendant_tips        = 5,
    num_cores                  = 1,
    shift_acceptance_threshold = -Inf,
    plot                       = FALSE,
    store_model_fit_history    = FALSE,
    method                     = "LL",
    verbose                    = FALSE,
    IC                         = "GIC",
    uncertaintyweights         = TRUE
  ))

  # Run 2: parallel IC weights (deterministic when num_cores = 1)
  set.seed(1414)
  res_par <- suppressWarnings(searchOptimalConfiguration(
    baseline_tree              = tr,
    trait_data                 = X,
    formula                    = "trait_data ~ 1",
    min_descendant_tips        = 5,
    num_cores                  = 1,
    shift_acceptance_threshold = -Inf,
    plot                       = FALSE,
    store_model_fit_history    = FALSE,
    method                     = "LL",
    verbose                    = FALSE,
    IC                         = "GIC",
    uncertaintyweights_par     = TRUE
  ))

  testthat::expect_true("ic_weights" %in% names(res_ser))
  testthat::expect_true("ic_weights" %in% names(res_par))
  testthat::expect_true(is.data.frame(res_ser$ic_weights))
  testthat::expect_true(is.data.frame(res_par$ic_weights))

  # If either is empty, both should be empty under fixed seed + num_cores=1
  if (nrow(res_ser$ic_weights) == 0L || nrow(res_par$ic_weights) == 0L) {
    testthat::expect_equal(nrow(res_ser$ic_weights), 0L)
    testthat::expect_equal(nrow(res_par$ic_weights), 0L)
    return(invisible(NULL))
  }

  cols <- c(
    "node", "ic_with_shift", "ic_without_shift", "delta_ic",
    "ic_weight_withshift", "ic_weight_withoutshift", "evidence_ratio"
  )

  testthat::expect_true(all(cols %in% names(res_ser$ic_weights)))
  testthat::expect_true(all(cols %in% names(res_par$ic_weights)))

  ser <- res_ser$ic_weights[, cols, drop = FALSE]
  par <- res_par$ic_weights[, cols, drop = FALSE]

  # order-independent compare
  ser <- ser[order(ser$node), , drop = FALSE]; rownames(ser) <- NULL
  par <- par[order(par$node), , drop = FALSE]; rownames(par) <- NULL

  testthat::expect_equal(ser$node, par$node)

  num_cols <- setdiff(cols, "node")
  for (nm in num_cols) {
    testthat::expect_equal(ser[[nm]], par[[nm]], tolerance = 1e-10)
  }
})

# ---- Test 15: CRAN-safety: do not write files to the working directory --------
test_that("searchOptimalConfiguration does not write files to the working directory", {
  skip_if_missing_deps()

  # Isolated working directory (base R only)
  wd <- tempfile("bifrost-wd-")
  dir.create(wd, recursive = TRUE)
  oldwd <- getwd()
  on.exit(setwd(oldwd), add = TRUE)
  setwd(wd)

  before <- list.files(".", recursive = TRUE, all.files = TRUE)

  set.seed(123)
  tr <- ape::rtree(30)
  X <- matrix(rnorm(30 * 2), ncol = 2)
  rownames(X) <- tr$tip.label

  res <- searchOptimalConfiguration(
    baseline_tree              = tr,
    trait_data                 = X,
    formula                    = "trait_data ~ 1",
    min_descendant_tips        = 10,
    num_cores                  = 1,
    shift_acceptance_threshold = 1e9,
    plot                       = FALSE,
    IC                         = "GIC",
    store_model_fit_history    = TRUE,  # key path we care about
    method                     = "LL",
    verbose                    = FALSE
  )

  after <- list.files(".", recursive = TRUE, all.files = TRUE)

  # Nothing new should appear in working directory
  testthat::expect_identical(after, before)

  # Optional sanity: if history is requested, it should be under tempdir(), not getwd()
  testthat::expect_true(
    !dir.exists(file.path(getwd(), "bifrost_fit_history"))
  )
})

# ---- Test 16: does not leak global bifrost.verbose option --------------------
test_that("searchOptimalConfiguration restores options(bifrost.verbose)", {
  skip_if_missing_deps()

  old_opt <- getOption("bifrost.verbose")
  options(bifrost.verbose = FALSE)
  on.exit(options(bifrost.verbose = old_opt), add = TRUE)

  set.seed(123)
  tr <- ape::rtree(25)
  X <- matrix(rnorm(25 * 2), ncol = 2)
  rownames(X) <- tr$tip.label

  res <- suppressMessages(searchOptimalConfiguration(
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
    verbose                    = TRUE   # triggers internal options() change
  ))

  testthat::expect_identical(getOption("bifrost.verbose"), FALSE)
})

# ---- Test 17: ic_weights schema is stable when requested ---------------------
test_that("ic_weights has stable schema when requested", {
  skip_if_missing_deps()

  set.seed(3)
  tr <- ape::rtree(25)
  X <- matrix(rnorm(25 * 2), ncol = 2)
  rownames(X) <- tr$tip.label

  res <- searchOptimalConfiguration(
    baseline_tree              = tr,
    trait_data                 = X,
    formula                    = "trait_data ~ 1",
    min_descendant_tips        = 10,
    num_cores                  = 1,
    shift_acceptance_threshold = 1e9, # likely no shifts
    plot                       = FALSE,
    store_model_fit_history    = FALSE,
    method                     = "LL",
    verbose                    = FALSE,
    IC                         = "GIC",
    uncertaintyweights_par     = TRUE
  )

  testthat::expect_true("ic_weights" %in% names(res))
  testthat::expect_true(is.data.frame(res$ic_weights))
  testthat::expect_true(all(c(
    "node","ic_with_shift","ic_without_shift","delta_ic",
    "ic_weight_withshift","ic_weight_withoutshift","evidence_ratio"
  ) %in% names(res$ic_weights)))
})

# ---- Test 18: model_fit_history ic_acceptance_matrix is well-formed ----------
test_that("model_fit_history ic_acceptance_matrix is well-formed", {
  skip_if_missing_deps()

  simdata <- load_simdata_fixture()
  built <- build_baseline_and_data(simdata)

  set.seed(123)
  res <- searchOptimalConfiguration(
    baseline_tree           = built$tree,
    trait_data              = built$X,
    formula                 = "trait_data ~ 1",
    min_descendant_tips     = 10,
    num_cores               = 1,
    shift_acceptance_threshold = 1e9,
    plot                    = FALSE,
    IC                      = "GIC",
    store_model_fit_history = TRUE,
    method                  = "LL",
    verbose                 = FALSE
  )

  testthat::expect_true(is.list(res$model_fit_history))
  testthat::expect_true("ic_acceptance_matrix" %in% names(res$model_fit_history))

  mat <- res$model_fit_history$ic_acceptance_matrix
  testthat::expect_true(is.matrix(mat))
  testthat::expect_equal(ncol(mat), 2L)

  # acceptance column should be coercible to logical without producing all NA
  acc <- mat[, 2]
  acc_lgl <- suppressWarnings(as.logical(acc))
  testthat::expect_false(all(is.na(acc_lgl)))
})

# ---- Test 19: no-shifts path yields a usable model_no_uncertainty -------------
test_that("no-shifts path yields a usable model_no_uncertainty", {
  skip_if_missing_deps()

  set.seed(999)
  tr <- ape::rtree(25)
  X <- matrix(rnorm(25 * 2), ncol = 2)
  rownames(X) <- tr$tip.label

  res <- searchOptimalConfiguration(
    baseline_tree              = tr,
    trait_data                 = X,
    formula                    = "trait_data ~ 1",
    min_descendant_tips        = 20,
    num_cores                  = 1,
    shift_acceptance_threshold = 1e9, # no shifts
    plot                       = FALSE,
    IC                         = "GIC",
    store_model_fit_history    = FALSE,
    method                     = "LL",
    verbose                    = FALSE
  )

  testthat::expect_equal(length(res$shift_nodes_no_uncertainty), 0L)
  # should not be NULL in a "nice" API; if it is today, this test will flag it.
  testthat::expect_true(!is.null(res$model_no_uncertainty))
})

# ---- Test 13 (NEW): warning handler path during shift evaluation --------------
test_that("searchOptimalConfiguration captures warnings from shift evaluation", {
  skip_if_missing_deps()

  ns <- asNamespace("bifrost")
  orig_fit <- get("fitMvglsAndExtractGIC.formula", envir = ns)

  testthat::local_mocked_bindings(
    fitMvglsAndExtractGIC.formula = function(formula, tree, trait_data, ...) {
      in_withCallingHandlers <- any(vapply(sys.calls(), function(cl) {
        is.call(cl) && is.name(cl[[1]]) && identical(as.character(cl[[1]]), "withCallingHandlers")
      }, logical(1)))

      if (in_withCallingHandlers) warning("forced warning from test")

      orig_fit(formula, tree, trait_data, ...)
    },
    .env = ns
  )

  set.seed(6)
  tr <- ape::rtree(25)
  X <- matrix(rnorm(25 * 2), ncol = 2)
  rownames(X) <- tr$tip.label

  res <- suppressWarnings(searchOptimalConfiguration(
    baseline_tree              = tr,
    trait_data                 = X,
    formula                    = "trait_data ~ 1",
    min_descendant_tips        = 5,
    num_cores                  = 1,
    shift_acceptance_threshold = 1e9,
    plot                       = FALSE,
    store_model_fit_history    = FALSE,
    method                     = "LL",
    verbose                    = FALSE,
    IC                         = "GIC"
  ))

  testthat::expect_true(!is.null(res$warnings))
  testthat::expect_true(any(grepl("Warning in evaluating shift at node", unlist(res$warnings))))
})

# ---- Test 14 (NEW): error handler + NA_real_ row in ic_acceptance_matrix ------
test_that("searchOptimalConfiguration records error entries in history and yields NA_real_ row", {
  skip_if_missing_deps()

  ns <- asNamespace("bifrost")

  testthat::local_mocked_bindings(
    addShiftToModel = function(tree, shift_node, shift_id) {
      list(tree = NULL, shift_id = shift_id + 1L)  # shifted_tree becomes NULL => fit errors
    },
    .env = ns
  )

  set.seed(7)
  tr <- ape::rtree(25)
  X <- matrix(rnorm(25 * 2), ncol = 2)
  rownames(X) <- tr$tip.label

  res <- suppressWarnings(searchOptimalConfiguration(
    baseline_tree              = tr,
    trait_data                 = X,
    formula                    = "trait_data ~ 1",
    min_descendant_tips        = 5,
    num_cores                  = 1,
    shift_acceptance_threshold = 1e9,
    plot                       = FALSE,
    store_model_fit_history    = TRUE,
    method                     = "LL",
    verbose                    = FALSE,
    IC                         = "GIC"
  ))

  testthat::expect_true(!is.null(res$model_fit_history$ic_acceptance_matrix))
  mat <- res$model_fit_history$ic_acceptance_matrix
  testthat::expect_true(any(is.na(mat[, 1]), na.rm = TRUE))

  testthat::expect_true(!is.null(res$warnings))
  testthat::expect_true(any(grepl("Error in evaluating shift at node", unlist(res$warnings))))
})

# ---- Test XX (NEW): cover .progress cat()/flush.console() branch (interactive only) ----
test_that("searchOptimalConfiguration uses cat() progress path in interactive RStudio plotting", {
  skip_if_missing_deps()
  testthat::skip_if_not(interactive())

  old_rstudio <- Sys.getenv("RSTUDIO", unset = NA_character_)
  Sys.setenv(RSTUDIO = "1")
  on.exit({
    if (is.na(old_rstudio)) Sys.unsetenv("RSTUDIO") else Sys.setenv(RSTUDIO = old_rstudio)
  }, add = TRUE)

  # Null device so plot=TRUE doesn't pop windows
  grDevices::pdf(NULL)
  on.exit(try(grDevices::dev.off(), silent = TRUE), add = TRUE)

  set.seed(123)
  tr <- ape::rtree(20)
  X <- matrix(rnorm(20 * 2), ncol = 2)
  rownames(X) <- tr$tip.label

  # Critical trick:
  # Set min_descendant_tips to Ntip so ONLY the root is eligible.
  # That makes candidate_trees_shifts empty => the main loop never runs,
  # so the plot code never calls getStates(shifted_tree,...).
  out <- testthat::capture_output({
    suppressWarnings(suppressMessages(searchOptimalConfiguration(
      baseline_tree              = tr,
      trait_data                 = X,
      formula                    = "trait_data ~ 1",
      min_descendant_tips        = ape::Ntip(tr),
      num_cores                  = 1,
      shift_acceptance_threshold = 1e9,
      plot                       = TRUE,
      IC                         = "GIC",
      store_model_fit_history    = FALSE,
      method                     = "LL",
      verbose                    = TRUE
    )))
  })

  txt <- paste(out, collapse = "\n")
  testthat::expect_true(grepl("Generating candidate shift models", txt))
})

# ---- Test XX (NEW): serial IC weights path executes BIC branch (mocked deterministic) ----
test_that("searchOptimalConfiguration serial ic_weights executes BIC branch", {
  skip_if_missing_deps()

  ns <- asNamespace("bifrost")

  # Deterministic decreasing BIC so shifts are accepted and weights run
  k <- 0L
  testthat::local_mocked_bindings(
    fitMvglsAndExtractBIC.formula = function(formula, tree, trait_data, ...) {
      k <<- k + 1L
      bic_val <- 1000 - 10 * k
      list(
        model = list(corrSt = list(phy = tree)),
        BIC = list(BIC = bic_val)
      )
    },
    # Make shift-removal safe & deterministic for weights loop
    removeShiftFromTree = function(tree, shift_node, stem = FALSE) tree,
    .env = ns
  )

  set.seed(999)
  tr <- ape::rtree(40)
  X <- matrix(rnorm(40 * 2), ncol = 2)
  rownames(X) <- tr$tip.label

  res <- suppressWarnings(searchOptimalConfiguration(
    baseline_tree              = tr,
    trait_data                 = X,
    formula                    = "trait_data ~ 1",
    min_descendant_tips        = 2,
    num_cores                  = 1,
    shift_acceptance_threshold = -Inf,   # accept shifts
    plot                       = FALSE,
    store_model_fit_history    = FALSE,
    verbose                    = FALSE,
    IC                         = "BIC",
    uncertaintyweights         = TRUE    # SERIAL weights path
  ))

  testthat::expect_true("ic_weights" %in% names(res))
  testthat::expect_true(is.data.frame(res$ic_weights))
  testthat::expect_true(nrow(res$ic_weights) >= 1L)
})

# ---- Test XX (NEW): warning handler path during shift evaluation (BIC) ----
test_that("searchOptimalConfiguration captures warnings from shift evaluation (BIC)", {
  skip_if_missing_deps()

  ns <- asNamespace("bifrost")
  orig_fit <- get("fitMvglsAndExtractBIC.formula", envir = ns)

  testthat::local_mocked_bindings(
    fitMvglsAndExtractBIC.formula = function(formula, tree, trait_data, ...) {
      in_withCallingHandlers <- any(vapply(sys.calls(), function(cl) {
        is.call(cl) && is.name(cl[[1]]) && identical(as.character(cl[[1]]), "withCallingHandlers")
      }, logical(1)))

      if (in_withCallingHandlers) warning("forced warning from test (BIC)")

      orig_fit(formula, tree, trait_data, ...)
    },
    .env = ns
  )

  set.seed(456)
  tr <- ape::rtree(25)
  X <- matrix(rnorm(25 * 2), ncol = 2)
  rownames(X) <- tr$tip.label

  res <- suppressWarnings(searchOptimalConfiguration(
    baseline_tree              = tr,
    trait_data                 = X,
    formula                    = "trait_data ~ 1",
    min_descendant_tips        = 5,
    num_cores                  = 1,
    shift_acceptance_threshold = 1e9,  # reject; still evaluates candidates and triggers handler
    plot                       = FALSE,
    store_model_fit_history    = FALSE,
    method                     = "LL",
    verbose                    = FALSE,
    IC                         = "BIC"
  ))

  testthat::expect_true(!is.null(res$warnings))
  testthat::expect_true(any(grepl("Warning in evaluating shift at node", unlist(res$warnings))))
})

