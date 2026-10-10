.search_thread_env <- c(
  "OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
  "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS"
)

test_that("thread limits depend on fitting workers independently of progress", {
  withr::local_envvar(c(
    OMP_NUM_THREADS = "4", OPENBLAS_NUM_THREADS = "3",
    MKL_NUM_THREADS = "2", VECLIB_MAXIMUM_THREADS = NA,
    NUMEXPR_NUM_THREADS = NA
  ))
  caller_threads <- c("4", "3", "2", NA_character_, NA_character_)
  caller_pid <- Sys.getpid()
  probe <- function(i) list(
    pid = Sys.getpid(),
    threads = unname(Sys.getenv(.search_thread_env, unset = NA_character_))
  )

  backends <- if (future::supportsMulticore()) c(TRUE, FALSE) else TRUE
  for (is_rstudio in backends) {
    for (progress in c(FALSE, TRUE)) {
      heartbeat <- if (progress) function() invisible(NULL) else NULL
      for (jobs in 1:2) {
        # One job needs one fitting worker even when two cores are requested.
        expected <- if (jobs == 1L) caller_threads else rep("1", 5L)
        for (cores in jobs:2) {
          result <- .bifrost_search_lapply(
            seq_len(jobs), probe, num_cores = cores, is_rstudio = is_rstudio,
            heartbeat = heartbeat
          )
          for (worker in result) expect_identical(worker$threads, expected)
          pids <- vapply(result, `[[`, integer(1), "pid")
          expect_length(unique(pids), jobs)
          if (progress || jobs > 1L) expect_false(any(pids == caller_pid))
          expect_identical(
            unname(Sys.getenv(.search_thread_env, unset = NA_character_)),
            caller_threads
          )
        }
      }
    }
  }
})

test_that("worker thread limits never change the calling process", {
  # Native thread counts cannot be queried portably. Catch either setter at
  # the library boundary without changing the test runner's thread settings.
  testthat::local_mocked_bindings(
    blas_set_num_threads = function(...) {
      stop("Attempted to change the caller's BLAS threads")
    },
    omp_set_num_threads = function(...) {
      stop("Attempted to change the caller's OpenMP threads")
    },
    .package = "RhpcBLASctl"
  )

  expect_invisible(.bifrost_search_limit_worker_threads(Sys.getpid()))
})

test_that("legacy Future plans are reinitialized after stage cleanup", {
  old_plan <- future::plan("list")
  on.exit(future::plan(old_plan), add = TRUE)
  future::plan(future::sequential)
  initialization_probe <- future::future(NA)
  active <- FALSE
  legacy <- function(..., workers = 1L) {
    active <<- TRUE
    initialization_probe
  }
  class(legacy) <- c("legacy", "future", "function")
  attr(legacy, "init") <- TRUE
  attr(legacy, "cleanup") <- function() active <<- FALSE
  future::plan(legacy)
  expect_true(active)
  expect_null(attr(future::plan("next"), "backend", exact = TRUE))
  expect_identical(attr(future::plan("next"), "init", exact = TRUE), "done")

  for (fail in c(FALSE, TRUE)) {
    work <- function() {
      if (fail) stop("synthetic fitting failure")
      42L
    }
    run <- function() .bifrost_search_with_future_plan(1L, TRUE, work)
    if (fail) expect_error(run(), "synthetic fitting failure")
    else expect_identical(run(), 42L)
    expect_true(active)
    expect_s3_class(future::plan("next"), "legacy")
  }
})

test_that("an existing multisession pool cannot bypass the thread policy", {
  withr::local_envvar(c(OMP_NUM_THREADS = "4", OPENBLAS_NUM_THREADS = "3"))
  old_plan <- future::plan("list")
  on.exit(future::plan(old_plan), add = TRUE)
  future::plan(future::sequential)
  connections_before <- nrow(showConnections())
  probe <- function(i) Sys.getenv("OPENBLAS_NUM_THREADS")
  for (initial_workers in c(1L, 2L)) {
    workers <- if (initial_workers == 1L) I(1L) else initial_workers
    future::plan(future::multisession, workers = workers)
    expect_identical(future::value(future::future(probe(1))), "3")
    for (progress in c(FALSE, TRUE)) {
      result <- .bifrost_search_lapply(
        1:2, probe, num_cores = 2L, is_rstudio = TRUE,
        heartbeat = if (progress) function() invisible(NULL) else NULL
      )
      expect_identical(result, list("1", "1"))
      expect_equal(future::nbrOfWorkers(), initial_workers)
      expect_identical(future::value(future::future(probe(1))), "3")
    }
    future::plan(future::sequential)
    expect_warning(gc(), NA)
    expect_equal(nrow(showConnections()), connections_before)
  }
})

test_that("disabled forking still runs fitting work outside the caller", {
  withr::local_options(future.fork.enable = FALSE)
  caller_pid <- Sys.getpid()
  for (cores in 1:2) {
    result <- .bifrost_search_lapply(
      seq_len(cores), function(i) Sys.getpid(), num_cores = cores,
      is_rstudio = FALSE, heartbeat = function() invisible(NULL)
    )
    expect_false(any(unlist(result) == caller_pid))
    expect_length(unique(unlist(result)), cores)
  }
})

test_that("worker counts are validated before being limited by job count", {
  for (jobs in 0:1) {
    for (cores in list(c(1, 2), numeric(), Inf, NA_real_, "2")) {
      expect_error(
        .bifrost_search_lapply(seq_len(jobs), identity, cores, is_rstudio = TRUE),
        "`num_cores` must be a single finite number", fixed = TRUE
      )
    }
  }
  expect_identical(
    .bifrost_search_lapply(integer(), identity, 2L, is_rstudio = TRUE), list()
  )
  withr::local_envvar(OMP_NUM_THREADS = "4")
  result <- .bifrost_search_lapply(
    1:2, function(i) Sys.getenv("OMP_NUM_THREADS"),
    num_cores = 1.5, is_rstudio = TRUE
  )
  expect_identical(result, list("4", "4"))
})

test_that("one-worker progress restores the caller plan after a failed fit", {
  withr::local_envvar(c(OMP_NUM_THREADS = "4", OPENBLAS_NUM_THREADS = NA))
  old_plan <- future::plan("list")
  on.exit(future::plan(old_plan), add = TRUE)
  future::plan(future::sequential)
  expect_error(
    .bifrost_search_lapply(
      1L, function(i) stop("synthetic serial failure"),
      num_cores = 1L, is_rstudio = TRUE,
      heartbeat = function() invisible(NULL)
    ),
    "synthetic serial failure"
  )
  expect_s3_class(future::plan("next"), "sequential")
  expect_identical(Sys.getenv("OMP_NUM_THREADS"), "4")
  expect_true(is.na(Sys.getenv("OPENBLAS_NUM_THREADS", unset = NA_character_)))
})

test_that("greedy search and serial weights preserve numerical threading", {
  withr::local_envvar(c(OMP_NUM_THREADS = "4", OPENBLAS_NUM_THREADS = "3"))
  set.seed(46)
  tree <- ape::rtree(10)
  baseline <- phytools::paintSubTree(tree, node = 11L, state = 0)
  candidate <- generatePaintedTrees(baseline, min_tips = 3)[2]
  caller_pid <- Sys.getpid()
  fit <- function(...) {
    # Observe settings at the fitting boundary, including in a real worker.
    if (!identical(Sys.getenv("OMP_NUM_THREADS"), "4") ||
        !identical(Sys.getenv("OPENBLAS_NUM_THREADS"), "3")) {
      stop("Numerical thread settings were changed during a serial fit")
    }
    # A fast fit can finish before the first heartbeat poll. Check its process
    # directly; deterministic polling coverage lives in test-search-spinner.R.
    if (!identical(Sys.getpid() != caller_pid, progress)) {
      stop("Fit ran in an unexpected process for the progress setting")
    }
    list(GIC = list(GIC = 90))
  }
  backends <- if (future::supportsMulticore()) c(TRUE, FALSE) else TRUE
  for (is_rstudio in backends) {
    for (progress in c(FALSE, TRUE)) {
      heartbeat <- if (progress) function() invisible(NULL) else NULL
      result <- .bifrost_search_forward(
        sorted_candidates = candidate, current_best_tree = baseline,
        current_best_ic = 100, shift_id = 0L, IC = "GIC",
        formula = trait_data ~ 1, trait_data = matrix(0, 10, 1),
        shift_acceptance_threshold = 5, store_model_fit_history = FALSE,
        sub_dir = NULL, plot = FALSE, verbose_log = function(...) NULL,
        heartbeat = heartbeat, is_rstudio = is_rstudio, fit = fit
      )
      expect_length(result$shift_vec, 1L)
      weights <- .bifrost_search_calculate_ic_weights(
        uncertaintyweights = TRUE, uncertaintyweights_par = FALSE,
        shift_vec = result$shift_vec,
        best_tree_no_uncertainty = result$best_tree_no_uncertainty,
        model_with_shift_no_uncertainty = result$model_with_shift_no_uncertainty,
        IC = "GIC", formula = trait_data ~ 1, trait_data = matrix(0, 10, 1),
        args_list = list(), num_cores = 1L, is_rstudio = is_rstudio,
        verbose_log = function(...) NULL, heartbeat = heartbeat, fit = fit
      )
      expect_equal(nrow(weights), 1L)
    }
  }
})
