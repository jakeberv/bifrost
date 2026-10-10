background_history_args <- function() {
  withr::local_seed(44)
  tree <- ape::rtree(12)
  baseline <- phytools::paintSubTree(tree, node = 13L, state = 0)
  candidates <- generatePaintedTrees(baseline, min_tips = 3)[-1][1:3]
  list(sorted_candidates = candidates, current_best_tree = baseline,
       current_best_ic = 100, shift_id = 0L, IC = "GIC",
       formula = trait_data ~ 1, trait_data = matrix(0, 12, 1),
       shift_acceptance_threshold = 5, store_model_fit_history = TRUE,
       plot = FALSE, verbose_log = function(...) invisible(NULL))
}

test_that("proposal history is saved in the fitting worker while the caller ticks", {
  withr::local_seed(82)
  backends <- c(TRUE, if (future::supportsMulticore()) FALSE)
  for (backend in backends) {
    # Private copies carry the instrumented writer into either worker backend.
    ns <- asNamespace("bifrost")
    env <- new.env(parent = ns)
    for (name in c(".bifrost_search_forward", ".bifrost_search_fit_proposal")) {
      fun <- get(name, ns)
      environment(fun) <- env
      env[[name]] <- fun
    }
    directory <- withr::local_tempdir()
    writing <- file.path(directory, "writing")
    finished <- file.path(directory, "finished")
    original_save <- .bifrost_search_save_history
    env$.bifrost_search_save_history <- function(...) {
      writeLines(as.character(Sys.getpid()), writing)
      Sys.sleep(0.35)
      original_save(...)
      file.create(finished)
      invisible(NULL)
    }
    args <- background_history_args()
    args$sorted_candidates <- args$sorted_candidates[1]
    args$sub_dir <- directory
    args$is_rstudio <- backend
    beats_during_save <- 0L
    args$heartbeat <- function() {
      if (file.exists(writing) && !file.exists(finished)) {
        beats_during_save <<- beats_during_save + 1L
      }
    }
    args$fit <- function(...) list(GIC = list(GIC = 90), pid = Sys.getpid())
    result <- do.call(env$.bifrost_search_forward, args)
    saved <- readRDS(file.path(directory, "iteration_2.rds"))
    expect_gt(beats_during_save, 0L)
    expect_false(identical(saved$model$pid, Sys.getpid()))
    expect_identical(as.integer(readLines(writing)), saved$model$pid)
    expect_identical(saved, result$model_fit_history)
    expect_identical(saved$status, "accepted")
  }
})

test_that("background history preserves accepted, rejected and failed proposals and RNG", {
  withr::local_seed(82)
  args <- background_history_args()
  # Identify each proposal by the newly assigned regime, not worker-local state.
  args$fit <- function(IC, formula, tree, ...) {
    draw <- runif(2)
    regime <- max(as.integer(colnames(tree$mapped.edge)))
    if (regime == 3L) stop("fit failed after drawing")
    if (regime == 2L) warning("fit warning")
    list(GIC = list(GIC = c(90, 88)[regime]), draw = draw)
  }
  run <- function(progress, backend) {
    directory <- withr::local_tempdir()
    set.seed(82)
    messages <- character()
    warnings <- character()
    result <- withCallingHandlers(do.call(.bifrost_search_forward, c(args, list(
      sub_dir = directory, heartbeat = if (progress) function() NULL else NULL,
      is_rstudio = backend, tick = function(message) messages <<- c(messages, message)
    ))), warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    })
    files <- list.files(directory, pattern = "\\.rds$", full.names = TRUE)
    list(result = result, seed = .Random.seed, kind = RNGkind(),
         history = lapply(files, readRDS), files = basename(files),
         messages = messages, warnings = warnings)
  }
  expected <- run(FALSE, TRUE)
  expect_identical(vapply(expected$history, `[[`, "", "status"),
                   c("accepted", "rejected", "error"))
  expect_identical(expected$files, paste0("iteration_", 2:4, ".rds"))
  for (history in c(TRUE, FALSE)) {
    args$store_model_fit_history <- history
    for (backend in c(TRUE, if (future::supportsMulticore()) FALSE)) {
      for (progress in c(FALSE, TRUE)) {
        actual <- run(progress, backend)
        expect_identical(actual$seed, expected$seed)
        expect_identical(actual$kind, expected$kind)
        expect_identical(actual$messages, expected$messages)
        expect_identical(actual$warnings, expected$warnings)
        if (history) {
          expect_identical(actual, expected)
        } else {
          expect_length(actual$history, 0L)
          expect_identical(actual$result$model_fit_history, list())
          actual$result$model_fit_history <- expected$result$model_fit_history
          expect_identical(actual$result, expected$result)
        }
      }
    }
  }
})

test_that("history I/O failures escape after the proposal tick without becoming fit failures", {
  args <- background_history_args()
  args$sorted_candidates <- args$sorted_candidates[1]
  args$fit <- function(...) { runif(1); list(GIC = list(GIC = 90)) }
  withr::local_seed(82)
  for (progress in c(FALSE, TRUE)) {
    # A file cannot be used as a history directory: exercise real saveRDS errors.
    directory <- withr::local_tempfile()
    writeLines("not a directory", directory)
    messages <- character()
    warnings <- character()
    set.seed(82)
    error <- withCallingHandlers(tryCatch(
      do.call(.bifrost_search_forward, c(args, list(
        sub_dir = directory, is_rstudio = TRUE,
        heartbeat = if (progress) function() NULL else NULL,
        tick = function(message) messages <<- c(messages, message)
      ))), error = identity
    ), warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    })
    expect_s3_class(error, "error")
    expect_match(conditionMessage(error), "cannot open")
    expect_length(messages, 1L)
    expect_match(messages, "accepted")
    expect_length(warnings, 1L)
    expect_false(grepl("evaluating shift", warnings, fixed = TRUE))
    if (!progress) expected_seed <- .Random.seed
    else expect_identical(.Random.seed, expected_seed)
  }
})

test_that("a failed history write retains the fitted value and ordered I/O conditions", {
  directory <- withr::local_tempdir()
  blocked <- file.path(directory, "not-a-directory")
  writeLines("blocked", blocked)
  fit <- function() list(GIC = list(GIC = 90), draw = runif(1))
  withr::local_seed(82)
  run <- function(destination) {
    set.seed(82)
    .bifrost_search_fit_proposal(
      fit, list(), "GIC", 100, 5,
      history = list(step = 1L, candidate_node = 14L, regime_id = "1"),
      sub_dir = destination
    )
  }
  successful <- run(directory)
  seed <- .Random.seed
  failed <- run(blocked)
  expect_identical(.Random.seed, seed)
  expect_identical(failed$value, successful$value)
  expect_identical(readRDS(file.path(directory, "iteration_2.rds")), successful$value)
  expect_length(successful$conditions, 0L)
  expect_length(failed$conditions, 2L)
  expect_s3_class(failed$conditions[[1L]], "warning")
  expect_s3_class(failed$conditions[[2L]], "error")
})
