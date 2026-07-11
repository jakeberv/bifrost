.collect_progress_handler <- function(events) {
  record <- function(config, state, progression, ...) {
    events$types <- c(events$types, progression$type)
    events$messages <- c(events$messages, paste(state$message, collapse = ""))
    events$amounts <- c(events$amounts, progression$amount)
    invisible(NULL)
  }

  progressr::make_progression_handler(
    "bifrost-test-collector",
    reporter = list(
      reset = function(...) invisible(NULL),
      hide = function(...) invisible(NULL),
      unhide = function(...) invisible(NULL),
      initiate = record,
      update = record,
      finish = record,
      interrupt = record
    ),
    enable = TRUE
  )
}

.recording_search_renderer <- function(events) {
  list(
    create = function(label, total) {
      id <- paste0("row-", length(events$created) + 1L)
      events$created <- c(events$created, label)
      list(id = id, label = label, total = total)
    },
    update = function(row, state, current, total, status, force = TRUE) {
      events$updates <- c(events$updates, list(list(
        id = row$id,
        state = state,
        current = current,
        total = total,
        status = status
      )))
      invisible(NULL)
    },
    output = function(row, text) {
      events$output <- c(events$output, text)
      invisible(NULL)
    },
    done = function(row, result) {
      events$done <- c(events$done, row$id)
      invisible(NULL)
    }
  )
}

test_that("search status rows appear progressively and finalize together", {
  events <- new.env(parent = emptyenv())
  events$created <- character()
  events$updates <- list()
  events$output <- character()
  events$done <- character()
  session <- .bifrost_search_progress_session(
    TRUE,
    .recording_search_renderer(events)
  )
  testthat::expect_false(session$has_rows())

  run_stage <- function(label) {
    .bifrost_search_run_stage(
      enabled = TRUE,
      steps = 1L,
      initial_message = label,
      work = function(tick) {
        tick(message = paste(label, "complete"))
        list(value = label, done = paste(label, "done"))
      },
      session = session
    )
  }

  run_stage("[1/3] Candidate scoring")
  testthat::expect_true(session$has_rows())
  testthat::expect_length(events$created, 1L)
  testthat::expect_length(events$done, 0L)
  session$output("verbose between stages")

  run_stage("[2/3] Greedy shift search")
  testthat::expect_length(events$created, 2L)
  testthat::expect_length(events$done, 0L)

  session$skip("[3/3] IC-weight re-estimation", "not requested")
  testthat::expect_length(events$created, 3L)
  testthat::expect_equal(tail(events$updates, 1L)[[1L]]$state, "skipped")
  testthat::expect_identical(session$rows(), c(
    "[1/3] Candidate scoring",
    "[2/3] Greedy shift search",
    "[3/3] IC-weight re-estimation"
  ))
  testthat::expect_length(events$done, 0L)

  session$finalize()
  session$finalize()
  testthat::expect_length(events$done, 3L)
  testthat::expect_identical(events$done, paste0("row-", 1:3))
  testthat::expect_identical(events$output, "verbose between stages")

  event_counts <- c(
    created = length(events$created),
    updates = length(events$updates),
    output = length(events$output),
    done = length(events$done)
  )
  session$update_stage(
    "[4/4] Too late",
    list(max_steps = 2L),
    list(step = 1L, message = "ignored update"),
    list(amount = 1),
    "active"
  )
  session$skip("[4/4] Too late", "already finalized")
  session$output("ignored output")
  session$finalize()
  testthat::expect_identical(c(
    created = length(events$created),
    updates = length(events$updates),
    output = length(events$output),
    done = length(events$done)
  ), event_counts)
})

test_that("disabled search status session never touches its renderer", {
  events <- new.env(parent = emptyenv())
  events$created <- character()
  events$updates <- list()
  events$output <- character()
  events$done <- character()
  session <- .bifrost_search_progress_session(
    FALSE,
    .recording_search_renderer(events)
  )
  handler <- session$handler("disabled stage")
  reporter <- get("reporter", envir = environment(handler))
  config <- list(max_steps = 2L)
  state <- list(step = 1L, message = "ignored")
  progression <- list(amount = 1)

  testthat::expect_false(session$has_rows())
  testthat::expect_identical(session$rows(), character())
  reporter$initiate(config, state, progression)
  reporter$update(config, state, progression)
  reporter$finish(config, state, progression)
  reporter$interrupt(
    config,
    state,
    structure(
      list(message = "ignored failure", call = NULL, amount = 1),
      class = c("simpleError", "error", "condition")
    )
  )
  session$update_stage(
    "disabled direct stage",
    config,
    state,
    progression,
    "active"
  )
  session$skip("disabled skipped stage", "not requested")
  session$output("disabled output")
  session$finalize()
  session$finalize()

  testthat::expect_false(session$has_rows())
  testthat::expect_identical(session$rows(), character())
  testthat::expect_length(events$created, 0L)
  testthat::expect_length(events$updates, 0L)
  testthat::expect_length(events$output, 0L)
  testthat::expect_length(events$done, 0L)
})

test_that("search status session marks failures before rethrowing them", {
  events <- new.env(parent = emptyenv())
  events$created <- character()
  events$updates <- list()
  events$output <- character()
  events$done <- character()
  session <- .bifrost_search_progress_session(
    TRUE,
    .recording_search_renderer(events)
  )

  testthat::expect_error(
    .bifrost_search_run_stage(
      enabled = TRUE,
      steps = 2L,
      initial_message = "failing stage",
      work = function(tick) {
        tick(message = "one complete")
        stop("synthetic session failure")
      },
      handler = session$handler("failing stage")
    ),
    "synthetic session failure"
  )

  testthat::expect_true(any(vapply(
    events$updates,
    function(event) identical(event$state, "failed"),
    logical(1)
  )))
  testthat::expect_length(events$done, 0L)
  session$finalize()
  session$finalize()
  testthat::expect_identical(events$done, "row-1")
})

test_that("search CLI renderer can persist and finalize a skipped row", {
  old_cli_dynamic <- options(cli.dynamic = TRUE)
  on.exit(options(old_cli_dynamic), add = TRUE)
  session <- .bifrost_search_progress_session(TRUE)

  testthat::expect_no_error(
    rendered <- utils::capture.output({
      session$skip("[3/3] IC-weight re-estimation", "not requested")
      session$finalize()
    }, type = "message")
  )
  testthat::expect_true(any(grepl("[3/3] IC-weight", rendered, fixed = TRUE)))
})

test_that("search CLI renderer finalizes persistent rows as separate lines", {
  old_cli_dynamic <- options(cli.dynamic = TRUE)
  on.exit(options(old_cli_dynamic), add = TRUE)
  session <- .bifrost_search_progress_session(TRUE)

  rendered <- utils::capture.output({
    session$skip("[1/2] First stage", "not requested")
    session$skip("[2/2] Second stage", "not requested")
    session$finalize()
  }, type = "message")

  testthat::expect_gte(length(rendered), 2L)
})

test_that("search CLI renderer grows a live multi-line stage stack", {
  old_cli_options <- options(
    cli.dynamic = TRUE,
    cli.num_colors = 1,
    cli.width = 200
  )
  on.exit(options(old_cli_options), add = TRUE)
  session <- .bifrost_search_progress_session(TRUE)
  on.exit(session$finalize(), add = TRUE)

  run_stage <- function(label) {
    .bifrost_search_run_stage(
      enabled = TRUE,
      steps = 1L,
      initial_message = label,
      work = function(tick) {
        tick(message = paste(label, "complete"))
        list(value = NULL, done = paste(label, "done"))
      },
      session = session
    )
  }

  rendered <- utils::capture.output({
    run_stage("[1/3] Candidate scoring")
    run_stage("[2/3] Greedy shift search")
  }, type = "message")

  testthat::expect_match(
    paste(rendered, collapse = "\n"),
    "\033[1A",
    fixed = TRUE
  )
})

test_that("search CLI renderer redraws active rows below verbose output", {
  old_cli_options <- options(cli.dynamic = TRUE, cli.width = 200)
  on.exit(options(old_cli_options), add = TRUE)
  session <- .bifrost_search_progress_session(TRUE)
  on.exit(session$finalize(), add = TRUE)

  .bifrost_search_run_stage(
    enabled = TRUE,
    steps = 1L,
    initial_message = "[1/3] Candidate scoring",
    work = function(tick) {
      tick(message = "[1/3] Candidate scoring complete")
      list(value = NULL, done = "[1/3] Candidate scoring done")
    },
    session = session
  )
  .bifrost_search_run_stage(
    enabled = TRUE,
    steps = 1L,
    initial_message = "[2/3] Greedy shift search",
    work = function(tick) {
      tick(message = "[2/3] Greedy shift search complete")
      list(value = NULL, done = "[2/3] Greedy shift search done")
    },
    session = session
  )

  rendered <- utils::capture.output(
    session$output("verbose detail"),
    type = "message"
  )

  plain <- gsub(
    "\r",
    "",
    cli::ansi_strip(paste(rendered, collapse = "\n")),
    fixed = TRUE
  )
  verbose_pos <- max(gregexpr("verbose detail", plain, fixed = TRUE)[[1L]])
  stage1_pos <- max(gregexpr("[1/3] Candidate scoring", plain, fixed = TRUE)[[1L]])
  stage2_pos <- max(gregexpr("[2/3] Greedy shift search", plain, fixed = TRUE)[[1L]])
  testthat::expect_gt(verbose_pos, 0L)
  testthat::expect_gt(stage1_pos, verbose_pos)
  testthat::expect_gt(stage2_pos, verbose_pos)
})

test_that("search progress stage runner reports lifecycle and returns work value", {
  events <- new.env(parent = emptyenv())
  events$types <- character()
  events$messages <- character()

  value <- .bifrost_search_run_stage(
    enabled = TRUE,
    steps = 2L,
    initial_message = "[1/3] Candidate scoring",
    work = function(tick) {
      tick(message = "candidate 1")
      tick(message = "candidate 2")
      list(value = 42L, done = "[1/3] Candidate scoring complete")
    },
    handler = .collect_progress_handler(events)
  )

  testthat::expect_identical(value, 42L)
  testthat::expect_identical(events$types[[1L]], "initiate")
  testthat::expect_equal(sum(events$types == "update"), 2L)
  testthat::expect_gte(sum(events$types == "finish"), 1L)
  testthat::expect_true(any(grepl("Candidate scoring complete", events$messages, fixed = TRUE)))
})

test_that("search progress stage runner owns and finalizes its default session", {
  events <- new.env(parent = emptyenv())
  events$constructed <- 0L
  events$created <- character()
  events$updates <- list()
  events$output <- character()
  events$done <- character()
  progress_session <- .bifrost_search_progress_session
  testthat::local_mocked_bindings(
    .bifrost_search_progress_session = function(enabled) {
      events$constructed <- events$constructed + 1L
      progress_session(enabled, .recording_search_renderer(events))
    }
  )

  value <- .bifrost_search_run_stage(
    enabled = TRUE,
    steps = 1L,
    initial_message = "owned stage",
    work = function(tick) {
      tick(message = "owned stage complete")
      list(value = 73L, done = "owned stage done")
    }
  )

  testthat::expect_identical(value, 73L)
  testthat::expect_identical(events$constructed, 1L)
  testthat::expect_identical(events$created, "owned stage")
  testthat::expect_identical(events$done, "row-1")
})

test_that("a heartbeat after the final completion remains valid until stage finish", {
  events <- new.env(parent = emptyenv())
  events$types <- character()
  events$messages <- character()
  events$amounts <- numeric()

  value <- .bifrost_search_run_stage(
    enabled = TRUE,
    steps = 1L,
    initial_message = "boundary stage",
    work = function(tick) {
      tick(message = "fit complete")
      tick(amount = 0, message = "collecting result")
      list(value = 11L, done = "boundary stage complete")
    },
    handler = .collect_progress_handler(events)
  )

  testthat::expect_identical(value, 11L)
  testthat::expect_identical(events$types[[1L]], "initiate")
  testthat::expect_equal(sum(events$types == "update"), 2L)
  testthat::expect_gte(sum(events$types == "finish"), 1L)
})

test_that("disabled search progress stage runs without creating a handler", {
  value <- .bifrost_search_run_stage(
    enabled = FALSE,
    steps = 1L,
    initial_message = "disabled",
    work = function(tick) {
      tick(message = "ignored")
      list(value = "result", done = "done")
    },
    handler = stop("handler must remain lazy when progress is disabled")
  )

  testthat::expect_identical(value, "result")
})

test_that("search progress stage preserves failure lifecycle and rethrows errors", {
  events <- new.env(parent = emptyenv())
  events$types <- character()
  events$messages <- character()

  testthat::expect_error(
    .bifrost_search_run_stage(
      enabled = TRUE,
      steps = 2L,
      initial_message = "failing stage",
      work = function(tick) {
        tick(message = "one complete")
        stop("synthetic stage failure")
      },
      handler = .collect_progress_handler(events)
    ),
    "synthetic stage failure"
  )

  testthat::expect_true("interrupt" %in% events$types)
  testthat::expect_false(any(events$types %in% c("finish", "shutdown")))
})

test_that("animated Future failures cancel remaining work and restore settings", {
  thread_vars <- c(
    "OMP_NUM_THREADS",
    "OPENBLAS_NUM_THREADS",
    "MKL_NUM_THREADS",
    "VECLIB_MAXIMUM_THREADS",
    "NUMEXPR_NUM_THREADS"
  )
  old_threads <- Sys.getenv(thread_vars, unset = NA_character_)
  old_plan <- future::plan("list")
  on.exit({
    future::plan(old_plan)
    for (nm in thread_vars) {
      value <- old_threads[[nm]]
      if (is.na(value)) {
        Sys.unsetenv(nm)
      } else {
        do.call(Sys.setenv, stats::setNames(list(value), nm))
      }
    }
  }, add = TRUE)

  future::plan(future::sequential)
  Sys.setenv(
    OMP_NUM_THREADS = "7",
    OPENBLAS_NUM_THREADS = "8",
    MKL_NUM_THREADS = "9",
    VECLIB_MAXIMUM_THREADS = "10",
    NUMEXPR_NUM_THREADS = "11"
  )

  started <- proc.time()[["elapsed"]]
  testthat::expect_error(
    .bifrost_search_future_lapply(
      1:2,
      function(i) {
        if (i == 1L) {
          Sys.sleep(0.15)
          stop("synthetic worker failure")
        }
        Sys.sleep(5)
        i
      },
      workers = 2L,
      is_rstudio_flag = TRUE,
      heartbeat = function() invisible(NULL),
      interval = 0.05
    ),
    "synthetic worker failure"
  )
  elapsed <- proc.time()[["elapsed"]] - started

  if (!identical(Sys.getenv("R_COVR"), "true")) {
    testthat::expect_lt(elapsed, 3)
  }
  testthat::expect_s3_class(future::plan("next"), "sequential")
  testthat::expect_identical(Sys.getenv("OMP_NUM_THREADS"), "7")
  testthat::expect_identical(Sys.getenv("OPENBLAS_NUM_THREADS"), "8")
  testthat::expect_identical(Sys.getenv("MKL_NUM_THREADS"), "9")
  testthat::expect_identical(Sys.getenv("VECLIB_MAXIMUM_THREADS"), "10")
  testthat::expect_identical(Sys.getenv("NUMEXPR_NUM_THREADS"), "11")
})

test_that("animated Future RNG is caller-safe and independent of chunk layout", {
  run_random_work <- function(workers) {
    set.seed(42)
    .bifrost_search_future_lapply(
      1:6,
      function(i) stats::runif(1),
      workers = workers,
      is_rstudio_flag = TRUE,
      heartbeat = function() invisible(NULL),
      interval = 0.01
    )
  }

  set.seed(101)
  seed_before <- .Random.seed
  deterministic <- .bifrost_search_future_lapply(
    1:3,
    identity,
    workers = 2L,
    is_rstudio_flag = TRUE,
    heartbeat = function() invisible(NULL),
    interval = 0.01
  )

  testthat::expect_identical(deterministic, as.list(1:3))
  testthat::expect_identical(.Random.seed, seed_before)
  testthat::expect_equal(run_random_work(1L), run_random_work(2L))
  testthat::expect_equal(run_random_work(2L), run_random_work(3L))
})

test_that("verbose message and plotting-style stdout coexist with stage progress", {
  events <- new.env(parent = emptyenv())
  events$types <- character()
  events$messages <- character()
  output <- character()

  messages <- testthat::capture_messages(
    output <- testthat::capture_output(
      value <- .bifrost_search_run_stage(
        enabled = TRUE,
        steps = 1L,
        initial_message = "plotting stage",
        work = function(tick) {
          cat("plotting-style verbose output\n")
          message("ordinary verbose message")
          tick(message = "fit complete")
          list(value = 7L, done = "plotting stage complete")
        },
        handler = .collect_progress_handler(events)
      )
    )
  )

  testthat::expect_identical(value, 7L)
  testthat::expect_match(paste(output, collapse = "\n"), "plotting-style verbose output", fixed = TRUE)
  testthat::expect_match(messages, "ordinary verbose message", fixed = TRUE)
  testthat::expect_true("finish" %in% events$types)
})

test_that("search CLI handler redraws zero-increment heartbeat events", {
  events <- new.env(parent = emptyenv())
  events$created <- character()
  events$updates <- list()
  events$output <- character()
  events$done <- character()
  session <- .bifrost_search_progress_session(
    TRUE,
    .recording_search_renderer(events)
  )
  handler <- session$handler("heartbeat stage")
  reporter <- get("reporter", envir = environment(handler))
  config <- list(max_steps = 3L)

  reporter$initiate(
    config,
    list(step = 0L, message = "starting"),
    list(amount = 1)
  )
  reporter$update(
    config,
    list(step = 1L, message = "first complete"),
    list(amount = 1)
  )
  reporter$update(
    config,
    list(step = 2L, message = "heartbeat"),
    list(amount = 0)
  )

  heartbeat <- tail(events$updates, 1L)[[1L]]
  testthat::expect_identical(heartbeat$current, 1L)
  testthat::expect_identical(heartbeat$total, 2L)
  testthat::expect_identical(heartbeat$status, "heartbeat")
})

test_that("non-dynamic search output stays bounded across heartbeats", {
  old_cli_options <- options(cli.dynamic = FALSE, cli.width = 200)
  on.exit(options(old_cli_options), add = TRUE)
  session <- .bifrost_search_progress_session(TRUE)
  handler <- session$handler("non-dynamic heartbeat stage")
  reporter <- get("reporter", envir = environment(handler))
  config <- list(max_steps = 3L)

  rendered <- utils::capture.output({
    reporter$initiate(
      config,
      list(step = 0L, message = "starting"),
      list(amount = 1)
    )
    reporter$update(
      config,
      list(step = 1L, message = "one complete"),
      list(amount = 1)
    )
    for (i in seq_len(8L)) {
      reporter$update(
        config,
        list(step = 2L, message = paste("heartbeat", i)),
        list(amount = 0)
      )
    }
    reporter$finish(
      config,
      list(step = 3L, message = "stage complete"),
      list(amount = 1)
    )
    session$finalize()
  }, type = "message")

  testthat::expect_lte(length(rendered), 5L)
  testthat::expect_true(any(grepl("stage complete", rendered, fixed = TRUE)))
})

test_that("skipped search stages are persistent only when progress is enabled", {
  shown <- testthat::capture_messages(
    .bifrost_search_report_skipped_stage(
      TRUE,
      "[3/3] IC-weight re-estimation",
      "not requested"
    )
  )
  hidden <- testthat::capture_messages(
    .bifrost_search_report_skipped_stage(
      FALSE,
      "[3/3] IC-weight re-estimation",
      "not requested"
    )
  )

  testthat::expect_match(shown, "[3/3] IC-weight re-estimation", fixed = TRUE)
  testthat::expect_match(shown, "skipped: not requested", fixed = TRUE)
  testthat::expect_length(hidden, 0L)
})

test_that("search progress is a default-on public argument", {
  testthat::expect_identical(
    formals(searchOptimalConfiguration)$progress,
    TRUE
  )
})

test_that("public search keeps positional dots while progress is keyword-only", {
  testthat::skip_if_not_installed("ape")
  testthat::skip_if_not_installed("phytools")

  set.seed(20260715)
  tree <- ape::rtree(6)
  traits <- matrix(stats::rnorm(12), ncol = 2)
  rownames(traits) <- tree$tip.label
  forwarded_methods <- character()
  enabled_sessions <- logical()
  events <- new.env(parent = emptyenv())
  events$created <- character()
  events$updates <- list()
  events$output <- character()
  events$done <- character()
  progress_session <- .bifrost_search_progress_session

  testthat::local_mocked_bindings(
    generatePaintedTrees = function(tree, min_tips, state = "shift") {
      list(tree)
    },
    fitMvglsAndExtractGIC.formula = function(formula, tree, trait_data, ...) {
      dots <- list(...)
      forwarded_methods <<- c(
        forwarded_methods,
        if (length(dots) == 0L) "<missing>" else dots[[1L]]
      )
      list(
        model = list(corrSt = list(phy = tree)),
        GIC = list(GIC = 100)
      )
    },
    .bifrost_search_forward = function(
        sorted_candidates,
        current_best_tree,
        current_best_ic,
        ...) {
      list(
        current_best_tree = current_best_tree,
        current_best_ic = current_best_ic,
        shift_vec = list(),
        shifts_no_uncertainty = integer(),
        model_with_shift_no_uncertainty = NULL,
        best_tree_no_uncertainty = NULL,
        warnings_list = list(),
        outcome_counts = c(accepted = 0L, rejected = 0L, error = 0L),
        model_fit_history = list()
      )
    },
    .bifrost_search_progress_session = function(enabled) {
      enabled_sessions <<- c(enabled_sessions, enabled)
      progress_session(enabled, .recording_search_renderer(events))
    }
  )

  positional_args <- list(
    baseline_tree = tree,
    trait_data = traits,
    formula = "trait_data ~ 1",
    min_descendant_tips = ape::Ntip(tree),
    num_cores = 1,
    ic_uncertainty_threshold = 1,
    shift_acceptance_threshold = 1e9,
    uncertaintyweights = FALSE,
    uncertaintyweights_par = FALSE,
    plot = FALSE,
    IC = "GIC",
    store_model_fit_history = FALSE,
    verbose = FALSE
  )
  positional_args[[length(positional_args) + 1L]] <- "positional-method"

  default_result <- suppressWarnings(do.call(
    searchOptimalConfiguration,
    positional_args
  ))
  positional_args$progress <- FALSE
  disabled_result <- suppressWarnings(do.call(
    searchOptimalConfiguration,
    positional_args
  ))

  testthat::expect_identical(
    tail(names(formals(searchOptimalConfiguration)), 2L),
    c("...", "progress")
  )
  testthat::expect_identical(
    forwarded_methods,
    rep("positional-method", 2L)
  )
  testthat::expect_true(default_result$user_input$progress)
  testthat::expect_false(disabled_result$user_input$progress)
  testthat::expect_identical(enabled_sessions, c(TRUE, FALSE))
  testthat::expect_identical(events$created, c(
    "[1/3] Candidate scoring",
    "[2/3] Greedy shift search",
    "[3/3] IC-weight re-estimation"
  ))
})

test_that("dual weight flags fail before public search work or progress rows", {
  testthat::skip_if_not_installed("ape")
  testthat::skip_if_not_installed("phytools")

  set.seed(20260716)
  tree <- ape::rtree(6)
  traits <- matrix(stats::rnorm(12), ncol = 2)
  rownames(traits) <- tree$tip.label
  candidate_calls <- 0L
  fit_calls <- 0L
  events <- new.env(parent = emptyenv())
  events$created <- character()
  events$updates <- list()
  events$output <- character()
  events$done <- character()
  progress_session <- .bifrost_search_progress_session

  testthat::local_mocked_bindings(
    generatePaintedTrees = function(tree, min_tips, state = "shift") {
      candidate_calls <<- candidate_calls + 1L
      list(tree)
    },
    fitMvglsAndExtractGIC.formula = function(formula, tree, trait_data, ...) {
      fit_calls <<- fit_calls + 1L
      list(
        model = list(corrSt = list(phy = tree)),
        GIC = list(GIC = 100)
      )
    },
    .bifrost_search_forward = function(
        sorted_candidates,
        current_best_tree,
        current_best_ic,
        ...) {
      list(
        current_best_tree = current_best_tree,
        current_best_ic = current_best_ic,
        shift_vec = list(),
        shifts_no_uncertainty = integer(),
        model_with_shift_no_uncertainty = NULL,
        best_tree_no_uncertainty = NULL,
        warnings_list = list(),
        outcome_counts = c(accepted = 0L, rejected = 0L, error = 0L),
        model_fit_history = list()
      )
    },
    .bifrost_search_progress_session = function(enabled) {
      progress_session(enabled, .recording_search_renderer(events))
    }
  )

  testthat::expect_error(
    searchOptimalConfiguration(
      baseline_tree = tree,
      trait_data = traits,
      min_descendant_tips = ape::Ntip(tree),
      num_cores = 1,
      uncertaintyweights = TRUE,
      uncertaintyweights_par = TRUE,
      plot = FALSE,
      IC = "GIC",
      store_model_fit_history = FALSE,
      verbose = FALSE,
      progress = TRUE,
      method = "LL"
    ),
    "Exactly one of uncertaintyweights or uncertaintyweights_par must be TRUE.",
    fixed = TRUE
  )

  testthat::expect_identical(candidate_calls, 0L)
  testthat::expect_identical(fit_calls, 0L)
  testthat::expect_length(events$created, 0L)
  testthat::expect_length(events$updates, 0L)
  testthat::expect_length(events$output, 0L)
  testthat::expect_length(events$done, 0L)
})

test_that("search progress helpers cover empty work and seedless callers", {
  random_env <- globalenv()
  had_seed <- exists(".Random.seed", envir = random_env, inherits = FALSE)
  old_seed <- if (had_seed) get(".Random.seed", envir = random_env) else NULL
  old_plan <- future::plan("list")
  on.exit({
    future::plan(old_plan)
    if (had_seed) {
      assign(".Random.seed", old_seed, envir = random_env)
    } else if (exists(".Random.seed", envir = random_env, inherits = FALSE)) {
      rm(".Random.seed", envir = random_env)
    }
  }, add = TRUE)

  if (exists(".Random.seed", envir = random_env, inherits = FALSE)) {
    rm(".Random.seed", envir = random_env)
  }
  seeds <- .bifrost_search_rng_seeds(2L)
  testthat::expect_length(seeds, 2L)
  testthat::expect_false(exists(".Random.seed", envir = random_env, inherits = FALSE))

  random_value <- .bifrost_search_with_rng_seed(
    seeds[[1L]],
    function() stats::runif(1L)
  )
  testthat::expect_length(random_value, 1L)
  testthat::expect_false(exists(".Random.seed", envir = random_env, inherits = FALSE))

  future::plan(future::sequential)
  testthat::expect_identical(.bifrost_search_await_futures(list()), list())
  testthat::expect_identical(
    .bifrost_search_await_work(function() 12L, interval = 0.01),
    12L
  )
  testthat::expect_identical(
    .bifrost_search_future_lapply(
      integer(),
      identity,
      workers = 1L,
      is_rstudio_flag = TRUE,
      heartbeat = function() invisible(NULL)
    ),
    list()
  )
})

test_that("print.bifrost_search formats captured threshold expressions", {
  object <- structure(
    list(
      user_input = list(
        shift_acceptance_threshold = quote(-Inf),
        progress = TRUE
      ),
      shift_nodes_no_uncertainty = integer(),
      warnings = character()
    ),
    class = c("bifrost_search", "list")
  )

  testthat::expect_no_error(
    output <- paste(testthat::capture_output(print(object)), collapse = "\n")
  )
  testthat::expect_match(output, "Threshold: -Inf", fixed = TRUE)
})

.fake_candidate_fit <- function(IC, formula, tree, trait_data, ...) {
  list(GIC = list(GIC = tree$ic))
}

test_that("candidate scoring ticks once per completed fit in serial and multisession modes", {
  candidate_trees <- list(
    "Node 11" = list(ic = 95),
    "Node 12" = list(ic = 90),
    "Node 13" = list(ic = 97)
  )

  run_scoring <- function(num_cores, is_rstudio) {
    events <- new.env(parent = emptyenv())
    events$types <- character()
    events$messages <- character()
    result <- .bifrost_search_run_stage(
      enabled = TRUE,
      steps = length(candidate_trees),
      initial_message = "[1/3] Candidate scoring",
      work = function(tick) {
        value <- .bifrost_search_score_candidates(
          candidate_trees_shifts = candidate_trees,
          baseline_ic = 100,
          IC = "GIC",
          formula = trait_data ~ 1,
          trait_data = matrix(0, nrow = 1),
          args_list = list(),
          num_cores = num_cores,
          is_rstudio = is_rstudio,
          tick = tick,
          fit = .fake_candidate_fit
        )
        list(value = value, done = "[1/3] Candidate scoring complete")
      },
      handler = .collect_progress_handler(events)
    )
    updates <- events$messages[events$types == "update"]
    list(result = result, messages = updates)
  }

  serial <- run_scoring(num_cores = 1L, is_rstudio = FALSE)
  parallel <- run_scoring(num_cores = 2L, is_rstudio = TRUE)

  testthat::expect_identical(names(serial$result$sorted_candidates), c("Node 12", "Node 11", "Node 13"))
  testthat::expect_length(serial$messages, 3L)
  testthat::expect_length(parallel$messages, 3L)
  testthat::expect_setequal(serial$messages, paste("[1/3] Candidate scoring - completed", names(candidate_trees)))
  testthat::expect_setequal(parallel$messages, serial$messages)
})

test_that("candidate scoring emits heartbeats without advancing completion", {
  candidate_trees <- list("Node 11" = list(ic = 95))
  events <- new.env(parent = emptyenv())
  events$types <- character()
  events$messages <- character()
  events$amounts <- numeric()

  slow_fit <- function(IC, formula, tree, trait_data, ...) {
    Sys.sleep(0.3)
    list(GIC = list(GIC = tree$ic))
  }

  .bifrost_search_run_stage(
    enabled = TRUE,
    steps = 1L,
    initial_message = "[1/3] Candidate scoring",
    work = function(tick) {
      value <- .bifrost_search_score_candidates(
        candidate_trees_shifts = candidate_trees,
        baseline_ic = 100,
        IC = "GIC",
        formula = trait_data ~ 1,
        trait_data = matrix(0, nrow = 1),
        args_list = list(),
        num_cores = 1L,
        is_rstudio = TRUE,
        tick = tick,
        heartbeat = function() {
          tick(amount = 0, message = "[1/3] Candidate scoring - fitting")
        },
        fit = slow_fit
      )
      list(value = value, done = "[1/3] Candidate scoring complete")
    },
    handler = .collect_progress_handler(events)
  )

  update_amounts <- events$amounts[events$types == "update"]
  testthat::expect_equal(sum(update_amounts > 0), 1L)
  testthat::expect_gte(sum(update_amounts == 0), 2L)
})

test_that("public search shares progress rows with verbose output", {
  testthat::skip_if_not_installed("ape")
  testthat::skip_if_not_installed("phytools")
  testthat::skip_if_not_installed("mvMORPH")

  events <- new.env(parent = emptyenv())
  events$created <- character()
  events$updates <- list()
  events$output <- character()
  events$done <- character()
  progress_session <- .bifrost_search_progress_session
  testthat::local_mocked_bindings(
    .bifrost_search_progress_session = function(enabled) {
      progress_session(enabled, .recording_search_renderer(events))
    }
  )

  set.seed(20260711)
  tree <- ape::rtree(10)
  traits <- matrix(stats::rnorm(20), ncol = 2)
  rownames(traits) <- tree$tip.label

  skipped <- testthat::capture_messages(
    result <- searchOptimalConfiguration(
      baseline_tree = tree,
      trait_data = traits,
      min_descendant_tips = 4,
      num_cores = 1,
      shift_acceptance_threshold = 1e9,
      uncertaintyweights = FALSE,
      uncertaintyweights_par = FALSE,
      plot = FALSE,
      IC = "GIC",
      store_model_fit_history = FALSE,
      verbose = TRUE,
      method = "LL"
    )
  )

  testthat::expect_s3_class(result, "bifrost_search")
  testthat::expect_identical(result$user_input$progress, TRUE)
  testthat::expect_identical(events$created, c(
    "[1/3] Candidate scoring",
    "[2/3] Greedy shift search",
    "[3/3] IC-weight re-estimation"
  ))
  testthat::expect_true(any(grepl(
    "Sorting and evaluating shifts",
    events$output,
    fixed = TRUE
  )))
  testthat::expect_identical(events$done, paste0("row-", 1:3))
  testthat::expect_true(any(grepl(
    "Generating candidate shift models",
    skipped,
    fixed = TRUE
  )))
})

test_that("greedy search ticks for accepted, rejected, and recoverable-error fits", {
  testthat::skip_if_not_installed("ape")
  testthat::skip_if_not_installed("phytools")

  set.seed(44)
  tree <- ape::rtree(12)
  baseline <- phytools::paintSubTree(
    tree,
    node = ape::Ntip(tree) + 1L,
    state = 0
  )
  candidates <- generatePaintedTrees(baseline, min_tips = 3)[-1]
  candidates <- candidates[seq_len(min(3L, length(candidates)))]
  testthat::skip_if(length(candidates) < 3L, "fixture needs three candidate nodes")

  fit_number <- 0L
  fake_fit <- function(IC, formula, tree, trait_data, ...) {
    fit_number <<- fit_number + 1L
    if (fit_number == 3L) stop("synthetic fit failure")
    ic <- c(90, 88)[fit_number]
    list(GIC = list(GIC = ic), model = list(corrSt = list(phy = tree)))
  }
  messages <- character()

  result <- suppressWarnings(.bifrost_search_forward(
    sorted_candidates = candidates,
    current_best_tree = baseline,
    current_best_ic = 100,
    shift_id = 0,
    IC = "GIC",
    formula = trait_data ~ 1,
    trait_data = matrix(0, nrow = ape::Ntip(tree), ncol = 1),
    shift_acceptance_threshold = 5,
    store_model_fit_history = FALSE,
    sub_dir = NULL,
    plot = FALSE,
    verbose_log = function(...) invisible(NULL),
    tick = function(...) messages <<- c(messages, list(...)$message),
    fit = fake_fit
  ))

  testthat::expect_length(messages, 3L)
  testthat::expect_match(messages[1], "accepted", fixed = TRUE)
  testthat::expect_match(messages[2], "rejected", fixed = TRUE)
  testthat::expect_match(messages[3], "error", fixed = TRUE)
  testthat::expect_length(result$shifts_no_uncertainty, 1L)
  testthat::expect_length(result$warnings_list, 1L)
})

test_that("greedy search emits heartbeats while one proposal fit is running", {
  testthat::skip_if_not_installed("ape")
  testthat::skip_if_not_installed("phytools")

  set.seed(46)
  tree <- ape::rtree(10)
  baseline <- phytools::paintSubTree(
    tree,
    node = ape::Ntip(tree) + 1L,
    state = 0
  )
  candidate <- generatePaintedTrees(baseline, min_tips = 3)[2]
  events <- new.env(parent = emptyenv())
  events$types <- character()
  events$messages <- character()
  events$amounts <- numeric()

  slow_fit <- function(IC, formula, tree, trait_data, ...) {
    Sys.sleep(0.3)
    list(GIC = list(GIC = 99))
  }

  .bifrost_search_run_stage(
    enabled = TRUE,
    steps = 1L,
    initial_message = "[2/3] Greedy shift search",
    work = function(tick) {
      value <- .bifrost_search_forward(
        sorted_candidates = candidate,
        current_best_tree = baseline,
        current_best_ic = 100,
        shift_id = 0,
        IC = "GIC",
        formula = trait_data ~ 1,
        trait_data = matrix(0, nrow = ape::Ntip(tree), ncol = 1),
        shift_acceptance_threshold = 5,
        store_model_fit_history = FALSE,
        sub_dir = NULL,
        plot = FALSE,
        verbose_log = function(...) invisible(NULL),
        tick = tick,
        heartbeat = function() {
          tick(amount = 0, message = "[2/3] Greedy shift search - fitting")
        },
        is_rstudio = TRUE,
        fit = slow_fit
      )
      list(value = value, done = "[2/3] Greedy shift search complete")
    },
    handler = .collect_progress_handler(events)
  )

  update_amounts <- events$amounts[events$types == "update"]
  testthat::expect_equal(sum(update_amounts > 0), 1L)
  testthat::expect_gte(sum(update_amounts == 0), 2L)
})

test_that("IC-weight re-estimation ticks once per completed serial and parallel refit", {
  testthat::skip_if_not_installed("ape")
  testthat::skip_if_not_installed("phytools")
  testthat::skip_if_not_installed("mvMORPH")

  set.seed(45)
  tree <- ape::rtree(12)
  shifted_tree <- phytools::paintSubTree(
    tree,
    node = ape::Ntip(tree) + 1L,
    state = 0
  )
  candidate_names <- names(generatePaintedTrees(shifted_tree, min_tips = 3))[-1]
  shift_nodes <- as.integer(sub("Node ", "", candidate_names[seq_len(2L)]))
  shift_id <- 0L
  for (node in shift_nodes) {
    added <- addShiftToModel(shifted_tree, node, shift_id)
    shifted_tree <- added$tree
    shift_id <- added$shift_id
  }

  fake_fit <- function(IC, formula, tree, trait_data, ...) {
    Sys.sleep(0.3)
    list(GIC = list(GIC = 90))
  }

  run_weights <- function(parallel) {
    events <- new.env(parent = emptyenv())
    events$types <- character()
    events$messages <- character()
    events$amounts <- numeric()
    value <- .bifrost_search_run_stage(
      enabled = TRUE,
      steps = length(shift_nodes),
      initial_message = "[3/3] IC-weight re-estimation",
      work = function(tick) {
        weights <- .bifrost_search_calculate_ic_weights(
          uncertaintyweights = !parallel,
          uncertaintyweights_par = parallel,
          shift_vec = as.list(shift_nodes),
          best_tree_no_uncertainty = shifted_tree,
          model_with_shift_no_uncertainty = list(GIC = list(GIC = 80)),
          IC = "GIC",
          formula = trait_data ~ 1,
          trait_data = matrix(0, nrow = ape::Ntip(tree), ncol = 1),
          args_list = list(),
          num_cores = if (parallel) 2L else 1L,
          is_rstudio = parallel,
          verbose_log = function(...) invisible(NULL),
          tick = tick,
          heartbeat = function() {
            tick(amount = 0, message = "[3/3] IC-weight re-estimation - fitting")
          },
          fit = fake_fit
        )
        list(
          value = weights,
          done = "[3/3] IC-weight re-estimation complete"
        )
      },
      handler = .collect_progress_handler(events)
    )
    list(
      value = value,
      messages = events$messages[events$types == "update" & events$amounts > 0],
      amounts = events$amounts[events$types == "update"]
    )
  }

  serial <- run_weights(FALSE)
  parallel <- run_weights(TRUE)

  expected <- paste("[3/3] IC-weight re-estimation - completed node", shift_nodes)
  testthat::expect_equal(nrow(serial$value), 2L)
  testthat::expect_equal(nrow(parallel$value), 2L)
  testthat::expect_setequal(serial$messages, expected)
  testthat::expect_setequal(parallel$messages, expected)
  testthat::expect_equal(sum(serial$amounts > 0), 2L)
  testthat::expect_equal(sum(parallel$amounts > 0), 2L)
  testthat::expect_gte(sum(serial$amounts == 0), 2L)
  testthat::expect_gte(sum(parallel$amounts == 0), 2L)
})

test_that("serial IC weights retain the direct no-heartbeat path", {
  testthat::skip_if_not_installed("ape")
  testthat::skip_if_not_installed("phytools")

  set.seed(47)
  tree <- ape::rtree(10)
  shifted_tree <- phytools::paintSubTree(
    tree,
    node = ape::Ntip(tree) + 1L,
    state = 0
  )
  shift_node <- as.integer(sub(
    "Node ",
    "",
    names(generatePaintedTrees(shifted_tree, min_tips = 3))[2L]
  ))
  shifted_tree <- addShiftToModel(shifted_tree, shift_node, 0L)$tree
  ticks <- character()

  weights <- .bifrost_search_calculate_ic_weights(
    uncertaintyweights = TRUE,
    uncertaintyweights_par = FALSE,
    shift_vec = list(shift_node),
    best_tree_no_uncertainty = shifted_tree,
    model_with_shift_no_uncertainty = list(GIC = list(GIC = 80)),
    IC = "GIC",
    formula = trait_data ~ 1,
    trait_data = matrix(0, nrow = ape::Ntip(tree), ncol = 1),
    args_list = list(),
    num_cores = 1L,
    is_rstudio = FALSE,
    verbose_log = function(...) invisible(NULL),
    tick = function(...) ticks <<- c(ticks, list(...)$message),
    heartbeat = NULL,
    fit = function(IC, formula, tree, trait_data, ...) {
      list(GIC = list(GIC = 90))
    }
  )

  testthat::expect_equal(nrow(weights), 1L)
  testthat::expect_length(ticks, 1L)
})

test_that("public search renders progress through accepted-shift weight re-estimation", {
  testthat::skip_if_not_installed("ape")
  testthat::skip_if_not_installed("phytools")
  testthat::skip_if_not_installed("mvMORPH")

  events <- new.env(parent = emptyenv())
  events$created <- character()
  events$updates <- list()
  events$output <- character()
  events$done <- character()
  fit_count <- 0L
  progress_session <- .bifrost_search_progress_session
  testthat::local_mocked_bindings(
    fitMvglsAndExtractGIC.formula = function(formula, tree, trait_data, ...) {
      fit_count <<- fit_count + 1L
      list(
        model = list(corrSt = list(phy = tree)),
        GIC = list(GIC = 1000 - fit_count)
      )
    },
    removeShiftFromTree = function(tree, shift_node, stem = FALSE) tree,
    .bifrost_search_progress_session = function(enabled) {
      progress_session(enabled, .recording_search_renderer(events))
    }
  )

  set.seed(20260714)
  tree <- ape::rtree(8)
  traits <- matrix(stats::rnorm(16), ncol = 2)
  rownames(traits) <- tree$tip.label

  result <- suppressWarnings(searchOptimalConfiguration(
    baseline_tree = tree,
    trait_data = traits,
    min_descendant_tips = 3,
    num_cores = 1,
    shift_acceptance_threshold = -Inf,
    uncertaintyweights = TRUE,
    uncertaintyweights_par = FALSE,
    plot = FALSE,
    IC = "GIC",
    store_model_fit_history = FALSE,
    verbose = FALSE,
    progress = TRUE,
    method = "LL"
  ))

  testthat::expect_gt(length(result$shift_nodes_no_uncertainty), 0L)
  testthat::expect_equal(
    nrow(result$ic_weights),
    length(result$shift_nodes_no_uncertainty)
  )
  testthat::expect_identical(events$created, c(
    "[1/3] Candidate scoring",
    "[2/3] Greedy shift search",
    "[3/3] IC-weight re-estimation"
  ))
  testthat::expect_identical(events$done, paste0("row-", 1:3))
})

test_that("requested IC-weight stage reports why zero-shift work is skipped", {
  testthat::skip_if_not_installed("ape")
  testthat::skip_if_not_installed("phytools")
  testthat::skip_if_not_installed("mvMORPH")

  set.seed(20260712)
  tree <- ape::rtree(10)
  traits <- matrix(stats::rnorm(20), ncol = 2)
  rownames(traits) <- tree$tip.label
  events <- new.env(parent = emptyenv())
  events$created <- character()
  events$updates <- list()
  events$output <- character()
  events$done <- character()
  progress_session <- .bifrost_search_progress_session
  testthat::local_mocked_bindings(
    .bifrost_search_progress_session = function(enabled) {
      progress_session(enabled, .recording_search_renderer(events))
    }
  )

  suppressWarnings(
    searchOptimalConfiguration(
      baseline_tree = tree,
      trait_data = traits,
      min_descendant_tips = 4,
      num_cores = 1,
      shift_acceptance_threshold = 1e9,
      uncertaintyweights = TRUE,
      uncertaintyweights_par = FALSE,
      plot = FALSE,
      IC = "GIC",
      store_model_fit_history = FALSE,
      verbose = FALSE,
      progress = TRUE,
      method = "LL"
    )
  )

  statuses <- vapply(events$updates, `[[`, character(1), "status")
  testthat::expect_true(any(grepl(
    "[3/3] IC-weight re-estimation - skipped: no accepted shifts",
    statuses,
    fixed = TRUE
  )))
  testthat::expect_identical(events$done, paste0("row-", 1:3))
})

test_that("zero-candidate searches leave all three skipped stage lines", {
  testthat::skip_if_not_installed("ape")
  testthat::skip_if_not_installed("phytools")
  testthat::skip_if_not_installed("mvMORPH")

  set.seed(20260713)
  tree <- ape::rtree(10)
  traits <- matrix(stats::rnorm(20), ncol = 2)
  rownames(traits) <- tree$tip.label
  events <- new.env(parent = emptyenv())
  events$created <- character()
  events$updates <- list()
  events$output <- character()
  events$done <- character()
  progress_session <- .bifrost_search_progress_session
  testthat::local_mocked_bindings(
    .bifrost_search_progress_session = function(enabled) {
      progress_session(enabled, .recording_search_renderer(events))
    }
  )

  suppressWarnings(
    searchOptimalConfiguration(
      baseline_tree = tree,
      trait_data = traits,
      min_descendant_tips = ape::Ntip(tree),
      num_cores = 1,
      shift_acceptance_threshold = 1e9,
      uncertaintyweights = FALSE,
      uncertaintyweights_par = FALSE,
      plot = FALSE,
      IC = "GIC",
      store_model_fit_history = FALSE,
      verbose = FALSE,
      progress = TRUE,
      method = "LL"
    )
  )

  testthat::expect_identical(events$created, c(
    "[1/3] Candidate scoring",
    "[2/3] Greedy shift search",
    "[3/3] IC-weight re-estimation"
  ))
  statuses <- vapply(events$updates, `[[`, character(1), "status")
  testthat::expect_true(any(grepl("[1/3] Candidate scoring - skipped", statuses, fixed = TRUE)))
  testthat::expect_true(any(grepl("[2/3] Greedy shift search - skipped", statuses, fixed = TRUE)))
  testthat::expect_true(any(grepl("[3/3] IC-weight re-estimation - skipped", statuses, fixed = TRUE)))
  testthat::expect_identical(events$done, paste0("row-", 1:3))
})
