.collect_progress_handler <- function(events) {
  record <- function(config, state, progression, ...) {
    events$types <- c(events$types, progression$type)
    events$messages <- c(events$messages, paste(state$message, collapse = ""))
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
  testthat::expect_identical(events$types, c("initiate", "update", "update", "finish"))
  testthat::expect_true(any(grepl("Candidate scoring complete", events$messages, fixed = TRUE)))
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

test_that("search progress handler is a persistent cli handler", {
  handler <- .bifrost_search_progress_handler()
  testthat::expect_s3_class(handler, "cli_progression_handler")
  testthat::expect_false(get("clear", envir = environment(handler)))
  testthat::expect_identical(get("enable_after", envir = environment(handler)), 0)
  testthat::expect_true(get("enable", envir = environment(handler)))
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

test_that("search emits candidate and greedy stages independently of verbose", {
  testthat::skip_if_not_installed("ape")
  testthat::skip_if_not_installed("phytools")
  testthat::skip_if_not_installed("mvMORPH")

  events <- new.env(parent = emptyenv())
  events$types <- character()
  events$messages <- character()
  testthat::local_mocked_bindings(
    .bifrost_search_progress_handler = function() .collect_progress_handler(events)
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
      verbose = FALSE,
      method = "LL"
    )
  )

  finish_messages <- events$messages[events$types == "finish"]
  testthat::expect_s3_class(result, "bifrost_search")
  testthat::expect_identical(result$user_input$progress, TRUE)
  testthat::expect_true(any(grepl("[1/3] Candidate scoring", finish_messages, fixed = TRUE)))
  testthat::expect_true(any(grepl("[2/3] Greedy shift search", finish_messages, fixed = TRUE)))
  testthat::expect_match(skipped, "[3/3] IC-weight re-estimation", fixed = TRUE)
  testthat::expect_false(any(grepl("Generating candidate shift models", skipped, fixed = TRUE)))
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
    list(GIC = list(GIC = 90))
  }

  run_weights <- function(parallel) {
    events <- new.env(parent = emptyenv())
    events$types <- character()
    events$messages <- character()
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
      messages = events$messages[events$types == "update"]
    )
  }

  serial <- run_weights(FALSE)
  parallel <- run_weights(TRUE)

  expected <- paste("[3/3] IC-weight re-estimation - completed node", shift_nodes)
  testthat::expect_equal(nrow(serial$value), 2L)
  testthat::expect_equal(nrow(parallel$value), 2L)
  testthat::expect_setequal(serial$messages, expected)
  testthat::expect_setequal(parallel$messages, expected)
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
  events$types <- character()
  events$messages <- character()
  testthat::local_mocked_bindings(
    .bifrost_search_progress_handler = function() .collect_progress_handler(events)
  )

  messages <- testthat::capture_messages(
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

  testthat::expect_true(any(grepl(
    "[3/3] IC-weight re-estimation - skipped: no accepted shifts",
    messages,
    fixed = TRUE
  )))
})

test_that("zero-candidate searches leave all three skipped stage lines", {
  testthat::skip_if_not_installed("ape")
  testthat::skip_if_not_installed("phytools")
  testthat::skip_if_not_installed("mvMORPH")

  set.seed(20260713)
  tree <- ape::rtree(10)
  traits <- matrix(stats::rnorm(20), ncol = 2)
  rownames(traits) <- tree$tip.label

  messages <- testthat::capture_messages(
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

  testthat::expect_true(any(grepl("[1/3] Candidate scoring - skipped", messages, fixed = TRUE)))
  testthat::expect_true(any(grepl("[2/3] Greedy shift search - skipped", messages, fixed = TRUE)))
  testthat::expect_true(any(grepl("[3/3] IC-weight re-estimation - skipped", messages, fixed = TRUE)))
})
