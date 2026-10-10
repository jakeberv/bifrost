legacy_start_cases <- function() {
  readRDS(test_path("fixtures", "mvgls-starts-1.2.1.rds"))$cases
}

legacy_case_args <- function(case) {
  args <- case$args
  args$formula <- stats::as.formula(args$formula)
  args
}

test_that("default starts reproduce official 1.2.1 across BM/BMM settings", {
  for (name in names(legacy_start_cases())) {
    case <- legacy_start_cases()[[name]]
    set.seed(71)
    before <- .Random.seed
    fit <- suppressWarnings(do.call(bifrost:::.bifrost_mvgls, legacy_case_args(case)))
    expect_equal(as.numeric(fit$start_values), case$start, tolerance = 1e-12, info = name)
    expect_identical(.Random.seed, before, info = name)
    expect_identical(attr(fit, "bifrost_initialization")$used, "legacy_1.2.1")
  }
})

test_that("native and supplied starts reach mvMORPH unchanged", {
  args <- legacy_case_args(legacy_start_cases()$singleton_hl)
  for (strategy in c("native", "legacy_1.2.1")) {
    supplied <- args
    if (strategy == "legacy_1.2.1") supplied$start <- c(.5, .2, .3, .01)
    set.seed(12)
    ref <- suppressWarnings(do.call(mvMORPH::mvgls, supplied))
    rng <- .Random.seed
    set.seed(12)
    fit <- suppressWarnings(do.call(bifrost:::.bifrost_mvgls,
                                    c(supplied, list(start_strategy = strategy))))
    expect_identical(fit$start_values, ref$start_values)
    expect_identical(fit$opt, ref$opt)
    expect_identical(.Random.seed, rng)
    expect_identical(attr(fit, "bifrost_initialization")$used,
                     if (strategy == "native") "native" else "supplied")
  }
})

test_that("starts are recalculated from each candidate's data", {
  args <- legacy_case_args(legacy_start_cases()$tipless_hl)
  first <- do.call(bifrost:::.bifrost_mvgls, args)
  args$data$Y <- args$data$Y * 10
  second <- do.call(bifrost:::.bifrost_mvgls, args)
  # Rate starts are square roots of variances; scaling Y tenfold scales them tenfold.
  expect_equal(as.numeric(second$start_values)[2:3],
               as.numeric(first$start_values)[2:3] * 10)
})

test_that("explicitly disabling the initialization grid is respected", {
  args <- legacy_case_args(legacy_start_cases()$bm_hl)
  args$grid.search <- FALSE
  ref <- do.call(mvMORPH::mvgls, args)
  fit <- do.call(bifrost:::.bifrost_mvgls, args)
  expect_identical(fit$start_values, ref$start_values)
  expect_identical(fit$opt, ref$opt)
  expect_identical(attr(fit, "bifrost_initialization")$used, "native_no_grid")
})

test_that("unknown strategies and unsupported legacy settings fail clearly", {
  args <- legacy_case_args(legacy_start_cases()$bm_hl)
  expect_error(do.call(bifrost:::.bifrost_mvgls,
                       c(args, list(start_strategy = "typo"))), "arg")
  args$method <- "EmpBayes"
  expect_error(do.call(bifrost:::.bifrost_mvgls, args), "start_strategy.*native")
})

test_that("every search stage uses the selected starting policy", {
  args <- legacy_case_args(legacy_start_cases()$bm_hl)
  Y <- args$data$Y[, 1:2]
  Y[1:8, ] <- Y[1:8, ] * 20
  ns <- asNamespace("bifrost")
  original <- get(".bifrost_mvgls", ns)
  seen <- list()
  local_rebind(".bifrost_mvgls", function(...) {
    fit <- original(...)
    seen[[length(seen) + 1L]] <<- attr(fit, "bifrost_initialization")
    fit
  }, ns)
  for (strategy in c("legacy_1.2.1", "native")) {
    seen <- list()
    result <- suppressWarnings(searchOptimalConfiguration(
      args$tree, Y, min_descendant_tips = 4, num_cores = 1,
      shift_acceptance_threshold = 0, uncertaintyweights = TRUE,
      method = "LL", progress = FALSE, start_strategy = strategy
    ))
    expect_gt(length(result$shift_nodes_no_uncertainty), 0L)
    expect_length(seen, 1L + result$num_candidates +
                    length(result$model_fit_history$fits) +
                    length(result$shift_nodes_no_uncertainty))
    expect_true(all(vapply(seen, function(x) identical(x$used, strategy), logical(1))))
    expect_identical(result$user_input$start_strategy, strategy)
    expect_identical(result$initialization,
                     attr(result$model_no_uncertainty, "bifrost_initialization"))
    expect_identical(result$initialization$mvMORPH_version,
                     as.character(utils::packageVersion("mvMORPH")))
  }
})

test_that("legacy preparation respects row alignment and a response override", {
  args <- legacy_case_args(legacy_start_cases()$formula_ll)
  ref <- do.call(bifrost:::.bifrost_mvgls, args)
  args$data <- args$data[rev(seq_len(nrow(args$data))), ]
  shuffled <- do.call(bifrost:::.bifrost_mvgls, args)
  expect_identical(shuffled$start_values, ref$start_values)
  expect_identical(shuffled$opt, ref$opt)
  args$response <- as.matrix(args$data[, c("y1", "y2")]) * 2
  overridden <- do.call(bifrost:::.bifrost_mvgls, args)
  args$response <- NULL
  args$data[, c("y1", "y2")] <- args$data[, c("y1", "y2")] * 2
  replaced <- do.call(bifrost:::.bifrost_mvgls, args)
  expect_identical(overridden$start_values, replaced$start_values)
  expect_identical(overridden$opt, replaced$opt)
})

test_that("native and explicit controls bypass the legacy backend entirely", {
  local_rebind(".bifrost_mvgls_start_backend", function() {
    stop("Legacy backend unavailable")
  }, asNamespace("bifrost"))
  args <- legacy_case_args(legacy_start_cases()$bm_hl)
  expect_error(do.call(bifrost:::.bifrost_mvgls, args), "Legacy backend unavailable")
  for (control in list(list(start_strategy = "native"),
                       list(start = legacy_start_cases()$bm_hl$start),
                       list(grid.search = FALSE), list(grid.search = 0))) {
    expect_s3_class(do.call(bifrost:::.bifrost_mvgls, c(args, control)), "mvgls")
  }
})

test_that("numeric logical controls retain mvgls coercion semantics", {
  case <- legacy_start_cases()$scaled_height
  args <- legacy_case_args(case)
  args$scale.height <- 1
  fit <- do.call(bifrost:::.bifrost_mvgls, args)
  expect_equal(as.numeric(fit$start_values), case$start, tolerance = 1e-12)
})
