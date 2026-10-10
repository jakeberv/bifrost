historical_start_cases <- function() {
  readRDS(test_path("fixtures", "mvgls-starts-historical.rds"))$cases
}

historical_case_args <- function(case) {
  args <- case$args
  args$formula <- stats::as.formula(args$formula)
  args
}

test_that("default starts reproduce historical starts and random draws across BM/BMM settings", {
  rng <- bifrost:::.bifrost_search_rng_state()
  withr::defer(bifrost:::.bifrost_search_restore_rng(rng))
  # Match the fixture's RNG configuration, including R-devel's new binomial kind.
  # R-devel warns about the historical binomial algorithm; these fits do not use it.
  suppressWarnings(RNGversion("4.4.2"))
  for (name in names(historical_start_cases())) {
    case <- historical_start_cases()[[name]]
    set.seed(71)
    expect_identical(.Random.seed, case$rng_before, info = name)
    fit <- suppressWarnings(do.call(bifrost:::.bifrost_mvgls, historical_case_args(case)))
    expect_equal(as.numeric(fit$start_values), case$start, tolerance = 1e-12, info = name)
    expect_identical(.Random.seed, case$rng_after, info = name)
    expect_identical(attr(fit, "bifrost_initialization")$used, "legacy_1.2.1")
  }
})

test_that("native and supplied starts reach mvMORPH unchanged", {
  args <- historical_case_args(historical_start_cases()$singleton_hl)
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
  args <- historical_case_args(historical_start_cases()$tipless_hl)
  first <- do.call(bifrost:::.bifrost_mvgls, args)
  args$data$Y <- args$data$Y * 10
  second <- do.call(bifrost:::.bifrost_mvgls, args)
  # Rate starts are square roots of variances; scaling Y tenfold scales them tenfold.
  expect_equal(as.numeric(second$start_values)[2:3],
               as.numeric(first$start_values)[2:3] * 10)
})

test_that("explicitly disabling the initialization grid is respected", {
  args <- historical_case_args(historical_start_cases()$bm_hl)
  args$grid.search <- FALSE
  ref <- do.call(mvMORPH::mvgls, args)
  fit <- do.call(bifrost:::.bifrost_mvgls, args)
  expect_identical(fit$start_values, ref$start_values)
  expect_identical(fit$opt, ref$opt)
  expect_identical(attr(fit, "bifrost_initialization")$used, "native_no_grid")
})

test_that("unknown strategies and unsupported historical settings fail clearly", {
  args <- historical_case_args(historical_start_cases()$bm_hl)
  expect_error(do.call(bifrost:::.bifrost_mvgls,
                       c(args, list(start_strategy = "typo"))), "arg")
  args$method <- "EmpBayes"
  expect_error(do.call(bifrost:::.bifrost_mvgls, args), "start_strategy.*native")
})

test_that("every search stage uses the selected starting policy", {
  args <- historical_case_args(historical_start_cases()$bm_hl)
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

test_that("historical preparation respects row alignment and a response override", {
  args <- historical_case_args(historical_start_cases()$formula_ll)
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

test_that("native and explicit controls bypass the historical backend entirely", {
  local_rebind(".bifrost_mvgls_start_backend", function() {
    stop("Historical backend unavailable")
  }, asNamespace("bifrost"))
  args <- historical_case_args(historical_start_cases()$bm_hl)
  expect_error(do.call(bifrost:::.bifrost_mvgls, args), "Historical backend unavailable")
  for (control in list(list(start_strategy = "native"),
                       list(start = historical_start_cases()$bm_hl$start),
                       list(grid.search = FALSE), list(grid.search = 0))) {
    expect_s3_class(do.call(bifrost:::.bifrost_mvgls, c(args, control)), "mvgls")
  }
})

test_that("numeric logical controls retain mvgls coercion semantics", {
  case <- historical_start_cases()$scaled_height
  args <- historical_case_args(case)
  args$scale.height <- 1
  fit <- do.call(bifrost:::.bifrost_mvgls, args)
  expect_equal(as.numeric(fit$start_values), case$start, tolerance = 1e-12)
})

test_that("historical initialization rejects incompatible trees, data, and methods", {
  args <- historical_case_args(historical_start_cases()$bm_hl)
  invalid <- list(
    list(change = list(tree = list()), error = "phylo tree"),
    list(change = list(model = "BMM", tree = structure(args$tree, class = "phylo")), error = "simmap"),
    list(change = list(penalty = "RidgeAlt"), error = "require RidgeArch"),
    list(change = list(data = list(Y = args$data$Y[-1, ])), error = "one row per tip"),
    list(change = list(data = list(Y = args$data$Y[, 1, drop = FALSE])),
         error = "multivariate datasets"),
    list(change = list(method = "LL"), error = "more variables than observations")
  )
  for (case in invalid) {
    invalid_args <- args
    invalid_args[names(case$change)] <- case$change
    expect_error(do.call(bifrost:::.bifrost_mvgls, invalid_args), case$error)
  }
})

test_that("a changed mvMORPH private interface fails with the native escape hatch", {
  args <- historical_case_args(historical_start_cases()$bm_hl)
  # Simulate an upstream signature change in memory, never in the installation.
  local_rebind(".setBounds", function(penalty) NULL, asNamespace("mvMORPH"))
  expect_error(do.call(bifrost:::.bifrost_mvgls, args), "incompatible.*native")
})

test_that("historical grid selection respects tolerances and rejects unusable scores", {
  corr <- list(model = "BM", nobs = 16L)
  grid_start <- function(penalty, tol, score) {
    bifrost:::.bifrost_mvgls_start_grid_historical(
      corr, "LOOCV", penalty, NULL, tol, score
    )
  }
  score <- function(par, ...) par[1L]^2
  # RidgeAlt candidates are log-transformed: tol = 1 leaves 10 as the minimum.
  expect_equal(grid_start("RidgeAlt", 1, score), c(log(10), 1))
  expect_error(grid_start("RidgeArch", 1, score), "excludes every.*grid point")
  expect_error(grid_start("RidgeArch", NULL, function(...) NA_real_),
               "No usable.*starting values")
})
