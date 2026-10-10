test_that("progress preserves serial and parallel candidate random draws", {
  withr::local_seed(51)
  backends <- if (future::supportsMulticore()) c(TRUE, FALSE) else TRUE
  draw <- function(i) list(sample = sample.int(100L, i), normal = rnorm(i))
  for (kind in c("Mersenne-Twister", "L'Ecuyer-CMRG")) {
    for (backend in backends) {
      for (cores in 1:2) {
        for (jobs in c(1L, 3L)) {
          items <- setNames(seq_len(jobs), paste0("candidate-", seq_len(jobs)))
          set.seed(51, kind = kind)
          plain <- .bifrost_search_lapply(items, draw, cores, backend)
          expected_seed <- .Random.seed
          expected_kind <- RNGkind()
          set.seed(51, kind = kind)
          displayed <- .bifrost_search_lapply(
            items, draw, cores, backend,
            heartbeat = function() invisible(NULL)
          )
          info <- paste(kind, backend, cores, jobs)
          expect_identical(displayed, plain, info = info)
          expect_identical(.Random.seed, expected_seed, info = info)
          expect_identical(RNGkind(), expected_kind, info = info)
        }
      }
    }
  }
})

test_that("background serial fits return the consumed RNG state even on failure", {
  withr::local_seed(61)
  for (fail in c(FALSE, TRUE)) {
    work <- function() {
      value <- sample.int(1000L, 7L)
      if (fail) stop("fit failed after drawing")
      value
    }
    set.seed(61)
    plain <- tryCatch(work(), error = identity)
    expected_seed <- .Random.seed
    set.seed(61)
    displayed <- tryCatch(.bifrost_search_with_future_plan(
      1L, TRUE, ensure_async = TRUE,
      work = function() .bifrost_search_await_work(work)
    ), error = identity)
    if (fail) {
      expect_identical(class(displayed), class(plain))
      expect_identical(conditionMessage(displayed), conditionMessage(plain))
    } else expect_identical(displayed, plain)
    expect_identical(.Random.seed, expected_seed)
  }
})

test_that("a background fit preserves an absent seed when it draws no randomness", {
  withr::local_seed(71)
  rm('.Random.seed', envir = globalenv())
  kind <- RNGkind()
  result <- .bifrost_search_with_future_plan(
    1L, TRUE, ensure_async = TRUE,
    work = function() .bifrost_search_await_work(function() 12L)
  )
  expect_identical(result, 12L)
  expect_false(exists('.Random.seed', envir = globalenv(), inherits = FALSE))
  expect_identical(RNGkind(), kind)
})

test_that("progress preserves a real singleton BMM fit and its next random draw", {
  # Multisession workers can load the current package during R CMD check, but
  # would load an older installation in a development load_all() session.
  installed <- file.exists(file.path(getNamespaceInfo("bifrost", "path"),
                                      "Meta", "package.rds"))
  backends <- c(if (future::supportsMulticore()) FALSE, if (installed) TRUE)
  skip_if(length(backends) == 0L, "current package unavailable to workers")
  case <- readRDS(test_path("fixtures", "mvgls-starts-historical.rds"))$cases$singleton_hl
  args <- case$args
  args$formula <- as.formula(args$formula)
  fit <- function() suppressWarnings(do.call(bifrost:::.bifrost_mvgls, args))
  withr::local_seed(71)
  for (backend in backends) {
    set.seed(71)
    plain <- fit()
    expected_next <- runif(1L)
    set.seed(71)
    displayed <- .bifrost_search_with_future_plan(
      1L, backend, ensure_async = TRUE,
      work = function() .bifrost_search_await_work(fit)
    )
    expect_identical(displayed$start_values, plain$start_values)
    expect_identical(displayed$opt, plain$opt)
    expect_identical(runif(1L), expected_next)
  }
})
