.load_comparison_script <- function() {
  script_path <- testthat::test_path(
    "..", "..", "tools", "compare-simulation-generators.R"
  )
  testthat::skip_if_not(
    file.exists(script_path),
    "Development-only comparison script is excluded from source tarballs"
  )

  script_env <- new.env(parent = globalenv())
  sys.source(script_path, envir = script_env)
  script_env
}

test_that("comparison script parses reproducible quick-run controls", {
  script_env <- .load_comparison_script()
  options <- script_env$.comparison_parse_args(c(
    "--replicates", "3",
    "--output-dir", "/tmp/comparison-output",
    "--workers", "2",
    "--diagnostic-draws", "25"
  ))

  testthat::expect_identical(options$replicates, 3L)
  testthat::expect_identical(options$output_dir, "/tmp/comparison-output")
  testthat::expect_identical(options$workers, 2L)
  testthat::expect_identical(options$diagnostic_draws, 25L)
  testthat::expect_identical(options$tree_tip_count, 200L)
  testthat::expect_identical(options$num_shifts, 5L)
})

test_that("comparison script rejects incomplete and invalid controls", {
  script_env <- .load_comparison_script()

  testthat::expect_error(
    script_env$.comparison_parse_args(c("--replicates", "0")),
    "replicates"
  )
  testthat::expect_error(
    script_env$.comparison_parse_args(c("--output-dir")),
    "requires a value"
  )
  testthat::expect_error(
    script_env$.comparison_parse_args(c("--unknown", "1")),
    "Unknown option"
  )
})

test_that("comparison script formats an empty warning log as character output", {
  script_env <- .load_comparison_script()

  testthat::expect_identical(
    script_env$.comparison_warning_lines(list()),
    character(0)
  )
  testthat::expect_identical(
    script_env$.comparison_warning_lines(list(a = c("first", "second"))),
    c("[a] first", "[a] second")
  )
})

test_that("comparison script uses public simulation generator vocabulary", {
  script_env <- .load_comparison_script()

  testthat::expect_identical(
    script_env$.comparison_generators,
    c("original", "empirical")
  )
  error_row <- script_env$.comparison_per_replicate(
    structure(list(error = "expected"), class = "comparison_error"),
    "empirical",
    "null"
  )
  testthat::expect_true("simulation_generator" %in% names(error_row))
  testthat::expect_false("simulation_design" %in% names(error_row))
})

test_that("comparison script excludes failed searches from inferred-shift means", {
  script_env <- .load_comparison_script()

  study <- list(
    simdata = list(
      list(covariance_matrix = diag(2L)),
      list(covariance_matrix = diag(2L))
    ),
    results = list(
      list(error = "fit failed", shift_nodes_no_uncertainty = seq_len(9L)),
      list(error = NULL, shift_nodes_no_uncertainty = seq_len(2L))
    )
  )
  rows <- script_env$.comparison_per_replicate(study, "original", "null")
  summary <- script_env$.comparison_summarize(rows)

  testthat::expect_identical(rows$status, c("error", "ok"))
  testthat::expect_true(is.na(rows$n_inferred_shifts[[1L]]))
  testthat::expect_identical(rows$n_inferred_shifts[[2L]], 2L)
  testthat::expect_identical(summary$mean_inferred_shifts, 2)
})
