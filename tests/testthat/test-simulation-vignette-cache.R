test_that("simulation vignette cache records the empirical generator", {
  cache_path <- testthat::test_path(
    "../../inst/extdata/simulation-study-cache/passerine_preview_tables.rds"
  )
  testthat::skip_if_not(file.exists(cache_path), "Simulation vignette cache not present")

  cache <- readRDS(cache_path)

  testthat::expect_identical(cache$schema_version, 2L)
  testthat::expect_identical(cache$provenance$simulation_generator, "empirical")
  testthat::expect_false("simulation_design" %in% names(cache$provenance))
  testthat::expect_identical(cache$provenance$n_replicates_per_setting, 100L)
  testthat::expect_identical(cache$provenance$tree_tip_count, 250L)
  testthat::expect_identical(
    cache$provenance$integration_power_range,
    c(0.5, 1.25)
  )
  testthat::expect_identical(
    cache$provenance$integration_exclude_range,
    c(0.8, 1.1)
  )

  testthat::expect_s3_class(cache$fixed_settings, "data.frame")
  testthat::expect_equal(nrow(cache$fixed_settings), 6L)
  testthat::expect_setequal(
    cache$fixed_settings$Scenario,
    c("Null", "Proportional", "Integration-rate")
  )
  testthat::expect_setequal(names(cache$grid_summary), c("gic", "bic"))
  testthat::expect_true(all(vapply(cache$grid_summary, nrow, integer(1L)) == 6L))

  serialized_text <- unlist(cache, recursive = TRUE, use.names = TRUE)
  serialized_text <- as.character(serialized_text)
  machine_specific <- grepl(
    "(^|[= ])(/Users/|/home/|/private/tmp/|[A-Za-z]:[/\\\\])",
    serialized_text
  )
  personal_email <- grepl(
    "[[:alnum:]._%+-]+@[[:alnum:].-]+\\.[[:alpha:]]{2,}",
    serialized_text,
    ignore.case = TRUE
  )
  testthat::expect_false(any(machine_specific))
  testthat::expect_false(any(personal_email))
})

test_that("empirical benchmark code pins every scenario explicitly", {
  part1_path <- testthat::test_path("../../vignettes/simulation-study-part-1.Rmd")
  part2_path <- testthat::test_path("../../vignettes/simulation-study-part-2.Rmd")
  testthat::skip_if_not(
    all(file.exists(c(part1_path, part2_path))),
    "Simulation vignette not present"
  )

  part1 <- paste(readLines(
    part1_path,
    warn = FALSE
  ), collapse = "\n")
  part2 <- paste(readLines(
    testthat::test_path("../../vignettes/simulation-study-part-2.Rmd"),
    warn = FALSE
  ), collapse = "\n")
  combined <- paste(part1, part2)

  testthat::expect_match(
    part1,
    'simulation_options = list\\(\\n    simulation_generator = "empirical"'
  )
  testthat::expect_match(
    part1,
    'null_simulation_options = list\\(\\n    simulation_generator = "empirical"'
  )
  occurrences <- gregexpr(
    'simulation_generator = "empirical"',
    combined,
    fixed = TRUE
  )[[1L]]
  testthat::expect_gte(sum(occurrences > 0L), 8L)
})

test_that("Part 1 preserves the original correlation color scale", {
  vignette_path <- testthat::test_path(
    "../../vignettes/simulation-study-part-1.Rmd"
  )
  testthat::skip_if_not(file.exists(vignette_path), "Simulation vignette not present")

  vignette_source <- readLines(vignette_path, warn = FALSE)

  testthat::expect_true(any(grepl(
    "limits = c(0, 0.75)",
    vignette_source,
    fixed = TRUE
  )))
  testthat::expect_false(any(grepl(
    "limits = c(-1, 1)",
    vignette_source,
    fixed = TRUE
  )))
  testthat::expect_false(any(grepl(
    "limits = c(0, 1)",
    vignette_source,
    fixed = TRUE
  )))
})

test_that("Part 1 demonstrates calibrated reduced integration with the original generator", {
  vignette_path <- testthat::test_path(
    "../../vignettes/simulation-study-part-1.Rmd"
  )
  testthat::skip_if_not(file.exists(vignette_path), "Simulation vignette not present")

  source <- paste(readLines(
    vignette_path,
    warn = FALSE
  ), collapse = "\n")

  testthat::expect_match(source, "Calibrating the manuscript generator")
  testthat::expect_match(source, 'simulation_generator = "original"', fixed = TRUE)
  testthat::expect_match(source, "original_calibration_grid", fixed = TRUE)
  testthat::expect_match(
    source,
    "scale_factor_range = c(0.35, 0.75)",
    fixed = TRUE
  )
  testthat::expect_match(source, "original_reduced_diagnostics", fixed = TRUE)
  testthat::expect_match(
    source,
    "does not preserve the empirical pairwise",
    fixed = TRUE
  )
})
