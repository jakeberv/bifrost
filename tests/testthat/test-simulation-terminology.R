test_that("manuscript-generator calibration is isolated as an experimental draft", {
  part1 <- testthat::test_path(
    "../../vignettes/simulation-study-part-1.Rmd"
  )
  draft <- testthat::test_path(
    "../../experimental/simulation-study-manuscript-generator-calibration.Rmd"
  )
  build_ignore <- testthat::test_path("../../.Rbuildignore")
  testthat::skip_if_not(
    file.exists(part1),
    "Part 1 vignette source is unavailable in installed tests"
  )
  testthat::skip_if_not(
    file.exists(draft),
    "Experimental calibration draft is unavailable in installed tests"
  )

  part1_text <- paste(readLines(part1, warn = FALSE), collapse = "\n")
  testthat::expect_false(grepl(
    "Optional: Empirical Calibration of the Manuscript Generator",
    part1_text,
    fixed = TRUE
  ))
  testthat::expect_false(grepl(
    "original_generator_calibration",
    part1_text,
    fixed = TRUE
  ))

  draft_text <- paste(readLines(draft, warn = FALSE), collapse = "\n")
  testthat::expect_match(
    draft_text,
    "Draft: Calibrating the Original Manuscript Generator",
    fixed = TRUE
  )
  testthat::expect_match(
    draft_text,
    'simulation_generator = "original"',
    fixed = TRUE
  )
  testthat::expect_false(grepl("VignetteIndexEntry", draft_text, fixed = TRUE))

  testthat::expect_true(file.exists(build_ignore))
  testthat::expect_true(any(readLines(build_ignore, warn = FALSE) == "^experimental$"))
})

test_that("simulation vignettes retain candidate and evaluability safeguards", {
  part1_path <- testthat::test_path(
    "../../vignettes/simulation-study-part-1.Rmd"
  )
  part2_path <- testthat::test_path(
    "../../vignettes/simulation-study-part-2.Rmd"
  )
  testthat::skip_if_not(
    all(file.exists(c(part1_path, part2_path))),
    "Simulation vignette sources are unavailable in installed tests"
  )

  part1 <- paste(readLines(part1_path, warn = FALSE), collapse = "\n")
  part2 <- paste(readLines(part2_path, warn = FALSE), collapse = "\n")

  testthat::expect_match(part1, "per_replicate$n_candidates", fixed = TRUE)
  testthat::expect_match(part2, "min_evaluable_fraction", fixed = TRUE)
})
