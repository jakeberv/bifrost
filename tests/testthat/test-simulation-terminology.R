test_that("simulation documentation uses empirically calibrated terminology", {
  paths <- c(
    testthat::test_path("../../NEWS.md"),
    testthat::test_path("../../R/simulation-template.R"),
    testthat::test_path("../../R/simulation-generators.R"),
    testthat::test_path("../../R/simulation-studies.R")
  )
  vignette_paths <- testthat::test_path(
    "../../vignettes",
    c("simulation-study-part-1.Rmd", "simulation-study-part-2.Rmd")
  )
  paths <- c(paths, vignette_paths)
  paths <- paths[file.exists(paths)]
  testthat::skip_if(
    length(paths) == 0L,
    "Simulation documentation sources are unavailable in installed tests"
  )
  text <- unlist(lapply(paths, readLines, warn = FALSE), use.names = FALSE)

  testthat::expect_false(any(grepl(
    "dataset[- ]matched",
    text,
    ignore.case = TRUE
  )))
  testthat::expect_true(any(grepl(
    "empirically calibrated",
    text,
    ignore.case = TRUE
  )))
})

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

  testthat::expect_true(file.exists(draft))
  if (file.exists(draft)) {
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
  }

  testthat::expect_true(file.exists(build_ignore))
  testthat::expect_true(any(readLines(build_ignore, warn = FALSE) == "^experimental$"))
})

reader_visible_r_chunks <- function(path) {
  text <- readLines(path, warn = FALSE)
  starts <- grep("^```\\{r(?:[ ,}])", text, perl = TRUE)

  chunks <- lapply(starts, function(start) {
    end <- which(text == "```" & seq_along(text) > start)[1L]
    header <- text[[start]]
    body <- text[seq.int(start + 1L, end - 1L)]
    body <- body[nzchar(trimws(body))]

    list(header = header, body = body)
  })

  chunks[!vapply(chunks, function(chunk) {
    grepl("(?:echo|include)\\s*=\\s*FALSE", chunk$header, perl = TRUE)
  }, logical(1L))]
}

test_that("reader-visible simulation chunks open with intent comments", {
  paths <- testthat::test_path(
    "../../vignettes",
    c("simulation-study-part-1.Rmd", "simulation-study-part-2.Rmd")
  )
  testthat::skip_if_not(
    all(file.exists(paths)),
    "Simulation vignette sources are unavailable in installed tests"
  )

  visible_chunks <- lapply(paths, reader_visible_r_chunks)
  testthat::expect_identical(lengths(visible_chunks), c(10L, 11L))

  for (chunks in visible_chunks) {
    missing_comments <- vapply(chunks, function(chunk) {
      length(chunk$body) == 0L || !grepl("^# ", chunk$body[[1L]])
    }, logical(1L))
    testthat::expect_equal(
      vapply(chunks[missing_comments], `[[`, character(1L), "header"),
      character()
    )
  }
})

test_that("Part 1 tables omit Evaluable without removing its safeguards", {
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

  testthat::expect_false(grepl('"Evaluable"', part1, fixed = TRUE))
  testthat::expect_match(part1, "unevaluable", fixed = TRUE)
  testthat::expect_match(part2, "min_evaluable_fraction", fixed = TRUE)
})
