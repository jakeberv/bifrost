paired_tuning_paths <- function() {
  root <- testthat::test_path("../..")
  list(
    exporter = file.path(root, "data-raw", "update_simulation_tuning_cache.R"),
    inputs = file.path(root, "local-cache", "paired-tuning-10-40-export-inputs"),
    baseline = file.path(
      root, "local-cache", "paired-tuning-10-40-export-inputs",
      "pre-part1-refresh-schema4.rds"
    ),
    artifact = file.path(
      root, "data-remote", "simulation-study-cache",
      "passerine_preview_tables.rds"
    )
  )
}

load_paired_tuning_helpers <- function() {
  paths <- paired_tuning_paths()
  testthat::skip_if_not(file.exists(paths$exporter), "paired-tuning exporter absent")
  env <- new.env(parent = globalenv())
  sys.source(paths$exporter, envir = env)
  env
}

write_paired_tuning_fixture <- function(cache, directory) {
  accounting <- cache$provenance$tuning$generation_accounting
  provenance <- cache$provenance$tuning[c(
    "package_version", "package_remote_sha", "package_github_sha1", "data_md5",
    "source_fingerprint", "runner_sha256", "template_traits", "template_tips",
    "template_covariance_sha256", "template_residual_df", "platform", "R",
    "runner_files_sha256"
  )]
  for (key in c("gic", "bic")) {
    saveRDS(list(
      schema_version = 1L,
      grid = cache$tuning_grids[[key]],
      policy = cache$provenance$tuning$policy,
      recommendation = list(),
      provenance = provenance,
      summary_audit = list()
    ), file.path(directory, paste0(key, "-tuning.rds")))
  }
  metric_suffixes <- c(
    "fuzzy_precision", "fuzzy_recall", "fuzzy_f1", "fuzzy_specificity",
    "fuzzy_fpr", "fuzzy_balanced_accuracy", "strict_precision",
    "strict_recall", "strict_f1", "strict_specificity", "strict_fpr",
    "strict_balanced_accuracy", "weighted_fuzzy_f1"
  )
  summaries <- do.call(rbind, lapply(cache$tuning_grids, function(grid) {
    summary <- grid$summary_table
    do.call(rbind, lapply(c("null", "proportional", "integration-rate"),
      function(scenario) {
        prefix <- if (scenario == "integration-rate") "correlation" else scenario
        out <- data.frame(
          config_id = paste(scenario, grid$IC, summary$setting_id, sep = "-"),
          scenario = scenario,
          setting_id = summary$setting_id,
          IC = grid$IC,
          shift_acceptance_threshold = summary$shift_acceptance_threshold,
          min_descendant_tips = summary$min_descendant_tips,
          study_seed = summary[[paste0(prefix, "_seed")]],
          completion_rate = summary[[paste0(prefix, "_completion_rate")]],
          mean_inferred_shifts = summary[[paste0(prefix, "_mean_inferred_shifts")]],
          evaluable_fraction = summary[[paste0(prefix, "_evaluable_fraction")]],
          n_attempted = summary[[paste0(prefix, "_n_attempted")]],
          n_generated = summary[[paste0(prefix, "_n_generated")]],
          n_placement_failed = summary[[paste0(prefix, "_n_placement_failed")]],
          n_search_attempted = summary[[paste0(prefix, "_n_search_attempted")]],
          n_evaluated = summary[[paste0(prefix, "_n_evaluated")]],
          stringsAsFactors = FALSE
        )
        out$mean_false_positive_rate <- NA_real_
        out$fraction_any_false_positive <- NA_real_
        planted <- cache$provenance$tuning$planted_clade_counts
        planted <- planted[planted$scenario == scenario, , drop = FALSE]
        out$n_true_shifts_10_19 <- planted$n_true_shifts_10_19
        out$n_true_shifts_20_40 <- planted$n_true_shifts_20_40
        for (metric in metric_suffixes) out[[metric]] <- NA_real_
        if (scenario == "null") {
          out$mean_false_positive_rate <- summary$null_mean_false_positive_rate
          out$fraction_any_false_positive <- summary$null_fraction_any_false_positive
        } else {
          for (metric in metric_suffixes) {
            out[[metric]] <- summary[[paste0(prefix, "_", metric)]]
          }
        }
        out
      }
    ))
  }))
  attempts <- data.frame(status = rep("generated", 1500L), stringsAsFactors = FALSE)
  saveRDS(list(
    selection_policy = cache$provenance$tuning$policy,
    source_fingerprint = provenance$source_fingerprint,
    design_fingerprint = cache$provenance$tuning$paired_design$design_fingerprint,
    all_attempts_accounted_for = TRUE,
    attempts = attempts,
    group_counts = accounting,
    summary_table = summaries,
    accounting_scope = cache$provenance$tuning$accounting_scope
  ), file.path(directory, "suite-summary.rds"))
  utils::write.csv(accounting, file.path(directory, "generation-accounting.csv"),
    row.names = FALSE)
}

test_that("paired-tuning exporter is sourceable without executing", {
  paths <- paired_tuning_paths()
  testthat::skip_if_not(file.exists(paths$exporter), "paired-tuning exporter absent")

  env <- new.env(parent = globalenv())
  testthat::expect_no_error(sys.source(paths$exporter, envir = env))
  testthat::expect_true(exists("build_paired_tuning_cache", envir = env))
  testthat::expect_true(exists("write_paired_tuning_cache", envir = env))
})

test_that("paired summaries derive fixed settings without changing Part 2", {
  paths <- paired_tuning_paths()
  testthat::skip_if_not(dir.exists(paths$inputs), "paired-tuning inputs absent")
  helpers <- load_paired_tuning_helpers()
  baseline <- readRDS(paths$baseline)
  updated <- helpers$build_paired_tuning_cache(paths$inputs)

  testthat::expect_identical(updated$schema_version, 4L)
  testthat::expect_identical(updated$tuning_grids, baseline$tuning_grids)
  testthat::expect_identical(updated$tuning, baseline$tuning)
  testthat::expect_identical(updated$grid_summary, baseline$grid_summary)
  testthat::expect_identical(updated$provenance$tuning, baseline$provenance$tuning)
  testthat::expect_identical(updated$fixed_settings,
    helpers$fixed_setting_rows(updated$tuning_grids))
  testthat::expect_identical(names(updated$tuning_grids), c("gic", "bic"))
  testthat::expect_true(all(vapply(
    updated$tuning_grids,
    inherits,
    logical(1L),
    what = "bifrost_search_tuning_grid"
  )))
  testthat::expect_true(all(vapply(
    updated$tuning_grids,
    function(grid) is.null(grid$studies) && isFALSE(grid$store_studies),
    logical(1L)
  )))

  policy <- list(
    max_false_positive_rate = 0.10,
    max_any_false_positive = 0.05,
    min_evaluable_fraction = 0.50,
    primary_metric = "fuzzy_balanced_accuracy",
    scenario_weights = c(proportional = 0.50, correlation = 0.50),
    tie_break = "conservative",
    allow_infeasible = FALSE
  )
  for (key in c("gic", "bic")) {
    source <- readRDS(file.path(paths$inputs, paste0(key, "-tuning.rds")))
    source_selection <- do.call(
      selectTunedSearchParameters,
      c(list(tuning_grid = source$grid), policy)
    )
    compact_selection <- do.call(
      selectTunedSearchParameters,
      c(list(tuning_grid = updated$tuning_grids[[key]]), policy)
    )
    testthat::expect_equal(
      compact_selection$selected_row,
      source_selection$selected_row
    )
  }

  selected <- updated$tuning$selected
  testthat::expect_equal(selected$Threshold, c(20, 10))
  testthat::expect_equal(selected$`Min clade`, c(10L, 10L))
  testthat::expect_true(all(c("Null any FP", "Score", "Status") %in%
    names(updated$tuning$gic)))
  testthat::expect_identical(
    updated$provenance$tuning$generation_accounting$n_generated,
    c(500L, 500L, 500L)
  )
  testthat::expect_identical(
    updated$provenance$tuning$total_evaluated_searches,
    18000L
  )
  testthat::expect_identical(
    updated$provenance$tuning$planted_clade_counts,
    data.frame(
      scenario = c("null", "proportional", "integration-rate"),
      n_true_shifts_10_19 = c(0L, 1411L, 1467L),
      n_true_shifts_20_40 = c(0L, 1089L, 1033L)
    )
  )

  required <- c(
    "null_mean_false_positive_rate", "null_fraction_any_false_positive",
    "null_evaluable_fraction", "proportional_fuzzy_balanced_accuracy",
    "correlation_fuzzy_balanced_accuracy"
  )
  testthat::expect_true(all(vapply(
    updated$grid_summary,
    function(x) all(is.finite(unlist(x[, required], use.names = FALSE))),
    logical(1L)
  )))
  serialized <- as.character(unlist(updated, recursive = TRUE, use.names = TRUE))
  testthat::expect_false(any(grepl(
    "(^|[= ])(/Users/|/home/|/private/tmp/|[A-Za-z]:[/\\\\])",
    serialized
  )))
})

test_that("tracked compact grids can drive an independent schema-4 export", {
  paths <- paired_tuning_paths()
  helpers <- load_paired_tuning_helpers()
  cache <- readRDS(paths$artifact)
  inputs <- withr::local_tempdir()
  write_paired_tuning_fixture(cache, inputs)

  rebuilt <- helpers$build_paired_tuning_cache(inputs)
  testthat::expect_no_error(helpers$validate_paired_tuning_cache(rebuilt))
  testthat::expect_identical(rebuilt$fixed_settings,
    helpers$fixed_setting_rows(rebuilt$tuning_grids))
  testthat::expect_equal(rebuilt$tuning$selected, cache$tuning$selected)
})

test_that("fixed-setting rows map independent metric values from shuffled settings", {
  helpers <- load_paired_tuning_helpers()
  # Every recovery metric has a distinct value, so swapping strict/fuzzy
  # columns or scenarios cannot pass by sharing the production mapping helper.
  values <- c(
    null_mean_false_positive_rate = .01, null_fraction_any_false_positive = .02,
    proportional_strict_precision = .11, proportional_strict_recall = .12,
    proportional_strict_f1 = .13, proportional_strict_specificity = .14,
    proportional_strict_fpr = .15, proportional_strict_balanced_accuracy = .16,
    proportional_fuzzy_precision = .17, proportional_fuzzy_recall = .18,
    proportional_fuzzy_f1 = .19, proportional_fuzzy_specificity = .20,
    proportional_fuzzy_fpr = .21, proportional_fuzzy_balanced_accuracy = .22,
    proportional_weighted_fuzzy_f1 = .23,
    correlation_strict_precision = .31, correlation_strict_recall = .32,
    correlation_strict_f1 = .33, correlation_strict_specificity = .34,
    correlation_strict_fpr = .35, correlation_strict_balanced_accuracy = .36,
    correlation_fuzzy_precision = .37, correlation_fuzzy_recall = .38,
    correlation_fuzzy_f1 = .39, correlation_fuzzy_specificity = .40,
    correlation_fuzzy_fpr = .41, correlation_fuzzy_balanced_accuracy = .42,
    correlation_weighted_fuzzy_f1 = .43,
    null_mean_inferred_shifts = .1, proportional_mean_inferred_shifts = 3,
    correlation_mean_inferred_shifts = 4
  )
  make_grid <- function(ic, offset, order) {
    summary <- data.frame(shift_acceptance_threshold = c(20, 10, 10),
                          min_descendant_tips = c(10L, 20L, 10L))
    for (name in names(values)) summary[[name]] <- c(-1, -2, values[[name]] + offset)
    for (scenario in c("null", "proportional", "correlation")) {
      summary[[paste0(scenario, "_evaluable_fraction")]] <- c(0, 0, 1)
      summary[[paste0(scenario, "_completion_rate")]] <- c(0, 0, 1)
      summary[[paste0(scenario, "_failure_rate")]] <- c(1, 1, 0)
    }
    list(IC = ic, summary_table = summary[order, , drop = FALSE])
  }
  grids <- list(gic = make_grid("GIC", 0, c(1, 3, 2)),
                bic = make_grid("BIC", .4, c(3, 2, 1)))
  expected <- data.frame(
    IC = rep(c("GIC", "BIC"), each = 3),
    Scenario = rep(c("Null", "Proportional", "Integration-rate"), 2),
    `Mean FP` = c(.01, NA, NA, .41, NA, NA),
    `Any FP` = c(.02, NA, NA, .42, NA, NA),
    `Strict precision` = c(NA, .11, .31, NA, .51, .71),
    `Strict recall` = c(NA, .12, .32, NA, .52, .72),
    `Strict F1` = c(NA, .13, .33, NA, .53, .73),
    `Strict specificity` = c(NA, .14, .34, NA, .54, .74),
    `Strict FPR` = c(NA, .15, .35, NA, .55, .75),
    `Strict balanced accuracy` = c(NA, .16, .36, NA, .56, .76),
    `Fuzzy precision` = c(NA, .17, .37, NA, .57, .77),
    `Fuzzy recall` = c(NA, .18, .38, NA, .58, .78),
    `Fuzzy F1` = c(NA, .19, .39, NA, .59, .79),
    `Fuzzy specificity` = c(NA, .20, .40, NA, .60, .80),
    `Fuzzy FPR` = c(NA, .21, .41, NA, .61, .81),
    `Fuzzy balanced accuracy` = c(NA, .22, .42, NA, .62, .82),
    `Weighted fuzzy F1` = c(NA, .23, .43, NA, .63, .83),
    `Mean shifts` = c(.1, 3, 4, .5, 3.4, 4.4),
    Evaluable = 1, Completion = 1, Failure = 0, check.names = FALSE
  )
  testthat::expect_equal(helpers$fixed_setting_rows(grids), expected)
})

test_that("cache writer round-trips valid data and replaces an existing file", {
  helpers <- load_paired_tuning_helpers()
  cache <- readRDS(paired_tuning_paths()$artifact)
  directory <- withr::local_tempdir()
  output <- file.path(directory, "nested", "cache.rds")
  helpers$write_paired_tuning_cache(cache, output)
  testthat::expect_identical(readRDS(output), cache)

  saveRDS(list(sentinel = "replace me"), output)
  helpers$write_paired_tuning_cache(cache, output)
  testthat::expect_identical(readRDS(output), cache)
  testthat::expect_identical(list.files(dirname(output), all.files = TRUE,
                                      no.. = TRUE), "cache.rds")
})

test_that("cache writer preserves existing bytes when validation fails", {
  helpers <- load_paired_tuning_helpers()
  cache <- readRDS(paired_tuning_paths()$artifact)
  directory <- withr::local_tempdir()
  output <- file.path(directory, "cache.rds")
  sentinel <- list(sentinel = "do not replace")
  saveRDS(sentinel, output)
  before <- digest::digest(file = output, algo = "sha256")
  cache$fixed_settings$`Fuzzy F1`[[2L]] <- .123
  testthat::expect_error(helpers$write_paired_tuning_cache(cache, output),
                         "fixed settings disagree")
  testthat::expect_identical(digest::digest(file = output, algo = "sha256"), before)
  testthat::expect_identical(readRDS(output), sentinel)
  testthat::expect_identical(list.files(directory, all.files = TRUE, no.. = TRUE),
                            "cache.rds")
})

test_that("completed-campaign counts cannot coexist with incomplete searches", {
  helpers <- load_paired_tuning_helpers()
  cache <- readRDS(paired_tuning_paths()$artifact)
  for (ic in c("GIC", "BIC")) {
    for (scenario in c("null", "proportional", "correlation")) {
      directory <- withr::local_tempdir()
      write_paired_tuning_fixture(cache, directory)
      path <- file.path(directory, paste0(tolower(ic), "-tuning.rds"))
      wrapper <- readRDS(path)
      wrapper$grid$summary_table[[paste0(scenario, "_completion_rate")]][1] <- .9
      wrapper$grid$summary_table[[paste0(scenario, "_failure_rate")]][1] <- .1
      saveRDS(wrapper, path)
      suite_path <- file.path(directory, "suite-summary.rds")
      suite <- readRDS(suite_path)
      suite_scenario <- if (scenario == "correlation") "integration-rate" else scenario
      row <- suite$summary_table$IC == ic &
        suite$summary_table$scenario == suite_scenario &
        suite$summary_table$setting_id == wrapper$grid$summary_table$setting_id[1]
      suite$summary_table$completion_rate[row] <- .9
      saveRDS(suite, suite_path)
      # Both sources agree on rates, but still claim 500/500 successful searches.
      testthat::expect_error(helpers$build_paired_tuning_cache(directory),
                             "completion rates contradict.*500/500")
    }
  }
})

test_that("cache validation rejects tampered fixed metrics and provenance", {
  paths <- paired_tuning_paths()
  helpers <- load_paired_tuning_helpers()
  cache <- readRDS(paths$artifact)

  tampered <- cache
  tampered$fixed_settings$`Fuzzy F1`[[2L]] <-
    tampered$fixed_settings$`Fuzzy F1`[[2L]] + 0.01
  testthat::expect_error(
    helpers$validate_paired_tuning_cache(tampered),
    "fixed settings disagree"
  )

  tampered <- cache
  tampered$provenance$fixed_settings$n_replicates_per_setting <- 100L
  testthat::expect_error(
    helpers$validate_paired_tuning_cache(tampered),
    "fixed-settings provenance disagrees"
  )
})

test_that("fixed-setting derivation requires one 10/10 setting per IC", {
  paths <- paired_tuning_paths()
  helpers <- load_paired_tuning_helpers()
  cache <- readRDS(paths$artifact)

  missing <- cache$tuning_grids
  missing$gic$summary_table$shift_acceptance_threshold[[1L]] <- 11
  testthat::expect_error(
    helpers$fixed_setting_rows(missing),
    "exactly one threshold-10/min-clade-10"
  )

  duplicate <- cache$tuning_grids
  duplicate$bic$summary_table$shift_acceptance_threshold[[2L]] <- 10
  testthat::expect_error(
    helpers$fixed_setting_rows(duplicate),
    "exactly one threshold-10/min-clade-10"
  )
})

test_that("paired exporter rejects contradictory generation and provenance inputs", {
  paths <- paired_tuning_paths()
  helpers <- load_paired_tuning_helpers()
  cache <- readRDS(paths$artifact)
  copied <- withr::local_tempdir()
  write_paired_tuning_fixture(cache, copied)
  accounting_path <- file.path(copied, "generation-accounting.csv")
  accounting <- utils::read.csv(accounting_path, stringsAsFactors = FALSE)
  accounting$n_generated[[2L]] <- 499L
  utils::write.csv(accounting, accounting_path, row.names = FALSE)
  testthat::expect_error(
    helpers$build_paired_tuning_cache(copied),
    "generation accounting"
  )

  copied <- withr::local_tempdir()
  write_paired_tuning_fixture(cache, copied)
  bic_path <- file.path(copied, "bic-tuning.rds")
  bic <- readRDS(bic_path)
  bic$provenance$source_fingerprint <- paste(rep("0", 64L), collapse = "")
  saveRDS(bic, bic_path)
  testthat::expect_error(
    helpers$build_paired_tuning_cache(copied),
    "matching paired.*provenance"
  )

  copied <- withr::local_tempdir()
  write_paired_tuning_fixture(cache, copied)
  gic_path <- file.path(copied, "gic-tuning.rds")
  gic <- readRDS(gic_path)
  gic$grid$summary_table$proportional_fuzzy_balanced_accuracy[[1L]] <- 0.5
  saveRDS(gic, gic_path)
  testthat::expect_error(
    helpers$build_paired_tuning_cache(copied),
    "does not correspond to the suite summary"
  )

  copied <- withr::local_tempdir()
  write_paired_tuning_fixture(cache, copied)
  gic_path <- file.path(copied, "gic-tuning.rds")
  gic <- readRDS(gic_path)
  gic$grid$grid$shift_acceptance_threshold[[1L]] <- 11
  saveRDS(gic, gic_path)
  testthat::expect_error(
    helpers$build_paired_tuning_cache(copied),
    "incomplete paired design"
  )

  copied <- withr::local_tempdir()
  write_paired_tuning_fixture(cache, copied)
  bic_path <- file.path(copied, "bic-tuning.rds")
  bic <- readRDS(bic_path)
  bic$grid$base_search_options$error <- TRUE
  saveRDS(bic, bic_path)
  testthat::expect_error(
    helpers$build_paired_tuning_cache(copied),
    "base search options"
  )

  copied <- withr::local_tempdir()
  write_paired_tuning_fixture(cache, copied)
  gic_path <- file.path(copied, "gic-tuning.rds")
  gic <- readRDS(gic_path)
  gic$grid$simulation_generators[["null"]] <- "original"
  saveRDS(gic, gic_path)
  testthat::expect_error(
    helpers$build_paired_tuning_cache(copied),
    "simulation generators"
  )

  copied <- withr::local_tempdir()
  write_paired_tuning_fixture(cache, copied)
  for (key in c("gic", "bic")) {
    path <- file.path(copied, paste0(key, "-tuning.rds"))
    wrapper <- readRDS(path)
    wrapper$provenance$package_remote_sha <- NULL
    saveRDS(wrapper, path)
  }
  testthat::expect_error(
    helpers$build_paired_tuning_cache(copied),
    "package commit"
  )

  copied <- withr::local_tempdir()
  write_paired_tuning_fixture(cache, copied)
  suite_path <- file.path(copied, "suite-summary.rds")
  suite <- readRDS(suite_path)
  suite$design_fingerprint <- NULL
  saveRDS(suite, suite_path)
  testthat::expect_error(
    helpers$build_paired_tuning_cache(copied),
    "design fingerprint"
  )

  copied <- withr::local_tempdir()
  write_paired_tuning_fixture(cache, copied)
  gic_path <- file.path(copied, "gic-tuning.rds")
  gic <- readRDS(gic_path)
  gic$grid$summary_table$null_failure_rate[[1L]] <- 0.25
  saveRDS(gic, gic_path)
  testthat::expect_error(
    helpers$build_paired_tuning_cache(copied),
    "failure rates"
  )

  copied <- withr::local_tempdir()
  write_paired_tuning_fixture(cache, copied)
  suite_path <- file.path(copied, "suite-summary.rds")
  suite <- readRDS(suite_path)
  suite$summary_table$n_true_shifts_10_19[[13L]] <- 1410L
  saveRDS(suite, suite_path)
  testthat::expect_error(
    helpers$build_paired_tuning_cache(copied),
    "planted-clade counts are inconsistent"
  )
})
