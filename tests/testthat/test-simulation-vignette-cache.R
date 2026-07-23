test_that("simulation vignette cache records the empirical generator", {
  cache_path <- testthat::test_path(
    "../../inst/extdata/simulation-study-cache/passerine_preview_tables.rds"
  )
  testthat::skip_if_not(file.exists(cache_path), "Simulation vignette cache not present")

  cache <- readRDS(cache_path)

  testthat::expect_identical(cache$schema_version, 3L)
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
  recovery_columns <- c(
    "Strict precision",
    "Strict recall",
    "Strict F1",
    "Strict specificity",
    "Strict FPR",
    "Strict balanced accuracy",
    "Fuzzy precision",
    "Fuzzy recall",
    "Fuzzy F1",
    "Fuzzy specificity",
    "Fuzzy FPR",
    "Fuzzy balanced accuracy"
  )
  has_recovery_columns <- all(recovery_columns %in% names(cache$fixed_settings))
  testthat::expect_true(has_recovery_columns)
  if (has_recovery_columns) {
    recovery_rows <- cache$fixed_settings[
      cache$fixed_settings$Scenario %in% c("Proportional", "Integration-rate"),
      recovery_columns,
      drop = FALSE
    ]
    testthat::expect_false(anyNA(recovery_rows))
  }

  tuning_metric_columns <- c(
    "Prop. Fuzzy balanced accuracy",
    "Integration Fuzzy balanced accuracy"
  )
  testthat::expect_true(all(vapply(
    cache$tuning[c("gic", "bic")],
    function(x) all(tuning_metric_columns %in% names(x)),
    logical(1L)
  )))
  selected_metric_columns <- c(tuning_metric_columns, "Score")
  has_selected_metrics <- all(selected_metric_columns %in%
    names(cache$tuning$selected))
  testthat::expect_true(has_selected_metrics)
  if (has_selected_metrics) {
    testthat::expect_false(anyNA(cache$tuning$selected[
      , selected_metric_columns,
      drop = FALSE
    ]))
    testthat::expect_equal(
      cache$tuning$selected$Score,
      rowMeans(cache$tuning$selected[, tuning_metric_columns, drop = FALSE])
    )
    for (ic in cache$tuning$selected$IC) {
      selected <- cache$tuning$selected[cache$tuning$selected$IC == ic, ]
      fixed <- cache$fixed_settings[cache$fixed_settings$IC == ic, ]
      testthat::expect_equal(
        fixed$`Fuzzy balanced accuracy`[
          fixed$Scenario == "Proportional"
        ],
        selected$`Prop. Fuzzy balanced accuracy`
      )
      testthat::expect_equal(
        fixed$`Fuzzy balanced accuracy`[
          fixed$Scenario == "Integration-rate"
        ],
        selected$`Integration Fuzzy balanced accuracy`
      )
    }
  }

  testthat::expect_setequal(names(cache$grid_summary), c("gic", "bic"))
  testthat::expect_true(all(vapply(cache$grid_summary, nrow, integer(1L)) == 6L))
  raw_metric_columns <- c(
    paste0("proportional_", c(
      "strict_specificity",
      "strict_fpr",
      "strict_balanced_accuracy",
      "fuzzy_specificity",
      "fuzzy_fpr",
      "fuzzy_balanced_accuracy"
    )),
    paste0("correlation_", c(
      "strict_specificity",
      "strict_fpr",
      "strict_balanced_accuracy",
      "fuzzy_specificity",
      "fuzzy_fpr",
      "fuzzy_balanced_accuracy"
    ))
  )
  testthat::expect_true(all(vapply(
    cache$grid_summary,
    function(x) all(raw_metric_columns %in% names(x)),
    logical(1L)
  )))

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

simulation_regeneration_paths <- function() {
  c(
    runner = testthat::test_path(
      "../../data-raw/run_simulation_study_vignette_grids.R"
    ),
    builder = testthat::test_path(
      "../../data-raw/build_simulation_study_vignette_cache.R"
    )
  )
}

skip_if_simulation_regeneration_scripts_absent <- function() {
  paths <- simulation_regeneration_paths()
  testthat::skip_if_not(
    all(file.exists(paths)),
    "Development-only simulation regeneration scripts are excluded from source tarballs"
  )
  paths
}

test_that("data-raw regeneration scripts are sourceable without executing", {
  paths <- skip_if_simulation_regeneration_scripts_absent()
  runner_env <- new.env(parent = globalenv())
  builder_env <- new.env(parent = globalenv())

  runner_result <- try(
    sys.source(paths[["runner"]], envir = runner_env),
    silent = TRUE
  )
  builder_result <- try(
    sys.source(paths[["builder"]], envir = builder_env),
    silent = TRUE
  )

  testthat::expect_false(inherits(runner_result, "try-error"))
  testthat::expect_false(inherits(builder_result, "try-error"))
  testthat::expect_true(exists("parse_cli_args", envir = runner_env))
  testthat::expect_true(exists("read_grid_run", envir = builder_env))
})

load_simulation_regeneration_helpers <- function() {
  paths <- skip_if_simulation_regeneration_scripts_absent()
  helper_env <- new.env(parent = globalenv())
  sys.source(paths[["runner"]], envir = helper_env)
  sys.source(paths[["builder"]], envir = helper_env)
  helper_env
}

make_tuning_grid_fixture <- function(ic, design, selected_setting = 1L) {
  settings <- expand.grid(
    shift_acceptance_threshold = design$shift_acceptance_thresholds,
    min_descendant_tips = design$min_descendant_tips_values,
    KEEP.OUT.ATTRS = FALSE
  )
  n_settings <- nrow(settings)
  scores <- rep(0.60, n_settings)
  scores[[selected_setting]] <- 0.90
  summary <- data.frame(
    setting_id = seq_len(n_settings),
    IC = rep(ic, n_settings),
    settings,
    null_seed = seq_len(n_settings),
    proportional_seed = seq_len(n_settings) + 100L,
    correlation_seed = seq_len(n_settings) + 200L,
    null_mean_false_positive_rate = rep(0, n_settings),
    null_fraction_any_false_positive = rep(0, n_settings),
    null_evaluable_fraction = rep(1, n_settings),
    null_completion_rate = rep(1, n_settings),
    null_failure_rate = rep(0, n_settings),
    null_mean_inferred_shifts = rep(0, n_settings),
    stringsAsFactors = FALSE
  )
  metric_suffixes <- c(
    "strict_precision", "strict_recall", "strict_f1",
    "strict_specificity", "strict_fpr", "strict_balanced_accuracy",
    "fuzzy_precision", "fuzzy_recall", "fuzzy_f1",
    "fuzzy_specificity", "fuzzy_fpr", "fuzzy_balanced_accuracy"
  )
  for (scenario in c("proportional", "correlation")) {
    for (suffix in metric_suffixes) {
      metric <- sub("^(strict|fuzzy)_", "", suffix)
      balanced_accuracy <- if (startsWith(suffix, "fuzzy_")) {
        scores
      } else {
        rep(0.70, n_settings)
      }
      specificity <- rep(0.95, n_settings)
      recall <- 2 * balanced_accuracy - specificity
      values <- switch(
        metric,
        precision = recall,
        recall = recall,
        f1 = recall,
        specificity = specificity,
        fpr = 1 - specificity,
        balanced_accuracy = balanced_accuracy
      )
      summary[[paste0(scenario, "_", suffix)]] <- values
    }
    summary[[paste0(scenario, "_evaluable_fraction")]] <- rep(1, n_settings)
    summary[[paste0(scenario, "_completion_rate")]] <- rep(1, n_settings)
    summary[[paste0(scenario, "_failure_rate")]] <- rep(0, n_settings)
    summary[[paste0(scenario, "_mean_inferred_shifts")]] <- rep(5, n_settings)
    summary[[paste0(scenario, "_weighted_fuzzy_f1")]] <- scores
  }

  grid <- list(
    IC = ic,
    grid = data.frame(setting_id = seq_len(n_settings), settings),
    summary_table = summary,
    simulation_generators = c(
      null = "empirical",
      proportional = "empirical",
      correlation = "empirical"
    ),
    base_search_options = design$base_search_options,
    null_replicates = design$null_replicates,
    recovery_replicates = design$recovery_replicates,
    fuzzy_distance = design$fuzzy_distance,
    weighted = design$weighted,
    store_studies = design$store_studies
  )
  class(grid) <- c("bifrost_search_tuning_grid", class(grid))
  grid
}

test_that("grid runner CLI parsing is explicit and deterministic", {
  helpers <- load_simulation_regeneration_helpers()
  required <- c("parse_cli_args", "simulation_grid_design")
  has_helpers <- all(vapply(required, exists, logical(1L), envir = helpers))
  testthat::expect_true(has_helpers)
  if (!has_helpers) {
    return(invisible(NULL))
  }

  cli <- helpers$parse_cli_args(
    c(
      "--mode=full",
      "--ic=GIC",
      "--cores=3",
      "--output-dir=fixture-output",
      "--overwrite"
    ),
    detect_cores = function() 8L
  )
  testthat::expect_identical(cli$mode, "full")
  testthat::expect_identical(cli$ic, "GIC")
  testthat::expect_identical(cli$cores, 3L)
  testthat::expect_identical(cli$output_dir, "fixture-output")
  testthat::expect_true(cli$overwrite)
  testthat::expect_identical(
    helpers$parse_cli_args(
      "--mode=full",
      detect_cores = function() 8L
    )$cores,
    6L
  )
  testthat::expect_identical(
    helpers$parse_cli_args(
      "--mode=smoke",
      detect_cores = function() 8L
    )$cores,
    1L
  )
  testthat::expect_error(helpers$parse_cli_args(character()), "--mode")
  testthat::expect_error(helpers$parse_cli_args("--mode=bad"), "smoke.*full")
  testthat::expect_error(
    helpers$parse_cli_args(c("--mode=smoke", "--cores=0")),
    "positive integer"
  )
  testthat::expect_error(
    helpers$parse_cli_args(c("--mode=smoke", "--unknown")),
    "Unknown argument"
  )
})

test_that("source provenance hashes relevant files and ignores unrelated dirt", {
  helpers <- load_simulation_regeneration_helpers()
  required <- c("collect_source_provenance", "relevant_source_files")
  has_helpers <- all(vapply(required, exists, logical(1L), envir = helpers))
  testthat::expect_true(has_helpers)
  if (!has_helpers) {
    return(invisible(NULL))
  }

  first <- tempfile(fileext = ".R")
  second <- tempfile(fileext = ".R")
  writeLines("first <- 1", first)
  writeLines("second <- 2", second)
  files <- c(
    "R/simulation-first.R" = first,
    "R/simulation-second.R" = second
  )
  clean <- helpers$collect_source_provenance(
    relevant_files = files,
    commit = "fixture-head",
    status_lines = " M _pkgdown.yml"
  )
  testthat::expect_identical(clean$commit, "fixture-head")
  testthat::expect_identical(names(clean$source_md5), names(files))
  testthat::expect_false(clean$source_dirty)

  changed <- clean
  writeLines("first <- 3", first)
  changed <- helpers$collect_source_provenance(
    relevant_files = files,
    commit = "fixture-head",
    status_lines = " M R/simulation-first.R"
  )
  testthat::expect_true(changed$source_dirty)
  testthat::expect_identical(
    changed$source_status,
    " M R/simulation-first.R"
  )
  testthat::expect_false(identical(
    clean$source_md5[["R/simulation-first.R"]],
    changed$source_md5[["R/simulation-first.R"]]
  ))

  actual_files <- helpers$relevant_source_files(
    testthat::test_path("../..")
  )
  testthat::expect_true(all(c(
    "data-raw/run_simulation_study_vignette_grids.R",
    "DESCRIPTION",
    "NAMESPACE",
    "inst/extdata/avian-skeleton/passerine_bodyplan_tree.tre",
    "inst/extdata/avian-skeleton/passerine_bodyplan_data.RDS"
  ) %in% names(actual_files)))
  testthat::expect_true(any(grepl("^R/.*\\.R$", names(actual_files))))
})

test_that("grid wrappers reject incomplete designs and mismatched provenance", {
  helpers <- load_simulation_regeneration_helpers()
  required <- c(
    "simulation_grid_design", "new_grid_wrapper", "validate_grid_wrapper"
  )
  has_helpers <- all(vapply(required, exists, logical(1L), envir = helpers))
  testthat::expect_true(has_helpers)
  if (!has_helpers) {
    return(invisible(NULL))
  }

  design <- helpers$simulation_grid_design("full")
  source_provenance <- list(
    commit = "fixture-head",
    source_md5 = c(
      "R/simulation.R" = "aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa"
    ),
    source_dirty = FALSE,
    source_status = character()
  )
  wrapper <- helpers$new_grid_wrapper(
    make_tuning_grid_fixture("GIC", design, selected_setting = 2L),
    design,
    source_provenance,
    setting_workers = 2L
  )
  testthat::expect_no_error(helpers$validate_grid_wrapper(
    wrapper,
    expected_ic = "GIC",
    expected_design = design,
    expected_source_provenance = source_provenance
  ))

  wrong_design <- wrapper
  wrong_design$provenance$design$tree_tip_count <- 249L
  testthat::expect_error(
    helpers$validate_grid_wrapper(
      wrong_design, "GIC", design, source_provenance
    ),
    "design"
  )
  wrong_settings <- wrapper
  wrong_settings$result$summary_table$shift_acceptance_threshold[[1L]] <- 99
  testthat::expect_error(
    helpers$validate_grid_wrapper(
      wrong_settings, "GIC", design, source_provenance
    ),
    "setting"
  )
  wrong_generator <- wrapper
  wrong_generator$result$simulation_generators[[1L]] <- "original"
  testthat::expect_error(
    helpers$validate_grid_wrapper(
      wrong_generator, "GIC", design, source_provenance
    ),
    "empirical generator"
  )
  missing_metric <- wrapper
  missing_metric$result$summary_table$proportional_strict_precision <- NULL
  testthat::expect_error(
    helpers$validate_grid_wrapper(
      missing_metric, "GIC", design, source_provenance
    ),
    "metrics"
  )
  out_of_range <- wrapper
  out_of_range$result$summary_table$proportional_strict_precision[[1L]] <- 999
  testthat::expect_error(
    helpers$validate_grid_wrapper(
      out_of_range, "GIC", design, source_provenance
    ),
    "between 0 and 1"
  )
  inconsistent_fpr <- wrapper
  inconsistent_fpr$result$summary_table$proportional_fuzzy_fpr[[1L]] <- 0.50
  testthat::expect_error(
    helpers$validate_grid_wrapper(
      inconsistent_fpr, "GIC", design, source_provenance
    ),
    "FPR.*specificity"
  )
  inconsistent_ba <- wrapper
  inconsistent_ba$result$summary_table$correlation_fuzzy_balanced_accuracy[[1L]] <- 0.10
  testthat::expect_error(
    helpers$validate_grid_wrapper(
      inconsistent_ba, "GIC", design, source_provenance
    ),
    "balanced accuracy"
  )
  invalid_completion <- wrapper
  invalid_completion$result$summary_table$null_completion_rate[[1L]] <- -7
  testthat::expect_error(
    helpers$validate_grid_wrapper(
      invalid_completion, "GIC", design, source_provenance
    ),
    "between 0 and 1"
  )
  wrong_source <- source_provenance
  wrong_source$source_md5[[1L]] <- "bbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb"
  testthat::expect_error(
    helpers$validate_grid_wrapper(
      wrapper, "GIC", design, wrong_source
    ),
    "provenance"
  )
})

test_that("cache builder rejects absent or malformed source provenance", {
  helpers <- load_simulation_regeneration_helpers()
  required <- c(
    "simulation_grid_design", "new_grid_wrapper", "read_grid_run"
  )
  has_helpers <- all(vapply(required, exists, logical(1L), envir = helpers))
  testthat::expect_true(has_helpers)
  if (!has_helpers) {
    return(invisible(NULL))
  }

  design <- helpers$simulation_grid_design("full")
  source_provenance <- list(
    commit = "fixture-head",
    source_md5 = c(
      "R/simulation.R" = "aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa"
    ),
    source_dirty = FALSE,
    source_status = character()
  )
  wrapper <- helpers$new_grid_wrapper(
    make_tuning_grid_fixture("GIC", design),
    design,
    source_provenance,
    setting_workers = 1L
  )
  variants <- list(
    omitted = within(wrapper, {
      provenance$commit <- NULL
      provenance$source_files <- NULL
      provenance$source_md5 <- NULL
      provenance$source_dirty <- NULL
      provenance$source_status <- NULL
    }),
    malformed_hash = within(wrapper, {
      provenance$source_md5[[1L]] <- "not-an-md5"
    }),
    unnamed_hash = within(wrapper, {
      provenance$source_md5 <- unname(provenance$source_md5)
    }),
    empty_inventory = within(wrapper, {
      provenance$source_files <- character()
      provenance$source_md5 <- stats::setNames(character(), character())
    }),
    dirty_without_status = within(wrapper, {
      provenance$source_dirty <- TRUE
      provenance$source_status <- character()
    })
  )

  for (name in names(variants)) {
    path <- tempfile(fileext = ".rds")
    saveRDS(variants[[name]], path)
    testthat::expect_error(
      helpers$read_grid_run(path, "GIC"),
      "source provenance",
      info = name
    )
  }
})

test_that("valid grid outputs resume while overwrite reruns", {
  helpers <- load_simulation_regeneration_helpers()
  required <- c(
    "simulation_grid_design", "new_grid_wrapper", "run_or_resume_grid"
  )
  has_helpers <- all(vapply(required, exists, logical(1L), envir = helpers))
  testthat::expect_true(has_helpers)
  if (!has_helpers) {
    return(invisible(NULL))
  }

  design <- helpers$simulation_grid_design("smoke")
  source_provenance <- list(
    commit = "fixture-head",
    source_md5 = c(
      "R/simulation.R" = "aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa"
    ),
    source_dirty = FALSE,
    source_status = character()
  )
  grid <- make_tuning_grid_fixture("BIC", design)
  wrapper <- helpers$new_grid_wrapper(
    grid,
    design,
    source_provenance,
    setting_workers = 1L
  )
  output_path <- tempfile(fileext = ".rds")
  saveRDS(wrapper, output_path)
  run_count <- 0L
  run_grid <- function() {
    run_count <<- run_count + 1L
    grid
  }

  resumed <- helpers$run_or_resume_grid(
    output_path,
    expected_ic = "BIC",
    design = design,
    source_provenance = source_provenance,
    setting_workers = 1L,
    overwrite = FALSE,
    run_grid = run_grid
  )
  testthat::expect_identical(resumed$status, "resumed")
  testthat::expect_identical(run_count, 0L)

  mismatched <- wrapper
  mismatched$provenance$source_md5[[1L]] <-
    "bbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb"
  saveRDS(mismatched, output_path)
  testthat::expect_error(
    helpers$run_or_resume_grid(
      output_path, "BIC", design, source_provenance, 1L, FALSE, run_grid
    ),
    "provenance"
  )
  overwritten <- helpers$run_or_resume_grid(
    output_path,
    expected_ic = "BIC",
    design = design,
    source_provenance = source_provenance,
    setting_workers = 1L,
    overwrite = TRUE,
    run_grid = run_grid
  )
  testthat::expect_identical(overwritten$status, "written")
  testthat::expect_identical(run_count, 1L)
  testthat::expect_identical(
    readRDS(output_path)$provenance$source_md5,
    source_provenance$source_md5
  )
})

test_that("cache builder uses selected fixture rows and validates before writing", {
  helpers <- load_simulation_regeneration_helpers()
  required <- c(
    "simulation_grid_design", "new_grid_wrapper", "build_simulation_cache",
    "validate_cache_object", "write_simulation_cache"
  )
  has_helpers <- all(vapply(required, exists, logical(1L), envir = helpers))
  testthat::expect_true(has_helpers)
  if (!has_helpers) {
    return(invisible(NULL))
  }

  design <- helpers$simulation_grid_design("full")
  source_files <- helpers$approved_generator_source_files()
  source_provenance <- list(
    commit = paste(rep("a", 40L), collapse = ""),
    source_md5 = stats::setNames(
      rep("aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa", length(source_files)),
      source_files
    ),
    source_dirty = FALSE,
    source_status = character()
  )
  gic_wrapper <- helpers$new_grid_wrapper(
    make_tuning_grid_fixture("GIC", design, selected_setting = 2L),
    design,
    source_provenance,
    setting_workers = 2L
  )
  bic_wrapper <- helpers$new_grid_wrapper(
    make_tuning_grid_fixture("BIC", design, selected_setting = 4L),
    design,
    source_provenance,
    setting_workers = 2L
  )
  gic_path <- tempfile(fileext = ".rds")
  bic_path <- tempfile(fileext = ".rds")
  saveRDS(gic_wrapper, gic_path)
  saveRDS(bic_wrapper, bic_path)

  cache <- helpers$build_simulation_cache(gic_path, bic_path)
  testthat::expect_no_error(helpers$validate_cache_object(cache))
  testthat::expect_identical(cache$tuning$selected$Threshold, c(20, 10))
  testthat::expect_identical(cache$tuning$selected$`Min clade`, c(10L, 20L))

  out_path <- tempfile(fileext = ".rds")
  saveRDS(list(sentinel = TRUE), out_path)
  helpers$write_simulation_cache(cache, out_path)
  testthat::expect_identical(readRDS(out_path)$schema_version, 3L)

  invalid <- cache
  invalid$tuning$selected$Score[[1L]] <- NA_real_
  saveRDS(list(sentinel = TRUE), out_path)
  testthat::expect_error(
    helpers$write_simulation_cache(invalid, out_path),
    "finite selected"
  )
  testthat::expect_identical(readRDS(out_path), list(sentinel = TRUE))

  inconsistent <- cache
  inconsistent$fixed_settings$`Fuzzy balanced accuracy`[
    inconsistent$fixed_settings$IC == "GIC" &
      inconsistent$fixed_settings$Scenario == "Proportional"
  ] <- 0
  testthat::expect_error(
    helpers$validate_cache_object(inconsistent),
    "selected.*fixed"
  )

  invalid_raw <- cache
  invalid_raw$grid_summary$gic$proportional_strict_precision[[1L]] <- 999
  testthat::expect_error(
    helpers$validate_cache_object(invalid_raw),
    "between 0 and 1"
  )

  invalid_tuning <- cache
  invalid_tuning$tuning$gic$`Null FP`[[1L]] <- 0.99
  testthat::expect_error(
    helpers$validate_cache_object(invalid_tuning),
    "tuning table.*raw grid"
  )

  invalid_design <- cache
  invalid_design$provenance$n_replicates_per_setting <- 1L
  testthat::expect_error(
    helpers$validate_cache_object(invalid_design),
    "design contract"
  )

  invalid_policy <- cache
  invalid_policy$provenance$selection$max_any_false_positive <- 0.99
  testthat::expect_error(
    helpers$validate_cache_object(invalid_policy),
    "selection policy"
  )

  missing_environment <- cache
  missing_environment$provenance$R <- NULL
  testthat::expect_error(
    helpers$validate_cache_object(missing_environment),
    "R and platform"
  )

  invalid_commit <- cache
  invalid_commit$provenance$generator_commit <- "not-a-commit"
  testthat::expect_error(
    helpers$validate_cache_object(invalid_commit),
    "generator source provenance"
  )

  invalid_inventory <- cache
  invalid_inventory$provenance$generator_source_md5 <- unname(
    invalid_inventory$provenance$generator_source_md5
  )
  testthat::expect_error(
    helpers$validate_cache_object(invalid_inventory),
    "generator source provenance"
  )

  invalid_generator_hash <- cache
  invalid_generator_hash$provenance$generator_source_md5[[1L]] <- "bad"
  testthat::expect_error(
    helpers$validate_cache_object(invalid_generator_hash),
    "generator source provenance"
  )

  invalid_grid_hashes <- cache
  invalid_grid_hashes$provenance$source_md5[[1L]] <- "bad"
  testthat::expect_error(
    helpers$validate_cache_object(invalid_grid_hashes),
    "input grid hashes"
  )

  misnamed_grid_hashes <- cache
  names(misnamed_grid_hashes$provenance$source_md5) <- c("GIC", "Other")
  testthat::expect_error(
    helpers$validate_cache_object(misnamed_grid_hashes),
    "input grid hashes"
  )

  wrong_winner <- cache
  alternate <- wrong_winner$grid_summary$gic[1L, , drop = FALSE]
  wrong_winner$tuning$selected[1L, ] <- helpers$build_selected_row(
    list(selected_row = transform(
      alternate,
      score = (
        proportional_fuzzy_balanced_accuracy +
          correlation_fuzzy_balanced_accuracy
      ) / 2
    ))
  )
  wrong_winner$fixed_settings[wrong_winner$fixed_settings$IC == "GIC", ] <-
    helpers$build_fixed_rows(list(selected_row = alternate))
  testthat::expect_error(
    helpers$validate_cache_object(wrong_winner),
    "selected.*ranking"
  )

  infeasible_gic <- gic_wrapper
  infeasible_gic$result$summary_table$null_mean_false_positive_rate <- 0.20
  saveRDS(infeasible_gic, gic_path)
  testthat::expect_error(
    suppressWarnings(helpers$build_simulation_cache(gic_path, bic_path)),
    "No feasible GIC setting"
  )
  saveRDS(gic_wrapper, gic_path)

  mismatched_environment <- bic_wrapper
  mismatched_environment$provenance$R <- "different R"
  saveRDS(mismatched_environment, bic_path)
  testthat::expect_error(
    helpers$build_simulation_cache(gic_path, bic_path),
    "same R and platform"
  )
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

test_that("simulation vignette sources use manuscript-aligned reporting", {
  part1_path <- testthat::test_path("../../vignettes/simulation-study-part-1.Rmd")
  part2_path <- testthat::test_path("../../vignettes/simulation-study-part-2.Rmd")
  testthat::skip_if_not(
    all(file.exists(c(part1_path, part2_path))),
    "Simulation vignettes not present"
  )

  part1 <- paste(readLines(part1_path, warn = FALSE), collapse = "\n")
  part2 <- paste(readLines(part2_path, warn = FALSE), collapse = "\n")

  for (source in list(part1, part2)) {
    testthat::expect_match(
      source,
      "identical(preview_tables$schema_version, 3L)",
      fixed = TRUE
    )
    testthat::expect_match(source, "rownames(x) <- NULL", fixed = TRUE)
  }
  testthat::expect_false(grepl(
    'primary_metric = "weighted_fuzzy_f1"',
    part2,
    fixed = TRUE
  ))
  testthat::expect_gte(
    lengths(regmatches(
      part2,
      gregexpr(
        'primary_metric = "fuzzy_balanced_accuracy"',
        part2,
        fixed = TRUE
      )
    )),
    4L
  )
  testthat::expect_gte(
    lengths(regmatches(
      part2,
      gregexpr(
        "scenario_weights = c(proportional = 0.50, correlation = 0.50)",
        part2,
        fixed = TRUE
      )
    )),
    4L
  )
  for (hard_stop in c(
    "!gic_preview_tuned$used_all_settings",
    "!bic_preview_tuned$used_all_settings",
    "stopifnot(!gic_tuned$used_all_settings)",
    "stopifnot(!bic_tuned$used_all_settings)"
  )) {
    testthat::expect_match(part2, hard_stop, fixed = TRUE)
  }

  testthat::expect_match(part1, "fixed_null_display", fixed = TRUE)
  testthat::expect_match(part1, "fixed_recovery_display", fixed = TRUE)
  for (label in c(
    "Fuzzy recall", "Fuzzy specificity", "Fuzzy F1",
    "Fuzzy balanced accuracy"
  )) {
    testthat::expect_match(part1, label, fixed = TRUE)
  }
  for (compact_label in c("Fuzzy rec.", "Fuzzy spec.", "Fuzzy BA")) {
    testthat::expect_match(part1, compact_label, fixed = TRUE)
  }
  testthat::expect_match(part1, "class imbalance", ignore.case = TRUE)
  testthat::expect_match(part1, "near misses", ignore.case = TRUE)
  testthat::expect_match(part1, "complements", ignore.case = TRUE)

  testthat::expect_match(part2, "add_combined_score", fixed = TRUE)
  testthat::expect_match(part2, "rowMeans", fixed = TRUE)
  testthat::expect_match(
    part2,
    "Prop. Fuzzy balanced accuracy",
    fixed = TRUE
  )
  testthat::expect_match(
    part2,
    "Integration Fuzzy balanced accuracy",
    fixed = TRUE
  )
  testthat::expect_match(part2, "Score", fixed = TRUE)
  testthat::expect_match(part2, "candidate-node universe", fixed = TRUE)
  for (safeguard in c(
    "max_false_positive_rate", "max_any_false_positive",
    "min_evaluable_fraction"
  )) {
    testthat::expect_match(part2, safeguard, fixed = TRUE)
  }
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

test_that("experimental draft preserves calibrated original-generator behavior", {
  part1_path <- testthat::test_path(
    "../../vignettes/simulation-study-part-1.Rmd"
  )
  draft_path <- testthat::test_path(
    "../../experimental/simulation-study-manuscript-generator-calibration.Rmd"
  )
  testthat::skip_if_not(file.exists(part1_path), "Simulation vignette not present")
  testthat::skip_if_not(
    file.exists(draft_path),
    "Experimental calibration draft is unavailable in installed tests"
  )

  source <- paste(readLines(
    draft_path,
    warn = FALSE
  ), collapse = "\n")

  testthat::expect_match(
    source,
    "Draft: Calibrating the Original Manuscript Generator",
    fixed = TRUE
  )
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
