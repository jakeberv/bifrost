paired_tuning_policy <- function() {
  list(
    max_false_positive_rate = 0.10,
    max_any_false_positive = 0.05,
    min_evaluable_fraction = 0.50,
    primary_metric = "fuzzy_balanced_accuracy",
    scenario_weights = c(proportional = 0.50, correlation = 0.50),
    tie_break = "conservative",
    allow_infeasible = FALSE
  )
}

approved_paired_tuning_identity <- function() {
  list(
    package_commit = "db18184ddb06a5019123647ce218f9b717f76e49",
    design_fingerprint =
      "c6327e5bbc59b6d5d948651d3142068ce3d732848e52851957d4fc0844fc518e",
    base_search_options = list(
      formula = "trait_data ~ 1",
      method = "LL",
      error = FALSE
    ),
    simulation_generators = c(
      null = "empirical",
      proportional = "empirical",
      correlation = "empirical"
    )
  )
}

paired_tuning_input_paths <- function(source_dir) {
  stats::setNames(
    file.path(source_dir, c(
      "gic-tuning.rds", "bic-tuning.rds", "suite-summary.rds",
      "generation-accounting.csv"
    )),
    c("gic", "bic", "suite", "accounting")
  )
}

sha256_file <- function(path) {
  if (!requireNamespace("digest", quietly = TRUE)) {
    stop("Exporting the paired tuning cache requires the digest package.")
  }
  digest::digest(file = path, algo = "sha256", serialize = FALSE)
}

required_paired_summary_columns <- function() {
  metric_suffixes <- c(
    "precision", "recall", "f1", "specificity", "fpr",
    "balanced_accuracy"
  )
  c(
    "setting_id", "IC", "shift_acceptance_threshold", "min_descendant_tips",
    "null_seed", "null_completion_rate", "null_failure_rate",
    "null_mean_inferred_shifts", "null_evaluable_fraction",
    "null_mean_false_positive_rate", "null_fraction_any_false_positive",
    unlist(lapply(c("proportional", "correlation"), function(scenario) c(
      paste0(scenario, "_seed"), paste0(scenario, "_completion_rate"),
      paste0(scenario, "_failure_rate"),
      paste0(scenario, "_mean_inferred_shifts"),
      paste0(scenario, "_evaluable_fraction"),
      paste0(scenario, "_fuzzy_", metric_suffixes),
      paste0(scenario, "_strict_", metric_suffixes),
      paste0(scenario, "_weighted_fuzzy_f1")
    )), use.names = FALSE),
    unlist(lapply(c("null", "proportional", "correlation"), function(scenario) {
      paste0(scenario, c(
        "_n_attempted", "_n_generated", "_n_placement_failed",
        "_n_search_attempted", "_n_evaluated"
      ))
    }), use.names = FALSE)
  )
}

validate_generation_accounting <- function(accounting) {
  expected <- data.frame(
    group_id = c(
      "null-n250", "proportional-n250-s5", "integration-rate-n250-s5"
    ),
    scenario = c("null", "proportional", "integration-rate"),
    n_attempted = c(500L, 500L, 500L),
    n_generated = c(500L, 500L, 500L),
    n_placement_failed = c(0L, 0L, 0L),
    stringsAsFactors = FALSE
  )
  rownames(accounting) <- NULL
  if (!identical(accounting, expected)) {
    stop("Source generation accounting does not match the completed campaign.")
  }
  invisible(accounting)
}

validate_paired_sources <- function(wrappers, suite, accounting) {
  identity <- approved_paired_tuning_identity()
  validate_generation_accounting(accounting)
  if (!identical(suite$group_counts, accounting) ||
      !isTRUE(suite$all_attempts_accounted_for) ||
      !identical(nrow(suite$attempts), 1500L) ||
      !identical(sum(suite$attempts$status == "generated"), 1500L) ||
      !identical(sum(suite$attempts$status == "placement_failed"), 0L)) {
    stop("Suite generation accounting contradicts the accounting export.")
  }
  if (!identical(suite$selection_policy, paired_tuning_policy())) {
    stop("Source selection policy does not match the approved policy.")
  }
  if (!is.character(suite$design_fingerprint) ||
      length(suite$design_fingerprint) != 1L ||
      !identical(suite$design_fingerprint, identity$design_fingerprint)) {
    stop("Source design fingerprint is missing or does not match the approved campaign.")
  }

  commits <- vapply(wrappers, function(x) {
    value <- x$provenance$package_remote_sha
    if (!is.character(value) || length(value) != 1L || is.na(value)) "" else value
  }, character(1L))
  if (!all(commits == identity$package_commit)) {
    stop("Source package commit is missing or does not match the approved campaign.")
  }

  fingerprints <- vapply(wrappers, function(x) x$provenance$source_fingerprint,
    character(1L)
  )
  designs <- vapply(wrappers, function(x) x$grid$paired_settings, logical(1L))
  if (!all(designs) || length(unique(fingerprints)) != 1L ||
      !identical(unname(fingerprints[[1L]]), suite$source_fingerprint) ||
      !identical(wrappers$gic$provenance, wrappers$bic$provenance) ||
      !identical(wrappers$gic$grid$study_seeds, wrappers$bic$grid$study_seeds)) {
    stop("GIC and BIC require matching paired design and provenance.")
  }

  required <- required_paired_summary_columns()
  expected_settings <- expand.grid(
    shift_acceptance_threshold = c(10, 20, 30),
    min_descendant_tips = c(10L, 20L),
    KEEP.OUT.ATTRS = FALSE
  )
  for (key in names(wrappers)) {
    wrapper <- wrappers[[key]]
    ic <- toupper(key)
    grid <- wrapper$grid
    if (!identical(wrapper$schema_version, 1L) ||
        !inherits(grid, "bifrost_search_tuning_grid") ||
        !identical(grid$IC, ic) || !identical(wrapper$policy, paired_tuning_policy()) ||
        !is.data.frame(grid$summary_table) || nrow(grid$summary_table) != 6L ||
        !all(required %in% names(grid$summary_table)) ||
        !identical(
          grid$grid[, c("shift_acceptance_threshold", "min_descendant_tips")],
          expected_settings
        ) ||
        !identical(grid$generation_accounting, accounting) ||
        !identical(grid$attempts_by_scenario, c(
          null = 500L, proportional = 500L, correlation = 500L
        )) || !identical(grid$replicates_by_scenario, c(
          null = 500L, proportional = 500L, correlation = 500L
        ))) {
      stop(ic, " source summary has incomplete paired design or accounting.")
    }
    summary <- grid$summary_table
    if (!identical(grid$base_search_options, identity$base_search_options)) {
      stop(ic, " base search options do not match the approved campaign.")
    }
    if (!identical(grid$simulation_generators, identity$simulation_generators)) {
      stop(ic, " simulation generators do not match the approved campaign.")
    }
    numeric_required <- setdiff(required, "IC")
    if (any(!is.finite(unlist(summary[, numeric_required], use.names = FALSE)))) {
      stop(ic, " source summary has non-finite required metrics.")
    }
    for (scenario in c("null", "proportional", "correlation")) {
      completion <- summary[[paste0(scenario, "_completion_rate")]]
      failure <- summary[[paste0(scenario, "_failure_rate")]]
      if (!isTRUE(all.equal(failure, 1 - completion, tolerance = 1e-12))) {
        stop(ic, " completion and failure rates are inconsistent.")
      }
      if (any(completion != 1) || any(failure != 0)) {
        stop(ic, " completion rates contradict the completed campaign's 500/500 search counts.")
      }
    }
    expected_counts <- list(null = c(500L, 500L, 0L, 500L, 500L),
      proportional = c(500L, 500L, 0L, 500L, 500L),
      correlation = c(500L, 500L, 0L, 500L, 500L))
    for (scenario in names(expected_counts)) {
      columns <- paste0(scenario, c(
        "_n_attempted", "_n_generated", "_n_placement_failed",
        "_n_search_attempted", "_n_evaluated"
      ))
      actual <- unique(summary[, columns, drop = FALSE])
      if (nrow(actual) != 1L || !identical(
        as.integer(actual[1L, ]), expected_counts[[scenario]]
      )) stop(ic, " source summary has contradictory generation accounting.")
    }

    metric_suffixes <- c(
      "fuzzy_precision", "fuzzy_recall", "fuzzy_f1", "fuzzy_specificity",
      "fuzzy_fpr", "fuzzy_balanced_accuracy", "strict_precision",
      "strict_recall", "strict_f1", "strict_specificity", "strict_fpr",
      "strict_balanced_accuracy", "weighted_fuzzy_f1"
    )
    for (scenario in c("null", "proportional", "integration-rate")) {
      prefix <- if (identical(scenario, "integration-rate")) "correlation" else scenario
      suite_rows <- suite$summary_table[
        suite$summary_table$IC == ic & suite$summary_table$scenario == scenario,
        , drop = FALSE
      ]
      suite_rows <- suite_rows[order(suite_rows$setting_id), , drop = FALSE]
      comparison <- list(
        setting_id = summary$setting_id,
        IC = summary$IC,
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
        n_evaluated = summary[[paste0(prefix, "_n_evaluated")]]
      )
      if (identical(scenario, "null")) {
        comparison$mean_false_positive_rate <- summary$null_mean_false_positive_rate
        comparison$fraction_any_false_positive <-
          summary$null_fraction_any_false_positive
      } else {
        for (metric in metric_suffixes) {
          comparison[[metric]] <- summary[[paste0(prefix, "_", metric)]]
        }
      }
      comparison <- as.data.frame(comparison, check.names = FALSE)
      if (nrow(suite_rows) != 6L || !all(names(comparison) %in% names(suite_rows)) ||
          !isTRUE(all.equal(
            suite_rows[, names(comparison), drop = FALSE], comparison,
            tolerance = 1e-12, check.attributes = FALSE
          ))) {
        stop(ic, " grid summary does not correspond to the suite summary.")
      }
    }
  }

  total <- sum(vapply(wrappers, function(wrapper) {
    summary <- wrapper$grid$summary_table
    sum(summary$null_n_evaluated + summary$proportional_n_evaluated +
      summary$correlation_n_evaluated)
  }, integer(1L)))
  config_keys <- unique(paste(
    suite$summary_table$IC,
    suite$summary_table$shift_acceptance_threshold,
    suite$summary_table$min_descendant_tips,
    sep = "/"
  ))
  if (!identical(total, 18000L) || nrow(suite$summary_table) != 36L ||
      length(config_keys) != 12L) {
    stop("Source summaries must represent 18,000 searches and twelve configurations.")
  }
  invisible(NULL)
}

derive_planted_clade_counts <- function(suite) {
  required <- c("scenario", "IC", "setting_id", "n_true_shifts_10_19",
    "n_true_shifts_20_40")
  if (!is.data.frame(suite$summary_table) ||
      !all(required %in% names(suite$summary_table))) {
    stop("Suite summary is missing planted-clade counts.")
  }
  counts <- suite$summary_table[, required, drop = FALSE]
  count_columns <- c("n_true_shifts_10_19", "n_true_shifts_20_40")
  if (anyNA(counts[, count_columns, drop = FALSE]) ||
      any(!is.finite(unlist(counts[, count_columns, drop = FALSE]))) ||
      any(unlist(counts[, count_columns, drop = FALSE]) < 0)) {
    stop("Suite planted-clade counts must be finite, non-negative, and complete.")
  }
  scenarios <- c("null", "proportional", "integration-rate")
  if (!setequal(unique(counts$scenario), scenarios)) {
    stop("Suite planted-clade counts have unexpected scenarios.")
  }
  out <- do.call(rbind, lapply(scenarios, function(scenario) {
    rows <- counts[counts$scenario == scenario, count_columns, drop = FALSE]
    unique_rows <- unique(rows)
    if (nrow(unique_rows) != 1L) {
      stop("Suite planted-clade counts are inconsistent across IC/settings for ",
        scenario, ".")
    }
    data.frame(scenario = scenario, unique_rows, row.names = NULL,
      check.names = FALSE)
  }))
  rownames(out) <- NULL
  out[, count_columns] <- lapply(out[, count_columns, drop = FALSE], as.integer)
  shifted <- out$scenario != "null"
  if (any(rowSums(out[shifted, count_columns, drop = FALSE]) != 2500L) ||
      any(unlist(out[!shifted, count_columns, drop = FALSE]) != 0L)) {
    stop("Suite planted-clade counts contradict the completed campaign.")
  }
  out
}

compact_paired_grid <- function(grid) {
  out <- grid[c(
    "IC", "grid", "summary_table", "base_search_options", "null_replicates",
    "simulation_generators", "fuzzy_distance", "weighted", "paired_settings",
    "store_studies", "study_seeds", "replicates_by_scenario",
    "attempts_by_scenario", "generation_accounting", "accounting_scope"
  )]
  out$studies <- NULL
  out$store_studies <- FALSE
  class(out) <- c("bifrost_search_tuning_grid", "list")
  out
}

paired_tuning_display <- function(grid, selection) {
  summary <- grid$summary_table
  feasible <- selection$feasible_table
  feasible_key <- paste(feasible$shift_acceptance_threshold,
    feasible$min_descendant_tips, sep = "/")
  key <- paste(summary$shift_acceptance_threshold, summary$min_descendant_tips,
    sep = "/")
  selected <- selection$selected_row
  selected_key <- paste(selected$shift_acceptance_threshold,
    selected$min_descendant_tips, sep = "/")
  score <- 0.5 * summary$proportional_fuzzy_balanced_accuracy +
    0.5 * summary$correlation_fuzzy_balanced_accuracy
  data.frame(
    Threshold = summary$shift_acceptance_threshold,
    `Min clade` = summary$min_descendant_tips,
    `Null FP` = summary$null_mean_false_positive_rate,
    `Null any FP` = summary$null_fraction_any_false_positive,
    `Prop. Fuzzy balanced accuracy` =
      summary$proportional_fuzzy_balanced_accuracy,
    `Integration Fuzzy balanced accuracy` =
      summary$correlation_fuzzy_balanced_accuracy,
    Score = score,
    Status = ifelse(key == selected_key, "Selected",
      ifelse(key %in% feasible_key, "Eligible", "Excluded")),
    check.names = FALSE
  )
}

paired_selected_row <- function(selection) {
  row <- selection$selected_row
  data.frame(
    IC = row$IC,
    Threshold = row$shift_acceptance_threshold,
    `Min clade` = row$min_descendant_tips,
    `Null FP` = row$null_mean_false_positive_rate,
    `Null any FP` = row$null_fraction_any_false_positive,
    `Prop. Fuzzy balanced accuracy` = row$proportional_fuzzy_balanced_accuracy,
    `Integration Fuzzy balanced accuracy` =
      row$correlation_fuzzy_balanced_accuracy,
    Score = row$score,
    Status = "Selected",
    check.names = FALSE
  )
}

sanitize_tuning_provenance <- function(provenance) {
  provenance[c(
    "package_version", "package_remote_sha", "package_github_sha1",
    "data_md5", "source_fingerprint", "runner_sha256", "template_traits",
    "template_tips", "template_covariance_sha256", "template_residual_df",
    "platform", "R", "runner_files_sha256"
  )]
}

fixed_setting_rows <- function(grids) {
  metric_names <- c(
    "Strict precision" = "strict_precision",
    "Strict recall" = "strict_recall",
    "Strict F1" = "strict_f1",
    "Strict specificity" = "strict_specificity",
    "Strict FPR" = "strict_fpr",
    "Strict balanced accuracy" = "strict_balanced_accuracy",
    "Fuzzy precision" = "fuzzy_precision",
    "Fuzzy recall" = "fuzzy_recall",
    "Fuzzy F1" = "fuzzy_f1",
    "Fuzzy specificity" = "fuzzy_specificity",
    "Fuzzy FPR" = "fuzzy_fpr",
    "Fuzzy balanced accuracy" = "fuzzy_balanced_accuracy"
  )
  rows <- do.call(rbind, lapply(grids, function(grid) {
    summary <- grid$summary_table
    selected <- summary[
      summary$shift_acceptance_threshold == 10 &
        summary$min_descendant_tips == 10L,
      , drop = FALSE
    ]
    if (nrow(selected) != 1L) {
      stop(grid$IC, " grid must contain exactly one threshold-10/min-clade-10 setting.")
    }
    do.call(rbind, lapply(c("null", "proportional", "correlation"), function(scenario) {
      label <- switch(scenario, null = "Null", proportional = "Proportional",
        correlation = "Integration-rate")
      out <- data.frame(
        IC = grid$IC,
        Scenario = label,
        `Mean FP` = if (scenario == "null") selected$null_mean_false_positive_rate else NA_real_,
        `Any FP` = if (scenario == "null") selected$null_fraction_any_false_positive else NA_real_,
        check.names = FALSE
      )
      for (display in names(metric_names)) {
        out[[display]] <- if (scenario == "null") NA_real_ else
          selected[[paste0(scenario, "_", metric_names[[display]])]]
      }
      out$`Weighted fuzzy F1` <- if (scenario == "null") NA_real_ else
        selected[[paste0(scenario, "_weighted_fuzzy_f1")]]
      out$`Mean shifts` <- selected[[paste0(scenario, "_mean_inferred_shifts")]]
      out$Evaluable <- selected[[paste0(scenario, "_evaluable_fraction")]]
      out$Completion <- selected[[paste0(scenario, "_completion_rate")]]
      out$Failure <- selected[[paste0(scenario, "_failure_rate")]]
      out
    }))
  }))
  rownames(rows) <- NULL
  rows
}

fixed_setting_provenance <- function(grids, suite, source_hashes) {
  list(
    simulation_generator = "empirical",
    n_replicates_per_setting = 500L,
    tree_tip_count = 250L,
    n_true_shifts = 5L,
    min_shift_tips = 10L,
    max_shift_tips = 40L,
    shift_acceptance_threshold = 10,
    min_descendant_tips = 10L,
    seed = 5L,
    calibration_error = TRUE,
    search_error = FALSE,
    integration_power_range = c(0.5, 1.25),
    integration_exclude_range = c(0.8, 1.1),
    package_commit = approved_paired_tuning_identity()$package_commit,
    design_fingerprint = suite$design_fingerprint,
    source_sha256 = source_hashes,
    metric_accounting_version = "candidate-node-aware-v1",
    accounting_scope = suite$accounting_scope
  )
}

validate_paired_tuning_cache <- function(cache) {
  if (!is.list(cache) || !identical(cache$schema_version, 4L) ||
      !is.list(cache$provenance) ||
      !identical(names(cache$provenance), c("fixed_settings", "tuning")) ||
      !is.data.frame(cache$fixed_settings) ||
      !identical(names(cache$tuning_grids), c("gic", "bic")) ||
      !identical(names(cache$grid_summary), c("gic", "bic")) ||
      !identical(names(cache$tuning), c("gic", "bic", "selected"))) {
    stop("Simulation vignette cache must have the schema-4 paired-tuning structure.")
  }
  if (!identical(cache$provenance$tuning$policy, paired_tuning_policy()) ||
      !identical(cache$provenance$tuning$total_evaluated_searches, 18000L) ||
      !is.data.frame(cache$provenance$tuning$planted_clade_counts)) {
    stop("Schema-4 tuning provenance is incomplete.")
  }
  validate_generation_accounting(cache$provenance$tuning$generation_accounting)
  expected_fixed <- fixed_setting_rows(cache$tuning_grids)
  if (!isTRUE(all.equal(cache$fixed_settings, expected_fixed, tolerance = 1e-12,
    check.attributes = FALSE))) {
    stop("Schema-4 fixed settings disagree with the threshold-10/min-clade-10 grid rows.")
  }
  fixed_provenance <- cache$provenance$fixed_settings
  tuning_provenance <- cache$provenance$tuning
  expected_fixed_provenance <- list(
    simulation_generator = tuning_provenance$simulation_generator,
    n_replicates_per_setting = 500L,
    tree_tip_count = tuning_provenance$paired_design$tree_tip_count,
    n_true_shifts = tuning_provenance$paired_design$n_true_shifts,
    min_shift_tips = tuning_provenance$paired_design$min_shift_tips,
    max_shift_tips = tuning_provenance$paired_design$max_shift_tips,
    shift_acceptance_threshold = 10,
    min_descendant_tips = 10L,
    seed = tuning_provenance$seed,
    calibration_error = tuning_provenance$paired_design$calibration_error,
    search_error = tuning_provenance$paired_design$search_error,
    integration_power_range = tuning_provenance$integration_power_range,
    integration_exclude_range = tuning_provenance$integration_exclude_range,
    package_commit = tuning_provenance$package_commit,
    design_fingerprint = tuning_provenance$paired_design$design_fingerprint,
    source_sha256 = tuning_provenance$source_sha256,
    metric_accounting_version = tuning_provenance$metric_accounting_version,
    accounting_scope = tuning_provenance$accounting_scope
  )
  if (!identical(fixed_provenance, expected_fixed_provenance)) {
    stop("Schema-4 fixed-settings provenance disagrees with tuning provenance.")
  }
  for (key in c("gic", "bic")) {
    grid <- cache$tuning_grids[[key]]
    selection <- do.call(selectTunedSearchParameters,
      c(list(tuning_grid = grid), paired_tuning_policy()))
    expected <- paired_selected_row(selection)
    actual <- cache$tuning$selected[cache$tuning$selected$IC == toupper(key), ,
      drop = FALSE]
    rownames(actual) <- NULL
    if (!isTRUE(all.equal(actual, expected, tolerance = 1e-12,
      check.attributes = FALSE)) ||
        !identical(cache$grid_summary[[key]], grid$summary_table) ||
        !isTRUE(all.equal(cache$tuning[[key]], paired_tuning_display(grid, selection),
          tolerance = 1e-12, check.attributes = FALSE))) {
      stop("Schema-4 displayed or selected tuning values disagree with the compact grid.")
    }
  }
  serialized <- as.character(unlist(cache, recursive = TRUE, use.names = TRUE))
  if (any(grepl("(^|[= ])(/Users/|/home/|/private/tmp/|[A-Za-z]:[/\\\\])",
    serialized))) stop("Schema-4 cache contains a private machine path.")
  invisible(cache)
}

build_paired_tuning_cache <- function(source_dir) {
  paths <- paired_tuning_input_paths(source_dir)
  missing <- paths[!file.exists(paths)]
  if (length(missing)) stop("Missing paired-tuning inputs: ", paste(basename(missing), collapse = ", "))
  wrappers <- list(gic = readRDS(paths[["gic"]]), bic = readRDS(paths[["bic"]]))
  suite <- readRDS(paths[["suite"]])
  accounting <- utils::read.csv(paths[["accounting"]], stringsAsFactors = FALSE)
  validate_paired_sources(wrappers, suite, accounting)
  planted_clade_counts <- derive_planted_clade_counts(suite)
  grids <- lapply(wrappers, function(x) compact_paired_grid(x$grid))
  selections <- lapply(grids, function(grid) do.call(
    selectTunedSearchParameters,
    c(list(tuning_grid = grid), paired_tuning_policy())
  ))
  source_hashes <- vapply(paths, sha256_file, character(1L))
  tuning_provenance <- c(sanitize_tuning_provenance(wrappers$gic$provenance), list(
    paired_design = list(
      paired_settings = TRUE,
      seed_scheme = "runSearchTuningGrid-paired-v1",
      tree_tip_count = 250L,
      response_traits = 12L,
      n_true_shifts = 5L,
      min_shift_tips = 10L,
      max_shift_tips = 40L,
      buffer = 3L,
      proportional_scale_factor_range = c(0.1, 2.0),
      proportional_exclude_range = c(0.5, 1.5),
      integration_power_range = c(0.5, 1.25),
      integration_exclude_range = c(0.8, 1.1),
      calibration_error = TRUE,
      search_error = FALSE,
      study_seeds = grids$gic$study_seeds,
      design_fingerprint = suite$design_fingerprint,
      attempts_by_scenario = grids$gic$attempts_by_scenario,
      replicates_by_scenario = grids$gic$replicates_by_scenario
    ),
    package_commit = wrappers$gic$provenance$package_remote_sha,
    seed = 5L,
    source_sha256 = source_hashes,
    metric_accounting_version = "candidate-node-aware-v1",
    accounting_scope = suite$accounting_scope,
    generation_accounting = accounting,
    total_evaluated_searches = 18000L,
    planted_clade_counts = planted_clade_counts,
    simulation_generator = "empirical",
    integration_power_range = c(0.5, 1.25),
    integration_exclude_range = c(0.8, 1.1),
    policy = paired_tuning_policy()
  ))
  cache <- list(
    schema_version = 4L,
    provenance = list(
      fixed_settings = fixed_setting_provenance(grids, suite, source_hashes),
      tuning = tuning_provenance
    ),
    fixed_settings = fixed_setting_rows(grids),
    tuning_grids = grids,
    tuning = list(
      gic = paired_tuning_display(grids$gic, selections$gic),
      bic = paired_tuning_display(grids$bic, selections$bic),
      selected = do.call(rbind, lapply(selections, paired_selected_row))
    ),
    grid_summary = lapply(grids, `[[`, "summary_table")
  )
  rownames(cache$tuning$selected) <- NULL
  validate_paired_tuning_cache(cache)
  cache
}

write_paired_tuning_cache <- function(cache, output_path) {
  validate_paired_tuning_cache(cache)
  dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
  temporary <- tempfile(paste0(".", basename(output_path), "-"),
    tmpdir = dirname(output_path), fileext = ".tmp")
  on.exit(unlink(temporary), add = TRUE)
  saveRDS(cache, temporary, compress = "xz", version = 3)
  validate_paired_tuning_cache(readRDS(temporary))
  if (!file.rename(temporary, output_path)) stop("Could not atomically replace ", output_path)
  invisible(output_path)
}

update_simulation_tuning_cache_main <- function(args = commandArgs(trailingOnly = TRUE)) {
  if (length(args) != 2L) stop(paste(
    "Usage: Rscript data-raw/update_simulation_tuning_cache.R",
    "<source-summary-dir> <schema-4-cache.rds>"
  ))
  if (!requireNamespace("pkgload", quietly = TRUE)) stop("Exporter requires pkgload.")
  pkgload::load_all(".", quiet = TRUE, export_all = FALSE, helpers = FALSE)
  cache <- build_paired_tuning_cache(args[[1L]])
  write_paired_tuning_cache(cache, args[[2L]])
}

if (sys.nframe() == 0L) update_simulation_tuning_cache_main()
