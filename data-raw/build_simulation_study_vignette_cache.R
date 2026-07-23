ensure_grid_runner_helpers <- function(root = ".") {
  required <- c(
    "simulation_grid_design",
    "validate_grid_wrapper",
    "validate_grid_summary_table",
    "normalize_source_provenance",
    "required_grid_summary_columns",
    "grid_recovery_metric_suffixes"
  )
  target_env <- environment(ensure_grid_runner_helpers)
  missing_helpers <- required[!vapply(
    required,
    exists,
    logical(1L),
    envir = target_env,
    inherits = TRUE
  )]
  if (length(missing_helpers) > 0L) {
    runner_path <- file.path(
      root,
      "data-raw",
      "run_simulation_study_vignette_grids.R"
    )
    if (!file.exists(runner_path)) {
      stop("Grid runner helpers not found: ", runner_path)
    }
    sys.source(runner_path, envir = target_env)
  }
  invisible(NULL)
}

wrapper_source_provenance <- function(wrapper) {
  normalize_source_provenance(wrapper$provenance)
}

read_grid_run <- function(path, expected_ic) {
  ensure_grid_runner_helpers()
  if (!file.exists(path)) {
    stop("Grid result not found: ", path)
  }
  wrapper <- readRDS(path)
  source_provenance <- wrapper_source_provenance(wrapper)
  validate_grid_wrapper(
    wrapper,
    expected_ic = expected_ic,
    expected_design = simulation_grid_design("full"),
    expected_source_provenance = source_provenance
  )
  list(
    wrapper = wrapper,
    grid = wrapper$result,
    path = path,
    source_provenance = source_provenance
  )
}

validate_matching_grid_provenance <- function(gic, bic) {
  if (!identical(gic$wrapper$provenance$R, bic$wrapper$provenance$R) ||
      !identical(
        gic$wrapper$provenance$platform,
        bic$wrapper$provenance$platform
      )) {
    stop("GIC and BIC source runs must use the same R and platform environment.")
  }
  if (!matching_source_provenance(
    gic$source_provenance,
    bic$source_provenance
  )) {
    stop("GIC and BIC source runs must use identical generator provenance.")
  }
  invisible(NULL)
}

select_grid_parameters <- function(grid_run) {
  selection <- selectTunedSearchParameters(
    grid_run$grid,
    max_false_positive_rate = 0.10,
    max_any_false_positive = 0.50,
    min_evaluable_fraction = 0.50,
    primary_metric = "fuzzy_balanced_accuracy",
    tie_break = "conservative",
    allow_infeasible = TRUE
  )
  if (isTRUE(selection$used_all_settings)) {
    stop(
      "No feasible ",
      grid_run$grid$IC,
      " setting passed the mandatory null and evaluability safeguards."
    )
  }
  selection
}

scenario_values <- function(row, metric_suffix) {
  c(
    NA_real_,
    unname(row[[paste0("proportional_", metric_suffix)]]),
    unname(row[[paste0("correlation_", metric_suffix)]])
  )
}

fixed_metric_columns <- function() {
  c(
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
}

build_fixed_rows <- function(selection) {
  row <- selection$selected_row
  metrics <- fixed_metric_columns()
  metric_values <- lapply(unname(metrics), function(metric_suffix) {
    scenario_values(row, metric_suffix)
  })
  names(metric_values) <- names(metrics)
  data.frame(
    list(
      IC = rep(row$IC, 3L),
      Scenario = c("Null", "Proportional", "Integration-rate"),
      `Mean FP` = c(row$null_mean_false_positive_rate, NA_real_, NA_real_),
      `Any FP` = c(row$null_fraction_any_false_positive, NA_real_, NA_real_)
    ),
    metric_values,
    list(
      `Weighted fuzzy F1` = c(
        NA_real_,
        row$proportional_weighted_fuzzy_f1,
        row$correlation_weighted_fuzzy_f1
      ),
      `Mean shifts` = c(
        row$null_mean_inferred_shifts,
        row$proportional_mean_inferred_shifts,
        row$correlation_mean_inferred_shifts
      ),
      Evaluable = c(
        row$null_evaluable_fraction,
        row$proportional_evaluable_fraction,
        row$correlation_evaluable_fraction
      ),
      Completion = c(
        row$null_completion_rate,
        row$proportional_completion_rate,
        row$correlation_completion_rate
      ),
      Failure = c(
        row$null_failure_rate,
        row$proportional_failure_rate,
        row$correlation_failure_rate
      )
    ),
    check.names = FALSE
  )
}

build_tuning_table <- function(grid_run) {
  summary <- grid_run$grid$summary_table
  data.frame(
    Threshold = summary$shift_acceptance_threshold,
    `Min clade` = summary$min_descendant_tips,
    `Null FP` = summary$null_mean_false_positive_rate,
    `Prop. Fuzzy balanced accuracy` =
      summary$proportional_fuzzy_balanced_accuracy,
    `Integration Fuzzy balanced accuracy` =
      summary$correlation_fuzzy_balanced_accuracy,
    check.names = FALSE
  )
}

build_selected_row <- function(selection) {
  row <- selection$selected_row
  data.frame(
    IC = row$IC,
    Threshold = row$shift_acceptance_threshold,
    `Min clade` = row$min_descendant_tips,
    `Null FP` = row$null_mean_false_positive_rate,
    `Prop. Fuzzy balanced accuracy` =
      row$proportional_fuzzy_balanced_accuracy,
    `Integration Fuzzy balanced accuracy` =
      row$correlation_fuzzy_balanced_accuracy,
    Score = row$score,
    check.names = FALSE
  )
}

approved_generator_source_files <- function() {
  c(
    "data-raw/run_simulation_study_vignette_grids.R",
    "DESCRIPTION",
    "NAMESPACE",
    "R/bifrost_search-methods.R",
    "R/formula-normalization.R",
    "R/icTrajectory.R",
    "R/lineage_rates.R",
    "R/plotting-utils.R",
    "R/rate-map-plot.R",
    "R/rate-map.R",
    "R/regime-integration.R",
    "R/searchOptimalConfiguration-helpers.R",
    "R/searchOptimalConfiguration.R",
    "R/shift-distributions.R",
    "R/simulation-generators.R",
    "R/simulation-helpers.R",
    "R/simulation-methods.R",
    "R/simulation-studies.R",
    "R/simulation-template.R",
    "R/simulation-tuning.R",
    "R/utils.R",
    "inst/extdata/avian-skeleton/passerine_bodyplan_tree.tre",
    "inst/extdata/avian-skeleton/passerine_bodyplan_data.RDS"
  )
}

validate_cache_provenance <- function(provenance) {
  if (!is.list(provenance)) {
    stop("Cache provenance must be a list.")
  }

  design_contract <- list(
    simulation_generator = "empirical",
    scenario_interpretation = paste(
      "Empirically calibrated spectral integration change with reciprocal",
      "marginal-variance scaling"
    ),
    original_generator_available = TRUE,
    n_replicates_per_setting = 100L,
    tree_tip_count = 250L,
    n_true_shifts = 5L,
    min_shift_tips = 10L,
    max_shift_tips = 20L,
    integration_power_range = c(0.5, 1.25),
    integration_exclude_range = c(0.8, 1.1),
    fuzzy_distance = 2L
  )
  design_matches <- vapply(names(design_contract), function(name) {
    identical(provenance[[name]], design_contract[[name]])
  }, logical(1L))
  if (!all(design_matches)) {
    stop("Cache provenance does not match the approved design contract.")
  }

  selection_policy <- list(
    primary_metric = "fuzzy_balanced_accuracy",
    scenario_weights = c(proportional = 0.5, correlation = 0.5),
    max_false_positive_rate = 0.10,
    max_any_false_positive = 0.50,
    min_evaluable_fraction = 0.50,
    tie_break = "conservative"
  )
  if (!identical(provenance$selection, selection_policy)) {
    stop("Cache provenance does not match the approved selection policy.")
  }

  nonempty_string <- function(x) {
    is.character(x) && length(x) == 1L && !is.na(x) && nzchar(trimws(x))
  }
  if (!nonempty_string(provenance$R) ||
      !nonempty_string(provenance$platform)) {
    stop("Cache provenance must record nonempty R and platform strings.")
  }

  generator_provenance <- tryCatch(
    normalize_source_provenance(list(
      commit = provenance$generator_commit,
      source_files = names(provenance$generator_source_md5),
      source_md5 = provenance$generator_source_md5,
      source_dirty = provenance$generator_source_dirty,
      source_status = provenance$generator_source_status
    )),
    error = identity
  )
  valid_commit <- !inherits(generator_provenance, "error") &&
    grepl("^[0-9a-f]{40}$", generator_provenance$commit)
  expected_source_files <- approved_generator_source_files()
  valid_inventory <- !inherits(generator_provenance, "error") &&
    identical(generator_provenance$source_files, expected_source_files)
  if (!valid_commit || !valid_inventory) {
    detail <- if (inherits(generator_provenance, "error")) {
      conditionMessage(generator_provenance)
    } else if (!valid_commit) {
      "commit must be a 40-character lowercase Git object ID"
    } else {
      "source inventory does not match the calibration generator"
    }
    stop("Cache generator source provenance is invalid: ", detail, ".")
  }

  input_hashes <- provenance$source_md5
  if (!is.character(input_hashes) ||
      !identical(names(input_hashes), c("GIC", "BIC")) ||
      anyNA(input_hashes) ||
      any(!grepl("^[0-9a-f]{32}$", input_hashes))) {
    stop("Cache input grid hashes must be named GIC/BIC lowercase MD5 values.")
  }
  invisible(provenance)
}

validate_cache_object <- function(cache) {
  ensure_grid_runner_helpers()
  if (!is.list(cache) || !identical(cache$schema_version, 3L)) {
    stop("Simulation vignette cache must use schema 3.")
  }
  if (!is.list(cache$provenance) ||
      !identical(
        cache$provenance$simulation_generator,
        "empirical"
      ) ||
      !identical(
        cache$provenance$selection$primary_metric,
        "fuzzy_balanced_accuracy"
      )) {
    stop("Cache provenance is incomplete or uses the wrong selection metric.")
  }
  validate_cache_provenance(cache$provenance)

  fixed <- cache$fixed_settings
  fixed_metrics <- fixed_metric_columns()
  required_fixed <- c(
    "IC", "Scenario", "Mean FP", "Any FP", names(fixed_metrics),
    "Weighted fuzzy F1", "Mean shifts", "Evaluable", "Completion", "Failure"
  )
  if (!is.data.frame(fixed) || nrow(fixed) != 6L ||
      !all(required_fixed %in% names(fixed)) ||
      !identical(as.character(fixed$IC), rep(c("GIC", "BIC"), each = 3L)) ||
      !identical(
        as.character(fixed$Scenario),
        rep(c("Null", "Proportional", "Integration-rate"), 2L)
      )) {
    stop("Cache fixed-settings table has an invalid structure.")
  }
  shifted <- fixed$Scenario != "Null"
  shifted_values <- unlist(
    fixed[shifted, names(fixed_metrics), drop = FALSE],
    use.names = FALSE
  )
  if (!all(is.finite(shifted_values))) {
    stop("Cache fixed-setting recovery metrics must be finite.")
  }

  tuning_names <- c(
    "Threshold", "Min clade", "Null FP",
    "Prop. Fuzzy balanced accuracy",
    "Integration Fuzzy balanced accuracy"
  )
  if (!is.list(cache$tuning) ||
      !identical(names(cache$tuning), c("gic", "bic", "selected"))) {
    stop("Cache tuning tables have an invalid structure.")
  }
  for (key in c("gic", "bic")) {
    table <- cache$tuning[[key]]
    if (!is.data.frame(table) || nrow(table) != 6L ||
        !identical(names(table), tuning_names) ||
        !all(is.finite(unlist(table, use.names = FALSE)))) {
      stop("Cache tuning table ", key, " is incomplete or non-finite.")
    }
  }

  selected <- cache$tuning$selected
  selected_names <- c("IC", tuning_names, "Score")
  if (!is.data.frame(selected) || nrow(selected) != 2L ||
      !identical(names(selected), selected_names) ||
      !identical(as.character(selected$IC), c("GIC", "BIC"))) {
    stop("Cache selected tuning table has an invalid structure.")
  }
  selected_numeric <- selected[
    , setdiff(selected_names, "IC"),
    drop = FALSE
  ]
  if (!all(is.finite(unlist(selected_numeric, use.names = FALSE)))) {
    stop("Cache must contain finite selected metrics and scores.")
  }
  weights <- cache$provenance$selection$scenario_weights
  if (!is.numeric(weights) ||
      !identical(names(weights), c("proportional", "correlation")) ||
      !identical(unname(weights), c(0.5, 0.5))) {
    stop("Cache selection weights are invalid.")
  }
  weights <- weights / sum(weights)
  expected_scores <-
    weights[[1L]] * selected$`Prop. Fuzzy balanced accuracy` +
    weights[[2L]] * selected$`Integration Fuzzy balanced accuracy`
  if (!isTRUE(all.equal(selected$Score, expected_scores, tolerance = 1e-12))) {
    stop("Cache selected scores are inconsistent with selected metrics.")
  }

  if (!is.list(cache$grid_summary) ||
      !identical(names(cache$grid_summary), c("gic", "bic"))) {
    stop("Cache raw grid summaries have an invalid structure.")
  }
  design <- simulation_grid_design("full")
  expected_settings <- expand.grid(
    shift_acceptance_threshold = design$shift_acceptance_thresholds,
    min_descendant_tips = design$min_descendant_tips_values,
    KEEP.OUT.ATTRS = FALSE
  )
  for (i in seq_len(nrow(selected))) {
    key <- tolower(selected$IC[[i]])
    raw <- cache$grid_summary[[key]]
    validate_grid_summary_table(raw, selected$IC[[i]], design)
    expected_tuning <- build_tuning_table(list(
      grid = list(summary_table = raw)
    ))
    actual_tuning <- cache$tuning[[key]]
    rownames(expected_tuning) <- NULL
    rownames(actual_tuning) <- NULL
    if (!isTRUE(all.equal(
      actual_tuning,
      expected_tuning,
      tolerance = 1e-12,
      check.attributes = FALSE
    ))) {
      stop("Cache tuning table does not match its raw grid summary.")
    }
    if (!is.data.frame(raw) || nrow(raw) != 6L ||
        !all(required_grid_summary_columns() %in% names(raw)) ||
        !identical(
          raw[, names(expected_settings), drop = FALSE],
          expected_settings
        )) {
      stop("Cache raw grid summary ", key, " is incomplete.")
    }
    tuning_grid <- list(
      IC = selected$IC[[i]],
      summary_table = raw,
      base_search_options = list()
    )
    class(tuning_grid) <- c("bifrost_search_tuning_grid", class(tuning_grid))
    expected_selection <- selectTunedSearchParameters(
      tuning_grid,
      max_false_positive_rate =
        cache$provenance$selection$max_false_positive_rate,
      max_any_false_positive =
        cache$provenance$selection$max_any_false_positive,
      min_evaluable_fraction =
        cache$provenance$selection$min_evaluable_fraction,
      primary_metric = cache$provenance$selection$primary_metric,
      scenario_weights = weights,
      tie_break = cache$provenance$selection$tie_break
    )
    if (isTRUE(expected_selection$used_all_settings)) {
      stop(
        "Cache raw grid has no feasible ",
        selected$IC[[i]],
        " setting under its mandatory safeguards."
      )
    }
    expected_selected <- build_selected_row(expected_selection)
    actual_selected <- selected[i, , drop = FALSE]
    rownames(expected_selected) <- NULL
    rownames(actual_selected) <- NULL
    if (!isTRUE(all.equal(
      actual_selected,
      expected_selected,
      tolerance = 1e-12,
      check.attributes = FALSE
    ))) {
      stop("Cache selected row does not match independent safeguard ranking.")
    }
    selected_index <-
      raw$shift_acceptance_threshold == selected$Threshold[[i]] &
      raw$min_descendant_tips == selected$`Min clade`[[i]]
    if (sum(selected_index) != 1L) {
      stop("Cache selected row is absent from its raw grid summary.")
    }
    raw_selected <- raw[selected_index, , drop = FALSE]
    if (!isTRUE(all.equal(
      selected$`Prop. Fuzzy balanced accuracy`[[i]],
      raw_selected$proportional_fuzzy_balanced_accuracy[[1L]],
      tolerance = 1e-12
    )) || !isTRUE(all.equal(
      selected$`Integration Fuzzy balanced accuracy`[[i]],
      raw_selected$correlation_fuzzy_balanced_accuracy[[1L]],
      tolerance = 1e-12
    ))) {
      stop("Cache selected metrics are inconsistent with raw grid summaries.")
    }

    fixed_ic <- fixed[fixed$IC == selected$IC[[i]], , drop = FALSE]
    expected_fixed <- build_fixed_rows(expected_selection)
    rownames(fixed_ic) <- NULL
    rownames(expected_fixed) <- NULL
    if (!isTRUE(all.equal(
      fixed_ic,
      expected_fixed,
      tolerance = 1e-12,
      check.attributes = FALSE
    ))) {
      stop("Cache selected and fixed rows do not match safeguard ranking.")
    }
    for (scenario in c("Proportional", "Integration-rate")) {
      prefix <- if (identical(scenario, "Proportional")) {
        "proportional_"
      } else {
        "correlation_"
      }
      fixed_row <- fixed_ic[fixed_ic$Scenario == scenario, , drop = FALSE]
      expected_fixed <- vapply(unname(fixed_metrics), function(suffix) {
        raw_selected[[paste0(prefix, suffix)]][[1L]]
      }, numeric(1L))
      actual_fixed <- unname(unlist(
        fixed_row[, names(fixed_metrics), drop = FALSE],
        use.names = FALSE
      ))
      if (!isTRUE(all.equal(
        actual_fixed,
        unname(expected_fixed),
        tolerance = 1e-12
      ))) {
        stop("Cache selected and fixed recovery metrics are inconsistent.")
      }
    }
  }
  invisible(cache)
}

build_simulation_cache <- function(gic_path, bic_path) {
  ensure_grid_runner_helpers()
  gic <- read_grid_run(gic_path, "GIC")
  bic <- read_grid_run(bic_path, "BIC")
  validate_matching_grid_provenance(gic, bic)
  gic_selection <- select_grid_parameters(gic)
  bic_selection <- select_grid_parameters(bic)
  source_provenance <- gic$source_provenance

  cache <- list(
    schema_version = 3L,
    provenance = list(
      simulation_generator = "empirical",
      scenario_interpretation = paste(
        "Empirically calibrated spectral integration change with reciprocal",
        "marginal-variance scaling"
      ),
      original_generator_available = TRUE,
      n_replicates_per_setting = 100L,
      tree_tip_count = 250L,
      n_true_shifts = 5L,
      min_shift_tips = 10L,
      max_shift_tips = 20L,
      integration_power_range = c(0.5, 1.25),
      integration_exclude_range = c(0.8, 1.1),
      fuzzy_distance = 2L,
      selection = list(
        primary_metric = "fuzzy_balanced_accuracy",
        scenario_weights = c(proportional = 0.5, correlation = 0.5),
        max_false_positive_rate = 0.10,
        max_any_false_positive = 0.50,
        min_evaluable_fraction = 0.50,
        tie_break = "conservative"
      ),
      generator_commit = source_provenance$commit,
      generator_source_md5 = source_provenance$source_md5,
      generator_source_dirty = source_provenance$source_dirty,
      generator_source_status = source_provenance$source_status,
      source_md5 = c(
        GIC = unname(tools::md5sum(gic$path)),
        BIC = unname(tools::md5sum(bic$path))
      ),
      R = gic$wrapper$provenance$R,
      platform = gic$wrapper$provenance$platform
    ),
    fixed_settings = rbind(
      build_fixed_rows(gic_selection),
      build_fixed_rows(bic_selection)
    ),
    tuning = list(
      gic = build_tuning_table(gic),
      bic = build_tuning_table(bic),
      selected = rbind(
        build_selected_row(gic_selection),
        build_selected_row(bic_selection)
      )
    ),
    grid_summary = list(
      gic = gic$grid$summary_table,
      bic = bic$grid$summary_table
    )
  )
  rownames(cache$fixed_settings) <- NULL
  rownames(cache$tuning$selected) <- NULL
  validate_cache_object(cache)
  cache
}

write_simulation_cache <- function(cache, out_path) {
  validate_cache_object(cache)
  dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)
  temp_path <- tempfile(
    pattern = paste0(".", basename(out_path), "-"),
    tmpdir = dirname(out_path),
    fileext = ".tmp"
  )
  on.exit(unlink(temp_path), add = TRUE)
  saveRDS(cache, temp_path, compress = "xz", version = 3)
  validate_cache_object(readRDS(temp_path))
  if (!file.rename(temp_path, out_path)) {
    stop("Could not atomically replace ", out_path)
  }
  message("Wrote ", out_path)
  invisible(out_path)
}

build_simulation_cache_main <- function(args = commandArgs(trailingOnly = TRUE)) {
  if (!file.exists("DESCRIPTION")) {
    stop("Run this script from the package root.")
  }
  if (length(args) != 2L) {
    stop(
      paste0(
        "Usage: Rscript data-raw/build_simulation_study_vignette_cache.R ",
        "<empirical-GIC-grid.rds> <empirical-BIC-grid.rds>"
      )
    )
  }
  if (!requireNamespace("pkgload", quietly = TRUE)) {
    stop("Building the vignette cache requires the suggested pkgload package.")
  }
  pkgload::load_all(".", quiet = TRUE, export_all = FALSE, helpers = FALSE)
  ensure_grid_runner_helpers()
  cache <- build_simulation_cache(args[[1L]], args[[2L]])
  out_path <- file.path(
    "inst",
    "extdata",
    "simulation-study-cache",
    "passerine_preview_tables.rds"
  )
  write_simulation_cache(cache, out_path)
}

if (sys.nframe() == 0L) {
  build_simulation_cache_main()
}
