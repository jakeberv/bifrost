if (!file.exists("DESCRIPTION")) {
  stop("Run this script from the package root.")
}

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) {
  stop(
    paste0(
      "Usage: Rscript data-raw/build_simulation_study_vignette_cache.R ",
      "<empirical-GIC-grid.rds> <empirical-BIC-grid.rds>"
    )
  )
}

read_grid_run <- function(path, expected_ic) {
  if (!file.exists(path)) {
    stop("Grid result not found: ", path)
  }
  wrapper <- readRDS(path)
  grid <- if (!is.null(wrapper$result)) wrapper$result else wrapper
  if (!identical(grid$IC, expected_ic)) {
    stop("Expected an ", expected_ic, " grid in ", path)
  }
  simulation_generators <- grid$simulation_generators
  if (!identical(unname(simulation_generators), rep("empirical", 3L))) {
    stop("All cached scenarios must explicitly use simulation_generator = 'empirical'.")
  }
  grid$simulation_generators <- simulation_generators
  grid$simulation_designs <- NULL
  if (!identical(grid$null_replicates, 100L) ||
      !identical(grid$recovery_replicates, 100L)) {
    stop("The vignette cache requires 100 replicates per grid setting.")
  }
  if (is.null(grid$summary_table) || nrow(grid$summary_table) != 6L) {
    stop("Each vignette grid must contain the six documented settings.")
  }
  list(wrapper = wrapper, grid = grid)
}

gic <- read_grid_run(args[[1L]], "GIC")
bic <- read_grid_run(args[[2L]], "BIC")

fixed_setting_row <- function(grid_run) {
  summary <- grid_run$grid$summary_table
  row <- summary[
    summary$shift_acceptance_threshold == 10 &
      summary$min_descendant_tips == 10,
    ,
    drop = FALSE
  ]
  if (nrow(row) != 1L) {
    stop("Each grid must contain exactly one threshold-10/min-clade-10 row.")
  }
  row
}

build_fixed_rows <- function(grid_run) {
  row <- fixed_setting_row(grid_run)
  data.frame(
    IC = rep(row$IC, 3L),
    Scenario = c("Null", "Proportional", "Integration-rate"),
    `Mean FP` = c(row$null_mean_false_positive_rate, NA_real_, NA_real_),
    `Any FP` = c(row$null_fraction_any_false_positive, NA_real_, NA_real_),
    `Strict F1` = c(
      NA_real_,
      row$proportional_strict_f1,
      row$correlation_strict_f1
    ),
    `Fuzzy F1` = c(
      NA_real_,
      row$proportional_fuzzy_f1,
      row$correlation_fuzzy_f1
    ),
    `Fuzzy recall` = c(
      NA_real_,
      row$proportional_fuzzy_recall,
      row$correlation_fuzzy_recall
    ),
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
    check.names = FALSE
  )
}

build_tuning_table <- function(grid_run) {
  summary <- grid_run$grid$summary_table
  data.frame(
    Threshold = summary$shift_acceptance_threshold,
    `Min clade` = summary$min_descendant_tips,
    `Null FP` = summary$null_mean_false_positive_rate,
    `Prop. Fuzzy F1` = summary$proportional_fuzzy_f1,
    `Integration Fuzzy F1` = summary$correlation_fuzzy_f1,
    check.names = FALSE
  )
}

build_selected_row <- function(grid_run) {
  row <- fixed_setting_row(grid_run)
  data.frame(
    IC = row$IC,
    Threshold = row$shift_acceptance_threshold,
    `Min clade` = row$min_descendant_tips,
    `Null FP` = row$null_mean_false_positive_rate,
    `Prop. Fuzzy F1` = row$proportional_fuzzy_f1,
    `Integration Fuzzy F1` = row$correlation_fuzzy_f1,
    check.names = FALSE
  )
}

commits <- unique(c(gic$wrapper$provenance$commit, bic$wrapper$provenance$commit))
if (length(commits) != 1L || is.na(commits) || !nzchar(commits)) {
  stop("GIC and BIC source runs must record the same generator commit.")
}

cache_obj <- list(
  schema_version = 2L,
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
    generator_commit = commits,
    source_md5 = c(
      GIC = unname(tools::md5sum(args[[1L]])),
      BIC = unname(tools::md5sum(args[[2L]]))
    ),
    R = gic$wrapper$provenance$R,
    platform = gic$wrapper$provenance$platform
  ),
  fixed_settings = rbind(
    build_fixed_rows(gic),
    build_fixed_rows(bic)
  ),
  tuning = list(
    gic = build_tuning_table(gic),
    bic = build_tuning_table(bic),
    selected = rbind(
      build_selected_row(gic),
      build_selected_row(bic)
    )
  ),
  grid_summary = list(
    gic = gic$grid$summary_table,
    bic = bic$grid$summary_table
  )
)
rownames(cache_obj$fixed_settings) <- NULL
rownames(cache_obj$tuning$selected) <- NULL

out_dir <- file.path("inst", "extdata", "simulation-study-cache")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
out_path <- file.path(out_dir, "passerine_preview_tables.rds")
saveRDS(cache_obj, out_path, compress = "xz", version = 3)
message("Wrote ", out_path)
