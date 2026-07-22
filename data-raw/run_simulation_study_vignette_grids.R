usage <- function() {
  paste(
    "Usage: Rscript data-raw/run_simulation_study_vignette_grids.R",
    "--mode=smoke|full [--ic=both|GIC|BIC] [--cores=N]",
    "[--output-dir=PATH] [--overwrite]"
  )
}

parse_cli_args <- function(args, mode = c("smoke", "full"),
                           detect_cores = function() {
                             parallel::detectCores(logical = FALSE)
                           }) {
  values <- list(
    mode = NULL,
    ic = "both",
    cores = NULL,
    output_dir = file.path(
      "local-cache",
      "simulation-study-vignette-grids"
    ),
    overwrite = FALSE
  )

  for (arg in args) {
    if (identical(arg, "--overwrite")) {
      values$overwrite <- TRUE
    } else if (startsWith(arg, "--mode=")) {
      values$mode <- sub("^--mode=", "", arg)
    } else if (startsWith(arg, "--ic=")) {
      values$ic <- sub("^--ic=", "", arg)
    } else if (startsWith(arg, "--cores=")) {
      values$cores <- suppressWarnings(as.integer(sub("^--cores=", "", arg)))
    } else if (startsWith(arg, "--output-dir=")) {
      values$output_dir <- sub("^--output-dir=", "", arg)
    } else if (arg %in% c("-h", "--help")) {
      message(usage())
      quit(save = "no", status = 0L)
    } else {
      stop("Unknown argument: ", arg, "\n", usage())
    }
  }

  if (is.null(values$mode)) {
    stop("--mode must be supplied explicitly.\n", usage())
  }
  values$mode <- match.arg(values$mode, mode)
  values$ic <- match.arg(values$ic, c("both", "GIC", "BIC"))
  if (!is.null(values$cores) &&
      (is.na(values$cores) || values$cores < 1L)) {
    stop("--cores must be a positive integer.")
  }
  if (!nzchar(values$output_dir)) {
    stop("--output-dir must not be empty.")
  }

  if (is.null(values$cores)) {
    available_cores <- detect_cores()
    if (length(available_cores) != 1L || is.na(available_cores)) {
      available_cores <- 1L
    }
    setting_count <- if (identical(values$mode, "full")) 6L else 1L
    values$cores <- max(1L, min(setting_count, available_cores - 1L))
  }
  values
}

simulation_grid_design <- function(mode = c("smoke", "full")) {
  mode <- match.arg(mode)
  full <- identical(mode, "full")
  list(
    mode = mode,
    shift_acceptance_thresholds = if (full) c(10, 20, 30) else 10,
    min_descendant_tips_values = if (full) c(10L, 20L) else 10L,
    tree_tip_count = 250L,
    null_replicates = if (full) 100L else 1L,
    recovery_replicates = if (full) 100L else 1L,
    seed = 5L,
    null_simulation_options = list(
      simulation_generator = "empirical"
    ),
    proportional_simulation_options = list(
      simulation_generator = "empirical",
      num_shifts = 5L,
      min_shift_tips = 10L,
      max_shift_tips = 20L,
      scale_factor_range = c(0.1, 2.0),
      exclude_range = c(0.5, 1.5),
      buffer = 3L
    ),
    correlation_simulation_options = list(
      simulation_generator = "empirical",
      integration_power_range = c(0.5, 1.25),
      integration_exclude_range = c(0.8, 1.1)
    ),
    base_search_options = list(
      formula = "trait_data ~ 1",
      method = "LL",
      error = TRUE,
      num_cores = 1L,
      uncertaintyweights_par = TRUE,
      plot = FALSE,
      store_model_fit_history = FALSE,
      progress = FALSE
    ),
    fuzzy_distance = 2L,
    weighted = TRUE,
    store_studies = FALSE
  )
}

grid_recovery_metric_suffixes <- function() {
  c(
    "strict_precision",
    "strict_recall",
    "strict_f1",
    "strict_specificity",
    "strict_fpr",
    "strict_balanced_accuracy",
    "fuzzy_precision",
    "fuzzy_recall",
    "fuzzy_f1",
    "fuzzy_specificity",
    "fuzzy_fpr",
    "fuzzy_balanced_accuracy"
  )
}

required_grid_summary_columns <- function() {
  metric_suffixes <- grid_recovery_metric_suffixes()
  c(
    "setting_id",
    "IC",
    "shift_acceptance_threshold",
    "min_descendant_tips",
    "null_seed",
    "proportional_seed",
    "correlation_seed",
    "null_mean_false_positive_rate",
    "null_fraction_any_false_positive",
    "null_evaluable_fraction",
    "null_completion_rate",
    "null_failure_rate",
    "null_mean_inferred_shifts",
    paste0("proportional_", metric_suffixes),
    "proportional_evaluable_fraction",
    "proportional_completion_rate",
    "proportional_failure_rate",
    "proportional_mean_inferred_shifts",
    "proportional_weighted_fuzzy_f1",
    paste0("correlation_", metric_suffixes),
    "correlation_evaluable_fraction",
    "correlation_completion_rate",
    "correlation_failure_rate",
    "correlation_mean_inferred_shifts",
    "correlation_weighted_fuzzy_f1"
  )
}

validate_grid_summary_table <- function(summary_table, expected_ic,
                                        expected_design,
                                        tolerance = 1e-10) {
  invalid <- function(detail) {
    stop("Invalid grid summary: ", detail, ".", call. = FALSE)
  }
  if (!is.data.frame(summary_table)) {
    invalid("summary_table must be a data frame")
  }
  missing_columns <- setdiff(
    required_grid_summary_columns(),
    names(summary_table)
  )
  if (length(missing_columns) > 0L) {
    invalid(paste0(
      "missing complete strict/fuzzy metrics: ",
      paste(missing_columns, collapse = ", ")
    ))
  }

  expected_settings <- expand.grid(
    shift_acceptance_threshold =
      expected_design$shift_acceptance_thresholds,
    min_descendant_tips = expected_design$min_descendant_tips_values,
    KEEP.OUT.ATTRS = FALSE
  )
  if (nrow(summary_table) != nrow(expected_settings) ||
      !identical(summary_table$setting_id, seq_len(nrow(expected_settings))) ||
      !identical(
        summary_table[
          , names(expected_settings),
          drop = FALSE
        ],
        expected_settings
      ) ||
      !is.character(summary_table$IC) ||
      !identical(summary_table$IC, rep(expected_ic, nrow(expected_settings)))) {
    invalid("rows do not match the exact requested IC and setting grid")
  }

  required_numeric <- setdiff(required_grid_summary_columns(), "IC")
  non_numeric <- required_numeric[!vapply(
    summary_table[required_numeric],
    is.numeric,
    logical(1L)
  )]
  if (length(non_numeric) > 0L) {
    invalid(paste0(
      "required columns must be numeric: ",
      paste(non_numeric, collapse = ", ")
    ))
  }
  non_finite <- required_numeric[!vapply(
    summary_table[required_numeric],
    function(x) all(is.finite(x)),
    logical(1L)
  )]
  if (length(non_finite) > 0L) {
    invalid(paste0(
      "required columns must be finite: ",
      paste(non_finite, collapse = ", ")
    ))
  }

  probability_columns <- c(
    "null_mean_false_positive_rate",
    "null_fraction_any_false_positive",
    "null_evaluable_fraction",
    "null_completion_rate",
    "null_failure_rate"
  )
  for (scenario in c("proportional", "correlation")) {
    probability_columns <- c(
      probability_columns,
      paste0(scenario, "_", grid_recovery_metric_suffixes()),
      paste0(scenario, c(
        "_evaluable_fraction",
        "_completion_rate",
        "_failure_rate",
        "_weighted_fuzzy_f1"
      ))
    )
  }
  outside_unit_interval <- probability_columns[!vapply(
    summary_table[probability_columns],
    function(x) all(x >= 0 & x <= 1),
    logical(1L)
  )]
  if (length(outside_unit_interval) > 0L) {
    invalid(paste0(
      "probability metrics must be between 0 and 1: ",
      paste(outside_unit_interval, collapse = ", ")
    ))
  }

  count_columns <- c(
    "null_mean_inferred_shifts",
    "proportional_mean_inferred_shifts",
    "correlation_mean_inferred_shifts"
  )
  if (any(unlist(summary_table[count_columns], use.names = FALSE) < 0)) {
    invalid("mean inferred-shift counts must be non-negative")
  }

  for (scenario in c("null", "proportional", "correlation")) {
    completion <- summary_table[[paste0(scenario, "_completion_rate")]]
    failure <- summary_table[[paste0(scenario, "_failure_rate")]]
    if (any(abs(completion + failure - 1) > tolerance)) {
      invalid(paste0(
        scenario,
        " completion and failure rates must sum to 1"
      ))
    }
  }

  for (scenario in c("proportional", "correlation")) {
    for (match in c("strict", "fuzzy")) {
      prefix <- paste0(scenario, "_", match, "_")
      precision <- summary_table[[paste0(prefix, "precision")]]
      recall <- summary_table[[paste0(prefix, "recall")]]
      f1 <- summary_table[[paste0(prefix, "f1")]]
      specificity <- summary_table[[paste0(prefix, "specificity")]]
      fpr <- summary_table[[paste0(prefix, "fpr")]]
      balanced_accuracy <-
        summary_table[[paste0(prefix, "balanced_accuracy")]]
      expected_f1 <- ifelse(
        precision + recall == 0,
        0,
        2 * precision * recall / (precision + recall)
      )
      if (any(abs(f1 - expected_f1) > tolerance)) {
        invalid(paste0(prefix, "F1 is inconsistent with precision and recall"))
      }
      if (any(abs(fpr - (1 - specificity)) > tolerance)) {
        invalid(paste0(prefix, "FPR is inconsistent with specificity"))
      }
      if (any(abs(
        balanced_accuracy - (recall + specificity) / 2
      ) > tolerance)) {
        invalid(paste0(
          prefix,
          "balanced accuracy is inconsistent with recall and specificity"
        ))
      }
    }
  }
  invisible(summary_table)
}

relevant_source_files <- function(root = ".") {
  r_files <- sort(list.files(
    file.path(root, "R"),
    pattern = "\\.R$",
    full.names = FALSE
  ))
  relative_paths <- c(
    "data-raw/run_simulation_study_vignette_grids.R",
    "DESCRIPTION",
    "NAMESPACE",
    file.path("R", r_files),
    file.path(
      "inst",
      "extdata",
      "avian-skeleton",
      c("passerine_bodyplan_tree.tre", "passerine_bodyplan_data.RDS")
    )
  )
  relative_paths <- gsub("^\\./", "", relative_paths)
  paths <- file.path(root, relative_paths)
  names(paths) <- relative_paths
  missing_paths <- names(paths)[!file.exists(paths)]
  if (length(missing_paths) > 0L) {
    stop(
      "Relevant simulation source files are missing: ",
      paste(missing_paths, collapse = ", ")
    )
  }
  paths
}

read_generator_commit <- function(root = ".") {
  commit <- tryCatch(
    system2(
      "git",
      c("-C", root, "rev-parse", "HEAD"),
      stdout = TRUE,
      stderr = FALSE
    ),
    error = function(e) character()
  )
  if (length(commit) != 1L || !nzchar(commit)) {
    stop("Could not determine the generator commit with git rev-parse HEAD.")
  }
  commit
}

read_worktree_status <- function(root = ".") {
  status <- tryCatch(
    system2(
      "git",
      c(
        "-C",
        root,
        "status",
        "--porcelain=v1",
        "--untracked-files=all"
      ),
      stdout = TRUE,
      stderr = FALSE
    ),
    error = function(e) character()
  )
  as.character(status)
}

status_entry_path <- function(status_line) {
  path <- substring(status_line, 4L)
  if (grepl(" -> ", path, fixed = TRUE)) {
    path <- sub("^.* -> ", "", path)
  }
  sub('^"(.*)"$', "\\1", path)
}

filter_relevant_status <- function(status_lines, relevant_paths) {
  if (length(status_lines) == 0L) {
    return(character())
  }
  status_paths <- vapply(status_lines, status_entry_path, character(1L))
  status_lines[status_paths %in% relevant_paths]
}

collect_source_provenance <- function(
    root = ".",
    relevant_files = relevant_source_files(root),
    commit = read_generator_commit(root),
    status_lines = read_worktree_status(root)) {
  if (is.null(names(relevant_files)) || any(!nzchar(names(relevant_files)))) {
    stop("relevant_files must be named with repository-relative paths.")
  }
  if (any(!file.exists(relevant_files))) {
    stop("Every relevant source file must exist before hashing.")
  }
  hashes <- unname(tools::md5sum(relevant_files))
  names(hashes) <- names(relevant_files)
  relevant_status <- filter_relevant_status(
    status_lines,
    names(relevant_files)
  )
  list(
    commit = commit,
    source_files = names(relevant_files),
    source_md5 = hashes,
    source_dirty = length(relevant_status) > 0L,
    source_status = unname(relevant_status)
  )
}

validate_source_provenance <- function(provenance) {
  invalid <- function(detail) {
    stop("Invalid source provenance: ", detail, ".", call. = FALSE)
  }
  if (!is.list(provenance)) {
    invalid("record must be a list")
  }
  if (!is.character(provenance$commit) ||
      length(provenance$commit) != 1L ||
      is.na(provenance$commit) ||
      !nzchar(trimws(provenance$commit))) {
    invalid("commit must be a nonempty character scalar")
  }

  source_files <- provenance$source_files
  if (!is.character(source_files) || length(source_files) < 1L ||
      anyNA(source_files) || any(!nzchar(source_files)) ||
      anyDuplicated(source_files) > 0L) {
    invalid("source_files must be a nonempty unique character inventory")
  }
  unsafe_paths <- startsWith(source_files, "/") |
    grepl("^[A-Za-z]:[/\\\\]", source_files) |
    vapply(strsplit(source_files, "/", fixed = TRUE), function(parts) {
      ".." %in% parts
    }, logical(1L))
  if (any(unsafe_paths)) {
    invalid("source_files must contain safe repository-relative paths")
  }

  hashes <- provenance$source_md5
  if (!is.character(hashes) || length(hashes) < 1L || anyNA(hashes) ||
      is.null(names(hashes)) || any(!nzchar(names(hashes))) ||
      !identical(names(hashes), source_files) ||
      any(!grepl("^[0-9a-f]{32}$", hashes))) {
    invalid("source_md5 must contain named lowercase MD5 hashes for every file")
  }
  if (!is.logical(provenance$source_dirty) ||
      length(provenance$source_dirty) != 1L ||
      is.na(provenance$source_dirty)) {
    invalid("source_dirty must be a non-missing logical scalar")
  }
  status <- provenance$source_status
  if (!is.character(status) || anyNA(status) || any(!nzchar(status))) {
    invalid("source_status must be a valid character vector")
  }
  if (!identical(provenance$source_dirty, length(status) > 0L)) {
    invalid("source_status must be nonempty exactly when source_dirty is TRUE")
  }
  if (length(status) > 0L) {
    if (any(nchar(status) < 4L)) {
      invalid("source_status entries must use git porcelain format")
    }
    status_paths <- vapply(status, status_entry_path, character(1L))
    if (any(!status_paths %in% source_files)) {
      invalid("source_status may reference only inventoried source files")
    }
  }
  invisible(provenance)
}

normalize_source_provenance <- function(provenance) {
  if (!is.list(provenance)) {
    stop("Invalid source provenance: record must be a list.", call. = FALSE)
  }
  source_files <- provenance$source_files
  if (is.null(source_files) && !is.null(names(provenance$source_md5))) {
    source_files <- names(provenance$source_md5)
  }
  normalized <- list(
    commit = provenance$commit,
    source_files = source_files,
    source_md5 = provenance$source_md5,
    source_dirty = provenance$source_dirty,
    source_status = provenance$source_status
  )
  validate_source_provenance(normalized)
  normalized
}

matching_source_provenance <- function(actual, expected) {
  actual <- normalize_source_provenance(actual)
  expected <- normalize_source_provenance(expected)
  identical(actual$commit, expected$commit) &&
    identical(actual$source_files, expected$source_files) &&
    identical(actual$source_md5, expected$source_md5) &&
    identical(actual$source_dirty, expected$source_dirty)
}

validate_grid_wrapper <- function(wrapper, expected_ic, expected_design,
                                  expected_source_provenance) {
  if (!is.list(wrapper) || !identical(wrapper$schema_version, 1L)) {
    stop("Grid output must be a schema-1 provenance wrapper.")
  }
  if (!is.list(wrapper$provenance) || is.null(wrapper$result)) {
    stop("Grid output must contain provenance and result components.")
  }
  if (!identical(wrapper$provenance$design, expected_design)) {
    stop("Grid output provenance does not match the complete requested design.")
  }
  if (!matching_source_provenance(
    wrapper$provenance,
    expected_source_provenance
  )) {
    stop("Grid output provenance does not match current source provenance.")
  }
  if (!is.character(wrapper$provenance$R) ||
      length(wrapper$provenance$R) != 1L ||
      !nzchar(wrapper$provenance$R) ||
      !is.character(wrapper$provenance$platform) ||
      length(wrapper$provenance$platform) != 1L ||
      !nzchar(wrapper$provenance$platform)) {
    stop("Grid output provenance must record R and platform.")
  }

  grid <- wrapper$result
  if (!inherits(grid, "bifrost_search_tuning_grid") ||
      !identical(grid$IC, expected_ic)) {
    stop("Grid output does not contain the expected ", expected_ic, " result.")
  }
  validate_grid_summary_table(
    grid$summary_table,
    expected_ic,
    expected_design
  )
  expected_settings <- expand.grid(
    shift_acceptance_threshold =
      expected_design$shift_acceptance_thresholds,
    min_descendant_tips = expected_design$min_descendant_tips_values,
    KEEP.OUT.ATTRS = FALSE
  )
  expected_grid <- data.frame(
    setting_id = seq_len(nrow(expected_settings)),
    expected_settings
  )
  actual_settings <- grid$summary_table[
    , c("shift_acceptance_threshold", "min_descendant_tips"),
    drop = FALSE
  ]
  if (!identical(actual_settings, expected_settings) ||
      !identical(grid$grid, expected_grid) ||
      !all(grid$summary_table$IC == expected_ic)) {
    stop("Grid output does not match the exact requested setting grid.")
  }
  if (!identical(grid$null_replicates, expected_design$null_replicates) ||
      !identical(
        grid$recovery_replicates,
        expected_design$recovery_replicates
      ) ||
      !identical(grid$fuzzy_distance, expected_design$fuzzy_distance) ||
      !identical(grid$weighted, expected_design$weighted) ||
      !identical(grid$store_studies, expected_design$store_studies) ||
      !identical(
        grid$base_search_options,
        expected_design$base_search_options
      )) {
    stop("Grid output does not match the requested study controls.")
  }
  if (!identical(unname(grid$simulation_generators), rep("empirical", 3L))) {
    stop("Every simulation scenario must use the empirical generator.")
  }
  invisible(wrapper)
}

new_grid_wrapper <- function(grid, design, source_provenance,
                             setting_workers,
                             started_at = Sys.time(),
                             finished_at = Sys.time()) {
  source_provenance <- normalize_source_provenance(source_provenance)
  list(
    schema_version = 1L,
    provenance = c(
      list(
        design = design,
        script = "data-raw/run_simulation_study_vignette_grids.R",
        started_utc = format(started_at, tz = "UTC", usetz = TRUE),
        finished_utc = format(finished_at, tz = "UTC", usetz = TRUE),
        elapsed_seconds = unname(as.numeric(difftime(
          finished_at,
          started_at,
          units = "secs"
        ))),
        R = R.version.string,
        platform = R.version$platform,
        setting_workers = as.integer(setting_workers)
      ),
      source_provenance
    ),
    result = grid
  )
}

write_grid_wrapper <- function(wrapper, path, overwrite = FALSE) {
  if (file.exists(path) && !overwrite) {
    stop("Grid output already exists; use --overwrite to replace it: ", path)
  }
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  temp_path <- tempfile(
    pattern = paste0(".", basename(path), "-"),
    tmpdir = dirname(path),
    fileext = ".tmp"
  )
  on.exit(unlink(temp_path), add = TRUE)
  saveRDS(wrapper, temp_path, compress = "xz", version = 3)
  if (!identical(readRDS(temp_path), wrapper)) {
    stop("Serialized grid wrapper did not round-trip before replacement.")
  }
  if (!file.rename(temp_path, path)) {
    stop("Could not atomically move the completed grid to ", path)
  }
  invisible(path)
}

run_or_resume_grid <- function(output_path, expected_ic, design,
                               source_provenance, setting_workers,
                               overwrite, run_grid) {
  if (file.exists(output_path) && !overwrite) {
    wrapper <- readRDS(output_path)
    validate_grid_wrapper(
      wrapper,
      expected_ic,
      design,
      source_provenance
    )
    message("Skipping existing validated output: ", output_path)
    return(list(status = "resumed", path = output_path, wrapper = wrapper))
  }

  started_at <- Sys.time()
  grid <- run_grid()
  finished_at <- Sys.time()
  wrapper <- new_grid_wrapper(
    grid,
    design,
    source_provenance,
    setting_workers,
    started_at,
    finished_at
  )
  validate_grid_wrapper(
    wrapper,
    expected_ic,
    design,
    source_provenance
  )
  write_grid_wrapper(wrapper, output_path, overwrite = overwrite)
  validate_grid_wrapper(
    readRDS(output_path),
    expected_ic,
    design,
    source_provenance
  )
  message("Wrote and validated ", output_path)
  list(status = "written", path = output_path, wrapper = wrapper)
}

build_bodyplan_template <- function() {
  tree_path <- file.path(
    "inst",
    "extdata",
    "avian-skeleton",
    "passerine_bodyplan_tree.tre"
  )
  trait_path <- file.path(
    "inst",
    "extdata",
    "avian-skeleton",
    "passerine_bodyplan_data.RDS"
  )
  if (!file.exists(tree_path) || !file.exists(trait_path)) {
    stop("The packaged passerine tree and body-plan data are required.")
  }

  bird_tree <- ape::read.tree(tree_path)
  bodyplan_data <- readRDS(trait_path)
  bodyplan_data <- as.matrix(
    bodyplan_data[bird_tree$tip.label, , drop = FALSE]
  )
  skeletal_cols <- setdiff(colnames(bodyplan_data), "vertnet_mass")
  bodyplan_data <- bodyplan_data[
    , c(skeletal_cols, "vertnet_mass"),
    drop = FALSE
  ]
  formula_str <- sprintf(
    "trait_data[, 1:%d] ~ trait_data[, %d]",
    length(skeletal_cols),
    ncol(bodyplan_data)
  )

  createSimulationTemplate(
    baseline_tree = bird_tree,
    trait_data = bodyplan_data,
    formula = formula_str,
    response_columns = seq_along(skeletal_cols),
    predictor_columns = ncol(bodyplan_data),
    method = "LL",
    error = TRUE
  )
}

run_simulation_grid_main <- function(args = commandArgs(trailingOnly = TRUE)) {
  if (!file.exists("DESCRIPTION")) {
    stop("Run this script from the package root.")
  }
  cli <- parse_cli_args(args)

  Sys.setenv(
    OMP_NUM_THREADS = "1",
    OPENBLAS_NUM_THREADS = "1",
    MKL_NUM_THREADS = "1",
    VECLIB_MAXIMUM_THREADS = "1",
    BLIS_NUM_THREADS = "1",
    NUMEXPR_NUM_THREADS = "1"
  )
  if (!requireNamespace("pkgload", quietly = TRUE)) {
    stop("Regenerating the grids requires the suggested pkgload package.")
  }
  pkgload::load_all(".", quiet = TRUE, export_all = FALSE, helpers = FALSE)

  design <- simulation_grid_design(cli$mode)
  source_provenance <- collect_source_provenance()
  message("Fitting the empirical passerine simulation template...")
  bodyplan_template <- build_bodyplan_template()
  ics <- if (identical(cli$ic, "both")) c("GIC", "BIC") else cli$ic

  for (ic in ics) {
    suffix <- if (identical(cli$mode, "smoke")) "-smoke" else ""
    output_path <- file.path(
      cli$output_dir,
      paste0("empirical-", ic, "-grid", suffix, ".rds")
    )
    run_grid <- function() {
      message(
        "Running ", cli$mode, " ", ic, " grid with ", cli$cores,
        " parallel setting worker(s)..."
      )
      runSearchTuningGrid(
        template = bodyplan_template,
        IC = ic,
        shift_acceptance_thresholds =
          design$shift_acceptance_thresholds,
        min_descendant_tips_values = design$min_descendant_tips_values,
        tree_tip_count = design$tree_tip_count,
        null_replicates = design$null_replicates,
        recovery_replicates = design$recovery_replicates,
        null_simulation_options = design$null_simulation_options,
        proportional_simulation_options =
          design$proportional_simulation_options,
        correlation_simulation_options =
          design$correlation_simulation_options,
        base_search_options = design$base_search_options,
        fuzzy_distance = design$fuzzy_distance,
        weighted = design$weighted,
        num_cores = cli$cores,
        seed = design$seed,
        store_studies = design$store_studies
      )
    }
    run_or_resume_grid(
      output_path,
      expected_ic = ic,
      design = design,
      source_provenance = source_provenance,
      setting_workers = cli$cores,
      overwrite = cli$overwrite,
      run_grid = run_grid
    )
  }
  invisible(NULL)
}

if (sys.nframe() == 0L) {
  run_simulation_grid_main()
}
