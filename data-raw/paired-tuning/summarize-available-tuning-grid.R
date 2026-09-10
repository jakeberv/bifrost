#!/usr/bin/env Rscript
# Analysis-only companion. It never generates datasets, fits models, or changes
# the source-fingerprinted campaign runner and its saved raw outputs.

collect_available_tuning <- function(e, groups, jobs, root, failure_log, namespace,
                                     provenance, design_fingerprint) {
  lines <- grep("^Failed ", readLines(failure_log, warn = FALSE), value = TRUE)
  pattern <- "^Failed ([a-z0-9-]+): Failed to place all requested shifts after 100 attempts[.]$"
  if (any(!grepl(pattern, lines))) stop("Unsupported failure in campaign log.")
  excluded <- sub(pattern, "\\1", lines)
  if (anyDuplicated(excluded)) stop("Duplicate placement failure in campaign log.")
  if (any(!excluded %in% jobs$job_id)) stop("Unknown placement-failure job ID.")
  expected <- vapply(seq_len(nrow(jobs)), function(i) e$job_output_path(root, "full", jobs[i, ]), character(1))
  observed <- list.files(file.path(root, "full", "replicates"), pattern = "[.]rds$",
                         recursive = TRUE, full.names = TRUE)
  if (length(setdiff(observed, expected))) stop("Unexpected replicate output.")
  present <- file.exists(expected)
  if (any(jobs$job_id[present] %in% excluded)) stop("Contradictory generated output and placement failure.")
  if (!setequal(jobs$job_id[!present], excluded)) stop("Unaccounted missing replicate outputs.")
  for (id in excluded) {
    group_id <- jobs$group_id[match(id, jobs$job_id)]
    checkpoint <- file.path(root, "full", "checkpoints", group_id, id)
    if (length(list.files(checkpoint, recursive = TRUE, all.files = TRUE, no.. = TRUE))) {
      stop("Placement failure has an existing checkpoint: ", id)
    }
  }
  attempts <- jobs[, c("job_id", "group_id", "replicate")]
  attempts$status <- ifelse(present, "generated", "placement_failed")
  attempts$reason <- ifelse(present, NA_character_, "Failed to place all requested shifts after 100 attempts.")
  summaries <- list(); diagnostics <- list(); group_counts <- list()
  for (i in seq_len(nrow(groups))) {
    group <- groups[i, ]; configs <- e$search_configs(group)
    indices <- which(jobs$group_id == group$id & present)
    if (!length(indices)) stop("No generated datasets for group: ", group$id)
    wrappers <- lapply(indices, function(j) {
      wrapper <- readRDS(expected[j])
      e$validate_job_wrapper(wrapper, jobs[j, ], group, configs, provenance, design_fingerprint)
      wrapper
    })
    attempted <- sum(jobs$group_id == group$id)
    generated <- length(wrappers)
    group_counts[[i]] <- data.frame(group_id = group$id, scenario = group$scenario,
      n_attempted = attempted, n_generated = generated, n_placement_failed = attempted - generated)
    for (j in seq_len(nrow(configs))) {
      config <- configs[j, ]
      summary <- e$summarize_config(group, config, wrappers, namespace, provenance,
                                     design_fingerprint, "full")
      summary$per_replicate$replicate <- jobs$replicate[indices]
      summary$per_replicate$job_id <- jobs$job_id[indices]
      summary$summary_row$n_attempted <- attempted
      summary$summary_row$n_generated <- generated
      summary$summary_row$n_placement_failed <- attempted - generated
      summary$summary_row$n_search_attempted <- generated
      summary$summary_row$n_evaluated <- summary$summary_row$n_evaluable_replicates
      summary$summary_row$generation_success_rate <- generated / attempted
      summary$accounting_scope <- "Search metrics are conditional on successful dataset generation."
      summaries[[config$config_id]] <- summary
      results <- lapply(wrappers, function(w) w$results[[config$config_id]])
      diag <- summary$per_replicate
      diag$scenario <- group$scenario
      diag$IC <- config$IC
      diag$shift_acceptance_threshold <- config$shift_acceptance_threshold
      diag$min_descendant_tips <- config$min_descendant_tips
      warning_text <- lapply(results, function(r) {
        x <- as.character(unlist(r$warnings, use.names = FALSE))
        x[!is.na(x) & nzchar(x)]
      })
      diag$warning_count <- lengths(warning_text)
      diag$warning_messages <- vapply(warning_text, paste, character(1), collapse = " | ")
      diag$nonfinite_ic <- vapply(results, function(r) {
        length(r$optimal_ic) != 1L || !is.finite(r$optimal_ic)
      }, logical(1))
      diagnostics[[config$config_id]] <- diag
    }
  }
  bind_rows <- function(rows) {
    columns <- unique(unlist(lapply(rows, names)))
    rows <- lapply(rows, function(row) {
      for (name in setdiff(columns, names(row))) row[[name]] <- NA
      row[, columns, drop = FALSE]
    })
    out <- do.call(rbind, rows); rownames(out) <- NULL; out
  }
  list(schema_version = 1L, suite = e$suite_name, mode = "full",
    all_groups_complete = all(present), all_attempts_accounted_for = TRUE,
    accounting_scope = "Search metrics are conditional on successful dataset generation; placement failures are excluded from search-performance denominators.",
    source_fingerprint = provenance$source_fingerprint, design_fingerprint = design_fingerprint,
    provenance = provenance, selection_policy = e$selection_policy,
    attempts = attempts, group_counts = do.call(rbind, group_counts),
    summary_table = bind_rows(lapply(summaries, `[[`, "summary_row")),
    summaries = summaries, diagnostics = bind_rows(diagnostics))
}

write_available_tuning <- function(result, e, namespace, destination, audit) {
  if (dir.exists(destination)) stop("Summary destination already exists; choose a new directory.")
  grids <- lapply(c("GIC", "BIC"), function(ic) {
    grid <- e$build_tuning_grid(result$summary_table, ic)
    grid$recovery_replicates <- NULL
    grid$replicates_by_scenario <- setNames(as.integer(result$group_counts$n_generated),
                                          c("null", "proportional", "correlation"))
    grid$attempts_by_scenario <- setNames(as.integer(result$group_counts$n_attempted),
                                         names(grid$replicates_by_scenario))
    grid$generation_accounting <- result$group_counts
    grid$accounting_scope <- result$accounting_scope
    for (scenario in c("null", "proportional", "integration-rate")) {
      prefix <- if (scenario == "integration-rate") "correlation" else scenario
      rows <- result$summary_table[result$summary_table$IC == ic & result$summary_table$scenario == scenario, ]
      rows <- rows[match(grid$summary_table$setting_id, rows$setting_id), ]
      for (field in c("n_attempted", "n_generated", "n_placement_failed", "n_search_attempted", "n_evaluated")) {
        grid$summary_table[[paste0(prefix, "_", field)]] <- rows[[field]]
      }
    }
    list(schema_version = 1L, grid = grid, policy = e$selection_policy,
         recommendation = e$select_grid(grid, namespace), provenance = result$provenance,
         summary_audit = audit)
  })
  names(grids) <- c("gic", "bic")
  dir.create(destination, recursive = TRUE)
  suite <- result; suite$summaries <- NULL; suite$diagnostics <- NULL
  suite$summary_audit <- audit
  e$atomic_save_rds(suite, file.path(destination, "suite-summary.rds"))
  e$atomic_write_csv(result$summary_table, file.path(destination, "suite-summary.csv"))
  e$atomic_write_csv(result$group_counts, file.path(destination, "generation-accounting.csv"))
  e$atomic_write_csv(result$attempts, file.path(destination, "attempts.csv"))
  e$atomic_write_csv(result$diagnostics, file.path(destination, "search-diagnostics.csv"))
  for (id in names(result$summaries)) e$atomic_save_rds(result$summaries[[id]], file.path(destination, paste0(id, ".rds")))
  for (ic in names(grids)) {
    e$atomic_save_rds(grids[[ic]], file.path(destination, paste0(ic, "-tuning.rds")))
    e$atomic_write_csv(grids[[ic]]$grid$summary_table, file.path(destination, paste0(ic, "-tuning.csv")))
    message(toupper(ic), ": ", grids[[ic]]$recommendation$status)
  }
  invisible(grids)
}

available_tuning_main <- function(args = commandArgs(TRUE)) {
  if (length(args) != 2L) stop("Usage: Rscript summarize-available-tuning-grid.R OUTPUT_ROOT ORIGINAL_STDERR_LOG")
  root <- normalizePath(args[1], mustWork = TRUE)
  log <- normalizePath(args[2], mustWork = TRUE)
  self <- normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))
  Sys.setenv(BIFROST_TUNING_RUNNER_PATH = file.path(dirname(self), "run-vignette-tuning-grid.R"))
  code <- new.env(parent = globalenv()); sys.source(Sys.getenv("BIFROST_TUNING_RUNNER_PATH"), code)
  e <- code$.tuning; e$active_output <- root; e$active_mode <- "full"
  destination <- file.path(root, "full", "summaries-available")
  if (dir.exists(destination)) stop("Summary destination already exists; refusing to overwrite.")
  e$set_single_threaded_math()
  ns <- e$load_development_package()
  for (file in c("provenance.rds", "calibration-template.rds")) {
    if (!file.exists(file.path(root, "full", file))) stop("Missing saved ", file, "; no fitting is permitted.")
  }
  raw_files <- sort(c(file.path(root, "full", c("provenance.rds", "calibration-template.rds")),
    list.files(file.path(root, "full", "replicates"), recursive = TRUE, full.names = TRUE),
    list.files(file.path(root, "full", "checkpoints"), recursive = TRUE, full.names = TRUE)))
  before <- tools::md5sum(raw_files)
  inputs <- e$resolve_data_inputs(NULL, ns)
  template <- e$build_template(inputs, ns) # loads the verified, already-saved template
  provenance <- e$collect_provenance(NULL, NULL, inputs, template, ns)
  groups <- e$build_groups("full"); jobs <- e$build_jobs(groups)
  design <- digest::digest(e$suite_contract(groups, "full"), algo = "sha256")
  result <- collect_available_tuning(e, groups, jobs, root, log, ns, provenance, design)
  stopifnot(identical(before, tools::md5sum(raw_files)))
  audit <- list(created_at = format(Sys.time(), tz = "UTC", usetz = TRUE),
    summarizer_sha256 = digest::digest(file = self, algo = "sha256", serialize = FALSE),
    failure_log = log, failure_log_sha256 = digest::digest(file = log, algo = "sha256", serialize = FALSE),
    input_md5 = before, R = R.version.string)
  write_available_tuning(result, e, ns, destination, audit)
  stopifnot(identical(before, tools::md5sum(raw_files)))
  print(result$group_counts)
  cat("Saved searches:", nrow(result$diagnostics), "\nFitting failures:", sum(result$diagnostics$status == "error"),
      "\nSearches with warnings:", sum(result$diagnostics$warning_count > 0),
      "\nNonfinite IC values:", sum(result$diagnostics$nonfinite_ic),
      "\nRaw files verified unchanged:", length(raw_files), "\nOutput:", destination, "\n")
}

if (sys.nframe() == 0L) available_tuning_main()
